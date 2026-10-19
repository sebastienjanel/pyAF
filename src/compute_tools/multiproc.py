# Copyright Michka Popoff (2011-2014) michkapopoff@gmail.com
# Copyright Antoine Dujardin (2016-2017) toine.dujardin@gmail.com
# Copyright Sébastien Janel (2024- ) sebastien.janel@cnrs.fr
#
# This software is a computer program whose purpose is to analyze force curves
# recorded with an atomic force microscope.
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software. You can use, modify
# and/ or redistribute the software under the terms of the CeCILL license as
# circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty and the software's author, the holder of the
# economic rights, and the successive licensors have only limited liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading, using, modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean that it is complicated to manipulate, and that also
# therefore means that it is reserved for developers and experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and, more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.

"""Computation of the stiffness and of the work on several processes.

A pool of worker processes is created at the first computation, and reused
for the next ones (starting a process costs about 1 s, the time to import
the computation modules). The processes are started with the "spawn"
method on all the systems (as on Windows and macOS): on Linux the default
"fork" method copies the running application (Qt) into the processes,
which is not safe.

Only the main process reads and writes the HDF5 file. The map is read in
blocks (see misc_tools.get_nbr_blocks), and each block is divided in small
parts, sent to the processes. The results are written when all the parts
of the block are computed. The parts are small enough for the processes
to share the work evenly (the cores do not all have the same speed).
"""

import atexit
import concurrent.futures
import math
import multiprocessing
import sys

import numpy
import tables

from .. import shared
from .. import widgets_list
from ..tools import misc_tools
from .workers import compute_columns

# Pool of processes, created at the first computation
_EXECUTOR = None
_EXECUTOR_WORKERS = 0

# Number of parts per process (more parts: better sharing of the work
# between the processes, but more communication)
PARTS_PER_PROCESS = 4

# Results of the curves, written in the HDF5 file (index in the result of
# the computation function of a curve)
STIFFNESS_RESULTS = {"pocs_indices": 2, "pocs_real": 3, "fits_poc": 4,
                     "topography": 6, "indentation": 8, "fit_quality": 9}
WORK_RESULTS = {"jocs1_indices": 2, "jocs2_indices": 3, "jocs1_real": 4,
                "jocs2_real": 5, "fits_joc": 6, "work": 7,
                "rupture_force1": 8}


def get_executor(nbr_workers):
    """Returns the pool of processes, created if needed."""
    global _EXECUTOR, _EXECUTOR_WORKERS

    if _EXECUTOR is not None and _EXECUTOR_WORKERS != nbr_workers:
        shutdown_pool()

    if _EXECUTOR is None:
        context = multiprocessing.get_context("spawn")
        _EXECUTOR = concurrent.futures.ProcessPoolExecutor(
            max_workers=nbr_workers, mp_context=context)
        _EXECUTOR_WORKERS = nbr_workers

    return _EXECUTOR


def shutdown_pool(wait=True):
    """Stops the processes of the pool (at the end of the application, or
    after an error)."""
    global _EXECUTOR, _EXECUTOR_WORKERS

    if _EXECUTOR is not None:
        _EXECUTOR.shutdown(wait=wait, cancel_futures=True)
        _EXECUTOR = None
        _EXECUTOR_WORKERS = 0


atexit.register(shutdown_pool)


def compute(mode, calc_id, params, path):
    """Computes the stiffness ("stiffness" mode) or the work ("work" mode)
    of the curves of a file.

    params contains the parameters of the computation (see the stiffness
    and work main.py files), path is the path of the HDF5 file (closed by
    the caller during the computation).

    Returns the list of the errors of the curves: ["fit_error" or
    "smoothing_error", "[col, row]"].
    """
    errors = []

    dt = shared.exp.list[calc_id]
    tilt_corr = dt.tilt_applied
    nbr_cols = dt.nbr_pixels_x
    nbr_rows = dt.nbr_pixels_y

    # The approach curves are needed for the stiffness, the retraction
    # curves for the work, and both for the tilt correction
    read_app = mode == "stiffness" or tilt_corr
    read_ret = mode == "work" or tilt_corr

    # Arrays of the whole map, sent by parts to the processes
    discarded = numpy.asarray(params["discarded_curves"])
    piezo_image = params.get("piezo_image", None)
    if isinstance(piezo_image, bool):
        piezo_image = None
    common = {key: value for key, value in params.items()
              if key not in ("discarded_curves", "piezo_image")}

    nbr_workers = max(1, min(shared.exp.cores, nbr_cols * nbr_rows))
    if sys.platform == "win32":
        # Maximum number of processes of a pool on Windows
        nbr_workers = min(nbr_workers, 61)
    executor = get_executor(nbr_workers) if nbr_workers > 1 else None

    h5_file = tables.open_file(path, "r+")
    try:
        st = "/data/_" + str(calc_id)
        if read_app:
            app_node = h5_file.get_node(st + "/curves/approach")
            app_pos_node = h5_file.get_node(
                st + "/positions/approach_positions")
        if read_ret:
            ret_node = h5_file.get_node(st + "/curves/retraction")
            ret_pos_node = h5_file.get_node(
                st + "/positions/retraction_positions")

        if mode == "stiffness":
            result_nodes = {name: h5_file.get_node(st + "/results", name)
                            for name in STIFFNESS_RESULTS}
            stiffness_node = h5_file.get_node(st + "/results",
                                              "stiffness_array")
        else:
            result_nodes = {name: h5_file.get_node(st + "/results", name)
                            for name in WORK_RESULTS}
        indices = STIFFNESS_RESULTS if mode == "stiffness" else WORK_RESULTS

        # Columns of the blocks
        nbr_blocks = misc_tools.get_nbr_blocks(nbr_cols, nbr_rows)
        cols_in_block = nbr_cols // nbr_blocks
        starts = [i * cols_in_block for i in range(nbr_blocks)]
        ends = starts[1:] + [nbr_cols]

        for block_start, block_end in zip(starts, ends):
            # Read the curves of the block
            app = app_pos = ret = ret_pos = None
            if read_app:
                app = app_node[block_start:block_end, :, :]
                app_pos = app_pos_node[block_start:block_end, :, :]
            if read_ret:
                ret = ret_node[block_start:block_end, :, :]
                ret_pos = ret_pos_node[block_start:block_end, :, :]

            # Results of the block, written at the end of the block
            buffers = {name: node[block_start:block_end]
                       for name, node in result_nodes.items()}
            stiffness = [None] * ((block_end - block_start) * nbr_rows)

            tasks = []
            for c0, c1, r0, r1 in split_block(
                    block_end - block_start, nbr_rows, nbr_workers):
                tasks.append({
                    "mode": mode,
                    "params": common,
                    "first_col": block_start + c0,
                    "first_row": r0,
                    "app": None if app is None else app[c0:c1, r0:r1],
                    "app_positions":
                        None if app is None else app_pos[c0:c1, r0:r1],
                    "ret": None if ret is None else ret[c0:c1, r0:r1],
                    "ret_positions":
                        None if ret is None else ret_pos[c0:c1, r0:r1],
                    "discarded": discarded[block_start + c0:
                                           block_start + c1, r0:r1],
                    "piezo_image":
                        None if piezo_image is None or app is None else
                        piezo_image[block_start + c0:block_start + c1,
                                    r0:r1]})

            for results in run_tasks(executor, tasks):
                for status, position, result in results:
                    if status != "continue":
                        errors.append([status, position])
                    if result is not None:
                        col = result[0] - block_start
                        row = result[1]
                        for name, index in indices.items():
                            buffers[name][col, row] = result[index]
                        if mode == "stiffness":
                            stiffness[col * nbr_rows + row] = result[5]
                widgets_list.widget_progressbar.update(len(results))

            widgets_list.widget_progressbar.set_label("Saving results")
            for name, node in result_nodes.items():
                node[block_start:block_end] = buffers[name]
            if mode == "stiffness":
                for values in stiffness:
                    stiffness_node.append(
                        [] if values is None else values)
            widgets_list.widget_progressbar.set_label("Computing")
    finally:
        h5_file.close()

    return errors


def split_block(nbr_cols, nbr_rows, nbr_workers):
    """Divides a block of curves in parts (col start, col end, row start,
    row end), for the processes.

    The parts contain whole columns if possible, else parts of columns
    (e.g. for a 4x5000 map).
    """
    nbr_curves = nbr_cols * nbr_rows
    if nbr_workers == 1:
        # Computed in the main process: parts only used for the progressbar
        size = max(1, nbr_curves // 100)
    else:
        size = max(1, math.ceil(nbr_curves / (nbr_workers *
                                              PARTS_PER_PROCESS)))

    parts = []
    if size >= nbr_rows:
        cols_per_part = size // nbr_rows
        for c0 in range(0, nbr_cols, cols_per_part):
            parts.append((c0, min(c0 + cols_per_part, nbr_cols),
                          0, nbr_rows))
    else:
        for c0 in range(nbr_cols):
            for r0 in range(0, nbr_rows, size):
                parts.append((c0, c0 + 1, r0, min(r0 + size, nbr_rows)))
    return parts


def run_tasks(executor, tasks):
    """Computes the tasks, and yields their results as they are finished.

    Without executor, the tasks are computed in the main process. If a
    task fails, the pool is stopped (it will be recreated at the next
    computation) and the error is raised.
    """
    if executor is None:
        for task in tasks:
            yield compute_columns(task)
        return

    futures = [executor.submit(compute_columns, task) for task in tasks]
    try:
        for future in concurrent.futures.as_completed(futures):
            yield future.result()
    except BaseException:
        # Error in a process (or process stopped): stop the computation
        for future in futures:
            future.cancel()
        shutdown_pool(wait=False)
        raise
