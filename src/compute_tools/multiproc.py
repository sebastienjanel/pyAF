# Copyright Michka Popoff (2011-2014) michkapopoff@gmail.com
# Copyright Antoine Dujardin (2016-2017) toine.dujardin@gmail.com
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

"""Tools for multiprocessing."""

import traceback
import numpy
import tables
import multiprocessing
import queue
from .. import shared
from .. import widgets_list
from ..tools import misc_tools
from .stiffness.stiffness_compute import main_compute_stiff
from .stiffness.stiffness_fit import main_fit_stiffness
from .work_and_rupture_force.work_and_rupture_force_compute \
    import main_compute_work


# pylint: disable=E1103

class FileAccess(multiprocessing.Process):
    """All access to the file goes through a single instance of this class.

    It contains several queues that are used to communicate with other
    processes.
    The read_queue is used for requests to read data from the HDF5 file.
    A list of result_queues is used to send data back to client processes.
    The write_queue is used for requests to modify the HDF5 file.
    One end of a pipe, shutdown is used to signal the process to terminate.

    Note : the exceptions in the child processes are fetched following the
    ideas provided by : http://stackoverflow.com/a/19951172
    """

    def __init__(self, path, queues, shutdown, data_id, mode, sizes,
                 tilt_corr):
        self.path = path
        self.read_q = queues[0]
        self.result_q = queues[1]
        self.write_q = queues[2]
        self.shutdown = shutdown
        self.use_block = False
        self.unique_id = data_id
        self.mode = mode
        self.sizes = sizes
        self.tilt_corr = tilt_corr

        super().__init__()

    def run(self):
        """Called when the Process is started.

        Will call the safe_run method. We catch the exceptions here to send
        them back to the main thread if needed.
        """
        # pylint: disable=W0703

        try:
            self.safe_run()
        except Exception:
            self.fetch_error(traceback.format_exc())
        return

    def fetch_error(self, e):
        """Catches the exceptions and sends them back to the main thread.

        We use the result queue for this. Use the first one for this as we
        can not know the proc_num here.
        """
        self.result_q[0].put(["exception", e])

    def safe_run(self):
        """Reads and writes the data from the hdf5 file."""
        # Open the hdf5 file
        h5_file = tables.open_file(self.path, "r+")
        st = "/data/_" + self.unique_id
        tf = h5_file

        nbr_rows = self.sizes[0]
        nbr_cols_in_block = self.sizes[1]
        shift_block = self.sizes[2]
        nbr_curves_in_block = nbr_rows * nbr_cols_in_block

        cols_start = shift_block
        cols_end = shift_block + nbr_cols_in_block

        app_curves = None
        ret_curves = None
        app_positions = None
        ret_positions = None

        # Get the nodes
        # Use numpy arrays to read the data in one block
        if self.mode == "stiffness" or self.tilt_corr:
            a1 = tf.get_node(st + "/curves/approach")
            a3 = tf.get_node(st + "/positions/approach_positions")
            app_curves = numpy.empty([nbr_cols_in_block, nbr_rows], list)
            app_curves = a1[cols_start:cols_end, :, :]
            app_positions = numpy.empty([nbr_cols_in_block, nbr_rows], list)
            app_positions = a3[cols_start:cols_end, :, :]

            pocs_indices = tf.get_node(st + "/results", "pocs_indices")
            pocs_real = tf.get_node(st + "/results", "pocs_real")
            topography = tf.get_node(st + "/results", "topography")
            fits_poc = tf.get_node(st + "/results", "fits_poc")
            indentation = tf.get_node(st + "/results", "indentation")
            stiffness_file = tf.get_node(st + "/results", "stiffness_array")
            stiffness_numpy = numpy.empty([nbr_curves_in_block], list)

        if self.mode == "work" or self.tilt_corr:
            a2 = tf.get_node(st + "/curves/retraction")
            a4 = tf.get_node(st + "/positions/retraction_positions")
            ret_curves = numpy.empty([nbr_cols_in_block, nbr_rows], list)
            ret_curves = a2[cols_start:cols_end, :, :]
            ret_positions = numpy.empty([nbr_cols_in_block, nbr_rows], list)
            ret_positions = a4[cols_start:cols_end, :, :]

            jocs1_indices = tf.get_node(st + "/results", "jocs1_indices")
            jocs2_indices = tf.get_node(st + "/results", "jocs2_indices")
            jocs1_real = tf.get_node(st + "/results", "jocs1_real")
            jocs2_real = tf.get_node(st + "/results", "jocs2_real")
            fits_joc = tf.get_node(st + "/results", "fits_joc")
            work = tf.get_node(st + "/results", "work")
            rupture_force1 = tf.get_node(st + "/results", "rupture_force1")
        count = 0
        another_loop = True

        while another_loop:
            # Check if the process has received the shutdown signal.
            if self.shutdown.poll():
                another_loop = False
                if count < nbr_curves_in_block:
                    # Be sure to quit the loop when the results are all saved
                    another_loop = True

            # Check for any data requests in the read_queue.
            try:
                proc_num = self.read_q.get(self.use_block)

                position = [app_positions, ret_positions]

                result_queue = self.result_q[proc_num]
                result_queue.put(([app_curves, ret_curves], position))

            except queue.Empty:
                pass

            # Check for any write requests in the write_queue.
            try:
                proc_num, res = self.write_q.get(self.use_block)

                if self.mode == "stiffness":
                    col = res[0]
                    row = res[1]

                    pocs_indices[col, row] = res[2]
                    pocs_real[col, row] = res[3]
                    fits_poc[col, row] = res[4]
                    topography[col, row] = res[6]
                    indentation[col, row] = res[8]

                    # Shift for col start in block
                    col = col - cols_start

                    stiffness_numpy[col * nbr_rows + row] = res[5]

                else:
                    jocs1_indices[res[0], res[1]] = res[2]
                    jocs2_indices[res[0], res[1]] = res[3]
                    jocs1_real[res[0], res[1]] = res[4]
                    jocs2_real[res[0], res[1]] = res[5]
                    fits_joc[res[0], res[1]] = res[6]
                    work[res[0], res[1]] = res[7]
                    rupture_force1[res[0], res[1]] = res[8]

                count += 1

            except queue.Empty:
                pass

        # Write to disk
        if self.mode == "stiffness":
            for index in range(nbr_curves_in_block):
                stiffness_file.append(stiffness_numpy[index])

        # Close the HDF5 file before shutting down
        h5_file.close()


def make_queues(
        num_processors,
        path,
        unique_id,
        mode,
        sizes,
        tilt_corr):
    """
    This function starts the FileAccess class instance and sets up all the
    queues used to communicate with it.

    Read queue : read the data from the hdf5 file.
    Write queue : write the data to the hdf5 file.
    Result queue(s) : pass data from read queue to write queue.
    Progress queue : updates the progressbar in the main thread.

    """

    # Create the queues
    read_queue = multiprocessing.Queue()
    write_queue = multiprocessing.Queue()
    progress_queue = multiprocessing.Queue()
    result_queue = [multiprocessing.Queue() for _ in range(num_processors)]

    # Signals for closing of the multiprocessing
    shutdown_recv, shutdown_send = multiprocessing.Pipe(False)

    queues = [read_queue, result_queue, write_queue, progress_queue]

    # File Access class
    file_access = FileAccess(
        path, queues, shutdown_recv, str(unique_id), mode, sizes, tilt_corr)
    file_access.start()

    return queues, shutdown_send, file_access


def compute(mode, calc_id, params, path):
    """Function called to start the computation.

    Will create the queus and go through the data for the computations.
    The mode defines the type if computation to do.
    """
    list_of_smoothing_errors = []

    dt = shared.exp.list[calc_id]
    tilt_corr = dt.tilt_applied

    nbr_pixels = dt.nbr_pixels_x * dt.nbr_pixels_y

    # Get the number of processes to spawn
    nbr_processes = misc_tools.get_processes(nbr_pixels)

    # Get the number of memory blocks to use
    nbr_of_blocks = misc_tools.get_nbr_blocks(dt.nbr_pixels_x, dt.nbr_pixels_y)
    total_nbr_cols = dt.nbr_pixels_x
    cols_in_block = total_nbr_cols // nbr_of_blocks
    last_block_nbr_cols = 0

    if total_nbr_cols % nbr_of_blocks != 0:
        last_block_nbr_cols = total_nbr_cols - \
                              (nbr_of_blocks - 1) * cols_in_block
    cols_sizes = []
    for block in range(nbr_of_blocks):
        cols_sizes.append(cols_in_block)
    if last_block_nbr_cols != 0:
        cols_sizes[-1] = last_block_nbr_cols

    for block in range(nbr_of_blocks):
        # Prepare cols
        nbr_cols_in_block = cols_sizes[block]

        # Check if there are too many cores for the given block size
        # In this case, reduce the number of processes
        # Should be a very rare case (can happen with a 4x5000 file)
        # Happens also with single curves
        if nbr_processes > nbr_cols_in_block:
            nbr_processes = nbr_cols_in_block
            if nbr_processes == 0:
                # Single curve case, use only one process
                # nbr_processes == 0 because nbr_cols_in_block = 1
                nbr_processes = 1

        # The last proc can contain less or more curves, in the case where
        # nbr_cols_in_block % nbr_processes != 0
        procs_nbr_cols = []
        modulo_cols_in_proc = nbr_cols_in_block % nbr_processes
        nbr = nbr_cols_in_block // nbr_processes
        for i in range(nbr_processes):
            procs_nbr_cols.append(nbr)

        # Example: 10x10 file with 4 processes (4 cores).
        # The algorithm will ask for 1 block, divided by 4,
        # so that nbr_cols_in_block = 10
        # In this case: modulo_cols_in_proc = 10 % 4 is != 0
        # and nbr = 2, so that procs_nbr_cols = [2, 2, 2, 2],
        # but we need [2, 2, 2, 4], which we will do below
        if modulo_cols_in_proc != 0:
            total = 0
            for i in range(nbr_processes):
                total += nbr
            # How many cols are missing so that we can add them
            diff = nbr_cols_in_block - total
            # Update the last position
            procs_nbr_cols[-1] += diff

        shift_block = 0
        for i in range(0, block, 1):
            shift_block += cols_sizes[i]

        # Create the queues
        sizes = [dt.nbr_pixels_y, nbr_cols_in_block, shift_block]
        queues, shutdown_send, file_access = make_queues(nbr_processes, path,
                                                         calc_id, mode,
                                                         sizes, tilt_corr)

        # Create the computations objects (one per core)
        processors = []
        for i in range(nbr_processes):
            shift_col = i * procs_nbr_cols[i - 1]
            processor = DataProcessor(mode, queues, nbr_processes, i,
                                      params, nbr_cols_in_block,
                                      procs_nbr_cols[i],
                                      shift_col,
                                      shift_block)
            processors.append(processor)

        # Start all the DataProcessor instances
        for processor in processors:
            processor.start()

        # Update the progressbar
        cont = True
        count = 0
        while cont:
            try:
                progress_queue = queues[3]
                res = progress_queue.get()
                # Get the progress. Depending on the first argument, store the
                # errors, continue or raise exception
                if res[0] == "smoothing_error":
                    # Store the error, the excecution will continue
                    list_of_smoothing_errors.append(res)
                elif res[0] == "exception":
                    # No cleanup for the moment, were screwed here ...
                    raise Exception(res[1])
                # Continue (update progress bar)
                widgets_list.widget_progressbar.update()
                if count + 1 == nbr_cols_in_block * dt.nbr_pixels_y:
                    cont = False
                count += 1
            except queue.Empty:
                pass
            except IOError:
                cont = False

        # Wait for all DataProcessor instances to finish
        for processor in processors:
            processor.join()

        widgets_list.widget_progressbar.set_label("Saving results ...")

        # Shut down the FileAccess instance
        shutdown_send.send(0)

        # Be sure the file_acces Process is closed and the data written
        # to the disk
        file_access.join()

        widgets_list.widget_progressbar.set_label("Getting stiffness")

    return list_of_smoothing_errors


class DataProcessor(multiprocessing.Process):
    """Processes the data (one per core).

    Computes the stiffness.
    """

    def __init__(self, mode, queues, nbr_procs, proc_num, dt,
                 nbr_cols_in_block, nbr_cols_in_proc, shift_col, shift_block):
        self.read_q = queues[0]
        self.result_queue = queues[1][proc_num]
        self.write_q = queues[2]
        self.progress_queue = queues[3]
        self.proc_num = proc_num
        self.nbr_procs = nbr_procs
        self.nbr_pixels_x = dt["nbr_pixels_x"]
        self.nbr_pixels_y = dt["nbr_pixels_y"]
        self.dt = dt
        self.nbr_cols_in_proc = nbr_cols_in_proc
        self.nbr_cols_in_block = nbr_cols_in_block
        self.shift_col = shift_col
        self.shift_block = shift_block
        self.mode = mode

        super().__init__()

    def run(self):
        """Called when the Process is started.

        Will call the safe_run method. We catch the exceptions here to send
        them back to the main thread if needed.
        """
        # pylint: disable=W0703

        try:
            self.safe_run()
        except Exception:
            self.fetch_error(traceback.format_exc())
        return

    def fetch_error(self, e):
        """Catches the exceptions and sends them back to the main thread."""
        self.progress_queue.put(["exception", e])

    def safe_run(self):
        """Actual part were the computing is done.

        Will go through a part of the data and compute stuff.
        """
        nbr_rows = self.nbr_pixels_y

        nbr_cols_in_proc = self.nbr_cols_in_proc
        nbr_of_curves_in_proc = nbr_cols_in_proc * nbr_rows

        # Shift col
        row = 0
        col = self.shift_block + self.shift_col
        limit = self.shift_block + self.shift_col

        # Read the data
        self.read_q.put((self.proc_num))
        # Get the data back from the read queue, will be written in the
        # results queue
        res = self.result_queue.get()
        if res[0] == "exception":
            # No cleanup here, this should be a big bug in the child process
            # so we are screwed anyway ...
            raise Exception(res[1])

        dt = {
            "approach_position": None,
            "retraction_position": None,
            "is_discarded": None,
            "piezo_image": None,
            "tilt_limit_1": self.dt["tilt_limit_1"],
            "tilt_limit_2": self.dt["tilt_limit_2"],
            "deflection_sensitivity": self.dt["deflection_sensitivity"],
            "spring_constant": self.dt["spring_constant"],
            "tilt_applied": self.dt["tilt_applied"],
            "sg_smoothing_enabled": self.dt["sg_smoothing_enabled"],
            "sg_smoothing_order": self.dt["sg_smoothing_order"],
            "sg_smoothing_width": self.dt["sg_smoothing_width"],
            "sg_smoothing_uniform": self.dt["sg_smoothing_uniform"],
            "stretch_applied_app": self.dt["stretch_applied_app"],
            "stretch_applied_ret": self.dt["stretch_applied_ret"],
            "stretch_app_lim1": self.dt["stretch_app_lim1"],
            "stretch_app_lim2": self.dt["stretch_app_lim2"],
            "stretch_len_app": self.dt["stretch_len_app"],
            "stretch_ret_lim1": self.dt["stretch_ret_lim1"],
            "stretch_ret_lim2": self.dt["stretch_ret_lim2"],
            "stretch_len_ret": self.dt["stretch_len_ret"]}

        if self.mode == "stiffness":
            fit_params = {
                "poc_skip_start": self.dt["poc_skip_start"],
                "poc_fit_length": self.dt["poc_fit_length"],
                "poc_refit_option": self.dt["poc_refit_option"],
                "poc_noise_multiplicator": self.dt["poc_noise_multiplicator"],
                "poc_refit_times": self.dt["poc_refit_times"],
                "trig_threshold": self.dt["trig_threshold"],
                "poisson_ratio": self.dt["poisson_ratio"],
                "tip_radius": self.dt["tip_radius"],
                "fit_range_type": self.dt["fit_range_type"]}

        elif self.mode == "work":
            fit_params = {
                "joc_skip_start": self.dt["joc_skip_start"],
                "joc_fit_length": self.dt["joc_fit_length"],
                "joc_refit_option": self.dt["joc_refit_option"],
                "joc_noise_multiplicator": self.dt["joc_noise_multiplicator"],
                "joc_refit_times": self.dt["joc_refit_times"]}

        if self.mode == "stiffness":
            o_params = {
                "indentation_start": self.dt["indentation_start"],
                "indentation_stop": self.dt["indentation_stop"],
                "indentation_step": self.dt["indentation_step"],
                "force_start": self.dt["force_start"],
                "force_stop": self.dt["force_stop"],
                "model_selected": self.dt["model_selected"],
                "tomography": self.dt["tomography"],
                "coeff": self.dt["coeff"]
            }

        curves = res[0]

        # Go trough a part of the data
        for _ in range(nbr_of_curves_in_proc):
            dt["is_discarded"] = self.dt["discarded_curves"][col][row]

            # Get the curves
            if curves[0] is not None:
                app = curves[0][col - self.shift_block][row]
                new_col_pos = col - self.shift_block
                dt["approach_position"] = res[1][0][new_col_pos][row]

                # Get the piezo image. In case of tilt correction during
                # work and rupture force computation curves[0] is not empty
                # In this case we have to check if the piezo image exists,
                # if not we do not read it (because we don't need it)
                # Note : I am comparing the types here, because if I compare
                # the array to False, I would need to specify which values to
                # compare in the array (any or all).
                # Case 1 : no piezo image : False != False => False
                # Case 2 : piezo image : type(array) != bool => True
                tp_false = type(False)
                if not isinstance(self.dt.get("piezo_image", False), tp_false):
                    dt["piezo_image"] = self.dt["piezo_image"][col][row]
            else:
                app = None

            if curves[1] is not None:
                ret = curves[1][col - self.shift_block][row]
                new_col_pos = col - self.shift_block
                dt["retraction_position"] = res[1][1][new_col_pos][row]
            else:
                ret = None

            # Compute the stiffness
            if self.mode == "stiffness":
                if self.dt["perform_fit"]:
                    result = main_fit_stiffness(col, row,
                                                app, ret,
                                                dt, fit_params, o_params)
                else:
                    result = main_compute_stiff(col, row,
                                                app, ret,
                                                dt, fit_params, o_params)
            elif self.mode == "work":
                result = main_compute_work(col, row,
                                           app, ret,
                                           dt, fit_params)

            # Write the data
            # print(result)
            self.write_q.put((self.proc_num, result))

            # Fetch smoothing errors and update progressbar
            if result is None:
                # Just confirm that this loop is over
                self.progress_queue.put((["continue", 1]))
            elif self.mode == "stiffness" and result[7]:
                val = "[" + str(col + 1) + ", " + str(row + 1) + "]"
                self.progress_queue.put((["smoothing_error", val]))
            else:
                # Just confirm that this loop is over
                self.progress_queue.put((["continue", 1]))

            # Update the position
            col = col + 1
            if col == limit + nbr_cols_in_proc:
                col = limit
                row += 1
