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

"""
This file lets you share data across the PYAF project (and the plugins).

This technique is explained in the python docs :
http://docs.python.org/2/faq/programming.html\
#how-do-i-share-global-variables-across-modules

The most important variable here is the **exp** variable. It will contain the
experiment object (from experiment.py), once PYAF has loaded a file. In exp,
there is a variable called list, containing **data** objets. Each data object
corresponds to a loaded file. You can then access the data of the i'th file by
using the following statement in your code :
.. code-block:: python

   import shared
   data = shared.exp.list[i]

The other variables contain data once computations have been done. They are
filled by the **tools/results_sorting.py** file.

The last variables, files_removed, is set to true once the remove button has
been used to delete a file. It will tell the saving algorithm in
load_and_save/save.py to trigger a special procedure to refactor the file
containing the data. It will then free the disk space previously used by the
removed file(s).

The layer_list contains all the layers for the 3D view.

The clipboard_layer_positions contains positions of layers (x, y, z) which can
be saved in the VTK layer.

VTK_first : see description in widgets_vtk/main.py file

"""

exp = None
single_data = None
groups_data = None
conditions_data = None
statistics_data = None  # Data frame containing data for analysis

single_values = None
groups_values = None
conditions_values = None

single_factors = None
groups_factors = None
conditions_factors = None

single_frequencies = None
groups_frequencies = None
conditions_frequencies = None

single_pdfs_x = None
single_pdfs_y = None
groups_pdfs_x = None
groups_pdfs_y = None
conditions_pdfs_x = None
conditions_pdfs_y = None

files_removed = False

layer_list = None

clipboard_layer_positions = None
VTK_first = True

# Statistical tests available
stat_tests = ["Two-sample T-test", "Paired T-test"]

# Statistial plots available
stat_plots = ["Boxplot", "Swarmplot", "Violinplot", "Superplot"]

significance = {"***": 0.001, "**": 0.01, "*": 0.05}

# Some basic colors used for the plots
# Corresponds to "Blue", "Red", "Green", "Orange", "Black"
colors_list = [
    [38/255.0, 0.0, 255/255.0, 1.0],
    [252/255.0, 5/255.0, 22/255.0, 1.0],
    [0.0, 255/255.0, 4/255.0, 1.0],
    [250/255.0, 187/255.0, 40/255.0, 1.0],
    [0.0, 0.0, 0.0, 1.0]]


# A list of ids to compute. Allows to recompute a part of data from a plugin
# Leave to None for normal usage. Else you can use it like this from your
# plugin : shared.force_calc_list = [5, 6, 7]
force_calc_list = None

# Element which is zoomed in on the curve in the compute tab
# Can be either an event on the retraction curve, the jump of contact,
# or the point of contact on the approach curve.
zoomed_in_element = 0

# Flag to know if an error has been triggered somewhere. In this case the
# flag is set to true, and when saving the file the user is told that his
# data may be corrupted
error_triggered = False
