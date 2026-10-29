#Development notes

##pyAF
*****************************
Version 1: Michka Popoff, PhD work
-> stiffness tomography, adhesion, events, GUI, 3D, histogram stats

Version 2: Antoine Dujardin, PhD work

Version 3: Javier Lopez-Alonso, PhD work
-> migration to python3 and fitting algorythm
-> export to r for data plotting

Version 4: Sébastien Janel and Nuno Duarte
-> Fixes for compatibilities with new libraries
-> Box plots
-> New fitting and BEC models

Version 5: Sébastien Janel
→ Updated python and librairies
→ Revamped interface (toggles, selectors, theme, autoscale, LUT, plotting)
→ Additional elasticity models
→ Updated Bottom effect corrections factors (Garcia)
→ New Ting viscoelasticity fit (Garcia)
→ New Microrheology (Alcaraz)
→ New Breakthrough analysis
→ Complete new statistics: violin plots + superplots
→ Revamped result export with JSON

##Checking the GUI

tests/check_gui.py loads your files, computes everything on them and uses
every control of the Data, Compute and Results tabs, before and after the
computations. It reports the errors with the action which triggered them:

    .venv/bin/python tests/check_gui.py file1 [file2 ...]

Use files of different types together (e.g. a JPK force map, a JPK QI map
and a Bruker file). It takes about 15 minutes for 3 files, and does not
modify your pyAF preferences.
