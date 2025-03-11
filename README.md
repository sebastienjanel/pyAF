#pyAF : An AFM Force Curve Analysis Software

pyAF is a tool to analyze force curves & force maps acquired by atomic force microscopes (Bruker and JPK).
Developed at CMPI lab (https://cmpi.cnrs.fr).

It computes:
- Elasticity with several models and bottom effect correction (BEC)
- Stiffness tomography
- Work of detachment and rupture force
- Unbinding events and loading rates
- Zero-force topography
- Regions of interest for substrate definition and stats
- Curve tilt correction and smoothing
- 3D map rendering
- Text export (to R), Statistics (UNPOLISHED)
- Multi files and multiprocessing

This software is governed by the CeCILL license under French law and abiding by the rules of distribution of free software. You can use, modify and/or redistribute the software under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL "http://www.cecill.info".

Main coding made by Michka Popoff (https://theses.fr/api/v1/document/2014LIL10220).
Modified by Antoine Dujardin, Simone Bovio, Javier Lopez-Alonso, Nuno Duarte, Sébastien Janel.

When using the program, cite the original article : Stiffness tomography of eukaryotic intracellular compartments by atomic force microscopy (Nanoscale, 2019,11, 10320-10328).
DOI	https://doi.org/10.1039/C8NR08955H

DISCLAIMER: pyAF is a lab-developed software, it is not a commercial software. Check as much as possible results provided vs other software (e.g. JPK DP) to rule out any no supported edges cases. 