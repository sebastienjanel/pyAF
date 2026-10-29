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


import numpy

coeff = None
# Force per unit modulus of the tabulated models (sphere, blunted tips),
# see contact_models.unit_force
unit = None

# Models for curve_fit #################################################################################################

def hertz_sphere(delta, delta0, e0):
    y = numpy.zeros(delta.shape)
    for i in range(len(y)):
        if delta[i] <= delta0:
            y[i] = 0
        else:
            y[i] = coeff * e0 * numpy.power((delta[i] - delta0), 1.5)

    return y


def sneddon_cone(delta, delta0, e0):
    y = numpy.zeros(delta.shape)
    for i in range(len(y)):
        if delta[i] <= delta0:
            y[i] = 0
        else:
            y[i] = coeff * e0 * numpy.power((delta[i] - delta0), 2)

    return y


def barber_billings_pyramid(delta, delta0, e0):
    y = numpy.zeros(delta.shape)
    for i in range(len(y)):
        if delta[i] <= delta0:
            y[i] = 0
        else:
            y[i] = coeff * e0 * numpy.power((delta[i] - delta0), 2)

    return y

def flat_punch(delta, delta0, e0):
    y = numpy.zeros(delta.shape)
    for i in range(len(y)):
        if delta[i] <= delta0:
            y[i] = 0
        else:
            y[i] = coeff * e0 * numpy.power((delta[i] - delta0), 1)

    return y

# Models for lmfit #####################################################################################################

def hertz_sphere_lmfit(params, delta, data):
    delta0 = params['delta0_fit'].value
    e0 = params['e0_fit'].value

    y = numpy.zeros(delta.shape)
    for i in range(len(y)):
        if delta[i] <= delta0:
            y[i] = 0
        else:
            y[i] = coeff * e0 * numpy.power((delta[i] - delta0), 1.5)

    return data - y

def sneddon_cone_lmfit(params, delta, data):
    delta0 = params['delta0_fit'].value
    e0 = params['e0_fit'].value

    y = numpy.zeros(delta.shape)
    for i in range(len(y)):
        if delta[i] <= delta0:
            y[i] = 0
        else:
            y[i] = coeff * e0 * numpy.power((delta[i] - delta0), 2)

    return data - y


def barber_billings_pyramid_lmfit(params, delta, data):
    delta0 = params['delta0_fit'].value
    e0 = params['e0_fit'].value

    y = numpy.zeros(delta.shape)
    for i in range(len(y)):
        if delta[i] <= delta0:
            y[i] = 0
        else:
            y[i] = coeff * e0 * numpy.power((delta[i] - delta0), 2)

    return data - y

def flat_punch_lmfit(params, delta, data):
    delta0 = params['delta0_fit'].value
    e0 = params['e0_fit'].value

    y = numpy.zeros(delta.shape)
    for i in range(len(y)):
        if delta[i] <= delta0:
            y[i] = 0
        else:
            y[i] = coeff * e0 * numpy.power((delta[i] - delta0), 1)

    return data - y


def tabulated(delta, delta0, e0):
    """Tabulated contact models (sphere, blunted cone and pyramid)."""
    return e0 * unit(numpy.asarray(delta) - delta0)


def tabulated_lmfit(params, delta, data):
    delta0 = params['delta0_fit'].value
    e0 = params['e0_fit'].value
    return data - tabulated(delta, delta0, e0)
