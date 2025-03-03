import numpy

coeff = None

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


def bilodeau_pyramid(delta, delta0, e0):
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


def bilodeau_pyramid_lmfit(params, delta, data):
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

"""
def hertz_sphere_lmfit(params, delta, data):
    delta0 = params['delta0_fit'].value
    e0 = params['e0_fit'].value
    y = hertz_sphere(delta, delta0, e0)
    return data - y
    
def sneddon_cone_lmfit(params, delta, data):
    delta0 = params['delta0_fit'].value
    e0 = params['e0_fit'].value
    y = sneddon_cone(delta, delta0, e0)
    return data - y

def bilodeau_pyramid_lmfit(params, delta, data):
    delta0 = params['delta0_fit'].value
    e0 = params['e0_fit'].value
    y = bilodeau_pyramid(delta, delta0, e0)
    return data - y
"""