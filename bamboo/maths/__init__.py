
import numpy

def finite(x):
    """Returns finite values of the given array"""
    return numpy.array(x)[numpy.isfinite(x)]

def round_no_fail(a, decimals=0):
    try:    return numpy.round(a, decimals)
    except: return None
