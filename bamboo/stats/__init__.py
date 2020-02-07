import numpy

def modified_z_scores(array):
    """Z-scores calculated using deviations from the median rather than the mean"""
    array = numpy.array(array)
    # Calculate deviations from the median
    medn = numpy.median(array)
    devs = array-medn
    # Calculate median of deviations
    mdev = numpy.median(numpy.abs(devs))
    return 0.6745*devs/mdev

def quartiles(array):
    return numpy.percentile(array, [25,75])

def iqr(array):
    q1, q3 = quartiles(array)
    return q3-q1

def iqr_outliers(array):
    array = numpy.array(array)
    q1,q3 = quartiles(array)
    iqr_a = q3-q1
    return (array<(q1-1.5*iqr_a)) | (array>(q3+1.5*iqr_a))
