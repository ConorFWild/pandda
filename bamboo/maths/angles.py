import numpy

def rad2deg(a):
    return a*180.0/numpy.pi
def deg2rad(a):
    return a*numpy.pi/180.0

def angle_difference(a1, a2, deg=True, abs_val=False):
    """Calculate the difference between two angles"""

    if deg is False:
        a1 = rad2deg(a1)
        a2 = rad2deg(a2)

    d = (a2-a1+180.0)%360.0-180.0

    if abs_val:
        d = numpy.abs(d)

    if deg is False:
        return deg2rad(d)
    else:
        return d

def normalise_modular_range(value, min, max):
    """Normalises values (e.g. angles) into the range {min, max}"""
    return numpy.mod(value-min, max-min)+min

def find_optimal_angle_window(angles, jump=5):
    """Find the optimal search window for the angles using a window of `jump`"""
    # The optimal shift for the std
    min_shift = min_std = 999
    # Shift the angles by an amount
    for shift in range(0, 360, jump):
        # Reset to {0,360}
        new_angles = numpy.mod(angles-shift, 360)
        # Calculate the std of this shift
        shift_std = numpy.std(new_angles)
        if shift_std < min_std:
            min_std = shift_std
            min_shift = shift
    return min_shift

def mean_std_angles(angles, interval=(-180,180)):
    """Calculate the mean and std of a set of angles. Returns mean in the range {0,360}"""
    # Find the optimal shift to avoid circularity (will fail for large stds)
    opt_shift = find_optimal_angle_window(angles)
    # Calculate the mean and standard deviation with the optimal shift
    ang_mean = numpy.mean(numpy.mod(angles-opt_shift, 360))+opt_shift
    ang_mean = normalise_modular_range(ang_mean, min=interval[0], max=interval[1])
    ang_std  = numpy.std(numpy.mod(angles-opt_shift, 360))
    return ang_mean, ang_std

