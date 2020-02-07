import numpy

from scitbx.array_family import flex

def scale_map_to_reference(ref_vals, vals, mask_idxs=None):
    """Scale given map to reference over the point_mask"""

    if mask_idxs is not None:
        s_ref_vals = ref_vals.select(mask_idxs)
        s_vals     = vals.select(mask_idxs)
    else:
        s_ref_vals = ref_vals
        s_vals     = vals

    # Fit ref = a*map + b
    # a, b
    scale, offset = numpy.polyfit(x=s_vals, y=s_ref_vals, deg=1)
#    print '\nScale: {}\tOffset: {}'.format(scale, offset)

    return (vals*scale) + offset
