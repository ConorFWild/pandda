
import numpy

from scitbx.array_family import flex

def calculate_bdc_subtracted_map(query_map_data, ref_map_data, bdc):
    """Calculate query_map_data - bdc*ref_map_data"""
    assert isinstance(ref_map_data, flex.double), 'ref_map_data must of type flex.double'
    assert isinstance(query_map_data, flex.double), 'query_map_data must of type flex.double'
    assert ref_map_data.all() == query_map_data.all(), 'ref_map_data and query_map_data are not the same size: {!s} != {!s}'.format(ref_map_data.all(), query_map_data.all())
    diff_map_data = query_map_data - ref_map_data*bdc
    return diff_map_data

def calculate_varying_bdc_correlations(ref_map_data, query_map_data, feature_idxs, reference_idxs=None,
                                        min_remain=0.0, max_remain=1.0, bdc_increment=0.01, verbose=False):
    """
    Estimate the background corrections of a feature in "query_map", defined by grid point indices "feature_idxs".
    Reference should be an equivalent map differing only by not having the feature present. An optional mask of reference_idxs can be given.
    """

    # assert isinstance(ref_map_data, flex.double), 'ref_map_data must of type flex.double'
    # assert isinstance(query_map_data, flex.double), 'query_map_data must of type flex.double'
    # assert isinstance(feature_idxs, flex.size_t), 'grid_point_selection must be of type flex.size_t'
    # assert isinstance(reference_idxs, flex.size_t), 'grid_point_selection must be of type flex.size_t'
    assert ref_map_data.all() == query_map_data.all(), 'ref_map_data and query_map_data are not the same size: {!s} != {!s}'.format(ref_map_data.all(), query_map_data.all())

    assert 1.0 >= max_remain > min_remain >= 0.0, 'Min ({}) and Max ({}) remain values are not valid'.format(min_remain, max_remain)

    # Extract the reference map values around the site of interest
    feature_vals_ref = ref_map_data.as_1d().select(feature_idxs)
    # print("Mean reference value is: {}".format(numpy.mean(feature_vals_ref.as_numpy_array())))
    # Extract the reference map values for the correlation correction
    if reference_idxs:
        reference_vals_ref = ref_map_data.as_1d().select(reference_idxs)
    else:
        reference_vals_ref = ref_map_data.as_1d()

    return_vals = []

    # Create list of fractions to subtract
    all_feature_remains = numpy.arange(min_remain, max_remain+bdc_increment, bdc_increment)

    # Iterate through different amounts of reference map subtraction to estimate feature correction
    for i_remain, feature_remain in enumerate(all_feature_remains):
        # Calculate how much of the reference map to take away
        bdc = 1.0 - feature_remain
        # Calculate the background-subtracted difference map
        diff_map_data = calculate_bdc_subtracted_map(ref_map_data=ref_map_data, query_map_data=query_map_data, bdc=bdc)
        # Extract feature idxs values
        feature_vals_diff = diff_map_data.as_1d().select(feature_idxs)
        # print("Feature reference value is: {}".format(numpy.mean(feature_vals_diff.as_numpy_array())))

        # Extract reference idxs values
        if reference_idxs:
            reference_vals_diff = diff_map_data.as_1d().select(reference_idxs)
        else:
            reference_vals_diff = diff_map_data.as_1d()
        # Calculate the correlations to the reference maps
        feature_corr   = flex.linear_correlation(feature_vals_diff,   feature_vals_ref,   subtract_mean=False)
        reference_corr = flex.linear_correlation(reference_vals_diff, reference_vals_ref, subtract_mean=False)

        return_vals.append((feature_remain, feature_corr.coefficient(), reference_corr.coefficient()))

    return zip(*return_vals)

def calculate_maximum_series_discrepancy(labels, series_1, series_2):
    """Calculate the point at which two series are maximally different"""
    assert len(series_1) == len(series_2)
    assert len(series_1) == len(labels)
    diffs = [series_1[i] - series_2[i] for i in range(len(series_1))]
    return labels[diffs.index(max(diffs))]

def estimate_event_background_correction(ref_map_data, query_map_data, feature_idxs, reference_idxs=None, min_remain=0.0, max_remain=1.0, bdc_increment=0.01, method='value', verbose=False):
    """Calculate an estimate for the correction of a defined event"""
    # Check method is valid
    assert method in ['gradient', 'value']
    # Calculate the correlations to the reference map for a range of corrections
    event_remains, event_corrs, global_corrs = calculate_varying_bdc_correlations(
        ref_map_data     = ref_map_data,
        query_map_data   = query_map_data,
        feature_idxs     = feature_idxs,
        reference_idxs   = reference_idxs,
        min_remain       = min_remain,
        max_remain       = max_remain,
        bdc_increment    = bdc_increment,
        verbose          = verbose)

    # Estimate the event corrections (1-bdc)
    if method == 'value':
        # Calculate the maximum discrepancy in the correlation values as an estimation of event correction
        remain_est = calculate_maximum_series_discrepancy(
            labels   = event_remains,
            series_1 = global_corrs,
            series_2 = event_corrs)
    elif method == 'gradient':
        # Calculate the maximum discrepancy in the correlation gradients as an estimation of event correction
        remain_est = calculate_maximum_series_discrepancy(
            labels   = event_remains[1:-1],
            series_1 = [event_corrs[i+1] - event_corrs[i-1] for i in range(1, len(event_corrs)-1)],
            series_2 = [global_corrs[i+1] - global_corrs[i-1] for i in range(1, len(global_corrs)-1)])

    return remain_est


