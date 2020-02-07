
from scitbx.array_family import flex
from scitbx.math import gamma_complete, gamma_incomplete
from scitbx.math import distributions

def test_significance_of_group_of_z_values(vals):
    """Test the significance of a group of z-values - assumes sum-squared is distributed as a chi-squared"""
    sum_sq = flex.sum_sq(flex.double(vals))
    n = len(vals)
    return chi_squared_test(sum_sq, n)

def chi_squared_test(sum_sq, n, max_iter=500):
    assert n < 344, 'Cannot test more than 343 values due to method of normalisig gamma function!'
    return 1.0-gamma_incomplete(a=n/2.0, x=sum_sq/2.0, max_iterations=max_iter)

def convert_pvalue_to_zscore(pval, two_tailed=True):
    """Convert a p-value to a z-score for a standard normal N(0,1)"""
    # If two-tailed test, need to halve the p-value
    if two_tailed:
        pval = pval/2.0
    # Create normal distribution to convert - N(0,1)
    nrm = distributions.normal_distribution()
    # Calculate the probability quantile (z-score) corresponding to 1-pval
    try:
        zsco = nrm.quantile(1.0-pval)
    except RuntimeError:
        # pval too small - return default (8.0 is max calculable)
        if pval < 0.00000000000000006:
            zsco = 8.2
        else:
            raise
    return zsco
