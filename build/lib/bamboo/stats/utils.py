import numpy

from libtbx.math_utils import iceil, ifloor

def interpolate(x,p1,v1,p2,v2):
    """Interpolate the values of p1 and p2 to the point at x"""
    if x==p1: return v1
    if x==p2: return v2
    return v2*(x-p1)/(p2-p1) + v1*(p2-x)/(p2-p1)

def calculate_minimum_redundancy(unsampled_size, max_sample_size):
    """Calculate the redundancy required to down-sample a vector of size unsampled_size to size max_sample_size"""
    min_redundancy = int(1+(unsampled_size-1)/max_sample_size)
    return min_redundancy

def resample_ordered_list_of_values(vals, redundancy=8):
    """resample a list of values with interpolation"""
    # Number of vals given
    num_inp_vals = len(vals)
    # Number of vals to be returned
    num_samp_vals = int(1+(num_inp_vals-1)/redundancy)
    # Sort in descending order
    ordered_vals = sorted(vals, reverse=True)
    sampled_vals = []

    if num_samp_vals==1:
        return [ordered_vals[0]]
    else:
        sample_dist = (num_inp_vals-1)/(num_samp_vals-1)
#        # Make sure it doesn't overrun the end of the array
#        while sample_dist*(num_samp_vals-1) > num_inp_vals-1:
#            sample_dist = 0.99999999*sample_dist

    # Resample points with interpolation
    for i_point in range(num_samp_vals-1):
        sample_index = sample_dist*i_point
        p1 = ifloor(sample_index)
        v1 = ordered_vals[p1]
        p2 = iceil(sample_index)
        v2 = ordered_vals[p2]
        sample_val = interpolate(x=sample_index, p1=p1, v1=v1, p2=p2, v2=v2)
        sampled_vals.append(sample_val)
    # Add the last point
    sampled_vals.append(ordered_vals[-1])

    assert len(sampled_vals) == num_samp_vals

    return sampled_vals

def calculate_roc_vals(correct_class, rank_vals):
    """Calculate the ROC Curve values for the ordered outcomes"""

    assert len(correct_class) == len(rank_vals)

    # Sort rank vals into descending order
    sorted_indices = sorted(range(len(rank_vals)), key=lambda k: rank_vals[k], reverse=True)
    # Use this to reorder
    sorted_correct_class = [correct_class[i] for i in sorted_indices]
    sorted_rank_vals =     [rank_vals[i] for i in sorted_indices]

    # Calculate the list of cutoffs
    cutoffs = numpy.zeros(len(rank_vals) + 1)
    cutoffs[0] =  max(sorted_rank_vals) + 1
    cutoffs[-1] = min(sorted_rank_vals) - 1
    cutoffs[1:-1] = (numpy.array(sorted_rank_vals[1:]) + numpy.array(sorted_rank_vals[:-1])) / 2.
    cutoffs = cutoffs.tolist()

    T_BOOL = [1.0 if t     else 0.0 for t in sorted_correct_class]
    N_BOOL = [1.0 if not t else 0.0 for t in sorted_correct_class]
    assert sum(T_BOOL) + sum(N_BOOL) == len(sorted_correct_class)

    # Iterate through the data and calculate TP, TN, FP, FN Rates
    POS_TOTAL = sum(T_BOOL)
    NEG_TOTAL = sum(N_BOOL)
    TOTAL = POS_TOTAL + NEG_TOTAL

    # INITIAL VALUES - STARTING FROM CLASSIFYING EVERYTHING AS NEGATIVE: NO POSITIVES, ALL NEGATIVES
    TPR = [0.0]; FPR = [0.0]; TNR = [1.0]; FNR = [1.0];
    PPV = [0.0]; NPV = [NEG_TOTAL/TOTAL];
    ACC = [NEG_TOTAL/TOTAL]
    TP  = [0.0]; FP  = [0.0]; TN  = [NEG_TOTAL]; FN  = [POS_TOTAL];

    for index in range(len(correct_class)):
        tp = sum(T_BOOL[0:index+1])
        fp = sum(N_BOOL[0:index+1])
        tn = sum(N_BOOL[index+1:] )
        fn = sum(T_BOOL[index+1:] )
        assert tp + fn == POS_TOTAL, 'POS TOTAL != TP + FN'
        assert tn + fp == NEG_TOTAL, 'NEG TOTAL != TN + FP'

#        print 'TP', tp, 'TN', tn, 'FP', fp, 'FN', fn

        TP.append(tp)
        FP.append(fp)
        TN.append(tn)
        FN.append(fn)

        # Fraction of positives classed as hits
        TPR.append(tp/(tp+fn))
        # Fraction of negatives classed as hits
        FPR.append(fp/(fp+tn))

        # Fraction of negatives classed as negs
        TNR.append(tn/(tn+fp))
        # Fraction of positives classed as negs
        FNR.append(fn/(fn+tp))

        # Fraction of hits that are positive
        try:                        PPV.append(tp/(tp+fp))
        except ZeroDivisionError:   PPV.append(1.0)
        # Fraction of negs that are negative
        try:                        NPV.append(tn/(tn+fn))
        except ZeroDivisionError:   NPV.append(1.0)

        # Fraction of classifications that are correct
        ACC.append((tp+tn)/TOTAL)

    SENS = TPR
    SPEC = TNR
    PREC = PPV

    ROC_RESULTS =   {
                    'VALS'  : cutoffs,
                    'TP'    : TP,       'FP'    : FP,
                    'TN'    : TN,       'FN'    : FN,
                    'TPR'   : TPR,      'FPR'   : FPR,
                    'TNR'   : TNR,      'FNR'   : FNR,
                    'PPV'   : PPV,      'NPV'   : NPV,
                    'SENS'  : SENS,     'SPEC'  : SPEC,
                    'PREC'  : PREC,     'ACC'   : ACC
                    }

    return ROC_RESULTS
