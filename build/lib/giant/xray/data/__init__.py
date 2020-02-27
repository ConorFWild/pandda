import cctbx.miller

import numpy

from scitbx.array_family import flex

from bamboo.common import Info
from giant.stats.optimisation import LinearScaling

def extract_structure_factors(mtz_object, ampl_label, phas_label):
    # Get the crystal symmetry from the amplitudes' crystal
    try:
        ampl_col = mtz_object.get_column(ampl_label)
    except:
        raise Exception('Amplitude column not found: {}. Have you specified the right columns?'.format(ampl_label))
    # Get the symmetry associated with the column
    crystal_symmetry = ampl_col.mtz_crystal().crystal_symmetry()
    # Extract amplitudes and phases
    try:
        sf_com = mtz_object.extract_complex(column_label_ampl=ampl_label, column_label_phi=phas_label)
    except:
        raise Exception('Could not extract structure factors - Amplitudes:{}, Phases:{}. Have you specified the right columns?'.format(ampl_label, phas_label))
    # Convert to miller array
    mill_set = cctbx.miller.set(crystal_symmetry=crystal_symmetry, indices=sf_com.indices)
    mill_sfs = mill_set.array(sf_com.data)
#    mill_sfs.set_observation_type_xray_amplitude()
    mill_sfs = mill_sfs.as_non_anomalous_array()
    assert mill_sfs.is_complex_array(), 'STRUCTURE FACTORS SHOULD BE COMPLEX?!'
    return mill_sfs

def estimate_wilson_b_factor(miller_array, low_res_cutoff=4.0):

    miller_array = miller_array.resolution_filter(d_max=low_res_cutoff).as_intensity_array()
    # print(miller_array.data().as_numpy_array())
    # print(miller_array.data())

    # Setup binner and extract radial averages
    binner = miller_array.setup_binner(auto_binning=True)
    binned = miller_array.wilson_plot(use_binning=True)
    #
    # print(binned.data[1:-1])
    # print(type(binned.data[1:-1]))
    # print(type(binned.data[1]))
    # print(binned.data[1:-1])
    # print([type(x) for x in binned.data[1:-1]])
    # print(type(1.0))
    binned_data = [float(x) if type(x) == float else 1.0 for x in binned.data[1:-1]]
    # print(flex.double(binned_data))

    # Convert to scale
    y_values = flex.log(flex.double(binned_data))
    x_values = flex.pow2(binner.bin_centers(1))
    # Check all values are valid
    # mask = flex.bool((True - numpy.isnan(list(y_values)) - numpy.isnan(list(x_values))).tolist())
    mask = flex.bool((True ^ numpy.isnan(list(y_values)) ^ numpy.isnan(list(x_values))).tolist())

    # Perform scaling
    scl = LinearScaling(x_values   = x_values.select(mask),
                        ref_values = y_values.select(mask))

    return -0.5*scl.optimised_values[1]

#def extract_structure_factors(mtz_object, ampl_label, phas_label):
#
#    # Extract matching miller arrays
#    match_arrs = [a for a in mtz_object.as_miller_arrays() if a.info().labels==[ampl_label, phas_label]]
#    if not match_arrs: raise Exception('Could not extract structure factors - Amplitudes:{}, Phases:{}. Have you specified the right columns?'.format(ampl_label, phas_label))
#    assert len(match_arrs) == 1
#    mill_arr = match_arrs[0]
#    assert mill_arr.is_complex_array(), 'STRUCTURE FACTORS SHOULD BE COMPLEX?!'
#    mill_arr = mill_arr.as_non_anomalous_array()
#    return mill_arr

def sort_column_labels(mtz_object, return_labels=True):
    """Sort columns by type (F, I, etc). If return_labels=True, return labels, else return column data"""

    col_labels = mtz_object.column_labels()
    col_types = mtz_object.column_types()

    lab_hash = {}

    # H: Miller Indices, F: Amplitudes, I: Integers, P: Phases, Q: Standard Deviations, J: Intensities
    mill = [col_labels[i] for i,t in enumerate(col_types) if t=='H']
    sfac = [col_labels[i] for i,t in enumerate(col_types) if t=='F']
    inty = [col_labels[i] for i,t in enumerate(col_types) if t=='J']
    sdev = [col_labels[i] for i,t in enumerate(col_types) if t=='Q']
    phas = [col_labels[i] for i,t in enumerate(col_types) if t=='P']
    ints = [col_labels[i] for i,t in enumerate(col_types) if t=='I']

    lab_hash['miller'] = mill

    # Find the main F col
    lab_hash['f'] = [s for s in sfac if (s in ['F','FP','FCTR','FOSC','FOBS'] \
                                     or s.startswith('F_') or s.startswith('FP_') \
                                     or s.upper().startswith('F-OBS') or s.upper().startswith('FOUT_'))]
    # Use these if no others
    if not lab_hash['f']: lab_hash['f'] = [s for s in sfac if s.startswith('FP')]
    # Find the main SIGF col
    lab_hash['sigf'] = [s for s in sdev if (s in ['SIGF','SIGFP','SIGFOBS'] \
                            or s.startswith('SIGF_') or s.startswith('SIGFP_') \
                            or s.upper().startswith('SIGF-OBS') or s.upper().startswith('SIGFOUT_'))]
    # Use these if no others
    if not lab_hash['sigf']: lab_hash['sigf'] = [s for s in sdev if s.startswith('SIGFP')]

    # Find the I cols
    lab_hash['i'] = [s for s in inty if s.startswith('I')]
    # Find the SIGI cols
    lab_hash['sigi'] = [s for s in sdev if s.startswith('SIGI')]

    # Find the F_calc
    lab_hash['calc_f'] = [s for s in sfac if (s=='FC' or s=='FMODEL' or s.startswith('FC_') or s.upper().startswith('F-MODEL'))]
    # Find the PHI_calcs
    lab_hash['calc_p'] = [p for p in phas if (p=='PHIC' or p=='PHIFMODEL' or p.startswith('PHIC_') or p.upper().startswith('PHIF-MODEL'))]

    # Find the main phase col
    lab_hash['phi'] = lab_hash['p_calc']

    # Find the RFree Flag
    lab_hash['r_flags'] = [r for r in ints if ('RFREE' in r.upper() or 'FREER' in r.upper() or r=='FREE' or r.upper().startswith('R-FREE'))]

    # 2FOFC cols
    wt_f_map_opts = ['2FOFCWT','FWT']
    wt_p_map_opts = ['PH2FOFCWT','PHWT','PHFWT']
    # FOFC cols
    wt_f_dif_opts = ['FOFCWT','DELFWT']
    wt_p_dif_opts = ['PHFOFCWT','DELPHWT','PHDELWT']
    # Record 2FOFC pair (composite map)
    lab_hash['2fofc_f'] = [s for s in sfac if s in wt_f_map_opts]
    lab_hash['2fofc_p'] = [p for p in phas if p in wt_p_map_opts]
    # Record FOFC pair (different map)
    lab_hash['fofc_f'] = [s for s in sfac if s in wt_f_dif_opts]
    lab_hash['fofc_p'] = [p for p in phas if p in wt_p_dif_opts]

    # XXX Will probably delete this later
    for k in lab_hash:
        for v in lab_hash[k]:
            col_labels.remove(v)
    lab_hash['unknown'] = col_labels

    return Info(lab_hash)
