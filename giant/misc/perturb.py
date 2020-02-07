
import os, sys, random

from iotbx import reflection_file_reader
from iotbx.reflection_file_utils import reflection_file_server

from cctbx.array_family import flex

def perturb_structure_factors_proportional_to_observation_error(mtzin, mtzout, seed=0):
    """Take the amplitudes from perturbation_array and use them to perturb the amplitudes in main_array (by gaussian noise)"""

    assert os.path.exists(mtzin), 'MTZ File does not exist! {!s}'.format(mtzin)
    assert not os.path.exists(mtzout), 'MTZ File already exists! {!s}'.format(mtzout)

    random.seed(seed)

    output_column_label = 'F_PERTURBED'

    # Create an mtz file object
    mtz_file = reflection_file_reader.any_reflection_file(mtzin)
    # Create file server to serve up the mtz file
    mtz_server = reflection_file_server(crystal_symmetry=None, force_symmetry=True, reflection_files=[mtz_file], err=sys.stderr)

    # At the moment hardcode a ccp4 label
    sfobs_label = 'F,SIGF'
    rfree_label = 'FreeR_flag'

    # Get the observed miller array, and the R-free flags
    sfobs_miller = mtz_server.get_miller_array(sfobs_label)
    rfree_miller = mtz_server.get_miller_array(rfree_label)

    # F column
    sfobs_amplts = sfobs_miller.data()
    # SIGF column
    sfobs_scales = sfobs_miller.sigmas()

    # Add to observed data
    perturbed_data = [f + e*random.gauss(0,1) for f,e in zip(sfobs_amplts, sfobs_scales)]
    # Zero any negative amplitudes
    perturbed_data = [f if f>0.0 else 0.0 for f in perturbed_data]

    # Generate new miller array
    perturbed_array = sfobs_miller.customized_copy(data=flex.double(perturbed_data))
    assert perturbed_array.is_real_array()

    print 'CORRELATION:', sfobs_miller.correlation(perturbed_array).coefficient()

    # Convert to mtz dataset
    perturbed_mtz_dataset = perturbed_array.as_mtz_dataset(output_column_label)
    # Copy across the R-Free
    perturbed_mtz_dataset.add_miller_array(rfree_miller, rfree_label)

    # Write out mtz file
    perturbed_mtz_object = perturbed_mtz_dataset.mtz_object()
    perturbed_mtz_object.write(file_name=mtzout)

    return output_column_label

def perturb_structure_factors_proportional_to_model_error(mtzin, mtzout, seed=0):
    """Perturb the observed structure factors proportional to the model error (FOBS-FCALC)"""

    random.seed(seed)

    assert os.path.exists(mtzin), 'MTZ File does not exist! {!s}'.format(mtzin)
    assert not os.path.exists(mtzout), 'MTZ File already exists! {!s}'.format(mtzout)

    output_column_label = 'F_PERTURBED'

    # Create an mtz file object
    mtz_file = reflection_file_reader.any_reflection_file(mtzin)
    # Create file server to serve up the mtz file
    mtz_server = reflection_file_server(crystal_symmetry=None, force_symmetry=True, reflection_files=[mtz_file], err=sys.stderr)

    # At the moment hardcode a ccp4 label
    sfobs_label = 'F,SIGF'
    sfdif_label = 'DELFWT,PHDELWT'
    rfree_label = 'FreeR_flag'

    # Get the observed miller array, and the R-free flags
    sfobs_miller = mtz_server.get_miller_array(sfobs_label)
    sfdif_miller = mtz_server.get_miller_array(sfdif_label).amplitudes()
    rfree_miller = mtz_server.get_miller_array(rfree_label)

    # F column
    sfobs_amplts = sfobs_miller.data()
    # FOBS-FCALC column
    sfobs_scales = sfdif_miller.data()

    # Add to observed data
    perturbed_data = [f + e*random.gauss(0,1) for f,e in zip(sfobs_amplts, sfobs_scales)]
    # Zero any negative amplitudes
    perturbed_data = [f if f>0.0 else 0.0 for f in perturbed_data]

    # Generate new miller array
    perturbed_array = sfobs_miller.customized_copy(data=flex.double(perturbed_data))
    assert perturbed_array.is_real_array()

    print 'CORRELATION:', sfobs_miller.correlation(perturbed_array).coefficient()

    # Convert to mtz dataset
    perturbed_mtz_dataset = perturbed_array.as_mtz_dataset(output_column_label)
    # Copy across the R-Free
    perturbed_mtz_dataset.add_miller_array(rfree_miller, rfree_label)

    # Write out mtz file
    perturbed_mtz_object = perturbed_mtz_dataset.mtz_object()
    perturbed_mtz_object.write(file_name=mtzout)

    return output_column_label

def calculate_rms_mtz(mtzfiles, mtzout):
    """Take a list of mtz files and calculate the rms of the structure factors"""

    pass



