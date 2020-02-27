import copy, os
import time

from collections import OrderedDict

import numpy
import pandas as pd
import joblib as jl

from libtbx import easy_mp
from libtbx.utils import Failure
from libtbx.utils import Sorry, Failure
from scitbx.array_family import flex

from bamboo.common.path import splice_ext
from bamboo.stats import modified_z_scores

from giant.xray.scaling import IsotropicBfactorScalingFactory
from giant.xray.data import estimate_wilson_b_factor

from multi_dataset_crystalography.functions import wrapper_run, DatasetAligner


class PanddaDataChecker:

    def __init__(self, structure_factors, low_resolution_completeness, all_data_are_valid_values):

        self.name = "PanddaDataChecker"


        self.structure_factors = structure_factors
        self.low_resolution_completeness = float(low_resolution_completeness)
        self.all_data_are_valid_values = bool(all_data_are_valid_values)

        self.rejections = {}


    def __call__(self, mcd, reference=None):
        """Check that the datasets are analysable (have the right mtz columns, etc)"""

        # Create dicts to hold the new metadata
        dataset_meta_column_lables_dict = {}
        dataset_meta_rejected_dict = {}

        # ==============================>
        # Extract structure factor column names
        # ==============================>
        # TODO: make handle lists
        sf_cols = [[str(x) for x in c.split(',')] for c in self.structure_factors]
        # Print the columns we're looking for
        # ==============================>
        # Validate columns
        # ==============================>
        # Record errors for all datasets without erroring (error at end)
        errors = []
        # Check for these columns in each dataset
        for dtag, dataset in mcd.datasets.items():
            # Check that the input files exist
            if not os.path.exists(dataset.model.filename):
                raise Sorry('Model file does not exist for dataset {}: {}'.format(dataset.tag, dataset.model.filename))
            if not os.path.exists(dataset.data.filename):
                raise Sorry('Data file does not exist for dataset {}: {}'.format(dataset.tag, dataset.data.filename))
            # Extract reflection data
            mtz_obj = dataset.data.mtz_object()
            # The structure factors to be used in this dataset
            dataset_sfs = None
            # ==============================>
            # Find which structure factors are present in this dataset
            # ==============================>
            for sf_pair in sf_cols:
                # Check that the data contains the appropriate column
                if mtz_obj.has_column(sf_pair[0]) and mtz_obj.has_column(sf_pair[1]):
                    dataset_sfs = sf_pair
                    break
            # Error if no columns are identified
            if dataset_sfs is None:
                errors.append(
                    'No matching structure factors were found in the reflection data for dataset {}. \n'.format(
                        dataset.tag) + \
                    'Looking for structure factors: \n\t{}\n'.format('\n\t'.join(map(' and '.join, sf_cols))) + \
                    'Structure factors in this dataset: \n\t{}\n'.format('\n\t'.join(mtz_obj.column_labels())) + \
                    'You may need to change the diffraction_data.structure_factors option.')
                continue
            # Store the column labels in the dataset object
            dataset.meta.column_labels = ','.join(dataset_sfs)
            dataset_meta_column_lables_dict[dataset.tag] = ','.join(dataset_sfs)
            # Report which are being used
            # ==============================>
            # Check the data for the selected columns
            # ==============================>
            for c in dataset_sfs:
                # Extract column from the dataset, and associated miller_set
                col = mtz_obj.get_column(c)
                ms_d = col.mtz_crystal().miller_set(anomalous_flag=False)
                # Get the boolean selection for valid reflections from the column
                values, valid = col.extract_values_and_selection_valid(-1).as_tuple()
                # List of reflections that have missing values
                zero_sel = flex.bool(valid.size(), True).set_selected(valid, False)
                # ==============================>
                # Check that all data are valid values
                # ==============================>
                if self.all_data_are_valid_values is True:
                    if sum(zero_sel) > 0:
                        zero_refl = '\n\t'.join([
                            'Reflection: ({:3d}, {:3d}, {:3d}) - resolution: {:6.3f}A - value: {}'.format(
                                i[0], i[1], i[2], d, v)
                            for (i, d), v in
                            zip(ms_d.select(zero_sel).d_spacings(), values.select(zero_sel))])
                        errors.append('Structure factor column "{}" in dataset {} has invalid reflections. \n'.format(c,
                                                                                                                      dataset.tag) + \
                                      '{} reflection(s) are set to N/A or zero. \n'.format(sum(zero_sel)) + \
                                      'You should populate the structure factors for these reflections with their estimated values -- this is normally performed automatically in refinement with phenix or refmac. \n' + \
                                      'Analysing maps with missing reflections (especially low resolution reflections!) will degrade the quality of the analysis. \n' + \
                                      'However, you can continue by setting checks.all_data_are_valid_values=None. \n' + \
                                      'Missing reflections (-1 indicates reflection set to N/A): \n\t{}'.format(
                                          zero_refl))
                        # continue
                # ==============================>
                # Check that the data is complete up until a certain resolution
                # ==============================>
                if self.low_resolution_completeness is not None:
                    # Find selection for the low resolution reflections
                    sel_low = ms_d.resolution_filter_selection(
                        d_min=self.low_resolution_completeness, d_max=999)
                    # Select these reflections in the dataset miller set
                    ms_d_low = ms_d.select(sel_low)
                    # Extract a complete set of miller indices up to cutoff resolution
                    ms_c_low = ms_d.complete_set(d_min_tolerance=0.0,
                                                 d_min=self.low_resolution_completeness,
                                                 d_max=999)
                    # Calculate the lone set of the two sets
                    lone_set = ms_c_low.lone_set(ms_d_low)
                    # Check that there are the same number of reflections in this set as the other set
                    if lone_set.size() != 0:
                        miss_refl = '\n\t'.join(
                            ['Reflection: ({:3d}, {:3d}, {:3d}) - resolution: {:6.3f}A'.format(i[0], i[1], i[2], d)
                             for (i, d) in lone_set.d_spacings()])
                        dataset.data.miller_arrays['scaled'] = None

                        self.rejections[dataset.tag] = "rejected - rmsd to reference"
                        dataset_meta_rejected_dict[dataset.tag] = "rejected - rmsd to reference"

                        continue

                    else:
                        dataset.meta.rejected = False
                        dataset_meta_rejected_dict[dataset.tag] = False

                    # Calculate overlap between low resolution set and valid set to ensure none missing from low resolution set
                    zero_sel_low = zero_sel.select(sel_low)
                    # Check if any low resolution reflections are invalid
                    if sum(zero_sel_low) > 0:
                        zero_refl = '\n\t'.join([
                            'Reflection: ({:3d}, {:3d}, {:3d}) - resolution: {:6.3f}A - value: {}'.format(
                                i[0], i[1], i[2], d, v)
                            for (i, d), v in zip(ms_d_low.select(zero_sel_low).d_spacings(),
                                                 values.select(sel_low).select(zero_sel_low))])
                        errors.append(
                            'Structure factor column "{}" in dataset {} has invalid reflections below {}A. \n'.format(c,
                                                                                                                      dataset.tag,
                                                                                                                      self.low_resolution_completeness) + \
                            '{} reflection(s) are set to N/A or zero. \n'.format(sum(zero_sel_low)) + \
                            'You should populate the structure factors for these reflections with their estimated values -- this is normally performed automatically in refinement with phenix or refmac. \n' + \
                            'Analysing maps with missing reflections (especially low resolution reflections!) will degrade the quality of the analysis. \n' + \
                            'You can ignore this -- though I really would not -- by setting checks.low_resolution_completeness to None. \n' + \
                            'Missing reflections (-1 indicates reflection set to N/A): \n\t{}'.format(zero_refl))
                        continue

        # # Collect rejected datasets
        # for dataset in datasets:
        #     if dataset.tag in self.rejections:
        #         pass
        #     else:
        #         self.rejections[dataset.tag] = False

        # checked_datasets = [dataset for dataset in datasets if self.rejections[dataset.tag] is False]

        # ==============================>
        # Create a new dataset from non-rejected datasets
        # ==============================>
        checked_datasets = {dtag: d
                            for dtag, d
                            in mcd.datasets.items()
                            if not (dtag in self.rejections)}
        new_dataset = mcd.new_from_datasets(datasets=checked_datasets)


        return new_dataset

    def log(self):
        log = OrderedDict()
        log["Rejected datasets"] = "Datasets rejected: \n{}\n".format("\n".join(["{}: {}".format(dtag, reason)
                                                                                 for dtag, reason
                                                                                 in self.rejections.items()]))
        return log

    def repr(self):
        repr = {}
        repr["structure_factors"] = self.structure_factors
        repr["low_resolution_completeness"] = self.low_resolution_completeness
        repr["all_data_are_valid_values"] = self.all_data_are_valid_values
        return repr


class PanddaDiffractionScaler:

    def __init__(self, apply_b_factor_scaling):

        self.name = "PanddaDiffractionScaler"

        self.rejections = {}

        self.apply_b_factor_scaling = bool(apply_b_factor_scaling)

    def __call__(self, mcd, reference_dataset=None):
        """Extract amplitudes and phases for creating map"""

        # ==============================>
        # Load input
        # ==============================>

        # ==============================>
        # Prepare objects for scaling
        # ==============================>
        # Extract reference miller array and prepare for scaling
        ref_dataset = reference_dataset
        ref_miller = ref_dataset.data.miller_arrays[ref_dataset.meta.column_labels]
        factory = IsotropicBfactorScalingFactory(reference_miller_array=ref_miller.as_intensity_array())

        scaling_dict = {}

        # ==============================>
        # Report
        # ==============================>
        for dtag, dataset in mcd.datasets.items():

            # Get the columns to be loaded for this dataset
            dataset_sfs = dataset.meta.column_labels
            # ==============================>
            # Load structure factors
            # ==============================>
            if dataset_sfs in dataset.data.miller_arrays.keys():
                ma_unscaled_com = dataset.data.miller_arrays[dataset_sfs]
            else:
                ma_unscaled_com = dataset.data.get_structure_factors(columns=dataset_sfs)
                # Store in the dataset object
                dataset.data.miller_arrays[dataset_sfs] = ma_unscaled_com

            # ==============================>
            # Reject sets with zero amplitude unscaled
            # TODO: is this reasonable? Breaks otherwise...
            # ==============================>
            if (ma_unscaled_com.as_amplitude_array().data().as_numpy_array() == 0).any():
                self.rejections[dataset.tag] = "rejected - zero amplitude reflection"
                continue
            # ==============================>
            # Extract intensities and phases
            # ==============================>
            assert ma_unscaled_com.is_complex_array()
            ma_unscaled_int = ma_unscaled_com.as_intensity_array()
            ma_unscaled_phs = ma_unscaled_com * (1.0 / ma_unscaled_com.as_amplitude_array().data())
            # ==============================>
            # Scale to the reference dataset
            # ==============================>
            if not (hasattr(dataset.data, 'scaling') and (dataset.data.scaling is not None)):
                try:
                    scaling = factory.calculate_scaling(miller_array=ma_unscaled_int)
                    dataset.data.scaling = scaling
                    scaling_dict[dataset.tag] = scaling
                except Exception as e:
                    self.rejections[dataset.tag] = "rejected - scaling failed"

                    continue

            # self.rejections[dataset.tag] = False

            # ==============================>
            # Select which data to use for analysis
            # ==============================>
            if not self.apply_b_factor_scaling:
                ma_scaled_com = ma_unscaled_com
            else:
                # Need to create a copy to preserve the x-values of the original scaling object
                new_scaling = copy.deepcopy(scaling)
                new_scaling.new_x_values(x_values=ma_unscaled_int.d_star_sq().data())
                ma_scaled_int = ma_unscaled_int.array(
                    data=new_scaling.transform(ma_unscaled_int.data())).set_observation_type_xray_intensity()
                ma_scaled_com = ma_unscaled_phs * ma_scaled_int.as_amplitude_array().data()
                # Check for nan values and set to zero
                ma_scaled_com.data().set_selected((ma_unscaled_int.data() == 0.0), 0 + 0j)

            assert ma_scaled_com.is_complex_array()

            # Add scaled array to dataset
            dataset.data.miller_arrays['scaled'] = ma_scaled_com

        # ==============================>
        # Create a new dataset from scaled datasets
        # ==============================>

        checked_datasets = {dtag: d
                            for dtag, d
                            in mcd.datasets.items()
                            if not (dtag in self.rejections)}
        new_dataset = mcd.new_from_datasets(datasets=checked_datasets)


        return new_dataset

    def log(self):
        log = OrderedDict()
        log["Rejected datasets"] = "Datasets rejected: \n{}\n".format("\n".join(["{}: {}".format(dtag, reason)
                                                                                 for dtag, reason
                                                                                 in self.rejections.items()]))
        return log

    def repr(self):
        repr = {}
        repr["rejections"] = self.rejections
        repr["apply_b_factor_scaling"] = self.apply_b_factor_scaling
        return repr


class PanddaDatasetFiltererWilsonRMSD:

    def __init__(self, max_wilson_plot_z_score, apply_b_factor_scaling, dataset_sfs):

        self.name = "PanddaDatasetFiltererWilsonRMSD"

        self.rejections_structural_deviation = {}
        self.rejections_wilson_RMSD = {}

        self.max_wilson_plot_z_score = float(max_wilson_plot_z_score)
        self.apply_b_factor_scaling = bool(apply_b_factor_scaling)
        self.dataset_sfs = dataset_sfs[0]

    def __call__(self, mcd, reference=None):

        records = []
        for dtag, dataset in mcd.datasets.items():
            records.append(self.make_wilson_record(dataset))

        wilson_table = pd.DataFrame(records).set_index("dtag")

        # ==============================>
        # Exclude from characterisation if poor quality diffraction
        # ==============================>
        for dtag, dataset in mcd.datasets.items():
            if (wilson_table.get_value(index=dataset.tag,
                                       col='scaled_wilson_rmsd_all_z') > self.max_wilson_plot_z_score) or \
                    (wilson_table.get_value(index=dataset.tag,
                                            col='scaled_wilson_rmsd_<4A_z') > self.max_wilson_plot_z_score) or \
                    (wilson_table.get_value(index=dataset.tag,
                                            col='scaled_wilson_rmsd_>4A_z') > self.max_wilson_plot_z_score) or \
                    (wilson_table.get_value(index=dataset.tag,
                                            col='scaled_wilson_ln_rmsd_z') > self.max_wilson_plot_z_score) or \
                    (wilson_table.get_value(index=dataset.tag,
                                            col='scaled_wilson_ln_dev_z') > self.max_wilson_plot_z_score):

                self.rejections_wilson_RMSD[dataset.tag] = 'exclude_from_characterisation'

        # ==============================>
        # Create new dataset
        # ==============================>
        new_datasets = {dtag: d
                        for dtag, d
                        in mcd.datasets.items()
                        if not (dtag in self.rejections_wilson_RMSD)}
        new_dataset = mcd.new_from_datasets(datasets=new_datasets)

        return new_dataset

    def make_wilson_record(self, dataset):
        # TODO: Move these to report or summary method

        scaling = dataset.data.scaling
        ma_unscaled_com = dataset.data.miller_arrays[self.dataset_sfs]
        ma_scaled_com = dataset.data.miller_arrays["scaled"]

        record = {"dtag": dataset.tag,
                  "scaling": scaling}

        # ==============================>
        # Record metrics for unscaled data
        # ==============================>
        # Select high resolution and low resolution separately
        high_res_sel = scaling.x_values > 1 / (4.0 ** 2)
        low_res_sel = scaling.x_values <= 1 / (4.0 ** 2)
        # Log unscaled rmsds values
        record["unscaled_wilson_rmsd_all"] = numpy.round(scaling.unscaled_rmsd, 3)

        record['unscaled_wilson_rmsd_<4A'] = numpy.round(
            scaling.rmsd_to_ref(values=scaling.scl_values, sel=high_res_sel), 3)
        record['unscaled_wilson_rmsd_>4A'] = numpy.round(
            scaling.rmsd_to_ref(values=scaling.scl_values, sel=low_res_sel), 3)
        # Log the scaled log-rmsd values
        record['unscaled_wilson_ln_rmsd'] = numpy.round(scaling.unscaled_ln_rmsd, 3)
        record['unscaled_wilson_ln_dev'] = numpy.round(scaling.unscaled_ln_dev, 3)
        # ==============================>
        # Report metrics for scaled data
        # ==============================>
        # Log the scaling
        record['applied_b_factor_scaling'] = numpy.round(scaling.scaling_b_factor, 3)
        # Log the scaled rmsd values
        record['scaled_wilson_rmsd_all'] = numpy.round(scaling.scaled_rmsd, 3)
        record['scaled_wilson_rmsd_<4A'] = numpy.round(scaling.rmsd_to_ref(values=scaling.out_values, sel=high_res_sel),
                                                       3)
        record['scaled_wilson_rmsd_>4A'] = numpy.round(scaling.rmsd_to_ref(values=scaling.out_values, sel=low_res_sel),
                                                       3)
        # Log the scaled log-rmsd values
        record['scaled_wilson_ln_rmsd'] = numpy.round(scaling.scaled_ln_rmsd, 3)
        record['scaled_wilson_ln_dev'] = numpy.round(scaling.scaled_ln_dev, 3)

        # ==============================>
        # Wilson B-factors
        # ==============================>
        dataset.meta.unscaled_wilson_b = estimate_wilson_b_factor(miller_array=ma_unscaled_com)
        dataset.meta.scaled_wilson_b = estimate_wilson_b_factor(miller_array=ma_scaled_com)

        record['unscaled_wilson_B'] = dataset.meta.unscaled_wilson_b
        record['scaled_wilson_B'] = dataset.meta.scaled_wilson_b

        # # ==============================>
        # # Update Z-score columns
        # # ==============================>
        record['unscaled_wilson_rmsd_all_z'] = modified_z_scores(record['unscaled_wilson_rmsd_all'])
        record['unscaled_wilson_rmsd_<4A_z'] = modified_z_scores(record['unscaled_wilson_rmsd_<4A'])
        record['unscaled_wilson_rmsd_>4A_z'] = modified_z_scores(record['unscaled_wilson_rmsd_>4A'])
        record['unscaled_wilson_ln_rmsd_z'] = modified_z_scores(record['unscaled_wilson_ln_rmsd'])
        record['unscaled_wilson_ln_dev_z'] = modified_z_scores(record['unscaled_wilson_ln_dev'])
        record['scaled_wilson_rmsd_all_z'] = modified_z_scores(record['scaled_wilson_rmsd_all'])
        record['scaled_wilson_rmsd_<4A_z'] = modified_z_scores(record['scaled_wilson_rmsd_<4A'])
        record['scaled_wilson_rmsd_>4A_z'] = modified_z_scores(record['scaled_wilson_rmsd_>4A'])
        record['scaled_wilson_ln_rmsd_z'] = modified_z_scores(record['scaled_wilson_ln_rmsd'])
        record['scaled_wilson_ln_dev_z'] = modified_z_scores(record['scaled_wilson_ln_dev'])

        return record

    def log(self):
        log = OrderedDict()
        log["Rejected datasets"] = "Datasets rejected: \n{}\n".format("\n".join(["{}: {}".format(dtag, reason)
                                                                                 for dtag, reason
                                                                                 in self.rejections_wilson_RMSD.items()]))
        return log

    def repr(self):
        repr = {}
        repr["rejections_structural_deviation"] = self.rejections_structural_deviation
        repr["rejections_wilson_RMSD"] = self.rejections_wilson_RMSD
        repr["max_wilson_plot_z_score"] = self.max_wilson_plot_z_score
        repr["apply_b_factor_scaling"] = self.apply_b_factor_scaling
        repr["dataset_sfs"] = self.dataset_sfs
        return repr


class PanddaDatasetFilterer:

    def __init__(self, similar_models_only, max_rfree, same_space_group_only):

        self.name = "PanddaDatasetFilterer"

        self.similar_models_only = similar_models_only
        self.max_rfree = max_rfree
        self.same_space_group_only = same_space_group_only

        self.rejections_structural_deviation = {}

    def __call__(self, mcd, reference_dataset=None):
        """Filter out the datasets which contain different protein models (i.e. protein length, sequence, etc)"""

        # ==============================>
        # Remove poor/incompatible datasets
        # ==============================>
        for dtag, dataset in mcd.datasets.items():
            if self.same_space_group_only and (
                    dataset.model.space_group.info().symbol_and_number() != reference_dataset.model.space_group.info().symbol_and_number()):
                self.rejections_structural_deviation[dataset.tag] = "rejected - different space group"

            elif dataset.model.input.get_r_rfree_sigma().r_free > self.max_rfree:
                self.rejections_structural_deviation[dataset.tag] = "rejected - rfree"

            elif self.similar_models_only and (
                    not dataset.model.hierarchy.is_similar_hierarchy(reference_dataset.model.hierarchy)):
                self.rejections_structural_deviation[dataset.tag] = "rejected - non-identical structures"


        # # ==============================>
        # # Create a new dataset
        # # ==============================>
        new_datasets = {dtag: d
                        for dtag, d
                        in mcd.datasets.items()
                        if not (dtag in self.rejections_structural_deviation)}
        new_dataset = mcd.new_from_datasets(datasets=new_datasets)

        return new_dataset

    def log(self):
        log = OrderedDict()
        log["Rejected datasets"] = "Datasets rejected: \n{}\n".format("\n".join(["{}: {}".format(dtag, reason)
                                                                                 for dtag, reason
                                                                                 in self.rejections_structural_deviation.items()]))
        return log

    def repr(self):
        repr = {}
        repr["similar_models_only"] = self.similar_models_only
        repr["max_rfree"] = self.max_rfree
        repr["same_space_group_only"] = self.same_space_group_only
        repr["rejections_structural_deviation"] = self.rejections_structural_deviation

        return repr


class PanddaDefaultStructureAligner:

    def __init__(self, method, cpus):
        self.name = "PanddaDefaultStructureAligner"
        self.method = method
        self.cpus = cpus
        self.alignments = OrderedDict()

    def __call__(self, mcd, reference_dataset=None):
        """Align each structure the reference structure"""

        assert self.method in ['local', 'global'], 'METHOD NOT DEFINED: {!s}'.format(self.method)

        # ==============================>
        # Select the datasets for alignment
        # ==============================>
        datasets_for_alignment = [d for dtag, d in mcd.datasets.items()]
        # ==============================>
        # Delete alignments
        # ==============================>
        for d in datasets_for_alignment:
            d.model.alignment = None

        # ==============================>
        # Generate the alignments for each structure
        # ==============================>
        arg_list = [DatasetAligner(model=d.model, other=reference_dataset.model, method=self.method, id=d.tag) for d
                    in datasets_for_alignment]
        dataset_alignments = jl.Parallel(n_jobs=self.cpus,
                                         verbose=15)(jl.delayed(wrapper_run)(arg)
                                                     for arg
                                                     in arg_list)
        # dataset_alignments = easy_mp.pool_map(func=wrapper_run, args=arg_list, processes=self.cpus)

        # ==============================>
        # Catch errors and print at end
        # ==============================>
        errors = []
        for dataset, alignment in zip(datasets_for_alignment, dataset_alignments):
            # If errored, print and record
            if isinstance(alignment, str):
                errors.append((dataset, alignment))
                self.alignments[dataset.tag] = False
                print(alignment)
                continue
            # Attach alignment to dataset
            assert dataset.tag == alignment.id
            self.alignments[alignment.id] = True
            dataset.model.alignment = alignment
            # Output an aligned copy of the structure
            aligned_struc = dataset.model.hierarchy.deep_copy()
            aligned_struc.atoms().set_xyz(
                dataset.model.alignment.nat2ref(coordinates=dataset.model.hierarchy.atoms().extract_xyz()))
            # TODO: restore this functionality?
            # aligned_struc.write_pdb_file(file_name=os.path.join(self.file_manager.get_dir('aligned_structures'),
            #                                                     '{!s}-aligned.pdb'.format(dataset.tag)))
            # aligned_struc.write_pdb_file(
            #     file_name=splice_ext(dataset.file_manager.get_file('aligned_model'), 'ref', position=-1))
            # Write alignment summary to log

        # ==============================>
        # Report Errors
        # ==============================>
        if errors:
            for dataset, err_msg in errors:
                print('Failed to align dataset {}'.format(dataset.tag))
                print(err_msg)
            raise Failure('Failed to align {} datasets. Error messages printed above.'.format(len(errors)))

        # ==============================>
        # Make new dataset
        # ==============================>
        new_datasets = {d.tag: d
                        for d
                        in datasets_for_alignment}
        new_dataset = mcd.new_from_datasets(datasets=new_datasets)

        return new_dataset

    def log(self):
        log = OrderedDict()
        log["aligned_datasets"] = "Datasets aligned: \n{}\n".format("\n".join(["{}".format(dtag)
                                                                               for dtag, val
                                                                               in self.alignments.items()
                                                                               if val is True]))
        log["not_aligned_datasets"] = "Datasets not aligned: \n{}\n".format("\n".join(["{}".format(dtag)
                                                                                       for dtag, val
                                                                                       in self.alignments.items()
                                                                                       if val is False]))
        return log

    def repr(self):
        repr = {}
        repr["method"] = self.method
        repr["cpus"] = self.cpus
        repr["alignments"] = self.alignments

        return repr


def same_spacegroup(dataset, reference_dataset):
    if dataset.model.space_group.info().symbol_and_number() != reference_dataset.model.space_group.info().symbol_and_number():
        return False
    else:
        return True


def filter_spacegroup(mcd, reference_dataset):
    accepted_items = list(filter(lambda item: same_spacegroup(item[1],
                                                              reference_dataset),
                                 mcd.datasets.items()))
    new_mcd = mcd.new_from_datasets(datasets=OrderedDict(accepted_items))
    return new_mcd


def reasonable_rfree(dataset, max_rfree):
    if dataset.model.input.get_r_rfree_sigma().r_free < max_rfree:
        return True
    else:
        return False


def filter_rfree(mcd, max_rfree):
    accepted_items = list(filter(lambda item: reasonable_rfree(item[1],
                                                              max_rfree),
                                 mcd.datasets.items()))
    new_mcd = mcd.new_from_datasets(datasets=OrderedDict(accepted_items))
    return new_mcd


def similar_models(dataset, reference_dataset):
    if dataset.model.hierarchy.is_similar_hierarchy(reference_dataset.model.hierarchy):
        return True
    else:
        return False


def filter_similar_models(mcd, reference_dataset):
    accepted_items = list(filter(lambda item: similar_models(item[1],
                                                              reference_dataset),
                                 mcd.datasets.items()))
    new_mcd = mcd.new_from_datasets(datasets=OrderedDict(accepted_items))
    return new_mcd


def filter_dataset(filter_func, mcd):
    accepted_items = list(filter(lambda item: filter_func(item[1]),
                                 mcd.datasets.items()))
    new_mcd = mcd.new_from_datasets(datasets=OrderedDict(accepted_items))
    return new_mcd