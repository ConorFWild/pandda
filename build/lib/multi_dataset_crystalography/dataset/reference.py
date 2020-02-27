import os, copy

from libtbx.utils import Sorry, Failure
import scitbx

from bamboo.common.path import easy_directory, rel_symlink

from giant.dataset import ModelAndData
from giant.structure.align import GlobalAlignment
from giant.structure.select import calphas, protein, sel_altloc, non_water


class PanddaReferenceDataset(ModelAndData):
    _origin_shift = None

    def __init__(self, model, data):
        super(PanddaReferenceDataset, self).__init__(model=model, data=data)

        self.set_origin_shift((0., 0., 0.))
        self.child = None

    def set_origin_shift(self, shift):
        """Creates an alignment corresponding to an origin shift"""

        self._origin_shift = shift
        r = scitbx.matrix.rec([1, 0, 0, 0, 1, 0, 0, 0, 1], (3, 3))
        t = scitbx.matrix.rec(shift, (3, 1))
        rt = scitbx.matrix.rt((r, t))
        self.model.alignment = GlobalAlignment(alignment_mx=rt,
                                               alignment_sites=calphas(self.model.hierarchy).atoms().extract_xyz(),
                                               id='ref')
        return self.model.alignment

    def origin_shift(self):
        return self._origin_shift

    def ref2grid(self, *args, **kwargs):
        return self.model.alignment.nat2ref(*args, **kwargs)

    def grid2ref(self, *args, **kwargs):
        return self.model.alignment.ref2nat(*args, **kwargs)

    def copy(self):
        return copy.deepcopy(self)


class DefaultReferenceGetter:

    def __init__(self, out_dir, reference_pdb_path, reference_mtz_path, reference_structure_factors, structure_factors):
        self.out_dir = out_dir
        self.reference_pdb_path = reference_pdb_path
        self.reference_mtz_path = reference_mtz_path
        self.reference_structure_factors = reference_structure_factors
        self.structure_factors = structure_factors

    def __call__(self, datasets):

        # Use given reference dataset, or select reference dataset
        if self.reference_pdb_path and self.reference_mtz_path:
            ref_pdb, ref_mtz = self.reference_pdb_path, self.reference_mtz_path
        else:
            ref_pdb, ref_mtz = self.select_reference_dataset(datasets,
                                                             method='resolution'
                                                             )

        # Load the reference dataset
        reference_dataset = self.load_reference_dataset(ref_pdb_path=ref_pdb,
                                                        ref_mtz_path=ref_mtz
                                                        )

        # return output_paths
        return reference_dataset

    def select_reference_dataset(self, datasets, method='resolution', max_rfree=0.4, min_resolution=5):
        """Select dataset to act as the reference - scaling, aligning etc"""

        assert method in ['resolution',
                          'rfree'], 'METHOD FOR SELECTING THE REFERENCE DATASET NOT RECOGNISED: {!s}'.format(method)

        # ==============================>
        # Get the potential reference datasets
        # ==============================>
        # filtered_datasets = self.datasets.mask(mask_name='valid - all')
        filtered_datasets = [d for dtag, d in datasets.items()]
        if not filtered_datasets: raise Failure(
            "Can't select a reference dataset - NO SUITABLE (NON-REJECTED) DATASETS REMAINING")
        # ==============================>
        # Select by either R-free or Resolution
        # ==============================>
        if method == 'rfree':
            # Get RFrees of datasets (set to dummy value of 999 if resolution is too high so that it is not selected)
            r_frees = [d.model.input.get_r_rfree_sigma().r_free if (
                    d.data.mtz_object().max_min_resolution()[1] < min_resolution) else 999
                       for d
                       in filtered_datasets]
            if len(r_frees) == 0: raise Exception(
                'NO DATASETS BELOW RESOLUTION CUTOFF {!s}A - CANNOT SELECT REFERENCE DATASET'.format(
                    min_resolution))
            ref_dataset_index = r_frees.index(min(r_frees))

        elif method == 'resolution':
            # Get Resolutions of datasets (set to dummy value of 999 if r-free is too high so that it is not selected)
            resolns = [d.data.mtz_object().max_min_resolution()[1] if (
                    d.model.input.get_r_rfree_sigma().r_free < max_rfree) else 999
                       for d
                       in filtered_datasets]
            if len(resolns) == 0: raise Exception(
                'NO DATASETS BELOW RFREE CUTOFF {!s} - CANNOT SELECT REFERENCE DATASET'.format(max_rfree))
            ref_dataset_index = resolns.index(min(resolns))
        # ==============================>
        # Report and return
        # ==============================>
        reference = filtered_datasets[ref_dataset_index]

        return reference.model.filename, reference.data.filename

    def load_reference_dataset(self, ref_pdb_path, ref_mtz_path):
        """Set the reference dataset, to which all other datasets will be aligned and scaled"""

        # ==============================>
        # Output links to reference files
        # ==============================>
        link_ref_pdb = ref_pdb_path
        link_ref_mtz = ref_mtz_path
        # ==============================>
        # Remove any old links to dataset
        # ==============================>
        if os.path.abspath(ref_pdb_path) != os.path.abspath(link_ref_pdb):
            if os.path.exists(link_ref_pdb): os.unlink(link_ref_pdb)
            if os.path.exists(link_ref_mtz): os.unlink(link_ref_mtz)
        # ==============================>
        # Create links to dataset
        # ==============================>
        if not os.path.exists(link_ref_pdb): rel_symlink(orig=ref_pdb_path, link=link_ref_pdb)
        if not os.path.exists(link_ref_mtz): rel_symlink(orig=ref_mtz_path, link=link_ref_mtz)
        # ==============================>
        # Create and set reference dataset
        # ==============================>

        # ref_dataset = PanddaReferenceDataset.from_file(
        #     model_filename=os.path.relpath(link_ref_pdb, start=self.out_dir),
        #     data_filename=os.path.relpath(link_ref_mtz, start=self.out_dir)).label(num=-1, tag='reference')
        ref_dataset = PanddaReferenceDataset.from_file(
            model_filename=str(link_ref_pdb),
            data_filename=str(link_ref_mtz)).label(num=-1, tag='reference')

        # ==============================>
        # Extract reference dataset SFs
        # ==============================>
        sf_cols = [[str(x) for x in self.reference_structure_factors.split(',')]] if self.reference_structure_factors \
            else [[str(x) for x in sf.split(',')] for sf in self.structure_factors]

        # Record when a pair is found
        dataset_sfs = None
        # Extract mtz object from the reference dataset
        mtz_obj = ref_dataset.data.mtz_object()
        # Iterate through possible structure factor pairs
        for sf_pair in sf_cols:
            # Check that the data contains the appropriate column
            if mtz_obj.has_column(sf_pair[0]) and mtz_obj.has_column(sf_pair[1]):
                dataset_sfs = ','.join(sf_pair)
                break
        # Raise error if no columns are identified
        if dataset_sfs is None:
            raise Sorry(
                'No matching structure factors were found in the reflection data for reference dataset. \n' + \
                'Looking for structure factors: \n\t{}\n'.format('\n\t'.join(map(' and '.join, sf_cols))) + \
                'Structure factors in this dataset: \n\t{}\n'.format('\n\t'.join(mtz_obj.column_labels())) + \
                'You may need to change the diffraction_data.structure_factors or the reference.structure_factors option.')
        # Store column labels for later
        ref_dataset.meta.column_labels = dataset_sfs
        # Load the diffraction data
        ref_dataset.data.miller_arrays[dataset_sfs] = ref_dataset.data.get_structure_factors(columns=dataset_sfs)

        return ref_dataset

    def __repr__(self):
        repr = {}
        return repr
