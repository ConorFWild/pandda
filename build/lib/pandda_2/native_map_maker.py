from giant.xray.maps.local_align import create_native_map

class NativeMapMaker(object):

    def __init__(self, dataset, map_obj, sites_mask, filename, outer_mask, grid_spacing):
        """
        The main object for comparing each dataset-map to the ensemble-maps.
        Constructed so that the class can be initialised and then called within a multiprocessing function.

        Usage:
            new_proc = DatasetProcessor(...)
            output   = new_proc.run()
        or:
            output   = DatasetProcessor.process(...)
        """

        self.data = (dataset, map_obj, sites_mask, filename)
        self.outer_mask = outer_mask
        self.grid_spacing = grid_spacing

    @classmethod
    def process(cls, dataset, map_obj, sites_mask, filename):
        """Process the dataset immediately and return output"""
        return cls(dataset=dataset, map_obj=map_obj, sites_mask=sites_mask, filename=filename).run()

    def run(self):
        """Process the dataset"""

        dataset, map_obj, sites_mask, filename = self.data

        native_map_data = create_native_map(
            native_crystal_symmetry=dataset.model.crystal_symmetry,
            native_sites=dataset.model.alignment.ref2nat(sites_mask),
            alignment=dataset.model.alignment,
            reference_map=map_obj.make_dense(),
            site_mask_radius=self.outer_mask,
            step=self.grid_spacing,
            filename=filename
        )

        return native_map_data
