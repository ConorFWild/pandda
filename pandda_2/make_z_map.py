from collections import OrderedDict
import pathlib as p

from pandda_2.native_map_maker import NativeMapMaker


class MakeZMap:
    def __init__(self):
        pass

    def __call__(self, xmap, truncated_dataset, z_map_path, statistical_model, grid):
        statistical_map = statistical_model.evaluate(xmap)  # sample)
        statistical_map_maker = PanddaNativeStatisticalMapMaker(grid,
                                                                grid.get_mask("protein")._max_dist,
                                                                grid.grid_spacing(),
                                                                )
        file_path = p.Path(z_map_path)
        if not file_path.exists():
            file_path_string = str(file_path)
            print("Outputing statistical map")
            statistical_map_maker(truncated_dataset,
                                  statistical_map,
                                  file_path_string,
                                  )

    def repr(self):
        repr = OrderedDict()
        return repr


class PanddaNativeStatisticalMapMaker:

    def __init__(self, grid, outer_mask, grid_spacing):
        self.grid = grid
        self.outer_mask = outer_mask
        self.grid_spacing = grid_spacing

    def __call__(self, dataset, ref_z_map, file_path):
        # ============================================================================================= #
        #####                                                                                       #####
        #               Generate native-aligned maps (in the crystallographic unit cell)                #
        #####                                                                                       #####
        # ============================================================================================= #
        # ============================================================================>
        # Make Event-map
        # ============================================================================>
        ref_z_map = ref_z_map.normalised_copy()
        map_maker = NativeMapMaker(dataset=dataset,
                                   map_obj=ref_z_map,
                                   sites_mask=self.grid.global_mask().sites_cart,
                                   filename=file_path,
                                   outer_mask=self.outer_mask,
                                   grid_spacing=self.grid_spacing,
                                   )

        map_maker.run()
