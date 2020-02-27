from collections import OrderedDict
import pathlib as p

from pandda_2.native_map_maker import NativeMapMaker


class MakeEventMap:
    def __init__(self):
        pass

    def __call__(self, xmap, truncated_dataset, ref_map, event, bdc, event_map_path, grid):
        # Calculate event maps and output
        event_map_maker = PanddaNativeEventMapMaker(grid,
                                                    grid.get_mask("protein")._max_dist,
                                                    grid.grid_spacing(),
                                                    )

        # gc.collect()
        print("Outputing event map")
        file_path = p.Path(event_map_path)

        if not file_path.exists():
            file_path_string = str(file_path)
            event_map_maker(truncated_dataset,
                            xmap,  # sample
                            ref_map,
                            event,
                            bdc,
                            str(file_path_string),
                            )

    def repr(self):
        repr = OrderedDict()
        return repr


class PanddaNativeEventMapMaker:

    def __init__(self, grid, outer_mask, grid_spacing):
        self.grid = grid
        self.outer_mask = outer_mask
        self.grid_spacing = grid_spacing

    def __call__(self, dataset, sample, ref_map, event, bdc, file_path):
        # ============================================================================================= #
        #####                                                                                       #####
        #               Generate native-aligned maps (in the crystallographic unit cell)                #
        #####                                                                                       #####
        # ============================================================================================= #
        # ============================================================================>
        # Make Event-map
        # ============================================================================>

        ref_event_map = (sample - (ref_map * bdc))
        map_maker = NativeMapMaker(dataset=dataset,
                                   map_obj=ref_event_map,
                                   sites_mask=self.grid.global_mask().sites_cart,
                                   filename=file_path,
                                   outer_mask=self.outer_mask,
                                   grid_spacing=self.grid_spacing
                                   )

        map_maker.run()
