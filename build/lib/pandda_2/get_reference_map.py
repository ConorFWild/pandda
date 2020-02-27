from collections import OrderedDict


class GetReferenceMap:
    def __init__(self):
        pass

    def __call__(self,
                 reference_map_getter,
                 reference,
                 map_resolution,
                 grid,
                 ):
        ref_map = reference_map_getter(reference,
                                       map_resolution,
                                       grid,
                                       )

        return ref_map

    def repr(self):
        repr = OrderedDict()
        return repr
