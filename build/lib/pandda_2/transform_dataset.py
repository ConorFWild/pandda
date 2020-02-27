from pandda_2.functional import curry_kwargs, chain


class TransformDataset:
    def __init__(self,
                 transform_data_check,
                 transform_scale_diffraction,
                 transform_filter_structure,
                 transform_filter_wilson,
                 transform_align,
                 ):
        self.transform_data_check = transform_data_check
        self.transform_scale_diffraction = transform_scale_diffraction
        self.transform_filter_structure = transform_filter_structure
        self.transform_filter_wilson = transform_filter_wilson
        self.transform_align = transform_align

    def __call__(self, dataset, reference):
        print("\tTransform: data check")
        dataset = self.transform_data_check(dataset,
                                            reference=reference,
                                            )

        print("\tTransform: scale diffraction")
        dataset = self.transform_scale_diffraction(dataset,
                                                   reference_dataset=reference,
                                                   )
        print("\tTransform: filtering structure")
        dataset = self.transform_filter_structure(dataset,
                                                  reference_dataset=reference,
                                                  )

        print("\tTransform: filtering Wilson")
        dataset = self.transform_filter_wilson(dataset,
                                               reference=reference,
                                               )

        print("\tTransform: aligning dataset")
        dataset = self.transform_align(dataset,
                                       reference_dataset=reference,
                                       )

        return dataset

    def __repr__(self):
        repr = {}
        repr["transform_data_check"] = self.transform_data_check.repr()
        repr["transform_scale_diffraction"] = self.transform_scale_diffraction.repr()
        repr["transform_filter_structure"] = self.transform_filter_structure.repr()
        repr["transform_filter_wilson"] = self.transform_filter_wilson.repr()
        repr["transform_align"] = self.transform_align.repr()

        return repr
