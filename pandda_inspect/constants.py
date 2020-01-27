
REJECT_MASK_NAMES    = [    'rejected - non-identical structures',
                            'rejected - rmsd to reference',
                            'rejected - wilson plot rmsd',
                            'rejected - different space group',
                            'rejected - isomorphism cutoff',
                            'rejected - data quality',
                            'rejected - rfree',
                            'rejected - unknown'  ]
FLAG_MASK_NAMES      = [    'rejected - total',
                            'noisy zmap',
                            'analysed',
                            'interesting',
                            'exclude_from_z_map_analysis',
                            'exclude_from_characterisation'    ]

DATASET_INFO_FIELDS  = [    'high_resolution',
                            'low_resolution',
                            'r_work',
                            'r_free',
                            'rmsd_to_reference',
                            'space_group',
                            'uc_a',     'uc_b',    'uc_c',
                            'uc_alpha', 'uc_beta', 'uc_gamma',
                            'uc_vol',
                            'unscaled_wilson_B',   'scaled_wilson_B',
                            'applied_b_factor_scaling',
                            'unscaled_wilson_rmsd_all',   'unscaled_wilson_rmsd_>4A',   'unscaled_wilson_rmsd_<4A',
                            'unscaled_wilson_rmsd_all_z', 'unscaled_wilson_rmsd_<4A_z', 'unscaled_wilson_rmsd_>4A_z',
                            'unscaled_wilson_ln_rmsd',    'unscaled_wilson_ln_rmsd_z',
                            'unscaled_wilson_ln_dev',     'unscaled_wilson_ln_dev_z',
                            'scaled_wilson_rmsd_all',     'scaled_wilson_rmsd_>4A',     'scaled_wilson_rmsd_<4A',
                            'scaled_wilson_rmsd_all_z',   'scaled_wilson_rmsd_<4A_z',   'scaled_wilson_rmsd_>4A_z',
                            'scaled_wilson_ln_rmsd',      'scaled_wilson_ln_rmsd_z',
                            'scaled_wilson_ln_dev',       'scaled_wilson_ln_dev_z',
                            'rejection_reason'                ]
DATASET_MAP_FIELDS   = [    'analysed_resolution',
                            'map_uncertainty',
                            'obs_map_mean',
                            'obs_map_rms',
                            'z_map_mean',
                            'z_map_std',
                            'z_map_skew',
                            'z_map_kurt'            ]
DATASET_EVENT_FIELDS = [    'site_idx',
                            '1-BDC',
                            'z_peak',
                            'z_mean',
                            'cluster_size',
                            'x','y','z',
                            'refx','refy','refz',
                            'global_correlation_to_average_map',
                            'local_correlation_to_average_map' ]
SITE_TABLE_FIELDS    = [    'centroid',
                            'num_events',
                            'nearest_residue 1',
                            'nearest_residue 2',
                            'nearest_residue 3',
                            'native_centroid',
                            'near_crystal_contacts'     ]

class PanddaMaskNames:
    reject_mask_names    = REJECT_MASK_NAMES
    flag_mask_names      = FLAG_MASK_NAMES
    all_mask_names       = flag_mask_names + reject_mask_names
    write_mask_names     = flag_mask_names

class PanddaTableFields:
    all_dataset_fields     = DATASET_INFO_FIELDS
    all_dataset_map_fields = DATASET_MAP_FIELDS
    all_event_fields       = DATASET_EVENT_FIELDS
    all_site_fields        = SITE_TABLE_FIELDS

########################################################################

class PanddaDatasetFilenames:
    # Dataset Information
    dataset_info           = '{!s}-info.csv'
    dataset_log            = '{!s}.log'
    dataset_pickle         = 'dataset.pickle'
    # Input data...
    input_model            = '{!s}-pandda-input.pdb'
    input_data             = '{!s}-pandda-input.mtz'
    # Structure Files...
    aligned_model          = '{!s}-aligned.pdb'
    symmetry_model         = '{!s}-aligned-sym-contacts.pdb'
    # Native (Output) Maps
    native_obs_map = '{!s}-observed.native.ccp4'
    native_z_map = '{!s}-z_map.native.ccp4'
    native_event_map = '{!s}-event_{!s}_1-BDC_{!s}_map.native.ccp4'
    native_average_map = '{!s}-ground-state-average-map.native.ccp4'
    # Modelled Structures...
    modelled_structure     = '{!s}-pandda-model.pdb'
    ensemble_structure     = '{!s}-ensemble-model.pdb'
    # Output Maps
    input_map              = '{!s}-input-map.ccp4'
    zstat_map              = '{!s}-z_map.ccp4'
    event_map              = '{!s}-event_{!s}_1-BDC_{!s}_map.ccp4'
    average_map            = '{!s}-ground-state-average-map.ccp4'
    average_diff_map       = '{!s}-difference-from-mean.ccp4'
    # Different Z-map types
    z_map_uncertainty      = '{!s}-z_map_uncertainty.ccp4'
    z_map_uncertainty_norm = '{!s}-z_map_uncertainty_normalised.ccp4'
    z_map_corrected        = '{!s}-z_map_adjusted.ccp4'
    z_map_corrected_norm   = '{!s}-z_map_adjusted_normalised.ccp4'
#    z_map_naive            = '{!s}-z_map_naive.ccp4'
#    z_map_naive_norm       = '{!s}-z_map_naive_normalised.ccp4'
    # Grid masks files
    grid_mask              = '{!s}-masked_grid.ccp4'
    z_map_mask             = '{!s}-z_map_search_mask.ccp4'
    high_z_mask            = '{!s}-z_map_high_z_mask.ccp4'
    # Ligands files
    ligand_coordinates     = '{!s}-ligand.pdb'
    ligand_restraints      = '{!s}-ligand.cif'
    ligand_image           = '{!s}-ligand.png'
    # Miscellaneous files
    z_peaks_csv            = '{!s}-z_map_peaks.csv'
    pymol_script           = 'load_script_for_pymol.py'
    ccp4mg_script          = 'load_script_for_ccp4mg.py'

class PanddaDatasetPNGFilenames:
    s_map_png                   = '{!s}-map_histogram.png'
    d_mean_map_png              = '{!s}-mean_map_difference_histogram.png'
    obs_qqplot_sorted_png       = '{!s}-mean_map_scatter_sorted.png'
    obs_qqplot_unsorted_png     = '{!s}-mean_map_scatter_unsorted.png'
#    z_map_naive_png             = '{!s}-z_map_histogram_naive.png'
#    z_map_naive_norm_png        = '{!s}-z_map_histogram_naive_normalised.png'
    z_map_uncertainty_png       = '{!s}-z_map_histogram_uncertainty.png'
    z_map_uncertainty_norm_png  = '{!s}-z_map_histogram_uncertainty_normalised.png'
    z_map_corrected_png         = '{!s}-z_map_histogram_adjusted.png'
    z_map_corrected_norm_png    = '{!s}-z_map_histogram_adjusted_normalised.png'
    z_map_qq_plot_png           = '{!s}-qq_plot_z_map.png'
    unc_qqplot_png              = '{!s}-qq_plot_map.png'
    bdc_est_png                 = '{!s}-event_{!s}_bdc_estimation.png'
    wilson_plot                 = '{!s}-wilson_plot.png'

class PanddaAnalyserFoldernames:
    pass

class PanddaAnalyserFilenames:

    event_info              = 'pandda_analyse_events.csv'
    site_info               = 'pandda_analyse_sites.csv'
    dataset_info            = 'all_datasets_info.csv'
    dataset_map_info        = 'all_datasets_info_maps.csv'
    dataset_combined_info   = 'all_datasets_info_combined.csv'
    dataset_masks           = 'all_datasets_info_masks.csv'

    mean_map                = '{!s}A-mean_map.ccp4'
    medn_map                = '{!s}A-medn_map.ccp4'
    stds_map                = '{!s}A-stds_map.ccp4'
    sadj_map                = '{!s}A-sadj_map.ccp4'
    skew_map                = '{!s}A-skew_map.ccp4'
    kurt_map                = '{!s}A-kurt_map.ccp4'
    bimo_map                = '{!s}A-bimo_map.ccp4'

    reference_structure     = 'reference.pdb'
    reference_dataset       = 'reference.mtz'
    reference_symmetry      = 'reference-symmetry.pdb'

    reference_grid_mask     = 'grid-voronoi-{}.ccp4'

class PanddaInspectorFilenames:

    event_info              = 'pandda_inspect_events.csv'
    site_info               = 'pandda_inspect_sites.csv'

class PanddaHtmlFilenames:

    initial_html            = 'pandda_initial.html'

    map_html                = 'pandda_map_{:.2f}A.html'

    analyse_html            = 'pandda_analyse.html'
    analyse_site_graph      = 'pandda_analyse_sites_graph.png'
    analyse_site_graph_mult = 'pandda_analyse_sites_graph_{!s}.png'
    pymol_sites_png_1       = 'pandda_analyse_sites_pymol_1.png'
    pymol_sites_png_2       = 'pandda_analyse_sites_pymol_2.png'
    pymol_sites_py          = 'pandda_analyse_sites_pymol.py'
    pymol_sites_pml         = 'pandda_analyse_sites_pymol.pml'

    inspect_html            = 'pandda_inspect.html'
    inspect_site_graph      = 'pandda_inspect_sites_graph.png'
    inspect_site_graph_mult = 'pandda_inspect_sites_graph_{!s}.png'
