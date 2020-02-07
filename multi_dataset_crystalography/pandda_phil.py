import libtbx.phil

pandda_phil = libtbx.phil.parse("""
pandda
{
    input
        .help = "File names"
    {
        data_dirs = None
            .type = str
        pdb_style = final.pdb
            .type = str
        mtz_style = None
            .type = str
        lig_style = *.cif
            .type = str
        max_new_datasets = 500
            .help = "Maximum number of new datasets that can be processed for an analysis -- used to break the analysis into chunks when processing large numbers of datasets"
            .type = int
        regex
            .help = "Advanced dataset labelling regexs"
        {
            dir_regex = None
                .type = str
            pdb_regex = None
                .type = str
            mtz_regex = None
                .type = str
        }
        filter
            .help = "Provide a dataset to filter input datasets against"
        {
            pdb = None
                .type = path
        }
        reference
            .help = "Manually define reference dataset to align+scale all other datasets to"
        {
            pdb = None
                .type = path
            mtz = None
                .type = path
            structure_factors = "FTW,PHWT"
                .type = str
        }
        flags
            .help = "Flags for individual datasets"
        {
            ground_state_datasets = None
                .help = "Define a set of known 'ground-state' (e.g. unbound) datasets. ONLY these datasets will be used for characterising the ground-state electron density."
                .type = str
            exclude_from_z_map_analysis = None
                .help = "Don't analyse these datasets, only use them to build the distributions - comma separated list of dataset tags"
                .type = str
            exclude_from_characterisation = None
                .help = "Don't use these datasets to characterise the ground-state electron density, only analyse them - comma separated list of dataset tags"
                .type = str
            only_datasets = None
                .help = "Only load these datasets, don't load any others - comma separated list of dataset tags"
                .type = str
            ignore_datasets = None
                .help = "Reject these datasets, don't even load them - comma separated list of dataset tags"
                .type = str
            reprocess_datasets = None
                .help = "Selection of existing datasets to reproces (treat as new datasets) - comma separated list of dataset tags. Setting this will set flags.existing_datasets=reload."
                .type = str
            test = None
                .help = "Only load these datasets, don't load any others - comma separated list of dataset tags"
                .type = str
            train = None
                .help = "Only load these datasets, don't load any others - comma separated list of dataset tags"
                .type = str
            not_test = None
                .help = "Only load these datasets, don't load any others - comma separated list of dataset tags"
                .type = str
            not_train = None
                .help = "Only load these datasets, don't load any others - comma separated list of dataset tags"
                .type = str
        }
    }
    output
        .help = "Output directory"
    {
        out_dir = './pandda'
            .type = path
        new_analysis_dir = False
            .help = "Create a new analysis directory to prevent overwriting previous results"
            .type = bool
        dataset_prefix = ''
            .help = "Prefix to be attached to each dataset name"
            .type = str
        maps
            .help = "Control which maps are output by the program"
        {
            write_z_maps = none *interesting all
                .help = "Output the z-maps in the native frame of datasets for selected datasets"
                .type = choice
            write_average_map = none *interesting all
                .help = "Output the average map in the native frame of datasets for selected datasets"
                .type = choice
            write_dataset_map = *none interesting all
                .help = "Output the analysed maps in the native frame of datasets for selected datasets"
                .type = choice
            write_statistical_maps = *none reference
                .help = "Output statistical maps in the native frame of datasets"
                .type = choice
        }
        pickling
            .help = "Pickle Settings"
            .expert_level = 1
        {
            pickle_complete_pandda = False
                .help = "Pickle the entire PANDDA object for reloading - takes a lot of space"
                .type = bool
            pickle_map_analysers = False
                .help = "Pickle map analysers"
                .type = bool
            pickle_dataset_maps = False
                .help = "Pickle the dataset maps as part of the dataset pickle object"
                .type = bool
        }
        developer
            .help = "Development/testing settings (Not needed by most users)"
            .expert_level = 3
        {
            debug = False
                .help = "Activate all developer flags"
                .type = bool
            write_grid_frame_masks = False
                .help = "Output a selection of maps (e.g. grid masks) in the aligned grid coordinate frame"
                .type = bool
            write_reference_frame_common_masks_and_maps = False
                .help = "Output the reference dataset maps and grid masks in the reference coordinate frame"
                .type = bool
            write_reference_frame_statistical_maps = False
                .help = "Output the statistical maps in the reference coordinate frame"
                .type = bool
            write_reference_frame_dataset_masks_and_maps = False
                .help = "Output maps for datasets in the reference coordinate frame"
                .type = bool
            write_reference_frame_all_dataset_z_map_types = False
                .help = "Output all possible types of Z-maps in the reference coordinate frame"
                .type = bool
        }
    }
    flags
        .help = "control which datasets are loaded and processed, and when statistical maps are calculated"
    {
        stages = *variation_analysis *z_map_analysis
            .help = "Which parts of the program should be turned on?
                        variation_analysis:  perform statistical density variation analysis (analyse differences between datasets).
                        z_map_analysis:      contrast datasets against characterised density (identify events in each dataset)."
            .type = choice(multi=True)
        existing_datasets = reprocess *reload ignore
            .help = "What to do with previously-analysed datasets? reprocess: old datasets are treated as new and processed fully. reload: events identified in old datasets will be included in results. ignore: ..."
            .type = choice
        recalculate_statistical_maps = yes no *extend
            .help = "Set whether statistical maps are re-used from previous runs. If No, it looks for existing statistical maps and uses those (reverts to Yes if none are found). If Extend is chosen, existing maps are used, but additional maps are calculated at high and low resolutions if the data extends beyond the current range."
            .type = choice
        density_analysis_for = all_resolutions *datasets
            .help = "Select which resolutions density is analysed at - characterise the density always? or only when there is a dataset to be analysed?"
            .type = choice
    }
    shortcuts
        .help = "Shortcuts to set sets of parameters to defaults"
    {
        run_in_single_dataset_mode = False
            .help = "Set the default parameters to allow the characterisation to be performed using a single dataset (no variation analysis)"
            .type = bool
    }
    params
        .help = "Algorithm Parameters"
    {
        analysis
            .help = "Settings to control the selection of datasets"
        {
            dynamic_res_limits = True
                .help = 'Allow the analysed resolution limits to change depending on the dataset resolution ranges'
                .type = bool
            high_res_upper_limit = 0.0
                .help = 'Highest resolution limit (maps are never calulcated above this limit)'
                .type = float
            high_res_lower_limit = 4.0
                .help = 'Lowest resolution limit (datasets below this are ignored)'
                .type = float
            high_res_increment = 0.05
                .help = 'Increment of resolution shell for map analysis'
                .type = float
        }
        diffraction_data
        {
            structure_factors = FWT,PHWT
                .type = strings
            checks
                .help = "checks on the mtz file data provided for each dataset"
            {
                all_data_are_valid_values = true
                    .help = "check that all reflections in the diffraction data have valid values (are not zero or n/a)"
                    .type = bool
                low_resolution_completeness = 4.0
                    .help = "check that diffraction data is 100% complete up to this resolution cutoff. missing reflections at low resolution may seriously degrade analysis quality. set to none to turn off this check."
                    .type = float
            }
            apply_b_factor_scaling = True
                .help = "Apply b-factor scaling to reflections? (reciprocal space)"
                .type = bool
        }
        alignment
            .help = "Settings to control the alignment of the structures"
        {
            method = global *local
                .help = "How should the structures be aligned? 'global' is fast, but requires high structural conservation, whereas local is slower and accounts for structural variation"
                .type = choice
        }
        pre_pandda
            .help = "setting to control pre-PanDDA dataset clustering"
        {
            do_pre_pandda = False
                .help = "Whether or not to do a pre-PanDDA clustering of datasets"
                .type = bool
        }
        
        
        filtering
            .help = "Settings to control when datasets are rejected from the analysis"
        {
            max_rmsd_to_reference = 1.5
                .help = "Reject datasets that have a calpha rmsd of greater than this to the reference (after global alignment)"
                .type = float
                .multiple = False
            max_rfree = 0.4
                .help = 'Maximum allowed rfree for a structure (datasets above this are rejected)'
                .type = float
                .multiple = False
            flags {
                same_space_group_only = True
                    .help = "Reject datasets that are a different spacegroup to the reference/filter dataset - NOT YET IMPLEMENTED - MUST BE SET TO TRUE"
                    .type = bool
                similar_models_only = False
                    .help = "Reject datasets that have a different model composition to the reference/filter dataset. All models must have the same number and identity of atoms."
                    .type = bool
            }
        }
        excluding
            .help = "Parameters to control when datasets are automatically excluded from characterisation"
        {
            max_wilson_plot_z_score = 5.0
                .help = "Maximum Z-score for RMSD of dataset diffraction wilson curve compared to the reference dataset. Z-scores calculated relative to the rest of the datasets."
                .type = float
                .multiple = False
        }
        maps
            .help = "Settings to control how maps are generated and analysed"
        {
            resolution_factor = 0.25
                .help = 'Sampling factor for fft-ing the maps'
                .type = float
            grid_spacing = 0.5
                .help = 'Spacing of the grid points in the sampled maps (A) - fixed across resolutions'
                .type = float
            padding = 3
                .help = "Padding around the edge of the maps (A)"
                .type = float
            density_scaling = none *sigma volume
                .help = "Apply scaling to electron density? (real-space)"
                .type = choice
        }
        masks
            .help = "Parameters to control the masking of grid points around the protein"
        {
            pdb = None
                .help = "A PDB to mask the grid against (if none provided, use reference dataset)"
                .type = str
                .multiple = False
            align_mask_to_reference = False
                .help = "If masks.pdb is supplied, does it require alignment to the reference structure?
                            If selecting a fragment of the structure, masks.pdb must already be aligned prior to running pandda (can't align fragments)."
                .type = bool
            selection_string = None
                .help = "Define a custom region of the protein to load electron density for (selection + outer_mask - inner_mask). If not defined, uses all protein atoms."
                .type = str
            outer_mask = 6
                .help = 'include region within outer_mask of the atoms defined by selection_string (or protein if selection_string undefined)'
                .type = float
            inner_mask = 1.8
                .help = 'exclude region within inner_mask of the atoms defined by selection_string (or protein if selection_string undefined)'
                .type = float
            inner_mask_symmetry = 3.0
                .help = 'exclude region within inner_mask of the (symmetry) atoms defined by selection_string (or protein if selection_string undefined)'
                .type = float
        }
        statistical_maps
            .help = "Settings to control the calculation of z-maps"
        {
            min_build_datasets = 40
                .help = 'Minimum number of datasets needed to build distributions'
                .type = int
            max_build_datasets = 60
                .help = 'Maximum number of datasets used to build distributions'
                .type = int
            average_map = *mean_map medn_map
                .help = 'Which statistical map should the uncertainties and Z-map be calculated from?
                            mean_map: mean map      (simple averaging)
                            medn_map: median map    (less sensitive to outliers)'
                .type = choice
            z_map_type = uncertainty *adjusted+uncertainty quantile
                .help = 'Type of Z-map to calculate'
                .type = choice
        }
        statistics
            .help = "Settings to control the statistical model in PanDDA"
        {
            numba = False
                .help = "Whether or not to use numba accelerated calculation of statistical model"
                .type = bool
        }
        
        z_map_analysis
            .help = "Settings to control the finding of blobs"
            .expert_level = 1
        {
            clustering_method = agglomerative_hierarchical* DBSCAN HDBSCAN 
                .help = "Method of clustering to use"
                .type = choice
            fofc_contouring = False
                .help = "Whether or not to use fofc contouring to find blob size"
                .type = bool
            contour_level = 2.5
                .help = 'Contour level when looking for blobs'
                .type = float
            negative_values = False
                .help = 'Look for large negative values as well?'
                .type = bool
            min_blob_volume = 10.0
                .help = 'Blob volume filter for detecting blobs'
                .type = float
            min_blob_z_peak = 3.0
                .help = 'Blob Z-peak filter for detecting blobs'
                .type = float
            all_clusters = False
                .help = "Record all clusters"
                .type = bool
            write_clusters = False
                .help = "Write the cluster indicies to hdf5 files"
                .type = bool
            masks {
                selection_string = None
                    .help = 'Define a further mask to define a region for z-map analysis (selection + outer_mask - inner_mask). If not defined, uses the same region as defined in params.masks.'
                    .type = str
                outer_mask = None
                    .help = 'include region within outer_mask of the atoms defined by selection_string (or protein if selection_string undefined)'
                    .type = float
                inner_mask = None
                    .help = 'exclude region within inner_mask of the atoms defined by selection_string'
                    .type = float
            }
            agglomerative_hierarchical {
                clustering_cutoff = 1.732
                    .help = "cutoff range for agglomerative clustering"
                    .type = float
            }
        }
        background_correction
            .help = "Parameters to control the estimation of feature background corrections"
            .expert_level = 3
        {
            min_bdc = 0.0
                .type = float
                .help = 'Minimum background correction estimate'
            max_bdc = 1.0
                .type = float
                .help = 'Maximum background correction estimate'
            increment = 0.01
                .type = float
                .help = 'Resolution of background correction estimation'
            output_multiplier = 1.0
                .type = float
                .help = 'Empirical multiplier to be applied to the contrast-estimated value of 1-BDC (truncated after multiplication to range 0-1)'
        }
    }
    
    processing
        .help = "Change how the PanDDA is calculated"
    {
        backend = seriel parallel_multiprocess parallel_joblib parallel_dask *distributed_dask
            .help = "What backend to use for PanDDA math"
            .type = choice
    }
    
    results
        .help = "Change ordering/filtering of the output data"
    {
        events {
            order_by = z_peak z_mean *cluster_size
                .help = "How should events be ordered within each site?"
                .type = choice
        }
        sites {
            order_by = *num_events
                .help = "How should sites be ordered?"
                .type = choice
        }
    }
    exit_flags
        .help = "Flags for terminating the program early"
    {
        dry_run = False
            .help = "Find input files on disk, but exit before loading any datasets."
            .type = bool
        setup_only = False
            .help = "Load datasets and create grid partition, but exit before analysis."
            .type = bool
        calculate_first_average_map_only = False
           .help = "Will calculate the highest resolution average map and then exit - used for initial reference modelling."
            .type = bool
    }
}
include scope giant.phil.settings_phil
""", process_includes=True)

