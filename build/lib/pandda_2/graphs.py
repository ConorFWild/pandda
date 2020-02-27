import os, sys

#################################
try:
    import matplotlib
    print("Setting backend in graphs")
    matplotlib.use('tkagg')
    matplotlib.interactive(False)
    from matplotlib import pyplot
    pyplot.style.use('ggplot')
except:
    pass
#################################

import numpy

from scitbx.array_family import flex

from bamboo.plot import bar, simple_histogram

def filter_nans(x):
    return [v for v in x if not numpy.isnan(v)]

def failure_graph(title):
    fig = pyplot.figure()
    pyplot.title('Failed to make {}'.format( title ))
    return fig

#################################

def get_wilson_plot_vals(miller_array):
    """Get the atuomatic wilson plot values for an array"""
    # Create binning - just use the auto-binning for each dataset
    binner = miller_array.setup_binner(auto_binning=True)
    # Create the wilson plot
    binned = miller_array.wilson_plot(use_binning=True)
    x_bin_cent = binner.bin_centers(1)
    y_bin_data = binned.data[1:-1]
    assert len(x_bin_cent) == len(y_bin_data)
    return numpy.power(x_bin_cent,2), numpy.log(y_bin_data)

#################################

def map_value_distribution(f_name, plot_vals, plot_normal=False):
    """Plot histogram of values, with optional normal distribution"""
    from scitbx.math.distributions import normal_distribution
    fig = pyplot.figure()
    pyplot.title('Distribution of map values')
    pyplot.hist(x=plot_vals, bins=30, normed=True)
    if plot_normal:
        # Plot the distribution for N(0,1)
        nd_t = normal_distribution()
        theor_x = numpy.linspace(-5,5,101)
        theor_y = [nd_t.pdf(x) for x in theor_x]
        pyplot.plot(theor_x, theor_y, c='k', ls='--', marker='o')
        # Plot the distribution for the observed distribution
        nd_o = normal_distribution(mean=numpy.mean(plot_vals), sd=numpy.std(plot_vals))
        obs_x = numpy.linspace(-5,5,101)
        obs_y = [nd_o.pdf(x) for x in obs_x]
        pyplot.plot(obs_x, obs_y, c='g', ls='-', marker='o')
    pyplot.xlabel('Map value')
    pyplot.ylabel('Density')
    #pyplot.tight_layout()
    pyplot.subplots_adjust()
    pyplot.savefig(f_name)
    pyplot.close(fig)

def qq_plot_against_normal(f_name, plot_vals):
    """Sort and plot list of values against expected quantiles from a normal distribution"""
    from scitbx.math.distributions import normal_distribution
    fig = pyplot.figure()
    pyplot.title('Q-Q plot for map values against normal distribution')
    expected_vals = normal_distribution().quantiles(len(plot_vals))
    pyplot.plot([min(expected_vals)-1, max(expected_vals)+1], [min(expected_vals)-1, max(expected_vals)+1], 'b--')
    pyplot.plot(sorted(plot_vals), expected_vals, 'go-')
    pyplot.xlabel('Observed quantiles')
    pyplot.ylabel('Theoretical quantiles')
    #pyplot.tight_layout()
    pyplot.subplots_adjust()
    pyplot.savefig(f_name)
    pyplot.close(fig)

def mean_obs_scatter(f_name, mean_vals, obs_vals):
    """Plot mean map against observed map values"""
    fig = pyplot.figure()
    pyplot.title('Mean map plotted against dataset map (unsorted)')
    pyplot.plot([-3, 10], [-3, 10], 'b--')
    pyplot.plot(mean_vals, obs_vals, 'go')
    pyplot.xlabel('Mean map value')
    pyplot.ylabel('Dataset map value')
    #pyplot.tight_layout()
    pyplot.subplots_adjust()
    pyplot.savefig(f_name)
    pyplot.close(fig)

def sorted_mean_obs_scatter(f_name, mean_vals, obs_vals):
    """Plot sorted mean map against sorted observed map values"""
    fig = pyplot.figure()
    pyplot.title('Mean map plotted against dataset map (sorted)')
    pyplot.plot([-3, 10], [-3, 10], 'b--')
    pyplot.plot(mean_vals, obs_vals, 'go')
    pyplot.xlabel('Mean map value')
    pyplot.ylabel('Dataset map value')
    #pyplot.tight_layout()
    pyplot.subplots_adjust()
    pyplot.savefig(f_name)
    pyplot.close(fig)

def uncertainty_qqplot(f_name, map_off, map_unc, q_cut, obs_diff, quantile):
    """Plot diff-mean map against normal quantiles, with uncertainty lines"""
    fig = pyplot.figure()
    pyplot.title('Q-Q plot for map values against normal distribution')
    pyplot.plot([map_off-5*map_unc, map_off+5*map_unc], [-5, 5], 'b--')
    pyplot.plot([-1, 1], [q_cut, q_cut], 'k-.')
    pyplot.plot([-1, 1], [-q_cut, -q_cut], 'k-.')
    pyplot.plot(obs_diff, quantile, 'go-')
    pyplot.xlabel('Difference from mean map')
    pyplot.ylabel('Theoretical Quantiles')
    #pyplot.tight_layout()
    pyplot.subplots_adjust()
    pyplot.savefig(f_name)
    pyplot.close(fig)

def write_occupancy_graph(f_name, x_values, global_values, local_values):
    """Write output graph from occupancy estimation"""

    # Get the x-value of the maximum difference
    diff_values = list(numpy.array(global_values) - numpy.array(local_values))
    max_x = x_values[diff_values.index(max(diff_values))]

    fig, (axis_1_1, axis_2_1) = pyplot.subplots(2, sharex=True)
    # 1st Plot - 1st Y-Axis
    line_1_1, = axis_1_1.plot(x_values, global_values, 'g--', label='Global corr.', linewidth=2)
    line_1_2, = axis_1_1.plot(x_values, local_values, 'k--', label='Local corr.', linewidth=2)
    axis_1_1.set_ylabel('Corr. to ground state', color='k', size=16)
    axis_1_1.set_ylim((-1, 1))
    # 2nd Plot - 1st Y-Axis
    line_2_1, = axis_2_1.plot(x_values, diff_values, 'b-', label='Difference', linewidth=2)
    axis_2_1.set_ylabel('Corr. difference', color='k', size=16)
    axis_2_1.set_xlabel('1-BDC', color='k', size=16)
    axis_2_1.set_ylim((min(diff_values)-0.2, max(diff_values)+0.2))
    # Plot line at the maximum
    line_2_2, = axis_2_1.plot([max_x,max_x],[-1,1], 'k-', linewidth=2)
    text_2_1 = axis_2_1.text(0.02+max_x, min(diff_values), 'BDC='+str(1-max_x), size=14)
    # Joint legend
    axis_1_1.legend(handles=[line_1_1, line_1_2, line_2_1], loc=4, fontsize=16)
    # Remove spacing between subplots
    pyplot.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
    pyplot.setp([axis_2_1.get_xticklabels()], visible=True)
    #pyplot.tight_layout()
    pyplot.subplots_adjust()
    pyplot.savefig(f_name)
    pyplot.close(fig)

#################################

def write_individual_dataset_plots(pandda, datasets):
    """Write individual wilson plots for each dataset"""

    pandda.log.heading('Generating individual dataset graphs')
    pandda.log('-> wilson plots against the reference dataset')

    ref_d = pandda.datasets.reference()
    ref_x, ref_y = get_wilson_plot_vals(miller_array=ref_d.data.miller_arrays[ref_d.meta.column_labels].as_amplitude_array())

    for d in datasets:
        fig = pyplot.figure()
        pyplot.title('Wilson plot for dataset {}'.format(d.tag))
        # Plot for the reference dataset
        pyplot.plot(ref_x, ref_y, 'k-', linewidth=2, label='reference')
        # Plot for this dataset (unscaled)
        x, y = get_wilson_plot_vals(miller_array=d.data.miller_arrays[d.meta.column_labels].as_amplitude_array())
        pyplot.plot(x, y, 'r-', linewidth=2, label=d.tag+' (unscaled)')
        # Plot for this dataset (scaled)
        x, y = get_wilson_plot_vals(miller_array=d.data.miller_arrays['scaled'].as_amplitude_array())
        pyplot.plot(x, y, 'b-', linewidth=2, label=d.tag+' (scaled)')
        # Plot settings
        pyplot.axes().set_xticklabels([numpy.round(t**-0.5,2) if t >0 else 0 for t in pyplot.axes().get_xticks()])
        pyplot.xlabel('resolution ($\AA$)')
        pyplot.ylabel('ln(mean amplitude)')
        #pyplot.tight_layout()
        pyplot.legend()
        pyplot.subplots_adjust()
        f_name = d.file_manager.get_file('wilson_plot_png')
        pyplot.savefig(f_name)
        pyplot.close(fig)
        pandda.log('\t{}'.format(f_name))

def write_dataset_summary_graphs(pandda):
    """Plot dataset summary graphs of resolution, unit cell variation, etc"""

    pandda.log.heading('Generating summary graphs for non-rejected datasets')

    # Filter the dataset to non-rejected datasets
    non_rejected_dsets = pandda.datasets.mask(mask_name='rejected - total', invert=True)
    non_rejected_dtags = [d.tag for d in non_rejected_dsets]

    d_info = pandda.tables.dataset_info.loc[non_rejected_dtags]
    n_bins = 30

    # ================================================>
    # Simple histograms of data from the dataset_info table
    # ================================================>
    pandda.log('-> crystal variables and data quality histograms')

    # ================================================>
    # High & Low Resolutions
    # ================================================>
    try:
        fig = pyplot.figure()
        pyplot.title('Resolution Histograms')
        pyplot.subplot(2, 1, 1)
        pyplot.hist(x=d_info['high_resolution'], bins=n_bins)
        pyplot.xlabel('High Resolution Limit ($\AA$)')
        pyplot.ylabel('Count')
        pyplot.subplot(2, 1, 2)
        pyplot.hist(x=d_info['low_resolution'], bins=n_bins)
        pyplot.xlabel('Low Resolution Limit ($\AA$)')
        pyplot.ylabel('Count')
        #pyplot.tight_layout()
        pyplot.subplots_adjust(hspace=0.3)
    except:
        fig = failure_graph(title='Resolution Histograms')

    f_name = pandda.file_manager.get_file('d_resolutions')
    pyplot.savefig(f_name)
    pyplot.close(fig)
    pandda.log('\t{}'.format(f_name))
    # ================================================>
    # R-factors
    # ================================================>
    try:
        fig = pyplot.figure()
        pyplot.title('R-Factor Histograms')
        pyplot.subplot(2, 1, 1)
        pyplot.hist(x=d_info['r_free'], bins=n_bins)
        pyplot.xlabel('R-Free')
        pyplot.ylabel('Count')
        pyplot.subplot(2, 1, 2)
        pyplot.hist(x=d_info['r_work'], bins=n_bins)
        pyplot.xlabel('R-Work')
        pyplot.ylabel('Count')
        #pyplot.tight_layout()
        pyplot.subplots_adjust(hspace=0.3)
    except:
        fig = failure_graph(title='R-Factor Histograms')

    f_name = pandda.file_manager.get_file('d_rfactors')
    pyplot.savefig(f_name)
    pyplot.close(fig)
    pandda.log('\t{}'.format(f_name))
    # ================================================>
    # RMSD to reference structure
    # ================================================>
    try:
        fig = pyplot.figure()
        pyplot.title('RMSDs to Reference Structure Histogram')
        pyplot.hist(x=filter_nans(d_info['rmsd_to_reference']), bins=n_bins)
        pyplot.xlabel('RMSD (A)')
        pyplot.ylabel('Count')
        #pyplot.tight_layout()
        pyplot.subplots_adjust()
    except:
        fig = failure_graph(title='RMSDs to Reference Structure Histogram')

    f_name = pandda.file_manager.get_file('d_global_rmsd_to_ref')
    pyplot.savefig(f_name)
    pyplot.close(fig)
    pandda.log('\t{}'.format(f_name))
    # ================================================>
    # Unit cell size
    # ================================================>
    try:
        fig = pyplot.figure()
        pyplot.title('Unit Cell Axis Variation')
        pyplot.subplot(3, 1, 1)
        pyplot.hist(x=d_info['uc_a'], bins=n_bins)
        pyplot.xlabel('A (A)')
        pyplot.subplot(3, 1, 2)
        pyplot.hist(x=d_info['uc_b'], bins=n_bins)
        pyplot.xlabel('B (A)')
        pyplot.subplot(3, 1, 3)
        pyplot.hist(x=d_info['uc_c'], bins=n_bins)
        pyplot.xlabel('C (A)')
        #pyplot.tight_layout()
        pyplot.subplots_adjust(hspace=0.6)
    except:
        fig = failure_graph(title='Unit Cell Axis Variation')

    f_name = pandda.file_manager.get_file('d_cell_axes')
    pyplot.savefig(f_name)
    pyplot.close(fig)
    pandda.log('\t{}'.format(f_name))
    # ================================================>
    # Unit cell angles
    # ================================================>
    try:
        fig = pyplot.figure()
        pyplot.title('Unit Cell Angle Variation')
        pyplot.subplot(3, 1, 1)
        pyplot.hist(x=d_info['uc_alpha'], bins=n_bins)
        pyplot.xlabel('Alpha')
        pyplot.subplot(3, 1, 2)
        pyplot.hist(x=d_info['uc_beta'], bins=n_bins)
        pyplot.xlabel('Beta')
        pyplot.subplot(3, 1, 3)
        pyplot.hist(x=d_info['uc_gamma'], bins=n_bins)
        pyplot.xlabel('Gamma')
        #pyplot.tight_layout()
        pyplot.subplots_adjust(hspace=0.6)
    except:
        fig = failure_graph(title='Unit Cell Angle Variation')

    f_name = pandda.file_manager.get_file('d_cell_angles')
    pyplot.savefig(f_name)
    pyplot.close(fig)
    pandda.log('\t{}'.format(f_name))
    # ================================================>
    # Unit cell volume
    # ================================================>
    try:
        fig = pyplot.figure()
        pyplot.title('Unit Cell Volume Variation')
        pyplot.hist(x=d_info['uc_vol']/1000.0, bins=n_bins)
        pyplot.xlabel('Volume ($10^3 A^3$)')
        pyplot.ylabel('Count')
        #pyplot.tight_layout()
        pyplot.subplots_adjust()
    except:
        fig = failure_graph(title='Unit Cell Volume Variation')

    f_name = pandda.file_manager.get_file('d_cell_volumes')
    pyplot.savefig(f_name)
    pyplot.close(fig)
    pandda.log('\t{}'.format(f_name))

    # ================================================>
    # Summary plots for the loaded diffracton data
    # ================================================>
    pandda.log('-> writing wilson plots of input and scaled data')

    min_x = numpy.nanmin([numpy.min(d_info['unscaled_wilson_rmsd_>4A']), numpy.min(d_info['scaled_wilson_rmsd_>4A'])])
    max_x = numpy.nanmax([numpy.max(d_info['unscaled_wilson_rmsd_>4A']), numpy.max(d_info['scaled_wilson_rmsd_>4A'])])
    min_y = numpy.nanmin([numpy.min(d_info['unscaled_wilson_rmsd_<4A']), numpy.min(d_info['scaled_wilson_rmsd_<4A'])])
    max_y = numpy.nanmax([numpy.max(d_info['unscaled_wilson_rmsd_<4A']), numpy.max(d_info['scaled_wilson_rmsd_<4A'])])

    # ================================================>
    # Wilson plots RMSDS for unscaled (input) data
    # ================================================>
    try:
        fig = pyplot.figure()
        pyplot.title('Wilson Plot RMSD to Reference (unscaled)')
        pyplot.scatter(x=d_info['unscaled_wilson_rmsd_>4A'], y=d_info['unscaled_wilson_rmsd_<4A'])
        pyplot.xlabel('RMSD to Reference (>4A)')
        pyplot.ylabel('RMSD to Reference (<4A)')
        pyplot.xlim(min_x, max_x)
        pyplot.ylim(min_y, max_y)
        #pyplot.tight_layout()
        pyplot.subplots_adjust()
    except:
        raise
        fig = failure_graph(title='Wilson Plot RMSD to Reference (unscaled)')

    f_name = pandda.file_manager.get_file('d_unscaled_wilson_rmsds')
    pyplot.savefig(f_name)
    pyplot.close(fig)
    pandda.log('\t{}'.format(f_name))

    # ================================================>
    # Wilson plots RMSDS for scaled data
    # ================================================>
    try:
        fig = pyplot.figure()
        pyplot.title('Wilson Plot RMSD to Reference (scaled)')
        pyplot.scatter(x=d_info['scaled_wilson_rmsd_>4A'], y=d_info['scaled_wilson_rmsd_<4A'])
        pyplot.xlabel('RMSD to Reference (>4A)')
        pyplot.ylabel('RMSD to Reference (<4A)')
        pyplot.xlim(min_x, max_x)
        pyplot.ylim(min_y, max_y)
        #pyplot.tight_layout()
        pyplot.subplots_adjust()
    except:
        raise
        fig = failure_graph(title='Wilson Plot RMSD to Reference (scaled)')

    f_name = pandda.file_manager.get_file('d_scaled_wilson_rmsds')
    pyplot.savefig(f_name)
    pyplot.close(fig)
    pandda.log('\t{}'.format(f_name))

    # ================================================>
    # Wilson plots of the unscaled (input) data
    # ================================================>
    fig = pyplot.figure()
    pyplot.title('Wilson plot for unscaled (input) structure factors')
    for d in non_rejected_dsets:
        x, y = get_wilson_plot_vals(miller_array=d.data.miller_arrays[d.meta.column_labels].as_amplitude_array())
        pyplot.plot(x, y, '-', linewidth=1)
    pandda.log(pyplot.axes().get_xticks())
    pyplot.axes().set_xticklabels([numpy.round(t**-0.5,2) if t > 0 else 0 for t in pyplot.axes().get_xticks()])
    pyplot.xlabel('resolution ($\AA$)')
    pyplot.ylabel('ln(mean amplitude)')
    #pyplot.tight_layout()
    pyplot.subplots_adjust()
    f_name = pandda.file_manager.get_file('d_unscaled_wilson_plots')
    pyplot.savefig(f_name)
    pyplot.close(fig)
    pandda.log('\t{}'.format(f_name))
    # ================================================>
    # Wilson plots of the scaled data (characterisation)
    # ================================================>
    fig, (axis_1, axis_2) = pyplot.subplots(2, sharex=True)
    pyplot.title('Wilson plot for scaled structure factors')
    for d in non_rejected_dsets:
        x, y = get_wilson_plot_vals(miller_array=d.data.miller_arrays['scaled'].as_amplitude_array())
        ax = axis_1 if pandda.datasets.all_masks().get_value(name='exclude_from_characterisation', id=d.tag)==False else axis_2
        ax.plot(x, y, '-', linewidth=1)
    axis_1.set_title('Datasets to be used for characterisation')
    axis_2.set_title('Datasets excluded from characterisation')
    axis_2.set_xticklabels([numpy.round(t**-0.5,2) if t >0 else 0 for t in axis_2.get_xticks()])
    axis_2.set_xlabel('resolution ($\AA$)')
    axis_1.set_ylabel('ln(mean amplitude)')
    axis_2.set_ylabel('ln(mean amplitude)')
    #pyplot.tight_layout()
    pyplot.subplots_adjust()
    f_name = pandda.file_manager.get_file('d_scaled_wilson_plots')
    pyplot.savefig(f_name)
    pyplot.close(fig)
    pandda.log('\t{}'.format(f_name))

    return None

def write_truncated_data_plots(pandda, resolution, datasets):
    """Write summaries of the truncated data"""

    reslns = [d.data.miller_arrays['truncated'].d_min() for d in datasets]
    min_res, max_res = min(reslns), max(reslns)
    pandda.log('After Truncation - Resolution Range: {!s}-{!s}'.format(min_res, max_res))
    pandda.log('')
    pandda.log('-> Writing truncated data plots')

    # ================================================>
    # Resolution ranges of the truncated data
    f_name = pandda.file_manager.get_file('truncated_res_hist').format(resolution)
    try:
        simple_histogram(filename = f_name,
                         data     = reslns,
                         title    = 'Truncated dataset resolutions',
                         x_lab    = 'Resolution (A)',
                         n_bins   = 15)
    except:
        fig = failure_graph(title='Truncated dataset resolutions')
        pyplot.savefig(f_name)
        pyplot.close(fig)
    pandda.log('\t{}'.format(f_name))

    # ================================================>
    # Wilson plots of the truncated data
    fig = pyplot.figure()
    pyplot.title('Wilson plot for truncated structure factors')
    for d in datasets:
        x, y = get_wilson_plot_vals(miller_array=d.data.miller_arrays['truncated'].as_amplitude_array())
        pyplot.plot(x, y, '-', linewidth=1)
    pyplot.axes().set_xticklabels([numpy.round(t**-0.5,2) if t >0 else 0 for t in pyplot.axes().get_xticks()])
    pyplot.xlabel('resolution ($\AA$)')
    pyplot.ylabel('ln(mean amplitude)')
    #pyplot.tight_layout()
    pyplot.subplots_adjust()
    f_name = pandda.file_manager.get_file('truncated_wilson_plot').format(resolution)
    pyplot.savefig(f_name)
    pyplot.close(fig)
    pandda.log('\t{}'.format(f_name))

def write_map_analyser_reference_dataset_graphs(pandda, map_analyser):

    # Plot the mean map against the reference map (unsorted)
    mean_obs_scatter(f_name      = pandda.file_manager.get_file('ref_v_mean_map_unsort').format(map_analyser.meta.resolution),
                     mean_vals   = map_analyser.statistical_maps.mean_map.get_map_data(sparse=True),
                     obs_vals    = pandda.datasets.reference().child.get_map_data(sparse=True))
    # Plot the mean map against the reference map (sorted)
    sorted_mean_obs_scatter(f_name    = pandda.file_manager.get_file('ref_v_mean_map_sort').format(map_analyser.meta.resolution),
                            mean_vals = sorted(map_analyser.statistical_maps.mean_map.get_map_data(sparse=True)),
                            obs_vals  = sorted(pandda.datasets.reference().child.get_map_data(sparse=True)))
    # Plot the reference map distribution
    map_value_distribution(f_name    = pandda.file_manager.get_file('ref_map_dist').format(map_analyser.meta.resolution),
                           plot_vals = pandda.datasets.reference().child.get_map_data(sparse=True))

def write_map_analyser_graphs(pandda, resolution, analysis_mask_name, building_mask_name):

    # Extract the statistical maps at this resolution
    statistical_maps = pandda.stat_maps.get(resolution)
    # Extract map values
    mean_map_vals = list(statistical_maps.mean_map.get_map_data(sparse=True))
    medn_map_vals = list(statistical_maps.medn_map.get_map_data(sparse=True))
    stds_map_vals = list(statistical_maps.stds_map.get_map_data(sparse=True))
    sadj_map_vals = list(statistical_maps.sadj_map.get_map_data(sparse=True))
    # Extract the dataset tags from the pandda object
    analysis_tags = [d.tag for d in pandda.datasets.mask(mask_name=analysis_mask_name)]
    building_tags = [d.tag for d in pandda.datasets.mask(mask_name=building_mask_name)]
    combined_tags = sorted(set(analysis_tags+building_tags))

    # Dataset Parameters
    d_info = pandda.tables.dataset_info
    m_info = pandda.tables.dataset_map_info

    # All datasets
    high_res = [d_info['high_resolution'][t] for t in combined_tags]
    low_res =  [d_info['low_resolution'][t]  for t in combined_tags]
    rfree =    [d_info['r_free'][t]          for t in combined_tags]
    rwork =    [d_info['r_work'][t]          for t in combined_tags]

    # All datasets
    map_uncties = [m_info['map_uncertainty'][t] for t in combined_tags]
    # Analysed datasets only
    z_map_mean  = [m_info['z_map_mean'][t] for t in analysis_tags]
    z_map_std   = [m_info['z_map_std'][t]  for t in analysis_tags]
    z_map_skew  = [m_info['z_map_skew'][t] for t in analysis_tags]
    z_map_kurt  = [m_info['z_map_kurt'][t] for t in analysis_tags]

    ########################################################

    n_bins = 30

    ########################################################

    pandda.log('=> Writing Statistical Map Distributions')

    ##################################
    # STATISTICAL MAP HISTOGRAMS
    ##################################
    fig = pyplot.figure()
    pyplot.title('Mean and Median Map Histograms')
    # MEAN MAP
    pyplot.subplot(2, 1, 1)
    pyplot.hist(x=mean_map_vals, bins=n_bins)
    pyplot.xlabel('Mean Map Values')
    pyplot.ylabel('Count')
    # MEDIAN MAP
    pyplot.subplot(2, 1, 2)
    pyplot.hist(x=medn_map_vals, bins=n_bins)
    pyplot.xlabel('Median Map Values')
    pyplot.ylabel('Count')
    #pyplot.tight_layout()
    pyplot.subplots_adjust(hspace=0.4)
    pyplot.savefig(pandda.file_manager.get_file('map_mean_median_hist').format(resolution))
    pyplot.close(fig)

    ##################################
    # MEAN Values v MEDIAN Values
    ##################################
    fig = pyplot.figure()
    pyplot.title('Mean v Median Map Values')
    pyplot.scatter(x=mean_map_vals, y=medn_map_vals)
    # Plot straight line between the min and max values
    min_val = min(mean_map_vals+medn_map_vals)
    max_val = max(mean_map_vals+medn_map_vals)
    pyplot.plot([min_val, max_val], [min_val, max_val], 'b--')
    # Axis labels
    pyplot.xlabel('Mean Map Value')
    pyplot.ylabel('Median Map Value')
    #pyplot.tight_layout()
    pyplot.subplots_adjust(hspace=0.4)
    pyplot.savefig(pandda.file_manager.get_file('map_mean_median_scat').format(resolution))
    pyplot.close(fig)

    ##################################
    # MEAN-MEDIAN DIFFERENCE HISTOGRAM
    ##################################
    fig = pyplot.figure()
    pyplot.title('Mean-Median Difference Histogram')
    pyplot.hist(x=numpy.array(flex.double(mean_map_vals)-flex.double(medn_map_vals)), bins=30, normed=True)
    pyplot.xlabel('Difference Map Value')
    pyplot.ylabel('Density')
    #pyplot.tight_layout()
    pyplot.subplots_adjust(hspace=0.4)
    pyplot.savefig(pandda.file_manager.get_file('map_mean_median_diff').format(resolution))
    pyplot.close(fig)

    ##################################
    # STATISTICAL MAP HISTOGRAMS
    ##################################
    fig = pyplot.figure()
    pyplot.title('Point Variation Map Histograms')
    # STANDARD DEVIATION MAPS
    pyplot.subplot(2, 1, 1)
    pyplot.hist(x=stds_map_vals, bins=n_bins)
    pyplot.xlabel('"Raw" Variation of Map Values')
    pyplot.ylabel('Count')
    # ADJUSTED STANDARD DEVIATION MAPS
    pyplot.subplot(2, 1, 2)
    pyplot.hist(x=sadj_map_vals, bins=n_bins)
    pyplot.xlabel('"Adjusted" Variation of Map Values')
    pyplot.ylabel('Count')
    #pyplot.tight_layout()
    pyplot.subplots_adjust(hspace=0.4)
    pyplot.savefig(pandda.file_manager.get_file('map_stds_sadj_hist').format(resolution))
    pyplot.close(fig)

    ##################################
    # STD Values v ADJ STD Values
    ##################################
    fig = pyplot.figure()
    pyplot.title('Raw v Adjusted Map Variation')
    pyplot.scatter(x=stds_map_vals, y=sadj_map_vals)
    # Plot straight line between the min and max values
    min_val = min(stds_map_vals+sadj_map_vals)
    max_val = max(stds_map_vals+sadj_map_vals)
    pyplot.plot([min_val, max_val], [min_val, max_val], 'b--')
    # Axis labels
    pyplot.xlabel('"Raw" Variation of Map Values')
    pyplot.ylabel('"Adjusted" Variation of Map Values')
    #pyplot.tight_layout()
    pyplot.subplots_adjust(hspace=0.4)
    pyplot.savefig(pandda.file_manager.get_file('map_stds_sadj_scat').format(resolution))
    pyplot.close(fig)

    ########################################################

    pandda.log('=> Writing Map Uncertainties')

    # MAP PARAMS
    fig = pyplot.figure()
    pyplot.title('Map Statistics Histrograms')
    # MAP UNCERTAINTIES
    pyplot.hist(x=map_uncties, bins=n_bins, range=(numpy.nanmin(map_uncties)-0.1,numpy.nanmax(map_uncties)+0.1))
    pyplot.xlabel('Dataset Map Uncertainties')
    pyplot.ylabel('Count')
    #pyplot.tight_layout()
    pyplot.subplots_adjust(hspace=0.4)
    pyplot.savefig(pandda.file_manager.get_file('dataset_unc_hist').format(resolution))
    pyplot.close(fig)

    ########################################################

    pandda.log('=> Scatter Plots')

    # MAP RESOLUTION V UNCERTAINTY
    fig = pyplot.figure()
    pyplot.title('High Resolution Limit v Map Uncertainty')
    pyplot.scatter(x=high_res, y=map_uncties)
    pyplot.xlabel('Dataset Resolution')
    pyplot.ylabel('Dataset Map Uncertainty')
    #pyplot.tight_layout()
    pyplot.subplots_adjust(hspace=0.4)
    pyplot.savefig(pandda.file_manager.get_file('dataset_res_unc_scat').format(resolution))
    pyplot.close(fig)

    # MAP RESOLUTION V UNCERTAINTY
    fig = pyplot.figure()
    pyplot.title('High Resolution Limit v R-Free')
    pyplot.scatter(x=high_res, y=rfree)
    pyplot.xlabel('Dataset Resolution')
    pyplot.ylabel('Dataset R-Free')
    #pyplot.tight_layout()
    pyplot.subplots_adjust(hspace=0.4)
    pyplot.savefig(pandda.file_manager.get_file('dataset_res_rfree_scat').format(resolution))
    pyplot.close(fig)

    # RFREE V UNCERTAINTY
    fig = pyplot.figure()
    pyplot.title('R-Free v Uncertainty')
    pyplot.scatter(x=rfree, y=map_uncties)
    pyplot.xlabel('Dataset R-Free')
    pyplot.ylabel('Dataset Map Uncertainty')
    #pyplot.tight_layout()
    pyplot.subplots_adjust(hspace=0.4)
    pyplot.savefig(pandda.file_manager.get_file('dataset_unc_rfree_scat').format(resolution))
    pyplot.close(fig)

    ########################################################

    pandda.log('=> Z-Map Distribution')

    # R-FACTORS
    fig = pyplot.figure()
    pyplot.title('Z-Map Statistics Histograms')
    # RFree
    pyplot.subplot(2, 1, 1)
    pyplot.hist(x=z_map_mean, bins=n_bins, range=(numpy.nanmin(z_map_mean)-0.1,numpy.nanmax(z_map_mean)+0.1))
    pyplot.xlabel('Z-Map Mean')
    pyplot.ylabel('Count')
    # RWork
    pyplot.subplot(2, 1, 2)
    pyplot.hist(x=z_map_std, bins=n_bins, range=(0, numpy.nanmax(z_map_std)+0.1))
    pyplot.xlabel('Z-Map Std')
    pyplot.ylabel('Count')
    #pyplot.tight_layout()
    pyplot.subplots_adjust(hspace=0.4)
    pyplot.savefig(pandda.file_manager.get_file('zmap_mean_sadj_hist').format(resolution))
    pyplot.close(fig)

    # Z-MAP SKEW V UNCERTAINTY
    fig = pyplot.figure()
    pyplot.title('Z-Map Normality Plots')
    pyplot.scatter(x=z_map_skew, y=z_map_kurt)
    pyplot.xlabel('Skew')
    pyplot.ylabel('Kurtosis')
    #pyplot.tight_layout()
    pyplot.subplots_adjust(hspace=0.4)
    pyplot.savefig(pandda.file_manager.get_file('zmap_skew_kurt_scat').format(resolution))
    pyplot.close(fig)

