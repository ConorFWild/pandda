import numpy

from bamboo.common import Report

from scitbx.math import scale_curves, approx_equal_relatively
from scitbx.array_family import flex
from scitbx.python_utils.robust_statistics import percentile
from mmtbx.scaling import absolute_scaling

from giant.stats.optimisation import ExponentialScaling

class IsotropicBfactorScalingFactory(object):
    """Adapted from mmtbx.scaling.relative_wilson"""

    def __init__(self, reference_miller_array, d_star_min_sq=0.0625, d_star_max_sq=2.0, n_scale=2000):
        """
        Calculate scaling between any miller array and a reference miller array.
        This does not assume that the two arrays come from the same crystal.
        Data will only be scaled over the common resolution range of the crystal.
        For best results, the two miller arrays should have high completeness
        (although this is not checked).
        """

        # Number of points for scaling
        self._npoints = n_scale

        # Convert to intensities and store
        self.ref_miller = reference_miller_array.as_intensity_array()
        self.ref_kernel = self._kernel_normalisation(miller_array=self.ref_miller)

        # Store d_star limits and then update due to reference limits
        self.d_star_min_sq, self.d_star_max_sq = d_star_min_sq, d_star_max_sq
        self.d_star_min_sq, self.d_star_max_sq = self._common_d_star_sq_range(d_star_sq=self.ref_kernel.d_star_sq_array)

        return None

    def _common_d_star_sq_range(self, d_star_sq):
        """Find the common range over current limits and those of supplied array"""

        d_star_min_sq = max(self.d_star_min_sq, d_star_sq.min_max_mean().min)
        d_star_max_sq = min(self.d_star_max_sq, d_star_sq.min_max_mean().max)
        return d_star_min_sq, d_star_max_sq

    def _kernel_normalisation(self, miller_array, n_bins=45, n_term=17):

        return absolute_scaling.kernel_normalisation(miller_array=miller_array,
                                                     auto_kernel=True,
                                                     n_bins=n_bins,
                                                     n_term=n_term)

    def calculate_scaling(self, miller_array, convergence_crit_perc=0.01, convergence_reject_perc=97.5, max_iter=20):
        """Calculate the scaling between two arrays"""

        assert convergence_reject_perc > 90.0

        # Convert to intensities and extract d_star_sq
        new_miller = miller_array.as_intensity_array()
        new_kernel = self._kernel_normalisation(miller_array=new_miller)
        # Calculate new range of d_star_sq
        d_star_sq_min, d_star_sq_max = self._common_d_star_sq_range(d_star_sq=new_kernel.d_star_sq_array)

        # Create interpolator for the two arrays (new and reference)
        interpolator = scale_curves.curve_interpolator(d_star_sq_min, d_star_sq_max, self._npoints)

        # Interpolate the two curves (use full range of the two array)
        new_itpl_d_star_sq, new_itpl_mean_I, dummy, dummy = interpolator.interpolate(
                                                                x_array=new_kernel.d_star_sq_array,
                                                                y_array=new_kernel.mean_I_array)
        ref_itpl_d_star_sq, ref_itpl_mean_I, dummy, dummy = interpolator.interpolate(
                                                                x_array=self.ref_kernel.d_star_sq_array,
                                                                y_array=self.ref_kernel.mean_I_array)

        # Initalise convergence loop - begin by scaling over all points
        selection = flex.bool(self._npoints, True)
        # Set initial scale factor to small value
        curr_b = 1e-6
        # Percent change between iterations - convergence when delta <convergence_criterion
        n_iter = 0
        # Report in case of error
        report = Report('Scaling log:', verbose=False)
        while n_iter < max_iter:
            report('---')
            report('ITER: '+str(n_iter))

            if selection.all_eq(False):
                print("Selection now empty, breaking")
                break

            # Run optimisation on the linear scaling
            lsc = ExponentialScaling(x_values = interpolator.target_x,
                                     ref_values = ref_itpl_mean_I,
                                     scl_values = new_itpl_mean_I,
                                     weights = selection.as_double())
            # Calculate scaling B-factor
            lsc.scaling_b_factor = -0.5 * list(lsc.optimised_values)[0]
            # Break if fitted to 0
            if approx_equal_relatively(0.0, lsc.scaling_b_factor, 1e-6):
                report('Scaling is approximately 0.0 - stopping')
                break
            # Calculate percentage change
            report('Curr/New: '+str(curr_b)+'\t'+str(lsc.scaling_b_factor))
            delta = abs((curr_b-lsc.scaling_b_factor)/curr_b)
            report('Delta: '+str(delta))
            if delta < convergence_crit_perc:
                report('Scaling has converged to within tolerance - stopping')
                break
            # Update selection
            report('Curr Selection Size: '+str(sum(selection)))
            ref_diffs = flex.log(lsc.ref_values)-flex.log(lsc.out_values)
            #abs_diffs = flex.abs(ref_diffs)
            sel_diffs = ref_diffs.select(selection)
            rej_val_high = numpy.percentile(sel_diffs, convergence_reject_perc)
            rej_val_low  = numpy.percentile(sel_diffs, 100.0-convergence_reject_perc)
            report('Percentile: '+str(convergence_reject_perc)+'\t<'+str(rej_val_low)+'\t>'+str(rej_val_high))
            selection.set_selected(ref_diffs>rej_val_high, False)
            selection.set_selected(ref_diffs<rej_val_low,  False)

            report('New Selection Size: '+str(sum(selection)))
            # Update loop params
            curr_b = lsc.scaling_b_factor
            n_iter += 1


        lsc.unscaled_ln_rmsd = (flex.log(lsc.ref_values)-flex.log(lsc.scl_values)).norm()/(lsc.ref_values.size()**0.5)
        lsc.scaled_ln_rmsd   = (flex.log(lsc.ref_values)-flex.log(lsc.out_values)).norm()/(lsc.ref_values.size()**0.5)

        lsc.unscaled_ln_dev = flex.sum(flex.abs(flex.log(lsc.ref_values)-flex.log(lsc.scl_values)))
        lsc.scaled_ln_dev   = flex.sum(flex.abs(flex.log(lsc.ref_values)-flex.log(lsc.out_values)))

        return lsc

