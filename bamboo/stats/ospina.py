from scipy.optimize import fsolve
import numpy

import time

def variable_sigma_d_log_likelihood(est_sigma, est_mu, obs_vals, obs_error):
    """Calculate the value of the differentiated log likelihood for the values of mu, sigma"""
    term1 = (obs_vals - est_mu)**2 / ((est_sigma**2 + obs_error**2)**2)
    term2 = 1 / (est_sigma**2 + obs_error**2)
    return numpy.sum(term1) - numpy.sum(term2)

def estimate_true_underlying_sd(obs_vals, obs_error, est_mu=None, est_sigma=1e-16, try_number=1):
    """Given a  set of observations `obs_vals` with estimated errors `obs_sigma`, estimate the mean and sd of the underlying distribution"""

    if try_number > 1:
        #print 'ITERATION {!s} - EST SIGMA: {!s}'.format(try_number, est_sigma)
        if try_number > 10:
            raise Exception('TOO MANY ITERATIONS IN OSPINA FUNCTION')

    obs_vals  = numpy.array(obs_vals)
    obs_error = numpy.array(obs_error)

    if not est_mu:
        # Calculate the mean of the sample - this is a good estimator of the true mean
        est_mu = numpy.mean(obs_vals)

    # Estimate the sigma of the underlying distribution
    answer = abs(fsolve(func=variable_sigma_d_log_likelihood, x0=est_sigma, args=(est_mu, obs_vals, obs_error)))[0]

    if answer > 2e40:
        est_sigma = est_sigma/1000
        answer = estimate_true_underlying_sd(obs_vals=obs_vals, obs_error=obs_error, est_mu=est_mu, est_sigma=est_sigma, try_number=try_number+1)

    return answer

if __name__ == '__main__':

    # True values we are trying to estimate
    true_mean = 1
    true_sd   = 0.11

    # Number of observations
    num_obs = 200

    # Number of times to run
    num_cycles = 10000

    guesses  = []

    start_t = time.time()

    for attempt in xrange(num_cycles):

#        print('============================>')

        # Sample the original distribution
        true_vals = true_mean + true_sd*numpy.random.randn(num_obs)
#        print('MEAN OF TRUE VALS: {!s} ({!s})'.format(round(numpy.mean(true_vals),3), true_mean))
#        print(' STD OF TRUE VALS: {!s} ({!s})'.format(round(numpy.std(true_vals),3), true_sd))

        # Create a random selection of sigmas for the different observations
        obs_error = numpy.abs(0.45 + 0.01*numpy.random.randn(num_obs))
#        print('MEAN OF OBS ERROR: {!s} ({!s})'.format(round(numpy.mean(obs_error),3), '2ish'))
#        print(' STD OF OBS ERROR: {!s} ({!s})'.format(round(numpy.std(obs_error),3), '1ish'))

        # Noise to be added to the true observations
        obs_noise = numpy.random.randn(num_obs)
#        print('MEAN OF OBS NOISE: {!s} ({!s})'.format(round(numpy.mean(obs_noise),3), 0))
#        print(' STD OF OBS NOISE: {!s} ({!s})'.format(round(numpy.std(obs_noise),3), 1))

        # Create fake data!
        obs_vals = true_vals + obs_error * obs_noise
#        print(' MEAN OF OBS VALS: {!s}'.format(round(numpy.mean(obs_vals),3)))
#        print('  STD OF OBS VALS: {!s}'.format(round(numpy.std(obs_vals),3)))

        out = estimate_true_underlying_sd(obs_vals=obs_vals, obs_error=obs_error, est_sigma=0.1)

        print(' ESTIMATED VALUES: {!s}'.format(round(out,5)))

        guesses.append(out)

    print('============================>')

    end_t = time.time()
    print('    TIME TAKEN TOTAL: {!s} (Seconds)'.format(int(end_t - start_t)))
    print('TIME TAKEN PER CYCLE: {!s} (Millionths)'.format(int(1000000*(end_t-start_t)/num_cycles)))

    print('============================>')

    print('MEAN OF ESTIMATES: {!s}'.format(round(numpy.mean(guesses),5)))
    print(' STD OF ESTIMATES: {!s}'.format(round(numpy.std(guesses),5)))
