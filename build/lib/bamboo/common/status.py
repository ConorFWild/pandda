from __future__ import print_function

import os, sys, glob

CLEAR_LINE = '\x1b[2K'
START_LINE = '\x1B[0E'

def status_bar(n, n_max):
    n_jump = max(1,int(n_max/100))
    assert n<=n_max, 'n must be less than or equal to n_max'
    if (n==n_max):
        print('\r>> 100 %')
    elif (n%n_jump==0):
        print('\r>> {:3d} %'.format(int(round(100.0*n/n_max,0))), end='')
    sys.stdout.flush()

def status_bar_2(n, n_max, width=50):
    assert n<=n_max, 'n must be less than or equal to n_max'
    if (n==n_max):
        print('\r>> 100 % [{}]'.format('='*width))
        return
    # Get a spacing so that at most 100 updates are performed
    n_jump = max(1,int(n_max/100.0))
    # Check if we're on this spacing
    if (n%n_jump==0):
        frac = (1.0*n)/n_max
        num_bars = int(width*frac)
        bars   = '='*num_bars
        blanks = ' '*(width-num_bars-1)
        perc = int(100.0*frac)
        print('\r>> {:3d} % [{}>{}]'.format(perc, bars, blanks), end='')
        sys.stdout.flush()

