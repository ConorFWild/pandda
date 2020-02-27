import os, sys, glob, re
import time

import numpy

from cctbx import maptbx
from scitbx.array_family import flex

from giant.grid import Grid


class AtomicMask(object):

    def __init__(self, parent, sites_cart, max_dist, min_dist=0):
        """Take a grid and calculate all grid points with a certain distance cutoff of any point in sites_cart"""

        if min_dist: assert max_dist > min_dist, 'Minimum Mask Distance must be smaller than Maximum Mask Distance'

        assert isinstance(parent, Grid)
        self.parent = parent
        self._max_dist = max_dist
        self._min_dist = min_dist
        self.sites_cart = sites_cart

        # Get copies of commonly-used grid objects
        self.indexer = parent.indexer()

        # Create the grid mask selections
        self._mask()

    def _mask(self):

        # Calculate the masked indices defined by max distance from protein atoms
        self._outer_mask_indices = maptbx.grid_indices_around_sites(unit_cell=self.parent.unit_cell(),
                                                                    fft_n_real=self.parent.grid_size(),
                                                                    fft_m_real=self.parent.grid_size(),
                                                                    # Masking is performed relative to the grid origin, so need to apply origin shift
                                                                    sites_cart=self.sites_cart - self.parent.cart_origin(),
                                                                    site_radii=flex.double(self.sites_cart.size(),
                                                                                           self._max_dist))
        self._outer_mask_binary = numpy.zeros(self.parent.grid_size_1d(), dtype=bool)
        self._outer_mask_binary.put(self._outer_mask_indices, True)

        # Calculate the masked indices defined by min distance from protein atoms
        self._inner_mask_indices = maptbx.grid_indices_around_sites(unit_cell=self.parent.unit_cell(),
                                                                    fft_n_real=self.parent.grid_size(),
                                                                    fft_m_real=self.parent.grid_size(),
                                                                    # Masking is performed relative to the grid origin, so need to apply origin shift
                                                                    sites_cart=self.sites_cart - self.parent.cart_origin(),
                                                                    site_radii=flex.double(self.sites_cart.size(),
                                                                                           self._min_dist))
        self._inner_mask_binary = numpy.zeros(self.parent.grid_size_1d(), dtype=bool)
        self._inner_mask_binary.put(self._inner_mask_indices, True)

        # Calculate the combination of these masks
        self._total_mask_binary = numpy.zeros(self.parent.grid_size_1d(), bool)
        self._total_mask_binary.put(self._outer_mask_indices, True)
        self._total_mask_binary.put(self._inner_mask_indices, False)
        self._total_mask_indices = numpy.where(self._total_mask_binary)[0]

    def total_mask(self):
        """Return the grid points allowed by the mask - combination of max_dist (allowed) and min_dist (rejected)"""
        for p in flex.nested_loop(self.parent.grid_size()):
            if self._total_mask_binary[self.indexer(p)]: yield p

    def outer_mask(self):
        """Get grid points allowed subject to max_dist"""
        for p in flex.nested_loop(self.parent.grid_size()):
            if self._outer_mask_binary[self.indexer(p)]: yield p

    def inner_mask(self):
        """Get grid points rejected subject to min_dist"""
        for p in flex.nested_loop(self.parent.grid_size()):
            if self._inner_mask_binary[self.indexer(p)]: yield p

    def total_mask_binary(self):
        return self._total_mask_binary

    def outer_mask_binary(self):
        return self._outer_mask_binary

    def inner_mask_binary(self):
        return self._inner_mask_binary

    def total_mask_indices(self):
        return self._total_mask_indices

    def outer_mask_indices(self):
        return self._outer_mask_indices

    def inner_mask_indices(self):
        return self._inner_mask_indices

    def total_size(self):
        """Returns the number of grid points in the mask"""
        return len(self._total_mask_indices)

    def outer_size(self):
        """Returns the number of grid points inside max_dist"""
        return len(self._outer_mask_indices)

    def inner_size(self):
        """Returns the number of grid points inside min_dist"""
        return len(self._inner_mask_indices)

    def extent(self):
        """Returns the minimum and maximum grid points in the mask"""
        return min(self.total_mask()), max(self.total_mask())

    def summary(self):
        return '\n'.join(['----------------------------------->>>',
                          'Atomic Mask Summary:',
                          'Total Mask Size (1D): {!s}'.format(self.total_size()),
                          'Outer Mask Size (1D): {!s}'.format(self.outer_size()),
                          'Inner Mask Size (1D): {!s}'.format(self.inner_size()),
                          'Masked Grid Min/Max: {!s}'.format(self.extent())
                          ])


class GridMask(AtomicMask):

    def __init__(self, parent, sites_cart, max_dist, min_dist=0):
        """Mask a grid (approximately) against a set of cartesian points, without allowing for lattice symmetry (rejects atoms outside of the unit cell)"""

        if min_dist: assert max_dist > min_dist, 'Minimum Mask Distance must be smaller than Maximum Mask Distance'

        # assert isinstance(parent, Grid)
        self.parent = parent
        self._max_dist = max_dist
        self._min_dist = min_dist

        # Get copies of commonly-used grid objects
        self.indexer = parent.indexer()

        # Masking is performed relative to the grid origin, so need to apply origin shift
        sites_cart = sites_cart - self.parent.cart_origin()

        # Filter out sites that are outside of the unit cell
        sites_frac = self.parent.unit_cell().fractionalize(sites_cart)
        ax1_len_1 = (sites_frac.dot((1.0, 0.0, 0.0)) < 1.000001).iselection()
        ax2_len_1 = (sites_frac.dot((0.0, 1.0, 0.0)) < 1.000001).iselection()
        ax3_len_1 = (sites_frac.dot((0.0, 0.0, 1.0)) < 1.000001).iselection()
        ax1_len_2 = (sites_frac.dot((1.0, 0.0, 0.0)) > -0.000001).iselection()
        ax2_len_2 = (sites_frac.dot((0.0, 1.0, 0.0)) > -0.000001).iselection()
        ax3_len_2 = (sites_frac.dot((0.0, 0.0, 1.0)) > -0.000001).iselection()
        allowed_idxs_1 = ax1_len_1.intersection(ax2_len_1).intersection(ax3_len_1)
        allowed_idxs_2 = ax1_len_2.intersection(ax2_len_2).intersection(ax3_len_2)
        allowed_idxs = allowed_idxs_1.intersection(allowed_idxs_2)
        if not allowed_idxs: raise Exception('No Sites valid for masking!')

        # Reverse the origin shift
        self.sites_cart = sites_cart.select(allowed_idxs) + self.parent.cart_origin()

        self._mask()
