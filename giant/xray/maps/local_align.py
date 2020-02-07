import copy, time

import cctbx.uctbx, cctbx.sgtbx, cctbx.maptbx
import scitbx.sparse

from libtbx.math_utils import iceil
from scitbx.array_family import flex

from giant.structure.align import *
from giant.xray.symmetry import get_crystal_contact_operators

def make_supercell(unit_cell, size=(3,3,3)):
    """Create an enlarged supercell composed of multiple unit cells"""
    assert isinstance(size, int) or (len(size) == 3), 'Size must be either single scale factor or 3 scale factors for (x,y,z)'
    if isinstance(size, int):   scales = [size]*3   + [1]*3
    else:                       scales = list(size) + [1]*3
    scaled = map(float, scales)
    assert len(scales) == 6, 'Ooops, something seems to have gone wrong...: size {}, scales {}'.format(size, scales)
    old_params = unit_cell.parameters()
    new_params = [a*b for a,b in zip(old_params, scales)]
    return cctbx.uctbx.unit_cell(new_params)

def calculate_offset_to_centre_grid(grid_centre, centre_on):
    return [a-b for a,b in zip(centre_on, grid_centre)]

def get_subset_of_grid_points(gridding, grid_indices):
    """Use a set of indices to mask the grid - returns masked grid points"""
    import numpy
    mask_binary = numpy.zeros(gridding.n_grid_points(), dtype=bool)
    mask_binary.put(grid_indices, True)
    grid = flex.grid(gridding.n_real())
    for i, p in enumerate(flex.nested_loop(gridding.n_real())):
        assert i == grid(p)
        if mask_binary[i]: yield p

def create_native_map(native_crystal_symmetry, native_sites, alignment, reference_map, site_mask_radius=6.0, step=0.5, filename=None, verbose=False):
    """
    Transform the reference-aligned map back to the native crystallographic frame
    native_sites            - defines region that map will be masked around
    native_crystal_symmetry - crystal symmetry
    reference_map           - basic_map object of the map in the reference frame
    alignment               - Alignment object used to map between the reference frame and the native frame
    site_mask_radius        - Define mask radius around native_sites
    step                    - grid sampling step
    """

    start=time.time()

    # Output unit cell and spacegroup
    native_unit_cell = native_crystal_symmetry.unit_cell()
    native_space_group = native_crystal_symmetry.space_group()

    # ===============================================================================>>>
    # Create a supercell containing the native_sites at the centre
    # ===============================================================================>>>
    # How many unit cells to include in the supercell
    box_frac_min_max = native_unit_cell.box_frac_around_sites(sites_cart=native_sites, buffer=site_mask_radius+1.0)
    supercell_size = tuple([iceil(ma-mi) for mi, ma in zip(*box_frac_min_max)])
    # create supercell in the native frame
    supercell = make_supercell(native_unit_cell, size=supercell_size)
    # Create gridding for the real unit cell based on fixed step size
    tmp_gridding = cctbx.maptbx.crystal_gridding(unit_cell=native_unit_cell, step=step)
    # Extract the n_real from this gridding to be applied to the real griddings
    uc_n_real = tmp_gridding.n_real()
    sc_n_real = tuple(flex.int(tmp_gridding.n_real())*flex.int(supercell_size))
    if verbose:
        print('Generating supercell of size {}. Unit cell grid size: {}. Supercell grid size: {}.'.format(supercell_size, uc_n_real, sc_n_real))
    # Get griddings based on the n_real determined above
    uc_gridding = cctbx.maptbx.crystal_gridding(unit_cell=native_unit_cell, pre_determined_n_real=uc_n_real)
    sc_gridding = cctbx.maptbx.crystal_gridding(unit_cell=supercell,        pre_determined_n_real=sc_n_real)

    # ===============================================================================>>>
    # Mask the supercell grid around the protein model
    # ===============================================================================>>>
    # calculate the origin of the supercell (centred on the protein model) - adding origin translates "grid frame" to "crystallographic frame"
    model_centroid = tuple((flex.double(native_sites.max()) + flex.double(native_sites.min()))/2.0)
    origin = calculate_offset_to_centre_grid(grid_centre=supercell.orthogonalize((0.5,0.5,0.5)), centre_on=model_centroid)
    # sample the map points near to the protein (transform the structure to be centre of the grid)
    masked_points_indices = cctbx.maptbx.grid_indices_around_sites(unit_cell  = supercell,
                                                                   fft_n_real = sc_gridding.n_real(),
                                                                   fft_m_real = sc_gridding.n_real(),
                                                                   sites_cart = native_sites - origin,
                                                                   site_radii = flex.double(native_sites.size(),site_mask_radius))
    # Create iterator over these points
    masked_points_grid_iter = get_subset_of_grid_points(gridding=sc_gridding, grid_indices=masked_points_indices)
    # Convert grid points to cartesian points
    g2c = cctbx.maptbx.grid2cart(sc_gridding.n_real(), supercell.orthogonalization_matrix())
    masked_points_cart = flex.vec3_double(map(g2c, masked_points_grid_iter)) + origin

#    from bamboo.pymol_utils.shapes import Sphere
#    points = ['from pymol import cmd','from pymol.cgo import *']
#    for i,p in enumerate(masked_points_cart):
#        if i%100==0: points.append(Sphere(p, 0.2).as_cmd('steve'))
#    points = '\n'.join(points)
#    with open(filename+'.pml.py', 'w') as fh: fh.write(points)

    # ===============================================================================>>>
    # Sample masked points from the reference-aligned map
    # ===============================================================================>>>
    # Transform points to the reference frame
    masked_points_transformed = alignment.nat2ref(masked_points_cart)
    # Sample the map at these points
    masked_values = reference_map.get_cart_values(masked_points_transformed)

    # ===============================================================================>>>
    # Create a native-aligned map in the supercell
    # ===============================================================================>>>
    # Create a map of the density
    sc_map_data = numpy.zeros(sc_gridding.n_grid_points(), dtype=numpy.float64)
    sc_map_data.put(masked_points_indices, masked_values)
    sc_map_data = flex.double(sc_map_data)
    sc_map_data.reshape(flex.grid(sc_gridding.n_real()))
    # Transform the points back to the native frame (simple origin shift)
    sc_map_data = cctbx.maptbx.rotate_translate_map(unit_cell          = supercell,
                                                    map_data           = sc_map_data,
                                                    rotation_matrix    = scitbx.matrix.rec([1,0,0,0,1,0,0,0,1], (3,3)).elems,
                                                    translation_vector = (-1.0*scitbx.matrix.rec(origin, (3,1))).elems    )

    # ===============================================================================>>>
    # Apply translations to populate all unit cells of supercell
    # ===============================================================================>>>
    # Create a copy to contain the combined map data
    combined_sc_map_data = copy.deepcopy(sc_map_data)
    # Apply unit cell translation operators
    for x,y,z in flex.nested_loop(supercell_size):
        if x==y==z==0: continue
        # Calculate the translation vector
        unit_cell_shift = native_unit_cell.orthogonalize((x,y,z))
        # Tranform the map
        rt_map_data = cctbx.maptbx.rotate_translate_map(unit_cell          = supercell,
                                                        map_data           = sc_map_data,
                                                        rotation_matrix    = scitbx.matrix.rec([1,0,0,0,1,0,0,0,1], (3,3)).elems,
                                                        translation_vector = scitbx.matrix.rec(unit_cell_shift, (3,1)).elems    )
        # Set any values that are filled in combined_sc_map_data to 0
        rt_map_data.set_selected(flex.abs(combined_sc_map_data) > 1e-6, 0.0)
        # Add values to combined_sc_map_data
        combined_sc_map_data = combined_sc_map_data + rt_map_data

    # ===============================================================================>>>
    # Select the first (on origin) unit cell of the supercell
    # ===============================================================================>>>
    # Get the indices for the first unit cell
    supercell_grid = flex.grid(sc_gridding.n_real())
    supercell_mask = map(supercell_grid, flex.nested_loop(uc_gridding.n_real()))
    # Extract the map data for those values (and reshape to the right size of the unit cell)
    uc_map_data = combined_sc_map_data.select(supercell_mask)
    uc_map_data.reshape(flex.grid(uc_gridding.n_real()))

    # ===============================================================================>>>
    # Apply symmetry operations to generate whole unit cell
    # ===============================================================================>>>
    # Create a copy to contain the combined map data
    combined_uc_map_data = copy.deepcopy(uc_map_data)
    # Apply all symmetry operations to unit cell data
    for sym_op in native_space_group.all_ops():
        if sym_op.as_xyz() == 'x,y,z': continue
        # Get the transformation matrix
        rt_mx = sym_op.as_rational().as_float()
        # Tranform the map
        rt_map_data = cctbx.maptbx.rotate_translate_map(unit_cell          = native_unit_cell,
                                                        map_data           = combined_uc_map_data,
                                                        rotation_matrix    = rt_mx.r.elems,
                                                        #translation_vector = native_unit_cell.orthogonalize((-1.0*rt_mx.t).elems)   )
                                                        translation_vector = native_unit_cell.orthogonalize(rt_mx.t.elems)   )
        # Set any values that are filled in combined_uc_map_data to 0
        rt_map_data.set_selected(flex.abs(combined_uc_map_data) > 1e-6, 0)
        # Add values to combined_uc_map_data
        combined_uc_map_data = combined_uc_map_data + rt_map_data

    # ===============================================================================>>>
    # Write output maps
    # ===============================================================================>>>

    if filename is not None:
        iotbx.ccp4_map.write_ccp4_map(  file_name   = filename,
                                        unit_cell   = native_unit_cell,
                                        space_group = native_space_group,
                                        map_data    = combined_uc_map_data,
                                        labels      = flex.std_string(['Map from pandda'])     )
#        iotbx.ccp4_map.write_ccp4_map(  file_name   = filename.replace('.ccp4','.supercell.ccp4'),
#                                        unit_cell   = supercell,
#                                        space_group = cctbx.sgtbx.space_group('P1'),
#                                        map_data    = sc_map_data,
#                                        labels      = flex.std_string(['Map from pandda'])     )

    return combined_uc_map_data

