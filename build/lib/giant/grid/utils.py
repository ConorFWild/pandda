import numpy
from scitbx.array_family import flex
from libtbx.math_utils import ifloor, iceil

def idx_to_grid(idx, grid_size):
    """Return the grid point for the index"""
    gp = []
    tot = idx
    for i,j in enumerate(grid_size):
        # Whats the multiplier on this row
        mult = int(numpy.prod(grid_size[i+1:]))
        # How many of this can we get
        val=tot//mult
        tot-=val*mult
        gp.append(val)
    return tuple(gp)

def calculate_grid_size(min_carts, max_carts, grid_spacing):
    """Calculate the number of points to be sampled for a box size and sampling distance. Returns the number of points to be sampled along each axis. Box may be larger than max_carts."""
    cart_size = tuple([float(max_c-min_c) for min_c, max_c in zip(min_carts, max_carts)])
    grid_size = tuple([iceil(c_size/(1.0*grid_spacing))+1 for c_size in cart_size])
    return grid_size

def get_grid_points_within_distance_cutoff_of_origin(grid_spacing, distance_cutoff):
    """Find all points on isotropic grid within distance_cutoff of the origin"""
    grid_index_cutoff = (1.0*distance_cutoff)/grid_spacing
    return get_grid_points_within_index_cutoff_of_origin(grid_index_cutoff=grid_index_cutoff)

def get_grid_points_within_index_cutoff_of_origin(grid_index_cutoff):
    """Find all points relative to the origin within a certain number of grid points"""

    # Round the max grid index down to the nearest int - outer bound on the x,y,z coords
    outer_bound_box = int(grid_index_cutoff)
    # Calculate r/sqrt(3) - inner bound on x,y,z coords
    inner_bound_box = grid_index_cutoff/numpy.math.sqrt(3)
    # Calculate r^2 - limiting sphere
    rad_sq = grid_index_cutoff**2
    # List of allowed grid indices
    grid_points = []

    for x,y,z in flex.nested_loop(begin=(-outer_bound_box, -outer_bound_box, -outer_bound_box),end=(outer_bound_box+1, outer_bound_box+1, outer_bound_box+1)):
        if (abs(x) <= inner_bound_box) and (abs(y) <= inner_bound_box) and (abs(z) <= inner_bound_box):
            grid_points.append((x, y, z))
        elif (x**2 + y**2 + z**2) <= rad_sq:
            grid_points.append((x, y, z))

    return grid_points

def combine_grid_point_and_grid_vectors(start_point, grid_vectors):
    """Take a list of grid vectors and add them to a particular grid point"""

    start_point = flex.int(start_point)
    grid_points = [tuple(start_point+flex.int(vec)) for vec in grid_vectors]
    return grid_points

def get_grid_points_within_distance_cutoff_of_cart_sites(cart_sites, grid_spacing, max_dist, min_dist=None):
    """Find all points on isotropic grid within distance cutoff of cartesian sites"""

    # Normalise the cutoff and the sites to make them grid-independent
    max_grid_dist = (1.0*max_dist)/grid_spacing
    if min_dist: min_grid_dist = (1.0*min_dist)/grid_spacing
    else:        min_grid_dist = None
    grid_sites = [tuple([1.0*c/grid_spacing for c in coords]) for coords in cart_sites]

    return get_grid_points_within_index_cutoff_of_grid_sites(grid_sites=grid_sites, max_grid_dist=max_grid_dist, min_grid_dist=min_grid_dist)

def get_grid_points_within_distance_cutoff_of_cart_sites_2(cart_sites, grid_spacing, max_dist, min_dist=None):
    """Find all points on isotropic grid within distance cutoff of cartesian sites"""

    # Normalise the cutoff and the sites to make them grid-independent
    max_grid_dist = (1.0*max_dist)/grid_spacing
    if min_dist: min_grid_dist = (1.0*min_dist)/grid_spacing
    else:        min_grid_dist = None
    grid_sites = [tuple([1.0*c/grid_spacing for c in coords]) for coords in cart_sites]

    return get_grid_points_within_index_cutoff_of_grid_sites_2(grid_sites=grid_sites, max_grid_dist=max_grid_dist, min_grid_dist=min_grid_dist)

def get_grid_points_within_index_cutoff_of_grid_sites(grid_sites, max_grid_dist, min_grid_dist=None):
    """Find all points on a grid within a certain number of grid points of grid sites (not necessarily integer sites)"""

    # Calculate the size of the grid we need to check over
    min_x = ifloor(min([s[0] for s in grid_sites]) - max_grid_dist)
    if min_x < 0: min_x=0
    min_y = ifloor(min([s[1] for s in grid_sites]) - max_grid_dist)
    if min_y < 0: min_y=0
    min_z = ifloor(min([s[2] for s in grid_sites]) - max_grid_dist)
    if min_z < 0: min_z=0
    max_x = iceil(max([s[0] for s in grid_sites]) + max_grid_dist)
    max_y = iceil(max([s[1] for s in grid_sites]) + max_grid_dist)
    max_z = iceil(max([s[2] for s in grid_sites]) + max_grid_dist)
    # Grid Extremities
    min_grid = (min_x, min_y, min_z)
    max_grid = (max_x+1, max_y+1, max_z+1)

    # Round the max grid distance down to the nearest int - outer bound on the dx,dy,dz values
    outer_bound_box = int(max_grid_dist)
    # Calculate r/sqrt(3) - inner bound on dx, dy, dz values
    inner_bound_box = max_grid_dist/numpy.math.sqrt(3)
    # Calculate r^2 - limiting sphere
    rad_sq = max_grid_dist**2

    # List of allowed grid indices
    outer_points = []
    inner_points = []

    # Iterate through and add valid points
    for gp in flex.nested_loop(min_grid, max_grid):
        for site in grid_sites:
            dx, dy, dz = [abs(p1-p2) for p1,p2 in zip(gp, site)]

            if (dx > outer_bound_box) or (dy > outer_bound_box) or (dz > outer_bound_box):
                continue
            elif (dx <= inner_bound_box) and (dy <= inner_bound_box) and (dz <= inner_bound_box):
                outer_points.append(gp)
                break
            elif (dx**2 + dy**2 + dz**2) <= rad_sq:
                outer_points.append(gp)
                break

    # Filter the grid points that are too close to the protein
    if min_grid_dist:
        # Round the min grid distance up to the nearest int - outer bound on the dx,dy,dz values
        outer_bound_box = int(min_grid_dist) + 1
        # Calculate r/sqrt(3) - inner bound on dx, dy, dz values
        inner_bound_box = min_grid_dist/numpy.math.sqrt(3)
        # Calculate r^2 - limiting sphere
        rad_sq = min_grid_dist**2

        # Iterate through and add valid points
        for gp in outer_points:
            for site in grid_sites:
                dx, dy, dz = [abs(p1-p2) for p1,p2 in zip(gp, site)]

                if (dx > outer_bound_box) or (dy > outer_bound_box) or (dz > outer_bound_box):
                    continue
                elif (dx <= inner_bound_box) and (dy <= inner_bound_box) and (dz <= inner_bound_box):
                    inner_points.append(gp)
                    break
                elif (dx**2 + dy**2 + dz**2) <= rad_sq:
                    inner_points.append(gp)
                    break

    if min_grid_dist:
        total_point = [gp for gp in outer_points if gp not in inner_points]
    else:
        total_points = outer_points

    return total_points, outer_points, inner_points

def get_grid_points_within_index_cutoff_of_grid_sites_2(grid_sites, max_grid_dist, min_grid_dist=None):
    """Find all points on a grid within a certain number of grid points of grid sites (not necessarily integer sites)"""

    # Find the size of the grid that we'll need
    max_x = iceil(max([s[0] for s in grid_sites]) + max_grid_dist)
    max_y = iceil(max([s[1] for s in grid_sites]) + max_grid_dist)
    max_z = iceil(max([s[2] for s in grid_sites]) + max_grid_dist)
    max_grid = (max_x+2, max_y+2, max_z+2)

    # Round the grid sites to the nearest grid point
    int_grid_sites = [tuple([int(round(x,0)) for x in site]) for site in grid_sites]

    # Grid objects
    grid_indexer = flex.grid(max_grid)
    grid_size = flex.product(flex.int(max_grid))

    # Outer mask
    outer_indices_mask = numpy.zeros(grid_size, dtype=int)

    # Find all of the grid vectors within max_dist of grid sites
    outer_grid_vectors = get_grid_points_within_index_cutoff_of_origin(grid_index_cutoff=max_grid_dist)
    outer_grid_points = []
    for site in int_grid_sites:
        outer_grid_points.extend(combine_grid_point_and_grid_vectors(site, grid_vectors=outer_grid_vectors))

    # They may overlap, so find unique grid_points
    outer_grid_points = list(set(outer_grid_points))
    # Filter those that are outside of the grid
    outer_grid_points = [gp for gp in outer_grid_points if not ( ([i+1 for i in range(3) if gp[i]<0]) or ([i+1 for i in range(3) if gp[i]>max_grid[i]]) )]
    # Map the grid points to the associated 1d index
    outer_grid_indices = [grid_indexer(gp) for gp in outer_grid_points]
    # Create a binary mask of the points
    [outer_indices_mask.put(i, 1) for i in outer_grid_indices]

    # Inner mask
    inner_indices_mask = numpy.zeros(grid_size, dtype=int)
    if min_grid_dist:
        # Find all of the grid vectors within min_dist of grid sites
        inner_grid_vectors = get_grid_points_within_index_cutoff_of_origin(grid_index_cutoff=min_grid_dist)
        inner_grid_points = []
        for site in int_grid_sites:
            inner_grid_points.extend(combine_grid_point_and_grid_vectors(site, grid_vectors=inner_grid_vectors))

        # They may overlap, so find unique grid_points
        inner_grid_points = list(set(inner_grid_points))
        # Filter those that are outside of the grid
        inner_grid_points = [gp for gp in inner_grid_points if not ( ([i+1 for i in range(3) if gp[i]<0]) or ([i+1 for i in range(3) if gp[i]>=max_grid[i]]) )]
        # Map the grid points to the associated 1d index
        inner_grid_indices = [grid_indexer(gp) for gp in inner_grid_points]
        # Create a binary mask of the points
        [inner_indices_mask.put(i, 1) for i in inner_grid_indices]

    # Convert from the mask back to grid points
    outer_points = [gp for i, gp in enumerate(flex.nested_loop(max_grid)) if outer_indices_mask[i]==1]
    inner_points = [gp for i, gp in enumerate(flex.nested_loop(max_grid)) if inner_indices_mask[i]==1]
    total_points = [gp for i, gp in enumerate(flex.nested_loop(max_grid)) if (inner_indices_mask[i]==0 and outer_indices_mask[i]==1)]

    return total_points, outer_points, inner_points



