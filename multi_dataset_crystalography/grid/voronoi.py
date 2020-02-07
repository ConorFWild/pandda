from __future__ import division

import numpy
import scipy.spatial

def get_voronoi_simplices(points, bounding_box=((0,0,0),(1,1,1))):
    """Get the corners of the voronoi cells"""
    vertices, links = get_voronoi_cell_borders(points)
    cells = assign_divider_lines_to_voronoi_cells(points, vertices, links)
    out_cells = []
    # Iterate through the voronoi cell (list of border lines)
    for i_cell, link_idxs in enumerate(cells):
        out_vertices = []
        # Iterate through the borders of the cell (list of vertex indices)
        for i_link, (i_v1,i_v2) in enumerate(link_idxs):
            # Normalise the bounding line to the bounding box
            v1, v2 = truncate_lines_to_bounding_box((vertices[i_v1],vertices[i_v2]))
            out_vertices.append((v1,v2))
        out_vertices = numpy.array(out_vertices)
    out_cells.append(out_vertices)
    return out_cells

def get_voronoi_cell_borders(points):
    """Create voronoi cells. Returns cell vertices and links between vertices"""
    # Create delaunay object for the points
    delauny = scipy.spatial.Delaunay(points)
    # Get the simplices from the list of vertex indices
    simplices = delauny.points[delauny.vertices]
    # Calculate the circumcircle centres of the simplices
    circum_centres = numpy.array(map(simplex_circumcircle, simplices))
    # Endpoints outside of the convex hull of the points
    extra_endpoints = []
    # Indices of the ends of the lines
    links = []
    # Enumerate through the simplices
    for i, simplex in enumerate(simplices):
        # Extract the circumcentre for this simplex
        circum_centre = circum_centres[i]
        # Enumerate through the neighbours of this simplex
        for j, neighbour in enumerate(delauny.neighbors[i]):
            if neighbour != -1:
                # If it's a valid neighbour, link these two circumcentres
                links.append((i, neighbour))
            else:
                # Make line on the edge of the simplex
                ps1 = simplex[(j+3)%4] - simplex[(j+2)%4]
                ps2 = simplex[(j+2)%4] - simplex[(j+1)%4]
                # Crossproduct with a vector in the
                ps = numpy.cross(ps1, ps2)
                # Find the middle of the edge of the simplex
                middle = (simplex[(j+3)%3] + simplex[(j+2)%3] + simplex[(j+1)%4]) * 0.5
                di = middle - simplex[j]
                ps /= numpy.linalg.norm(ps)
                di /= numpy.linalg.norm(di)
                # Extend the length of the vector
                if numpy.dot(di, ps) < 0.0:
                    ps *= -10.0
                else:
                    ps *= 10.0
                # Append to the list of endpoints
                extra_endpoints.append(circum_centre + ps)
                # Add the index where this point will be after concatenation of the two lists
                links.append((i, len(circum_centres) + len(extra_endpoints)-1))
    # Concatenate all of the points
    vertices = numpy.vstack((circum_centres, extra_endpoints))
    # Sort, tupleise, and make unique
    links = numpy.sort(links) # make (1,2) and (2,1) both (1,2)
    links = sorted(set([tuple(row) for row in links]))
    return vertices, links

def simplex_circumcircle(pts):
    """Create circumcircle from the triangle of points"""
    rows, cols = pts.shape
    A = numpy.bmat([[2 * numpy.dot(pts, pts.T), numpy.ones((rows, 1))], [numpy.ones((1, rows)), numpy.zeros((1, 1))]])
    b = numpy.hstack((numpy.sum(pts * pts, axis=1), numpy.ones((1))))
    x = numpy.linalg.solve(A,b)
    bary_coords = x[:-1]
    return numpy.sum(pts * numpy.tile(bary_coords.reshape((pts.shape[0], 1)), (1, pts.shape[1])), axis=0)

def assign_divider_lines_to_voronoi_cells(points, vertices, links):
    """Assign lines to the nearest points by the distance between the line midpoint and the point"""
    kd = scipy.spatial.KDTree(points)
    cells = [[]]*len(points)
    for i1,i2 in links:
        v1,v2 = vertices[i1], vertices[i2]
        mid = (v1+v2)/2
        _, (p1Idx,p2Idx) = kd.query(mid, 2)
        cells[p1Idx].append((i1,i2))
        cells[p2Idx].append((i1,i2))
    return cells

def truncate_lines_to_bounding_box(end_points, box_min=(0,0,0), box_max=(1,1,1)):
    """Given two points, move endpoint to be within bounding box keeping the joinging vector"""
    box_min = numpy.array(box_min)
    box_max = numpy.array(box_max)
    xyz1, xyz2 = numpy.array(end_points)
    vec = (xyz2-xyz1)
    # Check which way the vector points and reverse if necessary
    a_vec = numpy.abs(vec)
    s_vec = numpy.sign(vec)
    b_vec = (vec.dot([1,1,1]) > 0)

    # Calculate how much of each vector needs to be added or subtracted
    d_11 = (box_min-xyz1 > 0).astype(float) * (box_min-xyz1)
    d_21 = (box_min-xyz2 > 0).astype(float) * (box_min-xyz2)
    d_12 = (xyz1-box_max > 0).astype(float) * (xyz1-box_max)
    d_22 = (xyz2-box_max > 0).astype(float) * (xyz2-box_max)
    i_11 = numpy.where(d_11==max(d_11))[0][0]
    i_21 = numpy.where(d_21==max(d_21))[0][0]
    i_12 = numpy.where(d_12==max(d_12))[0][0]
    i_22 = numpy.where(d_22==max(d_22))[0][0]
    v_11 = d_11[i_11]*vec/vec[i_11]
    v_21 = d_21[i_21]*vec/vec[i_21]
    v_12 = d_12[i_12]*vec/vec[i_12]
    v_22 = d_22[i_22]*vec/vec[i_22]
    # Create output vectors
    new1 = xyz1-(1-b_vec)*v_12+b_vec*v_11
    new2 = xyz2-b_vec*v_22+(1-b_vec)*v_21
    print vec,vec
    print '-------------------'
    print d_11, d_21
    print d_12, d_22
    print '-------------------'
    print xyz1,xyz2
    print new1,new2
    print '==================='
    # Check the output is within the box
    assert sum((box_max-new1)<-1e-3) == 0
    assert sum((box_max-new2)<-1e-3) == 0
    assert sum((new1-box_min)<-1e-3) == 0
    assert sum((new2-box_min)<-1e-3) == 0
    return new1, new2

