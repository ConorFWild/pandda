
set normalize_ccp4_maps, off

# Load the aligned pdb structure
load {{ file_dict['aligned_structure'] }}, aligned_structure

# Load the observed density map
load {{ file_dict['sampled_map'] }}, sampled_map

# Load the difference z map
load {{ file_dict['z_map_corrected_normalised'] }}, z_map

# Load the mask of the large z-values
load {{ file_dict['high_z_mask'] }}, z_map_blobs

color atomic, aligned_structure

isomesh sampled_map_mesh, sampled_map, 1, (aligned_structure), 10
color blue, sampled_map_mesh
hide mesh, sampled_map_mesh

isomesh z_map_mesh_pos, z_map, 3, (aligned_structure), 10
color yellow, z_map_mesh_pos
hide mesh, z_map_mesh_pos

isomesh z_map_mesh_neg, z_map, -3, (aligned_structure), 10
color magenta, z_map_mesh_neg
hide mesh, z_map_mesh_neg

isomesh z_map_blobs_mesh, z_map_blobs, 1, (aligned_structure), 10
color grey, z_map_blobs_mesh

