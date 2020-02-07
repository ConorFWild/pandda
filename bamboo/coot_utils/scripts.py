"""Contains strings used to build up a coot script"""

COOT_exit = """coot_real_exit(1)"""

COOT_load_pdb = """read_pdb('<pdbin>')"""
COOT_load_mtz = """auto_read_make_and_draw_maps('<mtzin>')"""

COOT_real_space_refine = """regularize_zone(<imol>, <chain_id>, <resno1>, <resno2>, <altconf>)"""
