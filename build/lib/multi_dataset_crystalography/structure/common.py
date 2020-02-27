
def split_main_and_alt_conf_atom_groups(residue_group):
    return main_conf_atom_groups(residue_group), alt_conf_atom_groups(residue_group)

def main_conf_atom_groups(residue_group):
    ags = [ag for ag in residue_group.atom_groups() if ag.altloc=='']
    if len(ags) > 1: raise Exception('More than one main-conf atom group! {}'.format(residue_group.resid()))
    if len(ags) == 0: ags = None
    else:             ags = ags[0]
    return ags

def alt_conf_atom_groups(residue_group):
    ags = [ag for ag in residue_group.atom_groups() if ag.altloc]
    if len(ags) == 0: ags = None
    return ags
