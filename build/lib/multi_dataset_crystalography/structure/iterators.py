
def residues_via_conformers(pdb_hierarchy):
    """Iterate hierarchy and yield pure residues from conformers of chains (chain->conformer->residue)"""
    for chain in pdb_hierarchy.deep_copy().chains():
        for conf in chain.conformers():
            for res in conf.residues():
                yield res

def conformers_via_residue_groups(pdb_hierarchy):
    """Iterate hierarchy and yield pure conformers from residue groups (chain->residue_group->conf)"""
    for chain in pdb_hierarchy.deep_copy().chains():
        for rg in chain.residue_groups():
            for conf in rg.conformers():
                yield conf

def residue_groups_with_complete_set_of_conformers(pdb_hierarchy):
    all_confs = set(pdb_hierarchy.altloc_indices())
    all_confs.discard('')
    for rg in pdb_hierarchy.residue_groups():
        conf_altlocs = [c.altloc for c in rg.conformers()]
        if conf_altlocs == [''] or not all_confs.difference(conf_altlocs):
            yield rg

########################################################################################

def generate_residue_triplets(pdb_hierarchy):
    """Generate groups of three residues (this is an iterator)"""
    for chain in pdb_hierarchy.deep_copy().chains():
        # Skip if not protein
        if not chain.is_protein(): continue
        for conf in chain.conformers():
            # Reset for each conformer
            prev = current = next = None
            for res in conf.residues():
                # Shuffle all variables down one
                prev = current
                current = next
                next = res
                # Do we have three residues?
                if not (prev and current and next):
                    continue
                # Are they continuous residues
                if not (prev.resseq_as_int()+1 == current.resseq_as_int() == next.resseq_as_int()-1):
                    continue
                # Return the residue triplet
                yield (prev, current, next)
