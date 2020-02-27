import os, sys, copy
import numpy
import iotbx.pdb
from scitbx.array_family import flex
from libtbx.utils import Sorry, Failure

from giant.structure.formatting import Labeller

#####################################################
###        MULTIPLE HIERARCHIES FUNCTIONS         ###
#####################################################

def find_unused_chain_ids(hierarchies):
    unused_chain_ids = iotbx.pdb.systematic_chain_ids()
    current_chain_ids=[]; [current_chain_ids.extend([ch.id for ch in h.chains()]) for h in hierarchies]
    current_chain_ids=list(set(current_chain_ids))
    for curr_id in current_chain_ids:   unused_chain_ids.remove(curr_id)
    return unused_chain_ids

#####################################################
###             HIERARCHY FUNCTIONS               ###
#####################################################

def resolve_residue_id_clashes(fixed_hierarchy, moving_hierarchy, in_place=False, verbose=False):
    """Move residues in mov_hierarchiy to new chains if they have the same resid as a residue in fixed_hierarchy but different resnames"""
    if (not in_place): moving_hierarchy = moving_hierarchy.deep_copy()
    # Sort all residues (by chain then id) for the fixed hierarchy - chains > residue_ids > residue_groups (objects)
    fixed_dict = {c.id:{} for c in fixed_hierarchy.only_model().chains()}
    [fixed_dict.get(rg.parent().id).setdefault(rg.resid(), []).append(rg) for rg in fixed_hierarchy.residue_groups()]
    # Find the residues with clashing resids
    residues_to_update = []
    for rg_mov in moving_hierarchy.only_model().residue_groups():
        # Extract equivalent group for this residue
        rg_ref = fixed_dict.get(rg_mov.parent().id,{}).get(rg_mov.resid(),[])
        # Check that we've found only one (rather than silently erroring)
        if len(rg_ref) == 0:
            continue
        elif len(rg_ref) > 1:
            err_str =  "A residue is present more than once in the output hierarchy.\n"
            err_str += "> Moving residue: \n\t{} ({})\n".format(Labeller.format(rg_mov), ', '.join(rg_mov.unique_resnames()))
            err_str += '> Target residues: \n\t{}\n'.format('\n\t'.join(['{} ({})'.format(Labeller.format(r), ', '.join(r.unique_resnames())) for r in rg_ref]))
            err_str += "Each label combination can only be present once unless labelled with appropriate alternate conformations.\n"
            raise Sorry(err_str)
        rg_ref = rg_ref[0]
        # Check to see if the residue is the same type as in the reference structure
        if map(str.strip, rg_ref.unique_resnames()) == map(str.strip, rg_mov.unique_resnames()):
            # Same residue -- that's fine
            # TODO allow this if the resnames match OR if it's a protein residue? (allows mutations/reactions?)
            continue
        else:
            # Will need to be added to different chain
            if verbose: print 'Different residues with same id - changing chains: {} != {}'.format(list(rg_ref.unique_resnames()), list(rg_mov.unique_resnames()))
            residues_to_update.append(rg_mov)
    # Nothing to do -- return
    if not residues_to_update:
        return moving_hierarchy
    # New chain to add clashing groups to
    new_chain = None
    old_chain_prev = None
    new_chain_ids = find_unused_chain_ids(hierarchies=[fixed_hierarchy, moving_hierarchy])
    # Go through and transfer the residue groups to new chains
    for rg_mov in residues_to_update:
        old_chain = rg_mov.parent()
        # See if there is a residue with this id already present in the new chain
        if (not new_chain) or (rg_mov.resid() in new_chain.get_residue_ids()) or (old_chain_prev.id != old_chain.id):
            new_chain = iotbx.pdb.hierarchy.chain(id=new_chain_ids.pop(0))
            old_chain.parent().append_chain(new_chain)
        if verbose: print '{} - moving to chain {}'.format(Labeller.format(rg_mov), new_chain.id)
        # Remove from old chain and add to the new chain
        old_chain.remove_residue_group(rg_mov)
        new_chain.append_residue_group(rg_mov)
        # Keep track of old_chain of previous residue
        old_chain_prev = old_chain

    return moving_hierarchy

def transfer_residue_groups_from_other(acceptor_hierarchy, donor_hierarchy, in_place=False, verbose=False):
    """Transfer atom_groups from donor_hierarchy to matching residue_groups in acceptor_hierarchy, creating new chains and residue groups only where necessary"""
    if not in_place: acceptor_hierarchy = acceptor_hierarchy.deep_copy()
    # Sort all residues (by chain then id) for the acceptor hierarchy
    accept_model = acceptor_hierarchy.only_model()
    accept_dict = {c.id:{} for c in accept_model.chains()}
    [accept_dict.get(rg.parent().id).setdefault(rg.resid(), []).append(rg) for rg in accept_model.residue_groups()]
    # Dictionary to link matching chains (allows multiple chain As to be linked uniquely to multiple chain As)
    link_dict = {}
    # Residues that don't have a matching partner in the old hierarchy
    tricky_rgs = []
    # Iterate through donor chains
    for donor_ch in donor_hierarchy.only_model().chains():
        # If chain not in hierarchy, simply copy across
        if accept_dict.get(donor_ch.id, None) is None:
            if verbose: print 'Transferring whole chain:    {}'.format(Labeller.format(donor_ch))
            accept_model.append_chain(donor_ch.detached_copy())
            continue
        # Chain present, copy by residue_group
        for donor_rg in donor_ch.residue_groups():
            # Find equivalent residue groups in the other hierarchy
            accept_rg = accept_dict.get(donor_ch.id).get(donor_rg.resid(),[])
            if len(accept_rg) > 1:
                # Should only be one...
                raise Exception('More than one residue group in hierarchy with the same residue_id and chain_id')
            elif len(accept_rg) == 1:
                accept_rg = accept_rg[0]
                # Record the links between these chains
                link_dict.setdefault(donor_rg.parent(), accept_rg.parent())
                # Transfer atom groups to this residue_group
                if verbose: print 'Transferring atom groups:    {} > {}'.format(Labeller.format(donor_rg), Labeller.format(accept_rg))
                for donor_ag in donor_rg.atom_groups():
                    accept_rg.append_atom_group(donor_ag.detached_copy())
            else:
                # Have the possibility of multiple chains with the same id, so at the moment, store for later
                tricky_rgs.append(donor_rg)
    # Transfer residues that have chain matches, but don't have residue matches in the acceptor structures
    for donor_rg in tricky_rgs:
        # Get chain from link_dict
        accept_ch = link_dict.get(donor_rg.parent().id, None)
        # If the chain isn't linked:
        if accept_ch is None:
            # If there's only one chain with the same ID, choose this one
            possible_chains = [c for c in accept_model.chains() if c.id == donor_rg.parent().id]
            if len(possible_chains) == 1:
                accept_ch = possible_chains[0]
            else:
                raise Failure("Don't know how to transfer {} to the output model".format(Labeller.format(donor_rg)))
        # Simply append to chain
        if verbose: print 'Transferring residue group:  {} > {}'.format(Labeller.format(donor_rg), Labeller.format(accept_ch))
        accept_ch.append_residue_group(donor_rg.detached_copy())

    return acceptor_hierarchy

#def compare_hierarchies(reference_hierarchy, query_hierarchy):
#    """Find which residues are the same, different, or unique to each of the hierarchies"""
#
#    # Can't be bothered to type those names all the time
#    ref_h = reference_hierarchy; qry_h = query_hierarchy
#    # Iterate through the reference hierarchy and find which residues have changed/are unique/are the same
#    ref_only = []; ref_diff = []; ref_same = []
#
#    # Iterate through residues
#    for rg_ref in ref_h.residue_groups():
#        # Check that the residue groups only have one residue name in each
#        assert len(rg_ref.unique_resnames()) == 1, "Can't handle residues with more than one residue name"
#        # Extract label for this residue
#        resid = rg_ref.resid(); chainid = rg_ref.parent().id
#        # Extract same rg for qry
#        rg_qry = [rg for rg in qry_h.residue_groups() if rg.resid()==resid and rg.parent().id==chainid]
#
#        # MORE THAN ONE MATCHING RESIDUE -- ERROR?
#        if len(rg_qry) > 1:
#            raise Exception('MORE THAN ONE MATCHING')
#        # PRESENT ONLY IN REFERENCE
#        elif len(rg_qry) == 0:
#            ref_only.append((chainid, resid))
#        # PRESENT IN BOTH
#        elif len(rg_qry) == 1:
#            # Check to see if the residue has moved
#            rg_qry = rg_qry[0]
#            # Check that the residue name is the same
#            assert map(str.strip, rg_ref.unique_resnames()) == map(str.strip,rg_qry.unique_resnames()), 'Cannot Merge - Different Residues!: {} != {}'.format(list(rg_ref.unique_resnames()), list(rg_qry.unique_resnames()))
#            # Extract atoms and check to see if the same
#            rg_ref_ats = rg_ref.atoms()
#            rg_qry_ats = rg_qry.atoms()
#            # Check to see if the same length
#            if len(rg_ref_ats) != len(rg_qry_ats):
#                # There's been a change in the structure
#                ref_diff.append((chainid, resid))
#            else:
#                # Extract coordinates for the atoms
#                rg_ref_xyz = rg_ref_ats.extract_xyz()
#                rg_qry_xyz = rg_qry_ats.extract_xyz()
#                # Calculate the rmsd of the atoms (ASSUMING IN THE SAME ORDER)
#                rmsd = (rg_ref_xyz - rg_qry_xyz).norm()
#                # Check to see if the atoms have moved
#                if rmsd != 0.0:
#                    # There's been a change in the structure
#                    ref_diff.append((chainid, resid))
#                else:
#                    # Same between the two structures
#                    ref_same.append((chainid, resid))
#
#    return ref_only, ref_diff, ref_same
#
