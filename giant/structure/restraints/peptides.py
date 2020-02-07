import numpy

from scitbx.array_family import flex
from giant.structure.select import backbone, extract_backbone_atoms
from giant.structure.dihedrals import get_peptide_bond_type

#Lengths
IDEAL_CA_C = 1.52
IDEAL_C_N  = 1.33
IDEAL_N_CA = 1.45
# Angles
IDEAL_CA_C_N = 116
IDEAL_C_N_CA = 122

"""
 1 -  6        Record name     "LINK  "
 13 - 16        Atom            name1       Atom name.
 17             Character       altLoc1     Alternate location indicator.
 18 - 20        Residue name    resName1    Residue name.
 22             Character       chainID1    Chain identifier.
 23 - 26        Integer         resSeq1     Residue sequence number.
 27             AChar           iCode1      Insertion code.
 43 - 46        Atom            name2       Atom name.
 47             Character       altLoc2     Alternate location indicator.
 48 - 50        Residue name    resName2    Residue name.
 52             Character       chainID2    Chain identifier.
 53 - 56        Integer         resSeq2     Residue sequence number.
 57             AChar           iCode2      Insertion code.
 60 - 65        SymOP           sym1        Symmetry operator for 1st atom.
 67 - 72        SymOP           sym2        Symmetry operator for 2nd atom.

"""

def format_link_record(atom_1, atom_2, chain_id_1=None, chain_id_2=None, link_type=None):
    """Generate a LINK record for two atoms"""

    link_template = "LINK        {0.name:4}{0.altloc:1}{0.resname:3} {0.chain_id:1}{0.resseq:4}{0.icode:1}               {1.name:4}{1.altloc:1}{1.resname:3} {1.chain_id:1}{1.resseq:4}{1.icode:1}  {2:6} {3:6}{4:7}"

    a1 = atom_1.fetch_labels()
    a2 = atom_2.fetch_labels()
    if chain_id_1 is not None: a1.chain_id=chain_id_1
    if chain_id_2 is not None: a2.chain_id=chain_id_2
    if link_type is None: link_type = ''
    return link_template.format(a1,a2,'','',link_type)

def generate_set_of_alternate_conformer_peptide_links(hierarchy):
    """
    Identify discontinuous peptide chain conformations and generate link records to fix them:
    Returns (link_restraints, warnings)
        link_restraints - list of link records
        warnings        - list of warning strings
    """

    warnings = []
    link_restraints = []

    # Find any problems with chains
    break_list = identify_model_discontinuities(hierarchy)

    # Go through one chain at a time
    for chn_id in sorted(break_list.keys()):
        # Find any breaks or residues without connections
        atom_pairs, no_conns = suggest_continuity_connections(break_list[chn_id])
        # Create link records for breaks
        for a1,a2 in atom_pairs:
            #print 'Generate LINK for {} of {} to {} of {}'.format(a1.parent().altloc, a1.parent().parent().resid(), a2.parent().altloc, a2.parent().parent().resid())
            link_type = get_peptide_bond_type(prev_Ca=a1.parent().get_atom('CA'), prev_C=a1, curr_N=a2, curr_Ca=a2.parent().get_atom('CA'))
            lr = (a1.fetch_labels().detached_copy(), a2.fetch_labels().detached_copy(), chn_id, chn_id, link_type)
            link_restraints.append(lr)
        # Create warnings for no connections
        for a1,a2 in no_conns:
            m1 = m2 = '?'
            if a1:
                m1 = a1.parent().parent().resid()+'-'+a1.parent().altloc
                m2 = 'next residue'
            if a2:
                m1 = 'previous residue'
                m2 = a2.parent().parent().resid()+'-'+a2.parent().altloc
            warnings.append('Chain {}: {} has no link to {}'.format(chn_id, m1, m2))

    return link_restraints, warnings

def identify_model_discontinuities(hierarchy):
    """Find discontinuities in the protein backbone where alternate conformers do not match from residue to residue"""

    all_breaks = {}
    for chn in backbone(hierarchy).chains():
        if not chn.is_protein(): continue
        break_list = identify_chain_discontinuities(chain=chn)
        all_breaks[chn.id] = break_list
    return all_breaks

def identify_chain_discontinuities(chain):
    """Find residues that are discontinuous for conformers"""

    assert chain.is_protein()

    # Initialise empty variables
    prev = prev_c = None
    # List of chain breaks
    break_list = []

    # Iterate through residue groups
    for curr in chain.residue_groups():

        # Extract N & C for this residue (do atom-wise rather than using conformer objects)
        curr_n = []; curr_c = []
        for ag in curr.atom_groups():
            n = ag.get_atom('N')
            if n: curr_n.append(n)
            c = ag.get_atom('C')
            if c: curr_c.append(c)

        # Previous residue and linked to previous
        if prev_c and curr.link_to_previous:

            # Assume there's a break
            model_break = True

            altlocs_prev = set([a.parent().altloc for a in prev_c])
            altlocs_curr = set([a.parent().altloc for a in curr_n])

            if '' in altlocs_prev: assert len(altlocs_prev) == 1
            if '' in altlocs_curr: assert len(altlocs_curr) == 1

            # Skip if all main chain on either side
            if ({''} in [altlocs_prev, altlocs_curr]):
                model_break = False
            # Skip if altlocs match perfectly
            elif not altlocs_prev.symmetric_difference(altlocs_curr):
                model_break = False

            # Append to return list
            if model_break:
                break_list.append((prev, curr))

        # END OF LOOP - UPDATE PREV
        prev   = curr
        prev_c = curr_c

    return break_list

def suggest_continuity_connections(residue_group_pairs):
    """
    Generate a set of links that need to be made to make model continuous.
    Returns (connections_to_make, missing_connections)
        connections_to_make - a list of atom pairs (connections to make)
        missing_connections - a list of atoms with no connecting partner
    """

    connections_to_make = []
    missing_connections = []
    # Given a set of residue pairs, find connections between them based on likely connectivity
    for rg_1, rg_2 in residue_group_pairs:
        atom_pairs, no_links = suggest_residue_groups_connections(rg_prev=rg_1, rg_curr=rg_2, verbose=True)
        #print 'RESIDUE:', rg_1.resid(), '-', rg_2.resid()
        #print 'LINK:', ', '.join([a1.parent().altloc+'-'+a2.parent().altloc for a1,a2 in atom_pairs])
        #print 'NONE:', ', '.join([a1.parent().parent().resid()+'-'+a1.parent().altloc+':?' if a1 else '?:'+a2.parent().parent().resid()+'-'+a2.parent().altloc for a1,a2 in no_links])
        connections_to_make.extend(atom_pairs)
        missing_connections.extend(no_links)

    return connections_to_make, missing_connections

def suggest_residue_groups_connections(rg_prev, rg_curr, tol_dist=0.05, tol_ang=5, chi_cut=10, verbose=True):
    """Assess the type and number of chains breaks between a set of carbons and nitrogens of the backbone"""

    prev_Cs = [ag.get_atom('C') for ag in rg_prev.atom_groups() if ag.get_atom('C') is not None]
    curr_Ns = [ag.get_atom('N') for ag in rg_curr.atom_groups() if ag.get_atom('N') is not None]

    altlocs_prev_C = [a.parent().altloc for a in prev_Cs]
    altlocs_curr_N = [a.parent().altloc for a in curr_Ns]

    common_set = set(altlocs_prev_C).intersection(altlocs_curr_N)
    diff_set_prev_C = set(altlocs_prev_C).difference(altlocs_curr_N)
    diff_set_curr_N = set(altlocs_curr_N).difference(altlocs_prev_C)

    assert diff_set_prev_C.issubset(altlocs_prev_C)
    assert diff_set_curr_N.issubset(altlocs_curr_N)

#    if verbose:
#        print '========>'
#        print 'Common Set: {:10}'.format(','.join(common_set))
#        print 'Unique res-c-:   {:10}'.format(','.join(diff_set_prev_C))
#        print 'Unique -n-res:   {:10}'.format(','.join(diff_set_curr_N))

    # Create matrix of closeness to idealness
    rmsds = numpy.zeros((len(prev_Cs),len(curr_Ns)))

    # Calculate the RMS from the ideal bonds and angles
    for i_p, p_C in enumerate(prev_Cs):
        # Get the previous C-ALPHA
        p_Ca = p_C.parent().get_atom('CA')
        # Match to possible Nitrogens
        for i_c, c_N in enumerate(curr_Ns):
            # If same altloc, set to zero
            if p_C.parent().altloc == c_N.parent().altloc:
                rmsds[i_p,i_c] = 0
                continue
            # Get the current C-ALPHA
            c_Ca = c_N.parent().get_atom('CA')
            assert p_Ca and c_Ca
            # Calculate bond lengths and angles
            dist_pC_cN    = p_C.distance(c_N)
            ang_pCA_pC_cN = p_C.angle(p_Ca, c_N, deg=True)
            ang_pC_cN_cCA = c_N.angle(p_C, c_Ca, deg=True)
            # Calculate errors in the bond lengths and angles
            e_d  = (dist_pC_cN - IDEAL_C_N) / tol_dist
            e_a1 = (ang_pCA_pC_cN - IDEAL_CA_C_N) / tol_ang
            e_a2 = (ang_pC_cN_cCA - IDEAL_C_N_CA) / tol_ang
            #print dist_pC_cN, e_d, ang_pCA_pC_cN, e_a1, ang_pC_cN_cCA, e_a2
            # Store the RMS value
            rmsds[i_p,i_c] = numpy.sqrt(numpy.sum(numpy.square([e_d, e_a1, e_a2])))

    # Create matrix of restrain links
    links = rmsds < chi_cut
    # Trim unnecessary links
    links = _trim_link_matrix(link_mx=links,   set1=altlocs_prev_C, set2=altlocs_curr_N, common_set=common_set)
    links = _trim_link_matrix(link_mx=links.T, set1=altlocs_curr_N, set2=altlocs_prev_C, common_set=common_set).T

    # Create list of atoms that need to be connected
    atom_pairs = [(prev_Cs[i], curr_Ns[j]) for i,j in zip(*numpy.where(links))]

    no_links = [(prev_Cs[i],None) for i,r in enumerate(links) if sum(r) == 0] + \
               [(None,curr_Ns[i]) for i,r in enumerate(links.T) if sum(r) == 0]

    return atom_pairs, no_links

def _trim_link_matrix(link_mx, set1, set2, common_set=[]):
    # Decide which links to keep
    for i, s1 in enumerate(set1):
        # Flag if no connections
        if sum(link_mx[i,:]) == 0:
            continue
        # Keep if only one link
        if sum(link_mx[i,:]) == 1:
            continue
        # Remove links with non-common-set if possible
        # (Preferably leave common set unlinked - most likely to be continuous chain)
        if s1 in common_set:

            for j, s2 in enumerate(set2):
                # Always leave common_set links
                if s1 == s2:
                    continue
                # Do nothing if no link
                if not link_mx[i,j]:
                    continue
                # Remove link if target has more than one link
                if sum(link_mx[:,j]) > 1:
                    link_mx[i, j] = 0

    return link_mx

