
import numpy
import mmtbx.alignment

def align_sequences_default(seq_a, seq_b):
    """Align two sequences using a default setting"""
    if isinstance(seq_a, list): seq_a=''.join(seq_a)
    if isinstance(seq_b, list): seq_b=''.join(seq_b)
    # Align the sequences of the two chains
    return mmtbx.alignment.align(
                    seq_a=seq_a, seq_b=seq_b,
                    gap_opening_penalty = 20,
                    gap_extension_penalty = 2,
                    similarity_function = 'blosum50',
                    style = 'local').extract_alignment()

def pairwise_sequence_identity(seqs_1, seqs_2, min_alignment=0.90, seq_identity_threshold=None):
    """
    Align all sequences in seqs_1 to all sequences in seqs_2
    Returns an array giving the pairwise sequence identity between the sets.
    If seq_identity_threshold is given, the array contains boolean values indicating whether it is above that value
    If min_alignment is given:
        as a float of 0.0-1.0   - the alignment must cover this fraction of the shortest_sequence
        as an integer           - the alignment must cover this number of amino acids
    """
    min_ali_frac = min_ali_num = 0
    if isinstance(min_alignment, float):
        assert 0.0<=min_alignment<=1.0, 'min_alignment must be either an integer or between 0.0-1.0 or None'
        min_ali_frac = min_alignment
    elif isinstance(min_alignment, int):
        min_ali_num  = min_alignment
    else:
        assert min_alignment is None, 'min_alignment must be either an integer or between 0.0-1.0 or None'

    arr = numpy.zeros((len(seqs_1),len(seqs_2)), dtype=float)
    for i1, s1 in enumerate(seqs_1):
        if not s1: continue
        for i2, s2 in enumerate(seqs_2):
            if not s2: continue
            ali = align_sequences_default(s1,s2)
            align_num = len(ali.match_codes)
            if align_num >= (min_ali_num + min_ali_frac*min(len(s1),len(s2))):
                arr[i1,i2] = ali.calculate_sequence_identity()
    if seq_identity_threshold is not None:
        return (arr>seq_identity_threshold).astype(int)
    return arr

def pairwise_chain_sequence_identity(chains_1, chains_2, seq_identity_threshold=None):
    """
    Align all chains in chains_1 to all chains in chains_2
    Returns an array giving the pairwise sequence identity between the chains.
    If seq_identity_threshold is given, the array contains boolean values indicating whether it is above that value
    """
    seqs_1 = [''.join(c.as_sequence()) for c in chains_1]
    seqs_2 = [''.join(c.as_sequence()) for c in chains_2]
    return pairwise_sequence_identity(seqs_1, seqs_2, seq_identity_threshold=seq_identity_threshold)
