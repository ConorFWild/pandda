import unittest

import numpy

from giant.structure.sequence import *


class TestSequence(unittest.TestCase):


    def setUp(self):
        # Reference sequence
        self.seq_a = 'SMSVKKPKRDDSKDLALCSMILTEMETHEDAWPFLLPVNLKLVPGYKKVIKKPMDFSTIREKLSSGQYPNLETFALDVRLVFDNCETFNEDDSDIGRAGHNMRKYFEKKWTDTFK'
        # Truncated sequence
        self.seq_b = 'SMSVKKPKRDDSKDLALCSMILTEMETHEDAWPFLLPVNLKLVPGYKKVIKKPMDFSTIREKLSSGQYPNLETFALDVRLVFDNCETFNEDDSDIGRAGHNMRKYFEKKW'
        # Mutated sequence
        self.seq_c = self.seq_a.replace('VNL','AAA')
        # Spliced with unrelated (half overlaps)
        self.seq_d = 'SMSVKKPKRDDSKDLALCSMILTEMETHEDAWPFLLPVAQNPNCNIMIFHPTKEEFNDFDKYIAYMESQGAHRAGLAKIIPPKEWKARETYDNISEILIATPLQQVASGRAGVFT'
        # Unrelated sequence
        self.seq_e = 'AQNPNCNIMIFHPTKEEFNDFDKYIAYMESQGAHRAGLAKIIPPKEWKARETYDNISEILIATPLQQVASGRAGVFTQYHKKKKAMTVGEYRHLANSKKYQTPPHQNFEDLERKYWKNRIYNSPIYGADISGSLFDENTKQWNLG'

        self.set_1 = [self.seq_a, self.seq_c, self.seq_d, self.seq_e]
        self.set_2 = [self.seq_b, self.seq_c, self.seq_d]

    def test_align_sequences_default_identical(self):
        ali = align_sequences_default(self.seq_a, self.seq_a)
        self.assertEqual(ali.a, self.seq_a)
        self.assertEqual(ali.b, self.seq_a)
        self.assertEqual(ali.calculate_sequence_identity(), 1.0)
        self.assertEqual(ali.matches(), '|'*len(self.seq_a))

    def test_align_sequences_default_truncated(self):
        ali = align_sequences_default(self.seq_a, self.seq_b)
        self.assertEqual(ali.a, self.seq_b)
        self.assertEqual(ali.b, self.seq_b)
        self.assertEqual(ali.calculate_sequence_identity(), 1.0)
        self.assertEqual(ali.matches(), '|'*len(self.seq_b))

    def test_align_sequences_default_mutation(self):
        ali = align_sequences_default(self.seq_a, self.seq_c)
        self.assertEqual(ali.a, self.seq_a)
        self.assertEqual(ali.b, self.seq_c)
        self.assertEqual(ali.calculate_sequence_identity(), 0.9739130434782609)
        self.assertEqual(ali.matches(), '|'*37+' '*3+'|'*75)

    def test_pairwise_sequence_identity_frac_sequence_alignment(self):
        arr = pairwise_sequence_identity(
                    seqs_1=self.set_1, seqs_2=self.set_2,
                    min_alignment=0.30, seq_identity_threshold=None)
        self.assertTrue(numpy.array_equal(
                            arr.round(5),
                            numpy.array([[ 1.        , 112/115.0 , 1.      ],
                                         [ 107/110.0 , 1.        , 38/39.0 ],
                                         [ 1.        , 38/39.0   , 1.      ],
                                         [ 0.        , 0.        , 1.      ]]).round(5)))

        arr = pairwise_sequence_identity(
                    seqs_1=self.set_1, seqs_2=self.set_2,
                    min_alignment=0.0, seq_identity_threshold=None)
        self.assertTrue(numpy.array_equal(
                            arr.round(5),
                            numpy.array([[ 1.        , 112/115.0 , 1.      ],
                                         [ 107/110.0 , 1.        , 38/39.0 ],
                                         [ 1.        , 38/39.0   , 1.      ],
                                         [ 6/18.0    , 6/18.0    , 1.      ]]).round(5)))

    def test_pairwise_sequence_identity_num_sequence_alignment(self):
        arr = pairwise_sequence_identity(
                    seqs_1=self.set_1, seqs_2=self.set_2,
                    min_alignment=19, seq_identity_threshold=None)
        self.assertTrue(numpy.array_equal(
                            arr.round(5),
                            numpy.array([[ 1.        , 112/115.0 , 1.      ],
                                         [ 107/110.0 , 1.        , 38/39.0 ],
                                         [ 1.        , 38/39.0   , 1.      ],
                                         [ 0.        , 0.        , 1.      ]]).round(5)))

    def test_pairwise_sequence_identity_seq_identity_threshold(self):
        arr = pairwise_sequence_identity(
                    seqs_1=self.set_1, seqs_2=self.set_2,
                    min_alignment=19, seq_identity_threshold=0.99)
        self.assertTrue(numpy.array_equal(
                            arr,
                            numpy.array([[ 1 , 0 , 1 ],
                                         [ 0 , 1 , 0 ],
                                         [ 1 , 0 , 1 ],
                                         [ 0 , 0 , 1 ]])))
