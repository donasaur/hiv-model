from process.AlternativeSplicing import *
from mainaux.TestHelpers import *
from state.Proteins import *
import unittest
import numpy as np

class TestAlternativeSplicing(unittest.TestCase):
    
    # This method is run prior to each test automatically.
    def setUp(self):
        # Default state
        self.state = State()
        self.splice_proc = AlternativeSplicing(self.state)

        mRNA_state = self.state.get_state('mRNAs')
        protein_state = self.state.get_state('proteins')

        self.proteins_nuc = protein_state.proteins_nuc # proteins_nuc[Proteins.index['Rev']] is Rev count in nucleus
        self.proteins_cyt = protein_state.proteins_cyt # proteins_cyt[Proteins.index['Rev']] is Rev count in cytoplasm

        self.full_len_transcripts_nuc = mRNA_state.full_len_transcripts_nuc
        self.full_len_transcripts_cyt = mRNA_state.full_len_transcripts_cyt
        self.single_splice_transcript_nuc = mRNA_state.single_splice_transcript_nuc
        self.single_splice_transcript_cyt = mRNA_state.single_splice_transcript_cyt
        self.multi_splice_transcript_nuc = mRNA_state.multi_splice_transcript_nuc
        self.multi_splice_transcript_cyt = mRNA_state.multi_splice_transcript_cyt
        
        # S1 state
        self.s1_state = s1_state()
        self.s1_splice_proc = AlternativeSplicing(self.s1_state)

        mRNA_state = self.s1_state.get_state('mRNAs')
        protein_state = self.s1_state.get_state('proteins')

        self.s1_proteins_nuc = protein_state.proteins_nuc # proteins_nuc[Proteins.index['Rev']] is Rev count in nucleus
        self.s1_proteins_cyt = protein_state.proteins_cyt # proteins_cyt[Proteins.index['Rev']] is Rev count in cytoplasm

        self.s1_full_len_transcripts_nuc = mRNA_state.full_len_transcripts_nuc
        self.s1_full_len_transcripts_cyt = mRNA_state.full_len_transcripts_cyt
        self.s1_single_splice_transcript_nuc = mRNA_state.single_splice_transcript_nuc
        self.s1_single_splice_transcript_cyt = mRNA_state.single_splice_transcript_cyt
        self.s1_multi_splice_transcript_nuc = mRNA_state.multi_splice_transcript_nuc
        self.s1_multi_splice_transcript_cyt = mRNA_state.multi_splice_transcript_cyt

        # S2 state
        self.s2_state = s2_state()
        self.s2_splice_proc = AlternativeSplicing(self.s2_state)

        mRNA_state = self.s2_state.get_state('mRNAs')
        protein_state = self.s2_state.get_state('proteins')

        self.s2_proteins_nuc = protein_state.proteins_nuc # proteins_nuc[Proteins.index['Rev']] is Rev count in nucleus
        self.s2_proteins_cyt = protein_state.proteins_cyt # proteins_cyt[Proteins.index['Rev']] is Rev count in cytoplasm

        self.s2_full_len_transcripts_nuc = mRNA_state.full_len_transcripts_nuc
        self.s2_full_len_transcripts_cyt = mRNA_state.full_len_transcripts_cyt
        self.s2_single_splice_transcript_nuc = mRNA_state.single_splice_transcript_nuc
        self.s2_single_splice_transcript_cyt = mRNA_state.single_splice_transcript_cyt
        self.s2_multi_splice_transcript_nuc = mRNA_state.multi_splice_transcript_nuc
        self.s2_multi_splice_transcript_cyt = mRNA_state.multi_splice_transcript_cyt

        # S3 state
        self.s3_state = s3_state()
        self.s3_splice_proc = AlternativeSplicing(self.s3_state)

        mRNA_state = self.s3_state.get_state('mRNAs')
        protein_state = self.s3_state.get_state('proteins')

        self.s3_proteins_nuc = protein_state.proteins_nuc # proteins_nuc[Proteins.index['Rev']] is Rev count in nucleus
        self.s3_proteins_cyt = protein_state.proteins_cyt # proteins_cyt[Proteins.index['Rev']] is Rev count in cytoplasm

        self.s3_full_len_transcripts_nuc = mRNA_state.full_len_transcripts_nuc
        self.s3_full_len_transcripts_cyt = mRNA_state.full_len_transcripts_cyt
        self.s3_single_splice_transcript_nuc = mRNA_state.single_splice_transcript_nuc
        self.s3_single_splice_transcript_cyt = mRNA_state.single_splice_transcript_cyt
        self.s3_multi_splice_transcript_nuc = mRNA_state.multi_splice_transcript_nuc
        self.s3_multi_splice_transcript_cyt = mRNA_state.multi_splice_transcript_cyt

        # S4 state
        self.s4_state = s4_state()
        self.s4_splice_proc = AlternativeSplicing(self.s4_state)

        mRNA_state = self.s4_state.get_state('mRNAs')
        protein_state = self.s4_state.get_state('proteins')

        self.s4_proteins_nuc = protein_state.proteins_nuc # proteins_nuc[Proteins.index['Rev']] is Rev count in nucleus
        self.s4_proteins_cyt = protein_state.proteins_cyt # proteins_cyt[Proteins.index['Rev']] is Rev count in cytoplasm

        self.s4_full_len_transcripts_nuc = mRNA_state.full_len_transcripts_nuc
        self.s4_full_len_transcripts_cyt = mRNA_state.full_len_transcripts_cyt
        self.s4_single_splice_transcript_nuc = mRNA_state.single_splice_transcript_nuc
        self.s4_single_splice_transcript_cyt = mRNA_state.single_splice_transcript_cyt
        self.s4_multi_splice_transcript_nuc = mRNA_state.multi_splice_transcript_nuc
        self.s4_multi_splice_transcript_cyt = mRNA_state.multi_splice_transcript_cyt
        
        # Other setup stuff
        self.index = Proteins.index

    # Start out with a fixed amt of Rev total in the system (e.g. 3)
    # See if the same amount of Rev remains in the system after AlternativeSplicing is run a couple of times
    def test_Rev_balance(self):
        for i in range(30):
            prev_Rev_count = count_total_Rev(self.s1_state)
            self.s1_splice_proc.evolve_state(i)
            curr_Rev_count = count_total_Rev(self.s1_state)
            
            self.assertEqual(prev_Rev_count, curr_Rev_count)

    # See whether given enough time (by intensifying rate of transfer from one bin to another)
    # whether all balls will trickle down to m-tat, m-rev, and m-nef bins
    def test_trickle_to_leaves(self):
        self.s1_splice_proc.PROB_SPLICE_FULL_TO_SINGLE = 1
        self.s1_splice_proc.PROB_SPLICE_SINGLE_TO_MULTI = 1
        self.s1_splice_proc.PROB_VPR_THIRD_SPLICE = 1
        self.s1_splice_proc.PROB_VIF_THIRD_SPLICE = 1
        self.s1_multi_splice_transcript_nuc[[3,4,9,10,14,15]] = 0
        self.s1_single_splice_transcript_nuc[np.arange(4, 63, 7)] = 0
        self.s1_single_splice_transcript_nuc[np.arange(5, 63, 7)] = 0
        total_mRNA_count_in_nuc = count_total_mRNA(self.s1_state, True)[0]
        for i in range(30):
            self.s1_splice_proc.evolve_state(i)
        sum_of_terminal_bins = np.sum(self.s1_multi_splice_transcript_nuc[[1,2,5,7,8,11,12,13,16]])

        self.assertEqual(total_mRNA_count_in_nuc, sum_of_terminal_bins)


    # Start with fixed num of full-length mRNA
    # At every timestep before splicing process is run,
    # record initial abundance at a particular node
    # and make sure that sum of lower-tier abundances
    # is less than initial abundance
    def test_sum_child_nodes_smaller(self):
        mRNA_state = self.state.get_state('mRNAs')
        mRNA_state.full_len_transcripts_nuc = np.array([230, 100, 230, 400, 200, 300, 200, 300, 150])
        initial_mRNA_sum = np.sum(mRNA_state.full_len_transcripts_nuc)

        for i in range(30):
            self.splice_proc.evolve_state(i)

        self.assertTrue(np.sum(self.single_splice_transcript_nuc) + np.sum(self.multi_splice_transcript_nuc) <= initial_mRNA_sum)

        self.setUp()

        self.s1_multi_splice_transcript_nuc[[0,6,1,2,5,7,8,11,12,13,16]] = 0
        initial_mRNA_sum = np.sum(self.s1_full_len_transcripts_nuc) + np.sum(self.s1_single_splice_transcript_nuc)

        for i in range(30):
            self.s1_splice_proc.evolve_state(i)

        self.assertTrue(np.sum(self.s1_multi_splice_transcript_nuc[[0,6,1,2,5,7,8,11,12,13,16]]) <= initial_mRNA_sum)


    # Make sure that seeds placed at a particular ssMRNA
    # will trickle down to only descendants
    def test_single_splice_trickle(self):
        self.single_splice_transcript_nuc[np.arange(0,63,7)] = 100
        initial_mRNA_sum = count_total_mRNA(self.state, True)[0]

        for i in range(15):
            self.splice_proc.evolve_state(i)

        self.assertEqual(np.sum(self.multi_splice_transcript_nuc[[0,1,2,5]]) + np.sum(self.single_splice_transcript_nuc[np.arange(0,63,7)]), initial_mRNA_sum)

        self.setUp()

        self.single_splice_transcript_nuc[np.arange(1,63,7)] = 100
        initial_mRNA_sum = count_total_mRNA(self.state, True)[0]

        for i in range(15):
            self.splice_proc.evolve_state(i)

        self.assertEqual(np.sum(self.multi_splice_transcript_nuc[[6,7,8,11]]) + np.sum(self.single_splice_transcript_nuc[np.arange(1,63,7)]), initial_mRNA_sum)


        self.setUp()

        self.single_splice_transcript_nuc[np.arange(2,63,7)] = 100
        initial_mRNA_sum = count_total_mRNA(self.state, True)[0]

        for i in range(15):
            self.splice_proc.evolve_state(i)

        self.assertEqual(np.sum(self.multi_splice_transcript_nuc[12]) + np.sum(self.single_splice_transcript_nuc[np.arange(2,63,7)]), initial_mRNA_sum)


        self.setUp()

        self.single_splice_transcript_nuc[np.arange(3,63,7)] = 100
        initial_mRNA_sum = count_total_mRNA(self.state, True)[0]

        for i in range(15):
            self.splice_proc.evolve_state(i)

        self.assertEqual(np.sum(self.multi_splice_transcript_nuc[13]) + np.sum(self.single_splice_transcript_nuc[np.arange(3,63,7)]), initial_mRNA_sum)

        self.setUp()

        self.single_splice_transcript_nuc[np.arange(4,63,7)] = 100
        initial_mRNA_sum = count_total_mRNA(self.state, True)[0]

        for i in range(15):
            self.splice_proc.evolve_state(i)

        self.assertEqual(np.sum(self.multi_splice_transcript_nuc[16]) + np.sum(self.single_splice_transcript_nuc[np.arange(4,63,7)]), initial_mRNA_sum)


    # Make sure total number of mRNA in the system remains constant
    def test_total_abundance(self):
        for i in range(30):
            prev_mRNA_count_in_nuc = count_total_mRNA(self.s1_state, True)[0]
            self.s1_splice_proc.evolve_state(i)
            curr_mRNA_count_in_nuc = count_total_mRNA(self.s1_state, True)[0]

            self.assertEqual(prev_mRNA_count_in_nuc, curr_mRNA_count_in_nuc)

    # Test whether if under the appropriate rate constants whether or not
    # MRNA will get stuck at a certain level of the tree hierarchy
    def test_stuck_at_ssMRNA(self):
        self.s1_splice_proc.PROB_SPLICE_FULL_TO_SINGLE = 1
        self.s1_splice_proc.PROB_SPLICE_SINGLE_TO_MULTI = 0

        # Zero out bins
        mRNA_state = self.s1_state.get_state('mRNAs')
        mRNA_state.multi_splice_transcript_nuc = mRNA_state.multi_splice_transcript_nuc * 0

        total_mRNA_count_in_nuc = count_total_mRNA(self.s1_state, True)[0]

        for i in range(5):
            self.s1_splice_proc.evolve_state(i)
        sum_of_terminal_bins = np.sum(self.s1_single_splice_transcript_nuc)

        self.assertEqual(total_mRNA_count_in_nuc, sum_of_terminal_bins)
        # Don't need an assert equals 0 for the other bins because other tests already cover
        # mRNA balance + ensures no negative abundances

    def test_stuck_at_full_len(self):
        self.s1_splice_proc.PROB_SPLICE_FULL_TO_SINGLE = 0

        # Zero out bins
        mRNA_state = self.s1_state.get_state('mRNAs')
        mRNA_state.multi_splice_transcript_nuc = mRNA_state.multi_splice_transcript_nuc * 0
        mRNA_state.single_splice_transcript_nuc = mRNA_state.single_splice_transcript_nuc * 0

        total_mRNA_count_in_nuc = count_total_mRNA(self.s1_state, True)[0]

        for i in range(5):
            self.s1_splice_proc.evolve_state(i)
        sum_of_terminal_bins = np.sum(self.s1_full_len_transcripts_nuc)

        self.assertEqual(total_mRNA_count_in_nuc, sum_of_terminal_bins)

    # See how changing F values affects what full-length mRNA becomes
    def test_single_slice_partition(self):
        F1 = 1
        F2 = 0
        F3 = 0
        F4 = 0
        F5 = 0
        self.splice_proc.CUMULATIVE_F1_TO_F5 = np.cumsum([F1, F2, F3, F4, F5])
        self.splice_proc.PROB_SPLICE_SINGLE_TO_MULTI = 0

        self.full_len_transcripts_nuc[:] = 100
        initial_mRNA_sum = count_total_mRNA(self.state, True)[0]

        for i in range(30):
            self.splice_proc.evolve_state(i)

        self.assertEqual(np.sum(self.full_len_transcripts_nuc) + sum_ss_type(self.state, 0), initial_mRNA_sum)

        self.setUp()

        F1 = 0
        F2 = 1
        F3 = 0
        F4 = 0
        F5 = 0
        self.splice_proc.CUMULATIVE_F1_TO_F5 = np.cumsum([F1, F2, F3, F4, F5])
        self.splice_proc.PROB_SPLICE_SINGLE_TO_MULTI = 0

        self.full_len_transcripts_nuc[:] = 100
        initial_mRNA_sum = count_total_mRNA(self.state, True)[0]

        for i in range(30):
            self.splice_proc.evolve_state(i)

        self.assertEqual(np.sum(self.full_len_transcripts_nuc) + sum_ss_type(self.state, 1), initial_mRNA_sum)

        self.setUp()

        F1 = 0
        F2 = 0
        F3 = 1
        F4 = 0
        F5 = 0
        self.splice_proc.CUMULATIVE_F1_TO_F5 = np.cumsum([F1, F2, F3, F4, F5])
        self.splice_proc.PROB_SPLICE_SINGLE_TO_MULTI = 0

        self.full_len_transcripts_nuc[:] = 100
        initial_mRNA_sum = count_total_mRNA(self.state, True)[0]

        for i in range(30):
            self.splice_proc.evolve_state(i)

        self.assertEqual(np.sum(self.full_len_transcripts_nuc) + sum_ss_type(self.state, 2), initial_mRNA_sum)

        self.setUp()

        F1 = 0
        F2 = 0
        F3 = 0
        F4 = 1
        F5 = 0
        self.splice_proc.CUMULATIVE_F1_TO_F5 = np.cumsum([F1, F2, F3, F4, F5])
        self.splice_proc.PROB_SPLICE_SINGLE_TO_MULTI = 0

        self.full_len_transcripts_nuc[:] = 100
        initial_mRNA_sum = count_total_mRNA(self.state, True)[0]

        for i in range(30):
            self.splice_proc.evolve_state(i)

        self.assertEqual(np.sum(self.full_len_transcripts_nuc) + sum_ss_type(self.state, 3), initial_mRNA_sum)

        self.setUp()

        F1 = 0
        F2 = 0
        F3 = 0
        F4 = 0
        F5 = 1
        self.splice_proc.CUMULATIVE_F1_TO_F5 = np.cumsum([F1, F2, F3, F4, F5])
        self.splice_proc.PROB_SPLICE_SINGLE_TO_MULTI = 0

        self.full_len_transcripts_nuc[:] = 100
        initial_mRNA_sum = count_total_mRNA(self.state, True)[0]

        for i in range(30):
            self.splice_proc.evolve_state(i)

        self.assertEqual(np.sum(self.full_len_transcripts_nuc) + sum_ss_type(self.state, 6), initial_mRNA_sum)

    # Do not allow for AlternativeSplicing
    # Make sure that all the mRNA bins stay the same, and that there is no shuffling within a particular level
    def test_constant_tree(self):
        self.s1_splice_proc.PROB_SPLICE_FULL_TO_SINGLE = 0
        self.s1_splice_proc.PROB_SPLICE_SINGLE_TO_MULTI = 0
        self.s1_splice_proc.PROB_VPR_THIRD_SPLICE = 0
        self.s1_splice_proc.PROB_VIF_THIRD_SPLICE = 0

        prev_s1_full_len_transcripts_nuc = np.copy(self.s1_full_len_transcripts_nuc)
        prev_s1_single_splice_transcript_nuc = np.copy(self.s1_single_splice_transcript_nuc)
        prev_s1_multi_splice_transcript_nuc = np.copy(self.s1_multi_splice_transcript_nuc)

        for i in range(30):
            self.s1_splice_proc.evolve_state(i)

        self.assertTrue(np.array_equal(prev_s1_full_len_transcripts_nuc, self.s1_full_len_transcripts_nuc))
        self.assertTrue(np.array_equal(prev_s1_single_splice_transcript_nuc, self.s1_single_splice_transcript_nuc))
        self.assertTrue(np.array_equal(prev_s1_multi_splice_transcript_nuc, self.s1_multi_splice_transcript_nuc))

    # See how changing F values affects what m-vif, m-vpr becomes
    def test_cumulative_prob_partition(self):
        self.splice_proc.PROB_VPR_THIRD_SPLICE = 1
        self.splice_proc.PROB_VIF_THIRD_SPLICE = 1
        F1 = 0
        F2 = 0
        F3 = 1
        F4 = 0
        F5 = 0
        self.splice_proc.CUMULATIVE_F1_TO_F5 = np.cumsum([F1, F2, F3, F4, F5])
        self.splice_proc.CUMULATIVE_F3_TO_F5 = np.cumsum([(F3/(F3+F4+F5)), (F4/(F3+F4+F5)), (F5/(F3+F4+F5))])

        self.multi_splice_transcript_nuc[0] = 300
        self.multi_splice_transcript_nuc[6] = 300
        initial_mRNA_sum = count_total_mRNA(self.state, True)[0]

        for i in range(5):
            self.splice_proc.evolve_state(i)

        self.assertEqual(np.sum(self.multi_splice_transcript_nuc[[1,7]]), initial_mRNA_sum)

        self.setUp()

        self.splice_proc.PROB_VPR_THIRD_SPLICE = 1
        self.splice_proc.PROB_VIF_THIRD_SPLICE = 1
        F1 = 0
        F2 = 0
        F3 = 0
        F4 = 1
        F5 = 0
        self.splice_proc.CUMULATIVE_F1_TO_F5 = np.cumsum([F1, F2, F3, F4, F5])
        self.splice_proc.CUMULATIVE_F3_TO_F5 = np.cumsum([(F3/(F3+F4+F5)), (F4/(F3+F4+F5)), (F5/(F3+F4+F5))])

        self.multi_splice_transcript_nuc[0] = 300
        self.multi_splice_transcript_nuc[6] = 300
        initial_mRNA_sum = count_total_mRNA(self.state, True)[0]

        for i in range(5):
            self.splice_proc.evolve_state(i)

        self.assertEqual(np.sum(self.multi_splice_transcript_nuc[[2,8]]), initial_mRNA_sum)

        self.setUp()

        self.splice_proc.PROB_VPR_THIRD_SPLICE = 1
        self.splice_proc.PROB_VIF_THIRD_SPLICE = 1
        F1 = 0
        F2 = 0
        F3 = 0
        F4 = 0
        F5 = 1
        self.splice_proc.CUMULATIVE_F1_TO_F5 = np.cumsum([F1, F2, F3, F4, F5])
        self.splice_proc.CUMULATIVE_F3_TO_F5 = np.cumsum([(F3/(F3+F4+F5)), (F4/(F3+F4+F5)), (F5/(F3+F4+F5))])

        self.multi_splice_transcript_nuc[0] = 300
        self.multi_splice_transcript_nuc[6] = 300
        initial_mRNA_sum = count_total_mRNA(self.state, True)[0]

        for i in range(5):
            self.splice_proc.evolve_state(i)

        self.assertEqual(np.sum(self.multi_splice_transcript_nuc[[5,11]]), initial_mRNA_sum)        

    # Start out with zero of everything  
    def test_zero_abundance(self):
        for i in range(30):
            self.splice_proc.evolve_state(i)
            if (abundance_is_nonzero(self.state)):
                self.fail('Abundances should be zero')
        self.assertEqual(1,1)
    
    # Make sure abundance values are non-negative integers at all times
    def test_non_negative_integer_abundances(self):
        for i in range(30):
            self.s1_splice_proc.evolve_state(i)
            if abundance_is_negative(self.s1_state):
                self.fail('One or more abundances is a negative value.')
            elif abundance_is_not_integer(self.s1_state):
                self.fail('One or more abundances is not an integer.')
            else:
                self.assertEquals(1,1)
        
    # Use this method if I want to change any conditions
    # immediately after running a test.
    def tearDown(self):
        pass
    
if __name__ == '__main__':
    unittest.main()

