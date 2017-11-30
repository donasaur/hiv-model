from process.RevBinding import *
from mainaux.TestHelpers import *
from state.Proteins import *
import unittest
import numpy as np

class TestRevBinding(unittest.TestCase):
    
    # This method is run prior to each test automatically.
    def setUp(self):
        # Default state
        self.state = State()
        self.r_binding_proc = RevBinding(self.state)

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
        self.s1_r_binding_proc = RevBinding(self.s1_state)

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
        self.s2_r_binding_proc = RevBinding(self.s2_state)

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
        self.s3_r_binding_proc = RevBinding(self.s3_state)

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
        self.s4_r_binding_proc = RevBinding(self.s4_state)

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
        index = Proteins.index

    # Test Rev balance
    def test_Rev_balance(self):
        for i in range(30):
            prev_Rev_count = count_total_Rev(self.s1_state)
            self.s1_r_binding_proc.evolve_state(i)
            curr_Rev_count = count_total_Rev(self.s1_state)

            self.assertEqual(prev_Rev_count, curr_Rev_count)

        for i in range(30):
            prev_Rev_count = count_total_Rev(self.s2_state)
            self.s2_r_binding_proc.evolve_state(i)
            curr_Rev_count = count_total_Rev(self.s2_state)

            self.assertEqual(prev_Rev_count, curr_Rev_count)

        for i in range(30):
            prev_Rev_count = count_total_Rev(self.s3_state)
            self.s3_r_binding_proc.evolve_state(i)
            curr_Rev_count = count_total_Rev(self.s3_state)

            self.assertEqual(prev_Rev_count, curr_Rev_count)

        for i in range(30):
            prev_Rev_count = count_total_Rev(self.s4_state)
            self.s4_r_binding_proc.evolve_state(i)
            curr_Rev_count = count_total_Rev(self.s4_state)
        
            self.assertEqual(prev_Rev_count, curr_Rev_count)

    def test_ODE_discretizer(self):
        soln = np.zeros((60, 10))
        soln_last_row = soln[-1,:]
        
        # Test 1: Check a soln with only whole numbers
        # Make sure that prev_mRNA_abundance = new_mRNA_abundance (total # of mRNA does not change)
        # Make sure that total Rev count remains constant
        # 2210 mRNA total
        soln[-1,:] = [213, 123, 238, 238, 548, 128, 238, 192, 292, 129]

        # 2210 mRNA total
        prev_mRNA_abundance = np.array([507, 283, 271, 102, 182, 288, 173, 172, 232])
        prev_free_Rev_count = count_bound_Rev(soln_last_row[0:9]) - count_bound_Rev(prev_mRNA_abundance)

        new_mRNA_abundance, new_Rev_count = self.r_binding_proc.ODE_discretizer(soln, prev_mRNA_abundance, prev_free_Rev_count)

        self.assertEqual(prev_mRNA_abundance.sum(), new_mRNA_abundance.sum())
        self.assertEqual(count_bound_Rev(prev_mRNA_abundance) + prev_free_Rev_count, count_bound_Rev(new_mRNA_abundance) + new_Rev_count)

        # Test 2: Check a soln where the bins are the same in value as the prev mRNA abundance bins + free Rev count bin
        # Make sure that each bin of new_mRNA_abundance has the same value as the corresponding bin in new_mRNA_abundance
        soln[-1,:] = [213, 123, 238, 238, 548, 128, 238, 192, 292, 129]

        prev_mRNA_abundance = np.array([213, 123, 238, 238, 548, 128, 238, 192, 292])
        prev_free_Rev_count = 129

        new_mRNA_abundance, new_Rev_count = self.r_binding_proc.ODE_discretizer(soln, prev_mRNA_abundance, prev_free_Rev_count)

        self.assertTrue(np.array_equal(prev_mRNA_abundance, new_mRNA_abundance))
        self.assertEqual(prev_free_Rev_count, new_Rev_count)

        # Test 3: Check a soln where all the bins are (-) except one
        # Make sure that the solution generated has all bins have a value of zero except one
        # Make sure that the remaining bin has a value that is equal to the sum of all the Rev in the nucleus
        soln[-1,:] = [40, -10, 0, 0, 0, 0, 0, 0, 0, 100]
        
        mRNA_abund_intermediate = np.array([10, 0, 0, 0, 0, 0, 0, 0, 0])
        
        prev_mRNA_abundance = np.array([30, 0, 0, 0, 0, 0, 0, 0, 0])
        prev_free_Rev_count = 90
        
        new_mRNA_abundance, new_Rev_count = self.r_binding_proc.ODE_discretizer(soln, prev_mRNA_abundance, prev_free_Rev_count)
        self.assertTrue(np.array_equal(new_mRNA_abundance, np.array([30, 0, 0, 0, 0, 0, 0, 0, 0])))
        self.assertEqual(new_Rev_count, 90)

        # Test 4: Check the solution where all bins are a decimal except one
        # Expect Rounding from ODE_discretizer to inflate the mRNA count by two
        # and so the bin with mRNA count greater than 2 will have 2 mRNA removed
        # And then since there is more Rev in the system than the previous timestep
        # Expect to remove that Rev from the amount of free Rev in the system
        soln[-1, :] = [15, 0, 0.5, 0, 0.5, 0, 0, 0, 0, 20]

        prev_mRNA_abundance = np.array([16, 0, 0, 0, 0, 0, 0, 0, 0])
        prev_free_Rev_count = 23

        new_mRNA_abundance, new_Rev_count = self.r_binding_proc.ODE_discretizer(soln, prev_mRNA_abundance, prev_free_Rev_count)

        self.assertTrue(np.array_equal(new_mRNA_abundance, np.array([14, 0, 1, 0, 1, 0, 0, 0, 0])) or np.array_equal(new_mRNA_abundance, np.array([15, 0, 0, 0, 1, 0, 0, 0, 0])) or np.array_equal(new_mRNA_abundance, np.array([15, 0, 1, 0, 0, 0, 0, 0, 0])))
        self.assertTrue(new_Rev_count == 17 or new_Rev_count == 19 or new_Rev_count == 21)
        

    # Start out with zero of everything  
    def test_zero_abundance(self):
        for i in range(30):
            self.r_binding_proc.evolve_state(i)
            if (abundance_is_nonzero(self.state)):
                self.fail('Abundances should be zero')
        self.assertEqual(1,1)
    
    # Make sure abundance values are non-negative integers at all times
    def test_non_negative_integer_abundances(self):
        for i in range(30):
            self.s1_r_binding_proc.evolve_state(i)
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

