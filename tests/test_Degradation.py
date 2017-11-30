from process.Degradation import *
from mainaux.TestHelpers import *
from state.Proteins import *
import unittest
import numpy as np

class TestDegradation(unittest.TestCase):
    
    # This method is run prior to each test automatically.
    def setUp(self):
        """
        Old setup
        # Protein abundance arrays
        self.proteins_nuc = np.zeros(7)
        self.proteins_cyt = np.zeros(7)

        # mRNA abundance arrays
        self.full_len_transcript_nuc = np.zeros(9)
        self.full_len_transcript_cyt = np.zeros(9)
        self.single_splice_transcript_nuc = np.zeros(63)
        self.single_splice_transcript_cyt = np.zeros(63)
        self.multi_splice_transcript_nuc = np.zeros(17)
        self.multi_splice_transcript_cyt = np.zeros(17)
        """
        # Default state
        self.state = State()
        self.deg_process = Degradation(self.state)

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
        self.s1_deg_process = Degradation(self.s1_state)

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

    # Test Rev balance
    def test_Rev_balance(self):
        # Assume no Rev degradation
        self.s1_deg_process.PROB_PROTEIN_DEG_NUC = 0
        self.s1_deg_process.PROB_PROTEIN_DEG_CYT = 0

        prev_Rev_count = count_total_Rev(self.s1_state)
        self.s1_deg_process.evolve_state(5)
        curr_Rev_count = count_total_Rev(self.s1_state)
        
        self.assertEqual(prev_Rev_count, curr_Rev_count)        
    
    # Set degradation constants equal to 1
    # See if product # less than or equal to starting #
    def test_protein_degradation(self):
        # All proteins in the nucleus set to degrade (except the Rev bound to mRNA)
        self.s1_deg_process.PROB_PROTEIN_DEG_NUC = 1
        self.s1_deg_process.PROB_PROTEIN_DEG_CYT = 0
        self.s1_deg_process.evolve_state(5)
        self.assertTrue(self.s1_proteins_nuc.sum() == 0)
        self.assertTrue(self.s1_proteins_cyt.sum() != 0)
        
        self.setUp()
        
        # All proteins in the cytoplasm set to degrade (except the Rev bound to mRNA)
        self.s1_deg_process.PROB_PROTEIN_DEG_NUC = 0
        self.s1_deg_process.PROB_PROTEIN_DEG_CYT = 1
        self.s1_deg_process.evolve_state(5)
        self.assertTrue(self.s1_proteins_nuc.sum() != 0)
        self.assertTrue(self.s1_proteins_cyt.sum() == 0)
        
        self.setUp()        
        
        # All proteins in both the nucleus and the cytoplasm set to degrade
        self.s1_deg_process.PROB_PROTEIN_DEG_NUC = 0
        self.s1_deg_process.PROB_PROTEIN_DEG_CYT = 0
        self.s1_deg_process.PROB_mRNA_DEG = 0 # Otherwise, Rev will get released
        prev_protein_count = count_total_protein(self.s1_state)
        self.s1_deg_process.evolve_state(5)
        curr_protein_count = count_total_protein(self.s1_state)

        # Assuming no one protein's abundance increased during degradation step
        self.assertEqual(prev_protein_count, curr_protein_count)

    
    # Change around degradation constants of mRNAs
    def test_mRNA_degradation(self):
        
        # mRNA degradation rate set to 1
        self.s1_deg_process.PROB_mRNA_DEG = 1

        self.s1_deg_process.evolve_state(5)
        mRNA_count = count_total_mRNA(self.s1_state)
        self.assertTrue(mRNA_count == 0)
        
        self.setUp()
        
        # mRNA degradation rate set to 0
        self.s1_deg_process.PROB_mRNA_DEG = 0

        prev_mRNA_count = count_total_mRNA(self.s1_state)        
        self.s1_deg_process.evolve_state(5)
        curr_mRNA_count = count_total_mRNA(self.s1_state)
        self.assertEqual(prev_mRNA_count, curr_mRNA_count)
        
    
    # Start out with zero of everything  
    def test_zero_abundance(self):
        for i in range(30):
            self.deg_process.evolve_state(i)
            if (abundance_is_nonzero(self.state)):
                self.fail('Abundances should be zero')
        self.assertEqual(1,1)
    
    # Make sure abundance values are non-negative integers at all times
    def test_non_negative_integer_abundances(self):
        for i in range(30):
            self.s1_deg_process.evolve_state(i)
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
