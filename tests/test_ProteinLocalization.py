from process.ProteinLocalization import *
from mainaux.TestHelpers import *
from state.Proteins import *
import unittest
import numpy as np

class TestProteinLocalization(unittest.TestCase):
    
    def setUp(self):
        # Default state
        self.state = State()
        self.shuffling_process = ProteinLocalization(self.state)

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
        self.s1_shuffling_process = ProteinLocalization(self.s1_state)

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
    # mass balance test: Rev (C+N) before = Rev (C+N) after
    def test_Rev_balance(self):
        for i in range(30):
            prev_Rev_count = count_total_Rev(self.s1_state)
            self.s1_shuffling_process.evolve_state(i)
            curr_Rev_count = count_total_Rev(self.s1_state)
            
            self.assertEqual(prev_Rev_count, curr_Rev_count)
    
    # Tat (C+N) before = Tat (C+N) after
    def test_Tat_balance(self):
        for i in range(30):
            prev_Tat_count = count_total_Tat(self.s1_state)
            self.s1_shuffling_process.evolve_state(i)
            curr_Tat_count = count_total_Tat(self.s1_state)            

            self.assertEqual(prev_Tat_count, curr_Tat_count)
    
    def test_shuffling_rates(self):
        
        # Examine Rev shuffling rates
        # Probability of each Rev shuttling into nucleus is 1
        self.s1_shuffling_process.PROB_REV_SHUTTLING_IN = 1
        self.s1_shuffling_process.PROB_REV_SHUTTLING_OUT = 0

        expected_Rev_count_in_nuc = self.s1_proteins_nuc[Proteins.index['Rev']] + self.s1_proteins_cyt[Proteins.index['Rev']]
        self.s1_shuffling_process.evolve_state(5)
        
        #self.assertEqual(self.s1_proteins_nuc[Proteins.index['Rev']], expected_Rev_count_in_nuc)
        #self.assertEqual(self.s1_proteins_cyt[Proteins.index['Rev']], 0)
        
        self.setUp()        
        # Probability of each Rev shuttling out of nucleus is 1
        self.s1_shuffling_process.PROB_REV_SHUTTLING_IN = 0
        self.s1_shuffling_process.PROB_REV_SHUTTLING_OUT = 1
        
        expected_Rev_count_in_cyt = self.s1_proteins_nuc[Proteins.index['Rev']] + self.s1_proteins_cyt[Proteins.index['Rev']]        
        self.s1_shuffling_process.evolve_state(5)
        
        #self.assertEqual(self.s1_proteins_cyt[Proteins.index['Rev']], expected_Rev_count_in_cyt)
        #self.assertEqual(self.s1_proteins_nuc[Proteins.index['Rev']], 0)
        
        self.setUp()                
        # Probability of each Rev shuttling in/out of nucleus is 0
        self.s1_shuffling_process.PROB_REV_SHUTTLING_IN = 0
        self.s1_shuffling_process.PROB_REV_SHUTTLING_OUT = 0
        
        prev_Rev_in_nuc = self.s1_proteins_nuc[Proteins.index['Rev']]
        prev_Rev_in_cyt = self.s1_proteins_cyt[Proteins.index['Rev']]        
        self.s1_shuffling_process.evolve_state(5)
        
        self.assertEqual(self.s1_proteins_nuc[Proteins.index['Rev']], prev_Rev_in_nuc)
        self.assertEqual(self.s1_proteins_cyt[Proteins.index['Rev']], prev_Rev_in_cyt)
        
        self.setUp()        
        # Examine Tat shuffling rates
        # Probability of each Tat shuttling into nucleus is 1
        self.s1_shuffling_process.PROB_TAT_SHUTTLING_IN = 1
        self.s1_shuffling_process.PROB_TAT_SHUTTLING_OUT = 0
        
        expected_Tat_count_in_nuc = self.s1_proteins_nuc[Proteins.index['Tat']] + self.s1_proteins_cyt[Proteins.index['Tat']]        
        self.s1_shuffling_process.evolve_state(5)
        
        #self.assertEqual(self.s1_proteins_nuc[Proteins.index['Tat']], expected_Tat_count_in_nuc)
        #self.assertEqual(self.s1_proteins_cyt[Proteins.index['Tat']], 0)
        
        self.setUp()        
        # Probability of each Tat shuttling out of nucleus is 1
        self.s1_shuffling_process.PROB_TAT_SHUTTLING_IN = 0
        self.s1_shuffling_process.PROB_TAT_SHUTTLING_OUT = 1
        
        expected_Tat_count_in_cyt = self.s1_proteins_nuc[Proteins.index['Tat']] + self.s1_proteins_cyt[Proteins.index['Tat']]        
        self.s1_shuffling_process.evolve_state(5)
        
        #self.assertEqual(self.s1_proteins_cyt[Proteins.index['Tat']], expected_Tat_count_in_cyt)
        #self.assertEqual(self.s1_proteins_nuc[Proteins.index['Tat']], 0)
        
        self.setUp()        
        # Probability of each Tat shuttling in/out of nucleus is 0
        self.s1_shuffling_process.PROB_TAT_SHUTTLING_IN = 0
        self.s1_shuffling_process.PROB_TAT_SHUTTLING_OUT = 0
        
        prev_Tat_in_nuc = self.s1_proteins_nuc[Proteins.index['Tat']]
        prev_Tat_in_cyt = self.s1_proteins_cyt[Proteins.index['Tat']] 
        self.s1_shuffling_process.evolve_state(5)
        
        self.assertEqual(self.s1_proteins_nuc[Proteins.index['Tat']], prev_Tat_in_nuc)
        self.assertEqual(self.s1_proteins_cyt[Proteins.index['Tat']], prev_Tat_in_cyt)
        
    # Start out with zero of everything  
    def test_zero_abundance(self):
        for i in range(30):
            self.shuffling_process.evolve_state(i)
            if (abundance_is_nonzero(self.state)):
                self.fail('Abundances should be zero')
        self.assertEqual(1,1)
    
    # Make sure abundance values are non-negative integers at all times
    def test_non_negative_integer_abundances(self):
        # Modify S1 state to contain a low number of Rev/Tat molecules in cytoplasm
        self.s1_proteins_cyt[Proteins.index['Rev']] = 3
        self.s1_proteins_cyt[Proteins.index['Tat']] = 2
        
        for i in range(30):
            self.s1_shuffling_process.evolve_state(i)
            if abundance_is_negative(self.s1_state):
                self.fail('One or more abundances is a negative value.')
            elif abundance_is_not_integer(self.s1_state):
                self.fail('One or more abundances is not an integer.')
            else:
                self.assertEquals(1,1)
                
    def tearDown(self):
        pass
    
if __name__ == '__main__':
    unittest.main()
