from process.Translation import *
from mainaux.TestHelpers import *
from state.Proteins import *
import unittest
import numpy as np

class TestTranslation(unittest.TestCase):
    
    # This method is run prior to each test automatically.
    def setUp(self):
        # Default state
        self.state = State()
        self.trans_process = Translation(self.state)

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
        self.s1_trans_process = Translation(self.s1_state)

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

        # Other setup stuff
        self.index = Proteins.index
    
    # Makes sure that all the proteins of interest are being created during translation step
    # Test relies on large abundances for each mRNA
    def test_protein_creation(self):
        self.s1_proteins_cyt[self.s1_proteins_cyt>=0] = 0
        for i in range(30):
            self.s1_trans_process.evolve_state(i)
            self.assertFalse(self.s1_proteins_cyt.any == 0)

    # See that different proteins are translated from the correct respective mRNAs
    # Try zeroing the mRNAs that are translated to create certain proteins
    def test_Env_creation(self):
        self.s1_proteins_cyt[self.s1_proteins_cyt>=0] = 0
        for j in range(7):
            for i in np.arange(0+j, (7*self.s1_trans_process.MAX_REV_PER_TRANSCRIPT)+1+j, 7):
                if j >= 3:
                    self.s1_single_splice_transcript_cyt[i] = 0
        for i in range(30):
            self.s1_trans_process.evolve_state(i)
            Env_count_in_cyt = self.s1_proteins_cyt[self.index['Env']]
            self.assertEqual(Env_count_in_cyt, 0)
    
    def test_Vif_creation(self):
        self.s1_proteins_cyt[self.s1_proteins_cyt>=0] = 0
        self.s1_multi_splice_transcript_cyt[0] = 0
        for i in np.arange(0, 7*self.s1_trans_process.MAX_REV_PER_TRANSCRIPT+1, 7):
            self.s1_single_splice_transcript_cyt[i] = 0
        for i in range(30):
            self.s1_trans_process.evolve_state(i)
            Vif_count_in_cyt = self.s1_proteins_cyt[self.index['Vif']]
            self.assertEqual(Vif_count_in_cyt, 0)
        
    def test_Vpr_creation(self):
        self.s1_proteins_cyt[self.s1_proteins_cyt>=0] = 0
        self.s1_multi_splice_transcript_cyt[6] = 0
        for i in np.arange(0+1, 7*self.s1_trans_process.MAX_REV_PER_TRANSCRIPT+1+1, 7):
            self.s1_single_splice_transcript_cyt[i] = 0
        for i in range(30):
            self.s1_trans_process.evolve_state(i)
            Vpr_count_in_cyt = self.s1_proteins_cyt[self.index['Vpr']]
            self.assertEqual(Vpr_count_in_cyt, 0)

        
    def test_Tat_creation(self):
        self.s1_proteins_cyt[self.s1_proteins_cyt>=0] = 0
        self.s1_multi_splice_transcript_cyt[1] = 0
        self.s1_multi_splice_transcript_cyt[7] = 0
        self.s1_multi_splice_transcript_cyt[12] = 0
        for i in np.arange(0+2, 7*self.s1_trans_process.MAX_REV_PER_TRANSCRIPT+1+2, 7):
            self.s1_single_splice_transcript_cyt[i] = 0
        for i in range(30):
            self.s1_trans_process.evolve_state(i)
            Tat_count_in_cyt = self.s1_proteins_cyt[self.index['Tat']]
            self.assertEqual(Tat_count_in_cyt, 0)
        
    def test_Rev_creation(self):
        self.s1_proteins_cyt[self.s1_proteins_cyt>=0] = 0
        self.s1_multi_splice_transcript_cyt[2] = 0
        self.s1_multi_splice_transcript_cyt[3] = 0
        self.s1_multi_splice_transcript_cyt[4] = 0
        self.s1_multi_splice_transcript_cyt[8] = 0
        self.s1_multi_splice_transcript_cyt[9] = 0
        self.s1_multi_splice_transcript_cyt[10] = 0
        self.s1_multi_splice_transcript_cyt[13] = 0
        self.s1_multi_splice_transcript_cyt[14] = 0
        self.s1_multi_splice_transcript_cyt[15] = 0
        for i in range(30):
            self.s1_trans_process.evolve_state(i)
            Rev_count_in_cyt = self.s1_proteins_cyt[self.index['Rev']]
            self.assertEqual(Rev_count_in_cyt, 0)      
        
    def test_Nef_creation(self):
        self.s1_proteins_cyt[self.s1_proteins_cyt>=0] = 0
        self.s1_multi_splice_transcript_cyt[5] = 0
        self.s1_multi_splice_transcript_cyt[11] = 0
        self.s1_multi_splice_transcript_cyt[16] = 0
        for i in range(30):
            self.s1_trans_process.evolve_state(i)
            Nef_count_in_cyt = self.s1_proteins_cyt[self.index['Nef']]
            self.assertEqual(Nef_count_in_cyt, 0)        
        
    def test_GagProPol_creation(self):
        self.s1_proteins_cyt[self.s1_proteins_cyt>=0] = 0
        self.s1_full_len_transcripts_cyt[self.s1_full_len_transcripts_cyt>=0] = 0
        for i in range(30):
            self.s1_trans_process.evolve_state(i)
            GagProPol_count_in_cyt = self.s1_proteins_cyt[self.index['GagProPol']]
            self.assertEqual(GagProPol_count_in_cyt, 0)  
        
    # Vary the translation rate of mRNA
    def test_translation_rate(self):
        self.s1_trans_process.FREQ_TRANSLATION = 0
        self.s1_proteins_cyt[self.s1_proteins_cyt>=0] = 0
        for i in range(30):
            self.s1_trans_process.evolve_state(i)
            self.assertTrue(self.s1_proteins_cyt.sum() == 0)  
    
    # Start out with zero of everything  
    def test_zero_abundance(self):
        for i in range(30):
            self.trans_process.evolve_state(i)
            if (abundance_is_nonzero(self.state)):
                self.fail('Abundances should be zero')
        self.assertEqual(1,1)
    
    # Make sure abundance values are non-negative integers at all times
    def test_non_negative_integer_abundances(self):
        for i in range(30):
            self.s1_trans_process.evolve_state(i)
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

