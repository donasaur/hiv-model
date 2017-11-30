from process.MRNAExport import *
from mainaux.TestHelpers import *
from state.Proteins import *
import unittest
import numpy as np

class TestMRNAExport(unittest.TestCase):
    
    # This method is run prior to each test automatically.
    def setUp(self):
        # Default state
        self.state = State()
        self.mRNAexport_process = MRNAExport(self.state)

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
        self.s1_mRNAexport_process = MRNAExport(self.s1_state)

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
        
    # After mRNA export, the total # of mRNA should still be the same
    def test_mRNA_balance(self):
        for i in range(30):
            prev_mRNA_count = count_total_mRNA(self.s1_state)
            self.s1_mRNAexport_process.evolve_state(i)
            curr_mRNA_count = count_total_mRNA(self.s1_state)
            
            self.assertEqual(prev_mRNA_count, curr_mRNA_count)
    
    # After mRNA export, the full-length/single-spliced mRNA strands with
    # < MRNAExport.NUM_OF_REV_REQ_FOR_EXPORT should have the same abundance
    # in their respective locations (nucleus or cytoplasm)
    # By default, only 8-Rev bound single-spliced mRNA, full-length mRNA
    # should be exported
    def test_stationary_mRNA_export(self):
        MIN_REV_REQ = self.s1_mRNAexport_process.NUM_OF_REV_REQ_FOR_EXPORT        
        
        for i in range(30):
            prev_stationary_mRNA_in_nuc = self.s1_full_len_transcripts_nuc[0:MIN_REV_REQ].sum() + self.s1_single_splice_transcript_nuc[0:7*MIN_REV_REQ].sum()
            prev_stationary_mRNA_in_cyt = self.s1_full_len_transcripts_cyt[0:MIN_REV_REQ].sum() + self.s1_single_splice_transcript_cyt[0:7*MIN_REV_REQ].sum()
            
            self.s1_mRNAexport_process.evolve_state(i)
            
            curr_stationary_mRNA_in_nuc = self.s1_full_len_transcripts_nuc[0:MIN_REV_REQ].sum() + self.s1_single_splice_transcript_nuc[0:7*MIN_REV_REQ].sum()
            curr_stationary_mRNA_in_cyt = self.s1_full_len_transcripts_cyt[0:MIN_REV_REQ].sum() + self.s1_single_splice_transcript_cyt[0:7*MIN_REV_REQ].sum()
            
            self.assertEqual(prev_stationary_mRNA_in_nuc, curr_stationary_mRNA_in_nuc)
            self.assertEqual(prev_stationary_mRNA_in_cyt, curr_stationary_mRNA_in_cyt)

    def test_stationary_mRNA_export_deep_equals(self):
    	MIN_REV_REQ = self.s1_mRNAexport_process.NUM_OF_REV_REQ_FOR_EXPORT

    	for i in range(10):
    		prev_stationary_full_len_in_nuc = np.copy(self.s1_full_len_transcripts_nuc[0:MIN_REV_REQ])
    		prev_stationary_sing_splice_in_nuc = np.copy(self.s1_single_splice_transcript_nuc[0:7*MIN_REV_REQ])

    		self.s1_mRNAexport_process.evolve_state(i)

    		curr_stationary_full_len_in_nuc = np.copy(self.s1_full_len_transcripts_nuc[0:MIN_REV_REQ])
    		curr_stationary_sing_splice_in_nuc = np.copy(self.s1_single_splice_transcript_nuc[0:7*MIN_REV_REQ])

    		self.assertTrue(np.array_equal(prev_stationary_full_len_in_nuc,curr_stationary_full_len_in_nuc))
    		self.assertTrue(np.array_equal(prev_stationary_sing_splice_in_nuc,curr_stationary_sing_splice_in_nuc))
        
    # Check extremes of export rate constants
    # potential mRNA = mRNA strand with potential to be exported
    def test_export_rate(self):
        MIN_REV_REQ = self.s1_mRNAexport_process.NUM_OF_REV_REQ_FOR_EXPORT
        
        # Expect all mRNA strands with the potential to be exported
        # to get exported to cytoplasm
        self.s1_mRNAexport_process.PROB_REV_INDEP_EXPORT = 1
        self.s1_mRNAexport_process.PROB_REV_DEP_EXPORT = 1
        
        self.s1_mRNAexport_process.evolve_state(5)
        
        potential_mRNA_count_in_nuc = count_potential_mRNA(self.s1_state, MIN_REV_REQ)
        
        self.assertEquals(potential_mRNA_count_in_nuc, 0)
        
        self.setUp()
        
        # Export rate is zero
        self.s1_mRNAexport_process.PROB_REV_INDEP_EXPORT = 0
        self.s1_mRNAexport_process.PROB_REV_DEP_EXPORT = 0

        prev_potential_mRNA_count_in_nuc = count_potential_mRNA(self.s1_state, MIN_REV_REQ)
        self.s1_mRNAexport_process.evolve_state(5)
        curr_potential_mRNA_count_in_nuc = count_potential_mRNA(self.s1_state, MIN_REV_REQ)
        
        self.assertEquals(prev_potential_mRNA_count_in_nuc, curr_potential_mRNA_count_in_nuc)
        
        self.setUp()        
        
        self.s1_mRNAexport_process.PROB_REV_DEP_EXPORT = 0
        self.s1_mRNAexport_process.PROB_REV_INDEP_EXPORT = 1
        
        prev_pot_dependent_mRNA_in_nuc = self.s1_single_splice_transcript_nuc[7*MIN_REV_REQ:].sum() + self.s1_full_len_transcripts_nuc[MIN_REV_REQ:].sum()
        
        self.s1_mRNAexport_process.evolve_state(5)
        
        curr_pot_dependent_mRNA_in_nuc = self.s1_single_splice_transcript_nuc[7*MIN_REV_REQ:].sum() + self.s1_full_len_transcripts_nuc[MIN_REV_REQ:].sum()
        curr_pot_indep_mRNA_in_nuc = self.s1_multi_splice_transcript_nuc.sum()
        
        self.assertEquals(prev_pot_dependent_mRNA_in_nuc, curr_pot_dependent_mRNA_in_nuc)
        self.assertEquals(curr_pot_indep_mRNA_in_nuc, 0)
        
        self.setUp()        
        
        self.s1_mRNAexport_process.PROB_REV_INDEP_EXPORT = 0
        self.s1_mRNAexport_process.PROB_REV_DEP_EXPORT = 1
        
        prev_pot_indep_mRNA_in_nuc = self.s1_multi_splice_transcript_nuc.sum()

        self.s1_mRNAexport_process.evolve_state(5)
        
        curr_pot_dependent_mRNA_in_nuc = self.s1_single_splice_transcript_nuc[7*MIN_REV_REQ:].sum() + self.s1_full_len_transcripts_nuc[MIN_REV_REQ:].sum()
        curr_pot_indep_mRNA_in_nuc = self.s1_multi_splice_transcript_nuc.sum()
        
        self.assertEquals(prev_pot_indep_mRNA_in_nuc, curr_pot_indep_mRNA_in_nuc)
        self.assertEquals(curr_pot_dependent_mRNA_in_nuc, 0)
    
    # Change the minimum # of bounded Rev required for export
    def test_num_of_Rev_to_export(self):
        # All potential Rev-bound mRNAs to be exported to cytoplasm
        # will be exported at each timestep

        for i in np.arange(4,9):        
            self.setUp()
            self.s1_mRNAexport_process.NUM_OF_REV_REQ_FOR_EXPORT = i
            self.s1_mRNAexport_process.PROB_REV_DEP_EXPORT = 1
            MIN_REV_REQ = self.s1_mRNAexport_process.NUM_OF_REV_REQ_FOR_EXPORT
            
            self.s1_mRNAexport_process.evolve_state(5)
            pot_dependent_mRNA_in_nuc = self.s1_single_splice_transcript_nuc[7*MIN_REV_REQ:].sum() + self.s1_full_len_transcripts_nuc[MIN_REV_REQ:].sum()
            
            self.assertEquals(pot_dependent_mRNA_in_nuc, 0)
        
    # Start out with zero of everything  
    def test_zero_abundance(self):
        for i in range(30):
            self.mRNAexport_process.evolve_state(i)
            if (abundance_is_nonzero(self.state)):
                self.fail('Abundances should be zero')
        self.assertEqual(1,1)
    
    # Make sure abundance values are non-negative integers at all times
    def test_non_negative_integer_abundances(self):
        for i in range(30):
            self.s1_mRNAexport_process.evolve_state(i)
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

