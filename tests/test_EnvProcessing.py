from process.EnvProcessing import *
from mainaux.TestHelpers import *
from state.Proteins import *
import unittest
import numpy as np

class TestEnvProcessing(unittest.TestCase):
    TIMESTEP = 2300
    
    # This method is run prior to each test automatically.
    def setUp(self):
        # Default state
        self.state = State()
        self.env_proc = EnvProcessing(self.state)

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
        
        # sx state
        if self.TIMESTEP == 2300:
            self.sx_state = s7_state()
        elif self.TIMESTEP == 1500:
            self.sx_state = s8_state()
        self.sx_env_proc = EnvProcessing(self.sx_state)

        mRNA_state = self.sx_state.get_state('mRNAs')
        protein_state = self.sx_state.get_state('proteins')

        self.sx_proteins_nuc = protein_state.proteins_nuc # proteins_nuc[Proteins.index['Rev']] is Rev count in nucleus
        self.sx_proteins_cyt = protein_state.proteins_cyt # proteins_cyt[Proteins.index['Rev']] is Rev count in cytoplasm
        self.sx_proteins_mem = protein_state.proteins_mem
        self.sx_env_misc = protein_state.env_misc

        self.sx_proteins_cyt[Proteins.index['Vif']] = 200
        self.sx_proteins_cyt[Proteins.index['Gag']] = 200
        self.sx_proteins_mem[Proteins.index['Gag']] = 500
        self.sx_proteins_cyt[Proteins.index['Gag_dimers']] = 200
        self.sx_proteins_mem[Proteins.index['Gag_dimers']] = 500       

        self.sx_container = self.sx_state.get_state('viral_progeny_container')

        self.sx_full_len_transcripts_nuc = mRNA_state.full_len_transcripts_nuc
        self.sx_full_len_transcripts_cyt = mRNA_state.full_len_transcripts_cyt
        self.sx_single_splice_transcript_nuc = mRNA_state.single_splice_transcript_nuc
        self.sx_single_splice_transcript_cyt = mRNA_state.single_splice_transcript_cyt
        self.sx_multi_splice_transcript_nuc = mRNA_state.multi_splice_transcript_nuc
        self.sx_multi_splice_transcript_cyt = mRNA_state.multi_splice_transcript_cyt

        self.sx_full_len_transcripts_Gag_bound = mRNA_state.full_len_transcripts_Gag_bound
        self.sx_full_len_transcripts_Gag_bound = self.sx_full_len_transcripts_Gag_bound.reshape((2,2,2,2))
        self.sx_full_len_transcripts_Gag_bound_no_reshape = mRNA_state.full_len_transcripts_Gag_bound
        
        # Other setup stuff
        self.index = Proteins.index

    def test_active_Env_balance(self):

        prev_active_Env_count = count_active_Env(self.sx_state)
        for i in range(200):
            self.sx_env_proc.evolve_state(i)
        curr_active_Env_count = count_active_Env(self.sx_state)

        self.assertEqual(prev_active_Env_count, curr_active_Env_count)

    def test_inactive_Env_balance(self):
        prev_inactive_Env_count = count_inactive_Env(self.sx_state)
        for i in range(200):
            self.sx_env_proc.evolve_state(i)
        curr_inactive_Env_count = count_inactive_Env(self.sx_state)

        self.assertEqual(prev_inactive_Env_count, curr_inactive_Env_count)

    # Run evolve_state up to steps 1 - 9 being the final step
    def test_Env_balance_when_stopping_after_certain_steps(self):
        num_times = 1
        for i in range(9):
            final_step = i + 1
            for j in range(num_times): # go up to final_step num_times
                self.setUp()
                prev_Env_count = count_total_Env(self.sx_state)
                self.sx_env_proc.evolve_state(30 + j, final_step)
                curr_Env_count = count_total_Env(self.sx_state)
                self.assertEqual(prev_Env_count, curr_Env_count)
        
    def test_transfer_of_Env_to_ER(self):
        self.sx_env_proc.rate_of_ER_localization = 1
        prev_env_cyt_count = self.sx_env_misc['Env : cytoplasm']
        self.sx_env_proc.evolve_state(self.TIMESTEP, single_step=1)
        curr_env_cyt_count = self.sx_env_misc['Env : cytoplasm']
        self.assertEqual(curr_env_cyt_count, 0)
        self.assertEqual(self.sx_env_misc['Env : ER'], prev_env_cyt_count)

    def test_transfer_to_G1_bin(self):
        self.sx_env_proc.rate_oligosaccharyltransferase = 0
        prev_env_ER_count = self.sx_env_misc['Env : ER']
        prev_env_G1_count = self.sx_env_misc['Env : ER : G1']
        self.sx_env_proc.evolve_state(self.TIMESTEP, single_step=2)
        curr_env_ER_count = self.sx_env_misc['Env : ER']
        curr_env_G1_count = self.sx_env_misc['Env : ER : G1']
        self.assertEqual(curr_env_ER_count, prev_env_ER_count)
        self.assertEqual(curr_env_G1_count, prev_env_G1_count)

    def test_transfer_to_G2_bin(self):
        self.sx_env_proc.rate_glucosidaseI = 0
        prev_env_G1_count = self.sx_env_misc['Env : ER : G1']
        prev_env_G2_count = self.sx_env_misc['Env : ER : G2']
        self.sx_env_proc.evolve_state(self.TIMESTEP, single_step=3)
        curr_env_G1_count = self.sx_env_misc['Env : ER : G1']
        curr_env_G2_count = self.sx_env_misc['Env : ER : G2']
        self.assertEqual(curr_env_G1_count, prev_env_G1_count)
        self.assertEqual(curr_env_G2_count, prev_env_G2_count)

    # expect about a 0% error rate
    def test_G5_error_rate(self):
        self.sx_env_proc.prob_golgi_glycosylation_error = 0
        prev_G5_count = self.sx_env_misc['Env : Golgi : G5']
        prev_G5_error_count = self.sx_env_misc['Env : Golgi : G5 : error']
        self.sx_env_proc.evolve_state(self.TIMESTEP, single_step=6)
        curr_G5_count = self.sx_env_misc['Env : Golgi : G5']
        curr_G5_error_count = self.sx_env_misc['Env : Golgi : G5 : error']
        self.assertEqual(prev_G5_count, curr_G5_count)
        self.assertEqual(prev_G5_error_count, curr_G5_error_count)

    def test_trimerization(self):
        self.sx_env_proc.evolve_state(self.TIMESTEP, single_step=7)
        total_trimers = sum(self.sx_env_misc['Env : trimers'])
        self.assertAlmostEqual(float(self.sx_env_misc['Env : trimers'][0])/total_trimers, float(1)/8, places=1)
        self.assertAlmostEqual(float(self.sx_env_misc['Env : trimers'][1])/total_trimers, float(3)/8, places=1)
        self.assertAlmostEqual(float(self.sx_env_misc['Env : trimers'][2])/total_trimers, float(3)/8, places=1)
        self.assertAlmostEqual(float(self.sx_env_misc['Env : trimers'][3])/total_trimers, float(1)/8, places=1)

    def test_cleavage_ratio(self):
        self.sx_env_proc.rate_cleavage = 0
        prev_trimers_count = sum(self.sx_env_misc['Env : trimers'])
        prev_cleaved_count = sum(self.sx_env_misc['Env : trimers : cleaved'])
        self.sx_env_proc.evolve_state(self.TIMESTEP, single_step=8)
        curr_trimers_count = sum(self.sx_env_misc['Env : trimers'])
        curr_cleaved_count = sum(self.sx_env_misc['Env : trimers : cleaved'])
        self.assertEqual(prev_trimers_count, curr_trimers_count)
        self.assertEqual(prev_cleaved_count, curr_cleaved_count)

    def test_membrane_localization(self):
        self.sx_env_proc.rate_membrane_localization = 0
        prev_cleaved_count = sum(self.sx_env_misc['Env : trimers : cleaved'])
        prev_membrane_count = sum(self.sx_env_misc['Env : trimers : membrane'])
        self.sx_env_proc.evolve_state(self.TIMESTEP, single_step=9)
        curr_cleaved_count = sum(self.sx_env_misc['Env : trimers : cleaved'])
        curr_membrane_count = sum(self.sx_env_misc['Env : trimers : membrane'])
        self.assertEqual(prev_cleaved_count, curr_cleaved_count)
        self.assertEqual(prev_membrane_count, curr_membrane_count)

    def test_virons_incorporation(self):
        self.sx_env_proc.evolve_state(self.TIMESTEP, single_step=10)
        total_trimers = self.sx_container.count_total_num_of_Env_t()
        self.assertAlmostEqual(float(self.sx_container.count_num_of_Env_t(0))/total_trimers, float(1)/8, places=1)
        self.assertAlmostEqual(float(self.sx_container.count_num_of_Env_t(1))/total_trimers, float(3)/8, places=1)
        self.assertAlmostEqual(float(self.sx_container.count_num_of_Env_t(2))/total_trimers, float(3)/8, places=1)
        self.assertAlmostEqual(float(self.sx_container.count_num_of_Env_t(3))/total_trimers, float(1)/8, places=1)

        self.setUp()
        self.sx_env_proc.rate_viron_incorporation = 0
        prev_total_trimers = self.sx_container.count_total_num_of_Env_t()
        self.sx_env_proc.evolve_state(self.TIMESTEP, single_step=10)
        curr_total_trimers = self.sx_container.count_total_num_of_Env_t()
        self.assertEqual(prev_total_trimers, curr_total_trimers)

    # Start out with zero of everything  
    def test_zero_abundance(self):
        for i in range(30):
            self.env_proc.evolve_state(i)
            if (abundance_is_nonzero(self.state)):
                self.fail('Abundances should be zero')
        self.assertEqual(1,1)
    
    # Make sure abundance values are non-negative integers at all times
    def test_non_negative_integer_abundances(self):
        for i in range(30):
            self.sx_env_proc.evolve_state(i)
            if abundance_is_negative(self.sx_state):
                self.fail('One or more abundances is a negative value.')
            elif abundance_is_not_integer(self.sx_state):
                self.fail('One or more abundances is not an integer.')
            else:
                self.assertEquals(1,1)        

    # Use this method if I want to change any conditions
    # immediately after running a test.
    def tearDown(self):
        pass

class TestEnvProcessing1500(TestEnvProcessing):
    TIMESTEP = 1500       
    
if __name__ == '__main__':
    unittest.main()

