from process.Packaging import *
from mainaux.TestHelpers import *
from state.Proteins import *
from state.ViralProgeny import *
import unittest
import numpy as np

class TestPackaging(unittest.TestCase):
    TIMESTEP = 1500
    
    # This method is run prior to each test automatically.
    def setUp(self):
        # Default state
        self.state = State()
        self.pack_proc = Packaging(self.state)

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
        if self.TIMESTEP == 1500:
            self.sx_state = s9_state()
        elif self.TIMESTEP == 1200:
            self.sx_state = s10_state()                       
        self.sx_pack_proc = Packaging(self.sx_state)

        mRNA_state = self.sx_state.get_state('mRNAs')
        protein_state = self.sx_state.get_state('proteins')

        self.sx_proteins_nuc = protein_state.proteins_nuc # proteins_nuc[Proteins.index['Rev']] is Rev count in nucleus
        self.sx_proteins_cyt = protein_state.proteins_cyt # proteins_cyt[Proteins.index['Rev']] is Rev count in cytoplasm
        self.sx_proteins_mem = protein_state.proteins_mem
        self.sx_proteins_virion = protein_state.proteins_virion

        self.sx_proteins_cyt[Proteins.index['Vif']] = 2000     

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

    def test_total_Gag_count(self):
        prev_total_Gag_count = count_total_Gag(self.sx_state)
        for i in range(100):
            self.sx_proteins_cyt[6] += 1000
            prev_total_Gag_count += 1000
            self.sx_proteins_cyt[Proteins.index['Vif']] = 2000
            self.sx_pack_proc.evolve_state(self.TIMESTEP+i)
        curr_total_Gag_count = count_total_Gag(self.sx_state)
        self.assertEqual(prev_total_Gag_count, curr_total_Gag_count)

    def test_total_Vif_count(self):
        prev_total_Vif_count = count_total_Vif(self.sx_state)
        for i in range(100):
            self.sx_pack_proc.evolve_state(self.TIMESTEP+i)
        curr_total_Vif_count = count_total_Vif(self.sx_state)
        self.assertEqual(prev_total_Vif_count, curr_total_Vif_count)        

    def test_active_gag_mRNA(self):
        prev_active_mRNA_count = count_active_gag_mRNA(self.sx_state)
        for i in range(400):
            self.sx_pack_proc.evolve_state(i)
        curr_active_mRNA_count = count_active_gag_mRNA(self.sx_state)
        self.assertEqual(prev_active_mRNA_count, curr_active_mRNA_count)

    def test_total_gag_MRNA(self):
        prev_total_mRNA_count = count_total_gag_mRNA(self.sx_state)
        for i in range(400):
            self.sx_pack_proc.evolve_state(i)
        curr_total_mRNA_count = count_total_gag_mRNA(self.sx_state)
        self.assertEqual(prev_total_mRNA_count, curr_total_mRNA_count)

    def test_total_Rev_count(self):
        prev_total_Rev_count = count_total_Rev(self.sx_state)
        for i in range(400):
            self.sx_pack_proc.evolve_state(i)
        curr_total_Rev_count = count_total_Rev(self.sx_state)
        self.assertEqual(prev_total_Rev_count, curr_total_Rev_count)

    def test_protein_binding_rate(self):
        prev_progeny_Vif_count = self.sx_container.count_num_of_Vif()
        prev_progeny_GagProPol_count = self.sx_container.count_num_of_GagProPol()
        prev_progeny_Vpr_count = self.sx_container.count_num_of_Vpr()
        prev_progeny_Nef_count = self.sx_container.count_num_of_Nef()    

        self.sx_container.AVE_VIF_PER_VIRON = 0
        self.sx_container.AVE_GAGPROPOL_PER_VIRON = 0
        self.sx_container.AVE_VPR_PER_VIRON = 0
        self.sx_container.AVE_NEF_PER_VIRON = 0    
        for i in range(400):
            self.sx_pack_proc.evolve_state(self.TIMESTEP+i, "virion_growth")
            self.assertEqual(prev_progeny_Vif_count, self.sx_container.count_num_of_Vif())
            self.assertEqual(prev_progeny_GagProPol_count, self.sx_container.count_num_of_GagProPol())
            self.assertEqual(prev_progeny_Vpr_count, self.sx_container.count_num_of_Vpr())
            self.assertEqual(prev_progeny_Nef_count, self.sx_container.count_num_of_Nef())

        self.setUp()
        self.sx_container.PROB_RNA_NUCLEATE_TRANSLOCATION = 1
        self.sx_container.AVE_VIF_PER_VIRON = 10000
        self.sx_container.AVE_GAGPROPOL_PER_VIRON = 10000
        self.sx_container.AVE_VPR_PER_VIRON = 10000
        self.sx_container.AVE_NEF_PER_VIRON = 10000             
        self.sx_pack_proc.evolve_state(self.TIMESTEP+1, "virion_growth")
        if [progeny for progeny in self.sx_container.list_of_progeny if progeny.state == ViralProgeny.NUCLEATE_MEM or progeny.state == ViralProgeny.GROWING_VIRION or progeny.state == ViralProgeny.NUCLEATE_CYT]:
            if self.sx_proteins_cyt[Proteins.index['Vif']] and self.sx_proteins_cyt[Proteins.index['Gag']] and self.sx_proteins_cyt[Proteins.index['Gag_dimers']]:
                self.assertNotEqual(prev_progeny_Vif_count, self.sx_container.count_num_of_Vif())
            if self.sx_proteins_cyt[Proteins.index['GagProPol']] and self.sx_proteins_cyt[Proteins.index['Gag']] and self.sx_proteins_cyt[Proteins.index['Gag_dimers']]:
                self.assertNotEqual(prev_progeny_GagProPol_count, self.sx_container.count_num_of_GagProPol())
            if self.sx_proteins_cyt[Proteins.index['Vpr']] and self.sx_proteins_cyt[Proteins.index['Gag']] and self.sx_proteins_cyt[Proteins.index['Gag_dimers']]:
                self.assertNotEqual(prev_progeny_Vpr_count, self.sx_container.count_num_of_Vpr())
            if self.sx_proteins_cyt[Proteins.index['Nef']] and self.sx_proteins_cyt[Proteins.index['Gag']] and self.sx_proteins_cyt[Proteins.index['Gag_dimers']]:                
                self.assertNotEqual(prev_progeny_Nef_count, self.sx_container.count_num_of_Nef())

    def test_gagnc_diss_rate(self):
        # If dissociation rate is 0, then expect SL bound counts to never decrease
        self.sx_pack_proc.GAGNC_DISS_RATE = 0
        self.sx_proteins_cyt[Proteins.index['Gag']] = 0
        self.sx_pack_proc.BINDING_CONSTANT_SL1 = 1
        self.sx_pack_proc.BINDING_CONSTANT_SL2 = 1
        self.sx_pack_proc.BINDING_CONSTANT_SL3 = 1
        self.sx_pack_proc.BINDING_CONSTANT_SL4 = 1

        prev_num_of_gag_bound_mRNA = self.sx_full_len_transcripts_Gag_bound[(1,1,1,1)]
        for i in range(50):
            self.sx_pack_proc.evolve_state(self.TIMESTEP + i, "cell_cycle_arrest")

        # Consider last bin only
        curr_num_of_gag_bound_mRNA = self.sx_full_len_transcripts_Gag_bound[(1,1,1,1)]
        self.assertEqual(prev_num_of_gag_bound_mRNA, curr_num_of_gag_bound_mRNA)

        self.setUp()

        # Make all the Gag-bound mRNA lose their attached Gag!
        self.sx_pack_proc.GAGNC_DISS_RATE = 1
        self.sx_proteins_cyt[Proteins.index['Gag']] = 0
        self.sx_pack_proc.BINDING_CONSTANT_SL1 = 0
        self.sx_pack_proc.BINDING_CONSTANT_SL2 = 0
        self.sx_pack_proc.BINDING_CONSTANT_SL3 = 0
        self.sx_pack_proc.BINDING_CONSTANT_SL4 = 0

        prev_num_of_gag_bound_mRNA = np.sum(self.sx_full_len_transcripts_Gag_bound_no_reshape)   
        for i in range(50):
            self.sx_pack_proc.evolve_state(self.TIMESTEP + i, "cell_cycle_arrest")

        curr_num_of_gag_bound_mRNA = np.sum(self.sx_full_len_transcripts_Gag_bound_no_reshape)                    
        self.assertEqual(curr_num_of_gag_bound_mRNA, 0)

    # self.GAG_DIFFUSION_PROB
    def test_gag_diffusion(self):
        self.sx_pack_proc.GAG_DIFFUSION_PROB = 0
        # Expect proteins_cyt bin to stay the same
        prev_cyt_count = self.sx_proteins_cyt[Proteins.index['Gag']]
        for i in range(50):
            self.sx_pack_proc.evolve_state(self.TIMESTEP + i, ["Gag_and_Gag_dimer_diffusion"])
        curr_cyt_count = self.sx_proteins_cyt[Proteins.index['Gag']]
        self.assertEqual(prev_cyt_count, curr_cyt_count)

        self.setUp()

        # Expect proteins_cyt bin to zero out        
        self.sx_pack_proc.GAG_DIFFUSION_PROB = 1
        prev_cyt_count = self.sx_proteins_cyt[Proteins.index['Gag']]
        for i in range(50):
            self.sx_pack_proc.evolve_state(self.TIMESTEP + i, ["Gag_and_Gag_dimer_diffusion"])
        curr_cyt_count = self.sx_proteins_cyt[Proteins.index['Gag']]
        self.assertEqual(curr_cyt_count, 0)        

    # self.GAG_DIMER_DIFFUSION_PROB
    def test_gag_dimer_diffusion(self):
        self.sx_pack_proc.GAG_DIMER_DIFFUSION_PROB = 0
        prev_cyt_count = self.sx_proteins_cyt[Proteins.index['Gag_dimers']]
        for i in range(50):
            self.sx_pack_proc.evolve_state(self.TIMESTEP + i, ["Gag_and_Gag_dimer_diffusion"])
        curr_cyt_count = self.sx_proteins_cyt[Proteins.index['Gag_dimers']]
        self.assertEqual(prev_cyt_count, curr_cyt_count)

        self.setUp()

        self.sx_pack_proc.GAG_DIMER_DIFFUSION_PROB = 1
        prev_cyt_count = self.sx_proteins_cyt[Proteins.index['Gag_dimers']]
        for i in range(50):
            self.sx_pack_proc.evolve_state(self.TIMESTEP + i, ["Gag_and_Gag_dimer_diffusion"])
        curr_cyt_count = self.sx_proteins_cyt[Proteins.index['Gag_dimers']]
        self.assertEqual(curr_cyt_count, 0)

    # self.PROB_GAG_BOUND_RNA_DIMERS
    def test_dimerize_mRNA_transcripts(self):
        # No transcripts will dimerize
        self.sx_pack_proc.PROB_GAG_BOUND_RNA_DIMERS = 0
        old_array = np.copy(self.sx_full_len_transcripts_Gag_bound)
        for i in range(50):
            self.sx_pack_proc.evolve_state(self.TIMESTEP + i, ["dimerize_gagbound_transcripts"])
        self.assertTrue(np.array_equal(old_array, self.sx_full_len_transcripts_Gag_bound))

        self.setUp()
        self.sx_pack_proc.PROB_GAG_BOUND_RNA_DIMERS = 1
        for i in range(50):
            self.sx_pack_proc.evolve_state(self.TIMESTEP + i, ["dimerize_gagbound_transcripts"])
            for j in range(2):
                bin_of_interest = self.sx_full_len_transcripts_Gag_bound[(1,1,1,j)]
                self.assertTrue(bin_of_interest == 0 or bin_of_interest == 1)


    # self.PROB_RNA_NUCLEATE_TRANSLOCATION
    def test_nucleate_mem_update(self):
        self.sx_container.PROB_RNA_NUCLEATE_TRANSLOCATION = 0
        self.sx_container.NUCLEATE_DISS_RATE = 0
        prev_count = self.sx_container.progeny_state_count[ViralProgeny.NUCLEATE_CYT]
        for i in range(50):
            self.sx_pack_proc.evolve_state(self.TIMESTEP + i, ["virion_growth"])
        curr_count = self.sx_container.progeny_state_count[ViralProgeny.NUCLEATE_CYT]
        self.assertEqual(prev_count, curr_count)

        self.setUp()

        self.sx_container.PROB_RNA_NUCLEATE_TRANSLOCATION = 1
        self.sx_container.NUCLEATE_DISS_RATE = 0
        for i in range(50):
            self.sx_pack_proc.evolve_state(self.TIMESTEP + i, ["virion_growth"])
        curr_count = self.sx_container.progeny_state_count[ViralProgeny.NUCLEATE_CYT]
        self.assertEqual(curr_count, 0)

    # self.NUCLEATE_DISS_RATE
    def test_nucleate_mem_diss(self):
        self.sx_container.NUCLEATE_DISS_RATE = 0
        prev_count = self.sx_container.progeny_state_count[ViralProgeny.NUCLEATE_CYT]
        for i in range(50):
            self.sx_pack_proc.evolve_state(self.TIMESTEP + i, ["virion_growth"])
            curr_count = self.sx_container.progeny_state_count[ViralProgeny.NUCLEATE_CYT]
            self.assertTrue(curr_count <= prev_count)

        self.sx_container.NUCLEATE_DISS_RATE = 1
        for i in range(50):
            self.sx_pack_proc.evolve_state(self.TIMESTEP + i, ["virion_growth"])
            curr_count = self.sx_container.progeny_state_count[ViralProgeny.NUCLEATE_MEM]
            self.assertEqual(curr_count, 0)

    # Poisson
    # self.GAG_DIMER_DIFFUSION_FOLD_CHANGE
    def test_Gag_dimers_progeny_binding(self):
        # Gag_dimers bins should remain constant
        # rate_of_Gag_dim_cyt_binding, rate_of_Gag_dim_mem_binding = 0
        self.sx_container.GAG_DIMER_DIFFUSION_FOLD_CHANGE = 0
        prev_count = self.sx_proteins_cyt[Proteins.index['Gag_dimers']] + self.sx_proteins_mem[Proteins.index['Gag_dimers']]
        for i in range(50):
            self.sx_pack_proc.evolve_state(self.TIMESTEP + i, ["virion_growth"])
        curr_count = self.sx_proteins_cyt[Proteins.index['Gag_dimers']] + self.sx_proteins_mem[Proteins.index['Gag_dimers']]
        self.assertEqual(prev_count, curr_count)

    # Poisson
    # self.GAG_LATERAL_DIFFUSION_FOLD_CHANGE
    def test_Gag_progeny_binding(self):
        # No membrane Gag/Gag dimer binding with progeny
        self.sx_container.GAG_LATERAL_DIFFUSION_FOLD_CHANGE = 0

        prev_count = self.sx_proteins_mem[Proteins.index['Gag']] + self.sx_proteins_mem[Proteins.index['Gag_dimers']]
        for i in range(50):
            self.sx_pack_proc.evolve_state(self.TIMESTEP + i, ["virion_growth"])
        curr_count = self.sx_proteins_mem[Proteins.index['Gag']] + self.sx_proteins_mem[Proteins.index['Gag_dimers']]
        self.assertEqual(prev_count, curr_count)

    # Start out with zero of everything  
    def test_zero_abundance(self):
        for i in range(30):
            self.pack_proc.evolve_state(i)
            if (abundance_is_nonzero(self.state)):
                self.fail('Abundances should be zero')
        self.assertEqual(1,1)
    
    # Make sure abundance values are non-negative integers at all times
    def test_non_negative_integer_abundances(self):
        for i in range(30):
            self.sx_pack_proc.evolve_state(i)
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

class TestPackaging1200(TestPackaging):
    TIMESTEP = 1200
    
if __name__ == '__main__':
    unittest.main()

