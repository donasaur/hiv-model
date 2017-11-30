from process.Transcription import Transcription
from mainaux.TestHelpers import *
from state.Proteins import *
import unittest
import numpy as np

class TestTranscription(unittest.TestCase):
    
    # This method is run prior to each test automatically.
    def setUp(self):
        # Default state
        self.state = State()
        self.transcribe_proc = Transcription(self.state)

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
        self.s1_transcribe_proc = Transcription(self.s1_state)

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
        self.s2_transcribe_proc = Transcription(self.s2_state)

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
        self.s3_transcribe_proc = Transcription(self.s3_state)

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
        self.s4_transcribe_proc = Transcription(self.s4_state)

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

    # Change the value of PROMOTER_ON_RATE and PROMOTER_OFF_RATE
    # and see how it affects the value of promoter_activity accordingly
    def test_integration_site_param(self):
        self.s3_transcribe_proc.PROMOTER_ON_RATE = 0
        self.s3_transcribe_proc.PROMOTER_OFF_RATE = 1
        self.s3_transcribe_proc.THRESH_TAT_FEEDBACK = 1000
        self.s3_state.get_state('DNAs').promoter_activity = 0
        reaction_rate_state = self.s3_state.get_state('reaction_rates')
        mRNA_state = self.s3_state.get_state('mRNAs')
        reaction_rate_state.Tat_derived_transcription_rate = 0

        for i in range(30):
            old_full_len_bins = np.copy(mRNA_state.full_len_transcripts_nuc)
            self.s3_transcribe_proc.evolve_state(i)
            self.assertTrue(self.s3_state.get_state('DNAs').promoter_activity == 0)
            self.assertTrue(np.array_equal(old_full_len_bins, mRNA_state.full_len_transcripts_nuc)) # Bins should not change (no creation of transcript)

        self.setUp()

        self.s3_transcribe_proc.PROMOTER_ON_RATE = 1
        self.s3_transcribe_proc.PROMOTER_OFF_RATE = 0
        self.s3_transcribe_proc.THRESH_TAT_FEEDBACK = 1000        
        self.s3_state.get_state('DNAs').promoter_activity = 0
        reaction_rate_state = self.s3_state.get_state('reaction_rates')
        mRNA_state = self.s3_state.get_state('mRNAs')
        reaction_rate_state.Tat_derived_transcription_rate = 0

        for i in range(30):
            old_full_len_bins = np.copy(mRNA_state.full_len_transcripts_nuc)
            self.s3_transcribe_proc.evolve_state(i)
            self.assertTrue(self.s3_state.get_state('DNAs').promoter_activity == 1)
            self.assertTrue(mRNA_state.full_len_transcripts_nuc[0] >= old_full_len_bins[0])

        self.setUp()

        self.s3_transcribe_proc.PROMOTER_OFF_RATE = 0
        self.s3_transcribe_proc.PROMOTER_ON_RATE = 1
        self.s3_transcribe_proc.THRESH_TAT_FEEDBACK = 1000        
        self.s3_state.get_state('DNAs').promoter_activity = 1
        reaction_rate_state = self.s3_state.get_state('reaction_rates')
        mRNA_state = self.s3_state.get_state('mRNAs')
        reaction_rate_state.Tat_derived_transcription_rate = 0

        for i in range(30):
            old_full_len_bins = np.copy(mRNA_state.full_len_transcripts_nuc)
            self.s3_transcribe_proc.evolve_state(i)
            self.assertTrue(self.s3_state.get_state('DNAs').promoter_activity == 1)
            self.assertTrue(mRNA_state.full_len_transcripts_nuc[0] >= old_full_len_bins[0])

        self.setUp()

        self.s3_transcribe_proc.PROMOTER_OFF_RATE = 1
        self.s3_transcribe_proc.PROMOTER_ON_RATE = 0
        self.s3_transcribe_proc.THRESH_TAT_FEEDBACK = 1000   
        self.s3_state.get_state('DNAs').promoter_activity = 1
        reaction_rate_state = self.s3_state.get_state('reaction_rates')
        mRNA_state = self.s3_state.get_state('mRNAs')
        reaction_rate_state.Tat_derived_transcription_rate = 0

        for i in range(30):
            old_full_len_bins = np.copy(mRNA_state.full_len_transcripts_nuc)
            self.s3_transcribe_proc.evolve_state(i)
            self.assertTrue(self.s3_state.get_state('DNAs').promoter_activity == 0)
            self.assertTrue(np.array_equal(old_full_len_bins, mRNA_state.full_len_transcripts_nuc))  # Bins should not change (no creation of transcript)    

    # Makes sure the Rev count does not change during iterations of this process
    def test_Rev_balance(self):
        for i in range(30):
            prev_Rev_count = count_total_Rev(self.s2_state)
            self.s2_transcribe_proc.evolve_state(i)
            curr_Rev_count = count_total_Rev(self.s2_state)

            self.assertEqual(prev_Rev_count, curr_Rev_count)

        for i in range(30):
            prev_Rev_count = count_total_Rev(self.s3_state)
            self.s3_transcribe_proc.evolve_state(i)
            curr_Rev_count = count_total_Rev(self.s3_state)

            self.assertEqual(prev_Rev_count, curr_Rev_count)

        for i in range(30):
            prev_Rev_count = count_total_Rev(self.s4_state)
            self.s4_transcribe_proc.evolve_state(i)
            curr_Rev_count = count_total_Rev(self.s4_state)
        
            self.assertEqual(prev_Rev_count, curr_Rev_count)

    # Make sure that all other bins besides the 0-Rev-bound bin remains constant
    def test_unaffected_transcript_bal(self):
        for i in range(30):
            prev_transcript_count = count_total_mRNA(self.s3_state)
            prev_count_minus_affected_bin = prev_transcript_count - self.s3_full_len_transcripts_nuc[0]
            self.s3_transcribe_proc.evolve_state(i)
            curr_transcript_count = count_total_mRNA(self.s3_state)
            curr_count_minus_affected_bin = curr_transcript_count - self.s3_full_len_transcripts_nuc[0]
            self.assertEqual(prev_count_minus_affected_bin,curr_count_minus_affected_bin)

    # See how the value of the Tat_derived_transcription_rate (magnitude: 52) affects the
    # value of PROMOTER_ALWAYS_ON and the ACTUAL_TRANSCRIPTION_RATE
    def test_Tat_derived_rate(self):
        self.s3_transcribe_proc.BASAL_TRANSCRIPTION_RATE = 60
        self.s3_transcribe_proc.THRESH_TAT_FEEDBACK = 70

        self.s3_transcribe_proc.evolve_state(200)

        self.assertFalse(self.s3_transcribe_proc.PROMOTER_ALWAYS_ON)
        self.assertEqual(self.s3_transcribe_proc.ACTUAL_TRANSCRIPTION_RATE, self.s3_transcribe_proc.BASAL_TRANSCRIPTION_RATE)

        self.setUp()

        self.s3_transcribe_proc.BASAL_TRANSCRIPTION_RATE = 60
        self.s3_transcribe_proc.THRESH_TAT_FEEDBACK = 50

        self.s3_transcribe_proc.evolve_state(200)

        self.assertTrue(self.s3_transcribe_proc.PROMOTER_ALWAYS_ON)
        self.assertEqual(self.s3_transcribe_proc.ACTUAL_TRANSCRIPTION_RATE, self.s3_transcribe_proc.BASAL_TRANSCRIPTION_RATE)

        self.setUp()

        Tat_derived_transcription_rate = self.s3_state.get_state('reaction_rates').Tat_derived_transcription_rate

        self.s3_transcribe_proc.BASAL_TRANSCRIPTION_RATE = 50
        self.s3_transcribe_proc.THRESH_TAT_FEEDBACK = 60

        self.s3_transcribe_proc.evolve_state(200)

        self.assertFalse(self.s3_transcribe_proc.PROMOTER_ALWAYS_ON)
        self.assertEqual(self.s3_transcribe_proc.ACTUAL_TRANSCRIPTION_RATE, Tat_derived_transcription_rate)

        self.setUp()

        Tat_derived_transcription_rate = self.s3_state.get_state('reaction_rates').Tat_derived_transcription_rate

        self.s3_transcribe_proc.BASAL_TRANSCRIPTION_RATE = 40
        self.s3_transcribe_proc.THRESH_TAT_FEEDBACK = 50

        self.s3_transcribe_proc.evolve_state(200)

        self.assertTrue(self.s3_transcribe_proc.PROMOTER_ALWAYS_ON)
        self.assertEqual(self.s3_transcribe_proc.ACTUAL_TRANSCRIPTION_RATE, Tat_derived_transcription_rate)         

    # With state s3, make sure that the promoter is always on
    # The promoter should always be on when the Tat derived transcription rate is greater than the threshold rate
    # TODO: Check to make sure it's a valid assumption to assume that promoter is always ON
    # after derived Rate > threshold rate for the first time (even if Rate derived < threshold rate in later timestep)
    def test_promoter_always_on(self):
        DNA_state = self.s3_state.get_state('DNAs')

        for i in np.arange(200,230):
            self.s3_transcribe_proc.evolve_state(i)
            self.assertTrue(self.s3_transcribe_proc.PROMOTER_ALWAYS_ON)
            self.assertEqual(DNA_state.promoter_activity, 1)

    # Checks that the actual transcription rate is bounded by the MAX_TAT_ENHANCEMENT parameter appropriately
    def test_upper_limit_transcription_rate(self):
        self.setUp()
        self.s2_transcribe_proc.evolve_state(100)
        self.assertNotEquals(self.s2_transcribe_proc.ACTUAL_TRANSCRIPTION_RATE, self.s2_transcribe_proc.THRESH_TAT_FEEDBACK*self.s2_transcribe_proc.MAX_TAT_ENHANCEMENT)

        self.setUp()
        self.s3_transcribe_proc.evolve_state(200)
        self.assertEqual(self.s3_transcribe_proc.ACTUAL_TRANSCRIPTION_RATE, self.s3_transcribe_proc.THRESH_TAT_FEEDBACK*self.s3_transcribe_proc.MAX_TAT_ENHANCEMENT)

        self.setUp()
        self.s4_transcribe_proc.evolve_state(300)
        self.assertEqual(self.s4_transcribe_proc.ACTUAL_TRANSCRIPTION_RATE, self.s4_transcribe_proc.THRESH_TAT_FEEDBACK*self.s4_transcribe_proc.MAX_TAT_ENHANCEMENT)

    # Makes sure the increase in the 0-Rev-bound GagProPol bin corresponds to
    # the increase in the number of transcripts synthesized
    def test_num_transcripts_transcribed(self):
        mRNA_state = self.s3_state.get_state('mRNAs')
        for i in range(200, 230):
            prev_zero_Rev_bound_GagProPol_count = self.s3_full_len_transcripts_nuc[0]
            prev_transcript_synthesized_count = mRNA_state.transcripts_synthesized
            self.s3_transcribe_proc.evolve_state(i)
            curr_zero_Rev_bound_GagProPol_count = self.s3_full_len_transcripts_nuc[0]
            curr_transcript_synthesized_count = mRNA_state.transcripts_synthesized
            self.assertEqual(curr_zero_Rev_bound_GagProPol_count-prev_zero_Rev_bound_GagProPol_count, curr_transcript_synthesized_count-prev_transcript_synthesized_count)

    def test_promoter_on_zero_rate(self):
        reaction_rate_state = self.s3_state.get_state('reaction_rates')
        reaction_rate_state.Tat_derived_transcription_rate = 0.00000001
        self.s3_transcribe_proc.THRESH_TAT_FEEDBACK = 0
        self.s3_transcribe_proc.BASAL_TRANSCRIPTION_RATE = 0
        mRNA_state = self.s3_state.get_state('mRNAs')

        for i in range(30):
            old_full_len_bins = np.copy(mRNA_state.full_len_transcripts_nuc)
            self.s3_transcribe_proc.evolve_state(i)
            self.assertTrue(self.s3_state.get_state('DNAs').promoter_activity == 1)
            self.assertTrue(np.array_equal(old_full_len_bins, mRNA_state.full_len_transcripts_nuc))  # Bins should not change (no creation of transcript)    
        
    # Start out with zero of everything  
    def test_zero_abundance(self):
        self.transcribe_proc.BASAL_TRANSCRIPTION_RATE = 0

        for i in range(30):
            self.transcribe_proc.evolve_state(i)
            if (abundance_is_nonzero(self.state)):
                self.fail('Abundances should be zero')
        self.assertEqual(1,1)
    
    # Make sure abundance values are non-negative integers at all times
    def test_non_negative_integer_abundances(self):
        for i in range(30):
            self.s3_transcribe_proc.evolve_state(i)
            if abundance_is_negative(self.s3_state):
                self.fail('One or more abundances is a negative value.')
            elif abundance_is_not_integer(self.s3_state):
                self.fail('One or more abundances is not an integer.')
            else:
                self.assertEquals(1,1)
        
    # Use this method if I want to change any conditions
    # immediately after running a test.
    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()

