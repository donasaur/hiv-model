from process.TatFeedback import *
from mainaux.TestHelpers import *
from state.Proteins import *
import unittest
import numpy as np
import math

class TestTatFeedback(unittest.TestCase):
    
    # This method is run prior to each test automatically.
    def setUp(self):
        # Default state
        self.state = State()
        self.tfeedback_proc = TatFeedback(self.state)

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
        self.s1_tfeedback_proc = TatFeedback(self.s1_state)

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
        self.s2_tfeedback_proc = TatFeedback(self.s2_state)

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
        self.s3_tfeedback_proc = TatFeedback(self.s3_state)

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
        self.s4_tfeedback_proc = TatFeedback(self.s4_state)

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

    # Make sure that the number of Tat remains the same in the
    # system after evolve_state is called at each timestep
    # Make sure to consider pTEFb-Tat complex (deacetyl + acetyl)
    def test_Tat_balance(self):
        for i in range(30):
            prev_Tat_count = count_total_Tat(self.s1_state)
            self.s1_tfeedback_proc.evolve_state(i)
            curr_Tat_count = count_total_Tat(self.s1_state)

            self.assertEqual(prev_Tat_count, curr_Tat_count)

        for i in range(30):
            prev_Tat_count = count_total_Tat(self.s2_state)
            self.s2_tfeedback_proc.evolve_state(i)
            curr_Tat_count = count_total_Tat(self.s2_state)

            self.assertEqual(prev_Tat_count, curr_Tat_count)

        for i in range(30):
            prev_Tat_count = count_total_Tat(self.s3_state)
            self.s3_tfeedback_proc.evolve_state(i)
            curr_Tat_count = count_total_Tat(self.s3_state)

            self.assertEqual(prev_Tat_count, curr_Tat_count)

        for i in range(30):
            prev_Tat_count = count_total_Tat(self.s4_state)
            self.s4_tfeedback_proc.evolve_state(i)
            curr_Tat_count = count_total_Tat(self.s4_state)

            self.assertEqual(prev_Tat_count, curr_Tat_count)

    # Change pTEFb creation rate to zero
    def test_pTEFb_balance(self):
        self.s3_tfeedback_proc.pTEFb_DOUBLING_RATE = 0
        for i in range(30):
            prev_pTEFb_bal = count_total_pTEFb(self.s3_state)
            self.s3_tfeedback_proc.evolve_state(i)
            curr_pTEFb_bal = count_total_pTEFb(self.s3_state)

            self.assertEqual(prev_pTEFb_bal, curr_pTEFb_bal)

    # The following tests change the rate constants that control
    # the abundance of Tat in a variety of forms

    def test_free_Tat_rate(self):
        self.s1_tfeedback_proc.RATE_TAT_ACT_TRANSCRIPTION = 0
        self.s1_tfeedback_proc.RATE_TAT_pTEFb_BIND = 0

        for i in range(30):
            prev_free_Tat_in_nuc = self.s1_proteins_nuc[self.index['Tat']]
            self.s1_tfeedback_proc.evolve_state(i)
            curr_free_Tat_in_nuc = self.s1_proteins_nuc[self.index['Tat']]

            self.assertEqual(prev_free_Tat_in_nuc, curr_free_Tat_in_nuc)

    def test_free_pTEFb_rate(self):
        self.s3_tfeedback_proc.RATE_TAT_ACT_TRANSCRIPTION = 0
        self.s3_tfeedback_proc.RATE_TAT_pTEFb_BIND = 0

        for i in np.arange(201,231):
            prev_free_pTEFb = get_factors(self.s3_state)[1]
            self.s3_tfeedback_proc.evolve_state(i)
            curr_free_pTEFb = get_factors(self.s3_state)[1]

            self.assertEqual(math.floor(prev_free_pTEFb + np.around(get_factors(self.s3_state)[0]*np.exp(self.s3_tfeedback_proc.pTEFb_DOUBLING_RATE*(i+1))) - np.around(get_factors(self.s3_state)[0]*np.exp(self.s3_tfeedback_proc.pTEFb_DOUBLING_RATE*(i)))), math.floor(curr_free_pTEFb))

    def test_Tat_pTEFb_deacetyl_rate(self):
        self.s1_tfeedback_proc.RATE_TAT_pTEFb_BIND = 0
        self.s1_tfeedback_proc.RATE_TAT_pTEFb_ACETYL = 0

        for i in range(30):
            prev_free_ptd = get_factors(self.s1_state)[2]
            self.s1_tfeedback_proc.evolve_state(i)
            curr_free_ptd = get_factors(self.s1_state)[2]

            self.assertEqual(prev_free_ptd, curr_free_ptd)

    def test_Tat_pTEFb_acetyl_rate(self):
        self.s1_tfeedback_proc.RATE_TAT_pTEFb_ACETYL = 0
        self.s1_tfeedback_proc.RATE_TAT_ACT_TRANSCRIPTION = 0

        for i in range(30):
            prev_free_pta = get_factors(self.s1_state)[3]
            self.s1_tfeedback_proc.evolve_state(i)
            curr_free_pta = get_factors(self.s1_state)[3]

            self.assertEqual(prev_free_pta, curr_free_pta)

    def test_ODE_discretizer(self):
        soln = np.zeros((60, 5))
        soln_last_row = soln[-1,:]

        # Test 1: Check whole number solution
        # The output from the ODE discretizer should return results that are the same as the first 4 bins of the last row of soln
        # Also make sure that the # of Tat + # of pTEFb is preserved across the function call
        soln[-1,:] = [30, 160, 40, 50, 102]
        prev_abundances = np.array([70, 200, 20, 30])

        free_Tat, pTEFb_nuc, Tat_pTEFb_deacetyl, Tat_pTEFb_acetyl = self.tfeedback_proc.ODE_discretizer(soln, *prev_abundances)

        self.assertTrue(np.array_equal(np.array([free_Tat, pTEFb_nuc, Tat_pTEFb_deacetyl, Tat_pTEFb_acetyl]), np.array([30, 160, 40, 50])))
        self.assertEqual(prev_abundances[0] + prev_abundances[2] + prev_abundances[3], free_Tat + Tat_pTEFb_deacetyl + Tat_pTEFb_acetyl)
        self.assertEqual(prev_abundances[1] + prev_abundances[2] + prev_abundances[3], pTEFb_nuc + Tat_pTEFb_deacetyl + Tat_pTEFb_acetyl)

        # Test 2: Check a soln where the bins are the same in value as the prev abundance bins
        soln[-1,:] = [70, 200, 20, 30, 202]
        prev_abundances = np.array([70, 200, 20, 30])

        free_Tat, pTEFb_nuc, Tat_pTEFb_deacetyl, Tat_pTEFb_acetyl = self.tfeedback_proc.ODE_discretizer(soln, *prev_abundances)

        self.assertTrue(np.array_equal(np.array([free_Tat, pTEFb_nuc, Tat_pTEFb_deacetyl, Tat_pTEFb_acetyl]), prev_abundances))
        self.assertEqual(prev_abundances[0] + prev_abundances[2] + prev_abundances[3], free_Tat + Tat_pTEFb_deacetyl + Tat_pTEFb_acetyl)
        self.assertEqual(prev_abundances[1] + prev_abundances[2] + prev_abundances[3], pTEFb_nuc + Tat_pTEFb_deacetyl + Tat_pTEFb_acetyl)

        # Test 3: Let two of the three Tat bins be negative and estimate in advance how the ODE_discretizer function behaves in response
        soln[-1,:] = [-5, 200, -10, 135, 10]
        prev_abundances = np.array([70, 275, 20, 30])

        free_Tat, pTEFb_nuc, Tat_pTEFb_deacetyl, Tat_pTEFb_acetyl = self.tfeedback_proc.ODE_discretizer(soln, *prev_abundances)

        self.assertTrue(np.array_equal(np.array([free_Tat, pTEFb_nuc, Tat_pTEFb_deacetyl, Tat_pTEFb_acetyl]), np.array([0,205,0,120])))
        self.assertEqual(prev_abundances[0] + prev_abundances[2] + prev_abundances[3], free_Tat + Tat_pTEFb_deacetyl + Tat_pTEFb_acetyl)
        self.assertEqual(prev_abundances[1] + prev_abundances[2] + prev_abundances[3], pTEFb_nuc + Tat_pTEFb_deacetyl + Tat_pTEFb_acetyl)

        # Test 4: Rounding case
        soln[-1,:] = [5.5, 150.5, 10.5, 100, 30]
        prev_abundances = [5, 150, 10, 101]

        free_Tat, pTEFb_nuc, Tat_pTEFb_deacetyl, Tat_pTEFb_acetyl = self.tfeedback_proc.ODE_discretizer(soln, *prev_abundances)

        self.assertEqual(prev_abundances[0] + prev_abundances[2] + prev_abundances[3], free_Tat + Tat_pTEFb_deacetyl + Tat_pTEFb_acetyl)

        # Expect the output
        # (5, 150, 11, 100) or (6, 151, 10, 100) or (6, 151, 11, 99) from ODE discretizer
        self.assertTrue(free_Tat == 5 or free_Tat == 6)
        self.assertTrue(pTEFb_nuc == 150 or pTEFb_nuc == 151)
        self.assertTrue(Tat_pTEFb_deacetyl == 10 or Tat_pTEFb_deacetyl == 11)
        self.assertTrue(Tat_pTEFb_acetyl == 100 or Tat_pTEFb_deacetyl == 99)

        # Test 5: pTEFb testing

    # Start out with zero of everything  
    def test_zero_abundance(self):
        for i in range(30):
            self.tfeedback_proc.evolve_state(i)
            if (abundance_is_nonzero(self.state)):
                self.fail('Abundances should be zero')
        self.assertEqual(1,1)
    
    # Make sure abundance values are non-negative integers at all times
    def test_non_negative_integer_abundances(self):
        for i in range(30):
            self.s1_tfeedback_proc.evolve_state(i)
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

