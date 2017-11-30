# -*- coding: utf-8 -*-
from mainaux.State import State
from state.Proteins import Proteins
import numpy as np
import cPickle
import os

def load_state(batch_num, sim_num, timestep):
    path = os.getcwd() + "/samplestates/Batch" + str(batch_num) + "/Sim" + str(sim_num) + "/"
    new_filename = "sample" + str(timestep)
    fullpath = os.path.join(path, new_filename)
    fileHandler = open(fullpath, 'rb')
    rtn_state = cPickle.load(fileHandler)
    fileHandler.close()
    return rtn_state

def count_total_Rev(state):
    # Setup stuff
    full_splice_matrix = np.array([0,1,2,3,4,5,6,7,8])        
        
    subset1 = np.concatenate((0*np.ones(7), 1*np.ones(7), 2*np.ones(7)))
    subset2 = np.concatenate((3*np.ones(7), 4*np.ones(7), 5*np.ones(7)))
    subset3 = np.concatenate((6*np.ones(7), 7*np.ones(7), 8*np.ones(7)))
    single_splice_matrix = np.concatenate((subset1, subset2, subset3))    
    
    mRNA_state = state.get_state('mRNAs')
    protein_state = state.get_state('proteins')

    proteins_nuc = protein_state.proteins_nuc # proteins_nuc[Proteins.Proteins.index['Rev']] is Rev count in nucleus
    proteins_cyt = protein_state.proteins_cyt # proteins_cyt[Proteins.index['Rev']] is Rev count in cytoplasm

    full_len_transcripts_nuc = mRNA_state.full_len_transcripts_nuc
    full_len_transcripts_cyt = mRNA_state.full_len_transcripts_cyt
    single_splice_transcript_nuc = mRNA_state.single_splice_transcript_nuc
    single_splice_transcript_cyt = mRNA_state.single_splice_transcript_cyt
    multi_splice_transcript_nuc = mRNA_state.multi_splice_transcript_nuc
    multi_splice_transcript_cyt = mRNA_state.multi_splice_transcript_cyt

    Rev_count_in_nuc = proteins_nuc[Proteins.index['Rev']] + np.sum(np.dot(full_splice_matrix,full_len_transcripts_nuc)) + np.sum(np.dot(single_splice_matrix,single_splice_transcript_nuc))
    Rev_count_in_cyt = proteins_cyt[Proteins.index['Rev']] + np.sum(np.dot(full_splice_matrix,full_len_transcripts_cyt)) + np.sum(np.dot(single_splice_matrix,single_splice_transcript_cyt))
    Rev_count = Rev_count_in_nuc + Rev_count_in_cyt

    return Rev_count

def count_bound_Rev(mRNA_array):
    full_splice_matrix = np.array([0,1,2,3,4,5,6,7,8])
    return np.dot(full_splice_matrix,mRNA_array)

def get_factors(state):
    host_factor_state = state.get_state('host_factors')
    return [host_factor_state.pTEFb_nuc_init, host_factor_state.pTEFb_nuc, host_factor_state.Tat_pTEFb_deacetyl, host_factor_state.Tat_pTEFb_acetyl]

def count_total_pTEFb(state):
    return get_factors(state)[1] + get_factors(state)[2] + get_factors(state)[3]

def count_total_Tat(state):
    mRNA_state = state.get_state('mRNAs')
    protein_state = state.get_state('proteins')
    host_factor_state = state.get_state('host_factors')

    proteins_nuc = protein_state.proteins_nuc # proteins_nuc[Proteins.index['Rev']] is Rev count in nucleus
    proteins_cyt = protein_state.proteins_cyt # proteins_cyt[Proteins.index['Rev']] is Rev count in cytoplasm

    Tat_count_in_nuc = proteins_nuc[Proteins.index['Tat']]
    Tat_count_in_cyt = proteins_cyt[Proteins.index['Tat']]
    Tat_count = Tat_count_in_nuc + Tat_count_in_cyt + host_factor_state.Tat_pTEFb_deacetyl + host_factor_state.Tat_pTEFb_acetyl

    return Tat_count

def count_total_Gag(state):
    protein_state = state.get_state('proteins')
    mRNA_state = state.get_state('mRNAs')

    multiplicity_array = [0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4]
    total_Gag_from_mRNAs = np.dot(mRNA_state.full_len_transcripts_Gag_bound, multiplicity_array)

    proteins_nuc = protein_state.proteins_nuc # proteins_nuc[Proteins.index['Rev']] is Rev count in nucleus
    proteins_cyt = protein_state.proteins_cyt # proteins_cyt[Proteins.index['Rev']] is Rev count in cytoplasm
    proteins_mem = protein_state.proteins_mem
    proteins_virion = protein_state.proteins_virion

    sum_of_Gag = proteins_nuc[Proteins.index['Gag']] + proteins_nuc[Proteins.index['GagProPol']] + 2 * proteins_nuc[Proteins.index['Gag_dimers']]
    sum_of_Gag += proteins_cyt[Proteins.index['Gag']] + proteins_cyt[Proteins.index['GagProPol']] + 2 * proteins_cyt[Proteins.index['Gag_dimers']]
    sum_of_Gag += proteins_mem[Proteins.index['Gag']] + proteins_mem[Proteins.index['GagProPol']] + 2 * proteins_mem[Proteins.index['Gag_dimers']]
    sum_of_Gag += proteins_virion[Proteins.index['Gag']] + proteins_virion[Proteins.index['GagProPol']] + 2 * proteins_virion[Proteins.index['Gag_dimers']]
    sum_of_Gag += total_Gag_from_mRNAs
    return sum_of_Gag

def count_total_Vif(state):
    protein_state = state.get_state('proteins')
    
    proteins_nuc = protein_state.proteins_nuc # proteins_nuc[Proteins.index['Rev']] is Rev count in nucleus
    proteins_cyt = protein_state.proteins_cyt # proteins_cyt[Proteins.index['Rev']] is Rev count in cytoplasm
    proteins_mem = protein_state.proteins_mem
    proteins_virion = protein_state.proteins_virion

    sum_of_Vif = proteins_nuc[Proteins.index['Vif']] + proteins_cyt[Proteins.index['Vif']] + proteins_mem[Proteins.index['Vif']] + proteins_virion[Proteins.index['Vif']]

    return sum_of_Vif

def count_active_gag_mRNA(state):
    mRNA_state = state.get_state('mRNAs')
    container = state.get_state('viral_progeny_container')
    gag_mRNA_count = 0
    gag_mRNA_count += np.sum(mRNA_state.full_len_transcripts_cyt)
    gag_mRNA_count += np.sum(mRNA_state.full_len_transcripts_Gag_bound)
    gag_mRNA_count += 2*container.count_progeny()
    return gag_mRNA_count

def count_total_gag_mRNA(state):
    mRNA_state = state.get_state('mRNAs')
    container = state.get_state('viral_progeny_container')

    gag_mRNA_count = 0
    gag_mRNA_count += np.sum(mRNA_state.full_len_transcripts_cyt)
    gag_mRNA_count += np.sum(mRNA_state.full_len_transcripts_Gag_bound)
    gag_mRNA_count += np.sum(mRNA_state.full_len_transcripts_nuc)
    gag_mRNA_count += 2*container.count_progeny()
    return gag_mRNA_count

def count_total_mRNA(state, split=False):
    mRNA_state = state.get_state('mRNAs')

    full_len_transcripts_nuc = mRNA_state.full_len_transcripts_nuc
    full_len_transcripts_cyt = mRNA_state.full_len_transcripts_cyt
    single_splice_transcript_nuc = mRNA_state.single_splice_transcript_nuc
    single_splice_transcript_cyt = mRNA_state.single_splice_transcript_cyt
    multi_splice_transcript_nuc = mRNA_state.multi_splice_transcript_nuc
    multi_splice_transcript_cyt = mRNA_state.multi_splice_transcript_cyt

    mRNA_count_in_nucleus = full_len_transcripts_nuc.sum() + single_splice_transcript_nuc.sum() + multi_splice_transcript_nuc.sum()
    mRNA_count_in_cytoplasm = full_len_transcripts_cyt.sum() + single_splice_transcript_cyt.sum() + multi_splice_transcript_cyt.sum()
    mRNA_count = mRNA_count_in_nucleus + mRNA_count_in_cytoplasm

    if split == True:
        return mRNA_count_in_nucleus, mRNA_count_in_cytoplasm
    return mRNA_count

def sum_ss_type(state, i, split=False):
    mRNA_state = state.get_state('mRNAs')

    single_splice_transcript_nuc = mRNA_state.single_splice_transcript_nuc
    single_splice_transcript_cyt = mRNA_state.single_splice_transcript_cyt

    nuc_sum = np.sum(single_splice_transcript_nuc[np.arange(i, 63, 7)])
    cyt_sum = np.sum(single_splice_transcript_cyt[np.arange(i, 63, 7)])

    if split == True:
        return nuc_sum, cyt_sum
    return nuc_sum + cyt_sum

def count_potential_mRNA(state, min_Rev_req, location="nuc"):
    mRNA_state = state.get_state('mRNAs')

    full_len_transcripts_nuc = mRNA_state.full_len_transcripts_nuc
    full_len_transcripts_cyt = mRNA_state.full_len_transcripts_cyt
    single_splice_transcript_nuc = mRNA_state.single_splice_transcript_nuc
    single_splice_transcript_cyt = mRNA_state.single_splice_transcript_cyt
    multi_splice_transcript_nuc = mRNA_state.multi_splice_transcript_nuc
    multi_splice_transcript_cyt = mRNA_state.multi_splice_transcript_cyt

    if location == "nuc":
        potential_mRNA_count = full_len_transcripts_nuc[min_Rev_req:].sum() + single_splice_transcript_nuc[7*min_Rev_req:].sum() + multi_splice_transcript_nuc.sum()
    else:
        potential_mRNA_count = full_len_transcripts_cyt[min_Rev_req:].sum() + single_splice_transcript_cyt[7*min_Rev_req:].sum() + multi_splice_transcript_nuc.sum()

    return potential_mRNA_count

# Set the additional argument option to 1 if want the function to return a tuple of the form [protein_count_in_nucleus, protein_count_in_cytoplasm] instead
def count_total_protein(state, option=0):
    protein_state = state.get_state('proteins')

    proteins_nuc = protein_state.proteins_nuc # proteins_nuc[Proteins.index['Rev']] is Rev count in nucleus
    proteins_cyt = protein_state.proteins_cyt # proteins_cyt[Proteins.index['Rev']] is Rev count in cytoplasm

    protein_count_in_nucleus = proteins_nuc.sum()
    protein_count_in_cytoplasm = proteins_cyt.sum()
    protein_count = protein_count_in_nucleus + protein_count_in_cytoplasm

    if (option == 0):
        return protein_count
    else:
        return protein_count_in_nucleus, protein_count_in_cytoplasm

def count_active_Env(state):
    protein_state = state.get_state('proteins') 

    env_misc = protein_state.env_misc
    proteins_cyt = protein_state.proteins_cyt
    proteins_virion = protein_state.proteins_virion
    active_Env_count = 0
    for env_state, qty in protein_state.env_misc.iteritems():
        if 'trimers' in env_state:
            active_Env_count += 3 * sum(qty)
        else:
            active_Env_count += qty
    active_Env_count += proteins_cyt[Proteins.index['Env']]
    active_Env_count += proteins_virion[Proteins.index['Env']]
    return active_Env_count

def count_inactive_Env(state):
    protein_state = state.get_state('proteins') 

    proteins_nuc = protein_state.proteins_nuc
    proteins_mem = protein_state.proteins_mem

    inactive_Env_count = proteins_nuc[Proteins.index['Env']] + proteins_mem[Proteins.index['Env']]

    return inactive_Env_count

def count_total_Env(state):
    return count_active_Env(state) + count_inactive_Env(state)

def abundance_is_nonzero(state):
    mRNA_state = state.get_state('mRNAs')
    protein_state = state.get_state('proteins')

    proteins_nuc = protein_state.proteins_nuc # proteins_nuc[Proteins.index['Rev']] is Rev count in nucleus
    proteins_cyt = protein_state.proteins_cyt # proteins_cyt[Proteins.index['Rev']] is Rev count in cytoplasm

    full_len_transcripts_nuc = mRNA_state.full_len_transcripts_nuc
    full_len_transcripts_cyt = mRNA_state.full_len_transcripts_cyt
    single_splice_transcript_nuc = mRNA_state.single_splice_transcript_nuc
    single_splice_transcript_cyt = mRNA_state.single_splice_transcript_cyt
    multi_splice_transcript_nuc = mRNA_state.multi_splice_transcript_nuc
    multi_splice_transcript_cyt = mRNA_state.multi_splice_transcript_cyt

    nuc_abundance_is_nonzero = proteins_nuc.any() != 0 or full_len_transcripts_nuc.any() != 0 or single_splice_transcript_nuc.any() != 0 or multi_splice_transcript_nuc.any() != 0
    cyt_abundance_is_nonzero = proteins_cyt.any() != 0 or full_len_transcripts_cyt.any() != 0 or single_splice_transcript_cyt.any() != 0 or multi_splice_transcript_cyt.any() != 0
    abundance_is_nonzero_bool = nuc_abundance_is_nonzero or cyt_abundance_is_nonzero

    return abundance_is_nonzero_bool

def abundance_is_negative(state):
    mRNA_state = state.get_state('mRNAs')
    protein_state = state.get_state('proteins')

    proteins_nuc = protein_state.proteins_nuc # proteins_nuc[Proteins.index['Rev']] is Rev count in nucleus
    proteins_cyt = protein_state.proteins_cyt # proteins_cyt[Proteins.index['Rev']] is Rev count in cytoplasm

    full_len_transcripts_nuc = mRNA_state.full_len_transcripts_nuc
    full_len_transcripts_cyt = mRNA_state.full_len_transcripts_cyt
    single_splice_transcript_nuc = mRNA_state.single_splice_transcript_nuc
    single_splice_transcript_cyt = mRNA_state.single_splice_transcript_cyt
    multi_splice_transcript_nuc = mRNA_state.multi_splice_transcript_nuc
    multi_splice_transcript_cyt = mRNA_state.multi_splice_transcript_cyt

    abundance_is_negative_bool = (proteins_nuc.any() < 0 or full_len_transcripts_nuc.any() < 0 or single_splice_transcript_nuc.any() < 0 or multi_splice_transcript_nuc.any() < 0) or (proteins_cyt.any() < 0 or full_len_transcripts_cyt.any() < 0 or single_splice_transcript_cyt.any() < 0 or multi_splice_transcript_cyt.any() < 0)

    return abundance_is_negative_bool

def abundance_is_not_integer(state):
    mRNA_state = state.get_state('mRNAs')
    protein_state = state.get_state('proteins')

    proteins_nuc = protein_state.proteins_nuc # proteins_nuc[Proteins.index['Rev']] is Rev count in nucleus
    proteins_cyt = protein_state.proteins_cyt # proteins_cyt[Proteins.index['Rev']] is Rev count in cytoplasm

    full_len_transcripts_nuc = mRNA_state.full_len_transcripts_nuc
    full_len_transcripts_cyt = mRNA_state.full_len_transcripts_cyt
    single_splice_transcript_nuc = mRNA_state.single_splice_transcript_nuc
    single_splice_transcript_cyt = mRNA_state.single_splice_transcript_cyt
    multi_splice_transcript_nuc = mRNA_state.multi_splice_transcript_nuc
    multi_splice_transcript_cyt = mRNA_state.multi_splice_transcript_cyt

    abundance_is_not_integer_bool = (proteins_nuc.any() % 1 != 0 or full_len_transcripts_nuc.any() % 1 != 0 or single_splice_transcript_nuc.any() % 1 != 0 or multi_splice_transcript_nuc.any() % 1 != 0) or (proteins_cyt.any() % 1 != 0 or full_len_transcripts_cyt.any() % 1 != 0 or single_splice_transcript_cyt.any() % 1 != 0 or multi_splice_transcript_cyt.any() % 1 != 0)

    return abundance_is_not_integer_bool

# Returns an instance of a S1 state
def s1_state():
    state = State()

    mRNA_state = state.get_state('mRNAs')
    protein_state = state.get_state('proteins')
    host_factor_state = state.get_state('host_factors')

    proteins_nuc = protein_state.proteins_nuc # proteins_nuc[Proteins.index['Rev']] is Rev count in nucleus
    proteins_cyt = protein_state.proteins_cyt # proteins_cyt[Proteins.index['Rev']] is Rev count in cytoplasm

    full_len_transcripts_nuc = mRNA_state.full_len_transcripts_nuc
    full_len_transcripts_cyt = mRNA_state.full_len_transcripts_cyt
    single_splice_transcript_nuc = mRNA_state.single_splice_transcript_nuc
    single_splice_transcript_cyt = mRNA_state.single_splice_transcript_cyt
    multi_splice_transcript_nuc = mRNA_state.multi_splice_transcript_nuc
    multi_splice_transcript_cyt = mRNA_state.multi_splice_transcript_cyt
    
    # Initialize the amount of full-length mRNA with different amounts of Rev attached
    for i in range(9):
        full_len_transcripts_nuc[i] = 200*i
        full_len_transcripts_cyt[i] = 100*i

    # Initialize the amount of each single-spliced mRNA type
    for i in range(63):
        single_splice_transcript_nuc[i] = 1000*i
        single_splice_transcript_cyt[i] = 100*i

    # Initialize the amount of each multi-spliced mRNA type
    for i in range(17):
        multi_splice_transcript_nuc[i] = 300*i
        multi_splice_transcript_cyt[i] = 200*i

    # Protein count in the nucleus
    for i in range(7):
        proteins_nuc[i] = 100*i
        proteins_cyt[i] = 200*i

    host_factor_state.Tat_pTEFb_deacetyl = 120
    host_factor_state.Tat_pTEFb_acetyl = 130

    return state

# Snapshot of a state at timestep = 100
#TODO: Rewrite S2, S3, S4 values
def s2_state():
    state = State()

    mRNA_state = state.get_state('mRNAs')
    protein_state = state.get_state('proteins')
    host_factor_state = state.get_state('host_factors')
    DNA_state = state.get_state('DNAs')
    reaction_rate_state = state.get_state('reaction_rates')


    proteins_nuc = protein_state.proteins_nuc # proteins_nuc[Proteins.index['Rev']] is Rev count in nucleus
    proteins_cyt = protein_state.proteins_cyt # proteins_cyt[Proteins.index['Rev']] is Rev count in cytoplasm

    protein_state.proteins_nuc = np.array([0,0,18,0,23,1581,0]);
    protein_state.proteins_cyt = np.array([0,0,8,0,109,0,0]);
    mRNA_state.full_len_transcripts_nuc = np.array([137,60,9,1,1,0,0,0,0]);
    mRNA_state.full_len_transcripts_cyt = np.array([0,0,0,0,0,0,0,0,0]);
    mRNA_state.single_splice_transcript_nuc = np.array([0,1,11,7,0,0,65,0,0,0,0,0,0,13,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]);
    mRNA_state.single_splice_transcript_cyt = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]);
    mRNA_state.multi_splice_transcript_nuc = np.array([0,0,0,0,0,0,0,0,0,0,0,2,5,5,0,0,28]);
    mRNA_state.multi_splice_transcript_cyt = np.array([0,0,0,0,0,0,0,0,0,0,0,0,2,3,0,0,15]);

    host_factor_state.pTEFb_nuc = 417.0
    host_factor_state.Tat_pTEFb_deacetyl = 249.0
    host_factor_state.Tat_pTEFb_acetyl = 2.0

    DNA_state.promoter_activity = 1
    reaction_rate_state.Tat_derived_transcription_rate = 15.0
    reaction_rate_state.translation_suppressed = 0



    return state

# Snapshot of a state at timestep = 200
def s3_state():
    state = State()

    mRNA_state = state.get_state('mRNAs')
    protein_state = state.get_state('proteins')
    host_factor_state = state.get_state('host_factors')
    DNA_state = state.get_state('DNAs')
    reaction_rate_state = state.get_state('reaction_rates')    

    proteins_nuc = protein_state.proteins_nuc # proteins_nuc[Proteins.index['Rev']] is Rev count in nucleus
    proteins_cyt = protein_state.proteins_cyt # proteins_cyt[Proteins.index['Rev']] is Rev count in cytoplasm

    protein_state.proteins_nuc = np.array([968,1674,757,16101,1548,98665,35144]);
    protein_state.proteins_cyt = np.array([0,0,13256,0,7456,0,0]);
    mRNA_state.full_len_transcripts_nuc = np.array([0,0,0,0,0,0,2,27,347]);
    mRNA_state.full_len_transcripts_cyt = np.array([0,0,0,0,0,0,0,0,388]);
    mRNA_state.single_splice_transcript_nuc = np.array([0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,1,2,3,0,0,3,0,1,6,8,0,0,21,3,1,21,23,0,0,148]);
    mRNA_state.single_splice_transcript_cyt = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,13,0,0,194]);
    mRNA_state.multi_splice_transcript_nuc = np.array([3,0,0,0,0,3,4,1,3,0,0,6,35,51,0,0,245]);
    mRNA_state.multi_splice_transcript_cyt = np.array([5,0,0,0,0,5,10,1,1,0,0,4,79,99,0,0,589]);

    host_factor_state.pTEFb_nuc = 0.0
    host_factor_state.Tat_pTEFb_deacetyl = 882.0
    host_factor_state.Tat_pTEFb_acetyl = 9.0

    DNA_state.promoter_activity = 1
    reaction_rate_state.Tat_derived_transcription_rate = 52.0
    reaction_rate_state.translation_suppressed = 0    

    return state

# Snapshot of a state at timestep = 300
def s4_state():
    state = State()

    mRNA_state = state.get_state('mRNAs')
    protein_state = state.get_state('proteins')
    host_factor_state = state.get_state('host_factors')
    DNA_state = state.get_state('DNAs')
    reaction_rate_state = state.get_state('reaction_rates')    

    proteins_nuc = protein_state.proteins_nuc # proteins_nuc[Proteins.index['Rev']] is Rev count in nucleus
    proteins_cyt = protein_state.proteins_cyt # proteins_cyt[Proteins.index['Rev']] is Rev count in cytoplasm

    protein_state.proteins_nuc = np.array([4153,8159,1515,179134,7557,408734,368315]);
    protein_state.proteins_cyt = np.array([0,0,76355,0,58447,0,0]);
    mRNA_state.full_len_transcripts_nuc = np.array([0,0,0,0,0,0,0,3,323]);
    mRNA_state.full_len_transcripts_cyt = np.array([0,0,0,0,0,0,0,0,1234]);
    mRNA_state.single_splice_transcript_nuc = np.array([0,0,8,0,0,0,0,0,1,7,0,0,0,0,0,0,3,0,0,0,0,0,0,2,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,1,1,0,0,0,8,0,4,0,32,0,0,141]);
    mRNA_state.single_splice_transcript_cyt = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,9,55,88,0,0,493]);
    mRNA_state.multi_splice_transcript_nuc = np.array([2,0,0,0,0,1,7,1,0,0,0,0,17,21,0,0,154]);
    mRNA_state.multi_splice_transcript_cyt = np.array([7,0,2,0,0,7,17,2,4,0,0,12,132,188,0,0,977]);

    host_factor_state.pTEFb_nuc = 0.0
    host_factor_state.Tat_pTEFb_deacetyl = 1178.0
    host_factor_state.Tat_pTEFb_acetyl = 12.0

    DNA_state.promoter_activity = 1
    reaction_rate_state.Tat_derived_transcription_rate = 70.0
    reaction_rate_state.translation_suppressed = 0  

    return state

def s5_state():
    state = load_state("T1", 0, 340)
    return state

def s6_state():
    state = load_state("T1", 0, 240)
    return state

def s7_state():
    state = load_state("T2", 0, 2300)
    return state

def s8_state():
    state = load_state("T2", 0, 1500)
    return state

def s9_state():
    state = load_state("T3", 0, 1500)
    return state

def s10_state():
    state = load_state("T3", 0, 1200)
    return state

def site_locator():
    one_index_dict = {}
    zero_index_dict = {}
    zero_index_dict['SL1_site'] = []
    zero_index_dict['SL2_site'] = []
    zero_index_dict['SL3_site'] = []
    zero_index_dict['SL4_site'] = []
    one_index_dict['SL1_site'] = []
    one_index_dict['SL2_site'] = []
    one_index_dict['SL3_site'] = []
    one_index_dict['SL4_site'] = []            
    for i in range(2):
        for j in range(2):
            for k in range(2):
                for l in range(2):
                    if i == 1:
                        one_index_dict['SL1_site'] += [(i,j,k,l)] # Stores all tuples where i == 1 (should be 8)
                        zero_index_dict['SL1_site'] += [(0,j,k,l)] # Stores all tuples where i == 0 (should be 8)
                    if j == 1:
                        one_index_dict['SL2_site'] += [(i,j,k,l)]
                        zero_index_dict['SL2_site'] += [(i,0,k,l)]
                    if k == 1:
                        one_index_dict['SL3_site'] += [(i,j,k,l)]
                        zero_index_dict['SL3_site'] += [(i,j,0,l)]
                    if l == 1:
                        one_index_dict['SL4_site'] += [(i,j,k,l)]
                        zero_index_dict['SL4_site'] += [(i,j,k,0)]
    return one_index_dict, zero_index_dict
one_index_dict, zero_index_dict = site_locator()
