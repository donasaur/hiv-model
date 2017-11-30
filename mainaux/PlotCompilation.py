import numpy as np
import mainaux.PlotHelpers as PlotHelpers

class PlottingKey(object):
    def __init__(self, plotting_key_name):
        self.name = plotting_key_name
        self.dependencies = []
        self.generate_key_data = None # place one-arg function here
        self.plotting_options = []

    def __eq__(self, other):
        return self.key == other.key

def generate_plotting_keys(list_of_key_names):
    list_of_plotting_keys = []
    list_of_dependent_keys = []
    set_of_req_record_key_names = set()

    for plotting_key_name in list_of_key_names:
        plotting_key = PlottingKey(plotting_key_name)

        if plotting_key_name == 'promoter_activity':
            plotting_key.plotting_options = [['Blue'], ('promoter_activity',), 
                                    'Promoter activity', 'Promoter activity', 'off']    
                                   
        elif plotting_key_name == 'full_len_transcripts_nuc':
            plotting_key.plotting_options = [['Pink', 'Red', 'Orange', 'Yellow', 'LawnGreen', 'Green', 'Cyan', 'Blue', 'Purple'],
                                    ('X=0', '1', '2', '3', '4', '5', '6', '7', '8'), 
                                    'Full length transcript abundance in nucleus', 'Transcript abundance', 'on']
        elif plotting_key_name == 'multi_splice_transcript_nuc':
            plotting_key.plotting_options = [['Pink', 'DarkRed', 'Red', 'Red', 'Red', 'Orange', 'Gold', 'Yellow', 'LawnGreen', 'LawnGreen', 'LawnGreen', 'Green', 'OliveDrab', 'Cyan', 'Cyan', 'Cyan', 'Blue'],
                                    ('vifA1A7', 'tatA1A3A7', 'revA1A4aA7', 'revA1A4aA7', 'revA1A4aA7', 'nefA1A5A7', 'vprA2A7', 'tatA2A3A7', 'revA2A4aA7', 'revA2A4bA7', 'revA2A4cA7', 'nefA2A5A7', 'tatA3A7', 'revA4aA7', 'revA4bA7', 'revA4cA7', 'nefA5A7'), 
                                    'Multi-spliced transcript abundance in nucleus', 'Transcript abundance', 'on']
        elif plotting_key_name == 'full_len_transcripts_cyt':        
            plotting_key.plotting_options = [['Pink', 'Red', 'Orange', 'Yellow', 'LawnGreen', 'Green', 'Cyan', 'Blue', 'Purple'],
                                    ('X=0', '1', '2', '3', '4', '5', '6', '7', '8'), 
                                    'Full length transcript abundance in cytoplasm', 'Transcript abundance', 'on']
        elif plotting_key_name == 'multi_splice_transcript_cyt':        
            plotting_key.plotting_options = [['Pink', 'DarkRed', 'Red', 'Red', 'Red', 'Orange', 'Gold', 'Yellow', 'LawnGreen', 'LawnGreen', 'LawnGreen', 'Green', 'OliveDrab', 'Cyan', 'Cyan', 'Cyan', 'Blue'],
                                    ('vifA1A7', 'tatA1A3A7', 'revA1A4aA7', 'revA1A4aA7', 'revA1A4aA7', 'nefA1A5A7', 'vprA2A7', 'tatA2A3A7', 'revA2A4aA7', 'revA2A4bA7', 'revA2A4cA7', 'nefA2A5A7', 'tatA3A7', 'revA4aA7', 'revA4bA7', 'revA4cA7', 'nefA5A7'), 
                                    'Multi-spliced transcript abundance in cytoplasm', 'Transcript abundance', 'on']
        elif plotting_key_name == 'transcripts_synthesized':
            plotting_key.plotting_options = [['Blue'],
                                    ('Transcripts Synthesized',), 
                                    'Synthesized transcript abundance', 'Transcript abundance', 'on'] 
        elif plotting_key_name == 'proteins_nuc':
            plotting_key.plotting_options = [['Pink', 'Red', 'Orange', 'Yellow', 'Green', 'Blue', 'Purple'],
                                    ('Vif', 'Vpr', 'Tat', 'Env', 'Rev', 'Nef', 'Gag/Pol'), 
                                    'Proteins in the Nucleus', 'Protein abundance', 'on']
        elif plotting_key_name == 'proteins_cyt':        
            plotting_key.plotting_options = [['Pink', 'Red', 'Orange', 'Yellow', 'Green', 'Blue', 'Purple'],
                                    ('Vif', 'Vpr', 'Tat', 'Env', 'Rev', 'Nef', 'Gag/Pol'), 
                                    'Proteins in the Cytoplasm', 'Protein abundance', 'on']
        elif plotting_key_name == 'proteins_mem':        
            plotting_key.plotting_options = [['Pink', 'Red', 'Orange', 'Yellow', 'Green', 'Blue', 'Purple'],
                                    ('Vif', 'Vpr', 'Tat', 'Env', 'Rev', 'Nef', 'Gag/Pol'), 
                                    'Proteins in the Membrane', 'Protein abundance', 'on']
        elif plotting_key_name == 'proteins_virion':        
            plotting_key.plotting_options = [['Pink', 'Red', 'Orange', 'Yellow', 'Green', 'Blue', 'Purple'],
                                    ('Vif', 'Vpr', 'Tat', 'Env', 'Rev', 'Nef', 'Gag/Pol'), 
                                    'Virion Proteins', 'Protein abundance', 'on']

        elif plotting_key_name == 'Tat_derived_transcription_rate':
            plotting_key.plotting_options = [['Blue'],
                            ('Tat_derived_transcription_rate',), 
                            'Tat derived transcription rate', 'transcription_rate', 'on', False]

    # Start of special keys
        elif plotting_key_name == 'single_splice_transcript_nuc_1':
            plotting_key.dependencies = ['single_splice_transcript_nuc']

            def generate_key_data(input_dict):
                output_simulation_matrix = input_dict['single_splice_transcript_nuc'][np.arange(0,63,7)]
                return output_simulation_matrix

            plotting_key.plotting_options = [['Pink', 'Red', 'Orange', 'Yellow', 'LawnGreen', 'Green', 'Cyan', 'Blue', 'Purple'],
                                    [], 
                                    'Single spliced transcripts in Nucleus (vif)', 'vif', 'on']    
        elif plotting_key_name == 'single_splice_transcript_nuc_2':
            plotting_key.dependencies = ['single_splice_transcript_nuc']

            def generate_key_data(input_dict):
                output_simulation_matrix = input_dict['single_splice_transcript_nuc'][np.arange(1,63,7)]
                return output_simulation_matrix

            plotting_key.plotting_options = [['Pink', 'Red', 'Orange', 'Yellow', 'LawnGreen', 'Green', 'Cyan', 'Blue', 'Purple'],
                                    [], 
                                    'Single spliced transcripts in Nucleus (vpr)', 'vpr', 'on']   
        elif plotting_key_name == 'single_splice_transcript_nuc_3':
            plotting_key.dependencies = ['single_splice_transcript_nuc']

            def generate_key_data(input_dict):
                output_simulation_matrix = input_dict['single_splice_transcript_nuc'][np.arange(2,63,7)]
                return output_simulation_matrix

            plotting_key.plotting_options = [['Pink', 'Red', 'Orange', 'Yellow', 'LawnGreen', 'Green', 'Cyan', 'Blue', 'Purple'],
                                    [], 
                                    'Single spliced transcripts in Nucleus (tat)', 'tat', 'on']   
        elif plotting_key_name == 'single_splice_transcript_nuc_4':
            plotting_key.dependencies = ['single_splice_transcript_nuc']

            def generate_key_data(input_dict):
                output_simulation_matrix = input_dict['single_splice_transcript_nuc'][np.arange(3,63,7)]
                return output_simulation_matrix

            plotting_key.plotting_options = [['Pink', 'Red', 'Orange', 'Yellow', 'LawnGreen', 'Green', 'Cyan', 'Blue', 'Purple'],
                                    [], 
                                    'Single spliced transcripts in Nucleus (env4a)', 'env4a', 'on']   
        elif plotting_key_name == 'single_splice_transcript_nuc_5':
            plotting_key.dependencies = ['single_splice_transcript_nuc']

            def generate_key_data(input_dict):
                output_simulation_matrix = input_dict['single_splice_transcript_nuc'][np.arange(4,63,7)]
                return output_simulation_matrix

            plotting_key.plotting_options = [['Pink', 'Red', 'Orange', 'Yellow', 'LawnGreen', 'Green', 'Cyan', 'Blue', 'Purple'],
                                    [], 
                                    'Single spliced transcripts in Nucleus (env4b)', 'env4b', 'on']   
        elif plotting_key_name == 'single_splice_transcript_nuc_6':
            plotting_key.dependencies = ['single_splice_transcript_nuc']

            def generate_key_data(input_dict):
                output_simulation_matrix = input_dict['single_splice_transcript_nuc'][np.arange(5,63,7)]
                return output_simulation_matrix

            plotting_key.plotting_options = [['Pink', 'Red', 'Orange', 'Yellow', 'LawnGreen', 'Green', 'Cyan', 'Blue', 'Purple'],
                                    [], 
                                    'Single spliced transcripts in Nucleus (env4c)', 'env4c', 'on']   
        elif plotting_key_name == 'single_splice_transcript_nuc_7':
            plotting_key.dependencies = ['single_splice_transcript_nuc']

            def generate_key_data(input_dict):
                output_simulation_matrix = input_dict['single_splice_transcript_nuc'][np.arange(6,63,7)]
                return output_simulation_matrix

            plotting_key.plotting_options = [['Pink', 'Red', 'Orange', 'Yellow', 'LawnGreen', 'Green', 'Cyan', 'Blue', 'Purple'],
                                    [], 
                                    'Single spliced transcripts in Nucleus (env5)', 'env5', 'on']

        elif plotting_key_name == 'single_splice_transcript_cyt_1':
            plotting_key.dependencies = ['single_splice_transcript_cyt']

            def generate_key_data(input_dict):
                output_simulation_matrix = input_dict['single_splice_transcript_cyt'][np.arange(0,63,7)]
                return output_simulation_matrix

            plotting_key.plotting_options = [['Pink', 'Red', 'Orange', 'Yellow', 'LawnGreen', 'Green', 'Cyan', 'Blue', 'Purple'],
                                    [], 
                                    'Single spliced transcripts in Cytoplasm (vif)', 'vif', 'on']    
        elif plotting_key_name == 'single_splice_transcript_cyt_2':
            plotting_key.dependencies = ['single_splice_transcript_cyt']

            def generate_key_data(input_dict):
                output_simulation_matrix = input_dict['single_splice_transcript_cyt'][np.arange(1,63,7)]
                return output_simulation_matrix

            plotting_key.plotting_options = [['Pink', 'Red', 'Orange', 'Yellow', 'LawnGreen', 'Green', 'Cyan', 'Blue', 'Purple'],
                                    [], 
                                    'Single spliced transcripts in Cytoplasm (vpr)', 'vpr', 'on']   
        elif plotting_key_name == 'single_splice_transcript_cyt_3':
            plotting_key.dependencies = ['single_splice_transcript_cyt']

            def generate_key_data(input_dict):
                output_simulation_matrix = input_dict['single_splice_transcript_cyt'][np.arange(2,63,7)]
                return output_simulation_matrix

            plotting_key.plotting_options = [['Pink', 'Red', 'Orange', 'Yellow', 'LawnGreen', 'Green', 'Cyan', 'Blue', 'Purple'],
                                    [], 
                                    'Single spliced transcripts in Cytoplasm (tat)', 'tat', 'on']   
        elif plotting_key_name == 'single_splice_transcript_cyt_4':
            plotting_key.dependencies = ['single_splice_transcript_cyt']

            def generate_key_data(input_dict):
                output_simulation_matrix = input_dict['single_splice_transcript_cyt'][np.arange(3,63,7)]
                return output_simulation_matrix

            plotting_key.plotting_options = [['Pink', 'Red', 'Orange', 'Yellow', 'LawnGreen', 'Green', 'Cyan', 'Blue', 'Purple'],
                                    [], 
                                    'Single spliced transcripts in Cytoplasm (env4a)', 'env4a', 'on']   
        elif plotting_key_name == 'single_splice_transcript_cyt_5':
            plotting_key.dependencies = ['single_splice_transcript_cyt']

            def generate_key_data(input_dict):
                output_simulation_matrix = input_dict['single_splice_transcript_cyt'][np.arange(4,63,7)]
                return output_simulation_matrix

            plotting_key.plotting_options = [['Pink', 'Red', 'Orange', 'Yellow', 'LawnGreen', 'Green', 'Cyan', 'Blue', 'Purple'],
                                    [], 
                                    'Single spliced transcripts in Cytoplasm (env4b)', 'env4b', 'on']   
        elif plotting_key_name == 'single_splice_transcript_cyt_6':
            plotting_key.dependencies = ['single_splice_transcript_cyt']

            def generate_key_data(input_dict):
                output_simulation_matrix = input_dict['single_splice_transcript_cyt'][np.arange(5,63,7)]
                return output_simulation_matrix

            plotting_key.plotting_options = [['Pink', 'Red', 'Orange', 'Yellow', 'LawnGreen', 'Green', 'Cyan', 'Blue', 'Purple'],
                                    [], 
                                    'Single spliced transcripts in Cytoplasm (env4c)', 'env4c', 'on']   
        elif plotting_key_name == 'single_splice_transcript_cyt_7':
            plotting_key.dependencies = ['single_splice_transcript_cyt']

            def generate_key_data(input_dict):
                output_simulation_matrix = input_dict['single_splice_transcript_cyt'][np.arange(6,63,7)]
                return output_simulation_matrix

            plotting_key.plotting_options = [['Pink', 'Red', 'Orange', 'Yellow', 'LawnGreen', 'Green', 'Cyan', 'Blue', 'Purple'],
                                    [], 
                                    'Single spliced transcripts in Cytoplasm (env5)', 'env5', 'on']    

        elif plotting_key_name == 'total_proteins':
            plotting_key.dependencies = ['proteins_nuc', 'proteins_cyt', 'proteins_mem', 'proteins_virion']

            def generate_key_data(input_dict):
                output_simulation_matrix = input_dict['proteins_cyt'] + input_dict['proteins_nuc'] + input_dict['proteins_mem'] + input_dict['proteins_virion']
                return output_simulation_matrix

            plotting_key.plotting_options = [['Pink', 'Red', 'Orange', 'Yellow', 'Green', 'Blue', 'Purple', 'Cyan', 'Black'],
                                    ('Vif', 'Vpr', 'Tat', 'Env', 'Rev', 'Nef', 'Gag', 'Gag/Pro/Pol', 'Gag dimers'), 
                                    'Total Proteins in the Cell', 'Protein abundance in cell', 'on']

        elif plotting_key_name == 'total_env_various_places':
            plotting_key.dependencies = ['env_cyt']
            plotting_key.dependencies += ['env_ER', 'env_ER_G1', 'env_ER_G2', 'env_ER_G3', 'env_ER_G3_folded', 'env_ER_G4_folded']
            plotting_key.dependencies += ['env_Golgi', 'env_Golgi_G5', 'env_Golgi_G5_error', 'env_trimers', 'env_trimers_cleaved']
            plotting_key.dependencies += ['env_trimers_membrane']
            plotting_key.dependencies += ['total_num_of_virion_Env_trimer']

            def generate_key_data(input_dict):
                Env_cyt = input_dict['env_cyt']
                Env_endo_reticulum = input_dict['env_ER'] + input_dict['env_ER_G1'] + input_dict['env_ER_G2'] + input_dict['env_ER_G3'] + input_dict['env_ER_G3_folded'] + input_dict['env_ER_G4_folded']
                Env_golgi_body = input_dict['env_Golgi'] + input_dict['env_Golgi_G5'] + input_dict['env_Golgi_G5_error'] + 3*input_dict['env_trimers'].sum(axis=0) + 3*input_dict['env_trimers_cleaved'].sum(axis=0)
                Env_membrane = np.array([3 * input_dict['env_trimers_membrane'].sum(axis=0)])
                Env_viron = np.array([3 * input_dict['total_num_of_virion_Env_trimer'].sum(axis=0)])
                output_simulation_matrix = np.concatenate((Env_cyt, Env_endo_reticulum, Env_golgi_body, Env_membrane, Env_viron))
                return output_simulation_matrix

            plotting_key.plotting_options = [['Blue', 'Green', 'Red', 'Cyan', 'Purple'],
                                    ('cytoplasm', 'ER', 'golgi', 'membrane', 'virons'), 
                                    'Env compartment abundance', 'Number of Env (monomers)', 'off']

        elif plotting_key_name == 'all_ER_forms':
            plotting_key.dependencies = ['env_ER', 'env_ER_G1', 'env_ER_G2', 'env_ER_G3', 'env_ER_G3_folded', 'env_ER_G4_folded']

            def generate_key_data(input_dict):
                Env_ER = input_dict['env_ER']
                Env_ER_G1 = input_dict['env_ER_G1']
                Env_ER_G2 = input_dict['env_ER_G2']
                Env_ER_G3 = input_dict['env_ER_G3']
                Env_ER_G3_folded = input_dict['env_ER_G3_folded']
                Env_ER_G4_folded = input_dict['env_ER_G4_folded']
                output_simulation_matrix = np.concatenate((Env_ER, Env_ER_G1, Env_ER_G2, Env_ER_G3, Env_ER_G3_folded, Env_ER_G4_folded))
                return output_simulation_matrix

            plotting_key.plotting_options = [['Blue', 'Green', 'Red', 'Cyan', 'Purple', 'Yellow'],
                                    ('unglycosylated', 'post-Oligosaccharyltransferase', 'post-GlucosidaseI', 'post-GlucosidaseII', 'folded by calnexin', 'post-GlucosidaseII'), 
                                    'Env in different glycosylation states in the ER', 'Number of Env (monomers)', 'off']

        elif plotting_key_name == 'all_Golgi_forms':
            plotting_key.dependencies = ['env_Golgi', 'env_Golgi_G5', 'env_Golgi_G5_error', 'env_trimers', 'env_trimers_cleaved']

            def generate_key_data(input_dict):
                Env_Golgi = input_dict['env_Golgi']
                Env_Golgi_G5 = input_dict['env_Golgi_G5']
                Env_Golgi_G5_error = input_dict['env_Golgi_G5_error']
                Env_trimers = np.array([3 * input_dict['env_trimers'].sum(axis=0)])
                Env_trimers_cleaved = np.array([3 * input_dict['env_trimers_cleaved'].sum(axis=0)])
                output_simulation_matrix = np.concatenate((Env_Golgi, Env_Golgi_G5, Env_Golgi_G5_error, Env_trimers, Env_trimers_cleaved))
                return output_simulation_matrix

            plotting_key.plotting_options = [['Blue', 'Green', 'Red', 'Cyan', 'Purple'],
                                    ('post-ER', 'successful golgi glycos.', 'error in golgi glycos.', 'in trimer', 'in cleaved trimer'), 
                                    'Env in different states in the golgi', 'Number of Env (monomers)', 'off']

        elif plotting_key_name == 'proteins_in_nuc_and_cyt':
            plotting_key.dependencies = ['proteins_nuc', 'proteins_cyt']

            def generate_key_data(input_dict):
                output_simulation_matrix = input_dict['proteins_nuc']+input_dict['proteins_cyt']
                return output_simulation_matrix

            plotting_key.plotting_options = [['Pink', 'Red', 'Orange', 'Yellow', 'Green', 'Blue', 'Purple'],
                                    ('Vif', 'Vpr', 'Tat', 'Env', 'Rev', 'Nef', 'Gag/Pol'), 
                                    'Abundance of different proteins in nucleus and cytoplasm', 'Protein abundance', 'on']
        elif plotting_key_name == 'total proteins_nuc':
            plotting_key.dependencies = ['proteins_nuc']

            def generate_key_data(input_dict):
                output_simulation_matrix = np.array([input_dict['proteins_nuc'].sum(axis=0)])
                return output_simulation_matrix

            plotting_key.plotting_options = [['Black'],
                                    ('Total Proteins in Nuc',), 
                                    'Total Proteins in Nucleus', 'Protein abundance', 'on', True]
        elif plotting_key_name == 'total proteins_cyt':
            plotting_key.dependencies = ['proteins_cyt']

            def generate_key_data(input_dict):
                output_simulation_matrix = np.array([input_dict['proteins_cyt'].sum(axis=0)])
                return output_simulation_matrix

            plotting_key.plotting_options = [['Black'],
                                    ('Total Proteins in Cyt',), 
                                    'Total Proteins in Cytoplasm', 'Protein abundance', 'on', False]
        elif plotting_key_name == 'total single spliced mRNA nuc':
            plotting_key.dependencies = ['single_splice_transcript_nuc']

            def generate_key_data(input_dict):
                subset_single_splice_transcript_nuc = input_dict['single_splice_transcript_nuc']
                output_simulation_matrix = np.zeros([7,np.shape(subset_single_splice_transcript_nuc)[1]])
                for mRNA in range(7):
                    subset_single_splice_transcript_nuc_temp = subset_single_splice_transcript_nuc[np.arange(mRNA,63,7),:]
                    output_simulation_matrix[mRNA,:] = subset_single_splice_transcript_nuc_temp.sum(axis=0)    
                    return output_simulation_matrix

            plotting_key.plotting_options = [['Pink', 'Red', 'Orange', 'Yellow', 'Yellow', 'Yellow', 'Green'],
                                    ('Vif', 'Vpr', 'Tat', 'Env4', 'Env4', 'Env4', 'Env5'), 
                                    'Total single spliced mRNA nuc', 'mRNA abundance', 'on']
        elif plotting_key_name == 'total single spliced mRNA cyt':
            plotting_key.dependencies = ['single_splice_transcript_cyt']

            def generate_key_data(input_dict):
                subset_single_splice_transcript_cyt = input_dict['single_splice_transcript_cyt']
                output_simulation_matrix = np.zeros([7,np.shape(subset_single_splice_transcript_cyt)[1]])
                for mRNA in range(7):
                    subset_single_splice_transcript_cyt_temp = subset_single_splice_transcript_cyt[np.arange(mRNA,63,7),:]
                    output_simulation_matrix[mRNA,:] = subset_single_splice_transcript_cyt_temp.sum(axis=0)
                return output_simulation_matrix

            plotting_key.plotting_options = [['Pink', 'Red', 'Orange', 'Yellow', 'Yellow', 'Yellow', 'Green'],
                                    ('Vif', 'Vpr', 'Tat', 'Env4', 'Env4', 'Env4', 'Env5'), 
                                    'Total single spliced mRNA cyt', 'mRNA abundance', 'on']                                    
        elif plotting_key_name == 'total single spliced mRNA':
            plotting_key.dependencies = ['single_splice_transcript_nuc', 'single_splice_transcript_cyt']

            def generate_key_data(input_dict):
                subset_single_splice_transcript_nuc = input_dict['single_splice_transcript_nuc']
                subset_single_splice_transcript_cyt = input_dict['single_splice_transcript_cyt']
                output_simulation_matrix = np.zeros([7,np.shape(subset_single_splice_transcript_nuc)[1]])
                for mRNA in range(7):
                    subset_single_splice_transcript_nuc_temp = subset_single_splice_transcript_nuc[np.arange(mRNA,63,7),:]
                    subset_single_splice_transcript_cyt_temp = subset_single_splice_transcript_cyt[np.arange(mRNA,63,7),:]
                    output_simulation_matrix[mRNA,:] = subset_single_splice_transcript_nuc_temp.sum(axis=0) + subset_single_splice_transcript_cyt_temp.sum(axis=0)    
                    return output_simulation_matrix

            plotting_key.plotting_options = [['Pink', 'Red', 'Orange', 'Yellow', 'Yellow', 'Yellow', 'Green'],
                                    ('Vif', 'Vpr', 'Tat', 'Env4', 'Env4', 'Env4', 'Env5'), 
                                    'Total single spliced mRNA Nuc + Cyt', 'mRNA abundance', 'on']
        elif plotting_key_name == 'full single multi mRNA nuc':
            plotting_key.dependencies = ['single_splice_transcript_nuc', 'full_len_transcripts_nuc', 'multi_splice_transcript_nuc']

            def generate_key_data(input_dict):
                subset_single_splice_transcript_nuc = input_dict['single_splice_transcript_nuc']
                output_simulation_matrix = np.zeros([3,np.shape(subset_single_splice_transcript_nuc)[1]])
                #full
                output_simulation_matrix[0,:] = input_dict['full_len_transcripts_nuc'].sum(axis=0)
                 #Single
                output_simulation_matrix[1,:] = input_dict['single_splice_transcript_nuc'].sum(axis=0)
                 #multi
                output_simulation_matrix[2,:] = input_dict['multi_splice_transcript_nuc'].sum(axis=0)
                return output_simulation_matrix

            plotting_key.plotting_options = [['Red', 'Green', 'Blue'],
                                    ('full', 'single', 'multi'), 
                                    'full single multi mRNA nuc', 'mRNA abundance', 'on']
        elif plotting_key_name == 'full single multi mRNA cyt':
            plotting_key.dependencies = ['full_len_transcripts_cyt', 'single_splice_transcript_cyt', 'multi_splice_transcript_cyt']

            def generate_key_data(input_dict):
                subset_single_splice_transcript_cyt = input_dict['single_splice_transcript_cyt']
                output_simulation_matrix = np.zeros([3,np.shape(subset_single_splice_transcript_cyt)[1]])
                #full
                output_simulation_matrix[0,:] = input_dict['full_len_transcripts_cyt'].sum(axis=0)
                 #Single
                output_simulation_matrix[1,:] = input_dict['single_splice_transcript_cyt'].sum(axis=0)
                 #multi
                output_simulation_matrix[2,:] = input_dict['multi_splice_transcript_cyt'].sum(axis=0)
                return output_simulation_matrix

            plotting_key.plotting_options = [['Red', 'Green', 'Blue'],
                                    ('full', 'single', 'multi'), 
                                    'full single multi mRNA cyt', 'mRNA abundance', 'on']
        elif plotting_key_name == 'full single multi mRNA':
            plotting_key.dependencies = ['full_len_transcript_cyt', 'full_len_transcripts_nuc', 'single_splice_transcript_cyt', 'single_splice_transcript_nuc', 'multi_splice_transcript_cyt', 'multi_splice_transcript_nuc']

            def generate_key_data(input_dict):
                subset_single_splice_transcript_cyt = input_dict['single_splice_transcript_cyt']
                output_simulation_matrix = np.zeros([3,np.shape(subset_single_splice_transcript_cyt)[1]])
                #full
                output_simulation_matrix[0,:] = input_dict['full_len_transcripts_cyt'].sum(axis=0)+input_dict['full_len_transcripts_nuc'].sum(axis=0)
                 #Single
                output_simulation_matrix[1,:] = input_dict['single_splice_transcript_cyt'].sum(axis=0)+input_dict['single_splice_transcript_nuc'].sum(axis=0)
                 #multi
                output_simulation_matrix[2,:] = input_dict['multi_splice_transcript_cyt'].sum(axis=0)+input_dict['multi_splice_transcript_nuc'].sum(axis=0)
                return output_simulation_matrix

            plotting_key.plotting_options = [['Red', 'Green', 'Blue'],
                                    ('full', 'single', 'multi'), 
                                    'full single multi mRNA', 'mRNA abundance', 'on']
        elif plotting_key_name == 'mRNA by protein product nuc':
            plotting_key.dependencies = ['full_len_transcripts_nuc', 'single_splice_transcript_nuc', 'multi_splice_transcript_nuc']

            def generate_key_data(input_dict):
                subset_full_len_transcripts_nuc = input_dict['full_len_transcripts_nuc']
                subset_single_splice_transcript_nuc = input_dict['single_splice_transcript_nuc']
                subset_multi_splice_transcript_nuc = input_dict['multi_splice_transcript_nuc']           
                output_simulation_matrix = np.zeros([7,np.shape(subset_single_splice_transcript_nuc)[1]])
                #vif:
                output_simulation_matrix[0,:] = subset_single_splice_transcript_nuc[np.arange(0,63,7),:].sum(axis=0) + subset_multi_splice_transcript_nuc[0,:]
                #vpr:
                output_simulation_matrix[1,:] = subset_single_splice_transcript_nuc[np.arange(1,63,7),:].sum(axis=0) + subset_multi_splice_transcript_nuc[6,:]
                #tat:
                output_simulation_matrix[2,:] = subset_single_splice_transcript_nuc[np.arange(2,63,7),:].sum(axis=0) + subset_multi_splice_transcript_nuc[(1,7,12),:].sum(axis=0)
                #env:
                output_simulation_matrix[3,:] = subset_single_splice_transcript_nuc[np.arange(3,63,7),:].sum(axis=0) + subset_single_splice_transcript_nuc[np.arange(4,63,7),:].sum(axis=0) + subset_single_splice_transcript_nuc[np.arange(5,63,7),:].sum(axis=0) + subset_single_splice_transcript_nuc[np.arange(6,63,7),:].sum(axis=0)
                #rev:
                output_simulation_matrix[4,:] = subset_multi_splice_transcript_nuc[(2,3,4,8,9,10,13,14,15),:].sum(axis=0)
                #nef:
                output_simulation_matrix[5,:] = subset_multi_splice_transcript_nuc[(5,11,16),:].sum(axis=0)
                #gag/pol:
                output_simulation_matrix[6,:] = subset_full_len_transcripts_nuc.sum(axis=0)   
                return output_simulation_matrix

            plotting_key.plotting_options = [['Pink', 'Red', 'Orange', 'Yellow', 'Green', 'Blue', 'Purple'],
                                    ('vif', 'vpr', 'tat', 'env', 'rev', 'nef', 'gag/pol'), 
                                    'mRNA by protein product nuc', 'mRNA abundance', 'on']
        elif plotting_key_name == 'mRNA by protein product cyt':
            plotting_key.dependencies = ['full_len_transcripts_cyt', 'single_splice_transcript_cyt', 'multi_splice_transcript_cyt']

            def generate_key_data(input_dict):
                subset_full_len_transcripts_cyt = input_dict['full_len_transcripts_cyt']
                subset_single_splice_transcript_cyt = input_dict['single_splice_transcript_cyt']
                subset_multi_splice_transcript_cyt = input_dict['multi_splice_transcript_cyt']           
                output_simulation_matrix = np.zeros([7,np.shape(subset_single_splice_transcript_cyt)[1]])
                #vif:
                output_simulation_matrix[0,:] = subset_single_splice_transcript_cyt[np.arange(0,63,7),:].sum(axis=0) + subset_multi_splice_transcript_cyt[0,:]
                #vpr:
                output_simulation_matrix[1,:] = subset_single_splice_transcript_cyt[np.arange(1,63,7),:].sum(axis=0) + subset_multi_splice_transcript_cyt[6,:]
                #tat:
                output_simulation_matrix[2,:] = subset_single_splice_transcript_cyt[np.arange(2,63,7),:].sum(axis=0) + subset_multi_splice_transcript_cyt[(1,7,12),:].sum(axis=0)
                #env:
                output_simulation_matrix[3,:] = subset_single_splice_transcript_cyt[np.arange(3,63,7),:].sum(axis=0) + subset_single_splice_transcript_cyt[np.arange(4,63,7),:].sum(axis=0) + subset_single_splice_transcript_cyt[np.arange(5,63,7),:].sum(axis=0) + subset_single_splice_transcript_cyt[np.arange(6,63,7),:].sum(axis=0)
                #rev:
                output_simulation_matrix[4,:] = subset_multi_splice_transcript_cyt[(2,3,4,8,9,10,13,14,15),:].sum(axis=0)
                #nef:
                output_simulation_matrix[5,:] = subset_multi_splice_transcript_cyt[(5,11,16),:].sum(axis=0)
                #gag/pol:
                output_simulation_matrix[6,:] = subset_full_len_transcripts_cyt.sum(axis=0)   
                return output_simulation_matrix

            plotting_key.plotting_options = [['Pink', 'Red', 'Orange', 'Yellow', 'Green', 'Blue', 'Purple'],
                                    ('vif', 'vpr', 'tat', 'env', 'rev', 'nef', 'gag/pol'), 
                                    'mRNA by protein product cyt', 'mRNA abundance', 'on']
        elif plotting_key_name == 'mRNA by protein product':
            plotting_key.dependencies = ['full_len_transcript_cyt', 'full_len_transcripts_nuc', 'single_splice_transcript_cyt', 'single_splice_transcript_nuc', 'multi_splice_transcript_cyt', 'multi_splice_transcript_nuc']

            def generate_key_data(input_dict):
                subset_full_len_transcripts_nuc = input_dict['full_len_transcripts_nuc']
                subset_single_splice_transcript_nuc = input_dict['single_splice_transcript_nuc']
                subset_multi_splice_transcript_nuc = input_dict['multi_splice_transcript_nuc']
                subset_full_len_transcripts_cyt = input_dict['full_len_transcripts_cyt']
                subset_single_splice_transcript_cyt = input_dict['single_splice_transcript_cyt']
                subset_multi_splice_transcript_cyt = input_dict['multi_splice_transcript_cyt']           
                output_simulation_matrix = np.zeros([7,np.shape(subset_single_splice_transcript_nuc)[1]])
                #vif:
                output_simulation_matrix[0,:] = subset_single_splice_transcript_nuc[np.arange(0,63,7),:].sum(axis=0) + subset_multi_splice_transcript_nuc[0,:] + subset_single_splice_transcript_cyt[np.arange(0,63,7),:].sum(axis=0) + subset_multi_splice_transcript_cyt[0,:]
                #vpr:
                output_simulation_matrix[1,:] = subset_single_splice_transcript_nuc[np.arange(1,63,7),:].sum(axis=0) + subset_multi_splice_transcript_nuc[6,:] + subset_single_splice_transcript_cyt[np.arange(1,63,7),:].sum(axis=0) + subset_multi_splice_transcript_cyt[6,:]
                #tat:
                output_simulation_matrix[2,:] = subset_single_splice_transcript_nuc[np.arange(2,63,7),:].sum(axis=0) + subset_multi_splice_transcript_nuc[(1,7,12),:].sum(axis=0) + subset_single_splice_transcript_cyt[np.arange(2,63,7),:].sum(axis=0) + subset_multi_splice_transcript_cyt[(1,7,12),:].sum(axis=0)
                #env:
                output_simulation_matrix[3,:] = subset_single_splice_transcript_nuc[np.arange(3,63,7),:].sum(axis=0) + subset_single_splice_transcript_nuc[np.arange(4,63,7),:].sum(axis=0) + subset_single_splice_transcript_nuc[np.arange(5,63,7),:].sum(axis=0) + subset_single_splice_transcript_nuc[np.arange(6,63,7),:].sum(axis=0) + subset_single_splice_transcript_cyt[np.arange(3,63,7),:].sum(axis=0) + subset_single_splice_transcript_cyt[np.arange(4,63,7),:].sum(axis=0) + subset_single_splice_transcript_cyt[np.arange(5,63,7),:].sum(axis=0) + subset_single_splice_transcript_cyt[np.arange(6,63,7),:].sum(axis=0)
                #rev:
                output_simulation_matrix[4,:] = subset_multi_splice_transcript_nuc[(2,3,4,8,9,10,13,14,15),:].sum(axis=0) + subset_multi_splice_transcript_cyt[(2,3,4,8,9,10,13,14,15),:].sum(axis=0)
                #nef:
                output_simulation_matrix[5,:] = subset_multi_splice_transcript_nuc[(5,11,16),:].sum(axis=0) + subset_multi_splice_transcript_cyt[(5,11,16),:].sum(axis=0)
                #gag/pol:
                output_simulation_matrix[6,:] = subset_full_len_transcripts_nuc.sum(axis=0) + subset_full_len_transcripts_cyt.sum(axis=0)  
                return output_simulation_matrix

            plotting_key.plotting_options = [['Pink', 'Red', 'Orange', 'Yellow', 'Green', 'Blue', 'Purple'],
                                    ('vif', 'vpr', 'tat', 'env', 'rev', 'nef', 'gag/pol'), 
                                    'mRNA by protein product', 'mRNA abundance', 'on']
        elif plotting_key_name == 'num_virion_Gag_40_ge':
            plotting_key.plotting_options = [['Blue'], ('virion_num'), 
                                    'Number of Virion with >=40 percent expected num of Gag', 'Virion Quantity', 'off']
        elif plotting_key_name == 'progeny_count':
            plotting_key.plotting_options = [['Blue'], ('progeny_count',), 
                                    'Progeny Count', 'Progeny Count', 'off']
        elif plotting_key_name == 'total_num_of_virion_Gag':
            plotting_key.plotting_options = [['Blue'], ('total_num_of_virion_Gag',), 
                                    'Total Num of Virion Gag', 'Total Num of Virion Gag', 'off']

        elif plotting_key_name == 'total_num_of_virion_Vif':
            plotting_key.plotting_options = [['Blue'], ('total_num_of_virion_Vif',), 
                                    'Total Num of Virion Vif', 'Total Num of Virion Vif', 'off']        

        elif plotting_key_name == 'progeny_state_count':
            plotting_key.plotting_options = [['White', 'White', 'White', 'Blue', 'White'],
                                    ('NUCLEATE_CYT', 'NUCLEATE_MEM', 'GROWING_VIRION', 'VIRION_PREBUDDING', 'BUDDED_VIRION'), 
                                    'Counts of progeny in different states', 'Counts of progeny in different states', 'on']

        elif plotting_key_name == 'Prebudded_count':
            plotting_key.plotting_options = [['Blue'],
                                    ('VIRION_PREBUDDING'), 
                                    'Number of viral progeny', 'Number of viral progeny', 'on']

        elif plotting_key_name == 'state_of_diff_progeny':
            plotting_key.plotting_options = [['Blue']*10000, ('state_of_diff_progeny',), 
                                    'States of Different Progeny', 'States of Different Progeny', 'off', False]

        elif plotting_key_name == 'num_of_Gag_of_diff_progeny':
            plotting_key.plotting_options = [['Blue']*10000, ('num_of_Gag_of_diff_progeny',), 
                                    'Num of Gag of Different Progeny', 'Num of Gag of Different Progeny', 'off', False]

        elif plotting_key_name == 'num_of_Vif_of_diff_progeny':
            plotting_key.plotting_options = [['Blue']*10000, ('num_of_Vif_of_diff_progeny',), 
                                    'Num of Vif of Different Progeny', 'Num of Vif of Different Progeny', 'off', False]

        elif plotting_key_name == 'num_of_Env_t_of_diff_progeny':
            plotting_key.plotting_options = [['Blue']*10000, ('num_of_Env_t_of_diff_progeny',), 
                                    'Num of Env Trimer of Different Progeny', 'Num of Env Trimer of Different Progeny', 'off', False]

        elif plotting_key_name == 'total_num_of_virion_Env_trimer':
            plotting_key.plotting_options = [['Blue', 'Green', 'Red', 'Cyan'], ('0 correct', '1 correct', '2 correct', '3 correct'), 
                                    '# of virion trimers with 0,1,2, or 3 successfully glycosylated Env', 'Number of Env Trimers', 'off']                                        

        elif plotting_key_name == 'env_trimers':
            plotting_key.plotting_options = [['Blue', 'Green', 'Red', 'Cyan'], ('0 correct', '1 correct', '2 correct', '3 correct'), 
                                    '# of trimers with 0,1,2, or 3 successfully glycosylated Env', 'Number of Env Trimers', 'off']        

        elif plotting_key_name == 'env_trimers_cleaved':
            plotting_key.plotting_options = [['Blue', 'Green', 'Red', 'Cyan'], ('0 correct', '1 correct', '2 correct', '3 correct'), 
                                    '# of cleaved trimers with 0,1,2, or 3 successfully glycosylated Env', 'Number of Env Trimers', 'off']  

        elif plotting_key_name == 'env_trimers_membrane':
            plotting_key.plotting_options = [['Blue', 'Green', 'Red', 'Cyan'], ('0 correct', '1 correct', '2 correct', '3 correct'), 
                                    '# of membrane trimers with 0,1,2, or 3 successfully glycosylated Env', 'Number of Env Trimers', 'off']          

        try:
            plotting_key.generate_key_data = generate_key_data
        except NameError:
            pass

        list_of_plotting_keys.append(plotting_key)
        if len(plotting_key.dependencies) > 0:
            for dependency_name in plotting_key.dependencies:
                set_of_req_record_key_names.add(dependency_name)
            list_of_dependent_keys.append(plotting_key)
        else:
            set_of_req_record_key_names.add(plotting_key_name)

    return (list_of_plotting_keys, list_of_dependent_keys, set_of_req_record_key_names)
