from sklearn.ensemble import RandomForestRegressor
import pandas as pd
import pickle
pd.options.mode.chained_assignment = None
from GlycoPredict.sugar_scripts import featurize_dataset
from GlycoPredict.sugar_scripts import get_anomeric_configuration
from rxnmapper import RXNMapper
import numpy as np
import ast
import rxn

###Predict major anomer product, if minor product is formed, and anomeric ratio###
def predict(input_path, save_name = 'filename', save = True, map_rxn = True, only_ratio = False):
    ###Prepare input data
    df_input = pd.read_csv(input_path)

    ####Map reaction if not mapped###
    if map_rxn == True:
        rxnmapper = RXNMapper()
        for ind in df_input.index:
            rxn = df_input['rxnsmiles'][ind]
            res = rxnmapper.get_attention_guided_atom_maps([rxn])
            df_input['rxnsmiles'][ind] = res[0]['mapped_rxn']

    ###Prepare input data###
    df = featurize_dataset.create_featurized_data_set(df_input, remove_ano_stereo_pro = False, binary_conf = False, split_glyc = True, save = False, get_one_hot = False, only_c2 = True, file_name = 'file', remove_none_major_products = False, create_activator_column = False)

    if only_ratio == False:
        ###Prepare prediction columns###
        df['Predicted Major Conf'] = ""
        df['Predicted Major Ano'] = ""
        df['Predicted Minor Producted'] = ""
        df['Predicted Ratio (%pos)'] = ""

        ###Load models###
        major_anomer_predicter = pickle.load(open("models/major_anomer_RF.pickle", "rb"))
        minor_anomer_predicter = pickle.load(open("models/minor_anomer_RF.pickle", "rb"))
        ratio_neg_predicter = pickle.load(open("models/ratio_neg_RF.pickle", "rb"))
        ratio_pos_predicter = pickle.load(open("models/ratio_pos_RF.pickle", "rb"))

        ###Make predictions for all entries###
        for ind in df.index:
            x = []
            x.append(np.concatenate([df['donor_conf'][ind],[df['temperature_C'][ind]],
                                        [df['donor_type_ohe_furanose'][ind]],[df['donor_type_ohe_pyranose'][ind]],[df['acceptor_type_ohe_N'][ind]],[df['acceptor_type_ohe_O'][ind]],[df['acceptor_type_ohe_S'][ind]],
                                        df['fp donor'][ind],df['fp acceptor'][ind],df['solvent_fp'][ind],df['activator_fp'][ind]]))
            predicted_major_conf = major_anomer_predicter.predict(x)
            ano, ano_conf = get_anomeric_configuration.get_anomer(predicted_major_conf,df['donor_conf'][ind],df['donor_type'][ind])
            df['Predicted Major Conf'][ind] = ano_conf
            df['Predicted Major Ano'][ind] = ano

            predicted_minor_conf = minor_anomer_predicter.predict(x)
            df['Predicted Minor Producted'][ind] = predicted_minor_conf

            if predicted_minor_conf == 0 and predicted_major_conf == 1:
                df['Predicted Ratio (%pos)'][ind] = 100
            
            elif predicted_minor_conf == 0 and predicted_major_conf == 0:
                df['Predicted Ratio (%pos)'][ind] = 0
            
            elif predicted_minor_conf == 1:
                pred_neg = ratio_neg_predicter.predict(x)
                pred_pos = ratio_pos_predicter.predict(x)
                df['Predicted Ratio (%pos)'][ind] = ((100 - pred_neg) + pred_pos) / 2
    
    ###Option to only predict anomeric ratio###
    elif only_ratio == True:
        ###Add column###
        df['Predicted Ratio (%pos)'] = ""
        
        ###Load models###
        ratio_neg_predicter = pickle.load(open("models/ratio_neg_RF.pickle", "rb"))
        ratio_pos_predicter = pickle.load(open("models/ratio_pos_RF.pickle", "rb"))

        ###Prepare data###
        for ind in df.index:
            x = []
            x.append(np.concatenate([df['donor_conf'][ind],[df['temperature_C'][ind]],
                                        [df['donor_type_ohe_furanose'][ind]],[df['donor_type_ohe_pyranose'][ind]],[df['acceptor_type_ohe_N'][ind]],[df['acceptor_type_ohe_O'][ind]],[df['acceptor_type_ohe_S'][ind]],
                                        df['fp donor'][ind],df['fp acceptor'][ind],df['solvent_fp'][ind],df['activator_fp'][ind]]))
            ###Predict###
            pred_neg = ratio_neg_predicter.predict(x)
            pred_pos = ratio_pos_predicter.predict(x)
            df['Predicted Ratio (%pos)'][ind] = ((100 - pred_neg) + pred_pos) / 2

    ###Save predictions as CSV###
    if save == True:
        df.to_csv('results/'f'{save_name}_results.csv')
    
    return df
            


