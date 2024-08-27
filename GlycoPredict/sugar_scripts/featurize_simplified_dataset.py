import pandas as pd
import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import ast
from GlycoPredict.sugar_scripts.determine_sugar_configuration import get_sugar_conf
from GlycoPredict.helper_scripts.format_strings import format_str, fix_smiles
pd.options.mode.chained_assignment = None
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

###Function that create featurized data for the simplified model###

def create_featurized_data_set(df, binary_conf = False, save = False, get_one_hot = False,  file_name = 'file'):
    df = df.reset_index(drop=True)
    #Add new columns
    df['lg_fp'] = ""
    df['c2_fp'] = ""
    df['activator'] = ""
    df['activator_fp']=""
    df['solvent_fp'] = ""

    ohe_names = ['donor_type_ohe_furanose','donor_type_ohe_pyranose','acceptor_type_ohe_N','acceptor_type_ohe_O','acceptor_type_ohe_S']
    for name in ohe_names:
        if name in df.columns:
            df = df.drop(columns=[name])
    
    if get_one_hot == True:
        df['activator_list'] = ""

    ###Loop over df####

    for ind in df.index:

        donor_type = df['donor_type'][ind]
        acceptor_type = df['acceptor_type'][ind]

        ###Set missing temperatures to room temperature
        try:
            if np.isnan(float(df['temperature_C'][ind])) == True:
                df['temperature_C'][ind] = 20.123
        except:
            print(ind)

        #####Calculate fps#####
        ###Solvent###
        if pd.isna(df['solvent_smiles'][ind]) == False:
            solvents = (df['solvent_smiles'][ind])
            ###Get sane smiles for BF4- and PF6-
            solvents = solvents.replace('[F-][B+3]([F-])([F-])[F-]', 'F[B-](F)(F)F')
            solvents = solvents.replace('[F-][P+5]([F-])([F-])([F-])([F-])[F-]', '[F][P-]([F])([F])([F])([F])[F]')

            solvents =  ast.literal_eval(solvents)

            mf_solvents_sum = np.zeros(2048)
            for solvent in solvents:
                m = Chem.MolFromSmiles(solvent)
                try:
                    fp = np.array(Chem.AllChem.GetMorganFingerprintAsBitVect(m,radius=2, nBits=2048))
                except:
                    print('fingerprint not found for')
                    #donor_fp = None
                    print(solvent)
                    continue
                mf_solvents_sum = mf_solvents_sum + fp
            df['solvent_fp'][ind] = mf_solvents_sum.tolist()
        if pd.isna(df['solvent_smiles'][ind]) == True:
            df['solvent_fp'][ind] = np.zeros(2048).tolist()


        #df['activator'][ind] = df['activator_smiles'][ind]

        if pd.isna(df['activator'][ind]) == False:
            activators = (df['activator'][ind])

            m = Chem.MolFromSmiles(activators)
            try:
                activator_fp = np.array(Chem.AllChem.GetMorganFingerprintAsBitVect(m,radius=2, nBits=2048))
            except:
                print('fingerprint not found for activator:')
                print(activators)
                continue

            df['activator_fp'][ind] = activator_fp.tolist()
        if pd.isna(df['activator'][ind]) == True:
            df['activator_fp'][ind] = np.zeros(2048).tolist()         
            
            
        smiles_lg = df['lg'][ind]               
        lg = Chem.MolFromSmiles(smiles_lg)
        try:
            lg_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(lg,radius=2, nBits=2048)
        except:
            print('argument error lg. Removing:')
            print(ind)
            df = df.drop([ind])
            continue
        df['lg_fp'][ind] = np.array(lg_fp).tolist()
            
        smiles_c2 =df['C2_smiles'][ind]
        c2 = Chem.MolFromSmiles(smiles_c2)
        try:
            c2_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(c2,radius=2, nBits=2048)
        except:
            print('argument error c2. Removing:')
            print(ind)
            df = df.drop([ind])
            continue
        df['c2_fp'][ind] = np.array(c2_fp).tolist()

    if binary_conf == True:
        df['product_conf'] = df['product_conf'].map({-1: 0, 1: 1})
    
    ##OHE donor_type#
    df['donor_type_ohe'] = df['donor_type']
    df = pd.get_dummies(df, columns=['donor_type_ohe', ],dtype=int)
    print('Added OHE donor columns: 2')

    ##OHE acceptor_type#
    df['acceptor_type_ohe'] = df['acceptor_type']
    df = pd.get_dummies(df, columns=['acceptor_type_ohe', ],dtype=int)
    print('Added OHE acceptor columns: 3')

    ##Add OHE if not present i data set
    ohe_names = ['donor_type_ohe_furanose','donor_type_ohe_pyranose','acceptor_type_ohe_N','acceptor_type_ohe_O','acceptor_type_ohe_S']
    for name in ohe_names:
        if name not in df.columns:
            df[name] = 0

    
    if get_one_hot == True:

        ##OHE solvent##
        old_len = df.shape[1]
        df = df.join(df.solvent_cas.str.join('|').str.get_dummies())
        print('Added OHE solvent columns:'f"{df.shape[1]-old_len}")

        ##OHE activators#
        old_len = df.shape[1]
        df = df.join(df.activator_list.str.join('|').str.get_dummies())
        print('Added OHE activators columns:'f"{df.shape[1]-old_len}")


        ##OHE lg##
        df['OHE_lg'] = df['lg']
        old_len = df.shape[1]
        df = pd.get_dummies(df, columns=['OHE_lg', ],dtype=int)
        print('Added OHE lg columns:'f"{df.shape[1]-old_len}")

        
        ##OHE C-2##
        df['OHE_C2'] = df['C2_smiles']
        old_len = df.shape[1]
        df = pd.get_dummies(df, columns=['OHE_C2', ],dtype=int)
        print('Added OHE C-2 columns:'f"{df.shape[1]-old_len}")

    if save == True:
        df.to_csv(f"{file_name}.csv")
        
    return df
            

    
