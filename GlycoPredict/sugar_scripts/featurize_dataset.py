import pandas as pd
import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import ast
from rdkit.Chem import MCS
from sugar_scripts.determine_sugar_configuration import get_sugar_conf
from sugar_scripts.determine_product_configuration import get_product_conf
from sugar_scripts.get_glyc_parts import get_parts
from sugar_scripts.get_substituents import get_substituents
from sugar_scripts.identify_anomeric_carbon import get_anomeric_carbon_map_no
from sugar_scripts.identify_anomeric_carbon import remove_ano_stereo_product
from helper_scripts.format_strings import format_str, fix_smiles, fix_rxns_failing_in_chemprop
from helper_scripts.remove_none_major_anomer import remove_none_major_anomer
pd.options.mode.chained_assignment = None
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

###Function that creates featurized dataset###

def create_featurized_data_set(df, remove_ano_stereo_pro = False, binary_conf = False, split_glyc = False, save = False, get_one_hot = False, only_c2 = True, file_name = 'file', remove_none_major_products = False, create_activator_column = True):
    ###df with mapped rxn smiles and parts###
    df = df.reset_index(drop=True)
    #Add new columns
    df['donor_conf'] = ""
    #df['product_anomer'] = ""
    #df['product_conf'] = ""
    df['activator'] = ""
    df['activator_fp']=""
    df['solvent_fp'] = ""
    df['fp donor'] = ""
    df['fp acceptor'] = ""


    if split_glyc == True:
        df['donor'] = ""
        df['lg']=""
        df['acceptor']=""
        df['product']=""
        df['C2_smiles'] = ""
    
    #if get_one_hot == True:
    df['activator_list'] = ""
        
    
    if only_c2 == False:
        df['C3_smiles'] = ""
        df['C4_smiles'] = ""
        df['C5_smiles'] = ""

    ###Loop over df####

    for ind in df.index:

        donor_type = df['donor_type'][ind]
        acceptor_type = df['acceptor_type'][ind]

        ###Set missing temperatures to room temperature
        if np.isnan(df['temperature_C'][ind]) == True:
            df['temperature_C'][ind] = 20.123
    

        if split_glyc == True:
            rxn = df['rxnsmiles'][ind]
            lg, donor, acceptor,product = get_parts(rxn, donor_type, acceptor_type)
            ###Remove glycals, inverse, simultanuos deprotection,etc###
            if lg == 'NO LG':
                df = df.drop([ind])
                print('No LG removing: 'f'{ind}')
                continue
            elif lg == 'NO PSEUDOPRODUCT':
                lg, donor, acceptor, product = get_parts(rxn, donor_type, acceptor_type, remove_salts=True)
            if lg == 'NO PSEUDOPRODUCT':
                df = df.drop([ind])
                print('No pseudoproduct removing: 'f'{ind}')
                continue
            if lg != 'NO PSEUDOPRODUCT' and lg != 'NO LG':      
                df['lg'][ind]=lg
                df['donor'][ind]=Chem.MolToSmiles(donor)
                df['acceptor'][ind]=Chem.MolToSmiles(acceptor)
                df['product'][ind]=Chem.MolToSmiles(product)
                df['C2_smiles'][ind] = get_substituents(df['donor'][ind] ,df['lg'][ind], donor_type)

    
        ####Calculate vector####
        donor = df['donor'][ind]
        lg =df['lg'][ind]
        v = get_sugar_conf(donor,lg, donor_type)
        df['donor_conf'][ind] = v

        
        ###Remove anomeric configuration in product####
        if remove_ano_stereo_pro == True:
            ano = get_anomeric_carbon_map_no(df['rxnsmiles'][ind],df['lg'][ind][1:],df['donor'][ind])
            smiles = df['rxnsmiles'][ind]
            #smiles = fix_rxns_failing_in_chemprop(smiles)
            rxn_split = smiles.split('>>')
            new_product = remove_ano_stereo_product(rxn_split[1],ano)
            df['rxnsmiles'][ind] = rxn_split[0] + '>>' + new_product

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

        #if pd.isna(df['solvent_name'][ind]) == False:
            #df['solvent_name'][ind] = ast.literal_eval(df['solvent_name'][ind])
        #if pd.isna(df['solvent_cas'][ind]) == False:
            #df['solvent_cas'][ind] = ast.literal_eval(df['solvent_cas'][ind])


        ####Activators#####
        if create_activator_column == True:
            #Make activator column###
            try:
                reagents = ast.literal_eval(df['reagent_smiles'][ind])
            except ValueError:
                reagents = []
            try:
                cats = ast.literal_eval(df['catalyst_smiles'][ind])
            except  ValueError:
                cats = []
            activators = reagents + cats
            df['activator_list'][ind] = activators

            sm = str(activators)

            if len(sm) != 0:
                sm = format_str(sm)
                sm = fix_smiles(sm)
            
            #if pd.isna(sm) == False:
                #m = Chem.MolFromSmiles(sm)
                #if m == None:
                    #print(ind)
                    #print(sm)
                df['activator'][ind] = sm
        
        else:
            df['activator'][ind] = df['activator_smiles'][ind]

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
            
            
        smiles_donor =df['donor'][ind]               
        donor = Chem.MolFromSmiles(smiles_donor)
        try:
            donor_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(donor,radius=2, nBits=2048)
        except:
            print('argument error donor. Removing:')
            print(ind)
            df = df.drop([ind])
            continue
        df['fp donor'][ind] = np.array(donor_fp).tolist()
            
        smiles_acceptor =df['acceptor'][ind]
        acceptor = Chem.MolFromSmiles(smiles_acceptor)
        try:
            acceptor_fp = Chem.AllChem.GetMorganFingerprintAsBitVect(acceptor,radius=2, nBits=2048)
        except:
            print('argument error acceptor. Removing:')
                #acceptor_fp = None
            print(ind)
            df = df.drop([ind])
            continue
        df['fp acceptor'][ind] = np.array(acceptor_fp).tolist()


        if only_c2 == False:
            c2, c3, c4, c5 = get_substituents(df['donor'][ind] ,df['lg'][ind], donor_type, only_c2 = False)
            df['C3_smiles'][ind] = c3
            df['C4_smiles'][ind] = c4
            df['C5_smiles'][ind] = c5
    
    if remove_none_major_products == True:
        df = df.reset_index(drop=True)
        print(len(df))
        df = remove_none_major_anomer(df)
        print('left after removing minor')
        print(len(df))
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

        if only_c2 == False:

            df['OHE_C3'] = df['C3_smiles']
            old_len = df.shape[1]
            df = pd.get_dummies(df, columns=['OHE_C3', ],dtype=int)
            print('Added OHE C-3 columns'f"{df.shape[1]-old_len}")

            df['OHE_C4'] = df['C4_smiles']
            old_len = df.shape[1]
            df = pd.get_dummies(df, columns=['OHE_C4', ],dtype=int)
            print('Added OHE C-4 columns:'f"{df.shape[1]-old_len}")


            old_len = df.shape[1]
            df['OHE_C5'] = df['C5_smiles']
            df = pd.get_dummies(df, columns=['OHE_C5', ],dtype=int)
            print('Added OHE C-5 columns:'f"{df.shape[1]-old_len}")



        

    if save == True:
        df.to_csv(f"{file_name}.csv")
        
    return df
            

    
