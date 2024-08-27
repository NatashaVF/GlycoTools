import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import ast
###Encoded expert trees###


###Expert helper functions###

def check_if_D_or_L(donor_vector):

    ###Assign D or L according to LURD###
    if donor_vector[-1] != 0:
        last_stereocenter = donor_vector[-1] 
        if last_stereocenter == 1:
            sugar = 'D'

        elif last_stereocenter == -1:
            sugar = 'L'

        else:
            sugar = None
            print('D or L could not be determined')
       
    elif donor_vector[-1] == 0:
        for i in range(1,(len(donor_vector))):
            last_stereocenter = donor_vector[-i]
            if last_stereocenter != 0:
                break
        
        if last_stereocenter == 1:
            sugar = 'L'
        
        elif last_stereocenter == -1:
            sugar = 'D'

        else:
            sugar = None
            print('D or L could not be determined')
    
    return sugar

def check_if_c2_is_ester(lg, donor):
    #Create ester substructure with correct lg
    ester_sub_structure = Chem.MolFromSmiles('O=COC1CCCOC1-{}'.format(lg[1:]))
    return donor.HasSubstructMatch(ester_sub_structure)

def check_if_allo_gluco_gulo_galacto(sugar, donor_vector):
    if sugar == 'D':
        if donor_vector[1:] == [-1,-1,-1,1]: #D-allo
            result = True
        elif donor_vector[1:] == [-1,1,-1,1]: #D-gluco
            result = True
        elif donor_vector[1:] == [-1,-1,1,1]: #D-gulo
            result = True
        elif donor_vector[1:] == [-1,1,1,1]: #D-galacto
            result = True
        else:
            result = False

    if sugar == 'L':
        if donor_vector[1:] == [1,1,1,-1]: #D-allo
            result = True
        elif donor_vector[1:] == [1,-1,1,-1]: #D-gluco
            result = True
        elif donor_vector[1:] == [1,1,-1,-1]: #D-gulo
            result = True
        elif donor_vector[1:] == [1,-1,-1,-1]: #D-galacto
            result = True
        else:
            result = False
    
    return result

def check_if_mannose(donor_vector):
    if donor_vector[1:] == [1,1,-1,1]: #D-manno
        result = True
    elif donor_vector[1:] == [-1,-1,1,-1]: #L-manno
        result = True
    else:
        result = False

    return result

def check_for_benzylidene(donor,lg):
    benzylidene_sub_structure = Chem.MolFromSmarts('O1C(COC(C2=CC=CC=C2)O3)C3CCC1-{}'.format(lg[1:]))
    
    return donor.HasSubstructMatch(benzylidene_sub_structure)

def is_solvent_nitrile(solvent_smiles):
    if "C#N" in solvent_smiles:
        result = True
    else:
        result = False
    return result

def check_for_C2_sub(donor_vector,c2_sub,donor,lg):
    c2 = donor_vector[1]

    if c2 == 0 and c2_sub == '*[H]':
        O_sub = Chem.MolFromSmiles('OC1CCCOC1-{}'.format(lg[1:]))

        return donor.HasSubstructMatch(O_sub)
    
    else:
        return False

###Check for axial C-2###
def is_C2_axial(donor_vector):
    if donor_vector[1] == 1 and donor_vector[4] == 1:
        return True
    elif donor_vector[1] == -1 and donor_vector[4] == -1:
        return True
    elif donor_vector[1] == -1 and donor_vector[4] == 0 and donor_vector[3] == 1:
        return True
    elif donor_vector[1] == 1 and donor_vector[4] == 0 and donor_vector[3] == -1:
        return True
    else:
        return False

def check_for_ngp(lg, donor):
    #Create ester substructure with correct lg
    ngp_ester_structure = Chem.MolFromSmiles('O=COC1CCCOC1-{}'.format(lg[1:]))
    ngp_amide_structure = Chem.MolFromSmiles('O=CNC1CCCOC1-{}'.format(lg[1:]))
    if donor.HasSubstructMatch(ngp_ester_structure) == True:
        return True
    elif donor.HasSubstructMatch(ngp_amide_structure) == True:
        return True
    else:
        return False
    
###Expert 1###    
def expert1(donor, donor_vector, lg, solvent, get_path = False):
    ano_conf = None
    donor = Chem.MolFromSmarts(donor)
    donor.UpdatePropertyCache(strict=False)
    sugar = check_if_D_or_L(donor_vector)
    path = []
    ###Does the donor have an ester in the 2-position###
    result_1 = check_if_c2_is_ester(lg, donor)

    ###If True, check if allo, gluco, gulo, ###
    if result_1 == True:
        path.append('node2a')
        result_2a = check_if_allo_gluco_gulo_galacto(sugar,donor_vector)
        if sugar == 'L':
            if result_2a == True:
                path.append('leaf2beta')
                ano_conf = -1 #Beta for L
            elif result_2a == False:
                path.append('leaf2alpha')
                ano_conf = 1 #Alpha for L

        ##D-sugar
        else:
            if result_2a == True:
                path.append('leaf2beta')
                ano_conf = 1 #Beta for D
            elif result_2a == False:
                path.append('leaf2alpha')
                ano_conf = -1 #Alpha for D
    
    ###If ester in 2 position false check if mannose###
    else:
        result_2b = check_if_mannose(donor_vector)
        path.append('node2b')

        ###If mannoside check if donor has 4,6-benzylidene###
        if result_2b == True:
            path.append('node3a')
            result_3a = check_for_benzylidene(donor,lg)
            if sugar == 'L':
                if result_3a == True:
                    path.append('leaf3a_beta')
                    ano_conf = -1 #Beta for L
                elif result_3a == False:
                    path.append('leaf3a_alpha')
                    ano_conf = 1 #Alpha for L

            ##D-sugar
            else: 
                if result_3a == True:
                    ano_conf = 1 #Beta for D
                    path.append('leaf3a_beta')
                elif result_3a == False:
                    path.append('leaf3a_alpha')
                    ano_conf = -1 #Alpha for D
        
        ###If not mannoside check if solvent is a nitrile solvent###
        if result_2b == False:
            path.append('node3b')
            if pd.isna(solvent) == False:
                result_3b = is_solvent_nitrile(df['solvent_smile'][ind])
            else:
                result_3b = False ## False for no solvent

            if sugar == 'L':
                if result_3b == True:
                    path.append('leaf3b_beta')
                    ano_conf = -1 #Beta for L
                elif result_3b == False:
                    path.append('leaf3b_alpha')
                    ano_conf = 1 #Alpha for L

            ##D-sugar
            else: 
                if result_3b == True:
                    path.append('leaf3b_beta')
                    ano_conf = 1 #Beta for D
                elif result_3b == False:
                    path.append('leaf3b_alpha')
                    ano_conf = -1 #Alpha for D
    
    if get_path == True:
        return ano_conf, path
    else:
        return ano_conf

###Exper2 algorithm###
def expert2(donor, donor_vector, c2_sub, lg, solvent, get_path = False):
    ano_conf = None
    donor = Chem.MolFromSmarts(donor)
    donor.UpdatePropertyCache(strict=False)
    sugar = check_if_D_or_L(donor_vector)
    path = []

    ###Does the sugar have a C-2 substituent###
    result_1 = check_for_C2_sub(donor_vector, c2_sub, donor, lg)

    ###False equals 2-deoxy which equals axial###
    if result_1 == True:
        path.append('leaf1')
        if sugar == 'L':
            ano_conf = 1 #Alpha for L
        else:
            ano_conf = -1 #Alpha for D


    ###If C-2 substituent present###
    elif result_1 == False:
        path.append('node2')
        ###Check if 2-sub is axial###
        result_2 = is_C2_axial(donor_vector)

        ###If axial product configuration is axial###
        if result_2 == True:
            path.append('leaf2')
            if sugar == 'L':
                ano_conf = 1 #Alpha for L
            else:
                ano_conf = -1 #Alpha for D

        ###If C-2 not axial###
        if result_2 == False:
            path.append('node3')
            ###Check for NGP###
            result_3 = check_for_ngp(lg, donor)

            ###If true equatorial/beta###
            if result_3 == True:
                path.append('leaf3')
                if sugar == 'L':
                    ano_conf = -1 #beta for L
                else:
                    ano_conf = 1 #Beta for D

            elif result_3 == False:
                path.append('node4')
                if pd.isna(solvent) == False:
                    result_3 = is_solvent_nitrile(df['solvent_smile'][ind])
                else:
                        result_3 = False ## False for no solvent

                if sugar == 'L':
                    if result_3 == True:
                        path.append('leaf4a')
                        ano_conf = -1 #Beta for L
                    elif result_3 == False:
                        path.append('leaf4b')
                        ano_conf = 1 #Alpha for L

                ##D-sugar
                else: 
                    if result_3 == True:
                        path.append('leaf4a')
                        ano_conf = 1 #Beta for D
                    elif result_3 == False:
                        path.append('leaf4b')
                        ano_conf = -1 #Alpha for D
            
    if get_path == True:
        return ano_conf, path
    else:
        return ano_conf
                

