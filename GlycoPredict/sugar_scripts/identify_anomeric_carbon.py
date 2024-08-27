import pandas as pd
import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import ast
from rdkit.Chem import MCS

def identify_anomeric_carbon_idx(donor, match):
    for atom in donor.GetAtoms():
        if atom.GetAtomMapNum() == 0:
            continue
        if atom.GetIdx() in match:
            anomeric_carbon_idx = atom.GetAtomMapNum()

    return anomeric_carbon_idx

def get_anomeric_carbon_map_no(rxn,lg,donor):
    rxn = AllChem.ReactionFromSmarts(rxn)
    pattern_str = 'C' + lg
    pattern = Chem.MolFromSmarts(pattern_str)
    donor_old = Chem.MolFromSmarts(donor)
    donor_old.UpdatePropertyCache(strict=False)
    donor = Chem.Mol(rxn.GetReactantTemplate(0))
    donor.UpdatePropertyCache(strict=False)
    match = donor.GetSubstructMatch(pattern)
    if donor.HasSubstructMatch(donor_old) == True:
        pass
    else:
        donor = Chem.Mol(rxn.GetReactantTemplate(1))
        donor.UpdatePropertyCache(strict=False)
        match = donor.GetSubstructMatch(pattern)
    
    ano = identify_anomeric_carbon_idx(donor,match)
    return ano

def remove_ano_stereo_product(rxn,ano):
    rxn = rxn.replace('C@@H:'f"{ano}", 'CH:'f"{ano}")
    rxn = rxn.replace('C@H:'f"{ano}", 'CH:'f"{ano}")
    rxn = rxn.replace('C@@:'f"{ano}", 'C:'f"{ano}")
    rxn = rxn.replace('C@:'f"{ano}", 'C:'f"{ano}")
    return rxn

def inverse_anomeric_carbon(rxn,ano):
    if 'C@@H:'f"{ano}" in rxn:
        rxn = rxn.replace('C@@H:'f"{ano}", 'C@H:'f"{ano}")
    elif 'C@H:'f"{ano}" in rxn:
        rxn = rxn.replace('C@H:'f"{ano}", 'C@@H:'f"{ano}")
    elif 'C@@:'f"{ano}" in rxn:
        rxn = rxn.replace('C@@:'f"{ano}", 'C@:'f"{ano}")
    elif 'C@:'f"{ano}" in rxn:
        rxn = rxn.replace('C@:'f"{ano}", 'C@@:'f"{ano}")
    return rxn
