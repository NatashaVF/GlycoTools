from rdkit import Chem
from rdkit.Chem import MCS

####Function that gives the donor vector for pyranosides####
def get_sugar_conf(donor,lg):
    conf_vector = []
    donor = Chem.MolFromSmarts(donor)
    donor.UpdatePropertyCache(strict=False)
    ###Get amomeric configuration###
    if donor.HasSubstructMatch(Chem.MolFromSmarts('[C:2]1[C:3][C:4][C:5][O:6][C@H:1]1-{}'.format(lg[1:])),useChirality=True) == True:
        conf_vector.append(1)
        ano_smart = '[C@H:1]'
    elif donor.HasSubstructMatch(Chem.MolFromSmarts('[C:2]1[C:3][C:4][C:5][O:6][C@@H:1]1-{}'.format(lg[1:])),useChirality=True) == True:
        conf_vector.append(-1)
        ano_smart = '[C@@H:1]'
    else:
        conf_vector.append(0)
        ano_smart = "[C:1]"

    ###Get configuration at C-2###
    if donor.HasSubstructMatch(Chem.MolFromSmarts('[C,N,Se,Br,F,I,O,S][C@H:2]1[C:3][C:4][C:5][O:6]{}1-{}'.format(ano_smart, lg[1:])),useChirality=True) == True:
        conf_vector.append(1)
    elif donor.HasSubstructMatch(Chem.MolFromSmarts('[C,N,Se,Br,F,I,O,S][C@@H:2]1[C:3][C:4][C:5][O:6]{}1-{}'.format(ano_smart,lg[1:])),useChirality=True) == True:
        conf_vector.append(-1)
    else:
        conf_vector.append(0)

    ###Get configuration at C-3###
    if donor.HasSubstructMatch(Chem.MolFromSmarts('[C:2]1[C@@H:3]([C,N,Se,Br,F,I,O,S])[C:4][C:5][O:6]{}1-{}'.format(ano_smart, lg[1:])),useChirality=True) == True:
        conf_vector.append(1)
    elif donor.HasSubstructMatch(Chem.MolFromSmarts('[C:2]1[C@H:3]([C,N,Se,Br,F,I,O,S])[C:4][C:5][O:6]{}1-{}'.format(ano_smart,lg[1:])),useChirality=True) == True:
        conf_vector.append(-1)
    else:
        conf_vector.append(0)


    ###Get configuration at C-4###
    if donor.HasSubstructMatch(Chem.MolFromSmarts('[C:2]1[C:3][C@H:4]([C,N,Se,Br,F,I,O,S])[C:5][O:6]{}1-{}'.format(ano_smart,lg[1:])),useChirality=True) == True:
        conf_vector.append(-1)
    elif donor.HasSubstructMatch(Chem.MolFromSmarts('[C:2]1[C:3][C@@H:4]([C,N,Se,Br,F,I,O,S])[C:5][O:6]{}1-{}'.format(ano_smart, lg[1:])),useChirality=True) == True:
        conf_vector.append(1)
    else:
        conf_vector.append(0)

    ###Get configuration at C-5###
    if donor.HasSubstructMatch(Chem.MolFromSmarts('[C:2]1[C:3][C:4][C@@H:5]([C,N,Se,Br,F,I,O,S])[O:6]{}1-{}'.format(ano_smart,lg[1:])),useChirality=True) == True:
        conf_vector.append(1)
    elif donor.HasSubstructMatch(Chem.MolFromSmarts('[C:2]1[C:3][C:4][C@H:5]([C,N,Se,Br,F,I,O,S])[O:6]{}1-{}'.format(ano_smart, lg[1:])),useChirality=True) == True:
        conf_vector.append(-1)
    else:
        conf_vector.append(0)

    return conf_vector