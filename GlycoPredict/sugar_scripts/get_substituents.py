from rdkit import Chem
from rdkit.Chem import AllChem

def get_substituents(donor,lg, donor_type, only_c2 = True):
    if donor_type == 'pyranose':
        donor = Chem.MolFromSmarts(donor)
        donor.UpdatePropertyCache(strict=False)
        cleave_C2_template =AllChem.ReactionFromSmarts('[*:2222][C:222]1[C:333][C:444][C:555][O:666][C:111]1-{}>>[C:222]1[C:333][C:444][C:555][O:666][C:111]1-{}.[*:2222][#0]'.format(lg[1:],lg[1:]))
        cleaving_C2 = cleave_C2_template.RunReactants((donor,))
        if len(cleaving_C2 ) == 0:
            C2_sub = '*[H]'
        else:
            for i in range(len(cleaving_C2)):
                C2_sub = cleaving_C2[i][1]
                [a.SetAtomMapNum(0) for a in C2_sub.GetAtoms()]
                C2_sub = Chem.MolToSmiles(C2_sub)
                if C2_sub != '*[H]':
                    break


        if only_c2 == False:
            ##Cleave C3##
            cleave_C3_template =AllChem.ReactionFromSmarts('[C:222]1[C:333]([*:3333])[C:444][C:555][O:666][C:111]1-{}>>[C:222]1[C:333][C:444][C:555][O:666][C:111]1-{}.[*:3333][#0]'.format(lg[1:],lg[1:]))
            cleaving_C3 = cleave_C3_template.RunReactants((donor,))
            if len(cleaving_C3 ) == 0:
                C3_sub = '*[H]'
            else:
                for i in range(len(cleaving_C3)):
                    C3_sub = cleaving_C3[i][1]
                    [a.SetAtomMapNum(0) for a in C3_sub.GetAtoms()]
                    C3_sub = Chem.MolToSmiles(C3_sub)
                    if C3_sub != '*[H]':
                        break
            
            ##Cleave C4##
            cleave_C4_template =AllChem.ReactionFromSmarts('[C:222]1[C:333][C:444]([*:4444])[C:555][O:666][C:111]1-{}>>[C:222]1[C:333][C:444][C:555][O:666][C:111]1-{}.[*:4444][#0]'.format(lg[1:],lg[1:]))
            cleaving_C4 = cleave_C4_template.RunReactants((donor,))
            if len(cleaving_C4 ) == 0:
                C4_sub = '*[H]'
            else:
                for i in range(len(cleaving_C4)):
                    C4_sub = cleaving_C4[i][1]
                    [a.SetAtomMapNum(0) for a in C4_sub.GetAtoms()]
                    C4_sub = Chem.MolToSmiles(C4_sub)
                    if C4_sub != '*[H]':
                        break
            ##Cleave C5#
            cleave_C5_template =AllChem.ReactionFromSmarts('[C:222]1[C:333][C:444][C:555]([*:5555])[O:666][C:111]1-{}>>[C:222]1[C:333][C:444][C:555][O:666][C:111]1-{}.[*:5555][#0]'.format(lg[1:],lg[1:]))
            cleaving_C5 = cleave_C5_template.RunReactants((donor,))
            if len(cleaving_C5 ) == 0:
                C5_sub = '*[H]'
            else:
                for i in range(len(cleaving_C5)):
                    C5_sub = cleaving_C5[i][1]
                    [a.SetAtomMapNum(0) for a in C5_sub.GetAtoms()]
                    C5_sub = Chem.MolToSmiles(C5_sub)
                    if C5_sub != '*[H]':
                        break

        
    if donor_type == 'furanose':
        donor = Chem.MolFromSmarts(donor)
        donor.UpdatePropertyCache(strict=False)
        cleave_C2_template =AllChem.ReactionFromSmarts('[*:2222][C:222]1[C:333][C:444][O:666][C:111]1-{}>>[C:222]1[C:333][C:444][O:666][C:111]1-{}.[*:2222][#0]'.format(lg[1:],lg[1:]))
        cleaving_C2 = cleave_C2_template.RunReactants((donor,))
        if len(cleaving_C2 ) == 0:
            C2_sub = '*[H]'
        else:
            for i in range(len(cleaving_C2)):
                C2_sub = cleaving_C2[i][1]
                [a.SetAtomMapNum(0) for a in C2_sub.GetAtoms()]
                C2_sub = Chem.MolToSmiles(C2_sub)
                if C2_sub != '*[H]':
                    break


        if only_c2 == False:
            ##Cleave C3##
            cleave_C3_template =AllChem.ReactionFromSmarts('[C:222]1[C:333]([*:3333])[C:444][O:666][C:111]1-{}>>[C:222]1[C:333][C:444][O:666][C:111]1-{}.[*:3333][#0]'.format(lg[1:],lg[1:]))
            cleaving_C3 = cleave_C3_template.RunReactants((donor,))
            if len(cleaving_C3 ) == 0:
                C3_sub = '*[H]'
            else:
                for i in range(len(cleaving_C3)):
                    C3_sub = cleaving_C3[i][1]
                    [a.SetAtomMapNum(0) for a in C3_sub.GetAtoms()]
                    C3_sub = Chem.MolToSmiles(C3_sub)
                    if C3_sub != '*[H]':
                        break
            
            ##Cleave C4##
            cleave_C4_template =AllChem.ReactionFromSmarts('[C:222]1[C:333][C:444]([*:4444])[O:666][C:111]1-{}>>[C:222]1[C:333][C:444][O:666][C:111]1-{}.[*:4444][#0]'.format(lg[1:],lg[1:]))
            cleaving_C4 = cleave_C4_template.RunReactants((donor,))
            if len(cleaving_C4 ) == 0:
                C4_sub = '*[H]'
            else:
                for i in range(len(cleaving_C4)):
                    C4_sub = cleaving_C4[i][1]
                    [a.SetAtomMapNum(0) for a in C4_sub.GetAtoms()]
                    C4_sub = Chem.MolToSmiles(C4_sub)
                    if C4_sub != '*[H]':
                        break
            
            ##Set C-5 to None#
            C5_sub = 'None'


    if only_c2 == False:
        return C2_sub,C3_sub,C4_sub,C5_sub
    
    else:
        return C2_sub