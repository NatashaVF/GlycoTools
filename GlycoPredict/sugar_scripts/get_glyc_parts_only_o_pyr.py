from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.SaltRemover import SaltRemover

def get_parts(rxn, remove_salts = False):
    
    rxn_smart = rxn
    rxn = AllChem.ReactionFromSmarts(rxn)
    product = Chem.Mol(rxn.GetProductTemplate(0))
    product.UpdatePropertyCache(strict=False)
    if remove_salts == True:
        remover = SaltRemover()
        product = remover.StripMol(product)
    pyr_temp = Chem.MolFromSmiles('C1CCCCO1')
    fur_temp = Chem.MolFromSmiles('C1CCCO1')
    reactant_1 = Chem.Mol(rxn.GetReactantTemplate(0))
    reactant_1.UpdatePropertyCache(strict=False)
    reactant_2 = Chem.Mol(rxn.GetReactantTemplate(1))
    reactant_2.UpdatePropertyCache(strict=False)
    donor = None
    acceptor = None

    o_glyc_temp_1 = AllChem.ReactionFromSmarts("[*:8][C:1]1[C:2][C:3][C:4][C:5][O:6]1.[*:9][O:7]>>[*:9][O:7][C:1]1[C:2][C:3][C:4][C:5][O:6]1.[*:8][#0]")
    o_glyc_temp_2 = AllChem.ReactionFromSmarts("[*:9][O:7].[*:8][C:1]1[C:2][C:3][C:4][C:5][O:6]1>>[*:9][O:7][C:1]1[C:2][C:3][C:4][C:5][O:6]1.[*:8][#0]")
    o_glyc_temp_reverse = AllChem.ReactionFromSmarts("[*:9][O:7][C:1]1[C:2][C:3][C:4][C:5][O:6]1.[*:8][#0]>>[*:8][C:1]1[C:2][C:3][C:4][C:5][O:6]1.[*:9][O:7]")


    glyc = o_glyc_temp_1.RunReactants((reactant_1,reactant_2))

    if len(glyc) == 0:
        donor = None


    elif len(glyc) != 0:


        pseudo_product = None
        for i in range(len(glyc)):
            test_product = glyc[i][0]
            test_product.UpdatePropertyCache(strict=False)
            product.UpdatePropertyCache(strict=False)
            if test_product.HasSubstructMatch(product) == True:
                pseudo_product = glyc[i][0]
                pseudo_product.UpdatePropertyCache(strict=False)
                lg = glyc[i][1]
                break
        
        if pseudo_product != None:

            glyc_reverse = o_glyc_temp_reverse.RunReactants((pseudo_product,lg))
            pseudo_donor = glyc_reverse[0][0]
            pseudo_donor.UpdatePropertyCache(strict=False)

            ####Sanity check####
            if reactant_1.HasSubstructMatch(pseudo_donor):
                donor = reactant_1
                acceptor = reactant_2
                [a.SetAtomMapNum(0) for a in lg.GetAtoms()]
                lg_smile = AllChem.MolToSmiles(lg)

            elif reactant_2.HasSubstructMatch(pseudo_donor):
                donor = reactant_2
                acceptor = reactant_1
                [a.SetAtomMapNum(0) for a in lg.GetAtoms()]
                lg_smile = AllChem.MolToSmiles(lg)

            else:
                donor = None
                #print('match none')

        else: 
            donor = None
        
    if donor == None:
        #print('temp2')
        glyc = o_glyc_temp_2.RunReactants((reactant_1,reactant_2))
        if len(glyc) == 0:
                lg_smile = 'NO LG'
        else:        
            pseudo_product = None
            for i in range(len(glyc)):
                test_product = glyc[i][0]
                test_product.UpdatePropertyCache(strict=False)
                product.UpdatePropertyCache(strict=False)
                if test_product.HasSubstructMatch(product) == True:
                    pseudo_product = glyc[i][0]
                    pseudo_product.UpdatePropertyCache(strict=False)
                    lg = glyc[i][1]
                    break

            if pseudo_product != None:          
                glyc_reverse = o_glyc_temp_reverse.RunReactants((pseudo_product,lg))
                pseudo_donor = glyc_reverse[0][0]
                pseudo_donor.UpdatePropertyCache(strict=False)

                if reactant_1.HasSubstructMatch(pseudo_donor):
                    donor = reactant_1
                    acceptor = reactant_2
                    [a.SetAtomMapNum(0) for a in lg.GetAtoms()]
                    lg_smile = AllChem.MolToSmiles(lg)

                elif reactant_2.HasSubstructMatch(pseudo_donor):
                    donor = reactant_2
                    acceptor = reactant_1
                    [a.SetAtomMapNum(0) for a in lg.GetAtoms()]
                    lg_smile = AllChem.MolToSmiles(lg)

                else:
                    donor = None
                    lg_smile = 'NO DONOR'
                    #print('NO SANE DONOR FOUND')
            else:
                lg_smile = 'NO PSEUDOPRODUCT'

    
    return lg_smile, donor, acceptor, product

