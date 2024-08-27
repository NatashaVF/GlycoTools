from rdkit import Chem
from rdkit.Chem import AllChem

###Function that determines anomeric configuration both by Mills projection and translates it to the alpha/beta system###
def get_product_conf(donor,lg,product,donor_type, acceptor_type, get_anomer = False,donor_vector = None,):
    if donor_type == 'pyranose':
        donor = Chem.MolFromSmarts(donor)
        donor.UpdatePropertyCache(strict=False)
        product = Chem.MolFromSmarts(product)
        product.UpdatePropertyCache(strict=False)
        reaction_template = AllChem.ReactionFromSmarts('[C:222]1[C:333][C:444][C:555][O:666][C:111]1-{}.['f'{acceptor_type}'':777]>>[C:222]1[C:333][C:444][C:555][O:666][C@H:111]1['f'{acceptor_type}'':777].[C:222]1[C:333][C:444][C:555][O:666][C@@H:111]1['f'{acceptor_type}'':777]'.format(lg[1:]),useSmiles = True)
        pseudo_products = reaction_template.RunReactants((donor,Chem.MolFromSmarts('['f'{acceptor_type}'':1234]')))
        pseudo_product_up = pseudo_products[0][0]
        pseudo_product_down = pseudo_products[0][1]
        if product.HasSubstructMatch(pseudo_product_up,useChirality=True) == True:
            ano_conf = 1
        elif product.HasSubstructMatch(pseudo_product_down,useChirality=True) == True:
            ano_conf = -1
        else:
            ano_conf = 0

        if get_anomer == True:
            if donor_vector[-1] != 0:
                last_stereocenter = donor_vector[-1] 
                if last_stereocenter == ano_conf:
                    anomer = 'beta'
            
                elif last_stereocenter == -ano_conf:
                    anomer = 'alpha'
                
                else:
                    anomer = 'none'

            elif donor_vector[-1] == 0:
                for i in range(1,(len(donor_vector))):
                    last_stereocenter = donor_vector[-i]
                    if last_stereocenter != 0:
                        break
                
                if last_stereocenter == ano_conf:
                    anomer = 'alpha'
                
                elif last_stereocenter == -ano_conf:
                    anomer = 'beta'

                else:
                    anomer = 'could not be determined'

    
    elif donor_type == 'furanose':
        donor = Chem.MolFromSmarts(donor)
        donor.UpdatePropertyCache(strict=False)
        product = Chem.MolFromSmarts(product)
        product.UpdatePropertyCache(strict=False)
        reaction_template = AllChem.ReactionFromSmarts('[C:222]1[C:333][C:444][O:666][C:111]1-{}.['f'{acceptor_type}'':777]>>[C:222]1[C:333][C:444][O:666][C@H:111]1['f'{acceptor_type}'':777].[C:222]1[C:333][C:444][O:666][C@@H:111]1['f'{acceptor_type}'':777]'.format(lg[1:]),useSmiles = True)
        pseudo_products = reaction_template.RunReactants((donor,Chem.MolFromSmarts('['f'{acceptor_type}'':1234]')))
        pseudo_product_up = pseudo_products[0][0]
        pseudo_product_down = pseudo_products[0][1]
        if product.HasSubstructMatch(pseudo_product_up,useChirality=True) == True:
            ano_conf = 1
        elif product.HasSubstructMatch(pseudo_product_down,useChirality=True) == True:
            ano_conf = -1
        else:
            ano_conf = 0

        if get_anomer == True:
            if donor_vector[-2] != 0:
                last_stereocenter = donor_vector[-2] 
                if last_stereocenter == ano_conf:
                    anomer = 'beta'
            
                elif last_stereocenter == -ano_conf:
                    anomer = 'alpha'
                
                else:
                    anomer = 'none'

            elif donor_vector[-2] == 0:
                for i in range(1,(len(donor_vector)-1)):
                    last_stereocenter = donor_vector[-i+1]
                    if last_stereocenter != 0:
                        break
                
                if last_stereocenter == ano_conf:
                    anomer = 'alpha'
                
                elif last_stereocenter == -ano_conf:
                    anomer = 'beta'

                else:
                    anomer = 'could not be determined'
                
    if get_anomer == True:
        return ano_conf, anomer
        
    else:
        return ano_conf