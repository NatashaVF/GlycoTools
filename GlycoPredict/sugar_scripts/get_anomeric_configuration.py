####Function that coverts Mills configuration into alpha or beta####
def get_anomer(product_conf, donor_vector, donor_type):
    if product_conf == 0: #In case the product configuration has been normalized
        ano_conf = -1
    elif product_conf == -1:
        ano_conf = -1
    else:
        ano_conf = 1

    if donor_type == 'pyranose':
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
    
    return anomer, ano_conf
