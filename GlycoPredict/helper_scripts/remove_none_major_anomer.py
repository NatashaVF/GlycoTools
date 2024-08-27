import pandas as pd

###Add reactions with anomeric ratios to new dataset, and keep only major in current dataset###
def remove_none_major_anomer(df):   

    df = df.sort_values(by=['rxnid'])
    isomers = []
    i_skip = []
    df['remove'] = ""
    for i in range(len(df)-1):
        if len(i_skip) > 0:
            for k in range(len(i_skip)):
                if i == i_skip[k]:
                    pass
            
            i_skip = []
                    

        else:
            rxn_0 = df['rxnid'][i]
            no_of_duplicates = 0
            for j in range(1,4):
                if rxn_0 == df['rxnid'][i+j]:
                    no_of_duplicates = no_of_duplicates + 1
                else:
                    break
            if no_of_duplicates == 0:
                i_skip = []
                continue
            ##If two duplicates ##
            elif no_of_duplicates == 1:
                if df['yield'][i] > df['yield'][i+1]:
                    df['remove'][i+1] = 'Y'
                elif df['yield'][i] < df['yield'][i+1]:
                    df['remove'][i] = 'Y'
                else: 
                    df['remove'][i] = 'Inconclusive'
                    df['remove'][i+1] = 'Inconclusive'

                ##Check if anomers##
                #if df['product_conf'][i] != df['product_conf'][i]:
                isomers.append(df.loc[i])
                #print(df.iloc[i])
                isomers.append(df.loc[i+1])
                #print(df.iloc[i+1])
            
                i_skip = [i + 1]
                
            ##If 3 duplicates##
            elif no_of_duplicates == 2:
                if df['yield'][i] > df['yield'][i+1] and df['yield'][i] > df['yield'][i+2]:
                    df['remove'][i+1] = 'Y'
                    df['remove'][i+2] = 'Y'
                elif df['yield'][i+1] > df['yield'][i] and df['yield'][i+1] > df['yield'][i+2]:
                    df['remove'][i] = 'Y'
                    df['remove'][i+2] = 'Y'
                elif df['yield'][i+2] > df['yield'][i] and df['yield'][i+2] > df['yield'][i+1]:
                    df['remove'][i] = 'Y'
                    df['remove'][i+1] = 'Y'
                else: 
                    df['remove'][i] = 'Inconclusive'
                    df['remove'][i+1] = 'Inconclusive'
                    df['remove'][i+2] = 'Inconclusive'

                ##Append isomers##
                isomers.append(df.loc[i])
                isomers.append(df.loc[i+1])
                isomers.append(df.loc[i+2])

                i_skip = [i + 1, i + 2]

            ##If 4 duplicates##
            elif no_of_duplicates == 3:
                if df['yield'][i] > df['yield'][i+1] and df['yield'][i] > df['yield'][i+2] and df['yield'][i] > df['yield'][i+3]:
                    df['remove'][i+1] = 'Y'
                    df['remove'][i+2] = 'Y'
                    df['remove'][i+3] = 'Y'
                elif df['yield'][i+1] > df['yield'][i] and df['yield'][i+1] > df['yield'][i+2] and df['yield'][i+1] > df['yield'][i+3]:
                    df['remove'][i] = 'Y'
                    df['remove'][i+2] = 'Y'
                    df['remove'][i+3] = 'Y'
                elif df['yield'][i+2] > df['yield'][i] and df['yield'][i+2] > df['yield'][i+1] and df['yield'][i+2] > df['yield'][i+3]:
                    df['remove'][i] = 'Y'
                    df['remove'][i+1] = 'Y'
                    df['remove'][i+3] = 'Y'
                elif df['yield'][i+3] > df['yield'][i] and df['yield'][i+3] > df['yield'][i+1] and df['yield'][i+3] > df['yield'][i+2]:
                    df['remove'][i] = 'Y'
                    df['remove'][i+1] = 'Y'
                    df['remove'][i+3] = 'Y'
                else: 
                    df['remove'][i] = 'Inconclusive'
                    df['remove'][i+1] = 'Inconclusive'
                    df['remove'][i+2] = 'Inconclusive'
                    df['remove'][i+3] = 'Inconclusive'

                ##Append isomers##
                isomers.append(df.loc[i])
                isomers.append(df.loc[i+1])
                isomers.append(df.loc[i+2])
                isomers.append(df.loc[i+3])

                i_skip = [i + 1, i + 2, i + 3]

    df_iso = pd.DataFrame.from_dict(isomers)
    df_iso.to_csv('data/stereoisomers_cas_all.csv')

    df = df[~df.remove.str.contains('Y')]
    df = df[~df.remove.str.contains('Inconclusive')]
    
    df = df [~df.product_anomer.str.contains('could not be determined')]
    df = df[df['product_conf'] != 0]
    df = df.drop('remove', axis=1)
    df = df.reset_index(drop=True)
    return df


