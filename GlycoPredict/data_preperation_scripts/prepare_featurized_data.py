import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import ast
from GlycoPredict.helper_scripts.scaffold import scaffold_split

from sklearn import svm
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, confusion_matrix, precision_score, recall_score, roc_auc_score, ConfusionMatrixDisplay
from sklearn.model_selection import RandomizedSearchCV, train_test_split
from scipy.stats import randint

###Prepare fingerprint data###
def prepare_data_fp(df):
    X = []
    for ind in df.index:
        X.append(np.concatenate([ast.literal_eval(df['donor_conf'][ind]),[df['temperature_C'][ind]],ast.literal_eval(df['fp donor'][ind]),ast.literal_eval(df['fp acceptor'][ind]),ast.literal_eval(df['solvent_fp'][ind]),ast.literal_eval(df['activator_fp'][ind])]))
    y = np.array(df['product_conf'])

    ###Labels###
    X_labels = ['C-1','C-2','C-3','C-4','C-5', 'Temperature']
    for i in range(2048):
        X_labels.append('donor_fp_'f"{i+1}")
    for i in range(2048):
        X_labels.append('acceptor_fp_'f"{i+1}")
    for i in range(2048):
        X_labels.append('solvent_fp_'f"{i+1}")
    for i in range(2048):
        X_labels.append('activator_fp_'f"{i+1}")

    return X, y, X_labels

###Prepare ohe data###
def prepare_data_ohe(df):
    df = df

    df_ohe = df.loc[:,'C5_smiles':]
    X = []
    for ind in df.index:
        x = ast.literal_eval(df['donor_conf'][ind])
        x.append(df['temperature_C'][ind])
        for i in range(1,df_ohe.shape[1]):
            x.append(df_ohe.iloc[ind,i])
        X.append(x)


    y = np.array(df['product_conf'])

    ###Labels###
    X_labels = ['C-1','C-2','C-3','C-4','C-5', 'Temperature']
    for i in range(1,df_ohe.shape[1]):
        X_labels.append(df_ohe.columns[i])
    
    return X, y, X_labels



