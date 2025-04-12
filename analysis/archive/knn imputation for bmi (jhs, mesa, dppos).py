# load libraries
import numpy as np
import pandas as pd
from sklearn.impute import KNNImputer
import os
# load the data

# change the path to the location of the dataset. Zhongyu uses a MacBook so it cannot read (but can write into) the path from the OneDrive folder.
if os.getlogin()=='JGUO258':
  path_diabetes_subphenotypes_predictors_folder = "C:/Users/JGUO258/OneDrive - Emory/Papers/Predictors of Subphenotypes"

data_mi = pd.read_csv(path_diabetes_subphenotypes_predictors_folder + '/working/processed/8 cohorts/dsppre01_analytic df for bmi imputation.csv')

#select variables 
selected_variables = ['study_id','bmi', 'study','wave','age']

#drop missing values in the selected variables
data_mi = data_mi[selected_variables]

data_mi.head()



############### Do KNN Imputation with K = 5 by Study Sites #####################
columns_to_impute = ['bmi']

# Function to impute data for one study site
def impute_study_data(data, columns_to_impute, n_neighbors=5):
    imputer = KNNImputer(n_neighbors=n_neighbors)
    data[columns_to_impute] = imputer.fit_transform(data[columns_to_impute])
    return data
# Impute data for each study site
study_sites = data_mi['study'].unique()
imputed_datasets = []

for site in study_sites:
    site_data = data_mi[data_mi['study'] == site].copy()
    imputed_data = impute_study_data(site_data, columns_to_impute)
    imputed_datasets.append(imputed_data)

# Merge all imputed datasets back to one
imputed_data_merged = pd.concat(imputed_datasets)

# Check if there are any missing values
print(imputed_data_merged.isnull().sum())
print(data_mi.isnull().sum())
# compare the imputed data with the original data
data_mi.describe()
imputed_data_merged.describe()

# save the imputed data
imputed_data_merged.to_csv(path_diabetes_subphenotypes_predictors_folder + '/working/processed/8 cohorts/dsppre01_analytic df with imputed bmi.csv', index=False)

