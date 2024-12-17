import numpy as np
import pandas as pd
from sklearn.impute import KNNImputer
import os

if os.getlogin()=='JGUO258':
    path_diabetes_subphenotypes_predictors_folder = "C:/Users/JGUO258/OneDrive - Emory/Predictors of Subphenotypes"

data_mi = pd.read_csv(path_diabetes_subphenotypes_predictors_folder + '/working/processed/dsppre01_analytic df.csv')

selected_variables = ["study_id", "study", "visit", "year", "exam", "cluster_study_id", "cluster", "newdm",
                      "age", "dmagediag", "sbp", "dbp", "height", "wc", "bmi", "hba1c", "insulinf",
                     "glucosef", "glucose2h", "tgl", "hdlc", "ldlc", "serumcreatinine", "urinecreatinine",
                     "egfr", "apo_a", "apo_b", "uric_acid", "vldlc", "diagDays", "lab_StudyDays", 
                     "anthro_StudyDays", "hc", "triceps", "iliac", "abdominal", "medial", "ast", "alt",
                     "urinealbumin", "uacr", "weight", "homa2b", "homa2ir", "race_eth", "race", "female"]

#drop missing values in the selected variables
data_mi = data_mi[selected_variables]
data_mi.shape

columns_to_impute = ["age", "sbp", "dbp", "height", "wc", "bmi", "hba1c", "insulinf",
                     "glucosef", "glucose2h", "tgl", "hdlc", "ldlc", "serumcreatinine", 
                     "ast", "alt","weight", "homa2b", "homa2ir",
                     "apo_a", "apo_b", "uric_acid", "vldlc", "urinecreatinine","urinealbumin", "uacr",
                     "hc", "triceps", "iliac", "abdominal", "medial","female"]

def impute_data(data, n_neighbors=4):
    imputer = KNNImputer(n_neighbors=n_neighbors)
    data[columns_to_impute] = imputer.fit_transform(data[columns_to_impute])
    return data

# Impute data for the entire dataset
imputed_data = impute_data(data_mi)

# Check if there are any missing values
print(imputed_data.isnull().sum())
print(data_mi.isnull().sum())
# compare the imputed data with the original data
data_mi.describe()
imputed_data.describe()


imputed_data.to_csv(path_diabetes_subphenotypes_predictors_folder + '/working/processed/dsppre02_knn imputation.csv', index=False)

