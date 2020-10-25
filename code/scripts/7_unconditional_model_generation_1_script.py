############################################################################################
#                                                                                          #
#                      UNCONDITIONAL MODEL GENERATION (FORMATION 1)                        #
#                                                                                          #
############################################################################################



#Importing requiring libraries
import os
import gc
import feather
import pandas as pd
import numpy as np

#Setting Working Directory
path="/media/fouedjio/095e0241-4724-4a1f-84f8-9ddda0df2da9/fouedjio/3d_stochastic_implicit_modeling/"
os.chdir(path)
os.getcwd()

#Importing implicit trend models and residual fields
implicit_trend_model=feather.read_dataframe("./outputs/data/ppo/trend_sdf_model_1.feather")
residual_field= feather.read_dataframe("./outputs/data/unconditional_realizations/residual_field_1.feather")
sampling_order=feather.read_dataframe("./outputs/data/unconditional_realizations/sampling_order_1.feather")

#Building prior model realizations
nb_ppo_simu=len(implicit_trend_model.columns)
nb_simu=len(residual_field.columns)
nb_grid=len(residual_field)
unconditional_sdf_model= np.zeros((nb_simu, implicit_trend_model.shape[0]))
for l in range(nb_simu):
 unconditional_sdf_model[l]=implicit_trend_model.iloc[:,int(sampling_order['sample_simu'][l]-1)] +  residual_field.iloc[:,l]


del [[implicit_trend_model,residual_field]]
gc.collect()

feather.write_dataframe(pd.DataFrame(unconditional_sdf_model), './outputs/data/unconditional_realizations/unconditional_sdf_model_1.feather')
feather.write_dataframe(pd.DataFrame(np.delete(unconditional_sdf_model,range(nb_grid),1)), './outputs/data/unconditional_realizations/unconditional_sdf_data_1.feather')
