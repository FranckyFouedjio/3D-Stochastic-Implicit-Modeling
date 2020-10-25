############################################################################################
#                                                                                          #
#                       FAST EIGEN VALUE DECOMPOSITION FUNCTION                            #
#                                                                                          #
############################################################################################


import pandas as pd
import numpy as np
import gc
import feather


def full_evd_fast(path):
  
   X=np.asarray(feather.read_dataframe(path))
   center_X=X.mean(axis=0)
   scale_X=X.std(axis=0)
   X = (X -  center_X)/scale_X
   feather.write_dataframe(pd.DataFrame(X), './outputs/data/pca_unconditional_realizations/X.feather')
   Y=X.dot(X.transpose())
   del [[X]]
   gc.collect()

   eig_val, eig_vec = np.linalg.eig(Y)
   del [[Y]]
   gc.collect()
   
   feather.write_dataframe(pd.DataFrame(eig_vec), './outputs/data/pca_unconditional_realizations/vectors.feather')


   new_eig_vecs = []
   new_pc_scores = []

   X=np.asarray(feather.read_dataframe('./outputs/data/pca_unconditional_realizations/X.feather')).T
 

   for i in range(X.shape[1]):
  

      new_vec = X.dot(eig_vec[:,i:i+1])

      new_eig_vecs.append(new_vec[:,0]/np.linalg.norm(new_vec))
      
      del [[new_vec]]
      gc.collect()
   
   del [[X]]
   gc.collect()
   
   new_eig_vecs = np.asarray(new_eig_vecs).T
   gc.collect()
   
   X=np.asarray(feather.read_dataframe('./outputs/data/pca_unconditional_realizations/X.feather'))
   new_pc_scores= X.dot(new_eig_vecs)
   
   del [[X]]
   gc.collect()

   return new_eig_vecs, eig_val, center_X, scale_X, new_pc_scores


