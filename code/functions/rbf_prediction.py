
import numpy as numpy
import scipy as scipy

def rbf_predict(A,XYZ):
  return A(XYZ[:,0],XYZ[:,1],XYZ[:,2])
