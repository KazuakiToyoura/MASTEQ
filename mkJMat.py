import numpy as np
import scipy as sp
from scipy.sparse import csr_matrix, csc_matrix, coo_matrix, lil_matrix
np.set_printoptions(threshold=np.inf)

# Define functions for reading input files.
def mkjm(n_site,jmpSites,jmpFreq,jmpVec,waveNumber_fk): # Read jmpdata.csv
	jmpMat = lil_matrix((n_site,n_site),dtype=complex) # Create jump matrix in lil_matrix.
	for i in range(jmpSites.shape[0]):
		jmpMat[jmpSites[i][0]-1,jmpSites[i][1]-1]=jmpMat[jmpSites[i][0]-1,jmpSites[i][1]-1]+jmpFreq[i]*np.exp(1.0j*np.dot(waveNumber_fk,-jmpVec[i]))
		jmpMat[jmpSites[i][0]-1,jmpSites[i][0]-1]=jmpMat[jmpSites[i][0]-1,jmpSites[i][0]-1]-jmpFreq[i]
	
	# Convert lil_matrix to csr_matrix
	jmpMat = jmpMat.tocsr()

	return jmpMat 


