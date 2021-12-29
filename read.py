import numpy as np
np.set_printoptions(threshold=np.inf)

# Define functions for reading input files.
# Read jmpdata.csv
def read_jmp(jmpdata):
	data = np.loadtxt(jmpdata,dtype="str",delimiter=",",comments="#")
	jmpSites = data[:,:2].astype(int) # Initial & final sites for each jump
	jmpVec = data[:,2:5].astype(float) # Jump vector (cartesian) for each jump
	jmpFreq = data[:,5].astype(float) # Jump frequency for each jump [Hz]
	n_site = np.max(jmpSites) # Number of sites (eq. max site ID)
	return jmpSites, jmpVec, jmpFreq, n_site

