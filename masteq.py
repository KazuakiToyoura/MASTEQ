"""
ver.1.1.1, Jan. 2022.
MASTEQ is a numerical estimator of atomic diffusivity in crystals using a master equation approach. Only specifying the initial and final sites, jump vector, and jump frequency of every atomic jump in unitcell, the diffusion tensor is calculated under the independent particle approximation. See the following references for the detailed theoretical background. Note that any interaction between diffusion carriers cannot be treated until ver.1.0.2. After ver.1.1.0, site blocking can be treated by correcting jump frequencies based on the site occupancies by siteblk.py in tools directory.
This code was written by Kazuaki Toyoura. Mr. Takeo Fujii has great contribution to the code development with fruitful discussions.


* Kazuaki Toyoura, Takeo Fujii, Kenta Kanamori, Ichiro Takeuchi, Sampling strategy in efficient potential energy surface mapping for predicting atomic diffusivity in crystals by machine learning, Physical Review B 101, 184117 (2020).
* Kazuaki Toyoura, Takeo Fujii, Naoyuki Hatada, Donglin Han, Tetsuya Uda, Carrier–Carrier Interaction in Proton-Conducting Perovskites: Carrier Blocking vs Trap-Site Filling, The Journal of Physical Chemistry C 123, 26823-26830 (2019).
* Takeo Fujii, Kazuaki Toyoura, Tetsuya Uda, Shusuke Kasamatsu, Theoretical study on proton diffusivity in Y-doped BaZrO3 with realistic dopant configurations, Phys. Chem. Chem. Phys. 23, 5908-5918 (2021).

"""

# Import modules
import os
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
import scipy as sp
from scipy.sparse import csr_matrix, lil_matrix
from scipy.sparse.linalg import eigsh
import argparse
import datetime
np.set_printoptions(threshold=np.inf)

# Import defined functions in separated files
import read
import mkJMat

if __name__ == "__main__":
	# Start main program.
	########################################
	print("#"*30)
	date = datetime.datetime.now()
	print(date)
	print("Start program!")
	print("#"*30, "\n")
	########################################

	# Parse arguments.
	parser = argparse.ArgumentParser()
	parser.add_argument("--jmp", type=str, default="jmpdata.csv",\
			                help="Atomic jump data.\
											Initial and final sites, jump vector & frequency in the csv format.\
											New line for different jumps.\
											Both jumps in the opposite directions have to be specified.")
	parser.add_argument("--nonzero", type=str, default="1,2,3,4,5,6",\
			                help="Non-zero elements in diffusion tensor, separated by comma.\
											ex.) 1,2,3: Off-diagonal elements are zero.\
											D_tensor = [ [ D[1], D[6], D[5] ],\
	                                 [ D[6], D[2], D[4] ],\
										               [ D[5], D[4], D[3] ] ]")
	parser.add_argument("--factor", type=float, default=0.001,\
			                help="The scale factor of wave number vectors.\
									    This value have to be enough small\
											to correspond to diffusion scale.")
	args = parser.parse_args()
	jmp_file = args.jmp
	nonzero = (args.nonzero).split(",")
	nonzero = np.array(nonzero,dtype=int)-1
	factor = args.factor


	# Read jmpdata.csv.
	print("Read jmpdata.csv")
	jmpSites, jmpVec, jmpFreq, n_site = read.read_jmp(jmp_file)


	# Determine the magnitude of wave number vector |Q|.
	# Six independent Q vectors are defined.
	waveNumber = [np.array([1.0,0.0,0.0]),\
			          np.array([0.0,1.0,0.0]),\
								np.array([0.0,0.0,1.0]),\
								np.array([0.0,1.0,1.0]),\
								np.array([1.0,0.0,1.0]),\
			          np.array([1.0,1.0,0.0])]
	jmpDist_max=np.max(np.linalg.norm(jmpVec,ord=2,axis=1))
	factor_k = 2.0 * np.pi / jmpDist_max * factor
	waveNumber_fk = np.array(waveNumber) * factor_k


  # Start main loop. Estimate D for each waveNumber.
	print("Start main loop.","\n")
	D = np.zeros(6)
	for n in nonzero:
		print("D("+str(n+1)+")")
		print("-"*40)
		date = datetime.datetime.now()
		print(date)
		# Create jump matrix. jmpMat is csr_matrix in scipy.sparse
		print("Make jump matrix.")
		jmpMat = mkJMat.mkjm(n_site,jmpSites,jmpFreq,jmpVec,waveNumber_fk[n,:])
		
		# Get eigenvalues of jump matrix.
		# Use sp.linalg.eigsh (scipy is much faster for large systems than numpy.)
		print("Solve eigenvalue problem.")
		eigenVal, eigenVec = eigsh(jmpMat.T,which='SM',return_eigenvectors=True)
		eigenVal_sort = np.sort(eigenVal)[::-1] # eigenvalue closest to zero
		
		# Convert eigenvalue to diffusion tensor.
		D_cm2s = -eigenVal_sort[0]/factor_k**2/10.0**16
		if n < 3:
			D[n] = D_cm2s
		elif n == 3:
			D[n] = (D_cm2s-D[1]-D[2])/2
		elif n == 4:
			D[n] = (D_cm2s-D[2]-D[0])/2
		elif n == 5:
			D[n] = (D_cm2s-D[0]-D[1])/2
		
		print("      D [cm2/s]: {: .5E}".format(D[n]))
		
		# Check the accuracy of the estimated diffusivity
		if n <= 2 and D[n] < 0:
			print("Warning: Diffusivity is negative.")
			print("         Eigenvalues are inaccurate.")
		ratio_eigs = np.abs(eigenVal_sort[0]/eigenVal_sort[1])
		if ratio_eigs > 1.0e-2:
			print("Ratio of two smallest eigenvalues: "+str(np.round(ratio_eigs,5)))
			print("Warning: Two (or more) competing diffusion pathways.")
			print("         Diffusivity can be inaccurate.")
		
		print("-"*40,"\n")


	# Output diffusion tensor
	date = datetime.datetime.now()
	print("Estimated diffusvity")
	print("="*55)
	print(date)
	# Output estimated diffusion tensor
	print("Estimated diffusion tensor [cm2/s]:")
	D_tensor = np.array([ [ D[0], D[5], D[4] ],\
			                  [ D[5], D[1], D[3] ],\
				                [ D[4], D[3], D[2] ] ])
	for i in range(3):
		print("{: .5E}  {: .5E}  {: .5E}"\
				  .format(D_tensor[i,0],D_tensor[i,1],D_tensor[i,2]))
	
	# Output diagonalized diffusion tensor
	eigenVal_D, eigenVec_D = np.linalg.eig(np.array(D_tensor))
	# Sort the diffusivity
	eigenVal_D_sort = eigenVal_D[np.argsort(-eigenVal_D)]
	eigenVec_D_sort = eigenVec_D[:,np.argsort(-eigenVal_D)]
	# Sign conversion of eigenvectors
	for i in range(3):
		eigenVec_maxID = np.argmax(np.abs(eigenVec_D_sort[:,i]))
		if eigenVec_D_sort[eigenVec_maxID,i] < 0.0:
			eigenVec_D_sort[:,i] = -eigenVec_D_sort[:,i]
	print("\n"+"Diagonalized diffusion Tensor [cm2/s]:")
	print("(in descending diffusivity sequence)")
	for i in range(3):
		print("-- Direction "+str(i+1)+" --")
		print("D [cm2/s]: {: .5E}".format(eigenVal_D_sort[i]))
		print("Vector: {: .8E} {: .8E} {: .8E}\n"\
				  .format(eigenVec_D_sort[0,i],eigenVec_D_sort[1,i],eigenVec_D_sort[2,i]))
	print("="*55,"\n")

	# Finish main program.
	########################################
	print("#"*30)
	date = datetime.datetime.now()
	print(date)
	print("Finish program!")
	print("#"*30, "\n")
	########################################


