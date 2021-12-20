"""
This code is a numerical estimator of atomic diffusivity in crystals using a master equation approach. Only specifying the initial and final sites, jump vector, and jump frequency of every atomic jump in unitcell, the diffusion tensor is calculated under the indipendent particle approximation. See the following references for the detailed theoretical background.

* Kazuaki Toyoura, Takeo Fujii, Kenta Kanamori, Ichiro Takeuchi, Sampling strategy in efficient potential energy surface mapping for predicting atomic diffusivity in crystals by machine learning, Physical Review B 101, 184117 (2020).
* Kazuaki Toyoura, Takeo Fujii, Naoyuki Hatada, Donglin Han, Tetsuya Uda, Carrier–Carrier Interaction in Proton-Conducting Perovskites: Carrier Blocking vs Trap-Site Filling, The Journal of Physical Chemistry C 123, 26823-26830 (2019).

This code was written by Kazuaki Toyoura. Mr. Takeo Fujii has great contribution to the code development with fruitful discussions.

"""

# Import modules
import os
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
import argparse
import warnings
warnings.filterwarnings('ignore')
np.set_printoptions(threshold=np.inf)

# Import defined functions in separated files
import read
import mkJmpMat

if __name__ == "__main__":
	# Parse arguments.
	parser = argparse.ArgumentParser()
	parser.add_argument("--emig", type=str, default="emig.csv",\
			                help="Migration energy Emig for each atomic jump.\
											Initial and final sites, jump vector, Emig, and\
											vibrational prefactor in the csv format.\
											New line for different jumps.\
											Both jumps in the opposite directions have to be specified.")
	parser.add_argument("--lb_T", type=float, default=300.0,\
			                help="Lower bound of temperatures for making jmpdata.csv.")
	parser.add_argument("--ub_T", type=float, default=1000.0,\
			                help="Upper bound of temperatures for making jmpdata.csv.")
	parser.add_argument("--int_T", type=float, default=100.0,\
			                help="Temperature interval for making jmpdata.csv.")
	args = parser.parse_args()
	emig_file = args.emig
	lb_T = args.lb_T
	ub_T = args.ub_T
	int_T = args.int_T


	# Define physical constants
	KtoeV = 8.617333262e-5


	# Read jmpdata.csv.
	jmpSites, jmpVec, emig, freq0 = read.read_emig(emig_file)


	# Make temperature lists
	T = np.arange(lb_T,ub_T+1.0e-10,int_T,dtype=float)


	# Make and save jmpdata.csv for every temperature.
	for t in T:
		# Estimate jump frequencies from emig and v0.
		freq = freq0*np.exp(-emig/(t*KtoeV)) # frequency = v0 EXP(-Emig/kT)
		freq = freq.reshape([freq.shape[0],1])
		jmpdata = np.concatenate([jmpSites,jmpVec,freq],axis=1)
		# Save jmpdata in jmpdata.csv.
		np.savetxt("jmpdata_"+str(int(t))+"K.csv",jmpdata,delimiter=',',fmt=["%d","%d","%.10f","%.10f","%.10f","%.10e"],\
				       header='#initialSiteID,finalSiteID,s_x[ang.],s_y[ang.],s_z[ang.],freq[Hz]')


