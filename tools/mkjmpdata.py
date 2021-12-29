# Import modules
import os
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
import argparse
import warnings
warnings.filterwarnings('ignore')
np.set_printoptions(threshold=np.inf)


# Define function for reading emig.csv
def read_emig(emig_file):
	data = np.loadtxt(emig_file,dtype="str",delimiter=",",comments="#")
	jmpSites = data[:,:2].astype(int) # Initial & final sites for each jump
	jmpVec = data[:,2:5].astype(float) # Jump vector (cartesian) for each jump
	emig = data[:,5].astype(float) # Migration energy for each jump [eV]
	v0 = data[:,6].astype(float) # Jump frequency for each jump [THz]
	return jmpSites, jmpVec, emig, v0


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
	jmpSites, jmpVec, emig, freq0 = read_emig(emig_file)


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







