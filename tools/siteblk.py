# Import modules
import os
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
import argparse
import datetime
import warnings
import sys
warnings.filterwarnings('ignore')
np.set_printoptions(threshold=np.inf)

# Define physical constants
KtoeV = 8.617333262e-5


# Define functions
# For reading jmpdata.csv
def read_jmp(jmpdata_file):
	data = np.loadtxt(jmpdata_file,dtype="str",delimiter=",",comments="#")
	jmpSites = data[:,:2].astype(int) # Initial & final sites for each jump
	jmpVec = data[:,2:5].astype(float) # Jump vector (cartesian) for each jump
	jmpFreq = data[:,5].astype(float) # Jump frequency for each jump [Hz]
	n_site = np.max(jmpSites) # Number of sites (eq. max site ID)
	return jmpSites, jmpVec, jmpFreq, n_site


# For reading sitePE.csv
def read_sitePE(sitePE_file):
	sitePE = np.loadtxt(sitePE_file,dtype=float,delimiter=",",comments="#")
	return sitePE # Indexes start at zero.


# For reading excl.csv
def read_excl(excl_file):
	excl_sids = []
	with open(excl_file) as f:
		for line in f:
			excl_tmp = np.array(line.split(",")).astype(int)-1
			excl_sids.append(excl_tmp)
	return excl_sids # Indexes start at zero.


# For checking the number of carriers
def chk_n_c(n_c):
	if n_c <= 1:
		print("Warning! n = "+str(n_c))
		print("The number of diffusion carriers should be more than 1.")
		print("Abort this program!")
		sys.exit()


# For checking the number of sites in sitePE.csv and excl.csv
def chk_n_site(sitePE,excl_sids):
	if sitePE.shape[0] != n_site:
		print("Warning! n_sites in sitePE.csv is "+str(sitePE.shape[0]))
		print("The number of sites in sitePE.csv is not consistent with jmpdata.csv.")
		print("Abort this program!")
		sys.exit()
	if len(excl_sids) != n_site:
		print("Warning! n_sites in excl.csv is "+str(len(excl_sids)))
		print("The number of sites in excl.dat is not consistent with jmpdata.csv.")
		print("Abort this program!")
		sys.exit()


# For estimating Boltzmann factor at each site
def boltz(sid,T_K,sitePE,excl_sids,occupancies):
	unocc = 1 - occupancies[excl_sids[sid]] # Unoccupied probabilities in exclusive region
	p_unocc = np.prod(unocc) # Unoccupied probability of exclusive region
	# Modified Boltzmann factor at a given site.
	f_boltz = p_unocc * np.exp(-sitePE[sid]/(T_K*KtoeV)) 
	return f_boltz


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
			                help="Atomic jump data under independent-particle approximation.")
	parser.add_argument("--sitePE", type=str, default="sitePE.csv",\
			                help="Site energy data [eV] in the csv format.")
	parser.add_argument("--excl", type=str, default="excl.csv",\
			                help="Exclusive region data. New line for different site.")
	parser.add_argument("--T", type=float, default=1000.0,\
			                help="Temperature [K]")
	parser.add_argument("--n", type=int, default=1,\
			                help="Number of diffusion carriers in the system.")
	parser.add_argument("--prec", type=str, default="exact",\
			                help="Option for estimating site occupancy. exact OR rough")
	args = parser.parse_args()
	jmpdata_file = args.jmp
	sitePE_file = args.sitePE
	excl_file = args.excl
	T_K = args.T
	n_c = args.n
	prec = args.prec

	# Check the number of carriers.
	print("Checking the number of diffusion carriers")
	chk_n_c(n_c)


	# Read jmpdata.csv, sitePE.csv, and excl.csv
	print("Reading jmpdata.csv, sitePE.csv, and excl.csv")
	jmpSites, jmpVec, jmpFreq, n_site = read_jmp(jmpdata_file) # ndarrays except n_site
	sitePE = read_sitePE(sitePE_file) # 1d ndarray. Indexes start at zero.
	excl_sids = read_excl(excl_file) # list of ndarrays. Indexes start at zero.


	# Check the number of sites in sitePE.csv and excl.csv
	print("Checking the number of sites in sitePE.csv and excl.csv")
	chk_n_site(sitePE,excl_sids)


	# Estimate site occupancy
	"""
	Diffusion carriers are introduced in a given system one by one according to Boltzmann distribution. When prec = "rough", n_c carriers are put in the system. When prec = "exact", n_c-1 carriers are introduced in the system except the exclusive region of a focused site occupied by a carrier. The latter manner is time consuming if the number of sites is large (> 10000).In such cases, please try to use prec = "rough".
	"""
	print("Estimating site occupancies")
	# Define variable for site occupancies
	occupancies = np.zeros(n_site)

	# Set number of loops
	n_max = n_c-1
	if prec.lower() == "rough":
		n_max = n_c
	
	print("Start rough estimation")
	for n in range(n_max): # n denotes index of diffusion carrier
		# Estimate Boltzmann factors at all sites for nth diffusion carrier.
		f_boltz = np.zeros(n_site) # Boltzmann factor at each site
		for i in range(n_site): # i denotes site id.
			f_boltz[i] = boltz(i,T_K,sitePE,excl_sids,occupancies) # Boltzmann factor at site i
		
		# Site occupancies of nth diffusion carrier
		occupancies_n = f_boltz/np.sum(f_boltz)
		
		# Total site occupancies of 1st-nth carriers
		occupancies = occupancies + occupancies_n


	# Additional loops for prec == "exact"
	print("Start additional estimation for prec = exact")
	if prec.lower() != "rough":
		occupancies_each = [] # Site occupancies when site i is occupied.
		for i in range(n_site): # site i occupied by a carrier
			occupancies_each.append(occupancies.copy())
			occupancies_each[-1][excl_sids[i]] = 0.0 # Occupancies in excl regions are 0.
			occupancies_each[-1][i] = 1.0 # Occupancy at site i is 1.
			n_rest = np.sum(occupancies[excl_sids[i]]) # A gap of total occupancy
			n_rest_int = int(np.floor(n_rest)) # integer part
			n_rest_float = n_rest-n_rest_int # float part
			# Deficient charge carriers is introduced one by one
			for n in range(n_rest_int+1):
				f_boltz = np.zeros(n_site)
				for j in range(n_site):
					f_boltz[j] = boltz(j,T_K,sitePE,excl_sids,occupancies_each[-1])
				# Site occupancies of nth diffusion carrier
				occupancies_n = f_boltz/np.sum(f_boltz)
				# Total site occupancies of 1st-nth carriers
				if n == n_rest_int:
					occupancies_each[-1] = occupancies_each[-1] + occupancies_n*n_rest_float
				else:
					occupancies_each[-1] = occupancies_each[-1] + occupancies_n

	print("Correcting jump frequencies")
	# Correct jump frequencies by occupancies of exclusive regions
	jmpFreq_corr = np.zeros(jmpFreq.shape[0])
	for i in range(jmpFreq.shape[0]):
		sid_i = jmpSites[i,0]-1 # Index starts at zero.
		sid_f = jmpSites[i,1]-1 # Index starts at zero.
		
		if sid_f in excl_sids[sid_i]: # Final site is in exclusive region of initial site
			jmpFreq_corr[i] = jmpFreq[i]
		else:
			# exclusive sites of final site without the exclusive sites of initial site
			excl_sids_woi = list(set(excl_sids[sid_f].tolist())-set(excl_sids[sid_i].tolist()))
			
			# Unoccupied probability at each site in exclusive region of final site
			if prec.lower() == "rough":
				unocc = 1 - occupancies[excl_sids_woi]
			else:
				unocc = 1 - occupancies_each[sid_i][excl_sids_woi]

			# Unoccupied probability of exclusive region
			p_unocc = np.prod(unocc)
			# Corrected jump frequency
			jmpFreq_corr[i] = p_unocc * jmpFreq[i]

			
	# Save corrected jmpdata in corr_jmpdata.csv.
	print("Output corrected jump frequencies and their correction factors.")
	jmpFreq_corr = jmpFreq_corr.reshape([jmpFreq_corr.shape[0],1])
	jmpdata = np.concatenate([jmpSites,jmpVec,jmpFreq_corr],axis=1)
	np.savetxt("corr_"+jmpdata_file,jmpdata,delimiter=',',\
			       fmt=["%d","%d","%.10f","%.10f","%.10f","%.10e"],\
						 header='#initialSiteID,finalSiteID,s_x[ang.],s_y[ang.],s_z[ang.],freq[Hz]')

	# Save correction factors in corrFactors.dat
	corrFactors = (jmpFreq_corr[:,0]/jmpFreq).reshape([jmpFreq_corr.shape[0],1])
	np.savetxt("corrFactors.dat",corrFactors)
	

	# Finish main program.
	########################################
	print("#"*30)
	date = datetime.datetime.now()
	print(date)
	print("Finish program!")
	print("#"*30)
	########################################





