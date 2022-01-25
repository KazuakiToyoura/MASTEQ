# Import modules
import os
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
import argparse
import warnings
warnings.filterwarnings('ignore')
import pymatgen as mg
from pymatgen.io.cif import CifParser


# define function for reading emig.csv
def read_emig_noneq(emig_noneq_file):
	data = np.loadtxt(emig_noneq_file,dtype="str",delimiter=",",comments="#")
	# initial & final site IDs for each jump
	jmpSiteIDs = data[:,:2]
	# distance for each jump
	jmpDist = data[:,2].astype(float) 
	emig = data[:,3].astype(float) # migration energy for each jump [ev]
	v0 = data[:,4].astype(float) # jump frequency for each jump [thz]
	return jmpSiteIDs, jmpDist, emig, v0


np.set_printoptions(threshold=np.inf)
if __name__ == "__main__":
	# Parse arguments.
	parser = argparse.ArgumentParser()
	parser.add_argument("--cif", type=str, default='site.cif',\
			                help="Structure file name in the cif format\
											containing only sites of diffusion species\
											with symmetry information.")
	parser.add_argument("--emig_noneq", type=str, default="emig_noneq.csv",\
			                help="File name with info. on non-equivalent atomic jumps.\
											Initial and final sites, jump vector, Emig, and\
											vibrational prefactor in the csv format.\
											New line for non-equivalent atomic jumps.\
											Both jumps in the opposite directions have to\
											be specified if they are non-equivalent.\
											Sites have to be specified by ELEMENT SYMBOL!\
											e.g., H, He, Li, Be, B, etc.")
	parser.add_argument("--supercell", type=str, default='1,1,1',\
			                help="Supercell size for expanding unitcell.\
										  Since the inter-site distance is defined\
											as the 1NN distance with periodicity,\
											longer jump vectors than the half of lattice vectors\
											are not accepted in this program.\
											'3,3,3' is enough even for a small cell with few sites.")
	parser.add_argument("--prec", type=float, default=1.0e-4,\
			                help="Precision for fractional coordinate.\
											     Default: 1.0e-4")
	args = parser.parse_args()
	cif_file = args.cif
	emig_noneq_file = args.emig_noneq
	cell_size = np.array((args.supercell).split(","),dtype=int)
	prec = float(args.prec)


	# Read site.cif and emig_noneq.csv
	# site.cif
	sitedata = CifParser(cif_file) # parse cif file
	siteSt = sitedata.get_structures(primitive=False)[0] # structure class
	siteSt.make_supercell(cell_size) # make supercell
	lattVec = siteSt.lattice._matrix # lattice vectors
	list_siteType = list(set(siteSt.species)) # Element class in the list
	

	# emig_noneq.csv
	jmpSiteIDs, jmpDist, emig, v0 = read_emig_noneq(emig_noneq_file)


	# Mapping siteID vs. coordinateID
	sid2cid={} # key: element class, value: ndarray including coord ids.
	for i in list_siteType:
		cids = np.where(np.array(siteSt.species)==mg.Element(i))[0]
		sid2cid[i] = cids

	
	# Make distance matrix
	distMat = siteSt.distance_matrix


	# Output lines for emig.csv one by one
	f = open('emig.csv','w')
	f.write("#initialSiteID,finalSiteID,s_x[Ang.],s_y[Ang.],s_z[Ang.],emig[eV],v0[Hz]\n")
	for n in range(jmpSiteIDs.shape[0]): # Loop for jump types
		cids_i = sid2cid[mg.Element(jmpSiteIDs[n,0])] # coordinate ids of initial site type
		cids_j = sid2cid[mg.Element(jmpSiteIDs[n,1])] # cid candidates of final site type

		for i in cids_i: # Loop for initial sites
			# Search final sites with a specified distance
			cids_f = np.where((distMat[i,cids_j]>jmpDist[n]-prec) &\
					              (distMat[i,cids_j]<jmpDist[n]+prec))[0] + np.min(cids_j)
			
			for j in cids_f: # Loop for final sites
				jmpVec_frac = siteSt.frac_coords[i,:] - siteSt.frac_coords[j,:] # jump vector
				# jmpVec_frac is modified in the range -0.5 <= x <= 0.5 
				while np.any(jmpVec_frac<-0.5): 
					jmpVec_frac = np.where(jmpVec_frac<-0.5,jmpVec_frac+1.0,jmpVec_frac)
				while np.any(jmpVec_frac>0.5):
					jmpVec_frac = np.where(jmpVec_frac>0.5,jmpVec_frac-1.0,jmpVec_frac)

				# Convert fractional to cartesian coordinate
				jmpVec_xyz = np.dot(jmpVec_frac,lattVec)
				
				# Values with negligible magnitude are replaced by 0.0
				jmpVec_xyz = np.where((jmpVec_xyz>-1.0e-10)&(jmpVec_xyz<1.0e-10),\
						     1.0e-10,jmpVec_xyz)
				# Output
				f.write("{:d},{:d},{:.8f},{:.8f},{:.8f},{:.4f},{:.5e}\n"\
						    .format(i+1,j+1,jmpVec_xyz[0],jmpVec_xyz[1],jmpVec_xyz[2],emig[n],v0[n]))
	
	f.close()



