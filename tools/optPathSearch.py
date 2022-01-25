# Import modules
import os
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
import argparse
import datetime
from anytree import Node, RenderTree, PreOrderIter
import warnings
warnings.filterwarnings('ignore')
np.set_printoptions(threshold=np.inf)

# Import function from other files
from siteblk import read_sitePE
from mkjmpdata import read_emig 


# Define functions
# Function for making adjacent information of a node.
def mkadjs(node,nodes_smpl,adjIDs,lattVec_inv,sdlPEs_mat):
	adjs = set()
	sid = node.name[0]
	for i in adjIDs[sid]:
		Vec_from_gmin_xyz=node.Vec_from_gmin+jmpVecs_tnsr[sid,i,:]
		Vec_from_gmin_frac = np.dot(Vec_from_gmin_xyz,lattVec_inv)
		cell_ID = np.floor(Vec_from_gmin_frac).astype(int)
		adjs.add((i,cell_ID[0],cell_ID[1],cell_ID[2]))
	adjs = list( adjs - set([node.name for node in PreOrderIter(nodes_smpl[0])]) )
	PEs_sdl = []
	for i in range(len(adjs)):
		PEs_sdl.append(sdlPEs_mat[sid,adjs[i][0]])
	return adjs, PEs_sdl


if __name__ == "__main__":
	# Start main program.
	########################################
	print("#"*30)
	date = datetime.datetime.now()
	print(date)
	print("Start program!")
	print("#"*30, "\n")
	########################################

	# Parse arguments
	parser = argparse.ArgumentParser()
	parser.add_argument("--emig", type=str, default="emig.csv", \
			help="Migration energy Emig for each atomic jump.\
			Initial and final sites, jump vector, Emig, and\
			vibrational prefactor in the csv format.\
			New line for different jumps.\
			Both jumps in the opposite directions have to be specified.")
	parser.add_argument("--sitePE", type=str, default="sitePE.csv",\
			help="Site energy data [eV] in the csv format.")
	parser.add_argument("--a", type=str, default=None,\
			help="Lattice vector a in cartesian coordinate [Ang.]. e.g., 4.23621,0.0,0.0 ")
	parser.add_argument("--b", type=str, default=None,\
			help="Lattice vector b in cartesian coordinate [Ang.]. e.g., 0.0,4.23621,0.0 ")
	parser.add_argument("--c", type=str, default=None,\
			help="Lattice vector c in cartesian coordinate [Ang.]. e.g., 0.0,0.0,4.23621 ")
	parser.add_argument("--n_path", type=int, default=1,\
			help="Number of optimal paths explored in this program. Default: 1")
	parser.add_argument("--gmin", type=int, default=0,\
			help="Site ID of the global minimum point.\
			      If gmin = 0, it is searched in sitePE.csv. Default: 0")

	args = parser.parse_args()
	emig_file = args.emig
	sitePE_file = args.sitePE
	lv_abc = [args.a, args.b, args.c]
	n_path = args.n_path
	sid_gmin = args.gmin

	# Lattice vectors in ndarray
	lattVec = []
	for i in range(3):
		lattVec.append(lv_abc[i].split(","))
	lattVec = np.array(lattVec,dtype=float)
	lattVec_inv = np.linalg.inv(lattVec)


	# Read emig.csv and sitePE.csv
	print("Reading input files.\n")
	jmpSites, jmpVecs, Emigs, freqs0 = read_emig(emig_file)
	n_sites = np.max(jmpSites) # Number of sites
	jmpSites = jmpSites - 1 # jmpSites indexes start at 0.
	sitePEs = read_sitePE(sitePE_file) # 1d ndarray. Indexes start at 0.
	sitePEs = sitePEs - np.min(sitePEs) # PE vs. global minimum site


	# Make saddlePEs matrix and jmpVecs tensor.
	sdlPEs_mat = np.full([n_sites,n_sites],-1.0)
	jmpVecs_tnsr = np.zeros([n_sites,n_sites,3])
	for i in range(Emigs.shape[0]):
		sdlPEs_mat[jmpSites[i][0],jmpSites[i][1]] = Emigs[i]+sitePEs[jmpSites[i][0]]
		jmpVecs_tnsr[jmpSites[i][0],jmpSites[i][1],:] = jmpVecs[i,:]

	# Average the saddle-point PEs of jumps back and forth, just in case.
	sdlPEs_mat = (sdlPEs_mat+sdlPEs_mat.T)/2.0

	# Make adjacent site ID lists
	adjIDs = []
	for i in range(n_sites):
		adjIDs.append(np.where(sdlPEs_mat[i,:]>0.0)[0])
	

	# Start searching optimal paths.
	"""
	Using anyTree, grow the tree with the global min. as the root.
	
	Each node has several attributes as follows:
	name: (siteID,cellID_a,cellID_b,cellID_c)
	parent: Parent [Node class]
	adjs: List of adjacent sites [ (siteID,cellID_a,cellID_b,cellID_c), ... ]
	PEs_sdl: Saddle-point PEs [eV] along the jumps to all adjacent sites.
	Vec_from_gmin: Cartesian coordinate with reference to the global min.

	Add the adjacent site with the lowest saddle-point PEs in all nodes as a new node,
	until two nodes with the same site ID and different cell ID exist in the tree.
	If n_path > 1, adding nodes is iterated,
	until the number of linearly-independent optimal paths reaches n_path.
	"""
	print("Start searching optimal paths!")
	print("-"*30)
	#Get the global minimum point and set the point as the root of tree
	if sid_gmin == 0:
		sid_gmin = np.argmin(sitePEs) # site ID of the global min.
	else:
		sid_gmin = sid_gmin-1
	
	nodes_sid_smpl = np.array([[sid_gmin,0,0,0]]) # Site and cell ids of sampled sites.
	nodes_smpl = [] # List of sampled sites containing Node class.
	# Add the glibal min. as the root of tree
	nodes_smpl.append(Node((sid_gmin,0,0,0),parent=None,adjs=[],PEs_sdl=[],Vec_from_gmin=np.zeros(3)))
	# Add information of adjacent sites to the root as attributes
	nodes_smpl[0].adjs, nodes_smpl[0].PEs_sdl = mkadjs(nodes_smpl[0],nodes_smpl,adjIDs,lattVec_inv,sdlPEs_mat)


	# Iteration of adding nodes
	optPaths = [] # List for storing the informaton of found optimal paths
	optPathVecs = [] # List for storing the vectors of found optimal paths
	n=1
	while len(optPaths) < n_path: # Iterated until # of optPaths becomes n_path
		# Search the lowest sdlPE jump to the adjacent sites from all nodes.
		sdlPE_min = 1.0e+100
		nxt = [] # The information of the next site added to the tree.
		for i in range(len(nodes_smpl)):
			for j in range(len(nodes_smpl[i].adjs)):
				if nodes_smpl[i].PEs_sdl[j] < sdlPE_min:
					sdlPE_min = nodes_smpl[i].PEs_sdl[j]
					nxt = [i,j]
		if len(nxt) == 0:
			break

		# Add the final site of the selected jump as a new node in the tree
		sid_i = nodes_smpl[nxt[0]].name[0] # site ID of parent
		sid_f = nodes_smpl[nxt[0]].adjs[nxt[1]][0] # site ID of a new node
		# Add the new node
		nodes_smpl.append(Node(nodes_smpl[nxt[0]].adjs[nxt[1]],\
				              parent=nodes_smpl[nxt[0]],\
											adjs=[],PEs_sdl=[],\
				              Vec_from_gmin=nodes_smpl[nxt[0]].Vec_from_gmin+jmpVecs_tnsr[sid_i,sid_f,:]))
		
		# Delete the new node from the list of adjacent sites in the parent node.
		del nodes_smpl[nxt[0]].adjs[nxt[1]]
		del nodes_smpl[nxt[0]].PEs_sdl[nxt[1]]

		# Check whether another node with the same site ID already exists as the new node.
		# If it exists, store the information of the optimal path, and the branch is stopped.
		# Otherwise, the information of adjacent sites is added in the new node.
		if nodes_smpl[-1].name[0] in nodes_sid_smpl[:,0]:
			# Find the node with the same site ID, and calculate the translation vector (difference in cell ID).
			row_id = np.where(nodes_sid_smpl[:,0]==nodes_smpl[-1].name[0])[0][0]
			transVec = np.array(nodes_smpl[-1].name[1:]) - nodes_sid_smpl[row_id,1:]
			optPathVecs.append(transVec)
			# If the new optimal path is linearly independent to the others, it is stored, otherwise deleted.
			if np.linalg.matrix_rank(optPathVecs) == len(optPaths):
				del optPathVecs[-1]
			else:
				optPaths.append([row_id,len(nodes_smpl)-1])
				print("Iteration {:d}:\n  Optimal path {:d} is found!".format(n,len(optPaths)))
		else:
			nodes_smpl[-1].adjs,nodes_smpl[-1].PEs_sdl = mkadjs(nodes_smpl[-1],nodes_smpl,adjIDs,lattVec_inv,sdlPEs_mat)
		nodes_sid_smpl = np.concatenate([nodes_sid_smpl,np.array([list(nodes_smpl[-1].name)])])
		
		n=n+1


	# Output the information of the found optimal paths
	print("-"*30+"\n\nFound optimal paths\n"+("-"*30))
	for i in range(len(optPaths)): # Iterations according to the number of the found optimal paths
		print("Optimal path {:1d}:".format(i+1))
		# Output the direction of optimal path.
		print("  Direction: [ {: 1d} {: 1d} {: 1d} ]".\
				  format(optPathVecs[i][0],optPathVecs[i][1],optPathVecs[i][2]))
		path_i = [] # The partial path of initial node to root
		path_f = [] # The partial path of final node to root
		for node in nodes_smpl[optPaths[i][0]].iter_path_reverse():
			path_i.append(node.name)
		for node in nodes_smpl[optPaths[i][1]].iter_path_reverse():
			path_f.append(node.name)
		# Connect the two partial paths
		path_all = np.vstack([np.array(path_i[:-1]), np.array(path_f[::-1])])
		
		# Check The maximum PE along the optimal path and rate-limiting path.
		emig = 0.0
		rate_lim = 0
		for i in range(path_all.shape[0]-1):
			if sdlPEs_mat[path_all[i,0],path_all[i+1,0]] > emig:
				emig = sdlPEs_mat[path_all[i,0],path_all[i+1,0]]
				rate_lim = i
		# Output the potential barrier [eV].
		print("  Emig [eV]: {:.5f}".format(emig))

		# Output the site IDs along the optimal path.
		print("  Path (site ID, [cell ID]):")
		for i in range(path_all.shape[0]):
			print("    ({:6d},  [ {: 1d} {: 1d} {: 1d} ] )".\
					  format(path_all[i,0]+1,path_all[i,1],path_all[i,2],path_all[i,3]))
			if i == rate_lim:
				print("      -- Rate limiting --")
		
		print("-"*30)


	# Finish main program.
	########################################
	print("\n"+"#"*30)
	date = datetime.datetime.now()
	print(date)
	print("Finish program!")
	print("#"*30, "\n")
	########################################



