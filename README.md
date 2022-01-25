# MASTEQ
Numerical estimator of atomic diffusivity in crystals using a master equation approach

# Demo
Proton diffusivity in cubic-perovskite-structured BaZrO<sub>3</sub>.
![H-BZO](https://user-images.githubusercontent.com/93914342/146717021-419ab676-3871-48a1-8b86-b87be6da5965.png)

# Programming language
Python 3.x (installed through anaconda, https://www.anaconda.com/products/individual)

Imported modules in this code: numpy, scipy, argparse, copy, datetime, os, sys

<b><i>mkemig.py</i></b> in tools requires <b><i>pymatgen</i></b>, https://pymatgen.org/.

<b><i>optPathSearch.py</i></b> in tools requires <b><i>anytree</i></b>, https://anytree.readthedocs.io/en/latest/index.html.

# Installation
Just download all files in your preferred directory.

# Input file
Only an input file, <b>jmpdata.csv</b>, is required, which includes initial and final site IDs, jump vectors, and jump frequencies for every atomic jump in unitcell. Note that all atomic jumps with initial sites in the unitcell have to be specified including the atomic jumps across the periodic boundary, and that both jumps in the opposite directions have to be specified separately. New line for different atomic jumps. “#” denotes a comment line. Put the following items for every atomic jump separated by comma. (Site IDs have to be sequential numbers starting from 1.)

- Initial site ID [integer]: ID of the initial site for an atomic jump
- Final site ID [integer]: ID of the final site for an atomic jump
- x component of the jump vector [float]: The x component of jump vector in Å for an atomic jump
- y component of the jump vector [float]: The y component of jump vector in Å for an atomic jump
- z component of the jump vector [float]: The z component of jump vector in Å for an atomic jump
- Jump frequency [float]: Jump frequency in Hz for an atomic jump

(jmpdata.csv example)
```
#InitialSiteID,FinalSiteID,jmpVec_x[Ang.],jmpVec_y[Ang.],jmpVec_z[Ang.],frequency[Hz]
1,7,-1.2067005467,1.2067005467,0.0000000000,5.4961104742e+11
1,9,0.0000000000,-0.9114043880,-0.9114043880,1.3907136217e+12
```

# Usage
Type the following command at the directory with jmpdata.csv. The estimated diffusion tensor is printed on the standard output. Please redirect it to a file, e.g., stdout.
```
python [MASTEQ_DIR]/masteq.py --jmp jmpdata.csv --nonzero 1,2,3 --factor 1.0e-3 > stdout
```
# Options
- -h, --help: Help information. List of options in this code.

- --jmp: File name of atomic jump data in a given system. If it is the default name (jmpdata.csv), the argument is not necessary.

- --nonzero: Non-zero elements in diffusion tensor. Specified non-zero elements are estimated in this program. Specify element indexes separated by comma. Default is 1,2,3,4,5,6. Note that the independent elements are six in diffusion tensor, where
<i>D</i><sub><i>xx</i></sub> = <i>D</i><sub>1</sub>, <i>D</i><sub><i>yy</i></sub> = <i>D</i><sub>2</sub>, <i>D</i><sub><i>zz</i></sub> = <i>D</i><sub>3</sub>, <i>D</i><sub><i>yz</i></sub> = <i>D</i><sub><i>zy</i></sub> = <i>D</i><sub>4</sub>, <i>D</i><sub><i>zx</i></sub> = <i>D</i><sub><i>xz</i></sub> = <i>D</i>.<sub>5</sub>, <i>D</i><sub><i>xy</i></sub> = <i>D</i><sub><i>yz</i></sub> = <i>D</i><sub>6</sub>
For example, “--nonzero 1,2,3” means estimating only the diagonal elements.

- --factor: Scaling factor <i>f</i> for wave vector <b>Q</b>. This factor have to be enough small to correspond to large scale for atomic diffusion. The length of the shortest jump vector lmin is used for the standardized scale of reciprocal space, given by 2π/<i>l<sub>min</sub></i>. The magnitude of wave vector |<b>Q</b>| for estimating diffusion tensor is then (2π/<i>l<sub>min</sub></i>)×<i>f</i>. Default value is 1.0e-3.

# Note
This program assumes a single diffusion network without any other competing network in a given system. For example, if several parallel two-dimensional (2D) networks coexist in a given system, the estimated diffusivity reflects only the diffusivity of the fastest one. In such a case, the jump matrix has several eigenvalues close to zero, and this program prints a warning on the standard output. This warning is often output for anisotropic diffusivity at low temperatures, because a single network at high temperatures can substantially be separated into several networks at low temperatures. Generally, the differences in jump frequency between atomic jumps become larger rapidly with decreasing temperatures, which creates negligible paths for atomic jumps at low temperatures, leading to the separated networks.

The accuracy of estimated diffusivity corresponds to the accuracy of the minimum-magnitude eigenvalue in the jump matrix. The smaller scaling factor <i>f</i> for wave vector |<b>Q</b>| is better in terms of larger scale for atomic diffusion, but worse in terms of the accuracy of eigenvalues in the jump matrix. Please carefully check the temperature dependence of estimated diffusivity, which should be a smooth curve (approximately a straight line) in the Arrhenius plot. In the future, this program has a function for checking the diffusion network (site connectivity) in a given system.

# Examples
### H_BaZrO3
The proton diffusivity in the perfect crystal of BaZrO<sub>3</sub> [3] is estimated in the range of 200-1000 K under the independent-particle approximation. In the perfect crystal, a proton migrates two types of proton jumps, i.e., rotation around single oxide ions and hopping between adjacent oxide ions. The calculated migration energies of proton rotation and hopping are 0.17 and 0.25 eV, respectively. The input files (jmpdata.csv) at all temperatures can be generated at once by mkjmpdata.py. An input file for mkjmpdata.py, emig.csv, is also generated by mkemig.py. See the next section for the details.

### H_Y-dopedBaZrO3
The proton diffusivity in 30%Y-doped BaZrO<sub>3</sub> at 600 K is estimated with simple proton-proton interaction by site blocking [1,3]. Specifically, a supercell of 10×10×10 BaZrO<sub>3</sub> unitcells with randomly-distributed 300 Y dopants is employed, where more than one protons are not bonded to the same O ion. The input file, jumpdata.csv, is corrected using the thermal equilibrium site occupancy, which can be performed by siteblk.py. Note that the correlation between subsequent atomic jumps is not treated in this example. See the next section for the details.

# Useful tools
### mkjmpdata.py
This program makes multiple jmpdata.csv all at once in a specified temperature range, using the migration energy &Delta;<i>E</i><sub>mig</sub> and vibrational prefactor &nu;<sub>0</sub> for atomic jump. The jump frequency &nu; is estimated from &nu; = &nu;<sub>0</sub>exp(–&Delta;<i>E</i><sub>mig</sub>/<i>k</i><sub>B</sub><i>T</i>), where <i>k</i><sub>B</sub> and <i>T</i> are the Boltzmann constant and temperature. The input file (default name: emig.csv) have to be prepared in the csv format. The file format is almost same as jmpdata.csv, where &nu; [Hz], is replaced by &Delta;<i>E</i><sub>mig</sub> [eV] and &nu;<sub>0</sub> [Hz] separated by comma. The file example is as follows:

(emig.csv example)
```
#InitialSiteID,FinalSiteID,jmpVec_x[Ang.],jmpVec_y[Ang.],jmpVec_z[Ang.],DEmig[eV],v0[Hz]
1,7,-1.2067005467,1.2067005467,0.0000000000,0.25,1.0e+13
1,9,0.0000000000,-0.9114043880,-0.9114043880,0.17,1.0e+13
```

The temperature range [K] is specified by three parameters, i.e., lb_T (lower bound of temperature), ub_T (upper bound of temperature), and int_T (interval of temperature), which are specified as arguments.

- -h, --help: Help information. List of options in this code.
- --emig:	File name of migration energy [eV] and vibrational prefactor [Hz] for every atomic jump. If it is the default name (emig.csv), this argument is not necessary.
- --lb_T: Lower bound of the temperature range [K]. Default value is 300.
- --ub_T:		Upper bound of the temperature range [K]. Default value is 1000.
- --int_T:		Interval of the temperature step [K]. Default value is 100.

When typing the following command, jmpdata_300K.csv, jmpdata_400K.csv, and jmpdata_500K.csv are generated.
```
python [MASTEQ_DIR]/tools/mkjmpdata.py --emig emig.csv --lb_T 300 --ub_T 500 --int_T 100
```

### mkemig.py
This program makes emig.csv from only the information on crystallographically non-equivalent atomic jumps. Specifically, two files are required, site.cif with site information of diffusion carriers and emig_noneq.csv with jump information including the initial and final site types, jump distance [ang.], migration energy [eV], and vibrational prefactor [Hz]. All equivalent sites in the unitcell are generated using symmetry operations in site.cif, and the equivalent atomic jumps are then searched by the jump distance and the initial and final site types specified in emig_noneq.csv.

For reading site.cif by <i>pymatgen</i>, non-equivalent site types of diffusion carriers have to be distinguished by element symbols (H, He, Li, Be, …), whatever the diffusion carriers are. The example of site information line is shown below.

(Site information example in site.cif)
```
loop_
   _atom_site_label
   _atom_site_occupancy
   _atom_site_fract_x
   _atom_site_fract_y
   _atom_site_fract_z
   _atom_site_type_symbol
   H      1.0    0.000000  0.000000  0.000000    H
   He     1.0    0.500000  0.500000  0.500000    He
   Li     1.0    0.000000  0.500000  0.500000    Li
```

Concerning the file format of emig_noneq.csv, non-equivalent atomic jumps are specified by the initial and final site types (element symbols in site.cif), the jump distance [ang.], &Delta;<i>E</i><sub>mig</sub> [eV], and &nu;<sub>0</sub> [Hz], separated by comma. New line for different atomic jumps. Both jumps in the opposite directions have to be specified separately, if they are non-equivalent. “#” denotes a comment line. The file example is as follows:

(emig_noneq.csv example)
```
#InitialSiteType,FinalSiteType,jumpDistance[Ang.],DEmig[eV],v0[Hz]
H,Li,1.41421356,0.30,1.0e+13
Li,H,1.41421356,0.20,1.0e+13
He,Li,1.00000000,0.35,1.0e+13
Li,He,1.00000000,0.30,1.0e+13
```

Type the following command at the directory with site.cif and emig_noneq.csv, emig.csv is generated.
```
python [MASTEQ_DIR]/tools/mkemig.py --cif site.cif –emig_noneq emig_noneq.csv --supercell 1,1,1 --prec 1.0e-4
```
The options of this program are as follows:
- -h, --help: Help information. List of options in this code.
- --cif: File name of the cif file with site information. If it is the default name (site.cif), this argument is not necessary.
- -- emig_noneq: File name with information on non-equivalent atomic jumps. Initial and final site type, jump vector, Emig, and vibrational prefactor in the csv format. New line for non-equivalent atomic jumps. Both jumps in the opposite directions have to be specified if they are non-equivalent. Site types have to be specified by ELEMENT SYMBOL! e.g., H, He, Li, Be, B, etc. If it is the default name (emig_noneq.cif), this argument is not necessary.
- -- supercell:	Supercell size for expanding unitcell. Since the inter-site distance is defined as the 1NN distance with periodicity, longer jump vectors than the half of lattice vectors are not accepted in this program. '3,3,3' is enough even for a small cell with few sites. Default value: 1,1,1
- --prec: Precision used in symmetry operation and atomic jump search. Default value: 1.0e-4

### siteblk.py
This program corrects jmpdata.csv taking a simple carrier-carrier interaction into account, i.e., <b><i>site blocking</i></b>. <b><i>Site blocking</i></b> means that any site in a given system cannot be occupied by more than one diffusion carriers at the same time. The thermal equilibrium site occupancy can easily be estimated from Fermi distribution. In this program, the site blocking is extended by defining exclusive region around every site, which cannot be occupied by the other carriers as long as the focused site is occupied by a diffusion carrier. Thermal equilibrium site occupancies are used as site blocking probability, where the correlation between subsequent atomic jumps is not considered.

The three input files are required in this program as follows:
- jmpdata.csv: jmpdata.csv under independent-particle approximation without any carrier-carrier interaction.
- sitePE.csv: Potential energy (PE) at each site in a given system. All the PE values [eV] are written in a single row separated by comma, OR are written one by one in a single column.

  (sitePE.csv example)
  ```
  0.12,0.34,0.03,0.42, ...
  ```
- excl.csv: Exclusive region (exclusive sites) around every site, which is specified by site IDs (Site IDs have to be sequential numbers starting from 1.). Each line specifies the exclusive sites around the focused site with the same ID as the line. The exclusive sites have to be separated by comma in each row, including the focused site ID.

  (excl.csv example)
  ```
  1,2,3
  2,3,4
  2,3,5,6
  ...
  ```

Type the following command at the directory with jmpdata.csv, sitePE.csv, and excl.csv, and the corrected jmpdata.csv, corr_jmpdata.csv, is generated. The other output file, corrFactors.dat is also generated, in which the correction factor for each atomic jump is listed.
```
python [MASTEQ_DIR]/tools/siteblk.py --jmp jmpdata.csv --sitePE sitePE.csv --excl excl.csv --T 500 --n 300 --prec exact
```
The options of this program are as follows:
- -h, --help: Help information. List of options in this code.
- --jmp: File name of atomic jump data in a given system under the independent-particle approximation. Default name (jmpdata.csv) can be skipped.
- --sitePE: File name with the information of site PEs. Default name (sitePE.csv) can be skipped.
- --excl: File name with the information of exclusive regions. Default name (excl.csv) can be skipped.
- --T: Temperature [K]. Default value is 1000.
- --n: The number of diffusion carriers in the system [integer]. Default value is 1, meaning no correction.
- --prec: Precision for estimation of site occupancies [string]. 'exact' or 'rough'. Default value is 'exact'. If the number of sites in a given system is too large (> 10000), the occupancy estimation is time consuming. In such case, try to use prec = 'rough'.

### optPathSearch.py
This program searches the optimal path in a given crystal, which is defined as the lowest-energy path between two global minimum points separated by a lattice translation vector. In other words, the optimal path corresponds to the long-range migration pathway in the crystal. In general, there are three linearly-independent optimal paths. This program is based on dynamic programing. See Ref. [4] for the details of the algorithm. This program requires <b><i>emig.csv</i></b> and <b><i>sitePE.csv</i></b>. In addition, the three lattice vectors (<b>a</b>, <b>b</b>, and <b>c</b>) are also required, which should be specified in Cartesian coordinate [Å]. Type the following command at the directory with emig.csv and sitePE.csv, and the information of optimal paths is printed on the standard output. The direction of optimal path is shown by integer ratios of the lattice vectors. Note that the direction and trajectory of the optimal path could be only an example, particularly in a highly-symmetric crystal. For example, the diffusivity is isotropic in a cubic system, where any long-range path has the same potential barrier.
```
python [MASTEQ_DIR]/tools/optPathSearch.py --emig emig.csv --sitePE sitePE.csv --a 4.23621,0.0,0.0 --b 0.0,4.23621,0.0 --c 0.0,0.0,4.23621 --n_path 3 --gmin 0 > stdout_ops
```
The options of this program are as follows:
- -h, --help: Help information. List of options in this code.
- --emig: File name of migration energy [eV] and vibrational prefactor [Hz] for every atomic jump. If it is the default name (emig.csv), this argument is not necessary.
- --sitePE: File name with the information of site PEs. Default name (sitePE.csv) can be skipped.
- --a: Lattice vector <b>a</b> in Cartesian coordinate [Å]. e.g.)  4.23621,0.00000,0.00000
- --b: Lattice vector <b>b</b> in Cartesian coordinate [Å]. e.g.)  0.00000, 4.23621,0.00000
- --c: Lattice vector <b>c</b> in Cartesian coordinate [Å]. e.g.)  0.00000,0.00000,4.23621
- --n_path: Number of optimal paths explored in this program. Default value is 1.
- --gmin: Site ID of the global minimum point. If gmin = 0, it is searched in sitePE.csv. Default value is 0.

# Author
* Kazuaki Toyoura, PhD.
  Department of Mater. Sci. & Eng., Kyoto Univ.

# Collaborator
* Takeo Fujii, Mr.
  Department of Mater. Sci. & Eng., Kyoto Univ.

# References
[1] K. Toyoura, T. Fujii, N. Hatada, D. Han, T. Uda, Carrier–Carrier Interaction in Proton-Conducting Perovskites: Carrier Blocking vs Trap-Site Filling, The Journal of Physical Chemistry C 123, 26823-26830 (2019).

[2] K. Toyoura, T. Fujii, K. Kanamori, I. Takeuchi, Sampling strategy in efficient potential energy surface mapping for predicting atomic diffusivity in crystals by machine learning, Phys. Rev. B 101, 184117/1-11 (2020).

[3] T. Fujii, K. Toyoura, T. Uda, S. Kasamatsu, Theoretical study on proton diffusivity in Y-doped BaZrO<sub>3</sub> with realistic dopant configurations, Phys. Chem. Chem. Phys. 23, 5908-5918 (2021).

[4] K. Kanamori, K. Toyoura, J. Honda, K. Hattori, A. Seko, M. Karasuyama, K. Shitara, M. Shiga, A. Kuwabara, I. Takeuchi, Exploring a potential energy surface by machine learning for characterizing atomic transport, Phys. Rev. B 97, 125124/1-6 (2018).

# License
MASTEQ is under [BSD 3-Clause License](https://opensource.org/licenses/BSD-3-Clause)

