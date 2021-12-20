# MASTEQ
Numerical estimator of atomic diffusivity in crystals using a master equation approach

# Demo
Proton diffusivity in cubic-perovskite-structured BaZrO<sub>3</sub>.
![H-BZO](https://user-images.githubusercontent.com/93914342/146717021-419ab676-3871-48a1-8b86-b87be6da5965.png)

# Programming language
Python 3.x (installed through anaconda

Imported modules in this code: numpy, scipy, argparse, copy, datetime, os

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
- --h: Help information. List of options in this code.

- --jmp: File name of atomic jump data in a given system. If it is the default name (jmpdata.csv), the argument is not necessary.

- --nonzero: Non-zero elements in diffusion tensor. Specified non-zero elements are estimated in this program. Specify element indexes separated by comma. Default is 1,2,3,4,5,6. Note that the independent elements are six in diffusion tensor, where
<i>D</i><sub><i>xx</i></sub> = <i>D</i><sub>1</sub>, <i>D</i><sub><i>yy</i></sub> = <i>D</i><sub>2</sub>, <i>D</i><sub><i>zz</i></sub> = <i>D</i><sub>3</sub>, <i>D</i><sub><i>yz</i></sub> = <i>D</i><sub><i>zy</i></sub> = <i>D</i><sub>4</sub>, <i>D</i><sub><i>zx</i></sub> = <i>D</i><sub><i>xz</i></sub> = <i>D</i>.<sub>5</sub>, <i>D</i><sub><i>xy</i></sub> = <i>D</i><sub><i>yz</i></sub> = <i>D</i><sub>6</sub>
For example, “--nonzero 1,2,3” means estimating only the diagonal elements.

- --factor: Scaling factor <i>f</i> for wave vector <b>Q</b>. This factor have to be enough small to correspond to large scale for atomic diffusion. The length of the shortest jump vector lmin is used for the standardized scale of reciprocal space, given by 2π/<i>l<sub>min</sub></i>. The magnitude of wave vector |<b>Q</b>| for estimating diffusion tensor is then (2π/<i>l<sub>min</sub></i>)×<i>f</i>. Default value is 1.0e-3.

# Note
This program assumes a single diffusion network without any other competing network in a given system. For example, if several parallel two-dimensional (2D) networks coexist in a given system, the estimated diffusivity reflects only the diffusivity of the fastest one. In such a case, the jump matrix has several eigenvalues close to zero, and this program prints a warning on the standard output. This warning is often output for anisotropic diffusivity at low temperatures, because a single network at high temperatures can substantially be separated into several networks at low temperatures. Generally, the differences in jump frequency between atomic jumps become larger rapidly with decreasing temperatures, which creates negligible paths for atomic jumps at low temperatures, leading to the separated networks.

The accuracy of estimated diffusivity corresponds to the accuracy of the minimum-magnitude eigenvalue in the jump matrix. The smaller scaling factor <i>f</i> for wave vector |<b>Q</b>| is better in terms of larger scale for atomic diffusion, but worse in terms of the accuracy of eigenvalues in the jump matrix. Please carefully check the temperature dependence of estimated diffusivity, which should be a smooth curve (approximately a straight line) in the Arrhenius plot. In the future, this program has a function for checking the diffusion network (site connectivity) in a given system.

# Example
In example directory, the proton diffusivity in the perfect crystal of BaZrO3 [3] is estimated in the range of 200-1000 K. In the perfect crystal, a proton migrates two types of proton jumps, i.e., rotation around single oxide ions and hopping between adjacent oxide ions. The calculated migration energies of proton rotation and hopping are 0.17 and 0.25 eV, respectively. The input files (jmpdata.csv) at all temperatures can be generated all at once by mkjmpdata.py in the next section. 

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

- --h:		Help information. List of options in this code.
- --emig:	File name of migration energy [eV] and vibrational prefactor [Hz] for every atomic jump. If it is the default name (emig.csv), this argument is not necessary.
- --lb_T:		Lower bound of the temperature range [K]. Default value is 300.
- --ub_T:		Upper bound of the temperature range [K]. Default value is 1000.
- --int_T:		Interval of the temperature step [K]. Default value is 100.

When typing the following command, jmpdata_300K.csv, jmpdata_400K.csv, and jmpdata_500K.csv are generated.
```
python [MASTEQ_DIR]/mkjmpdata.py --emig emig.csv --lb_T 300 --ub_T 500 --int_T 100
```

# Author
* Kazuaki Toyoura, PhD.
  Department of Mater. Sci. & Eng., Kyoto Univ.

# Collaborator
* Takeo Fujii, Mr.
  Department of Mater. Sci. & Eng., Kyoto Univ.

# License
MASTEQ is under [BSD 3-Clause License](https://opensource.org/licenses/BSD-3-Clause)


