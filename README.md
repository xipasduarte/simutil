simutil
=======

This repository represents a collection of scripts and small utilitary programs to assist in computational chemistry simulations (using DL_POLY, GROMACS or Towhee). The tools are used to aid in the research carried out at Instituto Superior Técnico (IST), Universidade de Lisboa (UL), Portugal (PT), by the Molecular and Engineering Thermodynamics (MET) group of the Centre for Structural Chemistry (CQE).

_Every tool is designed to work with as little parameters and arguments as possible. Idealy the programs run without the need to provide any, the tool runs taking what it needs from the folder in which it is called._

##Utilities

* ```cfg2pdb``` - Convert a CONFIG file, from DL_POLY, into a PDB (Protein Data Bank) file to be used in RASMOL. (The output can, of course, be used elsewhere.)
* ```cutz``` - Cut a CONFIG file from DL_POLY. Removes all water molecules below a given z value. Used to shorten the layers tickness when studying surface tensions.
* ```dl_energy``` - Retrieve, in a table format, any of the properties present in the STATIS file. (Currently only the first 27 parameters are available.)
* ```gro2cfg``` - Convert a GROMACS ```.gro``` file into a DL_POLY CONFIG file.
* ```hdb``` - Determine the hydrogen bond distribution. No bonds, bonding with the Hydrogen atom, bonding with the Oxygen atom, both.
* ```surfangle``` - Determine the angle between molecules forming a monolayer on top of water, and the surface of the water. Can be easely fitted for other systems. (Uses the DL_POLY HISTORY file.)
* ```surfplane``` - Determine the cluster of molecules that correspond to the top surface in a simulation to obtain surface interactions/properties. (Uses the DL_POLY HISTORY file.)

## License
The simutil package is licensed under The MIT License (MIT).