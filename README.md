# LASTER
Crystal Plasticity Elastic Strain Extractor
## Description
LASTER (crystal pLASTicity Elastic strain extractoR) is a python tool to extract grain family resolved elastic strains from ABAQUS crystal plasticity simulations (that utilise the BRISTOL UMAT). LASTER utilises the materials.dat file produced by Dream3D to identify grains within the [111], [200], [220] and [311] grain families and outputs their strains into respective text files where each column is a grain within that family. An additional text file is also generated that contains simulation time, RVE average strain and RVE average stress in individual columns
## User Instructions
The odB and materials.dat file must be located in the same directory. The script and these aforementioned files do not have to be in the same directory. It is recommnded that this script is placed in the 'C:\Temp' folder for easy access as this is the folder ABAQUS command line defaults.

- Place 'LASTER.py' in the 'C:\Temp' directory and ensure the odB 
