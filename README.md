# LASTER
Crystal Plasticity Elastic Strain Extractor
## Description
LASTER (crystal pLASTicity Elastic strain extractoR) is a python tool to extract grain family resolved elastic strains from ABAQUS crystal plasticity simulations (that utilise the BRISTOL UMAT). LASTER utilises the materials.dat file produced by Dream3D to identify grains within the [111], [200], [220] and [311] grain families and outputs their strains into respective text files where each column is a grain within that family. An additional text file is also generated that contains simulation time, RVE average strain and RVE average stress in individual columns
## User Instructions
The odB and materials.dat file must be located in the same directory. The script and these aforementioned files do not have to be in the same directory. It is recommnded that this script is placed in the 'C:\Temp' folder for easy access as this is the folder ABAQUS command line defaults.

- Place 'LASTER.py' in the 'C:\Temp' directory.
- Ensure the odB and materials.dat files are in the same directory. They do not have to be in 'C:\Temp'.
- Open ABAQUS command. This should default to the 'C:\Temp' directory. If it does not, navigate to where you have placed LASTER.py.
- Enter command 'abaqus python LASTER.py'. The tool should initialise and begin asking for user inputs.

### User Inputs
- File Directory: Paste in the directory where your odB and materials file are located.
- odB File Name: As described, case sensitive. Do not include the '.odb' file extension.
- Integration Angle: Degrees. Grains within this angle will be classed as a certain grain family. For example, imagine you are trying to identify grains in your CP simulation that are in the [111] family. There will be little to none grains that line up exactly with the [111] definition hence a tolerance is defined. Selecting this angle as say 7.5 degrees will identify grains within a 15 degree solid angle to the [111] orientation as being [111] grains
- Loading Direction: Acceptable inputs are 'X', 'Y' or 'Z' - case sensitive. This defines the direction in which you have loaded your RVE.
- Strain Direction: Acceptable inputs are 'X', 'Y' or 'Z' - case sensitive. This defines which strain you want to extract. For example, if you enter 'Y' for the loading direction and 'Y' for the strain direction then the strain outputted will be the strain in the loading direction. If you enter 'Y' for the loading direction and 'X' for the strain direction you will be extracting a strain perpendicular to the loading direction - a transverse strain.
- Material File Name: As described, case sensitive. Do not include the '.dat' file extension.
- Step Name: Case sensitive. The step name in the ABAQUS CP simulation that you want to extract data from.

## Outputs
Providing the code has run successfully 5 text files should be outputted. The grain family resolved strains:
- G_111.txt
- G_200.txt
- G_220.txt
- G_311.txt

and the RVE averaged values:
- RVE_Avg_Values.txt
