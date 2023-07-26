import numpy as np
import math
from odbAccess import *
from numpy import linalg as LA
from itertools import islice
import os

print("    __    ___   _____________________ ")
print("   / /   /   | / ___/_  __/ ____/ __ \ ")
print("  / /   / /| | \__ \ / / / __/ / /_/ /")
print(" / /___/ ___ |___/ // / / /___/ _, _/") 
print("/_____/_/  |_/____//_/ /_____/_/ |_|")
print("")
print("LASTER  is a  tool   to extract  grain")
print("family   resolved    elastic   lattice")
print("strains from ABAQUS crystal plasticity")
print("simulations. This tool requires 7 user")
print("inputs,   specific details/formats can")
print("be  found in the    README  or  on the")
print("GitHub page.  This  tool  will  output")
print("lattice strains from the [111], [200],")
print("[220]  and  [311] into respective text")
print("files  where each  column  is a  grain")
print("within that family. An additional text")
print("file  will  also  be   generated  that")
print("contains simulation  time, RVE average")
print("strain   and  RVE   average stress  in")
print("respective columns.")
print("")
print("--------------------------------------")
print("")
cd = raw_input("Enter Directory:    ")
os.chdir(cd)
odB_Name = raw_input("odB File Name:      ")
odB_Name = odB_Name + ".odb"
Int_Ang = input("Integration Angle:  ")
LDir_Flag = raw_input("Loading Direction:  ")
SDir_Flag = raw_input("Strain Direction:   ")
Mat_Name = raw_input("Material File Name: ")
Mat_Name = Mat_Name + ".dat"
Step_Name = raw_input("Step Name:          ")
print("")
print("--------------------------------------")
print("")

Load_Dirs = {"X": [1,0,0], "Y": [0,1,0], "Z": [0,0,1]}
Strain_Dirs = {"X": 'UVARM1', "Y": 'UVARM2', "Z": 'UVARM3'}
Grain_Fams = {1: [(1,0,0), (-1,0,0), (0,1,0), (0,-1,0), (0,0,1), (0,0,-1)], 2: [(1,1,0), (1,0,1), (0,1,1), (1,-1,0), (1,0,-1), (0,-1,1)], 3: [(1,1,1), (-1,-1,-1), (-1,1,1), (1,-1,-1), (-1,1,-1), (1,-1,1), (-1,-1,1), (1,1,-1)], 4: [(3,1,1), (-1,-3,-1), (1,3,1), (-3,-1,-1), (1,1,3), (-1,-1,-3), (1,3,-1), (1,-3,1), (1,1,-3), (-1,1,3), (-1,1,-3), (-1,-3,1), (-1,3,-1), (-1,3,1), (-1,-1,-3), (1,-1,3), (3,1,-1), (-3,1,1), (1,-3,-1), (3,-1,1), (-3,1,-1), (1,-1,-3), (3,-1,-1), (-3,-1,1)]}

L_Dir = Load_Dirs.get(LDir_Flag, None)
S_Dir = Strain_Dirs.get(SDir_Flag, None)    

print("Reading Materials Data...")
materials_data = np.loadtxt(Mat_Name, dtype=float)
materials_data_raw = materials_data
grains_total = int(materials_data[-1][4])
columns_to_check = [1, 2, 3, 4]

materials_data = materials_data[:, columns_to_check]
materials_data = [list(x) for x in set(tuple(row) for row in materials_data)]

Grain_List = np.zeros((4,10000))
print("Identifying Grains in Four Families...")
for z in range(1,5):
	Grain_Fam_N = Grain_Fams.get(z, None)
	L_Dir = np.array(L_Dir)
	LDirXtal = np.zeros((grains_total, 3))
	for i in range(0, grains_total):
		phi1 = math.radians(materials_data[i][0])
		phi2 = math.radians(materials_data[i][2])
		phi  = math.radians(materials_data[i][1])
		Rz2 = [(math.cos(phi2),math.sin(phi2),0), (math.sin(phi2)*-1,math.cos(phi2),0), (0,0,1)]
		Rx = [(1,0,0), (0,math.cos(phi),math.sin(phi)), (0,math.sin(phi)*-1,math.cos(phi))]
		Rz1 = [(math.cos(phi1),math.sin(phi1),0), (math.sin(phi1)*-1,math.cos(phi1),0), (0,0,1)]
		Rz2 = np.array(Rz2)
		Rx = np.array(Rx)
		Rz1 = np.array(Rz1)
		g_xtal = Rz2.dot(Rx).dot(Rz1)
		LDirXtal[i] = g_xtal.dot(L_Dir.transpose())
	Ang = np.zeros((1,len(Grain_Fam_N)))
	Grain_Fam_N = np.array(Grain_Fam_N)
	AngF = np.zeros((grains_total, 1))
	for j in range (0, grains_total):
		for k in range(0, len(Grain_Fam_N)):
			Ang[0][k]= math.degrees(math.acos((Grain_Fam_N[k,:].dot(LDirXtal[j,:])/(LA.norm(LDirXtal[j,:])*LA.norm(Grain_Fam_N[k,:])))))
		AngF[j][0]=np.amin(Ang)
	AngF = [item for sublist in AngF for item in sublist]
	IPFPoints = []
	for i in range(0, grains_total):
		if AngF[i] < Int_Ang:
			IPFPoint = (materials_data[i][0], materials_data[i][1], materials_data[i][2], materials_data[i][3])
			IPFPoints.append(IPFPoint)
	for i in range(0,len(IPFPoints)):
		Grain_List[z-1][i] = IPFPoints[i][3]

Grain_List = [tuple(filter(lambda x: x != 0, arr)) for arr in Grain_List]

Grain_Els = [[] for _ in range(len(set(materials_data_raw[:, 4])))]
for x, y in zip(materials_data_raw[:, 0], materials_data_raw[:, 4]):
    Grain_Els[int(y) - 1].append(x)

Grain_Res_Els = [tuple([tuple(Grain_Els[int(x) - 1]) for x in sublist]) for sublist in Grain_List]

odb = openOdb(odB_Name)
FrameList = odb.steps[Step_Name].frames
NumFrames = len(FrameList)
E22_Data = []
print("Calculating Resolved Elastic Strains...")
for j in range(0, NumFrames):
    lastFrame = odb.steps[Step_Name].frames[j]
    E22 = lastFrame.fieldOutputs[S_Dir]
    E22 = np.copy(E22.bulkDataBlocks[0].data)
    k=0
    noEl = len(E22)/8
    #initialize StressEl:
    MisO = np.zeros([noEl,9])
    #Compress to just elements.
    for i in range(0,len(E22)/8):
        MisO[i,0] = np.mean(E22[k:k+7])
        k=k+8
    #Create list of only front elements.
    for z in range(0, len(Grain_Res_Els)):
        for l in range(0, len(Grain_Res_Els[z])):
            i=0
            MisOFront = np.zeros([len(Grain_Res_Els[z][l]),10])
            for el in Grain_Res_Els[z][l]:
                MisOFront[i,0] = el
                MisOFront[i,1:10] =  MisO[el-1,:]
                i=i+1
            StressFront = np.array(MisOFront[:,1])
            E22_Average = np.mean(StressFront)
            E22_Data.append(E22_Average)

E22_Data = np.array(E22_Data)

def slice_per(source, step):
    return [source[i::step] for i in range(step)]

Resolved_Grain = len(Grain_Res_Els[0]) + len(Grain_Res_Els[1]) + len(Grain_Res_Els[2]) + len(Grain_Res_Els[3])
E22_Data = slice_per(E22_Data, Resolved_Grain)
Resolved_Grain = len(Grain_Res_Els[0]), len(Grain_Res_Els[1]), len(Grain_Res_Els[2]), len(Grain_Res_Els[3])
E22_Data = iter(E22_Data)
G_Res = [list(islice(E22_Data, elem)) for elem in Resolved_Grain]

np.savetxt("G_200.txt", np.transpose(G_Res[0]))
np.savetxt("G_220.txt", np.transpose(G_Res[1]))
np.savetxt("G_111.txt", np.transpose(G_Res[2]))
np.savetxt("G_311.txt", np.transpose(G_Res[3]))

Strain_Dirs = {"X": 0, "Y": 1, "Z": 2}
S_Dir = Strain_Dirs.get(SDir_Flag, None)

print("Calculating RVE Stress...")
S22_Data = []
for j in range(0, NumFrames):
    lastFrame = odb.steps[Step_Name].frames[j]
    StressNode = lastFrame.fieldOutputs['S']
    StressNode = np.copy(StressNode.bulkDataBlocks[0].data)
    k=0
    noEl = len(StressNode)/8
    #initialize StressEl:
    StressEl = np.zeros((noEl))
    for i in range(0,len(StressNode)/8):
        StressEl[i] = np.mean(StressNode[k:k+7,S_Dir])
        k=k+8
    S22_Temp = np.mean(StressEl)
    S22_Data.append(S22_Temp)

S22_Data = np.array(S22_Data)
S22_Data = np.transpose(S22_Data)

print("Calculating RVE Strain...")
LE22_Data = []
FrameValues = []
for j in range(0, NumFrames):
    lastFrame = odb.steps[Step_Name].frames[j]
    FrameValue = lastFrame.frameValue
    StressNode = lastFrame.fieldOutputs['LE']
    StressNode = np.copy(StressNode.bulkDataBlocks[0].data)
    k=0
    noEl = len(StressNode)/8
    #initialize StressEl:
    StressEl = np.zeros((noEl))
    for i in range(0,len(StressNode)/8):
        StressEl[i] = np.mean(StressNode[k:k+7,S_Dir])
        k=k+8
    LE22_Temp = np.mean(StressEl)
    LE22_Data.append(LE22_Temp)
    FrameValues.append(FrameValue)

LE22_Data = np.array(LE22_Data)
LE22_Data = np.transpose(LE22_Data)

RVE_Avg_Values = np.column_stack((FrameValues, LE22_Data, S22_Data))
np.savetxt('RVE_Avg_Values.txt', RVE_Avg_Values)

print("")
print("Complete!")