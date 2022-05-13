#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  7 21:03:06 2021

@author: maria
"""
import math
import csv 
import os
import sys
from datetime import date
###############################################################################
# USER INPUT # 

working_dir = "/path/to/working/directory/" # Working directory. All python scripts should be here 
data_dir = "/path/to/data/folder/" # Directory for output files

seed = None #Random seed for reproducibility. Any number but 0 

#Select Nup, N or C terminus, and Model
nup = "Nup107"
term = "N"
model = "5A9Q"

# Number of NPCs to be simulated
n_input = 16

############################
# Variability parameters
mag = 1 # Magnitude of irregular variability. 0: No deformation, 15: Strong deformation

# Geometric variability
symmet_input = 8 # Symmetry of simulated NPCs

# Mean taken from input-model if "None". 
rnew = None # Mean radius [nm]
rsigma = None # Standard-deviation radius 

dnew = None # Mean ring distance [nm]
dsigma = None # Standard deviation ring distance

thetanew = None # Mean twist angle [rad]
thetasigma = None # standard deviation twist angle

elliptnew = None #approximate semiminor/semimajor axis ratio. 
elliptsigma = None #Standard deviation of semiminor/semimajor axis ratio

###########################
# Plotting parameters
Overviewplot = {"plot": True, "ellipse": False, "circle": False}
Detailplot2D = {"plot": True, "showforce" : False}
XYvsTime = False
Detailplot3D = False
Animate = False

###########################
# Export parameters 
MakeCSV = False # coordinates of simulated NPCs
featuresCSV = False # Features of whole NPCs
featuresCSV2 = False # Features of NPC components 

###############################################################################
# Change to working directory
working= os.environ.get("WORKING_DIRECTORY", working_dir)
if len(sys.argv) > 1: working = sys.argv[1]
os.chdir( working )

# import other scripts in workind directory 
import DeformNPC
import NPC_plotting
import Analyse_deformed

if ((str(seed) == "0") or (str(seed) == "0.0")): 
    sys.exit("Please pick a random seed other than 0 or 0.0.")


rvar = {"rnew" : rnew, 
"rsigma" : rsigma }

dvar = {"dnew" : dnew,
        "dsigma" : dsigma}


thetavar = {"thetanew" : thetanew,
            "thetasigma" : thetasigma}


elliptvar = {"elliptnew" : elliptnew,
             "elliptsigma" : elliptsigma}

multipleNPCs = DeformNPC.MultipleNPC(nup, term, model, n = n_input, mag = mag, 
                                     nConnect = 2, symmet = symmet_input, 
                                     elliptvar = elliptvar, rvar=rvar, 
                                     dvar = dvar, thetavar = thetavar, seed = seed)

NPCs = multipleNPCs["NPCs"]
zNPCs = multipleNPCs["zNPCs"]
mag = multipleNPCs["mag"]
sigmamult = multipleNPCs["sigmamult"]
z = multipleNPCs["z"]
r = multipleNPCs["temp_r"]
ringAngles = multipleNPCs["ringAngles"]
symmet = multipleNPCs["symmet"]
elliptical = multipleNPCs["elliptical"]
nConnect = multipleNPCs["nConnect"]
thetaold = multipleNPCs["thetaold"]
thetaoffset = multipleNPCs["theta_offset"]

NPCscoords = DeformNPC.MultipleNPCs_coord(NPCs, multipleNPCs["zNPCs"], symmet_input) 


ellipseWholeNPC = Analyse_deformed.Ellipses(NPCscoords)
ellipse_allrings = Analyse_deformed.Ellipses(NPCscoords, ringmode = "z" )
ellipse_CRNR = Analyse_deformed.Ellipses(NPCscoords, ringmode =  "CRNR")

circleWholeNPC = Analyse_deformed.Circles(NPCscoords)
circle_allrings = Analyse_deformed.Circles(NPCscoords, ringmode="z")
circle_CRNR = Analyse_deformed.Circles(NPCscoords, ringmode = "CRNR")    


########Plotting
if Overviewplot["plot"]: NPC_plotting.OverviewPlot(NPCscoords, NPCs, mag, r,  
                                                   ellipse = Overviewplot["ellipse"], 
                                                   circle = Overviewplot["circle"])

####### Detail plots of NPC with index NPCi
NPCi = 0
solution = NPCs[NPCi]
randfs = multipleNPCs["randfs"][NPCi]
fcoords = multipleNPCs["fcoords"][NPCi]


if Detailplot2D["plot"]: NPC_plotting.Plot2D(solution, z, symmet, nConnect, 
                                             forces = fcoords, showforce = Detailplot2D["showforce"], 
                                             legend = True, 
                    trajectory = True, colourcode=True, 
                   springs = True, anchorsprings = True)

if XYvsTime: NPC_plotting.XYoverTime(solution, symmet = 8, legend = True)

if Detailplot3D: NPC_plotting.Plot3D(solution, z, symmet, randfs, fcoords, plotforces = True, viewFrame = -1)

if Animate:
    #%matplotlib qt # run to display animated plot
    NPC_plotting.AnimatedScatter(solution, nConnect, symmet, r, randfs)
######################

def names():
    
    def makeNames(mean, sigma, meanName, sigmaName):
        meanStr = str(meanName) + "_" + str(mean) + "_" if mean else ""
        sigmaStr = str(sigmaName) + "_" + str(sigma) + "_" if sigma else ""  
        return meanStr + sigmaStr
    
    today = date.today()
    rStr = makeNames(rnew, rsigma, "r_mean", "r_sigma")
    dStr = makeNames(dnew, dsigma, "dist_mean", "dist_sigma")
    thetaStr = makeNames(thetanew, thetasigma, "twist_ang_mean", "twist_ang_sigma")
    elongStr = makeNames(elliptnew, elliptsigma, "elong_mean", "elong_sigma")
    symmetStr = "_symmet_" + str(symmet) + "_" if symmet != 8 else ""
    
    
    name = (str(term) + "_" + str(nup) + "_" + str(model) + "_x" + str(n_input) 
        + "_deformmag_" + str(mag) + "_" +  symmetStr + rStr + dStr + thetaStr + 
        elongStr + "seed_" + str(seed) + "_" + today.strftime("%d-%b-%Y"))
    
    return {"rStr": rStr, "dStr" : dStr, "thetaStr" : thetaStr, 
            "elongStr" : elongStr, "symmetStr" : symmetStr}, name


def MakeCSV(): 
    "Export CSV file"
    nameDict, name = names()
    
    with open(data_dir + name + '.csv', 'w', newline='') as csvfile:
            writer1 = csv.writer(csvfile, delimiter=',',
                                    quotechar='|', quoting=csv.QUOTE_MINIMAL)  
            writer1.writerow(["x", "y", "z", "particle"])
            for i in range(len(NPCscoords)):               
                #writer1.writerow([NPCscoords[i,0], NPCscoords[i,1], NPCscoords[i,2], NPCscoords[i,3].astype(int)])
                writer1.writerow(list(NPCscoords[i, :3]) + [NPCscoords[i,3].astype(int)])
     
    # corresponding metadata file             
    f = open(data_dir + name + '_metadata.txt', "x") #TODO: different directory
    
    f.write(str(term) + "_" + str(nup) + "_" + str(model) + "\n")
    
    f.write("random seed: " + str(seed) + "\n")
    f.write("\nmodel layout: ###########\n")
    
    #f.write("kr = 0.7 \n")   
    f.write("nConnect: " + str(nConnect) + "\n")
    
    f.write("\nNPC parameters: ##################")
    
    f.write("\nOrganic deformation parameters: ###########\n")
    
    f.write("deformation magnitude: " + str(mag) + "\n")
    f.write("sigma multiplier " + str(sigmamult) + "\n")
    
    f.write("\nGeometric NPC properties: ###########\n")
    
    f.write("symmetry: " + str(symmet) + "\n")
    f.write("radii per ring: " + str(r) + "\n")
    f.write("z: " + str(z) + "\n") 
    f.write("ring angles (rad): " + str(ringAngles) + "\n")
    f.write("twist angle input (rad): " + str(thetaold) + "\n")
    f.write("Nearest rotational unit, NR - CR (degrees): " + str(math.degrees(thetaoffset)) + "\n")
    f.write("\n")
    
    
    if nameDict["rStr"]: f.write("changed radii: " + nameDict["rStr"].replace("_", " ") + "\n")
    if nameDict["dStr"]: f.write("changed distance: " + nameDict["dStr"].replace("_", " ") + "\n")
    if nameDict["elongStr"]: f.write("minor/major axis input: " + nameDict["elongStr"].replace("_", " ") + "\n")
    if nameDict["thetaStr"]: f.write("twist angles (rad): " + nameDict["thetaStr"].replace("_", " ") + "\n")

    f.close()


def featuresCSV():
    nameDict, name = names()   
    
    c_WholeNPC = circleWholeNPC["circlefeatures"]

    c_name = ["c_x0", "c_xy", "c_radius", "c_squaresums"] 

    
    with open(data_dir + "NPC_features_" + name + '.csv', 'w', newline='') as csvfile: # TODO: change directory 
                writer2 = csv.writer(csvfile, delimiter=',',
                                        quotechar='|', quoting=csv.QUOTE_MINIMAL)  
                writer2.writerow(["NPC"] + c_name + list(ellipseWholeNPC[0].keys()))
                for i in range(n_input):               
                    writer2.writerow([i] + c_WholeNPC[i] + list(ellipseWholeNPC[i].values()))



def featuresCSV_rings():
    nameDict, name = names()   
    
    ellipse_CR = ellipse_CRNR[0]
    ellipse_NR = ellipse_CRNR[1]
    
    _, zname = Analyse_deformed.splitRingsByZ(NPCscoords, 0)
    zname = ["_z" + str(i) for i in zname]
    
    c_name = ["c_x0", "c_xy", "c_radius", "c_squaresums"] 
    
    c_name_rings = []#[None] * len(zname)
    e_name_rings = []
    
    c_name_CR = [i + "_CR" for i in c_name]
    c_name_NR = [i + "_NR" for i in c_name]
    
    e_name_CR = list(ellipse_CR[0].keys())
    e_name_NR = list(ellipse_NR[0].keys())
    
    def flatten(t): return [item for sublist in t for item in sublist]
        
    for i in range(len(zname)):
        c_name_rings.extend( [j + zname[i] for j in c_name])
        e_name_rings.extend(list(ellipse_allrings[0][i].keys()))

    
    with open(data_dir + "NPC_features_rings_" + name + '.csv', 'w', newline='') as csvfile: # TODO: change directory 
                writer2 = csv.writer(csvfile, delimiter=',',
                                        quotechar='|', quoting=csv.QUOTE_MINIMAL)  
                writer2.writerow(["NPC"] + c_name_rings + c_name_CR + c_name_NR 
                                         + e_name_rings + e_name_CR + e_name_NR)
                for i in range(n_input):               
                    writer2.writerow([i] 
                                     + flatten(circle_allrings["circleAll"][i])
                                     + circle_CRNR['circleCR'][i]
                                     + circle_CRNR['circleNR'][i]
                                     
                                     + flatten([list(j.values()) for j in ellipse_allrings[i]])
                                     + list(ellipse_CR[i].values())
                                     + list(ellipse_NR[i].values()))    
    


if MakeCSV: MakeCSV()
if featuresCSV: featuresCSV()
if featuresCSV2: featuresCSV_rings()

#os.system('spd-say "your program has finished"')

