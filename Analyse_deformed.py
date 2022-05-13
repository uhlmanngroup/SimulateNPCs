#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 17 18:13:18 2021

@author: maria
"""
import math 
import numpy as np 
import circle_fit as cf

def NPCcoords_split(NPCscoords, z, symmet):
    nRings = len(z)
    n = int(len(NPCscoords)/(symmet*nRings))
    return nRings, n, NPCscoords.reshape((n, nRings, symmet, 4))


def opposite_corner_radius(NPCscoords, z, symmet):   
    nRings, n, NPCscoordssplit = NPCcoords_split(NPCscoords, z, symmet)

    r_op = np.zeros((NPCscoordssplit[:,:,:,0]).shape)

    isodd = symmet%2
    xyopp = lambda npc, ring, opp : np.mean(
        NPCscoordssplit[npc, ring, (opp, (opp+isodd)%symmet), 0:2], axis = 0)
        
    for npc in range(n):
        for ring in range(nRings):
            for i in range(symmet):#(nodes):
                opp = math.floor((i + (symmet/2))%symmet) # = always i + symmet/2 for even symmet
                x, y = NPCscoordssplit[npc, ring, i, (0, 1)]
                xopp, yopp = xyopp(npc, ring, opp)
                r_op[npc, ring, i] = math.hypot(xopp - x, yopp - y)/2
        
    return r_op #rall

def centroid_r(NPCscoords, z, symmet):
    nRings, n, NPCscoordssplit = NPCcoords_split(NPCscoords, z, symmet)
    
    r_c = np.zeros((NPCscoordssplit[:,:,:,0]).shape)
    
    centroid = np.zeros((n, nRings, 2))
    
    for npc in range(n):
        for ring in range(nRings):
            centroid[npc, ring, :] = np.mean(NPCscoordssplit[npc,ring,:,0:2], axis = 0)
            for i in range(symmet):
                x, y = NPCscoordssplit[npc, ring, i, (0, 1)]
                r_c[npc, ring, i] = math.hypot(centroid[npc, ring, 0] - x, centroid[npc, ring, 1] - y)
    return r_c, centroid 
 
def analyse_opp_cor_radius(rall, z):
     nRings = len(z)
     n = len(rall)
     minr = []
     maxr = []
     for npc in range(n):
         for ring in range(nRings):
             minr.append(min(rall[npc][ring]))
             maxr.append(max(rall[npc][ring]))
     
     ratio = [minr[i]/maxr[i] for i in range(len(maxr))]
     meanratio = [np.mean(ratio[i:i+nRings]) for i in range(0, n*nRings, nRings)]
    
     return minr, maxr, ratio, meanratio
 
 
def splitrings(NPCscoords, i):
    NPC = NPCscoords[NPCscoords[:,3]==i]
    midplane = np.mean(NPC[:,2]) 
    CR = NPC[NPC[:,2] > midplane]
    NR = NPC[NPC[:,2] < midplane]
        
    return CR, NR

def splitRingsByZ(NPCscoords, i):
    NPC = NPCscoords[NPCscoords[:,3]==i]
    
    zs = list(set(NPC[:,2]))
    
    zrings = [None] * len(zs)

    for i in range(len(zs)):
        zrings[i] = NPC[NPC[:,2] == zs[i]]
    
    return zrings, zs


 
def Ellipses(NPCscoords, ringmode = None):
    """ringmode = 'CRNR', 'z' or None"""
    if str(ringmode) not in {"CRNR", "z", "None"}: raise ValueError("Wrong ringmode")
   
    
    n = int(NPCscoords[-1][-1] + 1) # total number of NPCs
    
    
    def squaresum(NPC_i, ellipsefeatures_i, z = None, ring = None):
        _, dist = align2ellipse(NPC_i, ellipsefeatures_i, z, ring)
        return sum([j**2 for j in dist])
    
    def roundfeatures(ellipsefeatures, z = None, ring = None):
        
        if z != None:
            suffix = "_z" + str(z) 
        elif ring != None: 
            suffix = "_" + str(ring)
        else: 
            suffix = ""
        
        
        ellipsefeatures["el_x0" + suffix] = round(ellipsefeatures["el_x0" + suffix], ndigits = 1)
        ellipsefeatures["el_y0" + suffix] = round(ellipsefeatures["el_y0" + suffix], ndigits = 1)
        ellipsefeatures["el_major" + suffix] = round(ellipsefeatures["el_major" + suffix], ndigits = 1)
        ellipsefeatures["el_minor" + suffix] = round(ellipsefeatures["el_minor" + suffix], ndigits = 1)
        ellipsefeatures["el_squaresum" + suffix] = round(ellipsefeatures["el_squaresum" + suffix], ndigits = 1)
        ellipsefeatures["el_q" + suffix] = round(ellipsefeatures["el_q" + suffix], ndigits = 2)
        ellipsefeatures["el_rot" + suffix] = round(ellipsefeatures["el_rot" + suffix], ndigits = 3)
        return ellipsefeatures

    
    if str(ringmode) == "None":
        ellipsefeatures = [None] * n
        squaresums = [None] * n
        
        for i in range(n): 
            NPC_i = NPCscoords[NPCscoords[:,3]==i]
            _, ellipsefeatures[i] = fitEllipse(NPC_i)
            
            squaresums[i] = squaresum(NPC_i, ellipsefeatures[i])
            ellipsefeatures[i].update({"el_squaresum" : squaresums[i]})
            #ellipsefeatures[i] = roundfeatures(ellipsefeatures[i])  TODO: uncomment
            
        return ellipsefeatures
        
    if ringmode == "CRNR":
        ellipseCR = [None] * n
        ellipseNR = [None] * n
        squaresumsCR = [None] * n
        squaresumsNR = [None] * n
        
        CRstr = "CR"
        NRstr = "NR"
        
        for i in range(n):
            CR, NR = splitrings(NPCscoords, i)
            _, ellipseCR[i] = fitEllipse(CR[CR[:,3]==i], ring = CRstr)
            squaresumsCR[i] = squaresum(CR[CR[:,3]==i], ellipseCR[i], ring = CRstr)
            ellipseCR[i].update({"el_squaresum_" + CRstr: squaresumsCR[i]})         
            ellipseCR[i] = roundfeatures(ellipseCR[i], ring = CRstr)
            
            _, ellipseNR[i] = fitEllipse(NR[NR[:,3]==i], ring = NRstr)     
            squaresumsNR[i] = squaresum(NR[NR[:,3]==i], ellipseNR[i], ring = NRstr)
            ellipseNR[i].update({"el_squaresum_" + NRstr: squaresumsNR[i]})
            ellipseNR[i] = roundfeatures(ellipseNR[i], ring = NRstr)
        
        return [ellipseCR, ellipseNR]
    
    
    if ringmode == "z": 
        ring0, _ = splitRingsByZ(NPCscoords, 0)
        nring = len(ring0)
        ellipseAll = [None] * n
       # squaresumsAll = [None] * n
        
        for i in range(n): 
            rings, z = splitRingsByZ(NPCscoords, i)
            ringsinNPC = [None] * nring
            squaresumsRing = [None] * nring
            
            for j in range(nring):
                #z = rings[j][0,2] # z coordinate
                zj = z[j] #TODO: Also for ringmode CRNR. Does not work for varying z. 
                _, ringsinNPC[j] = fitEllipse(rings[j], zj)

                squaresumsRing[j] = squaresum(rings[j], ringsinNPC[j], zj)
                
                ringsinNPC[j].update({"el_squaresum_z" + str(zj) : squaresumsRing[j]})
                ringsinNPC[j] = roundfeatures(ringsinNPC[j], zj)    
                

            ellipseAll[i] = ringsinNPC
            
        return ellipseAll
        

 
def fitEllipse(NPC, z = None, ring = None):

    X = np.array([NPC[:,0]]).T # X coordinates NPC
    Y = np.array([NPC[:,1]]).T # Y coordinates NPC
    
    # Formulate and solve the least squares problem ||A0x - b ||^2
    A0 = np.hstack([X**2, X * Y, Y**2, X, Y])
    b0 = np.ones_like(X)
    

    x = np.linalg.lstsq(A0, b0, rcond = None)[0].squeeze()

    A, B, C, D, E = x
    F = -1

    a = -np.sqrt(2*(A*E**2 + C * D**2 - B*D*E + (B**2 - 4*A*C)*F) *  
                 ((A+C) - np.sqrt((A-C)**2 + B**2))) /  (B**2 - 4*A*C)
    
    b = -np.sqrt(2*(A*E**2 + C * D**2 - B*D*E + (B**2 - 4*A*C)*F) *  
                 ((A+C) + np.sqrt((A-C)**2 + B**2))) /  (B**2 - 4*A*C)
    
    
    x0 = (2*C*D - B*E)/(B**2 - 4*A*C)
    y0 = (2*A*E - B*D)/(B**2 - 4*A*C)

    minor = min(a, b)
    major = max(a, b)


    
    if B != 0:
        rho = np.arctan(1/B * (C - A - np.sqrt((A-C)**2 + B**2)))
    elif A<C:
        rho = 0
    elif A>C:
        rho = 0.5*np.pi    
 
    if A<0:
        rho -= 0.5*np.pi
        if rho < -0.5 * np.pi: rho += np.pi 

    q = minor/major  
    
    
    x0key = "el_x0"
    y0key = "el_y0"
    minorkey = "el_minor"
    majorkey = "el_major"
    qkey = "el_q"
    rotkey = "el_rot"
    
    keyList = [x0key, y0key, minorkey, majorkey, qkey, rotkey]
    valueList = [x0, y0, minor, major, q, rho]
    

    if z != None: 
        keyList = [i + "_z" + str(z) for i in keyList]
    elif ring != None:
        keyList = [i + "_" + str(ring) for i in keyList]
        
    zip_iterator = zip(keyList, valueList)

    return x, dict(zip_iterator) #{x0key : x0, y0key : y0, minorkey : minor, majorkey : major, qkey : q, rotkey : rho}
            
     

def dist2ellipse(semi_major, semi_minor, p):  
    #https://github.com/0xfaded/ellipse_demo/issues/1
    # project all points onto 1st quadrant 
    px = abs(p[0])
    py = abs(p[1])
    
    # tx = 0.707 # = cos(pi/4) = sqrt(0.5)
    # ty = 0.707
    t = math.pi/4
    
    a = semi_major
    b = semi_minor
    
    for _ in range(0, 3): # 3 iterations 
        # x = a * tx
        # y = b * ty
        x = a * math.cos(t)
        y = b * math.sin(t)
    
        # ex = (a*a - b*b) * tx**3 / a
        # ey = (b*b - a*a) * ty**3 / b
        ex = (a*a - b*b) * math.cos(t)**3 / a
        ey = (b*b - a*a) * math.sin(t)**3 / b
    
        rx = x - ex
        ry = y - ey
    
        qx = px - ex
        qy = py - ey
    
        r = math.hypot(ry, rx)
        q = math.hypot(qy, qx)
    
        # tx = min(1, max(0, (qx * r / q + ex) / a))
        # ty = min(1, max(0, (qy * r / q + ey) / b))
        # t = math.hypot(ty, tx)
        # tx /= t 
        # ty /= t 
        
        delta_c = r * math.asin((rx*qy - ry*qx)/(r*q))
        delta_t = delta_c / math.sqrt(a*a + b*b - x*x - y*y)

        t += delta_t
        t = min(math.pi/2, max(0, t))
    
    #closestpoint = (math.copysign(a * tx, p[0]), math.copysign(b * ty, p[1]))
    closestpoint = (math.copysign(a * math.cos(t), p[0]), math.copysign(b * math.sin(t), p[1]))
    return closestpoint, np.linalg.norm(p - closestpoint)
 

def rotate(origin, point, angle):
    """
    Rotate a point counterclockwise by a given angle around a given origin.
    The angle should be given in radians.
    """
    ox, oy = origin
    px, py = point


    qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
    qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
    return qx, qy

def centre0(NPCrot):
    """centre around 0"""
    NPCrot[:,0] -= np.mean(NPCrot[:,0])
    NPCrot[:,1] -= np.mean(NPCrot[:,1]) 
    return NPCrot



def align2ellipse(NPC, elfeatures, z = None, ring = None):    
    #elfeatures = elfeatures[0] 
#    "x0" : x0, "y0" : y0, "minor" : minor, "major" : major, "q": q, "rho" : rho}

    if z != None:
        suffix = "_z" + str(z) 
    elif ring != None: 
        suffix = "_" + str(ring)
    else: 
        suffix = ""
        

    minor = elfeatures["el_minor" + suffix]
    major = elfeatures["el_major" + suffix]
    rho = elfeatures["el_rot" + suffix]
    x0 = elfeatures["el_x0" + suffix]
    y0 = elfeatures["el_y0" + suffix]

    nodes = len(NPC) 
    NPCrot = np.array(NPC)

    #rotate around centre of fitted ellipse so that the major axis aligns with 
    # the y axis 
    for i in range(nodes): 
        NPCrot[i][0:2] = rotate([x0, y0], NPC[i][0:2], -rho)
        
    NPCrot = centre0(NPCrot) #offset so that the centre is at 0 
  
    dist = [None] * nodes
    closestpoints = [None] * nodes    

    for i in range(nodes):
        closestpoints[i], dist[i] = dist2ellipse(major, minor, NPCrot[i][0:2])
        
    return closestpoints, dist
    


#########################################################################

def Circles(NPCscoords, ringmode = None):
        """ringmode is "CRNR", "z" or None"""
        if str(ringmode) not in {"CRNR", "z", "None"}: raise ValueError("Wrong ringmode")
        
        n = int(NPCscoords[-1][-1]+1)
        
        if str(ringmode) == "None":
            circlefeatures = [None] * n
            
            for i in range(n): 
                circlefeatures[i] = fitCircle(NPCscoords[NPCscoords[:,3]==i]) 
                
            return {"circlefeatures": circlefeatures}
            
        if ringmode == "CRNR": 
            circleCR = [None] * n
            circleNR = [None] * n
            for i in range(n):
                CR, NR = splitrings(NPCscoords, i)
                circleCR[i] = fitCircle(CR[CR[:,3]==i])
                circleNR[i] = fitCircle(NR[NR[:,3]==i])
            return {"circleCR": circleCR, "circleNR" : circleNR}
        
        if ringmode == "z":
            
            ring0, _ = splitRingsByZ(NPCscoords, 0)
            nring = len(ring0)
            circleAll = [None] * n
            
            for i in range(n):
                
                rings, _ = splitRingsByZ(NPCscoords, i)
                ringsinNPC = [None] * len(rings)
                for j in range(nring):
                    
                    ringsinNPC[j] = fitCircle(rings[j])
                circleAll[i] = ringsinNPC
            return {"circleAll" : circleAll}
        


def fitCircle(NPC):     
    xc,yc,r,residual = cf.least_squares_circle(NPC) #xcentre, ycentre, radius
    xc = round(xc,  ndigits=1)
    yc = round(yc, ndigits=1)
    #r = round(r, ndigits=1) TODO
    residual = round(residual, ndigits=1)
    return [xc, yc, r, residual] #[xc, yc, np.round(r, 2), residual] # TODO: change back
 
