#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 14:36:41 2021

@author: maria
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import circulant
from warnings import warn
from sklearn.gaussian_process.kernels import RBF
import random
import sys
from Nups_Info import SelectNup


class DeformNPC:
    def __init__(self, nConnect = 2, mag = 25, symmet = 8, nRings = 1, r = 0, ringAngles = 0, z = 0, elliptical = False, sigmamult = 1, seed = None):
        '''
        Models deformed NPCs using solve_ivp based on a simple, rotationally symmetric node and spring model with deforming forces applied in xy direcion. 
        ### Input ###
        nConnect: Number of connected neighbour nodes in clock-wise and anti-clockwise direction
        mag: Magnitude of deformation [nm]. Number represents 3 standard deviations -> 99.7 % of forces on a node lie within this range
        symmet: Symmetry of the NPC (Default 8 )
        nRings: Number of NPC rings 
        r: Radius of NPC rings assuming 8-fold symmetry. Must be a list of length nRings
        ringAngles: Rotational angle offset of NPC rings. Must be a list of length nRings
        z: z position of NPC rings. Stays unchanged. Must be a list of length nRings
        
        ### Relevant output ###
        solution: solve_ivp output for all NPC rings 
        fcoords: Coordinates of force vectors for all NPC rings
        initcoords: starting coordinates 
        randfs: list of magnitude of forces for all NPC rings 
        r: Radius of NPC rings corrected for symmetry (i.e. increasing with symmetry)
        z: z position of NPC rings
       
        '''
        self.seed = seed
        self.solution = [] 
        self.symmet = symmet # Rotational symmetry of the NPC
        self.initcoords = [] # starting coordinates of nodes    
        self.fcoords = []    # Coordinates of forces 
        self.z = z # z coordinates (stays unchanged)
        self.elliptical = elliptical
        self.geodesic = np.zeros(len(ringAngles))

        
        damp = 1 # damping :TODO
        kr = 0.7 # spring constant of anchor spring 
        
        tlast = 30
        tspan = [0,tlast]      
        #teval = np.arange(0, tlast, 0.4) # needed for animation function
        teval = None        
        
        Lrests = [] # circulant matrix of resting spring lengths
        Ks = [] # circulant matrix of spring constants 
        y0s = [] # initial coordinates and velocities per node

        if(nConnect > self.symmet/2):
            nConnect = int(np.floor(self.symmet/2))
            warn("Selected number of neighbours nConnect too large. nConnect has been changed to " + str(nConnect) + ".")
        
        if(len(r) != nRings or len(ringAngles) != nRings or len(z) != nRings):
            warn("r, ringAngles, and z must be of length nRings: " + str(nRings))

        self.geodesic = np.asarray(ringAngles) * np.mean(r) # calculate geodesic distance of corners (for 8-fold NPC)   
        self.ringAngles_corrected = np.zeros(len(ringAngles))
        
        #correct radius 
        mean_r_corrected = self.adjustRadius(np.mean(r))       
        diff_r = r - np.mean(r)
        r_corrected = diff_r + mean_r_corrected

            
        for i in range(nRings):  
            self.ringAngles_corrected[i] = self.geodesic[i]/np.mean(r_corrected)#TODO: r_corrected[i] # correct geodesic distance
            initcoord = self.initialcoords(r_corrected[i], ringAngle = self.ringAngles_corrected[i]) # TODO
            self.initcoords.append(initcoord) 
            Lrest = self.springlengths(initcoord)
            Lrests.append(Lrest) 
            Ks.append(self.springconstants(Lrest))
            y0s.append(np.concatenate((initcoord.flatten(), np.zeros(2 * self.symmet)))) 
 
        self.r = r_corrected    
        sigma = np.min(self.r)*sigmamult # sigma: free parameter for RBF kernel
        
        self.randfs = self.forcesMultivariateNorm(self.initcoords, r_corrected, mag, sigma = sigma, nRings = nRings) # generate random forces
        # Solve ODE, ring 1 - 4       
        for i in range(nRings):
            self.solution.append(solve_ivp(self.npc, tspan, y0s[i], t_eval=teval, 
                                           method='RK45', args=(r_corrected[i], Lrests[i], 
                                        Ks[i], kr, self.randfs[i], damp, nConnect)))   
            
            self.fcoords.append(self.initialcoords(r_corrected[i], self.ringAngles_corrected[i], self.randfs[i])) # TODO
        

            
    ### Methods ###


    def npc(self, t, y, r, Lrest, K, kr, randf, damp, nConnect):
        '''
        t: time points 
        y: values of the solution at t 
        r: radius of NPC ring
        Lrest: Circulant matrix of resting lengths of all springs 
        K: Circulant matrix of all radial spring constants 
        kr: Spring constants of anchor springs 
        randf: array of forces (length = symmet) to be applied in radial direction to each node 
        damp: Damping factor 
        nConnect: Number of connected neighbours in cw and ccw direction for each node 

        output: solutions at t. x and y components of positions and velocities of each node for each timestep 
        '''
        v = np.reshape(y[2*self.symmet:], (self.symmet, 2))
        x = np.reshape(y[:2*self.symmet], (self.symmet, 2))
        
        anc = np.array([0., 0.]) # coordinates of anchor node   
        F = np.zeros((self.symmet, 2)) # Forces
        

        for i in range(self.symmet): 
            F[i] = randf[i]*x[i] / np.linalg.norm(x[i] - anc) 

        allaccarray = np.zeros((self.symmet, 2)) # array for accelerations of node 0 - 7
        
        
        
        for i in range(self.symmet): # i indicates the reference node        
            accarray = np.array([0., 0.]) # initiate acceleration array for each node i 
            
            for j in [k for k in range(-nConnect, nConnect+1) if k != 0]: # j is neighbour nodes -nConnect to +nConnect relative to i, skipping 0 (0=i)            
                jnew = (i+j)%self.symmet 
                accarray += K[i][jnew]  * (x[i]-x[jnew])/np.linalg.norm(x[i]-x[jnew]) * (Lrest[i][jnew] - np.linalg.norm(x[i]-x[jnew]))
    
            accarray += kr * (x[i] - anc)/np.linalg.norm(x[i] - anc) * (r - np.linalg.norm(x[i] - anc)) #anchor
            accarray = F[i] + accarray - damp*v[i]  # external force and damping
            allaccarray[i] = accarray 
    
        dxdt = np.concatenate((v.flatten(), allaccarray.flatten()))                                                                
        return dxdt
    
    
    def adjustRadius(self, r8):
        """Adjusts radius r with symmetry. No adjustment is made when symmetry is 8. Radius is viewed
        as the length of the symmetric side of an isoceles triangle whose tip (angle alpha) points towards the 
        center of the NPC and whose base is the section between two neighbouring nodes at the circumference. Think slice of cake.
        # Input: 
        r8: radius of a default 8-fold symmetrical NPC
        
        ## Output: 
        radius of NPC with symmetry equal to symmet (rnew = r8 if symmet = 8)
        
        """
        alpha = 2*np.pi/self.symmet # Angle at the tip of triangular slice (pointing to center of NPC)
        theta = 0.5 * (np.pi - alpha) # Either angle at the base of (isosceles) triangular slice
        halfbase = r8 * np.sin(np.pi/8) # half the distance between two corners of an NPC ring (= half the base of triangular slice)
        return halfbase/np.cos(theta) # new radius
      
    def pol2cart(self, rho, phi):
        '''Transforms polar coordinates of a point (rho: radius, phi: angle) to 2D cartesian coordinates.
        '''
        x = rho * np.cos(phi)
        y = rho * np.sin(phi)
        return(x, y)

    def cart2pol(self, x, y): 
        rho = np.sqrt(x**2 + y**2)
        phi = np.arctan2(y, x)
        return(rho, phi)    
    
    def initialcoords(self, r, ringAngle = 0, forces = 0): 
        '''
        Generates cartesian coordinates of the NPC given radius and self.symmet 
        ## Input ##
        r: NPC Radius
        ringAngle: Angular offset of ring (default 0)
        forces (optional): input 1D forcevector to view coordinates of forces; used for visualisation

        ## Return values ##
        2D Cartesian coordinates 
        '''
        
        if(type(forces) != int):
            if (len(forces) != self.symmet):
                warn("forces must be 0 or an array with len(self.symmet")
                
        #correct ringAngle for varying symmetries #TODO: test
        # ringAngle = ((ringAngle*(2*np.pi/self.symmet))
        #             /(2*np.pi/8))
        #ringAngle = ringAngle * 8/self.symmet
        
        forces = forces * np.ones(self.symmet) # forces is 0 or np.array with len(self.symmet)
        rotAngle = 0.

        initcoord = np.zeros((self.symmet, 2))
        
        for i in range(self.symmet):
            initcoord[i, 0], initcoord[i, 1] = self.pol2cart(r + forces[i], 
                                                             rotAngle+ringAngle)
            rotAngle += 2*np.pi/self.symmet

        return initcoord
    
    
    def springlengths(self, initcoord): 
        '''Compute lengths of springs from coordinates and returns circulant matrix
        '''

        l = np.zeros(self.symmet)
        for i in range(self.symmet):
            l[i] = np.linalg.norm(initcoord[0, :] - initcoord[i, :])      
        self.L = circulant(l)
        return circulant(l)
    
    
    def springconstants(self, Lrest):
        "Returns circulant matrix of spring constants "
        
        Lscale = Lrest/Lrest[0][1] # scaled so that shortest edge is 1
        
        k = np.ones(int(np.floor(self.symmet/2)))
        if(self.symmet%2 == 0): #if symmet is even
            k[-1] = k[-1]/2 # springs that connect opposite corners will be double-counted. Spring constant thus halved 
            K = circulant(np.append(0, np.append(k, np.flip(k[:-1]))))
        else: #if symmet is odd 
            K = circulant(np.append(0, [k, np.flip(k)]))
        
        K = K/Lscale
        
        self.K = K
        return K
    
    
    def forcesMultivariateNorm(self, initcoords, r, mag, sigma, nRings = 1):
        '''
        Returns array of Forces that are later applied in radial direction to the NPC corners
        ## Input ## 
        *coordring: Initial coordinates of nodes for an arbitrary number of rings. 
        mag: Total magnitude of distortion. Default is 50. 
        r: radius of all NPC rings 
        nRings: Number of NPC rings (default 1)
        ## Returns ## 
        For each ring, an array of forces applied to each node
        '''

        nodesTotal = self.symmet*nRings # total number of nodes over all rings 
        #AllD = np.zeros((nodesTotal, nodesTotal))       
        
        zcoords = np.zeros(nodesTotal)
        for i in range(nodesTotal):
            zcoords[i]=self.z[int(np.floor(i/self.symmet))]
         
        # scale z to sigma
        zcoords = (zcoords/max(zcoords)) * sigma 
        allcoords = np.asarray(initcoords) # separated by NPC rings
        allcoords = allcoords.reshape(nodesTotal, 2) # not separated by NPC rings 
        
        if self.elliptical:
            Fe = self.ellipse(allcoords, self.elliptical)
        
        # project all coordinates onto one circle so that
        # varying radii keep their shape relative to each other better 
        phi = np.arctan2(allcoords[:,1],allcoords[:,0])
        allcoords[:,0] = np.mean(r)*np.cos(phi) 
        allcoords[:,1] = np.mean(r)*np.sin(phi)

        if (mag == 0):
            F = np.zeros(nodesTotal)
        else:    
            
            allcoords3D = np.concatenate((allcoords.T, [zcoords])).T
            #sigma = np.min(r)/2 # free parameter of RBF kernel that outputs covariance of forces on two nodes based on their euclidean distance
            cov = np.zeros((nodesTotal, nodesTotal))
    
            kernel = 1.0 * RBF(sigma)
            cov = kernel.__call__(X = allcoords3D)
    
            var = (mag/3)**2    # 3*SD to Var TODO
            cov = var * cov
            cov = cov + np.identity(nodesTotal)*1e-6 # jitter matrix to surely make it positive definite

            rng = np.random.default_rng(seed = self.seed) 
            u = rng.standard_normal(len(allcoords))
            L = np.linalg.cholesky(cov)
            F = L @ u
            

        
        if self.elliptical:
            F = F + Fe
        return np.split(F, nRings)
    
    def ellipse(self, allcoords, elliptical):
        p = [i*2*np.pi for i in self.r] # perimeter of NPC circle
        
        q = elliptical # ratio of a to b
    
        a = [((np.sqrt(2)*i)/(2*np.pi))/np.sqrt(1+q**2) for i in p] #length of axis a given perimeter of the ellipse is rougly p
        b = [q*i for i in a]
        
        nodes = len(allcoords)
        Fe = np.zeros(np.shape(allcoords)) # forces to deform NPC into ellipse 
  
        polarcoords = np.zeros(np.shape(allcoords))  
        n = np.zeros(nodes) # scaling factors n[i] that transform allcoords into ellipse coords
        ellipsecoords = np.zeros(np.shape(allcoords))

        # rotate NPC ring randomly between -180 and +180 deg
        for i in range(nodes):
            polarcoords[i] = self.cart2pol(allcoords[i,0], allcoords[i, 1])  
        
        rng = np.random.default_rng(seed = self.seed) 
        rotateby = rng.uniform(-np.pi, np.pi)    
        
        polarcoords[:,1] += rotateby
        
        for i in range(nodes):  
            ring = int(np.floor(i/self.symmet))
            allcoords[i] = self.pol2cart(polarcoords[i, 0], polarcoords[i, 1]) # rotated NPC ring back to cartesian coords
            n[i] = np.sqrt((a[ring]**2 * b[ring]**2)/
                           (b[ring]**2 * allcoords[i, 0]**2 + 
                            a[ring]**2 * allcoords[i, 1]**2)) # scaling factors n[i] that transform allcoords into ellipse coords
            ellipsecoords[i] = allcoords[i] * n[i]

        sign = np.zeros(nodes)
        difference = ellipsecoords - allcoords

        for i in range(nodes):
            sign[i] = np.sign(np.dot(difference[i], allcoords[i]))  # -1 if >90 deg, +1 if < 90 deg
       
        Fe = sign * np.linalg.norm(difference, axis = 1)
        return Fe


# Multiple NPCs 
def MultipleNPC(nup, term, model, n = 1, mag = 0, sigmamult = 0.5, nConnect = 4, symmet = 8, elliptvar = None, rvar = None, dvar = None, thetavar = None, seed = None):
    "Generate n deformed NPCs"
    
    if seed:
        random.seed(seed)
        seeds = random.sample(range(min(2**32, sys.maxsize)), n)#, n)
    else: seeds = np.full(n, None)
    
    
    r = SelectNup(nup, term, model).r
    rold = np.array(r)

    z = SelectNup(nup, term, model).z
    zold = np.array(z)
    
    ringAngles = SelectNup(nup, term, model).ringAngles
    ringAnglesOld = np.array(ringAngles)
    
    thetaold = SelectNup(nup, term, model).rotAng # NR - CR
    thetaoffset = SelectNup(nup, term, model).rotOffset #0, -45 or +45. Only works for 8-fold symmetry
    #sigma = min(r)*sigmamult
    NPCs = []
    zNPCs = []
    fcoords = []
    randfs = []
    nRings = len(z) # TODO: Problems for nRings = 1
    
    
    for i in range(n):
        r = Change_radius(rold,  rvar["rnew"], rvar["rsigma"], seeds[i])
        z = Change_dist(zold, dvar["dnew"], dvar["dsigma"], seeds[i])
        ringAngles, newTheta = Change_rotang(ringAnglesOld, thetaold, thetaoffset, zold, thetavar["thetanew"], thetavar["thetasigma"], seeds[i])
        elliptical = Change_ellipt(elliptvar["elliptnew"], elliptvar["elliptsigma"], seeds[i])
        

        deformNPC_temp = DeformNPC(nConnect, mag, 
                symmet = symmet, nRings = nRings, 
                r = r, ringAngles = ringAngles, 
                elliptical= elliptical, z = z, 
                sigmamult= sigmamult, seed = seeds[i])#TODO: remove


        fcoord = deformNPC_temp.fcoords # list
        tempsolution = deformNPC_temp.solution # list
        tempz = deformNPC_temp.z #np array
        randf = deformNPC_temp.randfs # list
        
        NPCs.append(tempsolution) 
        zNPCs.append(list(tempz))

        fcoords.append(fcoord)
        randfs.append(randf)

    ringAngles_corrected = deformNPC_temp.ringAngles_corrected

    multipleNPCs = {"NPCs" : NPCs, "zNPCs" : zNPCs, "mag": mag, 
     "sigmamult" : sigmamult, "z" : z, "fcoords" : fcoords,
     "temp_r" : deformNPC_temp.r, "ringAngles" : ringAngles, 
     "ringAngles_corrected" : ringAngles_corrected, 
     "symmet": symmet, "elliptical" : elliptical, 
     "randfs" : randfs, "nConnect" : nConnect, 
     "thetaold" :thetaold, "theta_offset" : thetaoffset}
        
    return multipleNPCs

def Change_radius(r, rnew = False, rsigma = False, seed = None):
    """
    r: old radii of all rings, np.array
    rnew: new mean radius or False
    sigma: standard deviation of Gaussion of which to sample new mean radii from. 
    Gaussian is centred around rnew, if rnew is provided, and around np.mean(r) otherwise
    """

    rmean = rnew if rnew else np.mean(r)

    if rsigma: 
        np.random.seed(seed)
        rmean = np.random.normal(rmean, rsigma)
    
    return (rmean/np.mean(r)) * r  


def Change_dist(zold, dnew = False, dsigma = False, seed = None):
    z = np.array(zold)
    
    midplane = np.mean(z) 
    dist = np.mean(z[z > midplane])-np.mean(z[z < midplane])
    
    dnew = dnew if dnew else dist
    
    if dsigma:
        np.random.seed(seed)
        dnew = np.random.normal(dnew, dsigma)
    
    ddif = dist - dnew
    
    z[z>midplane] = z[z>midplane]-ddif
    return z

def Change_rotang(ringAnglesOld, thetaold, thetaoffset, zold, thetanew = False, thetasigma = False, seed = None):
    midplane = np.mean(zold)
    ringAngles = np.array(ringAnglesOld)

    if not (str(thetanew) == "0" or str(thetanew) == "0.0"):
        thetanew = thetanew if thetanew else thetaold
    
    if thetasigma:
        np.random.seed(seed)
        thetanew = np.random.normal(thetanew, thetasigma)
    
    thetadif = thetaold - thetanew
    ringAngles[zold > midplane] += thetadif/2
    ringAngles[zold < midplane] -= thetadif/2
    
    NRmean = np.mean(ringAngles[zold < midplane])
    CRmean = np.mean(ringAngles[zold > midplane])
    
    newminTheta = NRmean - CRmean - thetaoffset
    return ringAngles, newminTheta

def Change_ellipt(elliptnew = False, elliptsigma = False, seed = None):
    
    elliptnew = elliptnew if elliptnew else 1.
    if elliptsigma:
        np.random.seed(seed)
        elliptnew = np.random.normal(elliptnew, elliptsigma) 
    return elliptnew        

def Sol2D(solution):
    """
    input: DeformNPC.solution for a given NPC ring
    output: 2D arrays of position and of velocity of nodes in a ring over time [timeframe, node, dimension (x or y)]
    """
    nFrames = len(solution.t)
    symmet = int(len(solution.y)/4) #/4 because of 2 dim times both velocity and position 
    pos1D = solution.y.T[: , :2*symmet] # positions over time
    vel1D = solution.y.T[: , 2*symmet:] # velocities over time
    pos2D = np.reshape(pos1D, (nFrames, symmet, 2))
    vel2D = np.reshape(vel1D, (nFrames, symmet, 2))
    return pos2D, vel2D


def Pos2D(solution):
    pos2D, vel2D = Sol2D(solution)
    return pos2D   


def MultipleNPCs_coord(NPCs, zNPCs, symmet):
    "Input is OdeResults for each NPC in a list. Output is just the final coordinates of each NPC"
    nRings = len(zNPCs[0])
    NPCscoord = np.zeros((symmet*nRings*len(NPCs), 4)) # number of nodes, dimensions + label
 
    i = 0

    for NPC in range(len(NPCs)):        
        for ring in range(nRings):
            for node in range(symmet):
                NPCscoord[i] = np.append(Pos2D(NPCs[NPC][ring])[-1, node], [zNPCs[NPC][ring], NPC])
                i += 1

    return NPCscoord
