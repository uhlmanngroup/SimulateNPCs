#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 16:21:05 2021

@author: maria
"""
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import seaborn as sns
from warnings import warn
import matplotlib.animation as animation
from matplotlib.patches import Circle
import DeformNPC
import Analyse_deformed
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap

Pos2D = DeformNPC.Pos2D
Sol2D = DeformNPC.Sol2D


def OffsetNPCs(NPCcoords, NPCs, maxr): 
    "Arrange NPCs on a grid by offsetting them a multiple of their radius maxr in x and y direction"
    
    NPCoffset = np.copy(NPCcoords)
    
    n = len(NPCs)
    # Determine the number of rows and columns needed. The last cells on the grid might stay empty
    ncols = math.ceil(np.sqrt(n)) 
    nrows = math.ceil(n/ncols)
    
    x = 0 #indexing x coordinate
    y = 1 # indexing y coordinate
    i = 0 # will get updated

    for row in range(ncols):
        for col in range(nrows):    
            if (i < n):
                NPCoffset[np.where(NPCoffset[:,3] == i), y]  += col*3*maxr # TODO: switch row and col back 
                NPCoffset[np.where(NPCoffset[:,3] == i), x]  += row*3*maxr          
                i += 1
    return NPCoffset


def OverviewPlot(NPCcoords, NPCs, mag, r, ellipse = False, circle = False):
    maxr = max(r)
    NPCoffset = OffsetNPCs(NPCcoords, NPCs, maxr)
    
    if len(NPCs) == 1: 
        markersize = 20
    else: markersize = 5

    # prepare to colourcode z position
    zs = NPCoffset[:,2]
    
    zcolour = []    
    zcolour.extend(ColourcodeZ(list(zs)))  
    mincolour = float(min(zcolour))
    maxcolour = float(max(zcolour))
    
    # plot
    plt.rcParams.update({'font.size': 50})
    fig, ax = plt.subplots(1, 1, figsize = (25, 25))
    ax.set_title("mag " + str(mag))
    ax.scatter(NPCoffset[:,0], NPCoffset[:,1], c = zcolour, s = np.min(r)*markersize)
    
    
    # fit circle and/or ellipse
    if ellipse:
        for i in range(len(NPCs)): plotEllipse(NPCoffset[NPCoffset[:,3]==i], ax, len(NPCs))
    
    if circle: 
        for i in range(len(NPCs)): 
            plotCircle(NPCoffset[NPCoffset[:,3]==i], ax, len(NPCs))
    
    # colourbar for z position 
    incr = (maxcolour - mincolour)/max(zs) # increment in colour 
    cmap = (ListedColormap( [str(i) for i in list(np.arange(mincolour, maxcolour + 0.5*incr, incr))])) 
    norm = plt.Normalize(min(zs), max(zs))
    fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), shrink = 0.7, label = "z [nm]", ticks = [min(zs), max(zs)])


    ax.set(xlabel = "x [nm]", ylabel = "y [nm]")
    ax.axis("scaled") 
    fig.tight_layout()


def plotEllipse(NPC, ax, n):
    
        X = NPC[:,0]
        Y = NPC[:,1]
        
        x, values = Analyse_deformed.fitEllipse(NPC)
        
        # Plot the least squares ellipse
        x_coord = np.linspace(min(X)-10, max(X)+10, 20)
        y_coord = np.linspace(min(Y)-10, max(Y)+10, 20)
        X_coord, Y_coord = np.meshgrid(x_coord, y_coord)
        Z_coord = x[0] * X_coord ** 2 + x[1] * X_coord * Y_coord + x[2] * Y_coord**2 + x[3] * X_coord + x[4] * Y_coord
        ax.contour(X_coord, Y_coord, Z_coord, levels=[1], colors=('m'), linewidths=4)
        
        # print infotext 
        if n <= 16:
            fontsize = 25
            if n <= 9: fontsize = 35
            if n <= 4: fontsize = 40
            info = ("a = " + str(round(values["el_major"], ndigits = 1)) + "\n" + "b/a = " + 
            str(round(values["el_q"], ndigits = 2)) + "\n"+ "rho = " + str(round(math.degrees(values["el_rot"]), ndigits = 1)))        
            ax.text(np.mean(NPC[:,0])-25, np.mean(NPC[:,1])-10, info, fontsize = fontsize)


def plotCircle(NPC, ax, n):
        xc, yc, r, _ = Analyse_deformed.fitCircle(NPC) # x and y position of centre, radius 
        c = Circle((xc,yc), radius = r, facecolor='none', edgecolor = "b", linewidth=4)
        ax.add_artist(c)



def ColourcodeZ(z, darkest = 0.1, brightest = 0.8):
    '''colourcode z, smaller values are shown darker'''
    return [str(i) for i in np.interp(z, (min(z), max(z)), (darkest, brightest))]

    
def Plot2D(solution, z, symmet, nConnect,  linestyle = "-", trajectory = True, 
           colourcode = True, springs = True, anchorsprings = True, markersize = 20, 
           forces = 0, showforce = False, legend = False): # TODO 
    '''
    solution: Output of solve_ivp
    symmet: number of nodes per ring
    nConnect: number of neighbours connected on each side per node
    linestyle (default: "-"): Linestyle in 1st plot 
    legend (default: False): Whether to show a legend in the 1st plot 
    colourcode (default: True): colourcodes trajectory with velocity
    colourbar (default: True): Plots colourbar in 2nd plot if True and if colourcode is True
    mainmarkercolor: Colour of nodes in 2nd plot 
    '''
        
    nRings = len(z)
    plt.rcParams.update({'font.size': 25})

    fig, ax = plt.subplots(1, 1, figsize = (12, 12))
    viewFrame = -1 # 0 is the first frame, -1 is the last frame  
    mainmarkercolor = ColourcodeZ(z)
    
    for i in range(nRings):
        
        nFrames = len(solution[i].t)# Nodes at last timestep
        pos2D, vel2D = Sol2D(solution[i])
        
        ax.plot(pos2D[viewFrame, :symmet, 0], pos2D[viewFrame,:symmet, 1], 
        linestyle = "", marker = "o", color="gray", markerfacecolor = mainmarkercolor[i], 
        markersize = markersize, zorder = 50, label = str(round(z[i], 1)))
        
        if (anchorsprings):
            # Anchor springs
            ax.plot([0,0], [0,0], marker = "o", color = "lightgray", markersize = 15)
            for j in range(symmet):
                ax.plot((pos2D[viewFrame, j, 0], 0), (pos2D[viewFrame, j, 1], 0),
                linestyle = ":", marker = "", color="lightgray")   
            
        # circumferential springs 
        if(springs):
            for ni in range(1, nConnect+1): # neighbours to connect to
                for j in range(symmet): # node to connect from 
                    ax.plot(pos2D[viewFrame, (j, (j+ni)%symmet), 0], pos2D[viewFrame, (j, (j+ni)%symmet), 1], 
                    linestyle = ":", marker = "", color="gray")#, linewidth = 5)
        
        # Trajectories 
        if (trajectory):
            if (colourcode): # Colourcoded trajectory
                ### colourcoding velocities
                normvel = np.zeros((nFrames, symmet)) #nFrames, node
                
                for j in range(symmet):
                    for frame in range(nFrames):
                        normvel[frame, j] = np.linalg.norm([vel2D[frame, j, 0], vel2D[frame, j, 1]])
                        
                norm = plt.Normalize(normvel.min(), normvel.max()) 
                
                #####trajectory colorcoded for velocity        
                for j in range(symmet):
                    points = pos2D[:, j, :].reshape(-1, 1, 2)
                    segments = np.concatenate([points[:-1],points[1:]], axis = 1)
                    lc = LineCollection(segments, cmap = 'plasma', norm=norm, zorder = 100)
                    lc.set_array(normvel[:, j])
                    line = ax.add_collection(lc) # TODO will only be saved for the last ring 
                       
            else: # monochrome trajectory
                for j in range(symmet):
                    ax.plot(pos2D[:, j, 0], pos2D[:, j, 1], color = "blue", linestyle = "-")
    
        ### Force vectors
        if(showforce and type(forces) != int):
            #forces2d = forces.reshape(symmet, 2)
            forces2d = forces[i]
            for j in range(symmet):
                ax.arrow(x = pos2D[0, j, 0], y = pos2D[0, j, 1], 
                             dx = (forces2d[j, 0] - pos2D[0, j, 0]), 
                             dy = (forces2d[j, 1] - pos2D[0, j, 1]),
                             width = 0.7, color="blue")   


    if(trajectory and colourcode and legend):
        fig.legend(bbox_to_anchor=(0.1,-0.025), loc="lower left", ncol = 4, title = "z [nm]")#loc="best")
        axcb = fig.colorbar(line, ax=ax, shrink = 0.7, ticks = [0])   
        axcb.set_label('velocity (a.u.)')
            
    ax.axis("scaled")
    ax.set(xlabel = "x (nm)", ylabel = "y (nm)")  
    plt.tight_layout()




def XYoverTime(solution, symmet , legend = False): #TODO: remove nRings
    '''x and y positions over time'''
    
    nRings = len(solution)
    
    # determin number of rows and colums in final figure. One plot per NPC ring
    l = 2-nRings%2
    rows, cols = sorted((int((nRings/l)), l))

    fig, ax = plt.subplots(rows, cols, figsize = (10*rows, 10*cols))
    palette = sns.color_palette("hsv", 2*symmet)

    for ring in range(nRings):
        for i in range(symmet):

            ax = ax.flatten()
            labelx = labely = None
            if ring == 0: 
                labelx ="x" + str(i) 
                labely = "y" + str(i)
                
            ax[ring].plot(solution[ring].t, Pos2D(solution[ring])[:, i, 0], label = labelx, linestyle = "-", color = palette[i*2])
            ax[ring].plot(solution[ring].t, Pos2D(solution[ring])[:, i, 1], label = labely, linestyle = "--", color = palette[i*2])
            ax[ring].set_title("ring "+ str(ring))
            ax[ring].set(xlabel = 't (a.u.)')
            ax[ring].set(ylabel = "change in x or y [nm]")
                                                                                                                  
    if(legend): fig.legend(bbox_to_anchor=(1,0.9), loc="upper left")
    plt.tight_layout()
    
    plt.show()
    



def Plot3D(solution, z, symmet, randfs, fcoords, plotforces = False, viewFrame = -1, colour = ["black", "gray"]):
    '''viewFrame: 0 is first frame, -1 is last frame'''
    fig = plt.figure(figsize = (15,15))
    ax = fig.add_subplot(111, projection='3d')

    linewidth = 3
    nRings = len(z)
    colour = ColourcodeZ(z)
    
    for ring in range(nRings):
        ax.scatter(Pos2D(solution[ring])[viewFrame, : ,0], Pos2D(solution[ring])[viewFrame, :,1], z[ring], s = 300, c = str(colour[ring]), linewidths = linewidth,  marker = "o")
        
    if plotforces:    
        for ring in range(nRings):
            for node in range(symmet):
                ax.quiver(Pos2D(solution[ring])[0, node, 0], Pos2D(solution[ring])[0, node ,1], z[ring], fcoords[ring][node][0], fcoords[ring][node][1], 0, length = randfs[ring][node], normalize = True, linewidth = linewidth , edgecolor = "blue") 
    
    ax.set_xlabel('x [nm]', labelpad = 30)
    ax.set_ylabel('y [nm]', labelpad = 30)
    ax.set_zlabel('z [nm]', labelpad = 40, fontsize = 40)

    ax.set_xticks([-25, 25])
    ax.set_yticks([-25, 25])
    ax.set_zticks([0, 50])
    
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    
    # color of edges that aren't axes 
    ax.xaxis.pane.set_edgecolor('w')
    ax.yaxis.pane.set_edgecolor('w')
    ax.zaxis.pane.set_edgecolor('w')

    ax.grid(False)
    plt.show()

#%matplotlib inline
    
class AnimatedScatter(object):
    """An animated scatter plot using matplotlib.animations.FuncAnimation."""
    def __init__(self, solution, nConnect, symmet, r, randfs):
        self.solution = solution
        #self.nRings = nRings
        self.nRings = len(solution)

        framenumbers = []
        
        for i in range(self.nRings): # check framenumbers are consistent for each ring
            framenumbers.append(len(self.solution[i].t)) 
        if (len(set(framenumbers)) != 1):
            warn("Number of timesteps for all ring must be the same in order to animate deformation.")
            return
        # if (nRings != 4):
        #     warn("Animation function works correctly only for 4 NPC rings at the moment.")
        #     return
        
        nframes = len(self.solution[0].t)
        
        self.nConnect = nConnect
        self.symmet = symmet
        self.xy = self.xydata()        
        self.stream = self.data_stream(self.xy)
        
        # Setup the figure and axes...
        self.axscale = 1.2 * (np.amax(randfs) + max(r))
        self.fig, self.ax = plt.subplots(figsize = (9, 10))
        plt.rcParams.update({'font.size': 20})
        
        # Then setup FuncAnimation.
        self.ani = animation.FuncAnimation(self.fig, self.update, interval=(5000/nframes), 
                                          init_func=self.setup_plot, blit=True)
        #HTML(self.ani.to_html5_video())
        #self.ani.save("Damping0.mp4", dpi = 250)
        plt.show()

    def xydata(self):
        xy = []
        for ring in range(self.nRings):
            xy.append(Pos2D(self.solution[ring])[:, np.append(np.arange(self.symmet), 0)])
        return xy
            
    def setup_plot(self):
        """Initial drawing of the scatter plot."""

        self.lines = []
        for i in range(int(self.nRings*2 + self.symmet*self.nRings*self.nConnect)):   #TODO code for 4 rings only!
            if (i <= 1): # 0, 1: lower rings
                self.lobj = self.ax.plot([], [], marker = "o", color = "gray", markerfacecolor = "black", linestyle = "", markersize = 15) 
            elif (i > 1 and i <=3): #2, 3 upper rings
                self.lobj = self.ax.plot([], [], marker = "o", color = "gray", markerfacecolor = "gray", linestyle = "", markersize = 15) 
            elif (i > 3 and i <= 7): #4, 5, 6, 7: 4 rings to anchor
                self.lobj = self.ax.plot([], [], marker = "", color = "orange", linestyle = "-", zorder = 0) # anchor
            else: # 8 - ? #all circumferential springs
                self.lobj = self.ax.plot([], [], marker = "", color = "blue", linestyle = "-")

            self.lines.append(self.lobj)
        
        # self.lines = []
        # for i in range(int(self.nRings*2 + self.symmet*self.nRings*self.nConnect)):   #TODO code for 4 rings only!
        #     if (i <= 0): # 0 lower ring
        #         self.lobj = self.ax.plot([], [], marker = "o", color = "gray", markerfacecolor = "black", linestyle = "", markersize = 15) 
        #     elif (i > 0 and i <= 1): #1: ring to anchor
        #         self.lobj = self.ax.plot([], [], marker = "", color = "orange", linestyle = "-", zorder = 0) # anchor
        #     else: # 8 - ? #all circumferential springs
        #         self.lobj = self.ax.plot([], [], marker = "", color = "blue", linestyle = "-")

        #     self.lines.append(self.lobj)
        
        
        self.ax.axis("scaled")
        self.ax.set(xlabel = "x (nm)", ylabel = "y (nm)")  
        self.ax.axis([-self.axscale, self.axscale, -self.axscale, self.axscale])
        
        return [self.lines[i][0] for i in range(int(self.nRings*2 + self.symmet*self.nRings*self.nConnect))]
        
    def data_stream(self, pos):
        x = np.zeros((self.symmet+1, self.nRings))
        y = np.zeros((self.symmet+1, self.nRings))
        while True: 
            for i in range(len(self.xy[0])):
                for ring in range(self.nRings):
                    x[:, ring] = self.xy[ring][i][:, 0]
                    y[:, ring] = self.xy[ring][i][:, 1]
                yield x, y
        
    def update(self, i):
        """Update the plot."""

        x, y = next(self.stream)
        
        xa = np.zeros((2*self.symmet, self.nRings))
        ya = np.zeros((2*self.symmet, self.nRings))     
        
        for ring in range(self.nRings):
            for i in range(1, 2*self.symmet, 2): 
                xa[i, ring] = x[int((i-1)/2), ring]
                ya[i, ring] = y[int((i-1)/2), ring]
        
        xlist = list(x.T) + list(xa.T)
        ylist = list(y.T) + list(ya.T)
                
        for lnum, self.line in enumerate(self.lines):
            if lnum >= len(xlist):
                break
            self.line[0].set_data(xlist[lnum], ylist[lnum]) 

        # TODO code  works only for 4 rings!
        count = len(xlist)
        for lnum in range(self.nRings):
            for ni in range(1, self.nConnect+1): # neighbours to connect to
                for i in range(self.symmet): # node to connect from 
                    self.lines[count][0].set_data((xlist[lnum][i], xlist[lnum][(i+ni)%self.symmet]), (ylist[lnum][i], ylist[lnum][(i+ni)%self.symmet])) 
                    count += 1
        return [self.lines[i][0] for i in range(int(self.nRings*2 + self.symmet*self.nRings*self.nConnect))]

if __name__ == '__main__':
#    a = AnimatedScatter(solution, nRings, nConnect, symmet, r, randfs)

    plt.show()
    
    
#XYoverTime(solution)
#Plotforces(fcoords, initcoords)
#Plot2D(solution, anchorsprings=False, radialsprings=False, trajectory=True, legend = False)    
#Plot3D(solution, z, symmet, viewFrame = -1)#, colour = ["black", "black", "gray", "gray"])

#zcolour = list(np.ones(len(NPCoffset[:,2])))

#zcolour = 
#NPCoffsetlist = list(NPCoffset[:,2])
#zcolour = NPCoffsetlist[0.1 if i<20 else 0.9 for i in NPCoffsetlist]
#zcolour = ["0.2" if i<20 else "0.7" for i in NPCoffsetlist]

# fig, ax = plt.subplots(1, 1, figsize = (32, 18))
# #ax.set_title("mag " + str(mag))
# ax.scatter(NPCoffset[:,0], NPCoffset[:,1], c = zcolour)#, cmap = 'copper')
# ax.set(xlabel = "nm", ylabel = "nm")
# ax.axis("scaled")