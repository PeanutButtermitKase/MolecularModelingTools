# Author @ Nydia R. Varela-Rosales
# Version v1 2021
# Description : analyze order
# Requires    : scipy

import sys
import os
import numpy as np
import matplotlib
import scipy.special
import cmath,math
mode    = sys.argv[1] # file coming the xyz of kind -->[index|x|y|z|radius]
NNsFile = sys.argv[2] # file coming from the output of setVoronoi

"""if the mode is "g" means look for all files in the directory with the termination """
if (mode=="g"):
    fileName = './'
    """else we are just looking for a specific file """
else:
    fileName = mode
    
filelist = []
class get_path(object):
    def __init__(self,filename):
    	self.filename = filename
    def __call__(self,filename):
    	filename = self.filename
    	filelist=os.listdir(filename)
    	for f in filelist[:]:
    		if not(f.endswith(".dat")):
    			filelist.remove(f)

    	return filelist
lis = get_path(filename)
list = lis.__call__(filename)

def countLinesFile(fileName):
	num_lines = sum(1 for line in open(fileName))
	return num_lines

############## GATHERING DATA #######################
def readerNNsPomelo(fileName,NNsFile,newFile):
    neighList = np.arange(1,20)
    colorList = ["0x00748A","0x3FFF00","0x3393ff","0x009767","0x7DFF00","0x00F30B","0xF0FFFF","0x00DB23","0x00AE50","0xBC8F8F","0x4169E1","0x8B4513","0xFA8072","0xF4A460","0x2E8B57","0xFFF5EE","0xA0522D","0xC0C0C0","0x1E90FF","0xB22222"]#,"0xFFFAF0","0x228B22","0xFF00FF","0xDCDCDC","0xF8F8FF","0xFFD700","0xDAA520","0x808080","0x008000","0xADFF2F","0xF0FFF0","0xFF69B4","0xCD5C5C","0x4B0082","0x87CEEB","0xFFFF33"]
    f = open(fileName,"r",encoding='utf-8-sig')    ## Open xyz format from pomelo
    f_nns = open(NNsFile,"r",encoding='utf-8-sig') ## Open setVoronoiNeighbors.dat format from pomelo
    f_nns1 = open(NNsFile,"r",encoding='utf-8-sig')
    fOut = open(newFile,"w")
    file_q4_q6 = open(newFile+"q4","w") ## This file will contain the information to plot q6 vs q4
    defaultColor = "0x00748A" # general coloring
    """First header expceted of the form:
    boundary_condition = periodic_cuboidal, infile = hs-16384_0.50-coord.dat, num_sph = 18, boxsx =100, boxsy = 51.58374 , boxsz = 51.583743
    """
    f.readline() #read first header
    f.readline() # read second header ---> # x y z r
    '''First header from "setVoronoiNeighbors.dat" expceted of the form:
    #1_particle label #2_number of neighbours'''
    f_nns.readline() ## read header
    radii = np.zeros(countLinesFile(fileName)) ## zeros vector
    x = np.zeros(countLinesFile(fileName)) ## zeros vector
    y = np.zeros(countLinesFile(fileName)) ## zeros vector
    z = np.zeros(countLinesFile(fileName)) ## zeros vector
    neighbors = np.zeros(countLinesFile(NNsFile))
    for i in range(0,countLinesFile(fileName)-3):
        numbers = [(text) for text in f.readline().strip().split(" ")]
        #radii[i] = numbers[3]
        x[i] = numbers[0] #first column
        y[i] = numbers[1] # second column
        z[i] = numbers[2]
        #neighbors[i] = numbers[4]
    for i in range(0,countLinesFile(NNsFile)-1):
        numbers = [(text) for text in f_nns.readline().rstrip("()").split(" ")]
        neighbors[i] = numbers[1]
    connections = []
    f_nns1.readline() # read Header
    for i in range(0,countLinesFile(NNsFile)-1):
        a = f_nns1.readline().split("(")
        try:
            b = a[1].replace(")","").strip().split(" ")
            connections.append(b)
        except IndexError:
            pass
        continue
    """ fix values that were not assigned neigbors"""
    listNN = []
    for i in range(0,len(connections)):
        listNN.append([])
        for n in connections[i]:
            if (n==""):
                t = 0
            else:
                t = n
            listNN[i].append(int(t))
    """Write posFile base on the proper nubmer of neihgbor identified"""
    for i in range(0,len(x)):
        for j in neighList:
            print(j,neighbors[i])
            if (int(neighbors[i])==j):
                #print(j,neighbors[i])
                fOut.write("sphere "+str(1)+" "+str(colorList[j])+" "+str(x[i])+" "+str(y[i])+" "+str(z[i])+" "+"\n")
            else: continue
    """Here we are going to calculate the q6 values to plot the histogram
        what is needed: [coordinates] [connections]
        """
  
    lmax = 6 # Max value to expand the spherical harmonics
    qlm_sum = [0+0j for t in range(2 * lmax + 1)]
    qlm_count = 0
    for i in range(0,len(x)-3):
        for j in listNN[i]:
            if i==j:
            	continue
            dx = x[i]-x[j]
            dy = y[i]-y[j]
            dz = z[i]-z[j]
            distance = (dx*dx) + (dy*dy) + (dz*dz)
            theta = np.arctan2(dy, dx)
            phi   = np.arccos( dz/np.sqrt(distance) )
            qlm = scipy.special.sph_harm(range(-lmax ,lmax +1),lmax , theta, phi)
            print(qlm)
   

newFile = mode.replace(".xyz",".pos")
readerNNsPomelo(mode,NNsFile,newFile)
