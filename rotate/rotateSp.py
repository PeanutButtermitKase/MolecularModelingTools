# Author @ Nydia R. Varela-Rosales
# Version v1 2022
# Description : rotates structure with reference to another structure
# Requires    : numpy

import pickle
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import sys,os,math,subprocess
from random import seed
from random import random
import statistics

fileIn      = sys.argv[1]
intputFiles = sys.argv[2] 
targetAngle = np.arange(1,30)
rotate      = False
# example structure "lastFrame_nvt_temp_0.2relaxTime100simTime_10000000initLattice_specificFileOutInjavistauValues_0.1ord_0.53333333.pos"

class readFileTreatment:
    def __init__(self,fileName,jumpLines,separationFactor,indexesOfInterest,castFloat=True):
        self.fileName = fileName
        self.jumpLines = jumpLines
        self.separationFactor = separationFactor
        self.indexesOfInterest = indexesOfInterest
        self.castFloat = castFloat
        self.thirdIndex = None
    def countLinesFile(self):
    	num_lines = sum(1 for line in open(self.fileName))
    	return num_lines
    def reader(self):
        atomPosList = np.zeros(self.countLinesFile()-self.jumpLines-2)
        list = []
        fileInput = open(self.fileName,"r")
        for header in range(self.jumpLines+1):
            fileInput.readline()
        for lines in range(0,self.countLinesFile()-self.jumpLines-2):
            pos = [(text) for text in fileInput.readline().strip().split(self.separationFactor)]
            if (self.castFloat == True):
                list.append([float(pos[self.indexesOfInterest[0]]),float(pos[self.indexesOfInterest[1]])])
            else:
                list.append([(pos[self.indexesOfInterest[0]]),(pos[self.indexesOfInterest[1]])])
        return list

def rotate(p, origin=(0, 0), degrees=0):
    angle = np.deg2rad(degrees)
    R = np.array([[np.cos(angle), -np.sin(angle)],
                  [np.sin(angle),  np.cos(angle)]])
    o = np.atleast_2d(origin)
    p = np.atleast_2d(p)
    return np.squeeze((R @ (p.T-o.T) + o.T).T)

def rotateInitialSnapshot(fileName,angleDegrees):
    fileSpirals  = readFileTreatment(intputFiles,6,"\t",[0,1],True) # in file to put hexagona lattice on
    instMethod   = readFileTreatment(fileName,6,"\t",[0,1],True)    #readFileTreatment(fileName,3," ",[3,4],True)  # in hexagonal lattice
    instMethod1  = instMethod.reader()
    fileSpirals1 = fileSpirals.reader()
    fileOutSin      = "rotated_"+str(angleDegrees)+".pos"
    fileOutSinMix   = intputFiles.replace(".pos","")+"rotated_"+str(angleDegrees)+"_mixed.pos"
    rotatedStructure      = open(fileOutSin,"w")
    rotatedStructureMixed = open(fileOutSinMix,"w")
    for i in range(len(instMethod1)):
        rot = rotate(instMethod1[i], origin=(0, 0), degrees=angleDegrees)
        rotatedStructure.write("sphere 1 ffff0000 "+str(rot[0])+" "+str(rot[1])+" 0"+"\n")
        rotatedStructureMixed.write("sphere 0.5 ff0960EA "+str(rot[0])+" "+str(rot[1])+" 0"+"\n")
    for i in range(len(fileSpirals1)):
        rot = rotate(fileSpirals1[i], origin=(0, 0), degrees=0)
        rotatedStructureMixed.write("sphere 0.25 fff21414 "+str(rot[0])+" "+str(rot[1])+" 0"+"\n")
    rotatedStructure.close()
    rotatedStructureMixed.close()
    return [fileOutSin,fileOutSinMix]

def getDiffraction(inF,out):
    pathInjavis = "/nydia/"
    subprocess.run(["java","-jar",pathInjavis+"injavis.jar",inF,"-A","Diffraction","Show",
    "-s","Diffraction","zoom","0.5",
    "-s","Diffraction","size","1024","-M","Diffraction",out,"-q"])

getDiffraction(intputFiles,intputFiles.replace(".pos","diffraction.png"))

for i in targetAngle:
    single,mixed = rotateInitialSnapshot(fileIn,round(i,2))
    getDiffraction(mixed,mixed.replace(".pos","diffraction.png"))

if rotate:
    angleList = np.arange(0.6,2.0,0.1)
    #print(angleList)
    for i in range(len(angleList)):
        rotateInitialSnapshot(fileIn,round(angleList[i],2))
