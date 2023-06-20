# Author @ Nydia R. Varela-Rosales
# Version v1 2022
# Description: cut structure
# Requires: numpy
import sys
import os
import numpy as np
import matplotlib

class readFileTreatment:
    def __init__(self,fileName,jumpLines,separationFactor,indexesOfInterest,castFloat=True):
        """read files in an arbitrary fashion manner
        return: list of the wished columns
        """
        self.fileName = fileName # file to read
        self.jumpLines = jumpLines # lines to ignore
        self.separationFactor = separationFactor # separation of lines
        self.indexesOfInterest = indexesOfInterest # which index do we want to retrieve
        self.castFloat = castFloat # cast string to float
    def countLinesFile(self):
    	num_lines = sum(1 for line in open(self.fileName))
    	return num_lines
    def reader(self):
        atomPosList = np.zeros(self.countLinesFile()-self.jumpLines-1)
        list = []
        fileInput = open(self.fileName,"r")
        for header in range(self.jumpLines):
            fileInput.readline()
        #print(self.indexesOfInterest)
        for lines in range(0,self.countLinesFile()-self.jumpLines-1):
            pos = [(text) for text in fileInput.readline().strip().split(self.separationFactor)]
            #print(pos)
            if (self.castFloat == True):
                list.append([float(pos[self.indexesOfInterest[0]]),float(pos[self.indexesOfInterest[1]]),float(pos[self.indexesOfInterest[2]])])
            else:
                list.append([(pos[self.indexesOfInterest[0]]),(pos[self.indexesOfInterest[1]])])
        #print(list)
        return list


############ LOOPING OVER ALL DATA FILES ################
filename='./'
filelist = []
class using_os(object):
    def __init__(self,filename):
    	self.filename = filename
    def __call__(self,filename):
    	filename = self.filename
    	filelist=os.listdir(filename)
    	for fichier in filelist[:]:
    		if not(fichier.endswith("6.pos")):
    			filelist.remove(fichier)

    	return filelist
lis = using_os(filename)
listData = lis.__call__(filename)

def countLinesFile(fileName):
	num_lines = sum(1 for line in open(fileName))
	return num_lines



def cleanCoordinates(coords,radius,choppedFile):
    choppedCoordinates = []
    choppedFile.write("//date: Wednesday, September 28, 2022, 10:19:52 AM\n"+
    "translation	0	0	402.860174\n"+
    "box	116.295715	116.295715	116.295715\n"+
    "shape "+f'"{"sphere 0.4 ff54dcad"}"'+"\n")
    for i in range(len(coords)):
        part = coords[i]
        origin = [0,0,0]
        circleEvaluation = np.sqrt((part[0] - origin[0]) ** 2 + (part[1] - origin[1]) ** 2 + (part[2] - origin[2]) ** 2 )
        if (circleEvaluation < radius ):
            choppedCoordinates.append(part)
            choppedFile.write(str(part[0])+"\t"+str(part[1])+"\t"+str(part[2])+"\n")
            choppedFile.flush()
    #choppedFile.close()
    return choppedCoordinates


def chop(inputFileName,box):
    choppedFile = open(inputFileName.replace(".pos","")+"chopped.pos","w")
### read file and set up trajectories as feature vector
    posFileName = inputFileName #ico_2604_001.pos"
    fileRead = readFileTreatment(posFileName,7,"\t",[0,1,2],True)
    featureVectorPre = fileRead.reader()
    featureVector = cleanCoordinates(featureVectorPre,box,choppedFile)

for i in range(len(listData)):
    chop(listData[i],7)
