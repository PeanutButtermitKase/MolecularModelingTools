# Author @ Nydia R. Varela-Rosales
# Version v1 2020
# Description: constructs: substrate, superstructure, points on substrate

import numpy as np
import matplotlib.pyplot as plt
import subprocess,sys
import math
import seaborn as sns

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
        print(self.indexesOfInterest)
        for lines in range(0,self.countLinesFile()-self.jumpLines-2):
            pos = [(text) for text in fileInput.readline().strip().split(self.separationFactor)]
            if (self.castFloat == True):
                list.append([float(pos[self.indexesOfInterest[0]]),float(pos[self.indexesOfInterest[1]])])
            else:
                list.append([(pos[self.indexesOfInterest[0]]),(pos[self.indexesOfInterest[1]])])
        return list
    def readerCostum(self):
        mainList = []
        if (len(self.indexesOfInterest) == 3):
            atomPosList = np.zeros(self.countLinesFile()-self.jumpLines-2)
            #list = []
            fileInput = open(self.fileName,"r")
            for header in range(self.jumpLines+1):
                fileInput.readline()
            print(self.indexesOfInterest)
            for lines in range(0,self.countLinesFile()-self.jumpLines-2):
                pos = [(text) for text in fileInput.readline().strip().split(self.separationFactor)]
                if (self.castFloat == True):
                    mainList.append([float(pos[self.indexesOfInterest[0]]),float(pos[self.indexesOfInterest[1]]),float(pos[self.indexesOfInterest[2]])])
                else:
                    mainList.append([(pos[self.indexesOfInterest[0]]),(pos[self.indexesOfInterest[1]]),(pos[self.indexesOfInterest[2]])])
        else:
            mainList = self.reader()
        return mainList
    def makeLists2params(self):
        """
        This function adjust the depth of the potential as:
        (val)*(epsilon)/(temperature) ==reversed==>
        (val)/(epsilon)/(temperature)
        """
        lst = self.readerCostum()
        x = np.zeros(len(lst))
        y = np.zeros(len(lst))
        #score = np.zeros(len(lst))
        for i in range(len(lst)):
            x[i] = lst[i][0]
            y[i] = lst[i][1] #round((lst[i][1])/((epsilon)/lst[i][0]),1) #lst[i][1]
           
        return x,y
    def adjustScalinginDepth(self,epsilon):
        """
        This function adjust the depth of the potential as:
        (val)*(epsilon)/(temperature) ==reversed==>
        (val)/(epsilon)/(temperature)
        """
        lst = self.readerCostum()
        x = np.zeros(len(lst))
        y = np.zeros(len(lst))
        score = np.zeros(len(lst))
        for i in range(len(lst)):
            x[i] = lst[i][0]
            y[i] = round((lst[i][1])/((epsilon)/lst[i][0]),1)

            score[i] = lst[i][2]

        return x,y,score

class gnuPlot(readFileTreatment):
    def __init__(self,fileName,jumpLines,separationFactor,indexesOfInterest,castFloat,interestPoint1,interestPoint2):
        super().__init__(fileName,jumpLines,separationFactor,indexesOfInterest,castFloat)
        """
        self.interestPoint1 and self.interestPoint2 keep the [name,number]
        example = self.interestPoint1 = ["amorphous",1], self.interestPoint2 = ["crystalline",0]
        """
        self.interestPoint1 = interestPoint1
        self.interestPoint2 = interestPoint2
        self.figureMatplotlibExt = ".svg"
        self.gnuplotFileName = "gnuplot_"+self.fileName.replace(".pos",".txt")
        self._externalGNUfile = "plotGNUcustomPhaseDiagram.sh"
    def makeFileforGNUposttreatment(self,epsilon):
        x,y,z = self.adjustScalinginDepth(epsilon)
        openFile = open(self.gnuplotFileName,"w")
        for n in range(len(x)):
            openFile.write(str(x[n])+" "+str(y[n])+" "+str(int(z[n]))+"\n")
        openFile.close()

    def plotPhasePlotLib(self,epsilon):
        """
        Gather information from readFileTreatment, using special case
        self.readerCostum()
        """
        x,y,z = self.adjustScalinginDepth(epsilon)
        plt.xlabel("Temperature")
        plt.ylabel("Potential depth [e/T]")
        plt.scatter(x,y,s=75, c=list(z), alpha=.5)
        #plt.legend()
        plt.savefig("phase_diagram_lowTemperature"+self.figureMatplotlibExt)
        #plt.show()
        return x,y,z
    def plotPhaseGNUplot(self,epsilon):
        """
        this part is necessary to take from a file
        since this is executed independently
        requires: 3 rows in the main file
        """
        self.makeFileforGNUposttreatment(epsilon) #make file
        # open file and use external gnuplot file call: "plotGNUcustomPhaseDiagram.sh"
        subprocess.run(["gnuplot","-e","filename='"+self.gnuplotFileName+"'",self._externalGNUfile])

class plotPotentialOverSubstrate(readFileTreatment):
    def __init__(self,fileName,jumpLines,separationFactor,indexesOfInterest,castFloat,period,lambdaVar,depth):
        super().__init__(fileName,jumpLines,separationFactor,indexesOfInterest,castFloat)
        self.period = period
        self.lambdaVar = lambdaVar
        self.depth = depth

    def customPeriodicPotentialHOOMD(self,rx,ry,eps,temp):
        """
        k = 2*pi/lambda
        angle60 = 2*pi/6
        xi = cos(angle*i)
        yi = sin(angle*i)
        ri = xi + yi
        V(r) = cos(k*ri)
        """
        scalingFactor = 2*math.pi/self.period#2*math.pi/self.lambdaVar#1/(2*math.pi*self.period)
        #lambdaParam = self.lambdaVar
        k = scalingFactor#(2*math.pi/self.periodicity)
        angle60 = (2*math.pi)/6
        Vtotal = []
        Vi = 0
        #phi = self.orderingParam
        #### rememebr the depth is scaled by temp
        ####round((order*(lepsilon/temperature)),1)
        potentialDepth = 1#round((self.depth*(eps/temp)),1)#self.depth
        for i in range(3):
            xi = np.cos(angle60*i)*rx
            yi = np.sin(angle60*i)*ry
            ri = xi + yi
            Vi += (np.cos(k*ri))
            #Vtotal.append(Vi)
        return potentialDepth*Vi
    def plotPhasePlotLib(self,eps,temp):
        """
        Gather information from readFileTreatment, using special case
        self.readerCostum()
        """
        self.figureMatplotlibExt = ".svg"
        maxSteps = 200#0
        x,y = self.makeLists2params()
        max_X = 6#0#np.max(x)
        min_X = -max_X #np.min(x)
        max_Y = max_X#np.max(y)
        min_Y = -max_X#np.min(y)
        ##### rewrite x and y #####
        newPos = []
        #newY = []
        for i in range(len(x)):
            #if ((x[i] >= -max_X) and (y[i] >= -max_X)): continue
            if (((x[i] <= max_X) and (y[i] <= max_X)) and ((x[i] >= -max_X) and (y[i] >= -max_X))): #continue
                newPos.append([x[i],y[i]])
            else: continue
        newPos = np.asarray(newPos)
        ### define substrate range
        # max_X = np.max(x)
        # min_X = np.min(x)
        # max_Y = np.max(y)
        # min_Y = np.min(y)
        xSubstrate = np.linspace(min_X,max_X,maxSteps)
        ySubstrate = np.linspace(min_Y,max_Y,maxSteps)
        ### compute substrate
        zSubstrate = []
        for i in range(len(xSubstrate)):
            zSubstrate.append(self.customPeriodicPotentialHOOMD(xSubstrate[i],ySubstrate[i],eps,temp))
        ######
        X, Y = np.meshgrid(xSubstrate, ySubstrate)
        zs = np.array(self.customPeriodicPotentialHOOMD(np.ravel(X), np.ravel(Y),eps,temp))
        plt.xlabel("x-domain")
        plt.ylabel("y-domain")
        import matplotlib.tri as tri
        #plt.scatter(xSubstrate,ySubstrate,c=zs,s=75, alpha=.5)
        paletteEnergy = sns.cubehelix_palette(start=.5, rot=-.75, n_colors=len(zs)).as_hex()
        # plt.scatter(X,Y,c=zs,cmap='binary') # plot substrate
        def generateNN(positions):
            # n = []
            # listNN = []
            color = "ffff0000"
            filePos = open("particles-substrate_"+"dist_"+str(max_X)+str(self.period)+".pos","w")
            for i in range(len(positions)):
                pos = positions
                filePos.write("sphere 1 "+color+" "+str(pos[i][0])+" "+str(pos[i][1])+" "+str(0)+"\n")
            filePos.close()
            # Create the Triangulation; no triangles so Delaunay triangulation created.
            # triang = tri.Triangulation(newPos[:,0], newPos[:,1])
            # fig1, ax1 = plt.figure()
            # plt.scatter(X,Y,c=zs,cmap='binary') # plot substrate
            # ax1.set_aspect('equal')
            # ax1.triplot(triang, 'bo-', lw=1)
            # ax1.set_title('triplot of Delaunay triangulation')
            #plt.show()
            #plt.savefig("connections.png")

        neigh = generateNN(newPos)
        triang = tri.Triangulation(newPos[:,0], newPos[:,1])
        fig1, ax1 = plt.subplots()
        plt.scatter(X,Y,c=zs,cmap='binary') # plot substrate
        plt.scatter(newPos[:,0],newPos[:,1],c="red",edgecolor="r",linewidths=7) # plot particles
        ax1.set_aspect('equal')
        ax1.triplot(triang, 'bo-', lw=1)
        ax1.set_title('triplot of Delaunay triangulation')

        plt.savefig("particles-substrate_"+"dist_"+str(max_X)+str(self.period)+self.figureMatplotlibExt)
        plt.show()

        return x,y

class generateSuperstructure:
    def __init__(self,posFile):
        self.posFile = posFile
        self.outFileName = self.posFile.replace(".pos","patterns.pos")
    def patternsInjavisFinder(self):
        subprocess.run(["java","-jar","/nishome/students/nydia/injavis.jar",self.posFile,"-A","RDF","Show",
        "-A","RDF","Detect Peak","-A","Network","Find Triangles",
        "-A","Edit","Color: Green","-A","Network","Find Squares",
        "-A","Edit","Color: Blue","-A","Network","Find Pentagons","-o",self.outFileName])
    def findTilings(self):
        """
        generates tilings from injavis
        using file with inputName
        """
        ##### genereta initial patterns (square,triangel,rhombs)
        self.patternsInjavisFinder()
        ##### genereta superstructure (shield tilings)
        sphereSize = 0.3 # size of the "particles" in new File
        #fileWithPos = self.outFileName#inputName # file containing spheres info
        self.newName = self.outFileName#"patterns"+inputName
        outName = self.newName.replace(".pos","centers.pos")
        outName1 = self.newName.replace(".pos","tilings.pos")
        # open new file
        newFile = open(outName,"w")
        # open file with particles as the centers of the triangles,pentagons,squares
        newPosFile = open(outName1,"w")
        ### write spheres information ##
        fileWithPos = open(self.outFileName,"r")
        for l in fileWithPos.readlines():
            if (l.startswith("poly")):
                continue
            else:
                 newFile.write(l) # wirte in the file of the centers
                 #newPosFile.write(l) # write also in the posFile
        for line in open(self.newName).readlines():
            if any( line.startswith(x) for x in ['box', 'eof'] ):
                print(line)
                newPosFile.write(line)
                #newFile.write(line)
            if line.startswith('poly'):
                lsp = line.split()
                color = lsp[1]
                vertices = int(lsp[2])
                coords   = [ float(x) for x in lsp[3:] ]
                center   = [ sum([ coords[id+iv*3] for iv in range(vertices) ])/vertices for id in range(3) ]
                #tmpString = ('poly3d {} '.format(vertices), end='')
                #print('poly3d {} '.format(vertices), end='')
                newFile.write('poly3d '+str(vertices) + " ")
                for iv in range(vertices):
                    for id in range(3):
                        #print('{:f} '.format( coords[id+3*iv]-center[id] ),end='')
                        newFile.write(str( coords[id+3*iv]-center[id] ) + " ")
                #print(color,end='')
                newFile.write(color+" ")
                #print(' {} {} {}   1 0 0 0'.format(*center))
                newFile.write(str(center[0])+" "+str(center[1])+" "+str(str(center[2]))+' 1 0 0 0'+"\n")
                """ normal posfile constructed from the centers of the particles """
                newPosFile.write("sphere "+str(sphereSize)+" "+color+" "+str(center[0])+" "+str(center[1])+" "+str(center[2])+"\n")
        newFile.close()
        newPosFile.close()
        return [outName,outName1]


#### init instance for phase diagram computation ####
# inputFile ="phaseDiagramSquareLatticeConsidered.txt"# "phaseDiagram.txt"
# lepsilon = 1.8
# instFile = gnuPlot(inputFile,1," ",[0,1,2],True,["amorphous",1],["crystalline",0])
# plotLines = instFile.plotPhaseGNUplot(lepsilon)
###### potential set-up #####
#inputFilePot_Sub = "./lastFrame/lastFrame_nvt_temp_0.31relaxTime10000potDepth_2.4period_2_importantRegionSquare.pos"
inputFilePatterns = sys.argv[1]#"lastFrame_nvt_temp_0.35relaxTime100simTime_1000000initLattice_simBigQuasitauValues_0.1ord_0.0_relaxation.pos"#"./lastFrame/lastFrame_nvt_temp_0.16relaxTime10000potDepth_0.1period_2_importantRegion.pos"
#inputFileSuperStruct = "./lastFrame/lastFrame_nvt_temp_0.16relaxTime10000potDepth_0.1period_2_importantRegionpatternstilings.pos"
#periodParam = eval(sys.argv[2])#4#12#16
lambdaParam = 1#(2*math.pi)/1
depthParam = 0.4#0.05#2.4
lepsilon = 1#1.8
temperature = 1#0.16#0.31
### 0.1 ofr patterns temp 0.16
### 2.4 for squares temp 0.31
#for i in range(1,5):
instCompariosnsubstrate = plotPotentialOverSubstrate(inputFilePatterns,4,"	",[0,1],True,depthParam,lambdaParam,depthParam)
instCompariosnsubstrate.plotPhasePlotLib(lepsilon,temperature)
##### genereta superSctruture #####
#instSuperstruct = generateSuperstructure(inputFilePatterns)
#instSuperstruct.findTilings()
#### super structure centers and background potential ####
#instCompariosnsubstrate = plotPotentialOverSubstrate(inputFileSuperStruct,1," ",[3,4],True,periodParam,lambdaParam,depthParam)
#instCompariosnsubstrate.plotPhasePlotLib(lepsilon,temperature)
