# Author @ Nydia R. Varela-Rosales
# Version v1 2020
# Description: generate tiles
import subprocess,sys
class generateSuperstructure:
    def __init__(self,posFile):
        self.posFile = posFile
        self.outFileName = self.posFile.replace(".pos","patterns.pos")
    def patternsInjavisFinder(self):
        pathToInjavis = "/nishome/students/nydia/injavis.jar"
        subprocess.run(["java","-jar",pathToInjavis,self.posFile,"-A","RDF","Show",
        "-A","RDF","Detect Peak","-A","Network","Find Triangles",
        "-A","Edit","Color: Green","-A","Network","Find Squares",
        "-A","Edit","Color: Blue","-A","Network","Find Pentagons","-o",self.outFileName])
   