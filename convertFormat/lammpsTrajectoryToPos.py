# Author @ Nydia R. Varela-Rosales
# Version v1 2023
# Description: convert lammps trajectory to pos file
# Requires: numpy, sys, os

import numpy as np
import sys, os

inputFile  = sys.argv[1]
output     = inputFile.replace(".lammpstrj",".pos")


def lammps_to_posFile(input_file, output_directory):
    # Create the output directory if it doesn't exist
    #os.makedirs(output_directory, exist_ok=True)
    with open(input_file, 'r') as lammps_file:
       lines   = lammps_file.readlines()
    fileSize   = len(lines)
    fileIn     = open(inputFile, 'r')
    frame      = 0
    atom_count = 0
    outputFile = open(output,"w")
    #fileIn.readlines()
    for line in range(0,fileSize):
        lineLocal =  fileIn.readline()
        if 'ITEM: TIMESTEP' in lineLocal:
            frame     += 1
            atom_count = 0
            outputFile.write(f'#[data] Frame \n')
            outputFile.write(f'{frame}\n')
            outputFile.write(f'#[done]\n')
        elif 'ITEM: BOX BOUNDS' in lineLocal:
            box1      = fileIn.readline().split()
            box2      = fileIn.readline().split()
            box3      = fileIn.readline().split()
            box1_low  = float(box1[0])
            box1_high = float(box1[1])
            box2_low  = float(box2[0])
            box2_high = float(box2[1])
            box3_low  = float(box3[0])
            box3_high = float(box3[1])
            if (box1_low > 0):
                bx = box1_low - box1_high
            if (box1_low < 0):
                bx = -box1_low + box1_high
            
            by = -box2_low + box2_high
            bz = -box3_low + box3_high
            outputFile.write(f'box {bx} {by} {bz}\n')
            outputFile.write(f'def B "sphere 1 ffff0000"\n')
            outputFile.write(f'def A "sphere 1 ff3366ff"\n')
        elif 'ITEM: ATOMS' in lineLocal:
            while ('ITEM: TIMESTEP' not in fileIn.readline()):
                atom_data = fileIn.readline().split()
                atom_id = atom_data[1]
                if (atom_id == "1"):
                    ID = "B"
                else:
                    ID = "A"
                x = round(float(atom_data[2]) - bx/2,4)
                y = round(float(atom_data[3]) - by/2,4)
                z = round(float(atom_data[4]) - bz/2,4)
                outputFile.write(f'{ID} {x} {y} {z}\n')
            outputFile.write(f'eof\n')
        outputFile.flush()

# Usage example
lammps_to_posFile(inputFile, output)