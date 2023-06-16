##### OPP SIMULATIONS WITH LOCAL DESCRIPTORS
##### written 2021 by Nydia R. Varela-Rosales
import sys
import pythia, struct
import seaborn as sns
import numpy as np
import math, random, shutil
import scipy.optimize as opt
import math
# parallel libraries needed
import freud
from freud import * #data, order, box
import freud.box
import multiprocessing
from multiprocessing import Pool
import subprocess
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R
from ase import Atoms,Atom
import pyscal.core as pc
import pyscal.crystal_structures as pcs
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import normalize
from sklearn.svm import SVC

### structural descriptors ###
def mID(n, d):
    M = np.zeros( (n,) * d )
    M[ tuple([np.arange(n)] * d) ] = 1
    return np.asarray(M)
# inertia tensor descriptor

def inertiaTensorDescriptor(positions):
    inertiaDescriptor = []
    n = []
    sc = 0
    for i in range(len(positions)):
        for j in range(len(positions)):
            if (i == j): continue
            a = np.asarray(positions[i])
            b = np.asarray(positions[j])
            diff = a-b
            diffR = b-a
            d = np.linalg.norm(a-b)
            if (math.sqrt(d)<=2.5):
                n.append([i,j])
                sc += np.dot(diff,diff)
            cc = np.kron(diff,diffR)
            s = sc*mID(1, 3) - cc
        inertiaDescriptor.append(s)
    return inertiaDescriptor

def euler_to_quaternion(yaw,pitch,roll):
    qx = np.sin(roll/2) * np.cos(pitch/2) * np.cos(yaw/2) - np.cos(roll/2) * np.sin(pitch/2) * np.sin(yaw/2)
    qy = np.cos(roll/2) * np.sin(pitch/2) * np.cos(yaw/2) + np.sin(roll/2) * np.cos(pitch/2) * np.sin(yaw/2)
    qz = np.cos(roll/2) * np.cos(pitch/2) * np.sin(yaw/2) - np.sin(roll/2) * np.sin(pitch/2) * np.cos(yaw/2)
    qw = np.cos(roll/2) * np.cos(pitch/2) * np.cos(yaw/2) + np.sin(roll/2) * np.sin(pitch/2) * np.sin(yaw/2)
    return np.quaternion(qw, qx, qy, qz)

def quaternion_to_euler(w,x,y,z):
    t0 = +2.0 * (w * x + y * z)
    t1 = +1.0 - 2.0 * (x * x + y * y)
    roll = math.atan2(t0, t1)
    t2 = +2.0 * (w * y - z * x)
    t2 = +1.0 if t2 > +1.0 else t2
    t2 = -1.0 if t2 < -1.0 else t2
    pitch = math.asin(t2)
    t3 = +2.0 * (w * z + x * y)
    t4 = +1.0 - 2.0 * (y * y + z * z)
    yaw = math.atan2(t3, t4)
    return [yaw, pitch, roll]

def opt(stepsSearch):
    """
    return: initial list for optimization
    """
    ListX = []
    ListY = []
    ListZ = []
    ListXY = []
    ListXZ = []
    ListZY = []
    allDirections = []
    for i in range(stepsSearch):
    	c = math.cos(i)
    	s = math.sin(i)
    	ListX.append([c, 0, s, 0])
    	ListY.append([c, 0, 0, s])
    	ListZ.append([c, s, 0, 0])
    	ListXY.append([c, 0, s, s])
    	ListXZ.append([c, s, 0, s])
    	ListZY.append([c, s, 0, s])
    	allDirections.append([c, 0, s, 0])
    	allDirections.append([c, 0, 0, s])
    	allDirections.append([c, s, 0, 0])
    	allDirections.append([c, 0, s, s])
    	allDirections.append([c, s, 0, s])
    	allDirections.append([c, s, 0, s])
    return allDirections

def multiply(quatIn, quatOut):
    """
    return: modified orientation
    """
    br = quatOut[0] * quatIn[0] - quatOut[1] * quatIn[1] - quatOut[2] * quatIn[2] - quatOut[3] * quatIn[3];
    bi = quatOut[0] * quatIn[1] + quatOut[1] * quatIn[0] + quatOut[2] * quatIn[3] - quatOut[3] * quatIn[2];
    bj = quatOut[0] * quatIn[2] - quatOut[1] * quatIn[3] + quatOut[2] * quatIn[0] + quatOut[3] * quatIn[1];
    quatIn[3] = quatOut[0] * quatIn[3] + quatOut[1] * quatIn[2] - quatOut[2] * quatIn[1] + quatOut[3] * quatIn[0];
    quatIn[0] = br
    quatIn[1] = bi
    quatIn[2] = bj
    return quatIn;

def generateNN(positions):
    n = []
    for i in range(len(positions)):
        for j in range(len(positions)):
            if (i == j): continue
            a = np.asarray(positions[i])
            b = np.asarray(positions[j])
            diff = a-b
            diffR = b-a
            d = np.linalg.norm(a-b)
            if (math.sqrt(d)<=2.5):
                n.append([i,j])
    return n
########### ROTATE POSITIONS ##################
def rotate(ar,rot):
    rotatedPos = []
    for i in range(len(ar)):
        rotatedPos.append(np.dot(rot,ar[i]))
    return np.asarray(rotatedPos)
####### basic list search for rotation ########
def rotateEnvironment(positions,quaternions):
    """
    note: this function will be call by inertia rotation
    """
    outRot = []
    for q in range(len(quaternions)):
        # loop over mini-clusters
        for i in range(len(positions)):
            rotated = rotate(positions[i],R.from_quat(quaternions[q]).as_matrix())
            outRot.append(rotated) # rotatate all positions by q
    return outRot # return the rotated positions
############## STEINHARDT ORDER PARAMETERS ###########################
# test order parameters
reported_Q6 = 0.57452416
def testQlinFCC():
    # test for FCC
    fcc_system = freud.data.UnitCell.fcc().generate_system(4, scale=2)
    r_max = 1.5
    a = []
    b = []
    for l in range (1,13):
        ql = freud.order.Steinhardt(l)
        a.append(l)
        ql_fcc = ql.compute(fcc_system, neighbors={'num_neighbors': 12}).particle_order
        b.append(ql_fcc)
    means = [np.mean([el for el in sublist if el > 0] or 0) for sublist in b]
    print(means)
    print("Obtained ",means)
    print("Obtained", sum(ql_fcc)/len(ql_fcc))
    print("Reported in literature",reported_Q6)
    print("Difference: ",(sum(ql_fcc)/len(ql_fcc))-reported_Q6)

def testRotationsInFCC(sphSteinhardtHelpers,pyscalHelpers):
    # test for FCC
    fcc_system = freud.data.UnitCell.fcc().generate_system(4, scale=2)
    pOutBox = fcc_system[0]
    mainSizeB = pOutBox.Lx
    pOut = fcc_system[1]
    nlistCustom = generateNN(pOut)
    rotationsGlobal = opt(3)
    descriptorsList = []
    if (sphSteinhardtHelpers == True):
        intsPythia = pythia.spherical_harmonics.steinhardt_q(pOutBox,pOut, neighbors=12, lmax=6, rmax_guess=2.0)
        instPythiaGSteinhardt = pythia.spherical_harmonics.system_average(pOutBox,pOut, neigh_min=4, neigh_max=12, lmax=4, negative_m=True, reference_frame='neighborhood', orientations=None, rmax_guess=1.0, noise_samples=0, noise_magnitude=0, nlist=None)
        instPythiaGSteinhardtRotated = pythia.spherical_harmonics.system_average(pOutBox,pOut, neigh_min=4, neigh_max=12, lmax=12, negative_m=True, reference_frame='neighborhood', orientations=rotationsGlobal[0], rmax_guess=1.0, noise_samples=0, noise_magnitude=0, nlist=None)
        intsPythiaGlobalSt = pythia.spherical_harmonics.abs_neighbor_average(pOutBox,pOut, neigh_min=4, neigh_max=4, lmax=12, negative_m=True, reference_frame='neighborhood', orientations=None, rmax_guess=1.0, noise_samples=0, noise_magnitude=0, nlist=None)
        instBispectrum = pythia.spherical_harmonics.bispectrum(pOutBox,pOut,neighbors = 12, lmax=6, rmax_guess=2.0)
        descriptorsList = [[intsPythia,"pythia.spherical_harmonics.steinhardt_q"],
                        [instPythiaGSteinhardt,"pythia.spherical_harmonics.system_average"],
                        [instPythiaGSteinhardtRotated,"pythia.spherical_harmonics.system_average rotated"],
                        [intsPythiaGlobalSt,"pythia.spherical_harmonics.abs_neighbor_average"],
                        [instBispectrum,"pythia.spherical_harmonics.bispectrum"]]
        for i in range(len(descriptorsList)):
            print(len(descriptorsList[i][0]),descriptorsList[i][1])
        ###### BISPECTRUM #####
        #Computes bispectrum invariants of particle local environments.
        #These are rotationally-invariant descriptions similar to a power spectrum of the spherical harmonics
        #(i.e. steinhardt order parameters), but retaining more information.
        #print(instBispectrum,len(instBispectrum),instBispectrum.shape,len(pOut))
        #print(instBispectrum[0],instBispectrum[0].shape)
        plotStruct = False
        if (plotStruct==True):
            paletteEnergy = sns.cubehelix_palette(start=.5, rot=-.75, n_colors=len(intsPythiaGlobalSt)).as_hex()
            fig = plt.figure()
            ax = fig.gca(projection='3d')
            xs = pOut[:,0]
            ys = pOut[:,1]
            zs = pOut[:,2]
            def float_to_hex(f):
                return hex(struct.unpack('<I', struct.pack('<f', f))[0])
            ax.set_xlabel('X axis')
            ax.set_ylabel('Y axis')
            ax.set_zlabel('Z axis')
            ax.scatter(xs, ys, zs)
            plt.show()
        if (pyscalHelpers==True):
            pyScalHelp = []
            testStruct = False
            sys = pc.System()
            if (testStruct==True):
                atoms, box = pcs.make_crystal('fcc', lattice_constant=3.147, repetitions=[4,4,4])
                sys.atoms = atoms
                sys.box = box
                sys.find_neighbors(method='cutoff', cutoff='adaptive')
            boxMat = [[mainSizeB,0,0],[0,mainSizeB,0],[0,0,mainSizeB]]
            sys.box = boxMat
            listA = []
            for i in range(len(pOut)):
                ains = pc.Atom(pos=pOut[i], id=0)
                listA.append(ains)
            sys.atoms = listA
            sys.find_neighbors(method='voronoi')
            coord = [atom.coordination for atom in sys.atoms]
            plotCoord = False
            if (plotCoord == True):
                nos, counts = np.unique(coord, return_counts=True)
                plt.bar(nos, counts, color="#AD1457")
                plt.ylabel("density")
                plt.xlabel("coordination number")
                plt.title("Cutoff method")
                plt.show()
                print(coord)
            sys.calculate_entropy(mainSizeB, averaged=True, local=True)
            sys.calculate_angularcriteria()
            atoms = sys.atoms
            angular = [atom.angular for atom in atoms]
            solid_entropy = [atom.entropy for atom in sys.atoms]
            solid_avg_entropy = [atom.avg_entropy for atom in sys.atoms]
            ########  FIND FIRST IF IT IS A SOLID ########
            sys.find_solids(bonds=6, threshold=0.5, avgthreshold=0.6, cluster=True)
            solids = [atom.solid for atom in sys.atoms]
            pyScalHelp = [[angular,"angles"],
            [solid_entropy,"solidEntropy"],
            [solid_avg_entropy,"solidEntropyAvg"],
            [solids,"solidBool"],
            [coord,"coordination"]
            ]
            features = {}
            for i in pyScalHelp:
                descriptorsList.append(i)
            for i in range(len(descriptorsList)):
                features[f"{descriptorsList[i][1]}"] = descriptorsList[i][0]
            ## strcutures
            structures = {}
            structures["fcc"] = fcc_system#[pOutBox, pOut]
            for name, (box, positions) in structures.items():
                print(name, "has", len(positions), "particles.")
            ##
            structure_features = {}
            for name, (box, positions) in structures.items():
                structure_features[name] = features#get_features(box, positions, name)

            structure_dfs = {}
            ##################### GET FEATURES ######################

            for i, structure in enumerate(structure_features):
                print(i, structure)
                # df = pd.DataFrame.from_dict(structure_features[structure])
                # df["class"] = i
                # structure_dfs[structure] = df
            #print(features)

    else:
        from rascal.representations import SphericalInvariants
        # define the parameters of the spherical expansion
        hypers = dict(soap_type="PowerSpectrum",
                    interaction_cutoff=3.0,
                    max_radial=8,
                    max_angular=6,
                    gaussian_sigma_constant=0.5,
                    gaussian_sigma_type="Constant",
                    cutoff_function_type="RadialScaling",
                    cutoff_smooth_width=0.5,
                    cutoff_function_parameters=
                            dict(
                                    rate=1,
                                    scale=3.5,
                                    exponent=4
                                ),
                    radial_basis="GTO",
                    normalize=True,
                    optimization=
                            dict(
                                    Spline=dict(
                                    accuracy=1.0e-05
                                    )
                                ),
                    compute_gradients=False
                    )
        newPos = []
        length = 2*pOutBox
        unitCell = [[length,length,length],[length,length,length],[length,length,length]]
        aProp = Atoms(cell=unitCell,pbc=[1, 1, 1])
        for i in range(len(pOut)):
            aProp.append(Atom("Au",tuple(pOut[i])))
        soap = SphericalInvariants(**hypers)
        instSOAP = soap.transform([aProp])
        soapFeatures = instSOAP.get_features(soap)
        print(len(soapFeatures[0]),(soapFeatures.shape[0]),len(pOut))
        print(soapFeatures)
        print(dir(instSOAP))
    return rotateEnvironment(pOut,rotationsGlobal)

def computeSVC():
    #https://scikit-learn.org/stable/modules/generated/sklearn.svm.SVC.html
    import numpy as np
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    X = np.array([[-1, -1], [-2, -1], [1, 1], [2, 1]])
    y = np.array([1, 1, 2, 2])
    from sklearn.svm import SVC
    clf = make_pipeline(StandardScaler(), SVC(gamma='auto'))
    out = clf.fit(X, y)
    print(clf.predict([[1, 1]]))
    #print(out)
    #Pipeline(steps=[('standardscaler', StandardScaler()),('svc', SVC(gamma='auto'))])
#testQlinFCC()
helpers = True
pyHelp = True
testRotationsInFCC(helpers,pyHelp)
computeSVC()
############################## TEST SECTION ################################
#runHoomd(0.6,7)
exit()
############## PARALLELIZATION SECTION VARYING just ONE parameter ##########
k_Value = [round(float(i),2) for i in list(np.arange(kTarget[0],kTarget[1],1))] # intervals to sample variable-n
# if runed please modify runHOOMDSimulations() to receive just one variable as input or set the second one to default
for i in range(len(k_Value)):
    jobs = []
    p = multiprocessing.Process(target=runAll, args=(k_Value[i],)) # assign threads
    jobs.append(p)
    p.start() # start threads
