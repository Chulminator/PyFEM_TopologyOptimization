# ---------------------------------------------------------------- 
# Written by: Chulmin Kweon  (Chulmin.Kweon@columbia.edu)
# ----------------------------------------------------------------       

import numpy as np
import sys
sys.path.append("./include")  # Add the include directory to the system path for importing modules
import time as ct

from ReadAnalysis import *  # Import functions related to reading analysis data
import WriteResult  # Module for writing output results

from LinearElasticity import *  # Import linear elasticity analysis functions
import FemBulk as Bulk  # Import bulk FEM-related functions
import importlib
import pandas as pd
from MeshHandler import *  # Import mesh handling functions
from DataStructure import *  # Import data structure definitions

# Set the maximum number of iterations for optimization
max_itr = 20  # Global maximum number of iterations for optimization
tic = ct.time()  # Start the timer

# Ensure that the input arguments are correct
# The input file should be provided for the analysis
assert len(sys.argv) == 2, "check the input value\n execution example:\n\t python3 MainTopologyOptimization.py ./input/PlateWithHole/PlateWithHole_113_93_plastic.inp"

# Model generation and setup
# Initialize program attributes
input_name = sys.argv[1]
Program = ProgramAttrib()

# Read input file to generate the model
# The input file defines the properties of the model such as geometry, material properties, boundary conditions, etc.
Model = ReadInput(input_name, Program)

# Print the initialization information
print("\n" + "=" * 41)
print("| FEM analysis is starting in 2 seconds |")
print("=" * 41)
ct.sleep(2)

##########################################################################################
# FEM Nonlinear Analysis
# Setup elements for the model based on the input configuration
ElementSetUp(Program, Model)

# Set up attributes for writing results based on model dimension (2D or 3D)
if Model.Dim == 2:
    WriteAttributes = np.zeros((Model.totalstep, 11))
elif Model.Dim == 3:
    WriteAttributes = np.zeros((Model.totalstep, 17))

# Apply boundary conditions
# Boundary conditions are imported from a specified file
if Model.BCFile == True:
    tmp = importlib.import_module(Model.BCFileDirectory)
    tmp.ApplyBC(Model)
else:
    assert False, "Other way to apply boundary is not implemented yet"

# Initialize step count and setup boundary conditions
Model.step = 0
Bulk.BC_SetUpAtStep(Model)
CalculateNodalValues(Model)  # Calculate the initial nodal values

# Write initial model state to output file
WriteResult.WriteVTU(Program, Model)

# Start topology optimization iterations
# The outer loop represents multiple optimization steps to update material distribution
for step in range(7):
    # Model step increment (FE analysis)
    Model.step = 1  # Set the finite element analysis step

    # Penalization parameter for material properties
    # It is updated with each iteration to adjust material interpolation for optimization
    Model.penal = 1. + 0.5 * step  # Penalization factor to enforce 0/1 density

    # Nested iterations for density update (using Optimality Criteria or OC method)
    for jj in range(150):
        Program.PrintCommand(str(step + 1) + " - step: " + str(Model.TOstep) + " | penalty: " + str(Model.penal), 0)
        Model.TOstep += 1

        # Sensitivity analysis
        # `zz` and `yy` represent the material densities
        zz = np.zeros(Model.NElem)  # Initialize the element density array
        yy = np.zeros(Model.NElem)  # Initialize the filtered density array

        # Extract material density values for each element
        for ind1, element in enumerate(Model.Element):
            zz[ind1] = element.zz
        # Apply filter to density values (smooth distribution)
        yy = Model.P @ zz

        # Initialize derivative arrays for sensitivity analysis
        dEEdyy = np.zeros(Model.NElem)  # Derivative of elastic energy (compliance) with respect to density
        dVdyy = np.ones(Model.NElem)  # Derivative of volume with respect to density

        # Calculate energy sensitivity with respect to filtered density
        dEEdyy = (1. - 1e-4) * Model.penal * yy ** (Model.penal - 1.)

        # Update material properties for each element
        for ind1, element in enumerate(Model.Element):
            element.yy = yy[ind1]
            element.EE = 1e-4 + (1 - 1e-4) * element.yy ** Model.penal

        # Reset elements for the next iteration
        ElementReset(Program, Model)
        Bulk.BC_SetUpAtStep(Model)

        # Finite Element Analysis
        # Perform finite element analysis to calculate displacements
        for ii in range(30):
            Model.GlobalNRStep = ii
            ConstructStiffness(Model)  # Construct global stiffness matrix
            Global_K = Bulk.AssembleageDisp(Model)  # Assemble the displacement stiffness matrix
            loss = Solve(Model, Global_K)  # Solve the FEM system for displacements
            Program.PrintCommand("itr " + str(ii + 1) + " loss: " + str(loss), 1)
            if loss < 1e-10 and ii > 0:  # Convergence check
                break

        # Objective function and Constraint function
        dfdEE = np.zeros(Model.NElem)  # Derivative of objective (energy) with respect to elastic modulus
        dfdV = np.zeros(Model.NElem)  # Derivative of objective with respect to volume
        dgdEE = np.zeros(Model.NElem)  # Derivative of volume constraint with respect to elastic modulus
        dgdV = np.zeros(Model.NElem)  # Derivative of volume constraint with respect to volume

        # Calculate the compliance (objective sfunction)
        energy = 0.
        for node in Model.Node:
            for ind1 in range(Model.Dim):
                energy += node.u[ind1] * node.F_ext[ind1]

        # Calculate sensitivities for each element
        for ind3, element in enumerate(Model.Element):
            u = np.zeros(element.NNode * element.Dim)
            for ind1, node in enumerate(element.Connectivity):
                for ind2 in range(element.Dim):
                    u[ind1 * element.Dim + ind2] = node.u[ind2]
            dfdEE[ind3] = -u @ element.Stiffness @ u  # Compliance sensitivity (negative stiffness influence)
            #dfdV[ind3] = 0.
            #dgdEE[ind3] = 0.
            dgdV[ind3] = element.Area / Model.Area  # Volume sensitivity

        # Volume constraint calculation
        g = sum(element.Area * element.yy for element in Model.Element) / Model.Area - Model.VolFrac

        # Compute design sensitivity
        dfdzz = Model.P.T @ (dfdEE * dEEdyy + dfdV * dVdyy)  # Sensitivity of objective wrt design variables
        dgdzz = Model.P.T @ (dgdEE * dEEdyy + dgdV * dVdyy)  # Sensitivity of constraint wrt design variables

        # Update density using Optimality Criteria (OC) method
        OCmove = 0.2  # Move limit for density update
        OCeta = 0.5  # Exponent for OC update

        move = OCmove * (1. - 0.)
        eta = OCeta

        l1 = 0
        l2 = 1e+6

        # Use bisection method to update densities (Optimality Criteria loop)
        while l2 - l1 > 1e-4:
            zznew = np.zeros(Model.NElem)

            lmid = 0.5 * (l1 + l2)
            B = -(dfdzz / dgdzz) / lmid

            for ind, element in enumerate(Model.Element):
                zzCnd = 0.0 + (zz[ind] - 0.0) * B[ind] ** eta
                element.zz = min(min(max(max(zzCnd, zz[ind] - move), 0), zz[ind] + move), 1)
                zznew[ind] = element.zz
            dzz = zznew - zz

            # Check volume constraint satisfaction
            if g + np.dot(dgdzz, dzz) > 0:
                l1 = lmid
            else:
                l2 = lmid

        # Print results for each iteration
        print("\tObjective:  ", energy)
        print("\tConstraint: ", g)
        print("\tChange:     ", np.max(np.abs(dzz)))
        print("="*30)

        # Check convergence criteria for design variables
        if np.max(np.abs(dzz)) < Model.toler:
            break

        # Update nodal values and write results to output
        CalculateNodalValues(Model)
        WriteResult.WriteVTU(Program, Model)

# Write final result and metadata
Model.totalstep = Model.TOstep
WriteResult.WritePVD(Program, Model)

# End of program, print total execution time
toc = ct.time() - tic
Program.PrintCommand(' ',0)
Program.PrintCommand("Elapsed CPU time: "+str(toc)+"[sec]",0)
Program.PrintCommand("file log is save in ./log/" +Program.title + ".dat",0)
Program.LogClose()
