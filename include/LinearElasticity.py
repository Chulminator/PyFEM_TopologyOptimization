# ----------------------------------------------------------------
# Written by: CK in promechanics.org at Columbia University
# ----------------------------------------------------------------
import numpy as np
import sys
import time as ct
import os.path
import math
from numpy import ix_
from FemBulk import *
# from util.tensor_operations import *
# from util.coordinate_transforms import *
# from util.EigCalc import *
import importlib
from scipy import linalg
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve


class GPstate:
    eps: np.array
    deps: np.array
    deps: np.array
    dstress: np.array
    lamda: float
    stress: np.array


def ElementSetUp(Program, Model):

    #ModelElement =
    Dim          = Model.Dim
    Model.Area   = 0.
    for Element in Model.Element:
        Edof = Element.Connectivity

        NodeTarget = []
        for Node in Edof:
            NodeTarget.append(Node.Coord)

        Init_tensor = []
        Init_tensor1 = []
        Init_tensor2 = []
        Init_tensor3 = []
        B_matrix  = []
        N_matrix  = []
        Jacc_Elem = []
        Area = 0.

        for gp in Element.GP:
            if Model.Dim   == 2:
                Jacc, B = Strain(Element, NodeTarget, gp[0], gp[1])
                N, dN   = ShapeFunction(Element, gp[0], gp[1])

            elif Model.Dim ==3:
                Jacc, B = Strain(Element, NodeTarget, gp[0], gp[1], gp[2])
                N, dN   = ShapeFunction(Element, gp[0], gp[1], gp[2])

            else:
                assert False, "Check Model.Dim"

            Area = Area + Jacc * gp[-1]
            B_matrix.append(B)
            N_matrix.append(N)
            Jacc_Elem.append(Jacc)
            if Element.Dim ==2:
                if Model.twoD == "planestrain":
                    tmp_tensor = np.zeros((3))
                elif Model.twoD == "planestress":
                    tmp_tensor = np.zeros((3))
                else:
                    assert False, "Check Model.twoD"
                tmp_tensor1 = np.zeros((3))
                tmp_tensor1[0] = Model.HSP
                tmp_tensor1[1] = Model.HSP
            elif Element.Dim ==3:
                tmp_tensor = np.zeros((6))
                tmp_tensor1 = np.zeros((6))
                tmp_tensor1[0] = Model.HSP
                tmp_tensor1[1] = Model.HSP
                tmp_tensor1[2] = Model.HSP
            else:
                assert False, "Check Fem.Dimension"
            Init_tensor.append(np.copy(tmp_tensor))
            Init_tensor1.append(0.0)
            Init_tensor2.append(tmp_tensor1)
            Init_tensor3.append(np.copy(tmp_tensor))

        Element.B_matrix = B_matrix
        Element.N_matrix = N_matrix
        Element.Jacc     = Jacc_Elem
        Element.Area     = Area
        Model.Area      += Area

        Element.GPstrain   = Init_tensor
        Element.GPlamda    = Init_tensor1
        Element.GPstress   = Init_tensor2


    ################################################################

    P            = np.zeros((Model.NElem,Model.NElem))

    for ind1, element1 in enumerate( Model.Element ):
        for ind2 in range(ind1, Model.NElem):
            element2 = Model.Element[ind2]
            tmp      = (element1.Centroid - element2.Centroid)
            distance = np.sqrt(np.dot( tmp.T, tmp))
            if distance < Model.R:
               P[ind1,ind2] = 1 - distance/Model.R
               P[ind2,ind1] = 1 - distance/Model.R
    P    = (P.T / np.sum(P,axis=1)).T
    Model.P  = csc_matrix(P)

    #print(np.sum(P,axis=1))
    #print(np.sum(P.T,axis=1))
    #exit(1)
    ################################################################
    # Or
    #Model.P            = np.eye(Model.NElem)
    ################################################################
    #print(np.sum(P,axis=0)) should be 1

    return



def ElementReset(Program, Model):

    #ModelElement =
    Dim          = Model.Dim
    for Element in Model.Element:
        Edof = Element.Connectivity

        NodeTarget = []
        for Node in Edof:
            NodeTarget.append(Node.Coord)

        Init_tensor = []

        for gp in Element.GP:

            if Element.Dim ==2:
                if Model.twoD == "planestrain":
                    tmp_tensor = np.zeros((3))
                elif Model.twoD == "planestress":
                    tmp_tensor = np.zeros((3))
                else:
                    assert False, "Check Model.twoD"
            elif Element.Dim ==3:
                tmp_tensor = np.zeros((6))
            else:
                assert False, "Check Fem.Dimension"
            Init_tensor.append(np.copy(tmp_tensor))

        Element.GPstrain   = Init_tensor

    return




def ConstructStiffness(Model):
    epsil = 1e-4
    Dim    = Model.Dim
    for node in Model.Node:
        node.F_int = np.zeros(Dim)

    for ind1, Element in enumerate(Model.Element):
        G_Edof   = Element.G_Edof
        E_K      = np.zeros([G_Edof, G_Edof])
        B        = Element.B_matrix
        N        = Element.N_matrix
        Jacc     = Element.Jacc
        GPstress = Element.GPstress
        GPeps    = Element.GPstrain
        GPF_int  = np.zeros(G_Edof)
        du       = np.zeros(G_Edof)

        for ind, Node in enumerate(Element.Connectivity):
            du[ind*Dim:(ind+1)*Dim] = Node.du

        for ind, gp in enumerate(Element.GP):
            B_GP          = B[ind]
            N_GP          = N[ind]
            Jacc_GP       = Jacc[ind]
            GPdata        = GPstate()
            GPdata.eps    = GPeps[ind]
            GPdata.deps   = B_GP @ du
            GPdata.stress = GPstress[ind]
            #GPdata.RHO    = Element.RHOhat
            
            Element.ConstitutiveModel.NextStep(GPdata)

            GPeps[ind]    = GPdata.eps
            GPstress[ind] =  GPdata.stress

            if Model.Dim == 2:
                GP_K = B_GP.T @ GPdata.D @ B_GP * Model.width * Jacc_GP * gp[-1]
                GPF_int += B_GP.T @ GPstress[ind] * Model.width * Jacc_GP * gp[-1] * Element.EE
            elif Model.Dim == 3:
                GP_K = B_GP.T @ GPdata.D @ B_GP * Jacc_GP * gp[-1]
                GPF_int += B_GP.T @ GPstress[ind] * Jacc_GP * gp[-1] * element.EE
            else:
                assert False, "Dimension error"

            E_K += GP_K

        for ind, Node in enumerate(Element.Connectivity):
            Node.F_int += GPF_int[ind*Dim:(ind+1)*Dim]

        Element.Stiffness = E_K
    return


def Solve(Model, Global_K):
    # https://mae.ufl.edu/nkim/egm6352/Chap2.pdf
    # Kim, N.H., 2014. Introduction to nonlinear finite element analysis. Springer Science & Business Media.
    # [[K11 K12],  [[d1],     [[Fout1],
    #  [K21 K22]]   [d2]]   =  [Fout2]]
    #
    # K12 * d2 + K11 * d1 = F1
    # d1 = K11 \ ( Fout1 - K12 * d2 )

    Dim = Model.Dim
    IndexBCN =[]
    IndexBCE =[]

    F_ext = np.zeros(Model.NNode * Dim)
    F_int = np.zeros(Model.NNode * Dim)
    u     = np.zeros(Model.NNode * Dim)
    du    = np.zeros(Model.NNode * Dim)
    for ii, node in enumerate(Model.Node):
        for jj in range(Dim):
            u[ii*Dim+jj]     = node.u[jj]
            F_ext[ii*Dim+jj] = node.F_ext[jj]
            F_int[ii*Dim+jj] = node.F_int[jj]

            if np.isnan( node.BC_E[jj] ):
                IndexBCN.append(ii*Dim+jj)
            else:
                IndexBCE.append(ii*Dim+jj)





    u1 = np.copy(u)

    #Sliced_K12 = Global_K[ix_(IndexBCN, IndexBCE)]
    Sliced_K11 = Global_K[ix_(IndexBCN, IndexBCN)]

    F_total = F_ext[IndexBCN] - F_int[IndexBCN]
    #F_total = F_ext[IndexBCN] - Sliced_K12 @ u[IndexBCE]


    Sliced_K11_csc = csc_matrix(Sliced_K11)
    du[IndexBCN] = spsolve(Sliced_K11_csc, F_total)


    u1[IndexBCN] = u[IndexBCN]  + du[IndexBCN]
    du[IndexBCE] = u1[IndexBCE] -  u[IndexBCE]

    u  = np.copy(u1)

    for ii, node in enumerate( Model.Node):
        for jj in range(Dim):
            node.u[jj]  = u[ii*Dim + jj]
            node.du[jj] = du[ii*Dim + jj]

    return np.linalg.norm(F_total)**2/(1.+np.linalg.norm(F_ext[IndexBCN])**2)



class ConstitutiveLaw():
    def __init__(self, Model, Element):
        self.twoD = Model.twoD
        self.HSP = Model.HSP
        self.Dim  = Model.Dim

        self.E       = Element.E
        self.nu      = Element.nu
        self.MatProp = Element.MatProp

    def Voigt2Tensor(self, voigt, flag='strain'):
        # https://www.comsol.com/blogs/what-is-the-difference-between-plane-stress-and-plane-strain
        Tensor = np.zeros((3,3))
        if flag.lower() == 'strain':
            if self.twoD == 'planestrain':
                Tensor[0,0] = voigt[0]
                Tensor[0,1] = voigt[2]*0.5
                Tensor[1,0] = voigt[2]*0.5
                Tensor[1,1] = voigt[1]
                #assert False, "this function is not available"

            elif self.twoD == 'planestress':
                nu = self.nu
                Tensor[0,0] = voigt[0]
                Tensor[0,1] = voigt[2]*0.5
                Tensor[1,0] = voigt[2]*0.5
                Tensor[1,1] = voigt[1]
                Tensor[2,2] = -nu/(1.0-nu)*(Tensor[0,0] + Tensor[1,1])
                #assert False, "this function is not available"
            else:
                if self.Dim == 3:
                    Tensor[0,0] = voigt[0]
                    Tensor[1,1] = voigt[1]
                    Tensor[2,2] = voigt[2]

                    Tensor[0,1] = voigt[3]*0.5 # xy
                    Tensor[1,0] = voigt[3]*0.5

                    Tensor[1,2] = voigt[4]*0.5 # yz
                    Tensor[2,1] = voigt[4]*0.5

                    Tensor[0,2] = voigt[5]*0.5 # zx
                    Tensor[2,0] = voigt[5]*0.5
                else:
                    print("self.twoD", self.twoD)
                    print("self.Dim",  self.Dim)
                    assert False,"Check Fem.twoD or Fem.Dimension"

        elif flag.lower() == 'stress':
            if self.twoD == 'planestrain':
                nu = self.nu
                Tensor[0,0] = voigt[0]
                Tensor[0,1] = voigt[2]
                Tensor[1,0] = voigt[2]
                Tensor[1,1] = voigt[1]
                Tensor[2,2] = nu*(Tensor[0,0] + Tensor[1,1])
                assert False, "this function is not available"

            elif self.twoD == 'planestress':
                nu = self.nu
                Tensor[0,0] = voigt[0]
                Tensor[0,1] = voigt[2]
                Tensor[1,0] = voigt[2]
                Tensor[1,1] = voigt[1]
                assert False, "this function is not available"
                
            else:
                if self.Dim == 3:
                    Tensor[0,0] = voigt[0]
                    Tensor[1,1] = voigt[1]
                    Tensor[2,2] = voigt[2]

                    Tensor[0,1] = voigt[3] # xy
                    Tensor[1,0] = voigt[3]

                    Tensor[1,2] = voigt[4] # yz
                    Tensor[2,1] = voigt[4]

                    Tensor[0,2] = voigt[5] # zx
                    Tensor[2,0] = voigt[5]
                else:
                    assert False,"Check Fem.twoD or Fem.Dimension"
        return Tensor


    def Tensor2Voigt(self, Tensor, flag='strain'):        
        if self.Dim ==2:
            voigt = np.zeros((3))
            if flag == 'strain':
                voigt[0] = Tensor[0,0]
                voigt[1] = Tensor[1,1]
                voigt[2] = Tensor[0,1] + Tensor[1,0]
            elif flag == 'stress':
                voigt[0] = Tensor[0,0]
                voigt[1] = Tensor[1,1]
                voigt[2] = Tensor[0,1]
            else:
                assert False ,"check the flag in Tensor2Voigt"
        elif self.Dim ==3:
            voigt = np.zeros((6))
            if flag == 'strain':
                voigt[0] = Tensor[0,0]
                voigt[1] = Tensor[1,1]
                voigt[2] = Tensor[2,2]
                voigt[3] = Tensor[0,1]+Tensor[1,0]
                voigt[4] = Tensor[1,2]+Tensor[2,1]
                voigt[5] = Tensor[2,0]+Tensor[0,2]

            elif flag == 'stress':
                voigt[0] = Tensor[0,0]
                voigt[1] = Tensor[1,1]
                voigt[2] = Tensor[2,2]
                voigt[3] = Tensor[0,1]
                voigt[4] = Tensor[1,2]
                voigt[5] = Tensor[2,0]
            else:
                assert False ,"check the flag in Tensor2Voigt"
        else:
            assert False,"Check the dimension"
        return voigt


    def InitialPressure(self):
        return np.diag([self.HSP, self.HSP, self.HSP])

    def NextStep(self, GPstate):
        epsil = 1e-4

        eps_init    = self.Voigt2Tensor( GPstate.eps )
        deps        = self.Voigt2Tensor( GPstate.deps )
        #RHO         = GPstate.RHO        # Damage parameter
        eps         = eps_init + deps

        #E   = self.E * (epsil + (1 - epsil) * RHO)
        E   = self.E
        nu  = self.nu

        K     = E / (3*(1.0-2.0*nu))
        mu    = E / (2.0*(1.0+nu))
        lam   = E*nu / ((1.0+nu)*(1.0-2.0*nu))

        # Calculate the stress
        I          = np.eye(3)
        sigma      = lam * np.trace(eps)*I + 2. * mu * eps

        # Calculate the stiffness
        if self.Dim == 2:
            if self.twoD == 'planestress':
                D = np.array([[1., nu, 0.],
                              [nu, 1., 0.],
                              [0., 0., (1.-nu)*0.5]])*E/(1.-nu**2)
            else:
                ddsdde = FthTensor2Voigt(lam*IxI + 2.*mu*II)
                idx = [0,1,3]
                D = ddsdde[ix_(idx, idx)]
        elif self.Dim ==3:
            ddsdde = FthTensor2Voigt(lam*IxI + 2.*mu*II)
            D = ddsdde
        else:
            assert False, "Dimension Check"

        GPstate.stress  = self.Tensor2Voigt(sigma,'stress')
        GPstate.eps     = self.Tensor2Voigt(eps)
        GPstate.D       = D
        return




