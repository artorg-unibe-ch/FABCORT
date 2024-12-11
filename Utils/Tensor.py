#%% !/usr/bin/env python3

Description = """
Class to perform tensor algebra
"""

__author__ = ['Mathieu Simon']
__date_created__ = '05-12-2024'
__date__ = '05-12-2024'
__license__ = 'MIT'
__version__ = '1.0'


import numpy as np

class Tensor():

    def __init__(self):
        pass

    def Dyadic(self, A, B):

        if A.size == 3:

            C = np.zeros((3,3))

            for i in range(3):
                for j in range(3):
                    C[i,j] = A[i]*B[j]

        elif A.size == 9:

            C = np.zeros((3,3,3,3))

            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        for l in range(3):
                            C[i,j,k,l] = A[i,j] * B[k,l]


        else:
            print('Matrices sizes mismatch')

        return C

    def Fabric(self, eValues, eVectors):

        M = np.zeros((3,3))
        for m in range(3):
            M += eValues[m] * self.Dyadic(eVectors[m], eVectors[m])

        return M

    def Ellipsoid(self, EigenValues, EigenVectors, NPoints=100):

        # New coordinate system
        Q = np.array(EigenVectors)

        # Build ellipsoid
        u = np.arange(0, 2 * np.pi + 2 * np.pi / NPoints, 2 * np.pi / NPoints)
        v = np.arange(0, np.pi + np.pi / NPoints, np.pi / NPoints)
        X = EigenValues[0] * np.outer(np.cos(u), np.sin(v))
        Y = EigenValues[1] * np.outer(np.sin(u), np.sin(v))
        Z = EigenValues[2] * np.outer(np.ones_like(u), np.cos(v))
        nNorm = np.zeros(X.shape)

        for i in range(len(X)):
            for j in range(len(X)):
                [X[i, j], Y[i, j], Z[i, j]] = np.dot([X[i, j], Y[i, j], Z[i, j]], Q)

        return X, Y, Z

    def ProjectEllipsoid(self, A, Normal):
        
        # Projection onto plane
        P = np.eye(3) - Tensor.Dyadic(self, Normal, Normal)
        B = P @ A @ P

        # Compute eigenvalues and vectors
        eVals, eVecs = np.linalg.eig(B)
        eVecs = eVecs.T

        # Keep only relevant values and vectors
        Indices = np.where(eVals != 0)[0]

        return eVals[Indices], eVecs[Indices]
    
    def Engineering2MandelNotation(self, A):

        B = np.zeros((6,6))
        for i in range(6):
            for j in range(6):
                if i < 3 and j >= 3:
                    B[i,j] = A[i,j] * np.sqrt(2)
                elif i >= 3 and j < 3:
                    B[i,j] = A[i,j] * np.sqrt(2)
                elif i >= 3 and j >= 3:
                    B[i,j] = A[i,j] * 2
                else:
                    B[i, j] = A[i, j]

        return B

    def Mandel2EngineeringNotation(self, A):

        B = np.zeros((6,6))

        for i in range(6):
            for j in range(6):

                if i < 3 and j >= 3:
                    B[i,j] = A[i,j] / np.sqrt(2)

                elif i >= 3 and j < 3:
                    B[i,j] = A[i,j] / np.sqrt(2)

                elif i >= 3 and j >= 3:
                    B[i,j] = A[i,j] / 2

                else:
                    B[i, j] = A[i, j]

        return B

    def CheckMinorSymmetry(self, A):
        MinorSymmetry = True
        for i in range(3):
            for j in range(3):
                PartialTensor = A[:,:, i, j]
                if PartialTensor[1, 0] == PartialTensor[0, 1] and PartialTensor[2, 0] == PartialTensor[0, 2] and PartialTensor[1, 2] == PartialTensor[2, 1]:
                    MinorSymmetry = True
                else:
                    MinorSymmetry = False
                    break

        if MinorSymmetry == True:
            for i in range(3):
                for j in range(3):
                    PartialTensor = np.squeeze(A[i, j,:,:])
                    if PartialTensor[1, 0] == PartialTensor[0, 1] and PartialTensor[2, 0] == PartialTensor[0, 2] and PartialTensor[1, 2] == PartialTensor[2, 1]:
                        MinorSymmetry = True
                    else:
                        MinorSymmetry = False
                        break

        return MinorSymmetry

    def IsoMorphism66_3333(self, A):

        # Check symmetry
        Symmetry = True
        for i in range(6):
            for j in range(6):
                if not A[i,j] == A[j,i]:
                    Symmetry = False
                    break
        if Symmetry == False:
            print('Matrix is not symmetric!')
            return

        B = np.zeros((3,3,3,3))

        # Build 4th tensor
        B[0, 0, 0, 0] = A[0, 0]
        B[1, 1, 0, 0] = A[1, 0]
        B[2, 2, 0, 0] = A[2, 0]
        B[1, 2, 0, 0] = A[3, 0] / np.sqrt(2)
        B[2, 0, 0, 0] = A[4, 0] / np.sqrt(2)
        B[0, 1, 0, 0] = A[5, 0] / np.sqrt(2)

        B[0, 0, 1, 1] = A[0, 1]
        B[1, 1, 1, 1] = A[1, 1]
        B[2, 2, 1, 1] = A[2, 1]
        B[1, 2, 1, 1] = A[3, 1] / np.sqrt(2)
        B[2, 0, 2, 1] = A[4, 1] / np.sqrt(2)
        B[0, 1, 2, 1] = A[5, 1] / np.sqrt(2)

        B[0, 0, 2, 2] = A[0, 2]
        B[1, 1, 2, 2] = A[1, 2]
        B[2, 2, 2, 2] = A[2, 2]
        B[1, 2, 2, 2] = A[3, 2] / np.sqrt(2)
        B[2, 0, 2, 2] = A[4, 2] / np.sqrt(2)
        B[0, 1, 2, 2] = A[5, 2] / np.sqrt(2)

        B[0, 0, 1, 2] = A[0, 3] / np.sqrt(2)
        B[1, 1, 1, 2] = A[1, 3] / np.sqrt(2)
        B[2, 2, 1, 2] = A[2, 3] / np.sqrt(2)
        B[1, 2, 1, 2] = A[3, 3] / 2
        B[2, 0, 1, 2] = A[4, 3] / 2
        B[0, 1, 1, 2] = A[5, 3] / 2

        B[0, 0, 2, 0] = A[0, 4] / np.sqrt(2)
        B[1, 1, 2, 0] = A[1, 4] / np.sqrt(2)
        B[2, 2, 2, 0] = A[2, 4] / np.sqrt(2)
        B[1, 2, 2, 0] = A[3, 4] / 2
        B[2, 0, 2, 0] = A[4, 4] / 2
        B[0, 1, 2, 0] = A[5, 4] / 2

        B[0, 0, 0, 1] = A[0, 5] / np.sqrt(2)
        B[1, 1, 0, 1] = A[1, 5] / np.sqrt(2)
        B[2, 2, 0, 1] = A[2, 5] / np.sqrt(2)
        B[1, 2, 0, 1] = A[3, 5] / 2
        B[2, 0, 0, 1] = A[4, 5] / 2
        B[0, 1, 0, 1] = A[5, 5] / 2



        # Add minor symmetries ijkl = ijlk and ijkl = jikl

        B[0, 0, 0, 0] = B[0, 0, 0, 0]
        B[0, 0, 0, 0] = B[0, 0, 0, 0]

        B[0, 0, 1, 0] = B[0, 0, 0, 1]
        B[0, 0, 0, 1] = B[0, 0, 0, 1]

        B[0, 0, 1, 1] = B[0, 0, 1, 1]
        B[0, 0, 1, 1] = B[0, 0, 1, 1]

        B[0, 0, 2, 1] = B[0, 0, 1, 2]
        B[0, 0, 1, 2] = B[0, 0, 1, 2]

        B[0, 0, 2, 2] = B[0, 0, 2, 2]
        B[0, 0, 2, 2] = B[0, 0, 2, 2]

        B[0, 0, 0, 2] = B[0, 0, 2, 0]
        B[0, 0, 2, 0] = B[0, 0, 2, 0]



        B[0, 1, 0, 0] = B[0, 1, 0, 0]
        B[1, 0, 0, 0] = B[0, 1, 0, 0]

        B[0, 1, 1, 0] = B[0, 1, 0, 1]
        B[1, 0, 0, 1] = B[0, 1, 0, 1]

        B[0, 1, 1, 1] = B[0, 1, 1, 1]
        B[1, 0, 1, 1] = B[0, 1, 1, 1]

        B[0, 1, 2, 1] = B[0, 1, 1, 2]
        B[1, 0, 1, 2] = B[0, 1, 1, 2]

        B[0, 1, 2, 2] = B[0, 1, 2, 2]
        B[1, 0, 2, 2] = B[0, 1, 2, 2]

        B[0, 1, 0, 2] = B[0, 1, 2, 0]
        B[1, 0, 2, 0] = B[0, 1, 2, 0]



        B[1, 1, 0, 0] = B[1, 1, 0, 0]
        B[1, 1, 0, 0] = B[1, 1, 0, 0]

        B[1, 1, 1, 0] = B[1, 1, 0, 1]
        B[1, 1, 0, 1] = B[1, 1, 0, 1]

        B[1, 1, 1, 1] = B[1, 1, 1, 1]
        B[1, 1, 1, 1] = B[1, 1, 1, 1]

        B[1, 1, 2, 1] = B[1, 1, 1, 2]
        B[1, 1, 1, 2] = B[1, 1, 1, 2]

        B[1, 1, 2, 2] = B[1, 1, 2, 2]
        B[1, 1, 2, 2] = B[1, 1, 2, 2]

        B[1, 1, 0, 2] = B[1, 1, 2, 0]
        B[1, 1, 2, 0] = B[1, 1, 2, 0]



        B[1, 2, 0, 0] = B[1, 2, 0, 0]
        B[2, 1, 0, 0] = B[1, 2, 0, 0]

        B[1, 2, 1, 0] = B[1, 2, 0, 1]
        B[2, 1, 0, 1] = B[1, 2, 0, 1]

        B[1, 2, 1, 1] = B[1, 2, 1, 1]
        B[2, 1, 1, 1] = B[1, 2, 1, 1]

        B[1, 2, 2, 1] = B[1, 2, 1, 2]
        B[2, 1, 1, 2] = B[1, 2, 1, 2]

        B[1, 2, 2, 2] = B[1, 2, 2, 2]
        B[2, 1, 2, 2] = B[1, 2, 2, 2]

        B[1, 2, 0, 2] = B[1, 2, 2, 0]
        B[2, 1, 2, 0] = B[1, 2, 2, 0]



        B[2, 2, 0, 0] = B[2, 2, 0, 0]
        B[2, 2, 0, 0] = B[2, 2, 0, 0]

        B[2, 2, 1, 0] = B[2, 2, 0, 1]
        B[2, 2, 0, 1] = B[2, 2, 0, 1]

        B[2, 2, 1, 1] = B[2, 2, 1, 1]
        B[2, 2, 1, 1] = B[2, 2, 1, 1]

        B[2, 2, 2, 1] = B[2, 2, 1, 2]
        B[2, 2, 1, 2] = B[2, 2, 1, 2]

        B[2, 2, 2, 2] = B[2, 2, 2, 2]
        B[2, 2, 2, 2] = B[2, 2, 2, 2]

        B[2, 2, 0, 2] = B[2, 2, 2, 0]
        B[2, 2, 2, 0] = B[2, 2, 2, 0]



        B[2, 0, 0, 0] = B[2, 0, 0, 0]
        B[0, 2, 0, 0] = B[2, 0, 0, 0]

        B[2, 0, 1, 0] = B[2, 0, 0, 1]
        B[0, 2, 0, 1] = B[2, 0, 0, 1]

        B[2, 0, 1, 1] = B[2, 0, 1, 1]
        B[0, 2, 1, 1] = B[2, 0, 1, 1]

        B[2, 0, 2, 1] = B[2, 0, 1, 2]
        B[0, 2, 1, 2] = B[2, 0, 1, 2]

        B[2, 0, 2, 2] = B[2, 0, 2, 2]
        B[0, 2, 2, 2] = B[2, 0, 2, 2]

        B[2, 0, 0, 2] = B[2, 0, 2, 0]
        B[0, 2, 2, 0] = B[2, 0, 2, 0]


        # Complete minor symmetries
        B[0, 2, 1, 0] = B[0, 2, 0, 1]
        B[0, 2, 0, 2] = B[0, 2, 2, 0]
        B[0, 2, 2, 1] = B[0, 2, 1, 2]

        B[1, 0, 1, 0] = B[1, 0, 0, 1]
        B[1, 0, 0, 2] = B[1, 0, 2, 0]
        B[1, 0, 2, 1] = B[1, 0, 1, 2]

        B[2, 1, 1, 0] = B[2, 1, 0, 1]
        B[2, 1, 0, 2] = B[2, 1, 2, 0]
        B[2, 1, 2, 1] = B[2, 1, 1, 2]


        # Add major symmetries ijkl = klij
        B[0, 1, 1, 1] = B[1, 1, 0, 1]
        B[1, 0, 1, 1] = B[1, 1, 1, 0]

        B[0, 2, 1, 1] = B[1, 1, 0, 2]
        B[2, 0, 1, 1] = B[1, 1, 2, 0]


        return B

    def IsoMorphism3333_66(self, A):

        if self.CheckMinorSymmetry == False:
            print('Tensor does not present minor symmetry')
        else:

            B = np.zeros((6,6))

            B[0, 0] = A[0, 0, 0, 0]
            B[0, 1] = A[0, 0, 1, 1]
            B[0, 2] = A[0, 0, 2, 2]
            B[0, 3] = np.sqrt(2) * A[0, 0, 1, 2]
            B[0, 4] = np.sqrt(2) * A[0, 0, 2, 0]
            B[0, 5] = np.sqrt(2) * A[0, 0, 0, 1]

            B[1, 0] = A[1, 1, 0, 0]
            B[1, 1] = A[1, 1, 1, 1]
            B[1, 2] = A[1, 1, 2, 2]
            B[1, 3] = np.sqrt(2) * A[1, 1, 1, 2]
            B[1, 4] = np.sqrt(2) * A[1, 1, 2, 0]
            B[1, 5] = np.sqrt(2) * A[1, 1, 0, 1]

            B[2, 0] = A[2, 2, 0, 0]
            B[2, 1] = A[2, 2, 1, 1]
            B[2, 2] = A[2, 2, 2, 2]
            B[2, 3] = np.sqrt(2) * A[2, 2, 1, 2]
            B[2, 4] = np.sqrt(2) * A[2, 2, 2, 0]
            B[2, 5] = np.sqrt(2) * A[2, 2, 0, 1]

            B[3, 0] = np.sqrt(2) * A[1, 2, 0, 0]
            B[3, 1] = np.sqrt(2) * A[1, 2, 1, 1]
            B[3, 2] = np.sqrt(2) * A[1, 2, 2, 2]
            B[3, 3] = 2 * A[1, 2, 1, 2]
            B[3, 4] = 2 * A[1, 2, 2, 0]
            B[3, 5] = 2 * A[1, 2, 0, 1]

            B[4, 0] = np.sqrt(2) * A[2, 0, 0, 0]
            B[4, 1] = np.sqrt(2) * A[2, 0, 1, 1]
            B[4, 2] = np.sqrt(2) * A[2, 0, 2, 2]
            B[4, 3] = 2 * A[2, 0, 1, 2]
            B[4, 4] = 2 * A[2, 0, 2, 0]
            B[4, 5] = 2 * A[2, 0, 0, 1]

            B[5, 0] = np.sqrt(2) * A[0, 1, 0, 0]
            B[5, 1] = np.sqrt(2) * A[0, 1, 1, 1]
            B[5, 2] = np.sqrt(2) * A[0, 1, 2, 2]
            B[5, 3] = 2 * A[0, 1, 1, 2]
            B[5, 4] = 2 * A[0, 1, 2, 0]
            B[5, 5] = 2 * A[0, 1, 0, 1]

            return B
        
    def Transform(self, A, OriginalBasis, NewBasis):

        # Build change of coordinate matrix
        O = OriginalBasis
        N = NewBasis

        Q = np.zeros((3,3))
        for i in range(3):
            for j in range(3):
                Q[i,j] = np.dot(O[i,:],N[j,:])

        if A.size == 36:
            A4 = self.IsoMorphism66_3333(A)

        elif A.size == 81 and A.shape == (3,3,3,3):
            A4 = A

        TransformedA = np.zeros((3, 3, 3, 3))
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        for m in range(3):
                            for n in range(3):
                                for o in range(3):
                                    for p in range(3):
                                        TransformedA[i, j, k, l] += Q[m,i]*Q[n,j]*Q[o,k]*Q[p,l] * A4[m, n, o, p]
        if A.size == 36:
            TransformedA = self.IsoMorphism3333_66(TransformedA)

        return TransformedA

