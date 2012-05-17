# carry the Walsh analysis information with individual
import numpy as np
import math
import random

import copy

class Struct:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

class Individual:
    def __init__( self, n, neigh=False ):
        self.bit = []
        self.fit = 0
        self.dim = n 
        if neigh == True :
            self.fitG = 0
        self.threshold = 1e-15

    def init(self):
        self.bit = []
        for j in range(self.dim):
            if random.random()<0.5:
                self.bit.append('0')
            else:
                self.bit.append('1')
    
    def initWal(self, model):
        """ 
        1. 
        compute the sum array for the first time, according to the initial solution
        2. 
        initialize a full matrix, Coincidence, C_ij, which is a symmetric one, and diagonal is empty (as we don't want to recompute SumArr), empty whenever i >= j.
        Two options for implementing C matrix
            a. precompute C matrix 
            b. construct it on the fly using dict()
        3. 
        initialize interaction list for each variable
        4.
        initialize a dict of interaction structure, where interactive bits and the index of WAS (walsh including sign)
        """
        self.model = model
        self.sumArr = np.zeros(self.dim)
        self.WAS = np.tile(Struct(arr = [], w = 0), len(model.w.keys()))# Walsh coefficients with sign, represented in Array
        self.lookup = dict()
        self.infectBit = dict()
        self.C = np.zeros((self.dim,self.dim)) # coincidence matrix

        for i in range(len(model.WA)):
            W = int(math.pow(-1, self.binCount(model.WA[i].arr, self.bit))) * model.WA[i].w
            self.WAS[i] = Struct(arr = model.WA[i].arr, w = W)
            comb = self.genComb(len(model.WA[i].arr))

            for j in model.WA[i].arr:
                self.sumArr[j] = self.sumArr[j] + W
                # add list of order >= 3 Walsh terms for the purpose of updating C matrix
                if len(model.WA[i].arr) >= 3:
                    if j not in self.infectBit: 
                        self.infectBit[j] = [Struct(arr=model.WA[i].arr, WI=i)]
                    else :
                        self.infectBit[j].append(Struct(arr=model.WA[i].arr, WI=i))

            for j in comb: # for each list in comb
                j0 = model.WA[i].arr[int(j[0])]
                j1 = model.WA[i].arr[int(j[1])]
                self.C[j0,j1] = self.C[j0,j1] + W

    def initSC(self):
        # compute the SC array
        self.SC = np.zeros(self.dim)
        self.Z = np.zeros(self.dim)
        self.orderC = np.zeros((self.dim,self.dim))

        for p in range(self.dim):

            phi = np.zeros(self.model.k+1)
            if p in self.model.Inter:
                for i in self.model.Inter[p].WI:
                    order = len(self.WAS[i].arr)
                    phi[order-1] = phi[order-1] + self.WAS[i].w

            self.Z[p] = self.sumArr[p]
            for i in range(1, self.model.k+1):
                if phi[i] != 0:
                    self.Z[p] = self.Z[p] + i * phi[i]

            self.SC[p] = self.sumArr[p] - 2/float(self.dim) * self.Z[p]
            #self.SC[i] = self.sumArr[i] - 2/float(self.dim) * self.compCsum(i)

        for i in range(len(self.WAS)):
            lenArr = len(self.WAS[i].arr)
            comb = self.genComb(lenArr)
            for j in comb:
                j0 = self.WAS[i].arr[int(j[0])]
                j1 = self.WAS[i].arr[int(j[1])]
                self.orderC[j0,j1] = self.orderC[j0,j1] + lenArr * self.WAS[i].w

    ## def compPSum(self,bitStr):
    ##     """
    ##     use Elementary Landscape Analysis to obtain the average of neighs of given
    ##     individual
    ##     """
    ##     p = np.zeros(self.model.k+1)
    ##     for k in self.model.w.keys():
    ##         oneC = k.count('1')
    ##         if  oneC !=0 :
    ##             p[oneC-1] = p[oneC-1] + self.model.w[k] * math.pow(-1,wal.bc(k,bitStr))
    ## #            else :
    ## #                p[0] = p[0] + self.model.w[k]

    ##     s = 0
    ##     for i in range(self.model.k+1):
    ##         s = s + (i+1)*p[i]
    ##     return s


    def compPhisum(self,p):
        """
        \varphi_{p,i}^{\prime}(x) = \Sigma_{order j terms, that touches bit p}
        """
        phi = np.zeros(self.model.k+1)
        if p in self.model.Inter:
            for i in self.model.Inter[p].WI:
                order = len(self.WAS[i].arr)
                phi[order-1] = phi[order-1] + self.WAS[i].w

        s = self.sumArr[p]
        for i in range(1, self.model.k+1):
            if phi[i] != 0:
                s = s + i * phi[i]
        return s

    def compCsum(self,p):
        """
        \sigma_{i=1}^{N} C_{ip}: be careful with C_{ii}, i \in N
        """
        s = 0

        for i in range(p):
            s = s + self.C[i,p]

        for i in range(p+1, self.dim):
            s = s + self.C[p,i]

        s = s + self.sumArr[p]

        return s

    def update(self, p):
        """
        By keeping track of coincidence matrix, 
        Cij stands for S_i(y_j) = S_i(x) - C_ij
        partially update the Sum Array and self.WAS, given the bit which is changed
        """
        self.sumArr[p] = - self.sumArr[p]
        
        if p in self.model.Inter:
            for i in self.model.Inter[p].arr:
                if i < p:
                    self.sumArr[i] = self.sumArr[i] - 2*self.C[i,p]
                    self.C[i,p] = - self.C[i,p]
                else:
                    self.sumArr[i] = self.sumArr[i] - 2*self.C[p,i]
                    self.C[p,i] = - self.C[p,i]

        # update the rest of elements in C matrix
        if p in self.infectBit.keys():
            for i in self.infectBit[p]:
                arr = i.arr[:]
                arr.remove(p)
                comb = self.genComb(len(arr))
                for k in range(len(comb)):
                    k0 = arr[int(comb[k][0])]
                    k1 = arr[int(comb[k][1])]
                    self.C[k0,k1] = self.C[k0,k1] - 2 * self.WAS[i.WI].w

    def updateImprS(self, p, minimize):
        self.improveA.remove(p)
        if p in self.model.Inter:
            for i in self.model.Inter[p].arr: 
                if (minimize == True and self.sumArr[i] > - self.threshold) or (minimize == False and self.sumArr[i] < self.threshold ):
                #if (minimize == True and self.sumArr[i] > self.threshold) or (minimize == False and self.sumArr[i]< - self.threshold ):
                    if i not in self.improveA:
                        self.improveA.append(i)
                elif i in self.improveA:
                    self.improveA.remove(i)

    def updatePertImprS(self, p, minimize):
        if p in self.model.Inter:
            for i in self.model.Inter[p].arr : 
                if (minimize == True and self.sumArr[i] > - self.threshold) or (minimize == False and self.sumArr[i]< self.threshold ):
                #if (minimize == True and self.sumArr[i] > self.threshold) or (minimize == False and self.sumArr[i]< - self.threshold ):
                    if i not in self.improveA:
                        self.improveA.append(i)
                elif i in self.improveA:
                    self.improveA.remove(i)
        
        # (minimize == True and self.sumArr[p] > self.threshold) or (minimize == False and self.sumArr[p]< - self.threshold ):
        if (minimize == True and self.sumArr[p] > - self.threshold) or (minimize == False and self.sumArr[p] < self.threshold ):
            if p not in self.improveA:
                self.improveA.append(p)
        elif p in self.improveA:
            self.improveA.remove(p)

    def updateFake(self, p):
        """
        The fake version, the updates are made to the mirror data structures
        By keeping track of coincidence matrix, 
        Cij stands for S_i(y_j) = S_i(x) - C_ij
        partially update the Sum Array and self.WAS, given the bit which is changed
        """
        self.sumArrFake[p] = - self.sumArrFake[p]
        
        if p in self.model.Inter:
            for i in self.model.Inter[p].arr:
                if i < p:
                    self.sumArrFake[i] = self.sumArrFake[i] - 2*self.Cfake[i,p]
                else:
                    self.sumArrFake[i] = self.sumArrFake[i] - 2*self.Cfake[p,i]
    def updateWAS(self,p):
        if p in self.model.Inter:
            for i in self.model.Inter[p].WI:
                self.WAS[i].w = - self.WAS[i].w


    def updateSC(self, p):
        self.SC[p] = - self.SC[p]
        self.Z[p] = - self.Z[p]

        #update Z array
        if p in self.model.Inter:
            for i in self.model.Inter[p].arr:
                if i < p :
                    self.Z[i] = self.Z[i]  - 2* self.orderC[i,p]
                    self.orderC[i,p] = - self.orderC[i,p]
                else :
                    self.Z[i] = self.Z[i]  - 2* self.orderC[p,i]
                    self.orderC[p,i] = - self.orderC[p,i]
                self.SC[i] = self.sumArr[i] - 2/float(self.dim) * self.Z[i]

        if p in self.infectBit.keys():
            for i in self.infectBit[p]:
                arr = copy.deepcopy(i.arr)
                arr.remove(p)
                lenArr = len(arr)
                comb = self.genComb(lenArr)
                for k in range(len(comb)):
                    k0 = arr[int(comb[k][0])]
                    k1 = arr[int(comb[k][1])]
                    self.orderC[k0,k1] = self.orderC[k0,k1] - 2 * (lenArr + 1)* self.WAS[i.WI].w

    def updateImprSC(self, p, minimize):
        self.improveA.remove(p)
        if p in self.model.Inter:
            for i in self.model.Inter[p].arr:
                # if (minimize == True and self.SC[i] > self.threshold) or (minimize == False and self.SC[i] < - self.threshold):
                if (minimize == True and self.SC[i] > - self.threshold) or (minimize == False and self.SC[i] < self.threshold):
                    if i not in self.improveA:
                        self.improveA.append(i)
                elif i in self.improveA:
                    self.improveA.remove(i)

    def updatePertImprSC(self, p, minimize):
        if p in self.model.Inter:
            for i in self.model.Inter[p].arr:
                # if (minimize == True and self.SC[i] > self.threshold) or (minimize == False and self.SC[i] < - self.threshold):
                if (minimize == True and self.SC[i] > - self.threshold) or (minimize == False and self.SC[i] < self.threshold):
                    if i not in self.improveA:
                        self.improveA.append(i)
                elif i in self.improveA:
                    self.improveA.remove(i)
        #if (minimize == True and self.SC[p] > self.threshold) or (minimize == False and self.SC[p] < - self.threshold):
        if (minimize == True and self.SC[p] > - self.threshold) or (minimize == False and self.SC[p] < self.threshold):
            if p not in self.improveA:
                self.improveA.append(p)
        elif p in self.improveA:
            self.improveA.remove(p)

    def genImproveS(self,minimize):
        """
        generate the index of best neigh according to sumArr only (surrogate of fitness)
        """
        # check improving move
        self.improveA = []
        for i in range(self.dim):
            # if (minimize == True and self.sumArr[i] > self.threshold) or (minimize == False and self.sumArr[i]< - self.threshold):
            if (minimize == True and self.sumArr[i] > - self.threshold) or (minimize == False and self.sumArr[i] < self.threshold):
                self.improveA.append(i) 

    def genFitNext(self,minimize):
        """
        generate the index of next improving neigh according to sumArr only (surrogate of fitness)
        """
        # check improving move
        improve = False
        self.improveA = []
        for i in range(self.dim):
            #if (minimize == True and self.sumArr[i] > self.threshold) or (minimize == False and self.sumArr[i] < - self.threshold):
            if (minimize == True and self.sumArr[i] > - self.threshold) or (minimize == False and self.sumArr[i] < self.threshold):
                self.improveA.append(i) 
                improve = True

        if improve == False:
            return False, None

        # randomly pick an improving move, which takes only constant time 
        bestI = random.choice(self.improveA)
                    
        return True, bestI

    def genFitBest2(self,minimize):
        """
        generate the index of best distance 2 neighborhoods according to sumArr only (surrogate of fitness)

        return: 1) whether there is an improving move in distance 2 neigh
                2) the index of best distance 1 neigh for taking the next move
                3) the number of evaluations consumed by this step
        """
        improve = False
        self.improveA = []
        neighImprove = []

        # checking the distance 1 neigh
        for i in range(self.dim):
            #if (minimize == True and self.sumArr[i] > self.threshold) or (minimize == False and self.sumArr[i] < -self.threshold):
            if (minimize == True and self.sumArr[i] > - self.threshold) or (minimize == False and self.sumArr[i] < self.threshold ):
                self.improveA.append(i) 
                neighImprove.append(Struct(index =[i], val = self.sumArr[i])) # add dist 1 neigh into consideration as well
                improve = True

#        print self.sumArr
#        print self.improveA

        # checking the distance 2 neigh, remember to preserve context
        for i in range(self.dim):
            self.mirrorParam() # everything is reset, pretending nothing happened
            self.updateFake(i)
            #self.updateWASfake(i)
            for j in [k for k in range(self.dim) if k!=i]:
                self.sumArrFake[j] = self.sumArrFake[j]+self.sumArr[i]
                if (minimize == True and self.sumArrFake[j] > - self.threshold) or (minimize == False and self.sumArrFake[j]< self.threshold):
                #if (minimize == True and self.sumArrFake[j] > self.threshold) or (minimize == False and self.sumArrFake[j]< -self.threshold):
                    neighImprove.append(Struct(index =[i,j], val = self.sumArrFake[j]))
                    improve = True

#        for i in range(len(neighImprove)):
#            print neighImprove[i].index, neighImprove[i].val

#        for i in neighImprove:
#            print i.index, i.val
               
        if improve == False:
            #return False, None, self.dim*self.dim
            return False, None, self.dim

        for i in range(len(neighImprove)):
            if i == 0:
                best = neighImprove[i].val
                bestI = neighImprove[i].index
            # if ( best<neighImprove[i].val - self.threshold and minimize == True) or ( best>neighImprove[i].val + self.threshold and minimize == False ): # seek for max S
            elif ( best<neighImprove[i].val + self.threshold and minimize == True) or ( best>neighImprove[i].val - self.threshold and minimize == False ): # seek for max S
                best = neighImprove[i].val
                bestI = neighImprove[i].index

        bestIlist = []
        for i in range(len(neighImprove)):
            #if abs(best - neighImprove[i].val) < self.threshold:
            if abs(best - neighImprove[i].val) < self.threshold:
                candI = neighImprove[i].index
                if candI not in bestIlist:
                    bestIlist.append(candI)

        #print 'bestIlist',bestIlist
        bestI = random.choice(bestIlist)
#        print 'bestList', bestIlist
#        print 'bestI', bestI
        if type(bestI) is int:
            # make a consistent interface
            bestI = [bestI]
                    
        #return True, bestI, self.dim*self.dim
        return True, bestI, self.dim

    def genFitBestsm(self,minimize):
        """
        generate the index of best neigh according to sumArr only (surrogate of fitness)
        """
        #TODO need to update the threshold
        # check improving move
        improve = False
        self.Buffer = []
        for i in range(self.dim):
            if (minimize == True and self.sumArr[i] > 0) or (minimize == False and self.sumArr[i]<0):
                self.Buffer.append(i) 
                improve = True

        if improve == False:
            return False, None

        for i in self.Buffer:
            if i == self.Buffer[0]:
                best = self.sumArr[i]
                bestI = i
            elif (best<self.sumArr[i] and minimize == True) or (best>self.sumArr[i] and minimize == False): # seek for max S
                best = self.sumArr[i]
                bestI = i

        self.P = [bestI]

        # produce buffer list (the independent improving set)
        for i in [ j for j in self.Buffer if j != bestI ]:
            if i not in self.Inter:
                self.P.append(i)
            else :
                inter = False
                for j in self.P:
                    if j in self.Inter[i].arr:
                        inter = True
                if inter == False:
                    self.P.append(i)

        return True, bestI

    def updateFitBestsm(self, minimize):
        #TODO need to update the threshold
        for i in self.P:
            self.Buffer.remove(i)
            if i in self.Inter:
                for j in self.Inter[i].arr:
                    if ((minimize == True and self.sumArr[j] > 0) or (minimize == False and self.sumArr[j]<0)): 
                        if j not in self.Buffer:
                            self.Buffer.append(j)
                    elif j in self.Buffer:
                        self.Buffer.remove(j)

        if not self.Buffer:
            return False, None
        #print improveA

        for i in self.Buffer:
            if i == self.Buffer[0]:
                best = self.sumArr[i]
                bestI = i
            elif (best<self.sumArr[i] and minimize == True) or (best>self.sumArr[i] and minimize == False): # seek for max S
                best = self.sumArr[i]
                bestI = i

        self.P = [bestI]

        for i in [ j for j in self.Buffer if j != bestI]:
            if i not in self.Inter:
                self.P.append(i)
            else :
                inter = False
                for j in self.P:
                    if j in self.Inter[i].arr:
                        inter = True
                if inter == False:
                    self.P.append(i)

        return True, bestI

    def updateFitBest(self, p, minimize):
        self.improveA.remove(p)

        if p in self.Inter:
            for i in self.Inter[p].arr: 
                # if (minimize == True and self.sumArr[i] > self.threshold) or (minimize == False and self.sumArr[i] < - self.threshold ):
                if (minimize == True and self.sumArr[i] > - self.threshold) or (minimize == False and self.sumArr[i]< self.threshold ):
                    if i not in self.improveA:
                        self.improveA.append(i)
                elif i in self.improveA:
                    self.improveA.remove(i)

        if not self.improveA:
            return False, None

        random.shuffle(self.improveA)

        for i in self.improveA:
            if i == self.improveA[0]:
                best = self.sumArr[i]
                bestI = i
            elif (best<self.sumArr[i] and minimize == True) or (best>self.sumArr[i] and minimize == False): 
                best = self.sumArr[i]
                bestI = i
                    
        return True, bestI

    def steepFitDesc(self, minimize):
        if not self.improveA:
            return False, None

        random.shuffle(self.improveA)

        for i in self.improveA:
            if i == self.improveA[0]:
                best = self.sumArr[i]
                bestI = i
            #elif (best<self.sumArr[i] - self.threshold and minimize == True) or (best>self.sumArr[i] + self.threshold and minimize == False): 
            elif (best<self.sumArr[i] + self.threshold and minimize == True) or (best > self.sumArr[i] - self.threshold and minimize == False): 
                best = self.sumArr[i]
                bestI = i
                    
        return True, bestI

    def steepMeanDesc(self, minimize):
        if not self.improveA:
            return False, None

        random.shuffle(self.improveA)

        # find the best value
        for i in self.improveA:
            if i == self.improveA[0]:
                best = self.SC[i]
                bestI = i
            # elif ( best<self.SC[i] - self.threshold and minimize == True ) or ( best>self.SC[i] + self.threshold and minimize == False ): # seek for max S
            elif ( best<self.SC[i] + self.threshold and minimize == True ) or ( best>self.SC[i] - self.threshold and minimize == False ): # seek for max S
                best = self.SC[i]
                bestI = i
        return True, bestI

    ## def updateFitNext(self, p, minimize):
    ##     """
    ##     find the next improving move by the similar update trick
    ##     """
    ##     self.improveA.remove(p)
    ##     if p in self.Inter:
    ##         for i in self.Inter[p].arr: 
    ##             #if (minimize == True and self.sumArr[i] > self.threshold) or (minimize == False and self.sumArr[i]< - self.threshold ):
    ##             if (minimize == True and self.sumArr[i] > - self.threshold) or (minimize == False and self.sumArr[i]< self.threshold ):
    ##                 if i not in self.improveA:
    ##                     self.improveA.append(i)
    ##             elif i in self.improveA:
    ##                 self.improveA.remove(i)

    ##     if not self.improveA:
    ##         return False, None, evalCount

    ##     # randomly pick an improving move, which takes only constant time 
    ##     bestI = random.choice(self.improveA)
                   
    ##     return True, bestI, evalCount

    def nextDesc(self):
        """
        find the next improving move by the similar update trick
        """
        if not self.improveA:
            return False, None

        # randomly pick an improving move, which takes only constant time 
        bestI = random.choice(self.improveA)
                   
        return True, bestI

    def updateFitBest2(self, P, minimize):
        """ 
        generate the index of best distance 2 neighborhoods according to sumArr only (surrogate of fitness), by performing partial updates
        """
        neighImprove = []
        evalCount = 0
        
        for p in P:
            if p in self.improveA:
                self.improveA.remove(p)
            if p in self.Inter:
                for i in self.Inter[p].arr: 
                    evalCount = evalCount + 1
                    # if (minimize == True and self.sumArr[i] > self.threshold) or (minimize == False and self.sumArr[i]<-self.threshold):
                    if (minimize == True and self.sumArr[i] > - self.threshold) or (minimize == False and self.sumArr[i]< self.threshold):
                        if i not in self.improveA:
                            self.improveA.append(i)
                    elif i in self.improveA:
                        self.improveA.remove(i)
        p = P[-1]

        for i in self.improveA:
            """ add distance 1 neigh under consideration """
            neighImprove.append(Struct(index=[i], val = self.sumArr[i]))

        # checking the distance 2 neigh, remember to preserve context
        for i in [k for k in range(self.dim) if k!=p]:
            self.mirrorParam() # everything is reset, pretending nothing happened
            self.updateFake(i)
            for j in [k for k in range(self.dim) if k!=i]:
                self.sumArrFake[j] = self.sumArrFake[j]+self.sumArr[i]
                # if (minimize == True and self.sumArrFake[j] > self.threshold) or (minimize == False and self.sumArrFake[j]<-self.threshold):
                if (minimize == True and self.sumArrFake[j] > - self.threshold) or (minimize == False and self.sumArrFake[j] < self.threshold):
                    neighImprove.append(Struct(index =[i,j], val = self.sumArrFake[j]))

#        for i in range(len(neighImprove)):
#            print neighImprove[i].index, neighImprove[i].val
        
        if not neighImprove:
            #return False, None, evalCount + self.dim * (self.dim-1)
            return False, None, evalCount 

        for i in range(len(neighImprove)):
            if i == 0:
                best = neighImprove[i].val
                bestI = neighImprove[i].index[0]
            # elif (best<neighImprove[i].val - self.threshold and minimize == True) or (best>neighImprove[i].val + self.threshold and minimize == False): # seek for max S
            elif (best<neighImprove[i].val + self.threshold and minimize == True) or (best>neighImprove[i].val - self.threshold and minimize == False): # seek for max S
                best = neighImprove[i].val
                bestI = neighImprove[i].index[0]

        bestIlist = []
        for i in range(len(neighImprove)):
            if abs(best - neighImprove[i].val) < self.threshold:
                candI = neighImprove[i].index
                if candI not in bestIlist:
                    bestIlist.append(candI)

        bestI = random.choice(bestIlist)
        if type(bestI) is int:
            # make a consistent interface
            bestI = [bestI]

        #return True, bestI, evalCount + self.dim * (self.dim-1)
        return True, bestI, evalCount 

    def genImproveSC(self,minimize):
        """
        generate the index of best neigh according to {S_p(X)-2/N \Sigma_{i=1}^{N}C_{ip}(X)} only (surrogate of fitness)
        """
        # check improving move 
        self.improveA = []
        for i in range(self.dim):
            # if (minimize == True and self.SC[i] > self.threshold) or (minimize == False and self.SC[i]< - self.threshold):
            if (minimize == True and self.SC[i] > - self.threshold) or (minimize == False and self.SC[i] < self.threshold):
                self.improveA.append(i)

    def genMeanBest(self,minimize):
        """
        generate the index of best neigh according to {S_p(X)-2/N \Sigma_{i=1}^{N}C_{ip}(X)} only (surrogate of fitness)
        """
        # check improving move 
        improve = False
        self.improveA = []
        for i in range(self.dim):
            # if (minimize == True and self.SC[i] > self.threshold) or (minimize == False and self.SC[i]< - self.threshold):
            if (minimize == True and self.SC[i] > - self.threshold) or (minimize == False and self.SC[i] < self.threshold):
                self.improveA.append(i)
                improve = True

        if improve == False:
            return False, None

        random.shuffle(self.improveA)

        # find the best value
        for i in self.improveA:
            if i == self.improveA[0]:
                best = self.SC[i]
                bestI = i
            # elif ( best<self.SC[i] - self.threshold and minimize == True ) or ( best>self.SC[i] + self.threshold and minimize == False ): # seek for max S
            elif ( best<self.SC[i] + self.threshold and minimize == True ) or ( best>self.SC[i] - self.threshold and minimize == False ): # seek for max S
                best = self.SC[i]
                bestI = i
        return True, bestI

    def genMeanNext(self,minimize):
        """
        generate the index of next improving neigh according to {S_p(X)-2/N \Sigma_{i=1}^{N}C_{ip}(X)} only (surrogate of fitness)
        """
        # check improving move 
        improve = False
        self.improveA = []
        for i in range(self.dim):
            # if (minimize == True and self.SC[i] > self.threshold) or (minimize == False and self.SC[i] < - self.threshold):
            if (minimize == True and self.SC[i] > - self.threshold) or (minimize == False and self.SC[i] < self.threshold):
                self.improveA.append(i)
                improve = True

        if improve == False:
            return False, None

        bestI = random.choice(self.improveA)
        return True, bestI

    def updateMeanBest(self, p, minimize):
        self.improveA.remove(p)
        if p in self.Inter:
            for i in self.Inter[p].arr:
                # if (minimize == True and self.SC[i] > self.threshold) or (minimize == False and self.SC[i]< - self.threshold):
                if (minimize == True and self.SC[i] > - self.threshold) or (minimize == False and self.SC[i] < self.threshold):
                    if i not in self.improveA:
                        self.improveA.append(i)
                elif i in self.improveA:
                    self.improveA.remove(i)

        if not self.improveA:
            return False, None

        random.shuffle(self.improveA)

        for i in self.improveA:
            if i == self.improveA[0]:
                best = self.SC[i]
                bestI = i
            # elif (best<self.SC[i] - self.threshold and minimize == True) or (best>self.SC[i]+self.threshold and minimize == False): # seek for max S
            elif (best<self.SC[i] + self.threshold and minimize == True) or (best>self.SC[i] - self.threshold and minimize == False): # seek for max S
                best = self.SC[i]
                bestI = i
                    
        return True, bestI

    def updateMeanNext(self, p, minimize):
        self.improveA.remove(p)
        if p in self.Inter:
            for i in self.Inter[p].arr:
                # if (minimize == True and self.SC[i] > self.threshold) or (minimize == False and self.SC[i]< - self.threshold):
                if (minimize == True and self.SC[i] > - self.threshold) or (minimize == False and self.SC[i]< self.threshold):
                    if i not in self.improveA:
                        self.improveA.append(i)
                elif i in self.improveA:
                    self.improveA.remove(i)

        if not self.improveA:
            return False, None

        bestI = random.choice(self.improveA)
        return True, bestI

    def binCount(self, arr, bit):
        """
        count the one bit of union self.model.WA[i].arr and bit
        """
        s = 0
        for i in arr:
            if bit[i] == '1':
                s = s + 1
        return s
    
    def genComb(self,N):
        """ 
        Generate C_k^0 sequence, index are stored, because they are more general, Implemented in an *incremental* fashion.
        """
        if N in self.lookup.keys(): # the key exists before
            return self.lookup[N]
        else : 
            comb =  []
            c = 0
            for i in range(N):
                for j in [ k for k in range(N) if k > i]:
                   arr = np.zeros(2)
                   arr[0] = i
                   arr[1] = j
                   comb.append(arr)
                   c = c + 1    
            self.lookup[N] = comb
            return comb
