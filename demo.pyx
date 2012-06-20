import nkLandscape as nk
import nkqLandscape as nkq
import WalshAnalysis as wal
import geneticAlgorithm as ga
import AutoCorrelation as ac
#import matplotlib.pyplot as plt
import LocalOptima as lo
import CHC as chc
import MAXSAT as mx
import LocalSearch as ls
import LocalOptima as lo
import tool as tl
import os
import numpy as np
import random
import math
import time
import pdb
import sys
import batch
import argparse

def main():
    """ command line options """
    # Usage: python demo.py [ComputeMethod] [NameOfProblem] [NameOfAlgorithm] [fit/mean/std] [overwrite] [I] [PopSize] [N] [K] [Q]
    parser = argparse.ArgumentParser(description='Walsh Local Search')
    #print argv
    #parser.parse_args(argv)

    parser.add_argument('-m', 
                        action="store", 
                        help="Computational Method",
                        dest="compMeth",
                        )
    parser.add_argument('-p', 
                        action="store", 
                        help="Problem",
                        dest="probName", 
                        )
    parser.add_argument('-a', 
                        action="store", 
                        help="Algorithm",
                        dest="algoName",
                        )
    parser.add_argument('-f', 
                        action="store", 
                        help="Evaluation Function",
                        dest="fitName",
                        )
    parser.add_argument('-i', 
                        action="store", 
                        help="Instance ID",
                        dest="inst",
                        type=int,
                        )
    parser.add_argument('-s', 
                        action="store", 
                        help="Population size",
                        dest="popSize",
                        default=1,
                        type=int,
                        )
    parser.add_argument('-n', 
                        action="store", 
                        help="Dimension",
                        dest="n",
                        type=int,
                        )
    parser.add_argument('-k', 
                        action="store", 
                        help="K",
                        dest="k",
                        type=int,
                        )
    parser.add_argument('-q', 
                        action="store", 
                        help="Q",
                        dest="q",
                        default=0,
                        type=int,
                        )
    parser.add_argument('-r', 
                        action="store", 
                        help="Random Seed",
                        dest="rseed",
                        default=0,
                        type=int,
                        )
    parser.add_argument('-w', 
                        action="store", 
                        help="Width of Beam",
                        dest="w",
                        default=1,
                        type=int,
                        )
    parser.add_argument('-t',
                        action="store",
                        help="number of true assignment in one subfunction",
                        dest="t",
                        default=-1,
                        type=int,
                        )
    parser.add_argument('-c',
                        action="store",
                        help="number of subfunctions",
                        dest="c",
                        type=int,
                        )

    opt = parser.parse_args()

    #tl.checkParam(sys.argv)
    opt.s = opt.popSize

    #rseed = 0
    nameOfDir = './result/'
    runtimeDir = './runtime/'
    waltimeDir = './walshtime/'
    traceDir = './trace/'
    prefixNK = './benchmark/NK/'
    prefixNKQ = './benchmark/NKQ/'

    random.seed(opt.rseed)

    # maxFit = 100
    maxFit = 1000 * opt.n

    runs = 30
    # runs = 1

    t = opt.t
    
    crossoverR = 0.8 # typically in (0.6, 0.9)
    mutationR = 1.0/float(opt.n) # typically between 1/popSize and 1/dim
    # for CHC 
    D = opt.n/4.0
    DR = 0.35
    M = 1


    if opt.algoName.find('LS') != -1:
        popSize = 1

    if opt.probName == 'SAT':
        """ with SAT, we are forced to set n to 100 """
        """ 
        TODO : 
            need to perform multiple runs for each instance 
        """

        model = mx.MAXSAT()
        res = []

        model.setInstance(opt.inst)
        print 'Instance', opt.inst

        if opt.algoName.find('LS') != -1:
            algo = ls.LocalSearch(model.compFit, maxFit, opt.n)
        elif opt.algoName.find('GA') != -1:
            algo = ga.GeneticAlgorithm( model.compFit, maxFit, opt.popSize, opt.n )
        elif opt.algoName.find('CHC') != -1:
            algo = chc.CHC()

        tAll = np.zeros(runs)

        for i in range(runs):
            start = os.times()[0]
            if opt.algoName.find('GA') != -1:
                res.append(algo.run(crossoverR, mutationR, opt.fitName, minimize = False))
            elif opt.algoName == 'LS':
                res.append(algo.run(opt.fitName, minimize = False, restart = False, beamWidth=opt.w))
            elif opt.algoName == 'rLS':
                res.append(algo.run(opt.fitName, minimize = False, restart = True, beamWidth=opt.w))
            elif opt.algoName.find('CHC') != -1:
                res.append(algo.run(model.compFit, maxFit, opt.popSize, opt.n, D, DR, M, opt.fitName, minimize = False))
            tAll[i] = os.times()[0] - start

        if opt.probName == 'SAT':
            nameOfF = nameOfDir+opt.probName+'-'+opt.algoName+'-F'+opt.fitName+'-M'+opt.compMeth+'-I'+str(opt.inst)+'-S'+str(opt.s)+'-N'+str(opt.n)+'.txt'
        f = open(nameOfF, 'w')
        for i in range(len(res)):
            if opt.fitName != 'fit':
                print >>f,"%g\t%g\t%g" % (res[i]['sol'], res[i]['fitG'], res[i]['nEvals'])
            else:
                print >>f,"%g\t%g" % (res[i]['sol'], res[i]['nEvals'])
        f.close()

        """ store runtime to files """
        if opt.probName == 'SAT':
            nameOfF = runtimeDir+opt.probName+'-'+opt.algoName+'-F'+opt.fitName+'-M'+opt.compMeth+'-I'+str(opt.inst)+'-S'+str(opt.s)+'-N'+str(opt.n)+'.txt'

        f = open(nameOfF, 'w')
        for i in range(runs):
            print >>f,"%g" % (tAll[i])
        f.close()

    else:
        res = []

        if opt.probName == 'NK':
            model = nk.NKLandscape(opt.n,opt.k,opt.c,prefixNK+'NK-N'+str(opt.n)+'-K'+str(opt.k)+'-C'+str(opt.c)+'-I'+str(opt.inst))
            #model = nk.NKLandscape(n,k)
        elif opt.probName == 'NKQ':
            model = nkq.NKQLandcape(opt.n, opt.k, opt.c, opt.q, opt.t, prefixNKQ+'NKQ-N'+str(opt.n)+'-K'+str(opt.k)+'-C'+str(opt.c)+'-I'+str(opt.inst)+'-Q'+str(opt.q)+'-T'+str(t))
            #model = nkq.NKQLandcape(n, k, q)

        if opt.compMeth == 'walWalk' or opt.compMeth == 'walRest' or opt.compMeth == 'supm' or opt.compMeth == 'bitImp' or opt.compMeth == 'walSearch' or opt.compMeth == 'checkOptWal' or opt.compMeth == 'checkHyper' or opt.compMeth == 'checkHyperRank' or opt.compMeth == 'hyperSearch' or opt.compMeth == 'hyperSqSearch' or opt.compMeth == 'hyperWalSearch' or opt.compMeth == 'walWalkNext' or opt.compMeth == 'walRestNext' or opt.compMeth == 'BeamWalkNext' or opt.compMeth=='BeamWalk':
            start = os.times()[0]
            # Walsh analysis
            w = model.WalshCofLinearLinklist()
            walTime = os.times()[0] - start

            start = os.times()[0]
            if opt.compMeth == 'checkHyper' or opt.compMeth == 'checkHyperRank' or opt.compMeth == 'hyperSearch':
                model.genHyperVote()
            elif opt.compMeth == 'hyperSqSearch':
                model.genHyperSqVote()
            elif opt.compMeth == 'hyperWalSearch':
                model.genHyperWalVote()
            hyperTime = os.times()[0] - start

            # count the number of interative bits
            # model.countInterBits()

            """ store runtime to files """
            if opt.probName == 'NKQ':
                nameOfF = waltimeDir+opt.probName+'-'+opt.algoName+'-F'+opt.fitName+'-M'+opt.compMeth+'-I'+str(opt.inst)+'-S'+str(opt.s)+'-N'+str(opt.n)+'-C'+str(opt.c)+'-K'+str(opt.k)+'-Q'+str(opt.q)+'-T'+str(t)+'.txt'
            elif opt.probName == 'NK':
                nameOfF = waltimeDir+opt.probName+'-'+opt.algoName+'-F'+opt.fitName+'-M'+opt.compMeth+'-I'+str(opt.inst)+'-S'+str(opt.s)+'-N'+str(opt.n)+'-C'+str(opt.c)+'-K'+str(opt.k)+'.txt'

            f = open(nameOfF, 'w')
            print >>f,"%g\t%g" % (walTime,hyperTime) 
            f.close()

        # if opt.fitName == 'fit':
        #     bit,fit = tl.compFit(model)
        # elif opt.fitName == 'mean':
        #     bit, fit = tl.compMean(model)

        # bitF,fitF = tl.compFit(model)
        # bitA,fitA = tl.compMean(model)
        # for i in zip(bitF,fitF, bitA, fitA):
        #     print i[0],'%g\t%g' %(i[1],i[3])

        
                
        # a = sorted(zip(bit,fit), key=lambda a_entry: a_entry[1]) 
        # print 'opti\n', a[0][0], a[0][1]
        # print
        
        # c = 0
        # for i in (zip(bit,fit)):
        #     print c,'\t', i[0], '%.2f' %(i[1])
        #     c = c + 1


        # for i in zip(bit,fit):
        #     print i[0],'%.3f' %(i[1])

        # print 'bit',bit
        # print 'fit',fit
        # print 'mean',np.mean(fit)
        # print 'w', w

        # numOpt = lo.localOpt(bit, fit)
        # print numOpt

        # lo.plateauCount(bit, fit, opt)

        if opt.algoName.find('LS') != -1:
            algo = ls.LocalSearch(model, maxFit, opt.n)
        elif opt.algoName.find('GA') != -1:
            algo = ga.GeneticAlgorithm( model.compFit, maxFit, opt.popSize, opt.n )
        elif opt.algoName.find('CHC') != -1:
            algo = chc.CHC()

        tAll = np.zeros(runs)
        for i in range(runs):
            start = os.times()[0]
            if opt.algoName.find('GA') != -1:
                res.append(algo.run(crossoverR, mutationR, opt.fitName))
            elif opt.algoName == 'LS':
                res.append(algo.run(opt.fitName, minimize = False, restart = False,compM = opt.compMeth, beamWidth=opt.w ))
            elif opt.algoName == 'rLS':
                res.append(algo.run(opt.fitName, minimize = False, restart = True,compM = opt.compMeth, beamWidth=opt.w))
            elif opt.algoName.find('CHC') != -1:
                res.append(algo.run(model.compFit, maxFit,  opt.popSize, opt.n, D, DR, M, opt.fitName))
            tAll[i] = os.times()[0] - start
            
    #    trace = res[0]['trace']
    #    for i in trace:
    ##        print 'Eval', i.fitEval, 'fit', i.fit
    #        print 'Eval', i.fitEval, 'fit', i.fit, 'fitG', i.fitG
    #
    #    plt.plot([i.fitEval for i in trace],[i.fit for i in trace],'.-')
    #    plt.plot([i.fitEval for i in trace],[i.fitG for i in trace],'.-')
    #    plt.show()

        """ store results to files """
        if opt.probName == 'NKQ':
            nameOfF = nameOfDir+opt.probName+'-'+opt.algoName+'-F'+opt.fitName+'-M'+opt.compMeth+'-I'+str(opt.inst)+'-S'+str(opt.s)+'-W'+str(opt.w)+'-N'+str(opt.n)+'-K'+str(opt.k)+'-C'+str(opt.c)+'-Q'+str(opt.q)+'-T'+str(t)+'.txt'
        elif opt.probName == 'NK':
            nameOfF = nameOfDir+opt.probName+'-'+opt.algoName+'-F'+opt.fitName+'-M'+opt.compMeth+'-I'+str(opt.inst)+'-S'+str(opt.s)+'-W'+str(opt.w)+'-N'+str(opt.n)+'-K'+str(opt.k)+'-C'+str(opt.c)+'.txt'

    #    """ print the mean over multiple runs """
    #    r = np.zeros(runs)
    #    for i in range(runs):
    #        r[i] = res[i]['sol']
    #    print np.mean(r)
        f = open(nameOfF, 'w')
        for i in range(runs):
            if opt.fitName != 'fit':
                print >>f,"%g\t%g\t%g" % (res[i]['sol'], res[i]['fitG'], res[i]['nEvals'])
            else:
                print >>f,"%g\t%g" % (res[i]['sol'], res[i]['nEvals'])
        f.close()
        print nameOfF

        """ store trace to files """
        if opt.probName == 'NKQ':
            nameOfF = traceDir+opt.probName+'-'+opt.algoName+'-F'+opt.fitName+'-M'+opt.compMeth+'-I'+str(opt.inst)+'-S'+str(opt.s)+'-W'+str(opt.w)+'-N'+str(opt.n)+'-K'+str(opt.k)+'-C'+str(opt.c)+'-Q'+str(opt.q)+'-T'+str(t)+'.txt'
        elif opt.probName == 'NK':
            nameOfF = traceDir+opt.probName+'-'+opt.algoName+'-F'+opt.fitName+'-M'+opt.compMeth+'-I'+str(opt.inst)+'-S'+str(opt.s)+'-W'+str(opt.w)+'-N'+str(opt.n)+'-K'+str(opt.k)+'-C'+str(opt.c)+'.txt'
        f = open(nameOfF, 'w')
        for i in range(runs):
              print >>f,"%g\t%g" % (res[i]['initC'], res[i]['updateC'])
        f.close()

        """ store runtime to files """
        if opt.probName == 'NKQ':
            nameOfF = runtimeDir+opt.probName+'-'+opt.algoName+'-F'+opt.fitName+'-M'+opt.compMeth+'-I'+str(opt.inst)+'-S'+str(opt.s)+'-W'+str(opt.w)+'-N'+str(opt.n)+'-K'+str(opt.k)+'-C'+str(opt.c)+'-Q'+str(opt.q)+'-T'+str(t)+'.txt'
        elif opt.probName == 'NK':
            nameOfF = runtimeDir+opt.probName+'-'+opt.algoName+'-F'+opt.fitName+'-M'+opt.compMeth+'-I'+str(opt.inst)+'-S'+str(opt.s)+'-W'+str(opt.w)+'-N'+str(opt.n)+'-K'+str(opt.k)+'-C'+str(opt.c)+'.txt'

        f = open(nameOfF, 'w')
        print >>f,"All\t\tinit\t\tdesc\t\tpert\t\tupdate\t\tupdatePert\t"
        for i in range(runs):
            print >>f,"%0.2e\t%0.2e\t%0.2e\t%0.2e\t%0.2e\t%0.2e" % (tAll[i], res[i]['init'],res[i]['descT'], res[i]['pertT'], res[i]['updateT'], res[i]['updatePertT'])
        f.close()

        print
