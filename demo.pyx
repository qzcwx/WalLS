import nkLandscape as nk
import nkqLandscape as nkq
import WalshAnalysis as wal
import geneticAlgorithm as ga
import AutoCorrelation as ac
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
import argparse

def main():
    """ command line options """
    # Usage: python demo.py [ComputeMethod] [NameOfProblem] [NameOfAlgorithm] [fit/mean/std] [overwrite] [I] [PopSize] [N] [K] [Q]
    parser = argparse.ArgumentParser(description='Walsh Local Search')

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
    parser.add_argument('-d',
                        action="store",
                        help="Radius of Hamming Spheres/Balls",
                        dest="radius",
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
    parser.add_argument('-v',
                        action="store",
                        help="CV ratio",
                        dest="v",
                        type=float,
                        default=0,
                        )
    parser.add_argument('-e',
                        action="store",
                        help="factor for MaxFit",
                        dest="e",
                        type=int,
                        default=1,
                        )
    parser.add_argument('-l',
                        action="store",
                        help="walk length",
                        dest="l",
                        type=int,
                        default=10,
                        )

    
    opt = parser.parse_args()

    if opt.v != 0:
        opt.c = int(opt.v*opt.n)

    #tl.checkParam(sys.argv)
    opt.s = opt.popSize
    
    nameOfDir = './result/'
    runtimeDir = './runtime/'
    waltimeDir = './walshtime/'
    traceDir = './trace/'
    prefixNK = './benchmark/NK/'
    prefixNKQ = './benchmark/NKQ/'

    # print 'radius', opt.radius
    
    # print opt.rseed
    random.seed(opt.rseed)

    # maxFit = opt.e * opt.n
    # runs = 20
    
    # runs = 10
    # maxFit = 100000                      # 100 k, submit for running experiment
    # maxFit = 1000000                      # 1 million
    
    # maxFit = 2000
    
    runs = 1
    maxFit = 100
    
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

            start = time.time()
            if opt.algoName.find('GA') != -1:
                res.append(algo.run(crossoverR, mutationR, opt.fitName, minimize = False))
            elif opt.algoName == 'LS':
                res.append(algo.run(opt.fitName, minimize = False, restart = False, beamWidth=opt.w))
            elif opt.algoName == 'rLS':
                res.append(algo.run(opt.fitName, minimize = False, restart = True, beamWidth=opt.w))
            elif opt.algoName.find('CHC') != -1:
                res.append(algo.run(model.compFit, maxFit, opt.popSize, opt.n, D, DR, M, opt.fitName, minimize = False))
            tAll[i] = time.time() - start

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
            model = nk.NKLandscape(opt.n,opt.k,opt.c,prefixNK+opt.probName+'-N'+str(opt.n)+'-K'+str(opt.k)+'-C'+str(opt.c)+'-I'+str(opt.inst))
        elif opt.probName == 'NKQ':
            model = nkq.NKQLandcape(opt.n, opt.k, opt.c, opt.q, opt.t, prefixNKQ+opt.probName+'-N'+str(opt.n)+'-K'+str(opt.k)+'-C'+str(opt.c)+'-I'+str(opt.inst)+'-Q'+str(opt.q)+'-T'+str(t))
        elif opt.probName == 'NonNK':
            model = nk.NonNKLandscape(opt.n,opt.k,opt.c,prefixNK+opt.probName+'-N'+str(opt.n)+'-K'+str(opt.k)+'-C'+str(opt.c)+'-I'+str(opt.inst))
            
        # Walsh analysis
        if opt.compMeth == 'walWalk' or opt.compMeth == 'walRest' or opt.compMeth == 'supm' or opt.compMeth == 'bitImp' or opt.compMeth == 'walSearch' or opt.compMeth == 'checkOptWal' or opt.compMeth == 'checkHyper' or opt.compMeth == 'checkHyperRank' or opt.compMeth == 'hyperSearch' or opt.compMeth == 'hyperSqSearch' or opt.compMeth == 'hyperWalSearch' or opt.compMeth == 'walWalkNext' or opt.compMeth == 'walRestNext' or opt.compMeth == 'BeamWalkNext' or opt.compMeth=='BeamWalk' or opt.compMeth == 'walRestFlip' or opt.compMeth == 'walRestU' or opt.compMeth == 'walWalkULen' or opt.compMeth == 'walRestUDist' or opt.compMeth == 'walWalkU' or opt.compMeth == 'walWalkUDist' or opt.compMeth == 'walRestNextU' or opt.compMeth == 'walTLONextU' or opt.compMeth == 'walHS' or opt.compMeth == 'walBS':
            start = time.time()
            w = model.WalshCofLinearLinklist()
            walTime = time.time() - start
            # print 'walsh trans\n', "%.4g" %(walTime)
            

            
            start = time.time()
            if opt.compMeth == 'checkHyper' or opt.compMeth == 'checkHyperRank' or opt.compMeth == 'hyperSearch':
                model.genHyperVote()
            elif opt.compMeth == 'hyperSqSearch':
                model.genHyperSqVote()
            elif opt.compMeth == 'hyperWalSearch':
                model.genHyperWalVote()
            hyperTime = time.time() - start

            # count the number of interative bits
            # model.countInterBits()

            """ store runtime to files """
            if opt.probName == 'NKQ' :
                nameOfF = waltimeDir+opt.probName+'-'+opt.algoName+'-F'+opt.fitName+'-M'+opt.compMeth+'-I'+str(opt.inst)+'-S'+str(opt.s)+'-N'+str(opt.n)+'-C'+str(opt.c)+'-K'+str(opt.k)+'-Q'+str(opt.q)+'-T'+str(t)+'.txt'
            elif opt.probName == 'NK' or opt.probName == 'NonNK':
                nameOfF = waltimeDir+opt.probName+'-'+opt.algoName+'-F'+opt.fitName+'-M'+opt.compMeth+'-I'+str(opt.inst)+'-S'+str(opt.s)+'-N'+str(opt.n)+'-C'+str(opt.c)+'-K'+str(opt.k)+'.txt'

            f = open(nameOfF, 'w')
            print >>f,"%g\t%g" % (walTime,hyperTime)
            f.close()
            
        # release memory of list f 
        model.f = []
        # print 'release memory of list f'
        
        if opt.compMeth == 'countFreq':
            freq = model.countFreqInFunc()
            # print sum(freq)

        # print 'enumerating all solutions'
        # if opt.fitName == 'fit':
        #     bit,fit = tl.compFit(model)
        # elif opt.fitName == 'mean':
        #     bit, fit = tl.compMean(model)
        # print max(fit)
        # bitF,fitF = tl.compFit(model)
        # for i in zip(bitF,fitF):
        #     print i[0],i[1]


        # bitA,fitA = tl.compMean(model)
        # for i in zip(bitF,fitF, bitA, fitA):
        #     print i[0],'%g\t%g' %(i[1],i[3])

        # a = sorted(zip(bit,fit), key=lambda a_entry: a_entry[1])
        # print 'opti\n', a[0][0], a[0][1]
        # print

        # c = 0
        # for i in (zip(bit,fit)):
        #     print c,'\t', i[0], '%.4f' %(i[1])
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
        # print
        for i in range(runs):
            # print 'run', i 
            start = time.time()
            if opt.algoName.find('GA') != -1:
                res.append(algo.run(crossoverR, mutationR, opt.fitName))
            elif opt.algoName == 'LS':
                res.append(algo.run(opt.fitName, minimize = False, restart = False,compM = opt.compMeth, beamWidth=opt.w ))
            elif opt.algoName == 'rLS':
                res.append(algo.run(opt.fitName, minimize = False, restart = True,compM = opt.compMeth, beamWidth=opt.w, walklen = opt.l, radius = opt.radius))
            elif opt.algoName.find('CHC') != -1:
                res.append(algo.run(model.compFit, maxFit,  opt.popSize, opt.n, D, DR, M, opt.fitName))
            tAll[i] = time.time() - start
            # print
            
        model.destFunc()

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
            nameOfF = nameOfDir+opt.probName+'-'+opt.algoName+'-F'+opt.fitName+'-M'+opt.compMeth+'-I'+str(opt.inst)+'-S'+str(opt.s)+'-W'+str(opt.w)+'-N'+str(opt.n)+'-K'+str(opt.k)+'-C'+str(opt.c)+'-Q'+str(opt.q)+'-T'+str(t)+'-E'+str(opt.e)+'-L'+str(opt.l)+'.txt'
        elif opt.probName == 'NK' or opt.probName == 'NonNK':
            nameOfF = nameOfDir+opt.probName+'-'+opt.algoName+'-F'+opt.fitName+'-M'+opt.compMeth+'-I'+str(opt.inst)+'-S'+str(opt.s)+'-W'+str(opt.w)+'-N'+str(opt.n)+'-K'+str(opt.k)+'-C'+str(opt.c)+'-E'+str(opt.e)+'-L'+str(opt.l)+'-D'+str(opt.radius)+'.txt'

        # print nameOfF     # print result
    #    """ print the mean over multiple runs """
    #    r = np.zeros(runs)
    #    for i in range(runs):
    #        r[i] = res[i]['sol']
    #    print np.mean(r)

        f = open(nameOfF, 'w')
        if opt.compMeth == 'countFreq':
            for i in freq:
                print >>f,"%d" %(i)
        else:
            for i in range(runs):
                if opt.fitName != 'fit':
                    print >>f,"%g\t%g\t%g" % (res[i]['sol'], res[i]['fitG'], res[i]['nEvals'])
                else:
                    # print "%g\t%g" % (res[i]['sol'], res[i]['nEvals'])
                    print >>f,"%g\t%g" % (res[i]['sol'], res[i]['nEvals'])
        f.close()

        # print nameOfF
        if opt.compMeth == 'countFreq':
            exit()
        
        if (opt.compMeth != 'bf'  and opt.compMeth!='partEval') and (opt.fitName == 'switchF' or opt.fitName == 'switchAvg' or opt.fitName == 'combF' or opt.fitName == 'combAvg'):
            """ store trace to files: 1. the number of descent steps """
            if opt.probName == 'NKQ':
                nameOfF = traceDir+opt.probName+'-'+opt.algoName+'-F'+opt.fitName+'-M'+opt.compMeth+'-I'+str(opt.inst)+'-S'+str(opt.s)+'-W'+str(opt.w)+'-N'+str(opt.n)+'-K'+str(opt.k)+'-C'+str(opt.c)+'-Q'+str(opt.q)+'-T'+str(t)+'-E'+str(opt.e)+'-L'+str(opt.l)+'-D'+str(opt.radius)+'.txt'
            elif opt.probName == 'NK' or opt.probName == 'NonNK':
                nameOfF = traceDir+opt.probName+'-'+opt.algoName+'-F'+opt.fitName+'-M'+opt.compMeth+'-I'+str(opt.inst)+'-S'+str(opt.s)+'-W'+str(opt.w)+'-N'+str(opt.n)+'-K'+str(opt.k)+'-C'+str(opt.c)+'-E'+str(opt.e)+'-L'+str(opt.l)+'-D'+str(opt.radius)+'.txt'
            f = open(nameOfF, 'w')
            for i in range(runs):
                print >>f,"%g\t%g\t%g" % (res[i]['platC'], res[i]['restC'], res[i]['updateC'])
            f.close()
        else:
            if opt.probName == 'NKQ':
                nameOfF = traceDir+opt.probName+'-'+opt.algoName+'-F'+opt.fitName+'-M'+opt.compMeth+'-I'+str(opt.inst)+'-S'+str(opt.s)+'-W'+str(opt.w)+'-N'+str(opt.n)+'-K'+str(opt.k)+'-C'+str(opt.c)+'-Q'+str(opt.q)+'-T'+str(t)+'-E'+str(opt.e)+'-L'+str(opt.l)+'-D'+str(opt.radius)+'.txt'
            elif opt.probName == 'NK' or opt.probName == 'NonNK':
                nameOfF = traceDir+opt.probName+'-'+opt.algoName+'-F'+opt.fitName+'-M'+opt.compMeth+'-I'+str(opt.inst)+'-S'+str(opt.s)+'-W'+str(opt.w)+'-N'+str(opt.n)+'-K'+str(opt.k)+'-C'+str(opt.c)+'-E'+str(opt.e)+'-L'+str(opt.l)+'-D'+str(opt.radius)+'.txt'

            if 'TLO' not in opt.compMeth:
                f = open(nameOfF, 'w')
                # print >>f,"initC\tupdateC\t"
                print >>f,"initC\tupdateC\tbackC\t"
                for i in range(runs):
                    print >>f,"%g\t%g\t%g\t" % (res[i]['initC'], res[i]['updateC'], res[i]['backC'])
                f.close()

        if opt.compMeth != 'bf' and opt.compMeth != 'partEval':
            """ store trace to files: the number of descent steps """

            for i in range(runs):
                if opt.probName == 'NKQ':
                    nameOfF = traceDir+opt.probName+'-'+opt.algoName+'-F'+opt.fitName+'-M'+opt.compMeth+'-I'+str(opt.inst)+'-S'+str(opt.s)+'-W'+str(opt.w)+'-N'+str(opt.n)+'-K'+str(opt.k)+'-C'+str(opt.c)+'-Q'+str(opt.q)+'-T'+str(t)+'-R'+str(i)+'-E'+str(opt.e)+'-L'+str(opt.l)+'-D'+str(opt.radius)+'.txt'
                elif opt.probName == 'NK' or opt.probName == 'NonNK':
                    nameOfF = traceDir+opt.probName+'-'+opt.algoName+'-F'+opt.fitName+'-M'+opt.compMeth+'-I'+str(opt.inst)+'-S'+str(opt.s)+'-W'+str(opt.w)+'-N'+str(opt.n)+'-K'+str(opt.k)+'-C'+str(opt.c)+'-R'+str(i)+'-E'+str(opt.e)+'-L'+str(opt.l)+'-D'+str(opt.radius)+'.txt'

                if 'TLO' not in opt.compMeth :
                    f = open(nameOfF, 'w')
                    print >>f,"step\tbest\t"
                    if (opt.fitName == 'fit'):
                        for j in zip(res[i]['traceEval'], res[i]['traceFit']):
                            # print >>f,"%g\t%g" % (j[0], j[1])
                            print >>f,"%g\t%g" % (j[0], j[1])
                    else :
                        for j in zip(res[i]['traceEval'], res[i]['traceFit'],res[i]['traceFitG']):
                            print >>f,"%g\t%g\t%g" % (j[0], j[1], j[2])
                    f.close()

        """ store runtime to files """
        if opt.probName == 'NKQ':
            nameOfF = runtimeDir+opt.probName+'-'+opt.algoName+'-F'+opt.fitName+'-M'+opt.compMeth+'-I'+str(opt.inst)+'-S'+str(opt.s)+'-W'+str(opt.w)+'-N'+str(opt.n)+'-K'+str(opt.k)+'-C'+str(opt.c)+'-Q'+str(opt.q)+'-T'+str(t)+'-E'+str(opt.e)+'-L'+str(opt.l)+'-D'+str(opt.radius)+'.txt'
        elif opt.probName == 'NK' or opt.probName == 'NonNK':
            nameOfF = runtimeDir+opt.probName+'-'+opt.algoName+'-F'+opt.fitName+'-M'+opt.compMeth+'-I'+str(opt.inst)+'-S'+str(opt.s)+'-W'+str(opt.w)+'-N'+str(opt.n)+'-K'+str(opt.k)+'-C'+str(opt.c)+'-E'+str(opt.e)+'-L'+str(opt.l)+'-D'+str(opt.radius)+'.txt'
            
        f = open(nameOfF, 'w')
        if 'TLO' in opt.compMeth:
            print >>f,"All\t\tinit\t\tdesc\t\tupdate"
            # print "All\t\tinit\t\tdesc\t\tupdate"
            for i in range(runs):
                print >>f,"%0.4e\t%0.4e\t%0.4e\t%0.4e" % (tAll[i], res[i]['init'],res[i]['descT'], res[i]['updateT']) 
                # print "%0.4e\t%0.4e\t%0.4e\t%0.4e" % (tAll[i], res[i]['init'],res[i]['descT'], res[i]['updateT']) 
        elif opt.compMeth != 'bf' and opt.compMeth != 'partEval':
            # tracking for walsh-based approach
            print >>f,"All\t\tinit\t\tdesc\t\tpert\t\tupdate\t\tupdatePert\t"
            for i in range(runs):
                print >>f,"%0.4e\t%0.4e\t%0.4e\t%0.4e\t%0.4e\t%0.4e" % (tAll[i], res[i]['init'],res[i]['descT'], res[i]['pertT'], res[i]['updateT'], res[i]['updatePertT'])
        else:
            # for the brute force approach
            print >>f,"All\tupdate"
            for i in range(runs):
                print >>f,"%0.4e\t%0.4e" % (tAll[i], res[i]['updateT'])


        f.close()

        if opt.compMeth == 'walRestFlip':
            """ store trace of bit-flips to files """
            for i in range(runs):
                if opt.probName == 'NKQ' :
                    nameOfF = runtimeDir+opt.probName+'-'+opt.algoName+'-F'+opt.fitName+'-M'+opt.compMeth+'-I'+str(opt.inst)+'-S'+str(opt.s)+'-W'+str(opt.w)+'-N'+str(opt.n)+'-K'+str(opt.k)+'-C'+str(opt.c)+'-Q'+str(opt.q)+'-T'+str(t)+'-E'+str(opt.e)+'-R'+str(i)+'-L'+str(opt.l)+'-D'+str(opt.radius)+'.txt'
                elif opt.probName == 'NK' or opt.probName == 'NonNK':
                    nameOfF = runtimeDir+opt.probName+'-'+opt.algoName+'-F'+opt.fitName+'-M'+opt.compMeth+'-I'+str(opt.inst)+'-S'+str(opt.s)+'-W'+str(opt.w)+'-N'+str(opt.n)+'-K'+str(opt.k)+'-C'+str(opt.c)+'-E'+str(opt.e)+'-R'+str(i)+'-L'+str(opt.l)+'-D'+str(opt.radius)+'.txt'
                f = open(nameOfF, 'w')
                for j in res[i]['traceFlip']:
                    print >>f, "%d" %(j)
                f.close()

