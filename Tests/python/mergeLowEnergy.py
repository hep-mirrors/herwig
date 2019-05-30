#! /usr/bin/env python
import yoda,glob,math,optparse
op = optparse.OptionParser(usage=__doc__)

opts, args = op.parse_args()

if(len(args)!=1) :
    print 'Must be one and only 1 name'
    quit()


for runType in ["NonPerturbative","Perturbative"]:
    outhistos={}
    for fileName in glob.glob("Rivet-LowEnergy-EE-%s-*.yoda" % runType):
        energy = float(fileName.split("-")[-1].strip(".yoda"))
        energyMeV = energy*1000.
        aos = yoda.read(fileName)
        for hpath,histo in aos.iteritems():
            if("/_" in hpath or "TMP" in hpath) : continue
            if(type(histo)==yoda.core.Histo1D) :
                outhistos[hpath] = histo
                continue
            # create histo if it doesn't exist
            elif(hpath not in outhistos) :
                outhistos[hpath] = yoda.core.Scatter2D(histo.path,
                                                       histo.title)
            matched = False
            for i in range(0,aos[hpath].numPoints) :
                x = aos[hpath].points[i].x
                delta=1e-5
                if("KLOE_2009_I797438" in hpath or "KLOE_2005_I655225" in hpath or
                   "FENICE_1994_I377833" in hpath):
                   x=math.sqrt(x)
                   delta=1e-3
                if(abs(x-energy)<1e-3*delta or abs(x-energyMeV)<delta) :
                    duplicate = False
                    for j in range(0,outhistos[hpath].numPoints) :
                        if(outhistos[hpath].points[j].x==aos[hpath].points[i].x) :
                            duplicate = True
                            break
                    if(not duplicate) :
                        outhistos[hpath].addPoint(aos[hpath].points[i])
                    matched = True
                    break
            if(matched) : continue
            for i in range(0,aos[hpath].numPoints) :
                xmin = aos[hpath].points[i].xMin
                xmax = aos[hpath].points[i].xMax
                if("KLOE_2009_I797438" in hpath or "KLOE_2005_I655225" in hpath or
                   "FENICE_1994_I377833" in hpath) :
                   xmin=math.sqrt(xmin)
                   xmax=math.sqrt(xmax)
                if((energy    > xmin and energy    < xmax) or
                   (energyMeV > xmin and energyMeV < xmax) ) :
                    duplicate = False
                    for j in range(0,outhistos[hpath].numPoints) :
                        if(outhistos[hpath].points[j].x==aos[hpath].points[i].x) :
                            duplicate = True
                            break
                    if(not duplicate) :
                        outhistos[hpath].addPoint(aos[hpath].points[i])
                    break
    if len(outhistos) == 0: continue
    yoda.writeYODA(outhistos,"LowEnergy-EE-%s-%s.yoda" % (runType,args[0]))

