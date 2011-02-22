import re
import os
import histogram
td_command = '/export/pc/bin/tdps'

def compareCrossSections(fname1,fname2) :
    f1 = open(fname1)
    test = f1.readline()
    while test:
        test = f1.readline()
        loc = test.find("Total (from unweighted events)")
        if ( loc >= 0 ) :
            temp1 = test.rsplit(")")
            exponent = temp1[2]
            temp2 = temp1[1].rsplit("(")
            errors = temp2[1]
            temp1 = temp2[0].rsplit("  ",1)
            mantisa = temp1[1]
            p = re.compile( '[0-9]')
            error2 = p.sub( '0', mantisa)
            newerror = error2[0:error2.rindex("0")] + errors
            cross1 = float(mantisa  + exponent)
            error1 = float(newerror + exponent)
    f1.close()
    f2 = open(fname2)
    test = f2.readline()
    while test:
        test = f2.readline()
        loc = test.find("Total (from unweighted events)")
        if ( loc >= 0 ) :
            temp1 = test.rsplit(")")
            exponent = temp1[2]
            temp2 = temp1[1].rsplit("(")
            errors = temp2[1]
            temp1 = temp2[0].rsplit("  ",1)
            mantisa = temp1[1]
            p = re.compile( '[0-9]')
            error2 = p.sub( '0', mantisa)
            newerror = error2[0:error2.rindex("0")] + errors
            cross2 = float(mantisa  + exponent)
            error2 = float(newerror + exponent)
    ratio = (cross1-cross2)**2/(error1**2+error2**2)
    return [cross1,error1,cross2,error2,ratio]

def compareLEPQuarks(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LEP-Quarks.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LEP-Quarks-QuarksTest.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    fo = open(fname,'w')
    diff = histogram.compareTopdrawFiles(fname1,fname2,fo,40)
    fo.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    diff[1] /= float(diff[0])
    ws = "|| Quarks || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" \
        % (output[0],output[1],output[2],output[3],output[4],diff[1],op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareLEPLeptons(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LEP-Leptons.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LEP-Leptons-LeptonsTest.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    fo = open(fname,'w')
    diff = histogram.compareTopdrawFiles(fname1,fname2,fo,28)
    fo.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    diff[1] /= float(diff[0])
    ws = "|| Leptons || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],diff[1],op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareLEPVV(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LEP-VV.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LEP-VV-VVTest.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    fo = open(fname,'w')
    diff = histogram.compareTopdrawFiles(fname1,fname2,fo,6)
    fo.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    diff[1] /= float(diff[0])
    ws = "|| VV || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],diff[1],op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareLEPVH(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LEP-VH.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LEP-VH-VHTest.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    fo = open(fname,'w')
    diff = histogram.compareTopdrawFiles(fname1,fname2,fo,10)
    fo.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    diff[1] /= float(diff[0])
    ws = "|| VH || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],diff[1],op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareLEPVBF(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LEP-VBF.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LEP-VBF-TestVBF.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    fo = open(fname,'w')
    diff = histogram.compareTopdrawFiles(fname1,fname2,fo,16)
    fo.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    diff[1] /= float(diff[0])
    ws = "|| VBF || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],diff[1],op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareCharmShapes(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LEP-BB.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LEP-BB-BELLECharm.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    oname="LEP-BB-Shapes.top"
    fo = open(oname,'w')
    diff1 = histogram.compareTopdrawFiles(fname1,fname2,fo,6,True)
    fname="LEP-BB-CLEOCharm.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    diff2 = histogram.compareTopdrawFiles(fname1,fname2,fo,4,True)
    fo.close()
    diff = [diff1[0]+diff2[0],diff1[1]+diff2[1]]
    op1 = plotLocation + "/" + oname
    op2 = op1.replace(".top",".ps")
    diff[1] /= float(diff[0])
    ws = "|| Charm || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],diff[1],op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +oname 
    os.system(tdstring)


def compareLEPShapes(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LEP-default.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    oname="LEP-Shapes.top"
    fo = open(oname,'w')
    totalDegree = 0
    totalChi = 0.
    # read the event shapes
    fname="LEP-default-LEPEvent.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    for i in range(0,16) :
        h1 = histogram.readHistogram(f1,True)
        h2 = histogram.readHistogram(f2,True)
        h1.write(fo,True,"BLACK",True,False,False,True)
        h2.write(fo,False,"RED",True,False,False,True)
        if(h1.hasData) : h1.writeData(fo)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    f1.close()
    f2.close()
    # read the identified particle spectra
    fname="LEP-default-LEPIdent.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    for i in range(0,54) :
        h1 = histogram.readHistogram(f1,True)
        h2 = histogram.readHistogram(f2,True)
        h1.write(fo,True,"BLACK",True,False,False,True)
        h2.write(fo,False,"RED",True,False,False,True)
        if(h1.hasData) : h1.writeData(fo)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    f1.close()
    f2.close()
    # single particle spectra
    fname="LEP-default-LEPSingle.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    for i in range(0,7) :
        h1 = histogram.readHistogram(f1,True)
        h2 = histogram.readHistogram(f2,True)
        h1.write(fo,True,"BLACK",True,False,False,True)
        h2.write(fo,False,"RED",True,False,False,True)
        if(h1.hasData) : h1.writeData(fo)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    f1.close()
    f2.close()
    # LEP jets
    fname="LEP-default-LEPJet.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    for i in range(0,4) :
        h1 = histogram.readHistogram(f1,True)
        h2 = histogram.readHistogram(f2,True)
        h1.write(fo,True,"BLACK",True,False,True,True)
        h2.write(fo,False,"RED",True,False,True,True)
        if(h1.hasData) : h1.writeData(fo)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    for i in range(4,13) :
        h1 = histogram.readLine(f1,True)
        h2 = histogram.readLine(f2,True)
        h1.write(fo,True,"BLACK",True,False,True,True)
        h2.write(fo,False,"RED",True,False,True,True)
        if(h1.hasData) : h1.writeData(fo)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    f1.close()
    f2.close()
    fname="LEP-default-LEPFourJet.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    for i in range(0,4) :
        h1 = histogram.readHistogram(f1,True)
        h2 = histogram.readHistogram(f2,True)
        h1.write(fo,True,"BLACK",True,False,False,False)
        h2.write(fo,False,"RED",True,False,False,False)
        if(h1.hasData) : h1.writeData(fo)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    f1.close()
    f2.close()
    fname="LEP-default-BFrag.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    for i in range(0,2) :
        h1 = histogram.readHistogram(f1,True)
        h2 = histogram.readHistogram(f2,True)
        h1.write(fo,True,"BLACK",True,False,False,False)
        h2.write(fo,False,"RED",True,False,False,False)
        if(h1.hasData) : h1.writeData(fo)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    f1.close()
    f2.close()
    fo.close()
    op1 = plotLocation + "/" + oname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| LEP || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +oname 
    os.system(tdstring)

def compareW(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LHC-W.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LHC-W-WTest.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(0,38) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        if(i>25 and i < 30) or (i>31 and i <36) :
            h1.write(fo,True,"BLACK",True,False,False,i%2==0)
            h2.write(fo,False,"RED",True,False,False,i%2==0)
            out  = h1.writeDifference(h2,fo,True,"RED")
            totalDegree += out[0]
            totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| W || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareZ(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LHC-Z.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LHC-Z-ZTest.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,39) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        if(i>24 and i < 27) or (i>30 and i <33) :
            h1.write(fo,True,"BLACK",True,False,False,i%2==0)
            h2.write(fo,False,"RED",True,False,False,i%2==0)
            out  = h1.writeDifference(h2,fo,True,"RED")
            totalDegree += out[0]
            totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| Z || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareWJet(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LHC-WJet.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LHC-WJet-WJetTest.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,25) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        logPlot = i==1 or i==5 or i==9 or i==13 or i==16 or i==19 or i==22
        h1.write(fo,True,"BLACK",True,False,False,logPlot)
        h2.write(fo,False,"RED",True,False,False,logPlot)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| WJet || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareZJet(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LHC-ZJet.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LHC-ZJet-ZJetTest.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,11) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        logPlot = i==1 or i==5 or i==8
        h1.write(fo,True,"BLACK",True,False,False,logPlot)
        h2.write(fo,False,"RED",True,False,False,logPlot)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| ZJet || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareHiggs(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LHC-Higgs.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LHC-Higgs-HiggsTest.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,9) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        if(i==5 or i==6) :
            logPlot = i==6
            h1.write(fo,True,"BLACK",True,False,False,logPlot)
            h2.write(fo,False,"RED",True,False,False,logPlot)
            out  = h1.writeDifference(h2,fo,True,"RED")
            totalDegree += out[0]
            totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| Higgs || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareHiggsJet(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LHC-HiggsJet.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LHC-HiggsJet-HiggsJetTest.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,9) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        logPlot = i%2==0
        h1.write(fo,True,"BLACK",True,False,False,logPlot)
        h2.write(fo,False,"RED",True,False,False,logPlot)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| HiggsJet || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)


def compareWWVBF(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LHC-WWVBF.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LHC-WWVBF-WWVBF-Test.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,11) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        logPlot = i==3 or i==4 or i==7 or i==8
        h1.write(fo,True,"BLACK",True,False,False,logPlot)
        h2.write(fo,False,"RED",True,False,False,logPlot)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| WWVBF || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareZZVBF(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LHC-ZZVBF.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LHC-ZZVBF-ZZVBF-Test.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,11) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        logPlot =  i==3 or i==4 or i==7 or i==8
        h1.write(fo,True,"BLACK",True,False,False,logPlot)
        h2.write(fo,False,"RED",True,False,False,logPlot)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| ZZVBF || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareVBF(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LHC-VBF.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LHC-VBF-VBF-Test.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,11) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        logPlot =  i==3 or i==4 or i==7 or i==8
        h1.write(fo,True,"BLACK",True,False,False,logPlot)
        h2.write(fo,False,"RED",True,False,False,logPlot)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| VBF || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareWW(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LHC-WW.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LHC-WW-WWTest.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,11) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        if(i<7 or i==10) :
            logPlot = i==1 or i==4 or i==10
            h1.write(fo,True,"BLACK",True,False,False,logPlot)
            h2.write(fo,False,"RED",True,False,False,logPlot)
            out  = h1.writeDifference(h2,fo,True,"RED")
            totalDegree += out[0]
            totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| WW || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareZZ(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LHC-ZZ.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LHC-ZZ-ZZTest.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,11) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        if(i>=7) :
            logPlot = i==1 or i==4 or i==7 or i==10
            h1.write(fo,True,"BLACK",True,False,False,logPlot)
            h2.write(fo,False,"RED",True,False,False,logPlot)
            out  = h1.writeDifference(h2,fo,True,"RED")
            totalDegree += out[0]
            totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| ZZ || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareWZ(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LHC-WZ.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LHC-WZ-WZTest.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,11) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        logPlot = i==1 or i==4 or i==7 or i==10
        h1.write(fo,True,"BLACK",True,False,False,logPlot)
        h2.write(fo,False,"RED",True,False,False,logPlot)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| WZ || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareWGamma(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LHC-WGamma.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LHC-WGamma-WGammaTest.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,18) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        if(i<9 or i>=13) :
            logPlot = i==1 or i==5 or i==13 or i==17
            h1.write(fo,True,"BLACK",True,False,False,logPlot)
            h2.write(fo,False,"RED",True,False,False,logPlot)
            out  = h1.writeDifference(h2,fo,True,"RED")
            totalDegree += out[0]
            totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| WGamma || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareZGamma(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LHC-ZGamma.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LHC-ZGamma-ZGammaTest.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,18) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        if(i>=9) :
            logPlot = i==9 or i==13 or i==17
            h1.write(fo,True,"BLACK",True,False,False,logPlot)
            h2.write(fo,False,"RED",True,False,False,logPlot)
            out  = h1.writeDifference(h2,fo,True,"RED")
            totalDegree += out[0]
            totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| ZGamma || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareGammaGamma(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LHC-GammaGamma.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LHC-GammaGamma-GammaGammaTest.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,14) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        if(i != 3 and i !=10 and i !=11 and i !=12) :
            logPlot = True
            h1.write(fo,True,"BLACK",True,False,False,logPlot)
            h2.write(fo,False,"RED",True,False,False,logPlot)
            out  = h1.writeDifference(h2,fo,True,"RED")
            totalDegree += out[0]
            totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| Gamma Gamma || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareGammaJet(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LHC-GammaJet.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LHC-GammaJet-GammaJetTest.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,11) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        logPlot = i%2==0 and not i==10
        h1.write(fo,True,"BLACK",True,False,False,logPlot)
        h2.write(fo,False,"RED",True,False,False,logPlot)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| Gamma Jet || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareWH(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LHC-WH.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LHC-WH-WHTest.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,29) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        logPlot = i==2 or i==6 or i==10 or i==14 or i==17 or i==20 or i==23 or i==26
        h1.write(fo,True,"BLACK",True,False,False,logPlot)
        h2.write(fo,False,"RED",True,False,False,logPlot)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| WH || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareZH(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LHC-ZH.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LHC-ZH-ZHTest.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,15) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        logPlot = i==2 or i==6 or i==9 or i==12
        h1.write(fo,True,"BLACK",True,False,False,logPlot)
        h2.write(fo,False,"RED",True,False,False,logPlot)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| ZH || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)


def compareQCD(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LHC-QCD.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LHC-QCD-QCDTest.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,6) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        logPlot = i==2 or i==4
        h1.write(fo,True,"BLACK",True,False,False,logPlot)
        h2.write(fo,False,"RED",True,False,False,logPlot)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| QCD || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareQCDFast(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LHC-QCDFast.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LHC-QCDFast-QCDTest.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,6) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        logPlot =  i==2 or i==4
        h1.write(fo,True,"BLACK",True,False,False,logPlot)
        h2.write(fo,False,"RED",True,False,False,logPlot)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| QCDFast || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def comparebbH(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LHC-bbH.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LHC-bbH-bbHTest.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,24) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        logPlot =  i==1 or i==4 or i==7 or i==10 or i==11 or i==13 or i==15 or i==17
        h1.write(fo,True,"BLACK",True,False,False,logPlot)
        h2.write(fo,False,"RED",True,False,False,logPlot)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| bbH || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def comparettH(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LHC-ttH.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LHC-ttH-ttHTest.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,24) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        logPlot = i==1 or i==4 or i==7 or i==10 or i==11 or i==13 or i==15  or i==17
        h1.write(fo,True,"BLACK",True,False,False,logPlot)
        h2.write(fo,False,"RED",True,False,False,logPlot)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| ttH || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareWShower(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LHC-WShower.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LHC-WShower-VTest.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,22) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        logPlot = i<=4
        h1.write(fo,True,"BLACK",True,False,False,logPlot)
        h2.write(fo,False,"RED",True,False,False,logPlot)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| W || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareZShower(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LHC-ZShower.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LHC-ZShower-VTest.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,22) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        logPlot = i<=4
        h1.write(fo,True,"BLACK",True,False,False,logPlot)
        h2.write(fo,False,"RED",True,False,False,logPlot)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| Z || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareWPowheg(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LHC-WShower-Powheg.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LHC-WShower-Powheg-VTest.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,22) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        logPlot = i<=4
        h1.write(fo,True,"BLACK",True,False,False,logPlot)
        h2.write(fo,False,"RED",True,False,False,logPlot)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| W-Powheg || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareZPowheg(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LHC-ZShower-Powheg.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LHC-ZShower-Powheg-VTest.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,22) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        logPlot = i<=4
        h1.write(fo,True,"BLACK",True,False,False,logPlot)
        h2.write(fo,False,"RED",True,False,False,logPlot)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| Z-Powheg || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareHJet(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LHC-HJet.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LHC-HJet-HTest.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,22) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        logPlot = i<=4
        h1.write(fo,True,"BLACK",True,False,False,logPlot)
        h2.write(fo,False,"RED",True,False,False,logPlot)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| H || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareHPowheg(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LHC-HJet-Powheg.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LHC-HJet-Powheg-HTest.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,22) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        logPlot = i<=4
        h1.write(fo,True,"BLACK",True,False,False,logPlot)
        h2.write(fo,False,"RED",True,False,False,logPlot)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| H-Powheg|| %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareWHJet(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LHC-WHJet.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LHC-WHJet-VHTest.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,38) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        logPlot = i<=8
        h1.write(fo,True,"BLACK",True,False,False,logPlot)
        h2.write(fo,False,"RED",True,False,False,logPlot)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| WH || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareZHJet(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LHC-ZHJet.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LHC-ZHJet-VHTest.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,38) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        logPlot = i<=8
        h1.write(fo,True,"BLACK",True,False,False,logPlot)
        h2.write(fo,False,"RED",True,False,False,logPlot)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| ZH || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareWHJetPowheg(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LHC-WHJet-Powheg.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LHC-WHJet-Powheg-VHTest.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,38) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        logPlot = i<=8
        h1.write(fo,True,"BLACK",True,False,False,logPlot)
        h2.write(fo,False,"RED",True,False,False,logPlot)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| WH-Powheg || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareZHJetPowheg(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LHC-ZHJet-Powheg.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LHC-ZHJet-Powheg-VHTest.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,38) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        logPlot = i<=8
        h1.write(fo,True,"BLACK",True,False,False,logPlot)
        h2.write(fo,False,"RED",True,False,False,logPlot)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| ZH-Powheg || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareNeutralCurrent(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="DIS-Neutral.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="DIS-Neutral-NeutralTest.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,7) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        logPlot = i!=6
        h1.write(fo,True,"BLACK",True,False,False,logPlot)
        h2.write(fo,False,"RED",True,False,False,logPlot)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| Neutral || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareChargedCurrent(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="DIS-Charged.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="DIS-Charged-ChargedTest.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,7) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        logPlot = i!=6
        h1.write(fo,True,"BLACK",True,False,False,logPlot)
        h2.write(fo,False,"RED",True,False,False,logPlot)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| Charged || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareGammaFF(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="Gamma-FF.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="Gamma-FF-TestFF.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,73) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        logPlot = False
        h1.write(fo,True,"BLACK",True,False,False,logPlot)
        h2.write(fo,False,"RED",True,False,False,logPlot)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| gamma gamma -> fermions || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareGammaWW(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="Gamma-WW.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="Gamma-WW-TestWW.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,9) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        logPlot = False
        h1.write(fo,True,"BLACK",True,False,False,logPlot)
        h2.write(fo,False,"RED",True,False,False,logPlot)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| gamma gamma -> WW || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareGammaP(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="Gamma-P.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="Gamma-P-TestP.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,6) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        logPlot = i==2 or i==4
        h1.write(fo,True,"BLACK",True,False,False,logPlot)
        h2.write(fo,False,"RED",True,False,False,logPlot)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| gamma hadron -> jets || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)
def compareTopDecay(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LEP-TopDecay.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LEP-TopDecay-TopDecay.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,2) :
        h1 = histogram.readHistogram(f1,False)
        h2 = histogram.readHistogram(f2,False)
        logPlot = i==2 or i==4
        h1.write(fo,True,"BLACK",True,False,False,logPlot)
        h2.write(fo,False,"RED",True,False,False,logPlot)
        out  = h1.writeDifference(h2,fo,True,"RED")
        totalDegree += out[0]
        totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| t tbar || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)

def compareTop(directory1,directory2,wiki,plotLocation) :
    # first compare the section sections
    fname="LHC-Top.out"
    fname1=directory1 + fname
    fname2=directory2 + fname
    output = compareCrossSections(fname1,fname2)
    # now compare the distributions
    fname="LHC-Top-TopTest.top"
    fname1=directory1 + fname
    fname2=directory2 + fname
    f1 = open(fname1)
    f2 = open(fname2)
    fo = open(fname,'w')
    totalDegree = 0
    totalChi = 0.
    for i in range(1,14) :
        nplot = 1
        if(i==1 or i==3 or i==5 or i==7 or i==9) : nplot = 2
        h1 = histogram.readHistogram(f1,False,nplot)
        h2 = histogram.readHistogram(f2,False,nplot)
        if ( i!=2 and i!=4 and i!=10) :
            logPlot = False
            if(nplot == 1 ) :
                h1.write(fo,True,"BLACK",True,False,False,logPlot)
                h2.write(fo,False,"RED",True,False,False,logPlot)
                out  = h1.writeDifference(h2,fo,True,"RED")
                totalDegree += out[0]
                totalChi += out[1]
            else :
                h1[0].write(fo,True,"BLACK",True,False,False,logPlot)
                h2[0].write(fo,False,"RED",True,False,False,logPlot)
                h1[1].write(fo,False,"BLACK DASHES",True,False,False,logPlot)
                h2[1].write(fo,False,"RED   DASHES",True,False,False,logPlot)
                out  = h1[0].writeDifference(h2[0],fo,True,"RED")
                totalDegree += out[0]
                totalChi += out[1]
                out  = h1[1].writeDifference(h2[1],fo,False,"RED DASHES")
                totalDegree += out[0]
                totalChi += out[1]
    fo.close()
    f1.close()
    f2.close()
    op1 = plotLocation + "/" + fname
    op2 = op1.replace(".top",".ps")
    totalChi /= float(totalDegree)
    ws = "|| ttbar || %s || %s || %s || %s || %s || %s || [%s top] || [%s ps] ||\n" % (output[0],output[1],output[2],output[3],output[4],totalChi,op1,op2)
    wiki.write(ws)
    tdstring = td_command + " " +fname 
    os.system(tdstring)



