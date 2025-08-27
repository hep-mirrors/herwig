import logging,yoda
def mergeHistos(histUE, histJet, pTMerge=10.):
    title=""
    path=""
    if hasattr(histUE, 'title'):
        title=histUE.title()
    if hasattr(histUE, 'path'):
        path=histUE.path()
    if(type(histUE)==yoda.core.Counter) :
        newHisto = yoda.core.Counter(path,title)
    elif(type(histUE)==yoda.core.Scatter2D) :
        newHisto = yoda.core.Scatter2D(path,title)
    elif(type(histUE)==yoda.core.Profile1D) :
        newHisto = yoda.core.Profile1D(path,title)
        for i in range(0,histUE.numBins()) :
            newHisto.addBin(histUE.bins()[i].xMin(),
                            histUE.bins()[i].xMax())
    elif(type(histUE)==yoda.core.Histo1D) :
        newHisto = yoda.core.HistoUE1D(path,title)
        for i in range(0,histUE.numBins()) :
            newHisto.addBin(histUE.bins()[i].xMin(),
                            histUE.bins()[i].xMax())
    elif(type(histUE)==yoda.core.BinnedEstimate1D) :
        newHisto = yoda.core.BinnedEstimate1D(histUE.xEdges(),path,title)
        if len(histUE.maskedBins()) > 0:
            newHisto.maskBins(histUE.maskedBins())
    else :
        logging.error("Histogram %s is of unknown type" % histUE.path())
        sys.exit(1)

    if type(newHisto)==yoda.core.Profile1D or type(newHisto)==yoda.core.Histo1D  :
        for i in range(0,newHisto.numBins()) :
            if newHisto.bins()[i].xMin() > pTMerge :
                newHisto.bins()[i] += histJet.bins()[i]
            else :
                newHisto.bins()[i] += histUE.bins()[i]
    elif type(newHisto)==yoda.core.BinnedEstimate1D :
        for b in newHisto.bins():
            idx = b.index()
            if b.xMin() > pTMerge :
                oldb = histJet.bin(idx)
            else :
                oldb = histUE.bin(idx)
            b.setVal(oldb.val())
            for key in  oldb.sources() :
                b.setErr(oldb.err(key),key)
    else :
        print(path,title,newHisto)
        print(histUE, histJet)
        quit()
    return newHisto
