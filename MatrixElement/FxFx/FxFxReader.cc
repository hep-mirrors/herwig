// -*- C++ -*-
//
// FxFxReader.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FxFxReader class.
//

#include "FxFxReader.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Command.h"
//#include "config.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDF/PartonExtractor.h"
#include "ThePEG/PDF/NoPDF.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/EventRecord/TmpTransform.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Handlers/XComb.h"
#include "ThePEG/Handlers/CascadeHandler.h"
#include "FxFxEventHandler.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Utilities/HoldFlag.h"
#include "ThePEG/Utilities/Debug.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"

using namespace ThePEG;

FxFxReader::FxFxReader(bool active)
  : theNEvents(0), position(0), reopened(0), theMaxScan(-1), scanning(false),
    isActive(active), theCacheFileName(""), doCutEarly(true),
    preweight(1.0), reweightPDF(false), doInitPDFs(false),
    theMaxMultCKKW(0), theMinMultCKKW(0), lastweight(1.0), maxFactor(1.0), optionalnpLO(0), optionalnpNLO(0),
    weightScale(1.0*picobarn), skipping(false), theMomentumTreatment(0),
    useWeightWarnings(true),theReOpenAllowed(true), theIncludeSpin(true) {}

FxFxReader::FxFxReader(const FxFxReader & x)
  : HandlerBase(x), LastXCombInfo<>(x), heprup(x.heprup), hepeup(x.hepeup),
    inData(x.inData), inPDF(x.inPDF), outPDF(x.outPDF),
    thePartonExtractor(x.thePartonExtractor), thePartonBins(x.thePartonBins),
    theXCombs(x.theXCombs), theCuts(x.theCuts),
    theNEvents(x.theNEvents), position(x.position), reopened(x.reopened),
    theMaxScan(x.theMaxScan), scanning(false),
    isActive(x.isActive),
    theCacheFileName(x.theCacheFileName), doCutEarly(x.doCutEarly),
    stats(x.stats), statmap(x.statmap),
    thePartonBinInstances(x.thePartonBinInstances),
    reweights(x.reweights), preweights(x.preweights),
    preweight(x.preweight), reweightPDF(x.reweightPDF),
    doInitPDFs(x.doInitPDFs),
    theMaxMultCKKW(x.theMaxMultCKKW), theMinMultCKKW(x.theMinMultCKKW),
    lastweight(x.lastweight), maxFactor(x.maxFactor),
    weightScale(x.weightScale), xSecWeights(x.xSecWeights),
    maxWeights(x.maxWeights), skipping(x.skipping),
    theMomentumTreatment(x.theMomentumTreatment),
    useWeightWarnings(x.useWeightWarnings),
    theReOpenAllowed(x.theReOpenAllowed),
    theIncludeSpin(x.theIncludeSpin) {}

FxFxReader::~FxFxReader() {}

void FxFxReader::doinitrun() {
  HandlerBase::doinitrun();
  stats.reset();
  for ( StatMap::iterator i = statmap.begin(); i != statmap.end(); ++i )
    i->second.reset();
  open();
  if ( cacheFileName().length() ) openReadCacheFile();
  position = 0;
  reopened = 0;
}

bool FxFxReader::preInitialize() const {
  if ( HandlerBase::preInitialize() ) return true;
  if ( doInitPDFs && ! ( inPDF.first && inPDF.second ) ) return true;
  return false;
}

void FxFxReader::doinit() {
  HandlerBase::doinit();
  open();
  close();
  if ( !heprup.IDBMUP.first || !heprup.IDBMUP.second )
    Throw<FxFxInitError>()
      << "No information about incoming particles were found in "
      << "FxFxReader '" << name() << "'." << Exception::warning;
  inData = make_pair(getParticleData(heprup.IDBMUP.first),
		     getParticleData(heprup.IDBMUP.second));
  if ( heprup.EBMUP.first <= 0.0 || heprup.EBMUP.second <= 0.0 )
    Throw<FxFxInitError>()
    << "No information about the energy of incoming particles were found in "
    << "FxFxReader '" << name() << "'." << Exception::warning;
  
  if ( doInitPDFs && ! ( inPDF.first && inPDF.second ) ) {
    initPDFs();
    if ( ! ( inPDF.first && inPDF.second ) ) Throw<InitException>()
      << "FxFxReader '" << name()
      << "' could not create PDFBase objects in pre-initialization."
      << Exception::warning;
  }
  else if ( !inPDF.first || !inPDF.second ) Throw<FxFxInitError>()
    << "No information about the PDFs of incoming particles were found in "
    << "FxFxReader '" << name() << "'." << Exception::warning;
}

void FxFxReader::initPDFs() {
  if ( inPDF.first && inPDF.second ) return;

  string remhname;
  if ( heprup.PDFSUP.first && !inPDF.first) {
    inPDF.first = dynamic_ptr_cast<PDFPtr>
      (generator()->preinitCreate("ThePEG::LHAPDF", fullName() + "/PDFA",
				  "ThePEGLHAPDF.so"));
    if ( !inPDF.first ) {
      Throw<InitException>()
	<< "FxFxReader '" << name() << "' could not use information "
	<< "about the PDFs used because the LHAPDF library was not properly "
	"defined." << Exception::warning;
      return;
    }
    remhname = fullName() + "/DummyRemH";
    generator()->preinitCreate("ThePEG::NoRemnants", remhname);
    generator()->preinitInterface(inPDF.first, "RemnantHandler",
				  "set", remhname);
    if ( heprup.PDFGUP.first > 0 && heprup.PDFGUP.first < 10 ) {
      ostringstream os;
      os << heprup.PDFGUP.first << " " << heprup.PDFSUP.first;
      generator()->preinitInterface(inPDF.first, "PDFLIBNumbers",
				    "set", os.str());
    } else {
      ostringstream os;
      os << heprup.PDFGUP.first*1000 + heprup.PDFSUP.first;
      generator()->preinitInterface(inPDF.first, "PDFNumber",
				    "set", os.str());
    }
    generator()->preinitInterface(inPDF.first, "RangeException",
				    "newdef", "Freeze");
  }

  if ( heprup.PDFSUP.second && !inPDF.second) {
    inPDF.second = dynamic_ptr_cast<PDFPtr>
      (generator()->preinitCreate("ThePEG::LHAPDF", fullName() + "/PDFB",
				  "ThePEGLHAPDF.so"));
    if ( !inPDF.second ) {
      Throw<InitException>()
	<< "FxFxReader '" << name() << "' could not use information "
	<< "about the PDFs used because the LHAPDF library was not properly "
	"defined." << Exception::warning;
      return;
    }
    if ( remhname == "" ) {
      remhname = fullName() + "/DummyRemH";
      generator()->preinitCreate("ThePEG::NoRemnants", remhname);
    }
    generator()->preinitInterface(inPDF.second, "RemnantHandler",
				  "set", remhname);      
    if ( heprup.PDFGUP.second > 0 && heprup.PDFGUP.second < 10 ) {
      ostringstream os;
      os << heprup.PDFGUP.second << " " << heprup.PDFSUP.second;
      generator()->preinitInterface(inPDF.second, "PDFLIBNumbers",
				    "set", os.str());
    } else {
      ostringstream os;
      os << heprup.PDFGUP.second*1000 + heprup.PDFSUP.second;
      generator()->preinitInterface(inPDF.second, "PDFNumber",
				    "set", os.str());
    }
    generator()->preinitInterface(inPDF.second, "RangeException",
				    "newdef", "Freeze");
  }
  
  if ( ! ( inPDF.first && inPDF.second ) ) Throw<InitException>()
    << "FxFxReader '" << name()
    << "' could not find information about the PDFs used."
    << Exception::warning;
}

void FxFxReader::initialize(FxFxEventHandler & eh) {
  Energy2 Smax = ZERO;
  double Y = 0.0;
  if ( !theCuts ) {
    theCuts = eh.cuts();
    if ( !theCuts ) Throw<FxFxInitError>()
      << "No Cuts object was assigned to the FxFxReader '"
      << name() << "' nor was one\nassigned to the controlling "
      << "FxFxEventHandler '" << eh.name() << "'.\nAt least one of them "
      << "needs to have a Cuts object." << Exception::runerror;

    Smax = cuts().SMax();
    Y = cuts().Y();
  }

  theCKKW = eh.CKKWHandler();

  if ( !partonExtractor() ) {
    thePartonExtractor = eh.partonExtractor();
    if ( !partonExtractor() )  Throw<FxFxInitError>()
      << "No PartonExtractor object was assigned to the FxFxReader '"
      << name() << "' nor was one\nassigned to the controlling "
      << "FxFxEventHandler '" << eh.name() << "'.\nAt least one of them "
      << "needs to have a PartonExtractor object." << Exception::runerror;
  }
  open();

  Energy emax = 2.0*sqrt(heprup.EBMUP.first*heprup.EBMUP.second)*GeV;
  theCuts->initialize(sqr(emax),
		      0.5*log(heprup.EBMUP.first/heprup.EBMUP.second));
  if ( Smax > ZERO && ( Smax != cuts().SMax() || Y != cuts().Y() ) )
    Throw<FxFxInitError>()
      << "The FxFxReader '" << name() << "' uses the same Cuts object "
      << "as another FxFxReader which has not got the same energies of "
      << "the colliding particles. For the generation to work properly "
      << "different FxFxReader object with different colliding particles "
      << "must be assigned different (although possibly identical) Cuts "
      << "objects." << Exception::warning;

  thePartonBins = partonExtractor()->getPartons(emax, inData, cuts());
  for ( int i = 0, N = partonBins().size(); i < N; ++i ) {
    theXCombs[partonBins()[i]] =
      new_ptr(XComb(emax, inData, &eh, partonExtractor(), CKKWHandler(),
		    partonBins()[i], theCuts));
    partonExtractor()->nDims(partonBins()[i]);
  }
  outPDF = make_pair(partonExtractor()->getPDF(inData.first),
		     partonExtractor()->getPDF(inData.second));


  close();

  if ( !heprup.IDWTUP && useWeightWarnings )
    Throw<FxFxInitError>()
      << "No information about the weighting scheme was found. The events "
      << "produced by FxFxReader " << name()
      << " may not be sampled correctly." << Exception::warning;

  if ( abs(heprup.IDWTUP) > 1 && !eh.weighted() && useWeightWarnings  )
    Throw<FxFxInitError>()
      << "FxFxReader " << name() << " has the IDWTUP flag set to "
      << heprup.IDWTUP << " which is not supported by this reader, the "
      << "produced events may not be sampled correctly. It is up to "
      << "sub-classes of FxFxReader to correctly convert to match IDWTUP "
      << "+/- 1. Will try to make intelligent guesses to get "
      << "correct statistics.\nIn most cases this should be sufficient. "
      << "Unset <interface>WeightWarnings</interface> to avoid this message,"
      << "or set <interface>Weighted</interface> to on."
      << Exception::warning;

  if ( heprup.IDWTUP != eh.weightOption() && abs(heprup.IDWTUP) < 3 &&
       useWeightWarnings  )
    Throw<FxFxInitError>()
      << "FxFxReader " << name() << " has the IDWTUP flag set to "
      << heprup.IDWTUP 
      << ", which does not correspond\nto the weight option "
      << eh.weightOption() << " set in "
      << "the FxFxEventHandler " << eh.name() << ".\n\n"
      << "Use the following handler setting instead:\n"
      << "  set " << eh.name() << ":WeightOption " << heprup.IDWTUP
      << "\nWill try to make intelligent guesses to get "
      << "correct statistics. In most cases this should be sufficient. "
      << "Unset <interface>WeightWarnings</interface> to avoid this message"
      << Exception::warning;

  scan();
  initStat();
}


long FxFxReader::scan() {
  
  open();

  // Shall we write the events to a cache file for fast reading? If so
  // we write to a temporary file if the caches events should be
  // randomized.
  if ( cacheFileName().length() ) openWriteCacheFile();

  // Keep track of the number of events scanned.
  long neve = 0;
  long cuteve = 0;
  bool negw = false;

  // If the open() has not already gotten information about subprocesses
  // and cross sections we have to scan through the events.
  if ( !heprup.NPRUP || cacheFile() || abs(heprup.IDWTUP) != 1 ) { // why scan if IDWTUP != 1?

    HoldFlag<> isScanning(scanning);

    double oldsum = 0.0;
    vector<int> lprup;
    vector<double> newmax;
    vector<long> oldeve;
    vector<long> neweve;
    vector<double> sumlprup;
    vector<double> sumsqlprup;  
    vector<long> nscanned;
    for ( int i = 0; ( maxScan() < 0 || i < maxScan() ) && readEvent(); ++i ) {
      if ( !checkPartonBin() ) Throw<FxFxInitError>()
        << "Found event in LesHouchesReader '" << name()
        << "' which cannot be handeled by the assigned PartonExtractor '"
        << partonExtractor()->name() << "'." << Exception::runerror;
      vector<int>::iterator idit =
        find(lprup.begin(), lprup.end(), hepeup.IDPRUP);
      int id = lprup.size();
      if ( idit == lprup.end() ) {
        lprup.push_back(hepeup.IDPRUP);
        newmax.push_back(0.0);
        neweve.push_back(0);
        oldeve.push_back(0);
        sumlprup.push_back(0.);
        sumsqlprup.push_back(0.);
        nscanned.push_back(0);
      } else {
        id = idit - lprup.begin();
      }
      ++neve;
      ++oldeve[id];
      oldsum += hepeup.XWGTUP;
      sumlprup[id] += hepeup.XWGTUP;
      sumsqlprup[id] += sqr(hepeup.XWGTUP);
      ++nscanned[id];
      if ( cacheFile() ) {
        if ( eventWeight() == 0.0 ) {
          ++cuteve;
          continue;
        }
        cacheEvent();
      }
      ++neweve[id];
      newmax[id] = max(newmax[id], abs(eventWeight()));
      if ( eventWeight() < 0.0 ) negw = true;
    } //end of scanning events
    xSecWeights.resize(oldeve.size(), 1.0);
    for ( int i = 0, N = oldeve.size(); i < N; ++i )
      if ( oldeve[i] ) xSecWeights[i] = double(neweve[i])/double(oldeve[i]);

    if ( maxScan() < 0 || neve > NEvents() ) NEvents(neve - cuteve);

    if ( lprup.size() == heprup.LPRUP.size() ) {
      for ( int id = 0, N = lprup.size(); id < N; ++id ) {
        vector<int>::iterator idit =
          find(heprup.LPRUP.begin(), heprup.LPRUP.end(), hepeup.IDPRUP);
        if ( idit == heprup.LPRUP.end() ) {
          Throw<FxFxInitError>()
            << "When scanning events, the LesHouschesReader '" << name()
            << "' found undeclared processes."  << Exception::warning;
          heprup.NPRUP = 0;
          break;
        }
        int idh = idit - heprup.LPRUP.begin();
        heprup.XMAXUP[idh] = newmax[id];
      } 
    }
    if ( heprup.NPRUP == 0 ) {
      // No heprup block was supplied or something went wrong.
      heprup.NPRUP = lprup.size();
      heprup.LPRUP.resize(lprup.size());
      heprup.XMAXUP.resize(lprup.size());
      for ( int id = 0, N = lprup.size(); id < N; ++id ) {
        heprup.LPRUP[id] = lprup[id];
        heprup.XMAXUP[id] = newmax[id];
      }
    }
    if ( abs(heprup.IDWTUP) != 1 ) {
      // Try to fix things if abs(heprup.IDWTUP) != 1.
      double sumxsec = 0.0;
      if(abs(heprup.IDWTUP)==3) {
	for ( int id = 0; id < heprup.NPRUP; ++id ) sumxsec += heprup.XSECUP[id];
      }
      else {
	for ( int id = 0; id < heprup.NPRUP; ++id )  {
	  //set the cross section directly from the event weights read
	  heprup.XSECUP[id] = sumlprup[id]/nscanned[id];
	  heprup.XERRUP[id] = (sumsqlprup[id]/nscanned[id] - sqr(sumlprup[id]/nscanned[id])) / nscanned[id];
	  if(fabs(heprup.XERRUP[id]) < 1e-10) { heprup.XERRUP[id] = 0; }
	  if(fabs(newmax[id]) < 1e-10) { newmax[id] = 0; }
	  if(heprup.XERRUP[id] < 0.) {
	    if( heprup.XERRUP[id]/(sumsqlprup[id]/nscanned[id])>-1e-10)
	      heprup.XERRUP[id] = 0.;
	    else {
	      Throw<FxFxInitError>()
		<< "Negative error when scanning events in LesHouschesReader '" << name()
		<< Exception::warning;
	      heprup.XERRUP[id] = 0.;
	    }
	  }
	  heprup.XERRUP[id] = sqrt( heprup.XERRUP[id] );
	  heprup.XMAXUP[id] = newmax[id];
	  //cout << "heprup.XMAXUP[id] = " << heprup.XMAXUP[id]  << endl;
	  sumxsec += heprup.XSECUP[id];
	}
      }
      //cout << "sumxsec = " << sumxsec << endl;
      //weightScale = picobarn*neve*sumxsec/oldsum;
      weightScale = 1.0*picobarn; // temporary fix?
      // cout << "weightscale = " << weightScale/picobarn << endl;
    }
  }

  if ( cacheFile() ) closeCacheFile();

  if ( negw ) heprup.IDWTUP = min(-abs(heprup.IDWTUP), -1);
 
  return neve;
}

void FxFxReader::setWeightScale(long) {}

void FxFxReader::initStat() {

  stats.reset();
  statmap.clear();
  if ( heprup.NPRUP <= 0 ) return;

  double sumx = 0.0;
  xSecWeights.resize(heprup.NPRUP, 1.0);
  maxWeights.clear();
  for ( int ip = 0; ip < heprup.NPRUP; ++ip ) {
    sumx = max(heprup.XMAXUP[ip]*xSecWeights[ip], sumx);
    statmap[heprup.LPRUP[ip]] =
      XSecStat(heprup.XMAXUP[ip]*weightScale*xSecWeights[ip]);
    maxWeights[heprup.LPRUP[ip]] = heprup.XMAXUP[ip];
  }
  stats.maxXSec(sumx*weightScale);
  maxFactor = 1.0;
}

void FxFxReader::increaseMaxXSec(CrossSection maxxsec) {
  for ( int i = 0; i < heprup.NPRUP; ++i )
    statmap[heprup.LPRUP[i]].maxXSec(statmap[heprup.LPRUP[i]].maxXSec()*
      maxxsec/stats.maxXSec());
  maxFactor *= maxxsec/stats.maxXSec();
  stats.maxXSec(maxxsec);
}

tXCombPtr FxFxReader::getXComb() {
  if ( lastXCombPtr() ) return lastXCombPtr();
  fillEvent();
  connectMothers();
  tcPBPair sel = createPartonBinInstances();
  tXCombPtr lastXC = xCombs()[sel];
  // clean up the old XComb object before switching to a new one
  if ( theLastXComb && theLastXComb != lastXC ) 
    theLastXComb->clean();
  theLastXComb = lastXC;
  lastXCombPtr()->subProcess(SubProPtr());
  lastXCombPtr()->setPartonBinInstances(partonBinInstances(),
					sqr(hepeup.SCALUP)*GeV2);
  lastXCombPtr()->lastAlphaS(hepeup.AQCDUP);
  lastXCombPtr()->lastAlphaEM(hepeup.AQEDUP);
  return lastXCombPtr();
}

tSubProPtr FxFxReader::getSubProcess() {
  getXComb();
  if ( subProcess() ) return subProcess();
  lastXCombPtr()->subProcess(new_ptr(SubProcess(lastPartons(), tCollPtr(), this)));
  subProcess()->setOutgoing(outgoing().begin(), outgoing().end());
  subProcess()->setIntermediates(intermediates().begin(),
				 intermediates().end());
  return subProcess();
}

void FxFxReader::fillEvent() {
  if ( !particleIndex.empty() ) return;
  particleIndex.clear();
  colourIndex.clear();
  colourIndex(0, tColinePtr());
  createParticles();
  createBeams();
}

void FxFxReader::reopen() {
  // If we didn't know how many events there were, we know now.
  if ( NEvents() <= 0 ) NEvents(position);
  ++reopened;
  // How large fraction of the events have we actually used? And how
  // large will we have used if we go through the file again?
  double frac = double(stats.attempts())/double(NEvents());
  if ( frac*double(reopened + 1)/double(reopened) > 1.0 &&
    NEvents() - stats.attempts() <
       generator()->N() - generator()->currentEventNumber() ) {
    if(theReOpenAllowed)
      generator()->logWarning(FxFxReopenWarning()
			      << "Reopening FxFxReader '" << name()
			      << "' after accessing " << stats.attempts() 
			      << " events out of "
			      << NEvents() << Exception::warning);
    else
      throw FxFxReopenWarning()
	<< "More events requested than available in FxFxReader "
	<< name() << Exception::runerror;
  }
  if ( cacheFile() ) {
    closeCacheFile();
    openReadCacheFile();
    if ( !uncacheEvent() ) Throw<FxFxReopenError>()
      << "Could not reopen FxFxReader '" << name()
      << "'." << Exception::runerror;
  } else {  
    close();
    open();
    if ( !readEvent() ) Throw<FxFxReopenError>()
      << "Could not reopen FxFxReader '" << name()
      << "'." << Exception::runerror;
  }
}

void FxFxReader::reset() {
  particleIndex.clear();
  colourIndex.clear();
  if ( theLastXComb ) theLastXComb->clean();
  theLastXComb = tXCombPtr();
}

bool FxFxReader::readEvent() {

  reset();

  if ( !doReadEvent() ) return false;

  // If we are just skipping event we do not need to reweight or do
  // anything fancy.
  if ( skipping ) return true;

  if ( cacheFile() && !scanning ) return true;

  // Reweight according to the re- and pre-weights objects in the
  // FxFxReader base class.
  lastweight = reweight();

  if ( !reweightPDF && !cutEarly() ) return true;
  // We should try to reweight the PDFs or make early cuts here.

  fillEvent();

  double x1 = incoming().first->momentum().plus()/
    beams().first->momentum().plus();

  if ( reweightPDF &&
       inPDF.first && outPDF.first && inPDF.first != outPDF.first ) {
    if ( hepeup.XPDWUP.first <= 0.0 )
      hepeup.XPDWUP.first =
	inPDF.first->xfx(inData.first, incoming().first->dataPtr(),
			 sqr(hepeup.SCALUP*GeV), x1);
    double xf = outPDF.first->xfx(inData.first, incoming().first->dataPtr(),
				  sqr(hepeup.SCALUP*GeV), x1);
    lastweight *= xf/hepeup.XPDWUP.first;
    hepeup.XPDWUP.first = xf;
  }

  double x2 = incoming().second->momentum().minus()/
    beams().second->momentum().minus();

  if ( reweightPDF &&
       inPDF.second && outPDF.second && inPDF.second != outPDF.second ) {
    if ( hepeup.XPDWUP.second <= 0.0 )
      hepeup.XPDWUP.second =
	inPDF.second->xfx(inData.second, incoming().second->dataPtr(),
			 sqr(hepeup.SCALUP*GeV), x2);
    double xf =
      outPDF.second->xfx(inData.second, incoming().second->dataPtr(),
			 sqr(hepeup.SCALUP*GeV), x2);
    lastweight *= xf/hepeup.XPDWUP.second;
    hepeup.XPDWUP.second = xf;
  }

  if ( cutEarly() ) {
    if ( !cuts().initSubProcess((incoming().first->momentum() +
				 incoming().second->momentum()).m2(),
				0.5*log(x1/x2)) ) lastweight = 0.0;
    tSubProPtr sub = getSubProcess();
    TmpTransform<tSubProPtr> tmp(sub, Utilities::getBoostToCM(sub->incoming()));
    if ( !cuts().passCuts(*sub) ) lastweight = 0.0;
  }

  return true;
}

double FxFxReader::getEvent() {
  if ( cacheFile() ) {
    if ( !uncacheEvent() ) reopen();
  } else {
    if ( !readEvent() ) reopen();
  }
  ++position;

  double max = maxWeights[hepeup.IDPRUP]*maxFactor;

  // cout << "maxFactor = " << maxFactor << " maxWeights[hepeup.IDPRUP] = " <<  maxWeights[hepeup.IDPRUP] << endl;
  // normalize all the weights to the max weight
  for(map<string,double>::iterator it=optionalWeights.begin();
      it!=optionalWeights.end();++it) {
    if(it->first!="ecom" && it->second!=-999 && it->second!=-111 && it->second!=-222 && it->second!=-333) {
      it->second = (max != 0.0) ? it->second/max : 0.0;
    }
  }
  //cout << "hepeup.XWGTUP = " << hepeup.XWGTUP << endl;
  //cout << "max = " << max << " " << eventWeight() << endl; 
  return max != 0.0? eventWeight()/max: 0.0;

}

void FxFxReader::skip(long n) {
  HoldFlag<> skipflag(skipping);
  while ( n-- ) getEvent();
}

double FxFxReader::reweight() {
  preweight = 1.0;
  if ( reweights.empty() && preweights.empty() &&
       !( CKKWHandler() && maxMultCKKW() > 0 && maxMultCKKW() > minMultCKKW() ) )
    return 1.0;
  fillEvent();
  getSubProcess();
  for ( int i = 0, N = preweights.size(); i < N; ++i ) {
    preweights[i]->setXComb(lastXCombPtr());
    preweight *= preweights[i]->weight();
  }
  double weight = preweight;
  for ( int i = 0, N = reweights.size(); i < N; ++i ) {
    reweights[i]->setXComb(lastXCombPtr());
    weight *= reweights[i]->weight();
  }

  // If we are caching events we do not want to do CKKW reweighting.
  if ( cacheFile() ) return weight;

  if ( CKKWHandler() && maxMultCKKW() > 0 && maxMultCKKW() > minMultCKKW() ) {
    CKKWHandler()->setXComb(lastXCombPtr());
    weight *= CKKWHandler()->reweightCKKW(minMultCKKW(), maxMultCKKW());
  }
  return weight;
}

bool FxFxReader::checkPartonBin() {

  // First find the positions of the incoming partons.
  pair< vector<int>, vector<int> > inc;
  for ( int i = 0; i < hepeup.NUP; ++i ) {
    if ( hepeup.ISTUP[i] == -9 ) {
      if ( inc.first.empty() ) inc.first.push_back(i);
      else if ( inc.second.empty() ) inc.second.push_back(i);
    }
    else if ( hepeup.ISTUP[i] == -1 ) {
      if ( inc.first.size() &&
	   hepeup.MOTHUP[i].first == inc.first.back() + 1 )
	inc.first.push_back(i);
      else if ( inc.second.size() &&
		hepeup.MOTHUP[i].first == inc.second.back() + 1 )
	inc.second.push_back(i);
      else if ( inc.first.empty() ) {
	inc.first.push_back(-1);
	inc.first.push_back(i);
      }
      else if ( inc.second.empty() ) {
	inc.second.push_back(-1);
	inc.second.push_back(i);
      }
      else if ( inc.first.size() <= inc.second.size() )
	inc.first.push_back(i);
      else
	inc.second.push_back(i);
    }
  }

  // Now store the corresponding id numbers
  pair< vector<long>, vector<long> > ids;
  ids.first.push_back(inc.first[0] < 0? heprup.IDBMUP.first:
		      hepeup.IDUP[inc.first[0]]);
  for ( int i = 1, N = inc.first.size(); i < N; ++i )
    ids.first.push_back(hepeup.IDUP[inc.first[i]]);
  ids.second.push_back(inc.second[0] < 0? heprup.IDBMUP.second:
		       hepeup.IDUP[inc.second[0]]);
  for ( int i = 1, N = inc.second.size(); i < N; ++i )
    ids.second.push_back(hepeup.IDUP[inc.second[i]]);

  // Find the correct pair of parton bins.
  PBPair pbp;
  for ( int i = 0, N = partonBins().size(); i < N; ++i ) {
    tcPBPtr curr = partonBins()[i].first;
    int icurr = inc.first.size() - 1;
    while ( curr && icurr >= 0 ) {
      if ( curr->parton()->id () != ids.first[icurr] ) break;
      curr = curr->incoming();
      --icurr;
    }
    if(!(!partonBins()[i].first->incoming() &&
	 !partonBins()[i].first->particle() &&  
	 partonBins()[i].first->parton()->id () == ids.first[0] &&
	 ( inc.first.size()==1 ||
	   (inc.first.size()==2 && ids.first[0]==ids.first[1]))) &&
       ( curr || icurr >= 0 ) ) continue;

    curr = partonBins()[i].second;
    icurr = inc.second.size() - 1;
    while ( curr && icurr >= 0 ) {
      if ( curr->parton()->id () != ids.second[icurr] ) break;
      curr = curr->incoming();
      --icurr;
    }
    if(!(!partonBins()[i].second->incoming() &&
	 !partonBins()[i].second->particle() &&  
	 partonBins()[i].second->parton()->id () == ids.second[0] &&
	 ( inc.second.size()==1 ||
	   (inc.second.size()==2 && ids.second[0]==ids.second[1]))) &&
       ( curr || icurr >= 0 ) ) continue;

    pbp = partonBins()[i];
  }

  // If we are only checking we return here.
  return ( pbp.first && pbp.second );

}

namespace {
  bool recursionNotNull(tcPBPtr bin, tcPPtr p) {
    while ( bin && p ) {
      if ( p->dataPtr() != bin->parton() ) break;
      bin = bin->incoming();
      p = p->parents().size()? p->parents()[0]: tPPtr();
    }
    return bin || p;
  }
}


tcPBPair FxFxReader::createPartonBinInstances() {
  tcPBPair sel;
  for ( int i = 0, N = partonBins().size(); i < N; ++i ) {
    tcPBPtr bin = partonBins()[i].first;
    tcPPtr p = incoming().first;
    if ( recursionNotNull(bin,p) ) continue;
    bin = partonBins()[i].second;
    p = incoming().second;
    if ( recursionNotNull(bin,p) ) continue;
    sel = partonBins()[i];
    break;
  }
  if ( !sel.first || !sel.second ) Throw<FxFxInconsistencyError>()
    << "Could not find appropriate PartonBin objects for event produced by "
    << "FxFxReader '" << name() << "'." << Exception::runerror;

  Direction<0> dir(true);
  thePartonBinInstances.first =
    new_ptr(PartonBinInstance(incoming().first, sel.first,
			      -sqr(hepeup.SCALUP*GeV)));
  if ( thePartonBinInstances.first->xi() > 1.00001 ) {
    Throw<FxFxInconsistencyError>()
      << "Found an event with momentum fraction larger than unity (x1="
      << thePartonBinInstances.first->xi()
      << "). The event will be skipped." << Exception::warning;
    throw Veto();
  }
  dir.reverse();
  thePartonBinInstances.second =
    new_ptr(PartonBinInstance(incoming().second, sel.second,
			      -sqr(hepeup.SCALUP*GeV)));

  if ( thePartonBinInstances.second->xi() > 1.00001 ) {
    Throw<FxFxInconsistencyError>()
      << "Found an event with momentum fraction larger than unity (x2="
      << thePartonBinInstances.second->xi()
      << "). The event will be skipped." << Exception::warning;
    throw Veto();
  }
  return sel;

}
		     

void FxFxReader::createParticles() {
  theBeams = PPair();
  theIncoming = PPair();
  theOutgoing = PVector();
  theIntermediates = PVector();
  for ( int i = 0, N = hepeup.IDUP.size(); i < N; ++i ) {
    if ( !hepeup.IDUP[i] ) continue;
    Lorentz5Momentum mom(hepeup.PUP[i][0]*GeV, hepeup.PUP[i][1]*GeV,
			 hepeup.PUP[i][2]*GeV, hepeup.PUP[i][3]*GeV,
			 hepeup.PUP[i][4]*GeV);
    if(theMomentumTreatment == 1)      mom.rescaleEnergy();
    else if(theMomentumTreatment == 2) mom.rescaleMass();
    PDPtr pd = getParticleData(hepeup.IDUP[i]);
    if (!pd) {
      Throw<FxFxInitError>()
	<< "FxFxReader '" << name() << "' found unknown particle ID "
	<< hepeup.IDUP[i]
	<< " in Les Houches common block structure.\n"
	<< "You need to define the new particle in an input file.\n"
	<< Exception::runerror;
    }
    if ( ! pd->coloured() 
	 && ( hepeup.ICOLUP[i].first != 0 || hepeup.ICOLUP[i].second != 0 ) ) {
      Throw<FxFxInconsistencyError>()
	<< "FxFxReader " << name() << ": " << pd->PDGName() 
	<< " is not a coloured particle.\nIt should not have "
	<< "(anti-)colour lines " << hepeup.ICOLUP[i].first
	<< ' ' << hepeup.ICOLUP[i].second
	<< " set; the event file needs to be fixed."
	<< Exception::runerror;
    }
    PPtr p = pd->produceParticle(mom);
    if(hepeup.ICOLUP[i].first>=0 && hepeup.ICOLUP[i].second >=0) {
      tColinePtr c = colourIndex(hepeup.ICOLUP[i].first);
      if ( c ) c->addColoured(p);
      c = colourIndex(hepeup.ICOLUP[i].second);
      if ( c ) c->addAntiColoured(p);
    }
    else {
      tColinePtr c1 = colourIndex(abs(hepeup.ICOLUP[i].first ));
      tColinePtr c2 = colourIndex(abs(hepeup.ICOLUP[i].second));
      if(pd->hasColour()) {
	c1->addColouredIndexed(p,1);
	c2->addColouredIndexed(p,2);
      }
      else {
	c1->addAntiColouredIndexed(p,1);
	c2->addAntiColouredIndexed(p,2);
      }
    }
    particleIndex(i + 1, p);
    switch ( hepeup.ISTUP[i] ) {
    case -9:
      if ( !theBeams.first ) theBeams.first = p;
      else if ( !theBeams.second ) theBeams.second = p;
      else Throw<FxFxInconsistencyError>()
	<< "To many incoming beam particles in the FxFxReader '"
	<< name() << "'." << Exception::runerror;
      break;
    case -1:
      if ( !theIncoming.first ) theIncoming.first = p;
      else if ( !theIncoming.second ) theIncoming.second = p;
      else if ( particleIndex(theIncoming.first) == hepeup.MOTHUP[i].first )
	theIncoming.first = p;
      else if ( particleIndex(theIncoming.second) == hepeup.MOTHUP[i].first )
	theIncoming.second = p;
      else Throw<FxFxInconsistencyError>()
	<< "To many incoming particles to hard subprocess in the "
	<< "FxFxReader '"	<< name() << "'." << Exception::runerror;
      p->scale(sqr(hepeup.SCALUP*GeV));
      break;
    case 1:
      theOutgoing.push_back(p);
      p->scale(sqr(hepeup.SCALUP*GeV));
      break;
    case -2:
    case 2:
    case 3:
      theIntermediates.push_back(p);
      break;
    default:
      Throw<FxFxInconsistencyError>()
	<< "Unknown status code (" << hepeup.ISTUP[i]
	<< ") in the FxFxReader '" << name() << "'."
	<< Exception::runerror;
    }

    // value 9 is defined as "Unknown or unpolarized particles"
    double spinup = hepeup.SPINUP[i];
    if ( abs(spinup - 9) < 1.0e-3 )
      spinup = 0.;
    if ( spinup < -1. || spinup > 1. ) {
      Throw<FxFxInconsistencyError>()
	<< "Polarization must be between -1 and 1, not "
	<< spinup << " as found in the "
	<< "FxFx event file.\nThe event file needs to be fixed." 
	<< Exception::runerror;
    }
    if( theIncludeSpin 
	&& abs(pd->id()) == ParticleID::tauminus 
	&& spinup !=0) {
      if(pd->iSpin() == PDT::Spin1Half ) {
	vector<Helicity::SpinorWaveFunction> wave;
	Helicity::SpinorWaveFunction(wave,p,Helicity::outgoing,true);
	RhoDMatrix rho(pd->iSpin(),true);
	rho(0,0) = 0.5*(1.-spinup);
	rho(1,1) = 0.5*(1.+spinup);
	p->spinInfo()->rhoMatrix() = rho;
	p->spinInfo()->  DMatrix() = rho;
      }
    }
  }
  // check the colour flows, and if necessary create any sources/sinks
  // hard process
  // get the particles in the hard process
  PVector external;
  for ( int i = 0, N = hepeup.IDUP.size(); i < N; ++i ) {
    unsigned int moth; 
    switch ( hepeup.ISTUP[i] ) {
    case -1:
      external.push_back(particleIndex.find(i+1));
      break;
    case 1: case 2: case 3:
      moth = hepeup.MOTHUP[i].first;
      if(moth!=0 && (hepeup.ISTUP[moth]==-1||hepeup.ISTUP[moth]==-2||
		   hepeup.ISTUP[moth]==-9))
	external.push_back(particleIndex.find(i+1));
      moth = hepeup.MOTHUP[i].second;
      if(moth!=0 && (hepeup.ISTUP[moth]==-1||hepeup.ISTUP[moth]==-2||
		     hepeup.ISTUP[moth]==-9))
	external.push_back(particleIndex.find(i+1));
      break;
    case -2: case -9: default: 
      break;
    }
  }
  // check the incoming/outgoing lines match
  vector<tColinePtr> unMatchedColour,unMatchedAntiColour;
  for(unsigned int ix=0;ix<external.size();++ix) {
    vector<tcColinePtr> 
      col  = external[ix]->colourInfo()->    colourLines();
    vector<tcColinePtr> 
      anti = external[ix]->colourInfo()->antiColourLines();
    if(hepeup.ISTUP[particleIndex(external[ix])-1]<0)
      swap(col,anti);
    if(!col.empty()) {
      for(unsigned int ic1=0;ic1<col.size();++ic1) {
	bool matched=false;
	for(unsigned int iy=0;iy<external.size();++iy) {
	  vector<tcColinePtr> col2;
	  if(hepeup.ISTUP[particleIndex(external[iy])-1]<0) {
	    if(external[iy]->colourInfo()->colourLines().empty()) continue;
	    col2 = external[iy]->colourInfo()->colourLines();
	  } 
	  else if(hepeup.ISTUP[particleIndex(external[iy])-1]>0) {
	    if(external[iy]->colourInfo()->antiColourLines().empty()) continue;
	    col2 = external[iy]->colourInfo()->antiColourLines();
	  }
	  for(unsigned int ic2=0;ic2<col2.size();++ic2) {
	    if(col[ic1]==col2[ic2]) {
	      matched=true;
	      break;
	    }
	  }
	  if(matched) break;
	}
	if(!matched) unMatchedColour.push_back(const_ptr_cast<tColinePtr>(col[ic1]));
      }
    }
    if(!anti.empty()) {
      for(unsigned int ic1=0;ic1<col.size();++ic1) {
	bool matched=false;
	for(unsigned int iy=0;iy<external.size();++iy) {
	  vector<tcColinePtr> anti2;
	  if(hepeup.ISTUP[particleIndex(external[iy])-1]<0) {
	    if(external[iy]->colourInfo()->colourLines().empty()) continue;
	    anti2 = external[iy]->colourInfo()->antiColourLines();
	  } 
	  else if(hepeup.ISTUP[particleIndex(external[iy])-1]>0) {
	    if(external[iy]->colourInfo()->antiColourLines().empty()) continue;
	    anti2 = external[iy]->colourInfo()->colourLines();
	  }
	  for(unsigned int ic2=0;ic2<anti2.size();++ic2) {
	    if(col[ic1]==anti2[ic2]) {
	      matched=true;
	      break;
	    }
	  }
	  if(matched) break;
	}
	if(!matched) unMatchedAntiColour.push_back(const_ptr_cast<tColinePtr>(anti[ic1]));
      }
    }
  }
  
  // might have source/sink
  if( unMatchedColour.size() + unMatchedAntiColour.size() != 0) {
    if(unMatchedColour.size() == 3 ) {
      unMatchedColour[0]->setSourceNeighbours(unMatchedColour[1],
					      unMatchedColour[2]);
    }
    else if(unMatchedColour.size() != 0 && ThePEG_DEBUG_LEVEL) {
      Throw<FxFxInconsistencyError>()
 	<< "FxFxReader '" << name() << "' found inconsistent colour "
 	<< "flow in Les Houches common block structure for hard process.\n"
 	<< hepeup << Exception::runerror;
    }
    if(unMatchedAntiColour.size() == 3 ) {
      unMatchedAntiColour[0]->setSinkNeighbours(unMatchedAntiColour[1],
						unMatchedAntiColour[2]);
    }
    else if(unMatchedAntiColour.size() != 0 && ThePEG_DEBUG_LEVEL) {
      Throw<FxFxInconsistencyError>()
 	<< "FxFxReader '" << name() << "' found inconsistent colour "
 	<< "flow in Les Houches common block structure for hard process.\n"
 	<< hepeup << Exception::runerror;
    }
  }
  
  // any subsequent decays
  for ( int i = 0, N = hepeup.IDUP.size(); i < N; ++i ) {
    if(hepeup.ISTUP[i] !=2 && hepeup.ISTUP[i] !=3) continue;
    PVector external;
    external.push_back(particleIndex.find(i+1));
    for ( int j = 0; j < N; ++j ) {
      if(hepeup.MOTHUP[j].first==i+1||  hepeup.MOTHUP[j].second==i+1)
	external.push_back(particleIndex.find(j+1));
    }
    // check the incoming/outgoing lines match
    vector<tColinePtr> unMatchedColour,unMatchedAntiColour;
    for(unsigned int ix=0;ix<external.size();++ix) {
      vector<tcColinePtr> 
	col  = external[ix]->colourInfo()->    colourLines();
      vector<tcColinePtr> 
	anti = external[ix]->colourInfo()->antiColourLines();
      if(ix==0) swap(col,anti);
      if(!col.empty()) {
	for(unsigned int ic1=0;ic1<col.size();++ic1) {
	  bool matched=false;
	  for(unsigned int iy=0;iy<external.size();++iy) {
	    if(iy==ix) continue;
	    vector<tcColinePtr> col2;
	    if(iy==0) {
	      if(external[iy]->colourInfo()->colourLines().empty()) continue;
	      col2 = external[iy]->colourInfo()->colourLines();
	    } 
	    else {
	      if(external[iy]->colourInfo()->antiColourLines().empty()) continue;
	      col2 = external[iy]->colourInfo()->antiColourLines();
	    }
	    for(unsigned int ic2=0;ic2<col2.size();++ic2) {
	      if(col[ic1]==col2[ic2]) {
		matched=true;
		break;
	      }
	    }
	    if(matched) break;
	  }
	  if(!matched) unMatchedColour.push_back(const_ptr_cast<tColinePtr>(col[ic1]));
	}
      }
      if(!anti.empty()) {
	for(unsigned int ic1=0;ic1<anti.size();++ic1) {
	  bool matched=false;
	  for(unsigned int iy=0;iy<external.size();++iy) {
	    if(iy==ix) continue;
	    vector<tcColinePtr> anti2;
	    if(iy==0) {
	      if(external[iy]->colourInfo()->antiColourLines().empty()) continue;
	      anti2 = external[iy]->colourInfo()->antiColourLines();
	    } 
	    else {
	      if(external[iy]->colourInfo()->colourLines().empty()) continue;
	      anti2 = external[iy]->colourInfo()->colourLines();
	    }
	    for(unsigned int ic2=0;ic2<anti2.size();++ic2) {
	      if(anti[ic1]==anti2[ic2]) {
		matched=true;
		break;
	      }
	    }
	    if(matched) break;
	  }
	  if(!matched) unMatchedAntiColour.push_back(const_ptr_cast<tColinePtr>(anti[ic1]));
	}
      }
    }
    // might have source/sink
    if( unMatchedColour.size() + unMatchedAntiColour.size() != 0) {
      if(unMatchedColour.size() == 3 ) {
	unMatchedColour[0]->setSourceNeighbours(unMatchedColour[1],
						unMatchedColour[2]);
      }
      else if(unMatchedColour.size() != 0 && ThePEG_DEBUG_LEVEL) {
	Throw<FxFxInconsistencyError>()
	  << "FxFxReader '" << name() << "' found inconsistent colour "
	  << "flow in Les Houches common block structure for decay of \n"
	  << *external[0] << "\n"
	  << hepeup << Exception::runerror;
      }
      if(unMatchedAntiColour.size() == 3 ) {
	unMatchedAntiColour[0]->setSinkNeighbours(unMatchedAntiColour[1],
						  unMatchedAntiColour[2]);
      }
      else if(unMatchedAntiColour.size() != 0 && ThePEG_DEBUG_LEVEL) {
	Throw<FxFxInconsistencyError>()
	  << "FxFxReader '" << name() << "' found inconsistent colour "
	  << "flow in Les Houches common block structure for decay of\n"
	  << *external[0] << "\n"
	  << hepeup << Exception::runerror;
      }
    }
  }
}

void FxFxReader::createBeams() {

  if ( !theBeams.first && dynamic_ptr_cast<Ptr<NoPDF>::tcp>(inPDF.first) ) {
    theBeams.first = theIncoming.first;
  }
  else if ( !theBeams.first ) {
    theBeams.first = getParticleData(heprup.IDBMUP.first)->produceParticle();
    double m = theBeams.first->mass()/GeV;
    theBeams.first->set5Momentum
      (Lorentz5Momentum(ZERO, ZERO,
			sqrt(sqr(heprup.EBMUP.first) - sqr(m))*GeV,
			heprup.EBMUP.first*GeV, m*GeV));
    hepeup.IDUP.push_back(heprup.IDBMUP.first);
    hepeup.ISTUP.push_back(-9);
    hepeup.MOTHUP.push_back(make_pair(0, 0));
    hepeup.ICOLUP.push_back(make_pair(0, 0));
    hepeup.VTIMUP.push_back(0.0);
    hepeup.SPINUP.push_back(0.0);
    particleIndex(hepeup.IDUP.size(), theBeams.first);
    hepeup.MOTHUP[particleIndex(theIncoming.first) - 1].first =
      hepeup.IDUP.size();
  }
  if ( !theBeams.second && dynamic_ptr_cast<Ptr<NoPDF>::tcp>(inPDF.second) ) {
    theBeams.second = theIncoming.second;
  }
  else if ( !theBeams.second ) {
    theBeams.second = getParticleData(heprup.IDBMUP.second)->produceParticle();
    double m = theBeams.second->mass()/GeV;
    theBeams.second->set5Momentum
      (Lorentz5Momentum(ZERO, ZERO,
			-sqrt(sqr(heprup.EBMUP.second) - sqr(m))*GeV,
			heprup.EBMUP.second*GeV, m*GeV));
    hepeup.IDUP.push_back(heprup.IDBMUP.second);
    hepeup.ISTUP.push_back(-9);
    hepeup.MOTHUP.push_back(make_pair(0, 0));
    hepeup.ICOLUP.push_back(make_pair(0, 0));
    hepeup.VTIMUP.push_back(0.0);
    hepeup.SPINUP.push_back(0.0);
    particleIndex(hepeup.IDUP.size(), theBeams.second);
    hepeup.MOTHUP[particleIndex(theIncoming.second) - 1].first =
      hepeup.IDUP.size();
  }
}

void FxFxReader::connectMothers() {
  const ObjectIndexer<long,Particle> & pi = particleIndex;
  for ( int i = 0, N = hepeup.IDUP.size(); i < N; ++i ) {
    if ( pi(hepeup.MOTHUP[i].first) ) 
      pi(hepeup.MOTHUP[i].first)->addChild(pi(i + 1));
    if ( pi(hepeup.MOTHUP[i].second) 
	 && hepeup.MOTHUP[i].second != hepeup.MOTHUP[i].first ) 
      pi(hepeup.MOTHUP[i].second)->addChild(pi(i + 1));
  }
}

void FxFxReader::openReadCacheFile() {
  if ( cacheFile() ) closeCacheFile();
  cacheFile().open(cacheFileName(), "r");
  position = 0;
}

void FxFxReader::openWriteCacheFile() {
  if ( cacheFile() ) closeCacheFile();
  cacheFile().open(cacheFileName(), "w");
}

void FxFxReader::closeCacheFile() {
  cacheFile().close();
}

void FxFxReader::cacheEvent() const {
  static vector<char> buff;
  cacheFile().write(&hepeup.NUP, sizeof(hepeup.NUP));
  buff.resize(eventSize(hepeup.NUP));
  char * pos = &buff[0];
  pos = mwrite(pos, hepeup.IDPRUP);
  pos = mwrite(pos, hepeup.XWGTUP);
  pos = mwrite(pos, hepeup.XPDWUP);
  pos = mwrite(pos, hepeup.SCALUP);
  pos = mwrite(pos, hepeup.AQEDUP);
  pos = mwrite(pos, hepeup.AQCDUP);
  pos = mwrite(pos, hepeup.IDUP[0], hepeup.NUP);
  pos = mwrite(pos, hepeup.ISTUP[0], hepeup.NUP);
  pos = mwrite(pos, hepeup.MOTHUP[0], hepeup.NUP);
  pos = mwrite(pos, hepeup.ICOLUP[0], hepeup.NUP);
  for ( int i = 0; i < hepeup.NUP; ++i )
    pos = mwrite(pos, hepeup.PUP[i][0], 5);
  pos = mwrite(pos, hepeup.VTIMUP[0], hepeup.NUP);
  pos = mwrite(pos, hepeup.SPINUP[0], hepeup.NUP);
  pos = mwrite(pos, lastweight);
  pos = mwrite(pos, optionalWeights);
  for(size_t ff = 0; ff < optionalWeightsNames.size(); ff++) {
    pos = mwrite(pos, optionalWeightsNames[ff]);
  }
  /*  for(int f = 0; f < optionalWeightsNames.size(); f++) {
    cout << "optionalWeightsNames = " << optionalWeightsNames[f] << endl;
    }*/
  pos = mwrite(pos, optionalnpLO);
  pos = mwrite(pos, optionalnpNLO);
  pos = mwrite(pos, preweight);
  cacheFile().write(&buff[0], buff.size(), 1);
}

bool FxFxReader::uncacheEvent() {
  reset();
  static vector<char> buff;
  if ( cacheFile().read(&hepeup.NUP, sizeof(hepeup.NUP)) != 1 )
    return false;
  buff.resize(eventSize(hepeup.NUP));
  if ( cacheFile().read(&buff[0], buff.size()) != 1 ) return false;
  const char * pos = &buff[0];
  pos = mread(pos, hepeup.IDPRUP);
  pos = mread(pos, hepeup.XWGTUP);
  pos = mread(pos, hepeup.XPDWUP);
  pos = mread(pos, hepeup.SCALUP);
  pos = mread(pos, hepeup.AQEDUP);
  pos = mread(pos, hepeup.AQCDUP);
  hepeup.IDUP.resize(hepeup.NUP);
  pos = mread(pos, hepeup.IDUP[0], hepeup.NUP);
  hepeup.ISTUP.resize(hepeup.NUP);
  pos = mread(pos, hepeup.ISTUP[0], hepeup.NUP);
  hepeup.MOTHUP.resize(hepeup.NUP);
  pos = mread(pos, hepeup.MOTHUP[0], hepeup.NUP);
  hepeup.ICOLUP.resize(hepeup.NUP);
  pos = mread(pos, hepeup.ICOLUP[0], hepeup.NUP);
  hepeup.PUP.resize(hepeup.NUP, vector<double>(5));
  for ( int i = 0; i < hepeup.NUP; ++i ) 
    pos = mread(pos, hepeup.PUP[i][0], 5);
  hepeup.VTIMUP.resize(hepeup.NUP);
  pos = mread(pos, hepeup.VTIMUP[0], hepeup.NUP);
  hepeup.SPINUP.resize(hepeup.NUP);
  pos = mread(pos, hepeup.SPINUP[0], hepeup.NUP);
  pos = mread(pos, lastweight);
  pos = mread(pos, optionalWeights);
  for(size_t ff = 0; ff < optionalWeightsNames.size(); ff++) {
    pos = mread(pos, optionalWeightsNames[ff]);
  }
  pos = mread(pos, optionalnpLO);
  pos = mread(pos, optionalnpNLO);
  pos = mread(pos, preweight);

  // If we are skipping, we do not have to do anything else.
  if ( skipping ) return true;

  if ( CKKWHandler() && maxMultCKKW() > 0 && maxMultCKKW() > minMultCKKW() ) {
    // The cached event has not been submitted to CKKW reweighting, so
    // we do that now.
    fillEvent();
    getSubProcess();
    CKKWHandler()->setXComb(lastXCombPtr());
    lastweight *= CKKWHandler()->reweightCKKW(minMultCKKW(), maxMultCKKW());
  }
  return true;
}

void FxFxReader::persistentOutput(PersistentOStream & os) const {
  os << heprup.IDBMUP << heprup.EBMUP << heprup.PDFGUP << heprup.PDFSUP
     << heprup.IDWTUP << heprup.NPRUP << heprup.XSECUP << heprup.XERRUP
     << heprup.XMAXUP << heprup.LPRUP << hepeup.NUP << hepeup.IDPRUP
     << hepeup.XWGTUP << hepeup.XPDWUP << hepeup.SCALUP << hepeup.AQEDUP
     << hepeup.AQCDUP << hepeup.IDUP << hepeup.ISTUP << hepeup.MOTHUP
     << hepeup.ICOLUP << hepeup.PUP << hepeup.VTIMUP << hepeup.SPINUP
     << inData << inPDF << outPDF << thePartonExtractor << theCKKW
     << thePartonBins << theXCombs << theCuts << theNEvents << position
     << reopened << theMaxScan << isActive
     << theCacheFileName << doCutEarly << stats << statmap
     << thePartonBinInstances
     << theBeams << theIncoming << theOutgoing << theIntermediates
     << reweights << preweights << preweight << reweightPDF << doInitPDFs
     << theLastXComb << theMaxMultCKKW << theMinMultCKKW << lastweight << optionalWeights << optionalnpLO << optionalnpNLO
     << maxFactor << ounit(weightScale, picobarn) << xSecWeights << maxWeights
     << theMomentumTreatment << useWeightWarnings << theReOpenAllowed
     << theIncludeSpin;
}

void FxFxReader::persistentInput(PersistentIStream & is, int) {
  if ( cacheFile() ) closeCacheFile();
  is >> heprup.IDBMUP >> heprup.EBMUP >> heprup.PDFGUP >> heprup.PDFSUP
     >> heprup.IDWTUP >> heprup.NPRUP >> heprup.XSECUP >> heprup.XERRUP
     >> heprup.XMAXUP >> heprup.LPRUP >> hepeup.NUP >> hepeup.IDPRUP
     >> hepeup.XWGTUP >> hepeup.XPDWUP >> hepeup.SCALUP >> hepeup.AQEDUP
     >> hepeup.AQCDUP >> hepeup.IDUP >> hepeup.ISTUP >> hepeup.MOTHUP
     >> hepeup.ICOLUP >> hepeup.PUP >> hepeup.VTIMUP >> hepeup.SPINUP
     >> inData >> inPDF >> outPDF >> thePartonExtractor >> theCKKW
     >> thePartonBins >> theXCombs >> theCuts >> theNEvents >> position
     >> reopened >> theMaxScan >> isActive
     >> theCacheFileName >> doCutEarly >> stats >> statmap
     >> thePartonBinInstances
     >> theBeams >> theIncoming >> theOutgoing >> theIntermediates
     >> reweights >> preweights >> preweight >> reweightPDF >> doInitPDFs
     >> theLastXComb >> theMaxMultCKKW >> theMinMultCKKW >> lastweight >> optionalWeights >> optionalnpLO >> optionalnpNLO
     >> maxFactor >> iunit(weightScale, picobarn) >> xSecWeights >> maxWeights
     >> theMomentumTreatment >> useWeightWarnings >> theReOpenAllowed
     >> theIncludeSpin;
}

AbstractClassDescription<FxFxReader>
FxFxReader::initFxFxReader;
// Definition of the static class description member.

void FxFxReader::setBeamA(long id) { heprup.IDBMUP.first = id; }
long FxFxReader::getBeamA() const { return heprup.IDBMUP.first; }
void FxFxReader::setBeamB(long id) { heprup.IDBMUP.second = id; }
long FxFxReader::getBeamB() const { return heprup.IDBMUP.second; }
void FxFxReader::setEBeamA(Energy e) { heprup.EBMUP.first = e/GeV; }
Energy FxFxReader::getEBeamA() const { return heprup.EBMUP.first*GeV; }
void FxFxReader::setEBeamB(Energy e) { heprup.EBMUP.second = e/GeV; }
Energy FxFxReader::getEBeamB() const { return heprup.EBMUP.second*GeV; }
void FxFxReader::setPDFA(PDFPtr pdf) { inPDF.first = pdf; }
PDFPtr FxFxReader::getPDFA() const { return inPDF.first; }
void FxFxReader::setPDFB(PDFPtr pdf) { inPDF.second = pdf; }
PDFPtr FxFxReader::getPDFB() const { return inPDF.second; }

void FxFxReader::Init() {

  static ClassDocumentation<FxFxReader> documentation
    ("ThePEG::FxFxReader is an abstract base class to be used "
     "for objects which reads event files or streams from matrix element "
     "generators.");

  static Parameter<FxFxReader,long> interfaceBeamA
    ("BeamA",
     "The PDG id of the incoming particle along the positive z-axis. "
     "If zero the corresponding information is to be deduced from the "
     "event stream/file.",
     0, 0, 0, 0,
     true, false, false,
     &FxFxReader::setBeamA,
     &FxFxReader::getBeamA, 0, 0, 0);

  static Parameter<FxFxReader,long> interfaceBeamB
    ("BeamB",
     "The PDG id of the incoming particle along the negative z-axis. "
     "If zero the corresponding information is to be deduced from the "
     "event stream/file.",
     0, 0, 0, 0,
     true, false, false,
     &FxFxReader::setBeamB,
     &FxFxReader::getBeamB, 0, 0, 0);

  static Parameter<FxFxReader,Energy> interfaceEBeamA
    ("EBeamA",
     "The energy of the incoming particle along the positive z-axis. "
     "If zero the corresponding information is to be deduced from the "
     "event stream/file.",
     0, GeV, ZERO, ZERO, 1000000000.0*GeV,
     true, false, true,
     &FxFxReader::setEBeamA, &FxFxReader::getEBeamA, 0, 0, 0);

  static Parameter<FxFxReader,Energy> interfaceEBeamB
    ("EBeamB",
     "The energy of the incoming particle along the negative z-axis. "
     "If zero the corresponding information is to be deduced from the "
     "event stream/file.",
     0, GeV, ZERO, ZERO, 1000000000.0*GeV,
     true, false, true,
     &FxFxReader::setEBeamB, &FxFxReader::getEBeamB, 0, 0, 0);

  static Reference<FxFxReader,PDFBase> interfacePDFA
    ("PDFA",
     "The PDF used for incoming particle along the positive z-axis. "
     "If null the corresponding information is to be deduced from the "
     "event stream/file.",
     0, true, false, true, true, false,
     &FxFxReader::setPDFA, &FxFxReader::getPDFA, 0);

  static Reference<FxFxReader,PDFBase> interfacePDFB
    ("PDFB",
     "The PDF used for incoming particle along the negative z-axis. "
     "If null the corresponding information is to be deduced from the "
     "event stream/file.",
     0, true, false, true, true, false,
     &FxFxReader::setPDFB, &FxFxReader::getPDFB, 0);

  static Parameter<FxFxReader,long> interfaceMaxScan
    ("MaxScan",
     "The maximum number of events to scan to obtain information about "
     "processes and cross section in the intialization.",
     &FxFxReader::theMaxScan, -1, 0, 0,
     true, false, false);

  static Parameter<FxFxReader,string> interfaceCacheFileName
    ("CacheFileName",
     "Name of file used to cache the events form the reader in a fast-readable "
     "form. If empty, no cache file will be generated.",
     &FxFxReader::theCacheFileName, "",
     true, false);
  interfaceCacheFileName.fileType();

  static Switch<FxFxReader,bool> interfaceCutEarly
    ("CutEarly",
     "Determines whether to apply cuts to events before converting to "
     "ThePEG format.",
     &FxFxReader::doCutEarly, true, true, false);
  static SwitchOption interfaceCutEarlyYes
    (interfaceCutEarly,
     "Yes",
     "Event are cut before converted.",
     true);
  static SwitchOption interfaceCutEarlyNo
    (interfaceCutEarly,
     "No",
     "Events are not cut before converted.",
     false);

  static Reference<FxFxReader,PartonExtractor> interfacePartonExtractor
    ("PartonExtractor",
     "The PartonExtractor object used to construct remnants. If no object is "
     "provided the FxFxEventHandler object must provide one instead.",
     &FxFxReader::thePartonExtractor, true, false, true, true, false);


  static Reference<FxFxReader,Cuts> interfaceCuts
    ("Cuts",
     "The Cuts object to be used for this reader. Note that these "
     "must not be looser cuts than those used in the actual generation. "
     "If no object is provided the FxFxEventHandler object must "
     "provide one instead.",
     &FxFxReader::theCuts, true, false, true, true, false);

  static RefVector<FxFxReader,ReweightBase> interfaceReweights
    ("Reweights",
     "A list of ThePEG::ReweightBase objects to modify this the weight of "
     "this reader.",
     &FxFxReader::reweights, 0, false, false, true, false);

  static RefVector<FxFxReader,ReweightBase> interfacePreweights
    ("Preweights",
     "A list of ThePEG::ReweightBase objects to bias the phase space for this "
     "reader without influencing the actual cross section.",
     &FxFxReader::preweights, 0, false, false, true, false);

  static Switch<FxFxReader,bool> interfaceReweightPDF
    ("ReweightPDF",
     "If the PDFs used in the generation for this reader is different "
     "from the ones assumed by the associated PartonExtractor object, "
     "should the events be reweighted to fit the latter?",
     &FxFxReader::reweightPDF, false, true, false);
  static SwitchOption interfaceReweightPDFNo
    (interfaceReweightPDF, "No", "The event weights are kept as they are.",
     false);
  static SwitchOption interfaceReweightPDFYes
    (interfaceReweightPDF,
     "Yes", "The events are reweighted.", true);

  static Switch<FxFxReader,bool> interfaceInitPDFs
    ("InitPDFs",
     "If no PDFs were specified in <interface>PDFA</interface> or "
     "<interface>PDFB</interface>for this reader, try to extract the "
     "information from the event file and assign the relevant PDFBase"
     "objects when the reader is initialized.",
     &FxFxReader::doInitPDFs, false, true, false);
  static SwitchOption interfaceInitPDFsYes
    (interfaceInitPDFs,
     "Yes",
     "Extract PDFs during initialization.",
     true);
  static SwitchOption interfaceInitPDFsNo
    (interfaceInitPDFs,
     "No",
     "Do not extract PDFs during initialization.",
     false);

  static Parameter<FxFxReader,int> interfaceMaxMultCKKW
    ("MaxMultCKKW",
     "If this reader is to be used (possibly together with others) for CKKW-"
     "reweighting and veto, this should give the multiplicity of outgoing "
     "particles in the highest multiplicity matrix element in the group. "
     "If set to zero, no CKKW procedure should be applied.",
     &FxFxReader::theMaxMultCKKW, 0, 0, 0,
     true, false, Interface::lowerlim);

  static Parameter<FxFxReader,int> interfaceMinMultCKKW
    ("MinMultCKKW",
     "If this reader is to be used (possibly together with others) for CKKW-"
     "reweighting and veto, this should give the multiplicity of outgoing "
     "particles in the lowest multiplicity matrix element in the group. If "
     "larger or equal to <interface>MaxMultCKKW</interface>, no CKKW "
     "procedure should be applied.",
     &FxFxReader::theMinMultCKKW, 0, 0, 0,
     true, false, Interface::lowerlim);


  static Switch<FxFxReader,unsigned int> interfaceMomentumTreatment
    ("MomentumTreatment",
     "Treatment of the momenta supplied by the interface",
     &FxFxReader::theMomentumTreatment, 0, false, false);
  static SwitchOption interfaceMomentumTreatmentAccept
    (interfaceMomentumTreatment,
     "Accept",
     "Just accept the momenta given",
     0);
  static SwitchOption interfaceMomentumTreatmentRescaleEnergy
    (interfaceMomentumTreatment,
     "RescaleEnergy",
     "Rescale the energy supplied so it is consistent with the mass",
     1);
  static SwitchOption interfaceMomentumTreatmentRescaleMass
    (interfaceMomentumTreatment,
     "RescaleMass",
     "Rescale the mass supplied so it is consistent with the"
     " energy and momentum",
     2);


  static Switch<FxFxReader,bool> interfaceWeightWarnings
    ("WeightWarnings",
     "Determines if warnings about possible weight incompatibilities should "
     "be issued when this reader is initialized.",
     &FxFxReader::useWeightWarnings, true, true, false);
  static SwitchOption interfaceWeightWarningsWarnAboutWeights
    (interfaceWeightWarnings,
     "WarnAboutWeights",
     "Warn about possible incompatibilities with the weight option in the "
     "Les Houches common block and the requested weight treatment.",
     true);
  static SwitchOption interfaceWeightWarningsDontWarnAboutWeights
    (interfaceWeightWarnings,
     "DontWarnAboutWeights",
     "Do not warn about possible incompatibilities with the weight option "
     "in the Les Houches common block and the requested weight treatment.",
     false);

  static Switch<FxFxReader,bool> interfaceAllowedTopReOpen
    ("AllowedToReOpen",
     "Can the file be reopened if more events are requested than the file contains?",
     &FxFxReader::theReOpenAllowed, true, false, false);
  static SwitchOption interfaceAllowedTopReOpenYes
    (interfaceAllowedTopReOpen,
     "Yes",
     "Allowed to reopen the file",
     true);
  static SwitchOption interfaceAllowedTopReOpenNo
    (interfaceAllowedTopReOpen,
     "No",
     "Not allowed to reopen the file",
     false);

  static Switch<FxFxReader,bool> interfaceIncludeSpin
    ("IncludeSpin",
     "Use the spin information present in the event file, for tau leptons"
     " only as this is the only case which makes any sense",
     &FxFxReader::theIncludeSpin, true, false, false);
  static SwitchOption interfaceIncludeSpinYes
    (interfaceIncludeSpin,
     "Yes",
     "Use the spin information",
     true);
  static SwitchOption interfaceIncludeSpinNo
    (interfaceIncludeSpin,
     "No",
     "Don't use the spin information",
     false);



  interfaceCuts.rank(8);
  interfacePartonExtractor.rank(7);
  interfaceBeamA.rank(5);
  interfaceBeamB.rank(4);
  interfaceEBeamA.rank(3);
  interfaceEBeamB.rank(2);
  interfaceMaxMultCKKW.rank(1.5);
  interfaceMinMultCKKW.rank(1.0);

  interfaceBeamA.setHasDefault(false);
  interfaceBeamB.setHasDefault(false);
  interfaceEBeamA.setHasDefault(false);
  interfaceEBeamB.setHasDefault(false);
  interfaceMaxMultCKKW.setHasDefault(false);
  interfaceMinMultCKKW.setHasDefault(false);

}

namespace ThePEG {

ostream & operator<<(ostream & os, const HEPEUP & h) {
  os << "<event>\n"
     << " " << setw(4) << h.NUP
     << " " << setw(6) << h.IDPRUP
     << " " << setw(14) << h.XWGTUP
     << " " << setw(14) << h.SCALUP
     << " " << setw(14) << h.AQEDUP
     << " " << setw(14) << h.AQCDUP << "\n";

  for ( int i = 0; i < h.NUP; ++i )
    os << " " << setw(8) << h.IDUP[i]
       << " " << setw(2) << h.ISTUP[i]
       << " " << setw(4) << h.MOTHUP[i].first
       << " " << setw(4) << h.MOTHUP[i].second
       << " " << setw(4) << h.ICOLUP[i].first
       << " " << setw(4) << h.ICOLUP[i].second
       << " " << setw(14) << h.PUP[i][0]
       << " " << setw(14) << h.PUP[i][1]
       << " " << setw(14) << h.PUP[i][2]
       << " " << setw(14) << h.PUP[i][3]
       << " " << setw(14) << h.PUP[i][4]
       << " " << setw(1) << h.VTIMUP[i]
       << " " << setw(1) << h.SPINUP[i] << std::endl;
  os << "</event>" << std::endl;
  return os;
}

}

