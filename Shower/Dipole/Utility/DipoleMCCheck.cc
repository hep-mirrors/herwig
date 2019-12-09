// -*- C++ -*-
//
// DipoleMCCheck.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DipoleMCCheck class.
//

#include "DipoleMCCheck.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

DipoleMCCheck::DipoleMCCheck() 
  : HandlerBase(),
    theHardPtBins(10),
    theEmitterXBins(5), theSpectatorXBins(5),
    thePtBins(100), theZBins(100) {
}

DipoleMCCheck::~DipoleMCCheck() {}

IBPtr DipoleMCCheck::clone() const {
  return new_ptr(*this);
}

IBPtr DipoleMCCheck::fullclone() const {
  return new_ptr(*this);
}

vector<double> DipoleMCCheck::makeLogBins(double xlow, double xup, unsigned int n) const {

  vector<double> res;

  double c = log10(xup/xlow) / (n-1.);

  for ( unsigned int k = 0; k < n; ++k )
    res.push_back(xlow*pow(10.0,k*c));

  return res;

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void DipoleMCCheck::doinitrun() {
  HandlerBase::doinitrun();

  vector<double> ptbins;
  double w = 0.5/theHardPtBins;

  for ( unsigned int k = 1; k <= theHardPtBins; ++k )
    ptbins.push_back(k*w);

  vector<double> xebins;
  if ( theEmitterXBins > 1 )
    xebins = makeLogBins(1e-7,1.0,theEmitterXBins);
  else
    xebins.push_back(1.0);

  vector<double> xsbins;
  if ( theSpectatorXBins > 1 )
    xsbins = makeLogBins(1e-7,1.0,theSpectatorXBins);
  else
    xsbins.push_back(1.0);

  for ( vector<double>::const_iterator xeit = xebins.begin();
	xeit != xebins.end(); ++xeit ) {
    map<double,
      map<double,
      pair<Ptr<Histogram>::ptr,Ptr<Histogram>::ptr>
      >
      > xebin;
    for ( vector<double>::const_iterator xsit = xsbins.begin();
	  xsit != xsbins.end(); ++xsit ) {
      map<double,
	pair<Ptr<Histogram>::ptr,Ptr<Histogram>::ptr>
	> xsbin;
      for ( vector<double>::const_iterator ptit = ptbins.begin();
	    ptit != ptbins.end(); ++ptit ) {
	pair<Ptr<Histogram>::ptr,Ptr<Histogram>::ptr> ptbin
	  (new_ptr(Histogram(0.0,0.5,thePtBins)),
	   new_ptr(Histogram(0.0,1.0,theZBins)));
	xsbin[*ptit] = ptbin;
      }
      xebin[*xsit] = xsbin;
    }
    histoMap[*xeit] = xebin;
  }

}

void DipoleMCCheck::dofinish() {
  HandlerBase::dofinish();

  map<double,
    map<double,
    map<double,
    pair<Ptr<Histogram>::ptr,Ptr<Histogram>::ptr>
    >
    >
    >::iterator xeit;

  map<double,
    map<double,
    pair<Ptr<Histogram>::ptr,Ptr<Histogram>::ptr>
    >
    >::iterator xsit;

  map<double,
    pair<Ptr<Histogram>::ptr,Ptr<Histogram>::ptr>
    >::iterator ptit;

  double xelow = 0.0;
  double xeup = 0.0;

  for ( xeit = histoMap.begin(); xeit != histoMap.end(); ++xeit ) {

    xeup = xelow + xeit->first;

    double xslow = 0.0;
    double xsup = 0.0;

    for ( xsit = xeit->second.begin(); xsit != xeit->second.end(); ++xsit ) {

      xsup = xslow + xsit->first;

      // open files here

      ostringstream ptFileName;
      ptFileName << name() << "_pt_"
		 << xelow << "_" << xeup << "_"
		 << xslow << "_" << xsup << ".dat";

      ofstream ptFile(ptFileName.str().c_str());

      ostringstream zFileName;
      zFileName << name() << "_z_"
		<< xelow << "_" << xeup << "_"
		<< xslow << "_" << xsup << ".dat";

      ofstream zFile(zFileName.str().c_str());

      double ptlow = 0.0;
      double ptup = 0.0;

      for ( ptit = xsit->second.begin(); ptit != xsit->second.end(); ++ptit ) {

	ptup = ptlow + ptit->first;

	// dump histos here

	ptFile << "#\n# " << ptlow << " < \\kappa < " << ptup << "\n#\n";
	zFile << "#\n# " << ptlow << " < \\kappa < " << ptup << "\n#\n";

	ptit->second.first->simpleOutput(ptFile,false);
	ptit->second.second->simpleOutput(zFile,false);

	ptFile << "\n\n\n";
	zFile << "\n\n\n";

	ptlow = ptup;

      }

      xslow = xsup;

    }

    xelow = xeup;

  }

}

void DipoleMCCheck::book(double xe,double xs,
			 Energy dScale,
			 Energy hardPt,
			 Energy pt, double z,
			 double weight) {

  map<double,
    map<double,
    map<double,
    pair<Ptr<Histogram>::ptr,Ptr<Histogram>::ptr>
    >
    >
    >::iterator xeit;

  map<double,
    map<double,
    pair<Ptr<Histogram>::ptr,Ptr<Histogram>::ptr>
    >
    >::iterator xsit;

  map<double,
    pair<Ptr<Histogram>::ptr,Ptr<Histogram>::ptr>
    >::iterator ptit;

  if ( theEmitterXBins == 1 || xe >= 1. )
    xeit = --histoMap.end();
  else
    xeit = histoMap.upper_bound(xe);

  if ( theSpectatorXBins == 1 || xs >= 1. )
    xsit = --(xeit->second.end());
  else
    xsit = xeit->second.upper_bound(xs);

  if ( theHardPtBins == 1 || hardPt/dScale >= 0.5 )
    ptit = --(xsit->second.end());
  else
    ptit = xsit->second.upper_bound(hardPt/dScale);

  ptit->second.first->addWeighted(pt/dScale,weight);
  ptit->second.second->addWeighted(z,weight);

}

void DipoleMCCheck::persistentOutput(PersistentOStream & os) const {
  os << theHardPtBins << theEmitterXBins
     << theSpectatorXBins << thePtBins
     << theZBins;
}

void DipoleMCCheck::persistentInput(PersistentIStream & is, int) {
  is >> theHardPtBins >> theEmitterXBins
     >> theSpectatorXBins >> thePtBins
     >> theZBins;
}

ClassDescription<DipoleMCCheck> DipoleMCCheck::initDipoleMCCheck;
// Definition of the static class description member.

void DipoleMCCheck::Init() {

  static ClassDocumentation<DipoleMCCheck> documentation
    ("DipoleMCCheck is used to perform checks for "
     "the dipole shower.");


  static Parameter<DipoleMCCheck,unsigned int> interfaceHardPtBins
    ("HardPtBins",
     "HardPtBins",
     &DipoleMCCheck::theHardPtBins, 10, 1, 0,
     false, false, Interface::lowerlim);

  static Parameter<DipoleMCCheck,unsigned int> interfaceEmitterXBins
    ("EmitterXBins",
     "EmitterXBins",
     &DipoleMCCheck::theEmitterXBins, 5, 1, 0,
     false, false, Interface::lowerlim);

  static Parameter<DipoleMCCheck,unsigned int> interfaceSpectatorXBins
    ("SpectatorXBins",
     "SpectatorXBins",
     &DipoleMCCheck::theSpectatorXBins, 5, 1, 0,
     false, false, Interface::lowerlim);

  static Parameter<DipoleMCCheck,unsigned int> interfacePtBins
    ("PtBins",
     "PtBins",
     &DipoleMCCheck::thePtBins, 100, 1, 0,
     false, false, Interface::lowerlim);

  static Parameter<DipoleMCCheck,unsigned int> interfaceZBins
    ("ZBins",
     "ZBins",
     &DipoleMCCheck::theZBins, 100, 1, 0,
     false, false, Interface::lowerlim);

}

