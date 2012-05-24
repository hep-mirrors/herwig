// -*- C++ -*-
//
// DipoleSplittingInfo.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DipoleIndex and DipoleSplittingInfo classes.
//

#include "DipoleSplittingInfo.h"
#include "ThePEG/Handlers/StandardXComb.h"

#include <iterator>

using std::ostream_iterator;

using namespace Herwig;

DipoleIndex::DipoleIndex()
  : theInitialStateEmitter(false), theInitialStateSpectator(false) {}

DipoleIndex::DipoleIndex(tcPDPtr newEmitter, tcPDPtr newSpectator,
			 const PDF& newEmitterPDF, const PDF& newSpectatorPDF)
  : theEmitterData(newEmitter), theInitialStateEmitter(newEmitterPDF.pdf()),
    theEmitterPDF(newEmitterPDF),
    theSpectatorData(newSpectator), theInitialStateSpectator(newSpectatorPDF.pdf()),
    theSpectatorPDF(newSpectatorPDF) {}


bool DipoleIndex::operator ==(const DipoleIndex& x) const {
  return
    theEmitterData == x.theEmitterData &&
    theInitialStateEmitter == x.theInitialStateEmitter &&
    theEmitterPDF == x.theEmitterPDF &&
    theSpectatorData == x.theSpectatorData &&
    theInitialStateSpectator == x.theInitialStateSpectator &&
    theSpectatorPDF == x.theSpectatorPDF;
}

bool DipoleIndex::operator <(const DipoleIndex& x) const {
  if ( theEmitterData == x.theEmitterData ) {
    if ( theInitialStateEmitter == x.theInitialStateEmitter ) {
      if ( theEmitterPDF == x.theEmitterPDF ) {
	if ( theSpectatorData == x.theSpectatorData ) {
	  if ( theInitialStateSpectator == x.theInitialStateSpectator ) {
	    return theSpectatorPDF < x.theSpectatorPDF;
	  }
	  return theInitialStateSpectator < x.theInitialStateSpectator;
	}
	return theSpectatorData < x.theSpectatorData;
      }
      return theEmitterPDF < x.theEmitterPDF;
    }
    return theInitialStateEmitter < x.theInitialStateEmitter;
  }
  return theEmitterData < x.theEmitterData;
}

void DipoleIndex::swap() {
  std::swap(theEmitterData,theSpectatorData);
  std::swap(theInitialStateEmitter,theInitialStateSpectator);
  std::swap(theEmitterPDF,theSpectatorPDF);
}

pair<DipoleIndex,DipoleIndex> DipoleIndex::split(tcPDPtr emm) const {

  DipoleIndex first(emitterData(),emm,emitterPDF(),PDF());
  DipoleIndex second(emm,spectatorData(),PDF(),spectatorPDF());

  return make_pair(first,second);

}

void DipoleIndex::print(ostream& os) const {
  os << "[" << emitterData()->PDGName();
  if ( emitterPDF().pdf() ) {
    os << "<-" << emitterPDF().particle()->PDGName()
       << "(" << emitterPDF().pdf() << ")";
  }
  os << "," << spectatorData()->PDGName();
  if ( spectatorPDF().pdf() ) {
    os << "<-" << spectatorPDF().particle()->PDGName()
       << "(" << spectatorPDF().pdf() << ")";
  }
  os << "]";
  os << flush;
}

DipoleSplittingInfo::DipoleSplittingInfo()
  : theConfiguration(false,false),
    theScale(0.0*GeV),
    theEmitterX(1.0), theSpectatorX(1.0), 
    theHardPt(0.0*GeV), theLastPt(0.0*GeV),
    theLastZ(0.0), theLastPhi(0.0), theLastEmitterZ(0.0),
    theLastSpectatorZ(0.0), theLastValue(0.0),
    theStoppedEvolving(false) {}

void DipoleSplittingInfo::print(ostream& os) const {

  os << "--- DipoleSplittingInfo --------------------------------------------------------\n";

  os << " index = " << theIndex << "\n";
  os << " configuration = (" << theConfiguration.first << "," << theConfiguration.second << ")\n"
     << " momentum fractions = [" << theEmitterX << "," << theSpectatorX << "]\n"
     << " generated starting from hard pt/GeV = " << theHardPt/GeV << "\n";

  if ( theEmitterData && theEmissionData && theSpectatorData ) {

    os << " splitting products = [(" << theEmitterData->PDGName()
       << "," << theEmissionData->PDGName() << "),"
       << theSpectatorData->PDGName() << "]\n";

  } else {

    os << " splitting products not available.\n";

  }

  if ( theSplittingKinematics ) {
    os << " kinematic variables associated to '" << theSplittingKinematics->name() << "':\n"
       << " scale = " << (theScale/GeV)
       << " pt/GeV = " << (theLastPt/GeV) << " z = " << theLastZ << " phi = " << theLastPhi << "\n"
       << " emitter z = " << theLastEmitterZ << " spectator z = " << theLastSpectatorZ << "\n"
       << " splitting kernel value = " << theLastValue << "\n"
       << " further parameters = ";
    copy(theLastSplittingParameters.begin(),theLastSplittingParameters.end(),ostream_iterator<double>(os," "));
    os << "\n the splitting " << (theStoppedEvolving ? "terminated " : "did not terminate ") << "the evolution\n";
  } else {
    os << " No kinematic variables have been generated yet.\n";
  }

  if ( theEmitter && theSpectator && theSplitEmitter && theSplitSpectator && theEmission ) {
    os << " the splitting has been performed:\n"
       << " emitter before emission:\n" << (*theEmitter)
       << " spectator before emission:\n" << (*theSpectator)
       << " emitter after emission:\n" << (*theSplitEmitter)
       << " emission:\n" << (*theEmission)
       << " spectator after emission:\n" << (*theSplitSpectator);
  } else {
    os << " the splitting has not yet been performed.\n";
  }

  os << "--------------------------------------------------------------------------------\n";

  os << flush;

}

