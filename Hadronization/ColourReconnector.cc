// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ColourReconnector class.
//

#include "ColourReconnector.h"
#include <ThePEG/Interface/ClassDocumentation.h>
#include <ThePEG/Interface/Parameter.h>
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/Repository/EventGenerator.h>

using namespace Herwig;

void ColourReconnector::persistentOutput(PersistentOStream & os) const {
  os << _ClReco << _PReco;
}

void ColourReconnector::persistentInput(PersistentIStream & is, int) {
  is >> _ClReco >> _PReco;
}

ClassDescription<ColourReconnector> ColourReconnector::initColourReconnector;
// Definition of the static class description member.


void ColourReconnector::Init() {

  static ClassDocumentation<ColourReconnector> documentation
    ("This class is responsible of the colour reconnection.");

  static Parameter<ColourReconnector,int>
    interfaceClReco ("ClReco","colour reconnection option",
                     &ColourReconnector::_ClReco, 0, 0, 0, 1,false,false,false);
  static Parameter<ColourReconnector,double>
    interfacePReco ("PReco","probability of colour reconnection",
                     &ColourReconnector::_PReco, 0, (1.0/9.0) , 0.0, 1.0,false,false,false);
  
}


void ColourReconnector::rearrange(EventHandler & ch, 
				  const StepPtr & pstep, 
				  ClusterVector & clusters) 
   throw(Veto, Stop, Exception){
  // Scan the particles in the Event record, and the "usual" clusters
  // stored in collecCluPtr.
  // If a new colour rearrangement is accepted, then  
  //       collecCluPtr.clear(); 
  // to get rid of the old clusters, and then add the new ones:
  //       collecCluPtr.insert( collecCluPtr.end(), clusterPointer );    
}






