// -*- C++ -*-
//
// GenericVVVVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GenericVVVVertex class.
//

#include "GenericVVVVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Command.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;  

GenericVVVVertex::GenericVVVVertex()
  :pids(ZERO),oas(0),oaew(0){
  orderInGs(0);
  orderInGem(0);
}

void GenericVVVVertex::doinit() {
  //PDG codes for particles at vertices
  addToList(pids[0],pids[1],pids[2]);
  orderInGs(oas);
  orderInGem(oaew);
  VVVVertex::doinit();
}



string GenericVVVVertex::dopids(string in) {
  vector<string> process = StringUtils::split(in);
  if ( process.size() != 3 )
    throw InitException() << "accepts only three particles.";

  for ( vector<string>::iterator p = process.begin();
	p != process.end(); ++p ) {
       int tmp;
       istringstream(*p) >> tmp;
       pids.push_back(tmp);
  }
  return "";
}





void GenericVVVVertex::persistentOutput(PersistentOStream & os) const {
  os << pids<<oas<<oaew;
}

void GenericVVVVertex::persistentInput(PersistentIStream & is, int) {
  is >> pids>>oas>>oaew;
}

ClassDescription<GenericVVVVertex> GenericVVVVertex::initGenericVVVVertex;
// Definition of the static class description member.

void GenericVVVVertex::Init() {
  
  static ClassDocumentation<GenericVVVVertex> documentation
    ("This class implements the v->v,v vertex");

  static Command<GenericVVVVertex> interfacepids
    ("pids",
     "Set the pids.",
     &GenericVVVVertex::dopids, false);

  static Parameter<GenericVVVVertex, int> interfaceOrderoas
    ("OrderInAlphaS",
     "The order in alpha_S",
     &GenericVVVVertex::oas, 2, 0, 0,
     false, false, Interface::lowerlim);
            
  static Parameter<GenericVVVVertex, int> interfaceOrderoaew
    ("OrderInAlphaEW",
     "The order in alpha_EW",
     &GenericVVVVertex::oaew, 2, 0, 0,
     false, false, Interface::lowerlim);
}

void GenericVVVVertex::setCoupling(Energy2, tcPDPtr part2, tcPDPtr part3, tcPDPtr part1) {
  assert(part1 && part2 && part3);
  assert(part1->id() == pids[0] &&
	 part2->id() == pids[1]  && part3->id() == pids[2] );
}

