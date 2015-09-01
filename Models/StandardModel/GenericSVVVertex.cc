// -*- C++ -*-
//
// GenericSVVVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GenericSVVVertex class.
//

#include "GenericSVVVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Command.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;  

GenericSVVVertex::GenericSVVVertex()
  :pids(ZERO),oas(0),oaew(0){
  orderInGs(0);
  orderInGem(0);
}

void GenericSVVVertex::doinit() {
  //PDG codes for particles at vertices
  addToList(pids[0],pids[1],pids[2]);
  orderInGs(oas);
  orderInGem(oaew);
  GeneralVVSVertex::doinit();
}



string GenericSVVVertex::dopids(string in) {
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





void GenericSVVVertex::persistentOutput(PersistentOStream & os) const {
  os << pids<<oas<<oaew;
}

void GenericSVVVertex::persistentInput(PersistentIStream & is, int) {
  is >> pids>>oas>>oaew;
}

ClassDescription<GenericSVVVertex> GenericSVVVertex::initGenericSVVVertex;
// Definition of the static class description member.

void GenericSVVVertex::Init() {
  
  static ClassDocumentation<GenericSVVVertex> documentation
    ("This class implements the s->v,v vertex");

  static Command<GenericSVVVertex> interfacepids
    ("pids",
     "Set the pids.",
     &GenericSVVVertex::dopids, false);

  static Parameter<GenericSVVVertex, int> interfaceOrderoas
    ("OrderInAlphaS",
     "The order in alpha_S",
     &GenericSVVVertex::oas, 2, 0, 0,
     false, false, Interface::lowerlim);
            
  static Parameter<GenericSVVVertex, int> interfaceOrderoaew
    ("OrderInAlphaEW",
     "The order in alpha_EW",
     &GenericSVVVertex::oaew, 2, 0, 0,
     false, false, Interface::lowerlim);
}

void GenericSVVVertex::setCoupling(Energy2, tcPDPtr part2, tcPDPtr part3, tcPDPtr part1) {
  assert(part1 && part2 && part3);
  assert(part1->id() == pids[0] &&
	 part2->id() == pids[1]  && part3->id() == pids[2] );
}

