// -*- C++ -*-
//
// NLOJetMEgg2ggg.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLOJetMEgg2ggg class.
//

#include "NLOJetMEgg2ggg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/NLOJetPhasespace.h"

using namespace Herwig;

NLOJetMEgg2ggg::NLOJetMEgg2ggg() {}

NLOJetMEgg2ggg::~NLOJetMEgg2ggg() {}

IBPtr NLOJetMEgg2ggg::clone() const {
  return new_ptr(*this);
}

IBPtr NLOJetMEgg2ggg::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void NLOJetMEgg2ggg::doinit() {
  NLOJetMEBase<0,2,0>::doinit();
}

void NLOJetMEgg2ggg::persistentOutput(PersistentOStream &) const {}

void NLOJetMEgg2ggg::persistentInput(PersistentIStream &, int) {}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<Herwig::NLOJetMEgg2ggg,Herwig::NLOJetMEBase<0,2,0> >
  describeNLOJetMEgg2ggg("Herwig::NLOJetMEgg2ggg", "HwMatchboxNLOJet.so");

void NLOJetMEgg2ggg::Init() {

  static ClassDocumentation<NLOJetMEgg2ggg> documentation
    ("NLOJetMEgg2ggg");

}


void NLOJetMEgg2ggg::doGetDiagrams() const {
  // get the particle data objects

  PDPtr g = getParticleData(ParticleID::g); 
  addSafe(new_ptr((Tree2toNDiagram(2), g, g, 1, g, 3, g, 4, g, 4, g, 3, g, -1)));
  addSafe(new_ptr((Tree2toNDiagram(2), g, g, 1, g, 3, g, 4, g, 3, g, 4, g, -2)));
  addSafe(new_ptr((Tree2toNDiagram(2), g, g, 1, g, 3, g, 3, g, 4, g, 4, g, -3)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 1, g, 3, g, 2, g, -4)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 1, g, 2, g, 3, g, -5)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 2, g, 1, g, 4, g, 4, g, -6)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 3, g, 1, g, 2, g, -7)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 2, g, 1, g, 3, g, -8)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 2, g, 4, g, 1, g, 4, g, -9)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 3, g, 2, g, 1, g, -10)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 2, g, 3, g, 1, g, -11)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 2, g, 4, g, 4, g, 1, g, -12)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 1, g, 2, g, 4, g, 4, g, -13)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 1, g, 4, g, 2, g, 4, g, -14)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 1, g, 4, g, 4, g, 2, g, -15)));
}


Selector <MEBase::DiagramIndex> NLOJetMEgg2ggg::diagrams(const DiagramVector & diags) const {
  // select the diagram
  matchboxAmplitude()->prepareAmplitudes(this);

  static const int cv1at0[]  = { -1, 3, 2, 1, 0, -998};
  static const int cv1at1[]  = { -1, 0, 3, 2, 1, -998};
  static const int cv1at2[]  = { -1, 3, 1, 2, 0, -998};
  static const int cv1at3[]  = { -1, 0, 3, 1, 2, -998};
  static const int cv1at4[]  = { -1, 2, 1, 3, 0, -998};
  static const int cv1at5[]  = { -1, 0, 2, 1, 3, -998};
  static const int cv1at6[]  = { -1, 1, 2, 3, 0, -998};
  static const int cv1at7[]  = { -1, 0, 1, 2, 3, -998};
  static const int cv2at0[]  = { -1, 2, 3, 1, 0, -998};
  static const int cv2at1[]  = { -1, 0, 2, 3, 1, -998};
  static const int cv2at2[]  = { -1, 2, 1, 3, 0, -998};
  static const int cv2at3[]  = { -1, 0, 2, 1, 3, -998};
  static const int cv2at4[]  = { -1, 3, 1, 2, 0, -998};
  static const int cv2at5[]  = { -1, 0, 3, 1, 2, -998};
  static const int cv2at6[]  = { -1, 1, 3, 2, 0, -998};
  static const int cv2at7[]  = { -1, 0, 1, 3, 2, -998};
  static const int cv3at0[]  = { -1, 1, 3, 2, 0, -998};
  static const int cv3at1[]  = { -1, 0, 1, 3, 2, -998};
  static const int cv3at2[]  = { -1, 1, 2, 3, 0, -998};
  static const int cv3at3[]  = { -1, 0, 1, 2, 3, -998};
  static const int cv3at4[]  = { -1, 3, 2, 1, 0, -998};
  static const int cv3at5[]  = { -1, 0, 3, 2, 1, -998};
  static const int cv3at6[]  = { -1, 2, 3, 1, 0, -998};
  static const int cv3at7[]  = { -1, 0, 2, 3, 1, -998};
  static const int cv4at0[]  = { -1, 3, 2, 0, 1, -998};
  static const int cv4at1[]  = { -1, 1, 3, 2, 0, -998};
  static const int cv4at2[]  = { -1, 3, 0, 2, 1, -998};
  static const int cv4at3[]  = { -1, 1, 3, 0, 2, -998};
  static const int cv4at4[]  = { -1, 2, 0, 3, 1, -998};
  static const int cv4at5[]  = { -1, 1, 2, 0, 3, -998};
  static const int cv4at6[]  = { -1, 0, 2, 3, 1, -998};
  static const int cv4at7[]  = { -1, 1, 0, 2, 3, -998};
  static const int cv5at0[]  = { -1, 2, 3, 0, 1, -998};
  static const int cv5at1[]  = { -1, 1, 2, 3, 0, -998};
  static const int cv5at2[]  = { -1, 2, 0, 3, 1, -998};
  static const int cv5at3[]  = { -1, 1, 2, 0, 3, -998};
  static const int cv5at4[]  = { -1, 3, 0, 2, 1, -998};
  static const int cv5at5[]  = { -1, 1, 3, 0, 2, -998};
  static const int cv5at6[]  = { -1, 0, 3, 2, 1, -998};
  static const int cv5at7[]  = { -1, 1, 0, 3, 2, -998};
  static const int cv6at0[]  = { -1, 0, 3, 2, 1, -998};
  static const int cv6at1[]  = { -1, 1, 0, 3, 2, -998};
  static const int cv6at2[]  = { -1, 0, 2, 3, 1, -998};
  static const int cv6at3[]  = { -1, 1, 0, 2, 3, -998};
  static const int cv6at4[]  = { -1, 3, 2, 0, 1, -998};
  static const int cv6at5[]  = { -1, 1, 3, 2, 0, -998};
  static const int cv6at6[]  = { -1, 2, 3, 0, 1, -998};
  static const int cv6at7[]  = { -1, 1, 2, 3, 0, -998};
  static const int cv7at0[]  = { -1, 3, 1, 0, 2, -998};
  static const int cv7at1[]  = { -1, 2, 3, 1, 0, -998};
  static const int cv7at2[]  = { -1, 3, 0, 1, 2, -998};
  static const int cv7at3[]  = { -1, 2, 3, 0, 1, -998};
  static const int cv7at4[]  = { -1, 1, 0, 3, 2, -998};
  static const int cv7at5[]  = { -1, 2, 1, 0, 3, -998};
  static const int cv7at6[]  = { -1, 0, 1, 3, 2, -998};
  static const int cv7at7[]  = { -1, 2, 0, 1, 3, -998};
  static const int cv8at0[]  = { -1, 1, 3, 0, 2, -998};
  static const int cv8at1[]  = { -1, 2, 1, 3, 0, -998};
  static const int cv8at2[]  = { -1, 1, 0, 3, 2, -998};
  static const int cv8at3[]  = { -1, 2, 1, 0, 3, -998};
  static const int cv8at4[]  = { -1, 3, 0, 1, 2, -998};
  static const int cv8at5[]  = { -1, 2, 3, 0, 1, -998};
  static const int cv8at6[]  = { -1, 0, 3, 1, 2, -998};
  static const int cv8at7[]  = { -1, 2, 0, 3, 1, -998};
  static const int cv9at0[]  = { -1, 0, 3, 1, 2, -998};
  static const int cv9at1[]  = { -1, 2, 0, 3, 1, -998};
  static const int cv9at2[]  = { -1, 0, 1, 3, 2, -998};
  static const int cv9at3[]  = { -1, 2, 0, 1, 3, -998};
  static const int cv9at4[]  = { -1, 3, 1, 0, 2, -998};
  static const int cv9at5[]  = { -1, 2, 3, 1, 0, -998};
  static const int cv9at6[]  = { -1, 1, 3, 0, 2, -998};
  static const int cv9at7[]  = { -1, 2, 1, 3, 0, -998};
  static const int cv10at0[]  = { -1, 2, 1, 0, 3, -998};
  static const int cv10at1[]  = { -1, 3, 2, 1, 0, -998};
  static const int cv10at2[]  = { -1, 2, 0, 1, 3, -998};
  static const int cv10at3[]  = { -1, 3, 2, 0, 1, -998};
  static const int cv10at4[]  = { -1, 1, 0, 2, 3, -998};
  static const int cv10at5[]  = { -1, 3, 1, 0, 2, -998};
  static const int cv10at6[]  = { -1, 0, 1, 2, 3, -998};
  static const int cv10at7[]  = { -1, 3, 0, 1, 2, -998};
  static const int cv11at0[]  = { -1, 1, 2, 0, 3, -998};
  static const int cv11at1[]  = { -1, 3, 1, 2, 0, -998};
  static const int cv11at2[]  = { -1, 1, 0, 2, 3, -998};
  static const int cv11at3[]  = { -1, 3, 1, 0, 2, -998};
  static const int cv11at4[]  = { -1, 2, 0, 1, 3, -998};
  static const int cv11at5[]  = { -1, 3, 2, 0, 1, -998};
  static const int cv11at6[]  = { -1, 0, 2, 1, 3, -998};
  static const int cv11at7[]  = { -1, 3, 0, 2, 1, -998};
  static const int cv12at0[]  = { -1, 0, 2, 1, 3, -998};
  static const int cv12at1[]  = { -1, 3, 0, 2, 1, -998};
  static const int cv12at2[]  = { -1, 0, 1, 2, 3, -998};
  static const int cv12at3[]  = { -1, 3, 0, 1, 2, -998};
  static const int cv12at4[]  = { -1, 2, 1, 0, 3, -998};
  static const int cv12at5[]  = { -1, 3, 2, 1, 0, -998};
  static const int cv12at6[]  = { -1, 1, 2, 0, 3, -998};
  static const int cv12at7[]  = { -1, 3, 1, 2, 0, -998};
  static const int cv13at0[]  = { -1, 3, 2, 1, 0, -998};
  static const int cv13at1[]  = { -1, 3, 2, 0, 1, -998};
  static const int cv13at2[]  = { -1, 2, 3, 1, 0, -998};
  static const int cv13at3[]  = { -1, 2, 3, 0, 1, -998};
  static const int cv13at4[]  = { -1, 1, 0, 3, 2, -998};
  static const int cv13at5[]  = { -1, 0, 1, 3, 2, -998};
  static const int cv13at6[]  = { -1, 1, 0, 2, 3, -998};
  static const int cv13at7[]  = { -1, 0, 1, 2, 3, -998};
  static const int cv14at0[]  = { -1, 3, 1, 2, 0, -998};
  static const int cv14at1[]  = { -1, 3, 1, 0, 2, -998};
  static const int cv14at2[]  = { -1, 1, 3, 2, 0, -998};
  static const int cv14at3[]  = { -1, 1, 3, 0, 2, -998};
  static const int cv14at4[]  = { -1, 2, 0, 3, 1, -998};
  static const int cv14at5[]  = { -1, 0, 2, 3, 1, -998};
  static const int cv14at6[]  = { -1, 2, 0, 1, 3, -998};
  static const int cv14at7[]  = { -1, 0, 2, 1, 3, -998};
  static const int cv15at0[]  = { -1, 2, 1, 3, 0, -998};
  static const int cv15at1[]  = { -1, 2, 1, 0, 3, -998};
  static const int cv15at2[]  = { -1, 1, 2, 3, 0, -998};
  static const int cv15at3[]  = { -1, 1, 2, 0, 3, -998};
  static const int cv15at4[]  = { -1, 3, 0, 2, 1, -998};
  static const int cv15at5[]  = { -1, 0, 3, 2, 1, -998};
  static const int cv15at6[]  = { -1, 3, 0, 1, 2, -998};
  static const int cv15at7[]  = { -1, 0, 3, 1, 2, -998};

  Selector <MEBase::DiagramIndex> sel;
  for(MEBase::DiagramIndex i=0; i < diags.size(); ++i){
    if( diags[i]->id() == -1 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv1at0, sizeof(cv1at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv1at1, sizeof(cv1at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv1at2, sizeof(cv1at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv1at3, sizeof(cv1at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv1at4, sizeof(cv1at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv1at5, sizeof(cv1at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv1at6, sizeof(cv1at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv1at7, sizeof(cv1at7)/sizeof(int)), i);   
    else if( diags[i]->id() == -2 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv2at0, sizeof(cv2at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv2at1, sizeof(cv2at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv2at2, sizeof(cv2at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv2at3, sizeof(cv2at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv2at4, sizeof(cv2at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv2at5, sizeof(cv2at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv2at6, sizeof(cv2at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv2at7, sizeof(cv2at7)/sizeof(int)), i);   
    else if( diags[i]->id() == -3 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv3at0, sizeof(cv3at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv3at1, sizeof(cv3at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv3at2, sizeof(cv3at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv3at3, sizeof(cv3at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv3at4, sizeof(cv3at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv3at5, sizeof(cv3at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv3at6, sizeof(cv3at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv3at7, sizeof(cv3at7)/sizeof(int)), i);   
    else if( diags[i]->id() == -4 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv4at0, sizeof(cv4at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv4at1, sizeof(cv4at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv4at2, sizeof(cv4at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv4at3, sizeof(cv4at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv4at4, sizeof(cv4at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv4at5, sizeof(cv4at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv4at6, sizeof(cv4at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv4at7, sizeof(cv4at7)/sizeof(int)), i);   
    else if( diags[i]->id() == -5 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv5at0, sizeof(cv5at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv5at1, sizeof(cv5at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv5at2, sizeof(cv5at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv5at3, sizeof(cv5at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv5at4, sizeof(cv5at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv5at5, sizeof(cv5at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv5at6, sizeof(cv5at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv5at7, sizeof(cv5at7)/sizeof(int)), i);   
    else if( diags[i]->id() == -6 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv6at0, sizeof(cv6at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv6at1, sizeof(cv6at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv6at2, sizeof(cv6at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv6at3, sizeof(cv6at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv6at4, sizeof(cv6at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv6at5, sizeof(cv6at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv6at6, sizeof(cv6at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv6at7, sizeof(cv6at7)/sizeof(int)), i);   
    else if( diags[i]->id() == -7 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv7at0, sizeof(cv7at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv7at1, sizeof(cv7at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv7at2, sizeof(cv7at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv7at3, sizeof(cv7at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv7at4, sizeof(cv7at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv7at5, sizeof(cv7at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv7at6, sizeof(cv7at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv7at7, sizeof(cv7at7)/sizeof(int)), i);   
    else if( diags[i]->id() == -8 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv8at0, sizeof(cv8at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv8at1, sizeof(cv8at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv8at2, sizeof(cv8at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv8at3, sizeof(cv8at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv8at4, sizeof(cv8at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv8at5, sizeof(cv8at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv8at6, sizeof(cv8at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv8at7, sizeof(cv8at7)/sizeof(int)), i);   
    else if( diags[i]->id() == -9 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv9at0, sizeof(cv9at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv9at1, sizeof(cv9at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv9at2, sizeof(cv9at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv9at3, sizeof(cv9at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv9at4, sizeof(cv9at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv9at5, sizeof(cv9at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv9at6, sizeof(cv9at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv9at7, sizeof(cv9at7)/sizeof(int)), i);   
    else if( diags[i]->id() == -10 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv10at0, sizeof(cv10at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv10at1, sizeof(cv10at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv10at2, sizeof(cv10at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv10at3, sizeof(cv10at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv10at4, sizeof(cv10at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv10at5, sizeof(cv10at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv10at6, sizeof(cv10at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv10at7, sizeof(cv10at7)/sizeof(int)), i);   
    else if( diags[i]->id() == -11 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv11at0, sizeof(cv11at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv11at1, sizeof(cv11at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv11at2, sizeof(cv11at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv11at3, sizeof(cv11at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv11at4, sizeof(cv11at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv11at5, sizeof(cv11at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv11at6, sizeof(cv11at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv11at7, sizeof(cv11at7)/sizeof(int)), i);   
    else if( diags[i]->id() == -12 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv12at0, sizeof(cv12at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv12at1, sizeof(cv12at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv12at2, sizeof(cv12at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv12at3, sizeof(cv12at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv12at4, sizeof(cv12at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv12at5, sizeof(cv12at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv12at6, sizeof(cv12at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv12at7, sizeof(cv12at7)/sizeof(int)), i);   
    else if( diags[i]->id() == -13 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv13at0, sizeof(cv13at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv13at1, sizeof(cv13at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv13at2, sizeof(cv13at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv13at3, sizeof(cv13at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv13at4, sizeof(cv13at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv13at5, sizeof(cv13at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv13at6, sizeof(cv13at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv13at7, sizeof(cv13at7)/sizeof(int)), i);   
    else if( diags[i]->id() == -14 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv14at0, sizeof(cv14at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv14at1, sizeof(cv14at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv14at2, sizeof(cv14at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv14at3, sizeof(cv14at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv14at4, sizeof(cv14at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv14at5, sizeof(cv14at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv14at6, sizeof(cv14at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv14at7, sizeof(cv14at7)/sizeof(int)), i);   
    else if( diags[i]->id() == -15 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv15at0, sizeof(cv15at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv15at1, sizeof(cv15at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv15at2, sizeof(cv15at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv15at3, sizeof(cv15at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv15at4, sizeof(cv15at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv15at5, sizeof(cv15at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv15at6, sizeof(cv15at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv15at7, sizeof(cv15at7)/sizeof(int)), i);   
  }
  return sel;
}


Selector<const ColourLines *> NLOJetMEgg2ggg::colourGeometries(tcDiagPtr diag) const {
  // colour lines for NLOJetMEgg2ggg

  static const ColourLines diag1[8] = { 
    ColourLines("1 3 7, -1 2, -2 -3 -4 -5, 5 -6, 6 4 -7"), 
    ColourLines("1 -2, -1 -3 -4 -5, 2 3 7, 5 -6, 6 4 -7"), 
    ColourLines("1 3 7, -1 2, -2 -3 -4 -6, 5 4 -7, -5 6"), 
    ColourLines("1 -2, -1 -3 -4 -6, 2 3 7, 5 4 -7, -5 6"), 
    ColourLines("1 3 4 6, -1 2, -2 -3 -7, 5 -6, -5 -4 7"), 
    ColourLines("1 -2, -1 -3 -7, 2 3 4 6, 5 -6, -5 -4 7"), 
    ColourLines("1 3 4 5, -1 2, -2 -3 -7, -5 6, -6 -4 7"), 
    ColourLines("1 -2, -1 -3 -7, 2 3 4 5, -5 6, -6 -4 7")
  }; 
  static const ColourLines diag2[8] = { 
    ColourLines("1 3 6, -1 2, -2 -3 -4 -5, 5 -7, -6 4 7"), 
    ColourLines("1 -2, -1 -3 -4 -5, 2 3 6, 5 -7, -6 4 7"), 
    ColourLines("1 3 6, -1 2, -2 -3 -4 -7, 5 4 -6, -5 7"), 
    ColourLines("1 -2, -1 -3 -4 -7, 2 3 6, 5 4 -6, -5 7"), 
    ColourLines("1 3 4 7, -1 2, -2 -3 -6, 5 -7, -5 -4 6"), 
    ColourLines("1 -2, -1 -3 -6, 2 3 4 7, 5 -7, -5 -4 6"), 
    ColourLines("1 3 4 5, -1 2, -2 -3 -6, -5 7, 6 -4 -7"), 
    ColourLines("1 -2, -1 -3 -6, 2 3 4 5, -5 7, 6 -4 -7")
  }; 
  static const ColourLines diag3[8] = { 
    ColourLines("1 3 5, -1 2, -2 -3 -4 -6, -5 4 7, 6 -7"), 
    ColourLines("1 -2, -1 -3 -4 -6, 2 3 5, -5 4 7, 6 -7"), 
    ColourLines("1 3 5, -1 2, -2 -3 -4 -7, -5 4 6, -6 7"), 
    ColourLines("1 -2, -1 -3 -4 -7, 2 3 5, -5 4 6, -6 7"), 
    ColourLines("1 3 4 7, -1 2, -2 -3 -5, 5 -4 -6, 6 -7"), 
    ColourLines("1 -2, -1 -3 -5, 2 3 4 7, 5 -4 -6, 6 -7"), 
    ColourLines("1 3 4 6, -1 2, -2 -3 -5, 5 -4 -7, -6 7"), 
    ColourLines("1 -2, -1 -3 -5, 2 3 4 6, 5 -4 -7, -6 7")
  }; 
  static const ColourLines diag4[8] = { 
    ColourLines("1 2 7, -1 -5, 4 -3 -2 5, -4 -6, 6 3 -7"), 
    ColourLines("1 5, -1 -2 -3 4, -4 -6, -5 2 7, 6 3 -7"), 
    ColourLines("1 2 7, -1 -5, 4 6, -4 3 -7, 5 -2 -3 -6"), 
    ColourLines("1 5, -1 -2 -3 -6, 4 6, -4 3 -7, -5 2 7"), 
    ColourLines("1 2 3 6, -1 -5, 4 -3 7, -4 -6, 5 -2 -7"), 
    ColourLines("1 5, -1 -2 -7, 4 -3 7, -4 -6, -5 2 3 6"), 
    ColourLines("1 2 3 -4, -1 -5, 4 6, 5 -2 -7, -6 -3 7"), 
    ColourLines("1 5, -1 -2 -7, 4 6, -4 3 2 -5, -6 -3 7")
  }; 
  static const ColourLines diag5[8] = { 
    ColourLines("1 2 6, -1 -5, 4 -3 -2 5, -4 -7, -6 3 7"), 
    ColourLines("1 5, -1 -2 -3 4, -4 -7, -5 2 6, -6 3 7"), 
    ColourLines("1 2 6, -1 -5, 4 7, -4 3 -6, 5 -2 -3 -7"), 
    ColourLines("1 5, -1 -2 -3 -7, 4 7, -4 3 -6, -5 2 6"), 
    ColourLines("1 2 3 7, -1 -5, 4 -3 6, -4 -7, 5 -2 -6"), 
    ColourLines("1 5, -1 -2 -6, 4 -3 6, -4 -7, -5 2 3 7"), 
    ColourLines("1 2 3 -4, -1 -5, 4 7, 5 -2 -6, 6 -3 -7"), 
    ColourLines("1 5, -1 -2 -6, 4 7, -4 3 2 -5, 6 -3 -7")
  }; 
  static const ColourLines diag6[8] = { 
    ColourLines("1 2 -3, -1 -5, 3 4 7, 5 -2 -4 -6, 6 -7"), 
    ColourLines("1 5, -1 -2 -4 -6, 3 4 7, -3 2 -5, 6 -7"), 
    ColourLines("1 2 -3, -1 -5, 3 4 6, 5 -2 -4 -7, -6 7"), 
    ColourLines("1 5, -1 -2 -4 -7, 3 4 6, -3 2 -5, -6 7"), 
    ColourLines("1 2 4 7, -1 -5, 3 -2 5, -3 -4 -6, 6 -7"), 
    ColourLines("1 5, -1 -2 3, -3 -4 -6, -5 2 4 7, 6 -7"), 
    ColourLines("1 2 4 6, -1 -5, 3 -2 5, -3 -4 -7, -6 7"), 
    ColourLines("1 5, -1 -2 3, -3 -4 -7, -5 2 4 6, -6 7")
  }; 
  static const ColourLines diag7[8] = { 
    ColourLines("1 2 7, -1 -6, 4 -3 -2 6, -4 -5, 5 3 -7"), 
    ColourLines("1 6, -1 -2 -3 4, -4 -5, 5 3 -7, -6 2 7"), 
    ColourLines("1 2 7, -1 -6, 4 5, -4 3 -7, -5 -3 -2 6"), 
    ColourLines("1 6, -1 -2 -3 -5, 4 5, -4 3 -7, -6 2 7"), 
    ColourLines("1 2 3 5, -1 -6, 4 -3 7, -4 -5, 6 -2 -7"), 
    ColourLines("1 6, -1 -2 -7, 4 -3 7, -4 -5, 5 3 2 -6"), 
    ColourLines("1 2 3 -4, -1 -6, 4 5, -5 -3 7, 6 -2 -7"), 
    ColourLines("1 6, -1 -2 -7, 4 5, -4 3 2 -6, -5 -3 7")
  }; 
  static const ColourLines diag8[8] = { 
    ColourLines("1 2 5, -1 -6, 4 -3 -2 6, -4 -7, -5 3 7"), 
    ColourLines("1 6, -1 -2 -3 4, -4 -7, 5 2 -6, -5 3 7"), 
    ColourLines("1 2 5, -1 -6, 4 7, -4 3 -5, 6 -2 -3 -7"), 
    ColourLines("1 6, -1 -2 -3 -7, 4 7, -4 3 -5, 5 2 -6"), 
    ColourLines("1 2 3 7, -1 -6, 4 -3 5, -4 -7, -5 -2 6"), 
    ColourLines("1 6, -1 -2 -5, 4 -3 5, -4 -7, -6 2 3 7"), 
    ColourLines("1 2 3 -4, -1 -6, 4 7, 5 -3 -7, -5 -2 6"), 
    ColourLines("1 6, -1 -2 -5, 4 7, -4 3 2 -6, 5 -3 -7")
  }; 
  static const ColourLines diag9[8] = { 
    ColourLines("1 2 -3, -1 -6, 3 4 7, 5 -7, -5 -4 -2 6"), 
    ColourLines("1 6, -1 -2 -4 -5, 3 4 7, -3 2 -6, 5 -7"), 
    ColourLines("1 2 -3, -1 -6, 3 4 5, -5 7, 6 -2 -4 -7"), 
    ColourLines("1 6, -1 -2 -4 -7, 3 4 5, -3 2 -6, -5 7"), 
    ColourLines("1 2 4 7, -1 -6, 3 -2 6, -3 -4 -5, 5 -7"), 
    ColourLines("1 6, -1 -2 3, -3 -4 -5, 5 -7, -6 2 4 7"), 
    ColourLines("1 2 4 5, -1 -6, 3 -2 6, -3 -4 -7, -5 7"), 
    ColourLines("1 6, -1 -2 3, -3 -4 -7, 5 4 2 -6, -5 7")
  }; 
  static const ColourLines diag10[8] = { 
    ColourLines("1 2 6, -1 -7, 4 -3 -2 7, -4 -5, 5 3 -6"), 
    ColourLines("1 7, -1 -2 -3 4, -4 -5, 5 3 -6, 6 2 -7"), 
    ColourLines("1 2 6, -1 -7, 4 5, -4 3 -6, -5 -3 -2 7"), 
    ColourLines("1 7, -1 -2 -3 -5, 4 5, -4 3 -6, 6 2 -7"), 
    ColourLines("1 2 3 5, -1 -7, 4 -3 6, -4 -5, -6 -2 7"), 
    ColourLines("1 7, -1 -2 -6, 4 -3 6, -4 -5, 5 3 2 -7"), 
    ColourLines("1 2 3 -4, -1 -7, 4 5, -5 -3 6, -6 -2 7"), 
    ColourLines("1 7, -1 -2 -6, 4 5, -4 3 2 -7, -5 -3 6")
  }; 
  static const ColourLines diag11[8] = { 
    ColourLines("1 2 5, -1 -7, 4 -3 -2 7, -4 -6, -5 3 6"), 
    ColourLines("1 7, -1 -2 -3 4, -4 -6, 5 2 -7, -5 3 6"), 
    ColourLines("1 2 5, -1 -7, 4 6, -4 3 -5, -6 -3 -2 7"), 
    ColourLines("1 7, -1 -2 -3 -6, 4 6, -4 3 -5, 5 2 -7"), 
    ColourLines("1 2 3 6, -1 -7, 4 -3 5, -4 -6, -5 -2 7"), 
    ColourLines("1 7, -1 -2 -5, 4 -3 5, -4 -6, 6 3 2 -7"), 
    ColourLines("1 2 3 -4, -1 -7, 4 6, 5 -3 -6, -5 -2 7"), 
    ColourLines("1 7, -1 -2 -5, 4 6, -4 3 2 -7, 5 -3 -6")
  }; 
  static const ColourLines diag12[8] = { 
    ColourLines("1 2 -3, -1 -7, 3 4 6, 5 -6, -5 -4 -2 7"), 
    ColourLines("1 7, -1 -2 -4 -5, 3 4 6, -3 2 -7, 5 -6"), 
    ColourLines("1 2 -3, -1 -7, 3 4 5, -5 6, -6 -4 -2 7"), 
    ColourLines("1 7, -1 -2 -4 -6, 3 4 5, -3 2 -7, -5 6"), 
    ColourLines("1 2 4 6, -1 -7, 3 -2 7, -3 -4 -5, 5 -6"), 
    ColourLines("1 7, -1 -2 3, -3 -4 -5, 5 -6, 6 4 2 -7"), 
    ColourLines("1 2 4 5, -1 -7, 3 -2 7, -3 -4 -6, -5 6"), 
    ColourLines("1 7, -1 -2 3, -3 -4 -6, 5 4 2 -7, -5 6")
  }; 
  static const ColourLines diag13[8] = { 
    ColourLines("1 4 7, -1 -2 3, -3 -5, 5 2 -4 -6, 6 -7"), 
    ColourLines("1 4 7, -1 -2 -5, 3 5, -3 2 -4 -6, 6 -7"), 
    ColourLines("1 4 6, -1 -2 3, -3 -5, 5 2 -4 -7, -6 7"), 
    ColourLines("1 4 6, -1 -2 -5, 3 5, -3 2 -4 -7, -6 7"), 
    ColourLines("1 2 5, -1 -4 -6, 3 -2 4 7, -3 -5, 6 -7"), 
    ColourLines("1 2 -3, -1 -4 -6, 3 5, -5 -2 4 7, 6 -7"), 
    ColourLines("1 2 5, -1 -4 -7, 3 -2 4 6, -3 -5, -6 7"), 
    ColourLines("1 2 -3, -1 -4 -7, 3 5, -5 -2 4 6, -6 7")
  }; 
  static const ColourLines diag14[8] = { 
    ColourLines("1 4 7, -1 -2 3, -3 -6, 5 -7, -5 -4 2 6"), 
    ColourLines("1 4 7, -1 -2 -6, 3 6, -3 2 -4 -5, 5 -7"), 
    ColourLines("1 4 5, -1 -2 3, -3 -6, -5 7, 6 2 -4 -7"), 
    ColourLines("1 4 5, -1 -2 -6, 3 6, -3 2 -4 -7, -5 7"), 
    ColourLines("1 2 6, -1 -4 -5, 3 -2 4 7, -3 -6, 5 -7"), 
    ColourLines("1 2 -3, -1 -4 -5, 3 6, 5 -7, -6 -2 4 7"), 
    ColourLines("1 2 6, -1 -4 -7, 3 -2 4 5, -3 -6, -5 7"), 
    ColourLines("1 2 -3, -1 -4 -7, 3 6, 5 4 -2 -6, -5 7")
  }; 
  static const ColourLines diag15[8] = { 
    ColourLines("1 4 6, -1 -2 3, -3 -7, 5 -6, -5 -4 2 7"), 
    ColourLines("1 4 6, -1 -2 -7, 3 7, -3 2 -4 -5, 5 -6"), 
    ColourLines("1 4 5, -1 -2 3, -3 -7, -5 6, -6 -4 2 7"), 
    ColourLines("1 4 5, -1 -2 -7, 3 7, -3 2 -4 -6, -5 6"), 
    ColourLines("1 2 7, -1 -4 -5, 3 -2 4 6, -3 -7, 5 -6"), 
    ColourLines("1 2 -3, -1 -4 -5, 3 7, 5 -6, 6 4 -2 -7"), 
    ColourLines("1 2 7, -1 -4 -6, 3 -2 4 5, -3 -7, -5 6"), 
    ColourLines("1 2 -3, -1 -4 -6, 3 7, 5 4 -2 -7, -5 6")
  }; 

  static const int cv1at0[]  = { -1, 3, 2, 1, 0, -998};
  static const int cv1at1[]  = { -1, 0, 3, 2, 1, -998};
  static const int cv1at2[]  = { -1, 3, 1, 2, 0, -998};
  static const int cv1at3[]  = { -1, 0, 3, 1, 2, -998};
  static const int cv1at4[]  = { -1, 2, 1, 3, 0, -998};
  static const int cv1at5[]  = { -1, 0, 2, 1, 3, -998};
  static const int cv1at6[]  = { -1, 1, 2, 3, 0, -998};
  static const int cv1at7[]  = { -1, 0, 1, 2, 3, -998};
  static const int cv2at0[]  = { -1, 2, 3, 1, 0, -998};
  static const int cv2at1[]  = { -1, 0, 2, 3, 1, -998};
  static const int cv2at2[]  = { -1, 2, 1, 3, 0, -998};
  static const int cv2at3[]  = { -1, 0, 2, 1, 3, -998};
  static const int cv2at4[]  = { -1, 3, 1, 2, 0, -998};
  static const int cv2at5[]  = { -1, 0, 3, 1, 2, -998};
  static const int cv2at6[]  = { -1, 1, 3, 2, 0, -998};
  static const int cv2at7[]  = { -1, 0, 1, 3, 2, -998};
  static const int cv3at0[]  = { -1, 1, 3, 2, 0, -998};
  static const int cv3at1[]  = { -1, 0, 1, 3, 2, -998};
  static const int cv3at2[]  = { -1, 1, 2, 3, 0, -998};
  static const int cv3at3[]  = { -1, 0, 1, 2, 3, -998};
  static const int cv3at4[]  = { -1, 3, 2, 1, 0, -998};
  static const int cv3at5[]  = { -1, 0, 3, 2, 1, -998};
  static const int cv3at6[]  = { -1, 2, 3, 1, 0, -998};
  static const int cv3at7[]  = { -1, 0, 2, 3, 1, -998};
  static const int cv4at0[]  = { -1, 3, 2, 0, 1, -998};
  static const int cv4at1[]  = { -1, 1, 3, 2, 0, -998};
  static const int cv4at2[]  = { -1, 3, 0, 2, 1, -998};
  static const int cv4at3[]  = { -1, 1, 3, 0, 2, -998};
  static const int cv4at4[]  = { -1, 2, 0, 3, 1, -998};
  static const int cv4at5[]  = { -1, 1, 2, 0, 3, -998};
  static const int cv4at6[]  = { -1, 0, 2, 3, 1, -998};
  static const int cv4at7[]  = { -1, 1, 0, 2, 3, -998};
  static const int cv5at0[]  = { -1, 2, 3, 0, 1, -998};
  static const int cv5at1[]  = { -1, 1, 2, 3, 0, -998};
  static const int cv5at2[]  = { -1, 2, 0, 3, 1, -998};
  static const int cv5at3[]  = { -1, 1, 2, 0, 3, -998};
  static const int cv5at4[]  = { -1, 3, 0, 2, 1, -998};
  static const int cv5at5[]  = { -1, 1, 3, 0, 2, -998};
  static const int cv5at6[]  = { -1, 0, 3, 2, 1, -998};
  static const int cv5at7[]  = { -1, 1, 0, 3, 2, -998};
  static const int cv6at0[]  = { -1, 0, 3, 2, 1, -998};
  static const int cv6at1[]  = { -1, 1, 0, 3, 2, -998};
  static const int cv6at2[]  = { -1, 0, 2, 3, 1, -998};
  static const int cv6at3[]  = { -1, 1, 0, 2, 3, -998};
  static const int cv6at4[]  = { -1, 3, 2, 0, 1, -998};
  static const int cv6at5[]  = { -1, 1, 3, 2, 0, -998};
  static const int cv6at6[]  = { -1, 2, 3, 0, 1, -998};
  static const int cv6at7[]  = { -1, 1, 2, 3, 0, -998};
  static const int cv7at0[]  = { -1, 3, 1, 0, 2, -998};
  static const int cv7at1[]  = { -1, 2, 3, 1, 0, -998};
  static const int cv7at2[]  = { -1, 3, 0, 1, 2, -998};
  static const int cv7at3[]  = { -1, 2, 3, 0, 1, -998};
  static const int cv7at4[]  = { -1, 1, 0, 3, 2, -998};
  static const int cv7at5[]  = { -1, 2, 1, 0, 3, -998};
  static const int cv7at6[]  = { -1, 0, 1, 3, 2, -998};
  static const int cv7at7[]  = { -1, 2, 0, 1, 3, -998};
  static const int cv8at0[]  = { -1, 1, 3, 0, 2, -998};
  static const int cv8at1[]  = { -1, 2, 1, 3, 0, -998};
  static const int cv8at2[]  = { -1, 1, 0, 3, 2, -998};
  static const int cv8at3[]  = { -1, 2, 1, 0, 3, -998};
  static const int cv8at4[]  = { -1, 3, 0, 1, 2, -998};
  static const int cv8at5[]  = { -1, 2, 3, 0, 1, -998};
  static const int cv8at6[]  = { -1, 0, 3, 1, 2, -998};
  static const int cv8at7[]  = { -1, 2, 0, 3, 1, -998};
  static const int cv9at0[]  = { -1, 0, 3, 1, 2, -998};
  static const int cv9at1[]  = { -1, 2, 0, 3, 1, -998};
  static const int cv9at2[]  = { -1, 0, 1, 3, 2, -998};
  static const int cv9at3[]  = { -1, 2, 0, 1, 3, -998};
  static const int cv9at4[]  = { -1, 3, 1, 0, 2, -998};
  static const int cv9at5[]  = { -1, 2, 3, 1, 0, -998};
  static const int cv9at6[]  = { -1, 1, 3, 0, 2, -998};
  static const int cv9at7[]  = { -1, 2, 1, 3, 0, -998};
  static const int cv10at0[]  = { -1, 2, 1, 0, 3, -998};
  static const int cv10at1[]  = { -1, 3, 2, 1, 0, -998};
  static const int cv10at2[]  = { -1, 2, 0, 1, 3, -998};
  static const int cv10at3[]  = { -1, 3, 2, 0, 1, -998};
  static const int cv10at4[]  = { -1, 1, 0, 2, 3, -998};
  static const int cv10at5[]  = { -1, 3, 1, 0, 2, -998};
  static const int cv10at6[]  = { -1, 0, 1, 2, 3, -998};
  static const int cv10at7[]  = { -1, 3, 0, 1, 2, -998};
  static const int cv11at0[]  = { -1, 1, 2, 0, 3, -998};
  static const int cv11at1[]  = { -1, 3, 1, 2, 0, -998};
  static const int cv11at2[]  = { -1, 1, 0, 2, 3, -998};
  static const int cv11at3[]  = { -1, 3, 1, 0, 2, -998};
  static const int cv11at4[]  = { -1, 2, 0, 1, 3, -998};
  static const int cv11at5[]  = { -1, 3, 2, 0, 1, -998};
  static const int cv11at6[]  = { -1, 0, 2, 1, 3, -998};
  static const int cv11at7[]  = { -1, 3, 0, 2, 1, -998};
  static const int cv12at0[]  = { -1, 0, 2, 1, 3, -998};
  static const int cv12at1[]  = { -1, 3, 0, 2, 1, -998};
  static const int cv12at2[]  = { -1, 0, 1, 2, 3, -998};
  static const int cv12at3[]  = { -1, 3, 0, 1, 2, -998};
  static const int cv12at4[]  = { -1, 2, 1, 0, 3, -998};
  static const int cv12at5[]  = { -1, 3, 2, 1, 0, -998};
  static const int cv12at6[]  = { -1, 1, 2, 0, 3, -998};
  static const int cv12at7[]  = { -1, 3, 1, 2, 0, -998};
  static const int cv13at0[]  = { -1, 3, 2, 1, 0, -998};
  static const int cv13at1[]  = { -1, 3, 2, 0, 1, -998};
  static const int cv13at2[]  = { -1, 2, 3, 1, 0, -998};
  static const int cv13at3[]  = { -1, 2, 3, 0, 1, -998};
  static const int cv13at4[]  = { -1, 1, 0, 3, 2, -998};
  static const int cv13at5[]  = { -1, 0, 1, 3, 2, -998};
  static const int cv13at6[]  = { -1, 1, 0, 2, 3, -998};
  static const int cv13at7[]  = { -1, 0, 1, 2, 3, -998};
  static const int cv14at0[]  = { -1, 3, 1, 2, 0, -998};
  static const int cv14at1[]  = { -1, 3, 1, 0, 2, -998};
  static const int cv14at2[]  = { -1, 1, 3, 2, 0, -998};
  static const int cv14at3[]  = { -1, 1, 3, 0, 2, -998};
  static const int cv14at4[]  = { -1, 2, 0, 3, 1, -998};
  static const int cv14at5[]  = { -1, 0, 2, 3, 1, -998};
  static const int cv14at6[]  = { -1, 2, 0, 1, 3, -998};
  static const int cv14at7[]  = { -1, 0, 2, 1, 3, -998};
  static const int cv15at0[]  = { -1, 2, 1, 3, 0, -998};
  static const int cv15at1[]  = { -1, 2, 1, 0, 3, -998};
  static const int cv15at2[]  = { -1, 1, 2, 3, 0, -998};
  static const int cv15at3[]  = { -1, 1, 2, 0, 3, -998};
  static const int cv15at4[]  = { -1, 3, 0, 2, 1, -998};
  static const int cv15at5[]  = { -1, 0, 3, 2, 1, -998};
  static const int cv15at6[]  = { -1, 3, 0, 1, 2, -998};
  static const int cv15at7[]  = { -1, 0, 3, 1, 2, -998};

  Selector <const ColourLines *> sel;

  if( diag->id() == -1 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv1at0, sizeof(cv1at0)/sizeof(int)),  &(diag1[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv1at1, sizeof(cv1at1)/sizeof(int)),  &(diag1[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv1at2, sizeof(cv1at2)/sizeof(int)),  &(diag1[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv1at3, sizeof(cv1at3)/sizeof(int)),  &(diag1[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv1at4, sizeof(cv1at4)/sizeof(int)),  &(diag1[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv1at5, sizeof(cv1at5)/sizeof(int)),  &(diag1[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv1at6, sizeof(cv1at6)/sizeof(int)),  &(diag1[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv1at7, sizeof(cv1at7)/sizeof(int)),  &(diag1[7]) );
  } 
  else if( diag->id() == -2 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv2at0, sizeof(cv2at0)/sizeof(int)),  &(diag2[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv2at1, sizeof(cv2at1)/sizeof(int)),  &(diag2[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv2at2, sizeof(cv2at2)/sizeof(int)),  &(diag2[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv2at3, sizeof(cv2at3)/sizeof(int)),  &(diag2[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv2at4, sizeof(cv2at4)/sizeof(int)),  &(diag2[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv2at5, sizeof(cv2at5)/sizeof(int)),  &(diag2[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv2at6, sizeof(cv2at6)/sizeof(int)),  &(diag2[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv2at7, sizeof(cv2at7)/sizeof(int)),  &(diag2[7]) );
  } 
  else if( diag->id() == -3 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv3at0, sizeof(cv3at0)/sizeof(int)),  &(diag3[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv3at1, sizeof(cv3at1)/sizeof(int)),  &(diag3[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv3at2, sizeof(cv3at2)/sizeof(int)),  &(diag3[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv3at3, sizeof(cv3at3)/sizeof(int)),  &(diag3[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv3at4, sizeof(cv3at4)/sizeof(int)),  &(diag3[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv3at5, sizeof(cv3at5)/sizeof(int)),  &(diag3[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv3at6, sizeof(cv3at6)/sizeof(int)),  &(diag3[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv3at7, sizeof(cv3at7)/sizeof(int)),  &(diag3[7]) );
  } 
  else if( diag->id() == -4 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv4at0, sizeof(cv4at0)/sizeof(int)),  &(diag4[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv4at1, sizeof(cv4at1)/sizeof(int)),  &(diag4[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv4at2, sizeof(cv4at2)/sizeof(int)),  &(diag4[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv4at3, sizeof(cv4at3)/sizeof(int)),  &(diag4[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv4at4, sizeof(cv4at4)/sizeof(int)),  &(diag4[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv4at5, sizeof(cv4at5)/sizeof(int)),  &(diag4[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv4at6, sizeof(cv4at6)/sizeof(int)),  &(diag4[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv4at7, sizeof(cv4at7)/sizeof(int)),  &(diag4[7]) );
  } 
  else if( diag->id() == -5 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv5at0, sizeof(cv5at0)/sizeof(int)),  &(diag5[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv5at1, sizeof(cv5at1)/sizeof(int)),  &(diag5[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv5at2, sizeof(cv5at2)/sizeof(int)),  &(diag5[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv5at3, sizeof(cv5at3)/sizeof(int)),  &(diag5[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv5at4, sizeof(cv5at4)/sizeof(int)),  &(diag5[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv5at5, sizeof(cv5at5)/sizeof(int)),  &(diag5[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv5at6, sizeof(cv5at6)/sizeof(int)),  &(diag5[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv5at7, sizeof(cv5at7)/sizeof(int)),  &(diag5[7]) );
  } 
  else if( diag->id() == -6 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv6at0, sizeof(cv6at0)/sizeof(int)),  &(diag6[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv6at1, sizeof(cv6at1)/sizeof(int)),  &(diag6[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv6at2, sizeof(cv6at2)/sizeof(int)),  &(diag6[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv6at3, sizeof(cv6at3)/sizeof(int)),  &(diag6[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv6at4, sizeof(cv6at4)/sizeof(int)),  &(diag6[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv6at5, sizeof(cv6at5)/sizeof(int)),  &(diag6[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv6at6, sizeof(cv6at6)/sizeof(int)),  &(diag6[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv6at7, sizeof(cv6at7)/sizeof(int)),  &(diag6[7]) );
  } 
  else if( diag->id() == -7 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv7at0, sizeof(cv7at0)/sizeof(int)),  &(diag7[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv7at1, sizeof(cv7at1)/sizeof(int)),  &(diag7[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv7at2, sizeof(cv7at2)/sizeof(int)),  &(diag7[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv7at3, sizeof(cv7at3)/sizeof(int)),  &(diag7[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv7at4, sizeof(cv7at4)/sizeof(int)),  &(diag7[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv7at5, sizeof(cv7at5)/sizeof(int)),  &(diag7[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv7at6, sizeof(cv7at6)/sizeof(int)),  &(diag7[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv7at7, sizeof(cv7at7)/sizeof(int)),  &(diag7[7]) );
  } 
  else if( diag->id() == -8 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv8at0, sizeof(cv8at0)/sizeof(int)),  &(diag8[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv8at1, sizeof(cv8at1)/sizeof(int)),  &(diag8[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv8at2, sizeof(cv8at2)/sizeof(int)),  &(diag8[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv8at3, sizeof(cv8at3)/sizeof(int)),  &(diag8[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv8at4, sizeof(cv8at4)/sizeof(int)),  &(diag8[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv8at5, sizeof(cv8at5)/sizeof(int)),  &(diag8[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv8at6, sizeof(cv8at6)/sizeof(int)),  &(diag8[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv8at7, sizeof(cv8at7)/sizeof(int)),  &(diag8[7]) );
  } 
  else if( diag->id() == -9 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv9at0, sizeof(cv9at0)/sizeof(int)),  &(diag9[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv9at1, sizeof(cv9at1)/sizeof(int)),  &(diag9[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv9at2, sizeof(cv9at2)/sizeof(int)),  &(diag9[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv9at3, sizeof(cv9at3)/sizeof(int)),  &(diag9[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv9at4, sizeof(cv9at4)/sizeof(int)),  &(diag9[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv9at5, sizeof(cv9at5)/sizeof(int)),  &(diag9[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv9at6, sizeof(cv9at6)/sizeof(int)),  &(diag9[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv9at7, sizeof(cv9at7)/sizeof(int)),  &(diag9[7]) );
  } 
  else if( diag->id() == -10 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv10at0, sizeof(cv10at0)/sizeof(int)),  &(diag10[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv10at1, sizeof(cv10at1)/sizeof(int)),  &(diag10[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv10at2, sizeof(cv10at2)/sizeof(int)),  &(diag10[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv10at3, sizeof(cv10at3)/sizeof(int)),  &(diag10[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv10at4, sizeof(cv10at4)/sizeof(int)),  &(diag10[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv10at5, sizeof(cv10at5)/sizeof(int)),  &(diag10[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv10at6, sizeof(cv10at6)/sizeof(int)),  &(diag10[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv10at7, sizeof(cv10at7)/sizeof(int)),  &(diag10[7]) );
  } 
  else if( diag->id() == -11 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv11at0, sizeof(cv11at0)/sizeof(int)),  &(diag11[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv11at1, sizeof(cv11at1)/sizeof(int)),  &(diag11[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv11at2, sizeof(cv11at2)/sizeof(int)),  &(diag11[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv11at3, sizeof(cv11at3)/sizeof(int)),  &(diag11[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv11at4, sizeof(cv11at4)/sizeof(int)),  &(diag11[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv11at5, sizeof(cv11at5)/sizeof(int)),  &(diag11[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv11at6, sizeof(cv11at6)/sizeof(int)),  &(diag11[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv11at7, sizeof(cv11at7)/sizeof(int)),  &(diag11[7]) );
  } 
  else if( diag->id() == -12 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv12at0, sizeof(cv12at0)/sizeof(int)),  &(diag12[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv12at1, sizeof(cv12at1)/sizeof(int)),  &(diag12[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv12at2, sizeof(cv12at2)/sizeof(int)),  &(diag12[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv12at3, sizeof(cv12at3)/sizeof(int)),  &(diag12[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv12at4, sizeof(cv12at4)/sizeof(int)),  &(diag12[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv12at5, sizeof(cv12at5)/sizeof(int)),  &(diag12[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv12at6, sizeof(cv12at6)/sizeof(int)),  &(diag12[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv12at7, sizeof(cv12at7)/sizeof(int)),  &(diag12[7]) );
  } 
  else if( diag->id() == -13 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv13at0, sizeof(cv13at0)/sizeof(int)),  &(diag13[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv13at1, sizeof(cv13at1)/sizeof(int)),  &(diag13[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv13at2, sizeof(cv13at2)/sizeof(int)),  &(diag13[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv13at3, sizeof(cv13at3)/sizeof(int)),  &(diag13[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv13at4, sizeof(cv13at4)/sizeof(int)),  &(diag13[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv13at5, sizeof(cv13at5)/sizeof(int)),  &(diag13[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv13at6, sizeof(cv13at6)/sizeof(int)),  &(diag13[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv13at7, sizeof(cv13at7)/sizeof(int)),  &(diag13[7]) );
  } 
  else if( diag->id() == -14 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv14at0, sizeof(cv14at0)/sizeof(int)),  &(diag14[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv14at1, sizeof(cv14at1)/sizeof(int)),  &(diag14[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv14at2, sizeof(cv14at2)/sizeof(int)),  &(diag14[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv14at3, sizeof(cv14at3)/sizeof(int)),  &(diag14[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv14at4, sizeof(cv14at4)/sizeof(int)),  &(diag14[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv14at5, sizeof(cv14at5)/sizeof(int)),  &(diag14[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv14at6, sizeof(cv14at6)/sizeof(int)),  &(diag14[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv14at7, sizeof(cv14at7)/sizeof(int)),  &(diag14[7]) );
  } 
  else if( diag->id() == -15 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv15at0, sizeof(cv15at0)/sizeof(int)),  &(diag15[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv15at1, sizeof(cv15at1)/sizeof(int)),  &(diag15[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv15at2, sizeof(cv15at2)/sizeof(int)),  &(diag15[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv15at3, sizeof(cv15at3)/sizeof(int)),  &(diag15[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv15at4, sizeof(cv15at4)/sizeof(int)),  &(diag15[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv15at5, sizeof(cv15at5)/sizeof(int)),  &(diag15[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv15at6, sizeof(cv15at6)/sizeof(int)),  &(diag15[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv15at7, sizeof(cv15at7)/sizeof(int)),  &(diag15[7]) );
  } 
  return sel;
}
