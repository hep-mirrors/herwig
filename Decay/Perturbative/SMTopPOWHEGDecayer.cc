// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMTopPOWHEGDecayer class.
//

#include "SMTopPOWHEGDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

SMTopPOWHEGDecayer::SMTopPOWHEGDecayer() {}

IBPtr SMTopPOWHEGDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr SMTopPOWHEGDecayer::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void SMTopPOWHEGDecayer::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void SMTopPOWHEGDecayer::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<SMTopPOWHEGDecayer,SMTopDecayer>
  describeHerwigSMTopPOWHEGDecayer("Herwig::SMTopPOWHEGDecayer", "HwPertrubativeDecay.so");

void SMTopPOWHEGDecayer::Init() {

  static ClassDocumentation<SMTopPOWHEGDecayer> documentation
    ("There is no documentation for the SMTopPOWHEGDecayer class");

}

HardTreePtr SMTopPOWHEGDecayer::generateHardest(ShowerTreePtr) {
  cerr << "in generate hardest\n";

  // herwig stuff to get t, b and W from shower tree

  // call Alix's function to calculate hardest emission

  unsigned int npoint=100000;
  ofstream file("dalitz.top");
  for(unsigned int ix=0;ix<npoint;++ix) {
    vector<Lorentz5Momentum> momenta = hardMomenta();
    if(momenta.size()==4) { 
      double x_g = 2.*momenta[3].e()/momenta[0].mass();
      double x_w = 2.*momenta[2].e()/momenta[0].mass();
      file << x_g << "\t" << x_w << "\n";
    }
  }
  file.close();



  // Herwig stuff to put it in the hard tree and return it

  exit(1);
}

vector<Lorentz5Momentum>  SMTopPOWHEGDecayer::hardMomenta() {


  return vector<Lorentz5Momentum>();
}
