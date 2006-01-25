// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerParticle class.
//

#include "ShowerParticle.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ShowerIndex.h"

#define ShowerCast(a) (dynamic_ptr_cast<ShowerParticlePtr>(a))
#define tShowerCast(a) (dynamic_ptr_cast<tShowerParticlePtr>(a))

using namespace Herwig;
using namespace ThePEG;

ClassDescription<ShowerParticle> ShowerParticle::initShowerParticle;

void ShowerParticle::printInfo() {
  cout << this->data().id() << " " 
       << momentum().x() << " " << momentum().y() << " " 
       << momentum().z() << " " <<momentum().t() << " ";
  Lorentz5Momentum dum = sumParentsMomenta() - momentum(); 
  cout << dum.x() << " " << dum.y() << " " 
       << dum.z() << " " << dum.t() << " " << children().size() << "\t";
  if( parents().size() ) {
    cout << "<" << parents()[0]
	 << "," << parents()[0]->data().PDGName() << "> ";
  } else {
    cout << "<> ";
  }
  cout << "[" << this << "," << this->data().PDGName() << "] {"; 
  if ( children().size() > 0 ) {
    for (ParticleVector::const_iterator cit = children().begin(); 
	 cit != children().end(); ++cit) {
      cout << "[" << (*cit) << "," << (*cit)->data().PDGName() << "]";
    }
    cout << "} ";
  } else {
    cout << "} ";
  }
  cout << endl; 
  if( parents().size() ) {
    cout << "col " << parents()[0]->data().PDGName() << " > ";
  } else {
    cout << " > ";
  }
  cout << this->data().PDGName() << " > ("; 
  if ( children().size() > 0 ) {
    for (ParticleVector::const_iterator jit = children().begin(); 
	 jit != children().end(); ++jit) {
      cout << (*jit)->data().PDGName() << " ";
    }
    cout << ") [";
  } else {
    cout << " ) [";
  }
  cout << colourLine() << ", " 
       << antiColourLine() << "] ["
       << incomingColour() << ", " 
       << incomingAntiColour() << "] ["
       << outgoingColour() << ", " 
       << outgoingAntiColour() << "] ["
       << colourNeighbour() << ", " 
       << antiColourNeighbour() << "]" << endl; 
}

ShowerParticle::~ShowerParticle() { }

void ShowerParticle::deepPrintInfo() {
  printInfo();
  for (ParticleVector::const_iterator cit = children().begin(); 
       cit != children().end(); ++cit) {
    ShowerParticlePtr s = ShowerCast(*cit);
    if(s) s->deepPrintInfo();   
  }
}

Lorentz5Momentum ShowerParticle::sumParentsMomenta() {
  Lorentz5Momentum dum; 
  if (parents().size() == 1) {
    tShowerParticlePtr s = tShowerCast(parents()[0]);
    if (s) dum = s->sumParentsMomenta();
  }
  return dum += momentum(); 
}


tShowerParticleVector ShowerParticle::getFSChildren() {
  tShowerParticleVector fs;   
  if ( this->children().empty() ) {
    fs.push_back( this ); 
  } else {
    tShowerParticleVector yet; 
    yet.push_back( this ); 
    tShowerParticlePtr dum; 
    while (! yet.empty() ) {
      dum = yet.back(); 
      yet.pop_back(); 
      for (ParticleVector::const_iterator cit = dum->children().begin(); 
	   cit != dum->children().end(); ++cit) {
	if ( ShowerCast(*cit)->isReconstructionFixedPoint() ) {
	  fs.push_back( ShowerCast(*cit) ); 
	} else {
	  yet.push_back( ShowerCast(*cit) ); 
	}
      }
    } 
  }
  return fs; 
}	


void ShowerParticle::setShowerKinematics(const ShoKinPtr inputShowerKinematics) {
  _showerKinematics = inputShowerKinematics;
}

