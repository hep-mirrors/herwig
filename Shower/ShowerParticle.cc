// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerParticle class.
//

#include "ShowerParticle.h"
#include "Pythia7/EventRecord/Particle.h"
#include "ShowerIndex.h"

#define ShowerCast(a) (dynamic_ptr_cast<ShowerParticlePtr>(a))
#define tShowerCast(a) (dynamic_ptr_cast<tShowerParticlePtr>(a))

using namespace Herwig;
using namespace Pythia7;

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

void ShowerParticle::addChildrenEvtRec(const tStepPtr sp) {
  tPPtr dum; 
  tParticleVector yet; 
  tParticleVector addCh; 
  yet.push_back( dynamic_ptr_cast<tPPtr>(this) );  
  while( !yet.empty() ) { 
    dum = yet.back(); 
    yet.pop_back(); 
    for ( ParticleVector::const_iterator cit = dum->children().begin(); 
	  cit != dum->children().end(); ++cit ) {      
      yet.push_back( dynamic_ptr_cast<tPPtr>(*cit) ); 
      addCh.push_back( *cit ); 
    }
    while( !addCh.empty() ) {       
      sp->addDecayNoCol(dum, addCh.back());    
      // sp->addDecayProduct(dum, dynamic_ptr_cast<tPPtr>(addCh.back()));    
      addCh.pop_back(); 
    }
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


//bool ShowerParticle::addChildren(const tCollecShoParPtr & inputChildren) {
//bool isOK = true;
//if (! children().empty() ) {
    // When the current particle has already children (this would happen,
    // under normal conditions without bugs, if that particle has decayed 
    // before showering) then find between the (new)  inputChildren
    // which one has the same id as the current (decaying) particle:
    // this should be present because it represents the emitting
    // particle after the emission. If you can't find it then some
    // bugs is lurking somewhere; otherwise, transfer the original
    // children (decay products) of the current particle to the
    // one with the same identity in the (new)  inputChildren.
//  long idParticle = data().id();
//  tShowerParticlePtr partAfterEmission = tShowerParticlePtr();
//  for ( tShowerParticleVector::const_iterator cit = inputChildren.begin();
//  cit != inputChildren.end(); ++cit ) {
//    if ( (*cit)->data().id() == idParticle ) {
//partAfterEmission = *cit;
//    } 
//  }    
//  if ( partAfterEmission  &&  partAfterEmission->children().empty() ) {
//    for ( ShowerParticleVector::const_iterator cit = children().begin();
//    cit != children().end(); ++cit ) {
//partAfterEmission->addChild( *cit );
//    }    
//    _children.clear();
// } else {  
//    isOK = false;  // something wrong.
//  }
//}
//if ( isOK ) {
//  // Add normally the (new) children one by one, setting their parent pointers.
//  for ( tShowerParticleVector::const_iterator cit = inputChildren.begin();
//  cit != inputChildren.end(); ++cit ) {
//    addChild( *cit );
//  }
//}
//return isOK;
//}

