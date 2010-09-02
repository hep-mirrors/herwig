// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FortranShowerKinematics class.
//

#include "FortranShowerKinematics.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"

using namespace Herwig;

void FortranShowerKinematics::
updateChildren(const tShowerParticlePtr theParent, 
	       const ShowerParticleVector theChildren,
	       bool angularOrder) const {
  if(theChildren.size() != 2)
    throw Exception() << "FortranShowerKinematics1to2::updateChildren() " 
		      << "Warning! too many children!" << Exception::eventerror;
  // get the interaction type
  // if parent not from the shower set up
  if(theParent->showerParameters().empty()) {
    // DGRELL is resize intended here? The argument may not be 
    // used if there is a value already there.
    theParent->showerParameters().resize(1,1.);
    theParent->showerVariables().resize(5,0.*MeV);
    theParent->showerVariables()[0]=theParent->evolutionScale();
  }
  Energy dqtilde = scale();
  double dz = z(); 
//   double dphi = phi();
  theChildren[0]->showerVariables().resize(5,0.*MeV);
  theChildren[1]->showerVariables().resize(5,0.*MeV);
  // note that 1st child gets z, 2nd gets (1-z) by our convention.
  if(angularOrder) {
    theChildren[0]->setEvolutionScale(dz*dqtilde);
    theChildren[1]->setEvolutionScale((1.-dz)*dqtilde);
  }
  else {
    theChildren[0]->setEvolutionScale(dqtilde);
    theChildren[1]->setEvolutionScale(dqtilde);
  }
  // calculate and set xi
  double xi = sqr(dqtilde/theParent->showerVariables()[0]);
  theChildren[0]->showerParameters().resize(1,xi);
  theChildren[1]->showerParameters().resize(1,xi);
  // set energy
  theChildren[0]->showerVariables()[0]=    dz *theParent->showerVariables()[0];
  theChildren[1]->showerVariables()[0]=(1.-dz)*theParent->showerVariables()[0];
  // set up the colour connections
  splittingFn()->colourConnection(theParent,theChildren[0],theChildren[1],false);
  // make the products children of the parent
  theParent->addChild(theChildren[0]);
  theParent->addChild(theChildren[1]);
}

void FortranShowerKinematics::updateParent(const tShowerParticlePtr , 
					   const ShowerParticleVector ) const {
  throw Exception() << "FortranShowerKinematics::updateParent "
		    << " not implemented" << Exception::runerror;
}

void FortranShowerKinematics::updateLast(const tShowerParticlePtr ) const {
  throw Exception() << "FortranShowerKinematics::updateLast "
		    << " not implemented" << Exception::runerror;
}

void FortranShowerKinematics::reconstructChildren(const tShowerParticlePtr , 
						  const ShowerParticleVector ) const {
  throw Exception() << "FortranShowerKinematics::reconstructChildren "
		    << " not implemented" << Exception::runerror;
}

void FortranShowerKinematics::reconstructParent(const tShowerParticlePtr theParent, 
		       const ParticleVector theChildren) const {
  // calculate the mass of the parent
  tShowerParticlePtr c1=dynamic_ptr_cast<tShowerParticlePtr>(theChildren[0]);
  tShowerParticlePtr c2=dynamic_ptr_cast<tShowerParticlePtr>(theChildren[1]);
  Energy2 exi=c1->evolutionScale()*c2->evolutionScale();
  Energy2 m2=exi+sqr(theChildren[0]->mass())+sqr(theChildren[1]->mass());
  // calculate three-momentum of the parent
  Energy2 pisq=sqr(theParent->showerVariables()[0])-m2;
  theParent->showerVariables()[4]=sqrt(pisq);
  // Compute daughter' transverse and longitudinal momenta
  Energy2 pjpk=c1->showerVariables()[4]*c2->showerVariables()[4];
  Energy2 ejek=c1->showerVariables()[0]*c2->showerVariables()[0]-exi;
  Energy2 ptsq=(pjpk+ejek)*(pjpk-ejek)/pisq;
  c1->showerVariables()[3]= sqrt(ptsq);
  c2->showerVariables()[3]=-sqrt(ptsq);
  c1->showerVariables()[4]=sqrt(sqr(c1->showerVariables()[4])-ptsq);
  c2->showerVariables()[4]=theParent->showerVariables()[4]-c1->showerVariables()[4];
  theParent->set5Momentum(Energy(sqrt(m2)));
  cerr << "testing A " << theParent->PDGName() << " -> " 
       << c1->PDGName() << " " 
       << c2->PDGName() << "\n";
  cerr << "testing branching" << theParent->mass()/MeV << " " << c1->showerVariables()[3]/MeV
       << "\n";
}

void FortranShowerKinematics::reconstructLast(const tShowerParticlePtr theLast,
					      unsigned int ) const {
  Energy mass=theLast->dataPtr()->constituentMass();
  theLast->set5Momentum(mass);
  Energy2 p2=sqr(theLast->showerVariables()[0])-sqr(mass);
  theLast->showerVariables()[4]=sqrt(p2);
  cerr << "testing last " << theLast->mass()/MeV << "\n";
}

void FortranShowerKinematics::initialize(ShowerParticle & ) {
}

vector<Lorentz5Momentum> FortranShowerKinematics::getBasis() const {
  throw Exception() << "The Fortran shower doesn't use an explicit basis so "
		    << "FortranShowerKinematics::getBasis() should not be called" 
		    << Exception::runerror;
}
