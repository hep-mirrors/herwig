// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GravitonMassGenerator class.
//

#include "GravitonMassGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/PDT/GenericWidthGenerator.h"
#include "ADDModel.h"

using namespace Herwig;

GravitonMassGenerator::GravitonMassGenerator()
  : prefactor_(0.), delta_(2), md_(1000.*GeV), mMin_(MeV) 
{}

IBPtr GravitonMassGenerator::clone() const {
  return new_ptr(*this);
}

IBPtr GravitonMassGenerator::fullclone() const {
  return new_ptr(*this);
}

void GravitonMassGenerator::persistentOutput(PersistentOStream & os) const {
  os << prefactor_ << delta_ << ounit(md_,GeV) << ounit(mMin_,GeV);
}

void GravitonMassGenerator::persistentInput(PersistentIStream & is, int) {
  is >> prefactor_ >> delta_ >> iunit(md_,GeV) >> iunit(mMin_,GeV);
}

ClassDescription<GravitonMassGenerator> 
GravitonMassGenerator::initGravitonMassGenerator;
// Definition of the static class description member.

void GravitonMassGenerator::Init() {

  static ClassDocumentation<GravitonMassGenerator> documentation
    ("The GravitonMassGenerator class generates the mass for external gravitions "
     "in the ADD model.");

  static Parameter<GravitonMassGenerator,Energy> interfaceMinimumMass
    ("MinimumMass",
     "Minimum gravition mass to avoid numerical problems",
     &GravitonMassGenerator::mMin_, GeV, MeV, MeV, GeV,
     false, false, Interface::limited);

}

void GravitonMassGenerator::doinit() {
  GenericMassGenerator::doinit();
  tcHwADDPtr hwADD = dynamic_ptr_cast<tcHwADDPtr>(generator()->standardModel());
  if(!hwADD) 
    throw Exception() << "Must have ADDModel in GravitonMassGenerator::doinit()"
		      << Exception::runerror;
  delta_ = hwADD->delta();
  md_ =  hwADD->MD();
  // calculate the prefactor
  prefactor_ = sqr(hwADD->MPlanckBar()/md_);
  // even no of dimensions
  if(delta_%2==0) {
    unsigned int n = delta_/2;
    prefactor_ *= 2.*pow(Constants::pi,int(n));
    for(unsigned int ix=1;ix<n;++ix) {
      prefactor_ /= double(ix);
    }
  }
  // odd number of dimensions
  else {
    unsigned int n = (delta_-1)/2;
    prefactor_ *= 2.*pow(Constants::pi,int(n));
    for(unsigned int ix=0;ix<n;++ix) {
      prefactor_ /= double(ix)+0.5;
    }
  }
}

Energy GravitonMassGenerator::mass(double & wgt, const ParticleData & ,
				   const Energy low,const Energy upp, int,
				   double r) const {
  Energy low2 = max(mMin_,low);
  double rlow = pow(double(low2/md_),int(delta_))/double(delta_);
  double rupp = pow(double(upp /md_),int(delta_))/double(delta_);
  double rho = rlow + (rupp-rlow)*r;
  wgt = (rupp-rlow)*prefactor_;
  return pow(double(delta_)*rho,1./double(delta_))*md_;
}
