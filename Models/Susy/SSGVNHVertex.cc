// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSGVNHVertex class.
//

#include "SSGVNHVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "MSSM.h"

using namespace Herwig;

SSGVNHVertex::SSGVNHVertex() : sa_(0.), sb_(0.), ca_(0.), cb_(0.),
			       MPlanck_(2.4e18*GeV) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::SINGLET);
}

IBPtr SSGVNHVertex::clone() const {
  return new_ptr(*this);
}

IBPtr SSGVNHVertex::fullclone() const {
  return new_ptr(*this);
}

void SSGVNHVertex::persistentOutput(PersistentOStream & os) const {
  os << sa_ << sb_ << ca_ << cb_ << nmix_ << ounit(MPlanck_,GeV);
}

void SSGVNHVertex::persistentInput(PersistentIStream & is, int) {
  is >> sa_ >> sb_ >> ca_ >> cb_ >> nmix_ >> iunit(MPlanck_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SSGVNHVertex,Helicity::RFSVertex>
describeHerwigSSGVNHVertex("Herwig::SSGVNHVertex", "HwSusy.so");

void SSGVNHVertex::Init() {

  static ClassDocumentation<SSGVNHVertex> documentation
    ("The SSGVNHVertex class implments the coupling of the Higgs"
     " bosons to a gravitino and a neutralino");

}

void SSGVNHVertex::doinit() {
  long neu[4] = {ParticleID::SUSY_chi_10, ParticleID::SUSY_chi_20,
		 ParticleID::SUSY_chi_30, ParticleID::SUSY_chi_40};
  long higgs[3] =  {25, 35, 36};
  for(unsigned int i = 0; i < 3; ++i) {
    for(unsigned int j = 0; j < 4; ++j) {
      addToList(ParticleID::SUSY_Gravitino, neu[j], higgs[i]);
    }
  }
  RFSVertex::doinit();
  tMSSMPtr model = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
  if( !model )
    throw InitException() 
      << "SSGVNHVertex::doinit() - The pointer to the MSSM object is null!"
      << Exception::abortnow;
  double tanb = model->tanBeta();
  sb_ = tanb/sqrt(1. + sqr(tanb));
  cb_ = sqrt( 1. - sqr(sb_) );
  sa_ = sin(model->higgsMixingAngle());
  ca_ = sqrt(1. - sqr(sa_));
  nmix_ = model->neutralinoMix();
  MPlanck_ = model->MPlanck();
}

void SSGVNHVertex::setCoupling(Energy2 ,
#ifndef NDEBUG
			       tcPDPtr part1,
#else
			       tcPDPtr,
#endif
			       tcPDPtr part2,tcPDPtr part3) {
  assert(part1->id()==ParticleID::SUSY_Gravitino);
  assert(part3->iSpin()==PDT::Spin0);
  unsigned int neut = part2->id() - ParticleID::SUSY_chi_10;
  if(neut>1) neut = ( neut == 13 ) ? 3 : 2;
  int hid = part3->id();
  Complex coup;
  switch(hid) {
  case ParticleID::h0 :
    left (1.);
    right(1.);
    coup = -(*nmix_)(neut,2)*sa_+(*nmix_)(neut,3)*ca_;
    break;
  case ParticleID::H0 :
    left (1.);
    right(1.);
    coup =  (*nmix_)(neut,2)*ca_+(*nmix_)(neut,3)*sa_;
    break;
  case ParticleID::A0 :
    left (Complex(0.,-1.));
    right(Complex(0., 1.));
    coup =  (*nmix_)(neut,2)*sb_+(*nmix_)(neut,3)*cb_;
    break;
  default :
    assert(false);
  }
  norm(coup/MPlanck_*UnitRemoval::E);
}
