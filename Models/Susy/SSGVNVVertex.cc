// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSGVNVVertex class.
//

#include "SSGVNVVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "MSSM.h"

using namespace Herwig;

SSGVNVVertex::SSGVNVVertex() : sw_(0.), cw_(0.), sb_(0.), cb_(0.),
			       mz_(91.1876*GeV), MPlanck_(2.4e18*GeV) {
  orderInGem(1);
  orderInGs(0);
}

IBPtr SSGVNVVertex::clone() const {
  return new_ptr(*this);
}

IBPtr SSGVNVVertex::fullclone() const {
  return new_ptr(*this);
}

void SSGVNVVertex::persistentOutput(PersistentOStream & os) const {
  os << sw_ << cw_ << sb_ << cb_ << ounit(mz_,GeV) << nmix_ << ounit(MPlanck_,GeV);
}

void SSGVNVVertex::persistentInput(PersistentIStream & is, int) {
  is >> sw_ >> cw_ >> sb_ >> cb_ >> iunit(mz_,GeV) >> nmix_ >> iunit(MPlanck_,GeV);
}

ClassDescription<SSGVNVVertex> SSGVNVVertex::initSSGVNVVertex;
// Definition of the static class description member.

void SSGVNVVertex::Init() {

  static ClassDocumentation<SSGVNVVertex> documentation
    ("The SSGVNVVertex class implements the coupling of the gravitino"
     " to the neutralino and a photon or Z boson, or the gluino and gluon.");

}

void SSGVNVVertex::doinit() {
  long neu[4] = {ParticleID::SUSY_chi_10, ParticleID::SUSY_chi_20,
		 ParticleID::SUSY_chi_30, ParticleID::SUSY_chi_40};
  for(unsigned int j = 0; j < 4; ++j) {
    addToList(ParticleID::SUSY_Gravitino, neu[j], ParticleID::gamma);
    addToList(ParticleID::SUSY_Gravitino, neu[j], ParticleID::Z0);
  }
  addToList(ParticleID::SUSY_Gravitino, ParticleID::SUSY_g, ParticleID::g);
  RFVVertex::doinit();
  tMSSMPtr model = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
  if( !model )
    throw InitException() 
      << "SSGVNVVertex::doinit() - The pointer to the MSSM object is null!"
      << Exception::abortnow;
  double tanb = model->tanBeta();
  sb_ = tanb/sqrt(1. + sqr(tanb));
  cb_ = sqrt( 1. - sqr(sb_) );
  sw_ = sqrt(sin2ThetaW());
  cw_ = sqrt(1. - sin2ThetaW());
  nmix_ = model->neutralinoMix();
  MPlanck_ = model->MPlanck();
}

void SSGVNVVertex::setCoupling(Energy2 ,
#ifndef NDEBUG
			       tcPDPtr part1,
#else
			       tcPDPtr,
#endif
			       tcPDPtr part2,tcPDPtr part3) {
  assert(part1->id()==ParticleID::SUSY_Gravitino);
  assert(part3->iSpin()==PDT::Spin1);
  unsigned int neut = part2->id() - ParticleID::SUSY_chi_10;
  if(neut>1) neut = ( neut == 13 ) ? 3 : 2;
  int bid = part3->id();
  Complex coup[2];
  vector<Complex> lV,rV;
  switch(bid) {
  case ParticleID::gamma :
    coup[0] = (*nmix_)(neut,0)*cw_+(*nmix_)(neut,1)*sw_;
    lV.push_back(-coup[0]*part2->mass()*UnitRemoval::InvE);
    lV.push_back( coup[0]);
    lV.push_back(0.);
    rV=lV;
    break;
  case ParticleID::Z0 :
    coup[0] = -(*nmix_)(neut,0)*sw_+(*nmix_)(neut,1)*cw_;
    coup[1] = -(*nmix_)(neut,2)*cb_+(*nmix_)(neut,3)*sb_;
    lV.push_back((-coup[0]*part2->mass()-mz_*coup[1])*UnitRemoval::InvE);
    lV.push_back( coup[0]);
    lV.push_back(0.);
    rV=lV;
    rV[0] = Complex((-coup[0]*part2->mass()+mz_*coup[1])*UnitRemoval::InvE);
    break;
  case ParticleID::g :
    lV.push_back(-double(part2->mass()*UnitRemoval::InvE));
    lV.push_back(1.);
    lV.push_back(0.);
    rV=lV;
    break;
  default :
    assert(false);
  }
  left (lV);
  right(rV);
  norm(double(1./MPlanck_*UnitRemoval::E));
}
