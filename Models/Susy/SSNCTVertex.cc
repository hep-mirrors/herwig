// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSNCTVertex class.
//

#include "SSNCTVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Models/StandardModel/StandardCKM.h"

using namespace Herwig;

SSNCTVertex::SSNCTVertex() : MX_(2.e16*GeV), 
			     sw_(0.), cw_(0.), mw_(ZERO), 
			     sb_(0.), cb_(0.), q2last_(), couplast_(0.),
			     leftlast_(0.), rightlast_(0.), idlast_(0),
			     epsilon_(0.) {
  orderInGem(1);
  orderInGs(0);
}

IBPtr SSNCTVertex::clone() const {
  return new_ptr(*this);
}

IBPtr SSNCTVertex::fullclone() const {
  return new_ptr(*this);
}

void SSNCTVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(MX_,GeV) << nmix_ << sw_ << cw_ << ounit(mw_,GeV)
     << sb_ << cb_ << ounit(q2last_,GeV2) << couplast_
     << leftlast_ << rightlast_ << idlast_ << epsilon_;
}

void SSNCTVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(MX_,GeV) >> nmix_ >> sw_ >> cw_ >> iunit(mw_,GeV)
     >> sb_ >> cb_ >> iunit(q2last_,GeV2) >> couplast_
     >> leftlast_ >> rightlast_ >> idlast_ >> epsilon_;
}

ClassDescription<SSNCTVertex> SSNCTVertex::initSSNCTVertex;
// Definition of the static class description member.

void SSNCTVertex::Init() {

  static ClassDocumentation<SSNCTVertex> documentation
    ("The SSNCTVertex class implements the flavour violating"
     " coupling of the top quark to a charm quark and a neutralino");

  static Parameter<SSNCTVertex,Energy> interfaceMX
    ("MX",
     "Unification scale for the loop",
     &SSNCTVertex::MX_, GeV, 2.e16*GeV, 1.e14*GeV, 1.e20*GeV,
     false, false, Interface::limited);

}

void SSNCTVertex::doinit() {
  long neut[5] = {1000022, 1000023, 1000025, 1000035, 1000045};
  for(unsigned int nl = 0; nl < 5; ++nl) {
    addToList( neut[nl],  4, -1000006 );
    addToList( neut[nl], -4,  1000006 );
  }
  FFSVertex::doinit();
  // get the MSSM
  MSSMPtr model = dynamic_ptr_cast<MSSMPtr>(generator()->standardModel());
  if(!model)
    throw InitException() << "SSNCTVertex::doinit() - "
			  << "The model pointer is null."
			  << Exception::abortnow;
  // standard SUSY couplings
  // neutralino mixing
  nmix_ = model->neutralinoMix();
  if(!nmix_) throw InitException() << "SSNCTVertex::doinit() "
				   << "The neutralino mixing matrix pointer is null."
				   << Exception::abortnow;
  sw_ = sqrt(sin2ThetaW());
  mw_ = getParticleData(24)->mass();
  double tb = model->tanBeta();
  cw_ = sqrt(1. - sqr(sw_));
  sb_ = tb/sqrt(1 + sqr(tb));
  cb_ = sqrt(1 - sqr(sb_));
  // susy breaking scale
  Energy mSUSY = 
    sqrt(max(sqr(getParticleData(ParticleID::Z0)->mass()),
	     model->Mq3L()*model->MtR()));
  // CKM factor
  ThePEG::Ptr<Herwig::StandardCKM>::transient_const_pointer 
    CKMptr = ThePEG::dynamic_ptr_cast< ThePEG::Ptr<Herwig::StandardCKM>::
    transient_const_pointer>(model->CKM());
  if(!CKMptr)
    throw Exception() << "Must have access to the Herwig::StandardCKM object"
		      << "for the CKM matrix in SSNCTVertex::doinit()"
		      << Exception::runerror;
  vector< vector<Complex > > CKM;
  CKM = CKMptr->getUnsquaredMatrix(generator()->standardModel()->families());
  // SM masses
  Energy mb = getParticleData(ParticleID::b)->mass();
  Energy mt = getParticleData(ParticleID::t)->mass();
  // squark masses
  Energy mt1 = getParticleData(1000006)->mass();
  Energy mcL = getParticleData(1000004)->mass();
  // mixing parameter
  Complex pre = sqr(weakCoupling(sqr(mSUSY)))/16./sqr(Constants::pi)*
    log(MX_/mSUSY)*sqr(double(mb/mw_))/sqr(cb_)*conj(CKM[2][2])*CKM[1][2];
  complex<Energy2> deltaL = -0.5*pre*(sqr(model->Mq2L())+sqr(model->Mq3L()) +
				      2.*model->Mh12()+2.*sqr(model->MbR()) +
				      2.*real(     model->bottomTrilinear()*
					      conj(model->bottomTrilinear())));
  complex<Energy2> deltaR =  pre*mt*conj(model->bottomTrilinear());
  if(abs(mt1-mcL)/abs(mt1+mcL)<1e-10) {
    epsilon_ = 0.;
  }
  else {
    epsilon_ = (deltaL*(*model->stopMix())(0,0)-
		deltaR*(*model->stopMix())(0,1))/(sqr(mt1)-sqr(mcL));
  }
}


void SSNCTVertex::setCoupling(Energy2 q2,tcPDPtr part1,
			      tcPDPtr part2,tcPDPtr ) {
  long ism(abs(part1->id())), ineut(abs(part2->id()));
  tcPDPtr smfermion = part1;
  if( ism / 1000000 == 1 )  {
    swap( ism, ineut);
    smfermion = part2;
  }
  if(q2!= q2last_ || couplast_==0.) {
    couplast_ = -sqrt(2)*weakCoupling(q2);
    q2last_=q2;
  }
  norm(couplast_);
  if( ineut != idlast_) {
    idlast_ = ineut;
    // determine neutralino
    unsigned nl(0);
    switch( ineut ) {
    case 1000022 : nl = 0;
      break;
    case 1000023 : nl = 1;
      break;
    case 1000025 : nl = 2;
      break;
    case 1000035 : nl = 3;
      break;
    case 1000045 : nl = 4;
      break;
    default : assert(false);
    }
    // common primed neutralino matrices
    Complex n2prime = (*nmix_)(nl,1)*cw_ - (*nmix_)(nl,0)*sw_;
    Complex n1prime = (*nmix_)(nl,0)*cw_ + (*nmix_)(nl,1)*sw_;
    tcPDPtr smf = getParticleData(ism);
    double qf = smf->charge()/eplus;
    //Complex bracketl = qf*sw_*( conj(n1prime) - sw_*conj(n2prime)/cw_ );
    double lambda(0.);
    //neutralino mixing element
    Complex nlf(0.);
    lambda = -0.5 + qf*sqr(sw_);
    nlf = (*nmix_)(nl,3);
    Complex bracketr = sw_*qf*n1prime - n2prime*lambda/cw_;
    leftlast_  = 0.;
    rightlast_ = epsilon_*bracketr;
  }
  //determine the helicity order of the vertex
  if( smfermion->id() < 0 ) {
    left(conj(rightlast_));
    right(conj(leftlast_));
  }
  else {
    left(leftlast_);
    right(rightlast_);
  }
}
