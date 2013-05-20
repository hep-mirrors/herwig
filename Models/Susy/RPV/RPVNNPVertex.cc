// -*- C++ -*-
//
// RPVNNPVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RPVNNPVertex class.
//

#include "RPVNNPVertex.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Models/Susy/MixingMatrix.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig++/Looptools/clooptools.h"
#include "Herwig++/Utilities/Maths.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

namespace {

  unsigned int neutralinoIndex(long id) {
    if(id> 1000000)
      return id<1000025 ? id-1000022 : (id-1000005)/10;
    else if(abs(id)<=16) 
      return (abs(id)-4)/2;
    else
      return id-13;
  }
  
  unsigned int charginoIndex(long id) {
    return abs(id)>1000000 ? (abs(id)-1000024)/13 : (abs(id)-7)/2;
  }

}

RPVNNPVertex::RPVNNPVertex() : _includeOnShell(false),
			     sw_(0.), cw_(0.), id1Last_(0), 
			     id2Last_(0), q2Last_(ZERO), coupLast_(0.),
			     leftLast_(ZERO), rightLast_(ZERO) {
  orderInGem(3);
  orderInGs(0);
}

void RPVNNPVertex::doinit() {
  tMSSMPtr model = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
  if(!model)
    throw InitException() << "RPVNNPVertex::doinit() - "
			  << "The model pointer is null."
			  << Exception::abortnow;
  // neutralino and chargino mixing matrix
  Nmix_ = model->neutralinoMix();
  Umix_ = model->charginoUMix();
  Vmix_ = model->charginoVMix();
  if(!Nmix_ || !Umix_ || ! Vmix_)
    throw InitException() << "RPVNNPVertex::doinit - The neutralino "
			  << "mixing matrix pointer is null." 
			  << Exception::abortnow;
  Looptools::ltini();
  vector<long> ineu(4);
  ineu[0] = 1000022; ineu[1] = 1000023;
  ineu[2] = 1000025; ineu[3] = 1000035;
  if(Nmix_->size().first==7) {
    if(model->majoranaNeutrinos()) {
      ineu.push_back(17);
      ineu.push_back(18);
      ineu.push_back(19);
    }
    else {
      ineu.push_back(12);
      ineu.push_back(14);
      ineu.push_back(16);
    }
  }
  for(unsigned int i = 0; i < ineu.size(); ++i) {
    for(unsigned int j = 0; j < ineu.size(); ++j) {
      addToList(ineu[i], ineu[j], 22);
    }
  }
  GeneralFFVVertex::doinit();














  
  sw_ = sqrt(sin2ThetaW());
  cw_ = sqrt(1 - sw_*sw_);
  _mw = getParticleData(ParticleID::Wplus)->mass();
  double tb = model->tanBeta();
  sb_ = tb/sqrt(1 + sqr(tb));
  cb_ = sqrt(1 - sqr(sb_));
  stop_ = model->stopMix();
  sbot_ = model->sbottomMix();
  stau_ = model->stauMix();
  Looptools::ltexi();
}

void RPVNNPVertex::dofinish() {
  Looptools::ltexi();
  GeneralFFVVertex::dofinish();
}

void RPVNNPVertex::doinitrun() {
  Looptools::ltini();
  GeneralFFVVertex::doinitrun();
}

void RPVNNPVertex::persistentOutput(PersistentOStream & os) const {
  os << sw_ << cw_ << Nmix_ << ounit(_mw,GeV) << sb_ << cb_
     << stop_ << sbot_ << stau_ << Umix_ << Vmix_ << _includeOnShell;
}

void RPVNNPVertex::persistentInput(PersistentIStream & is, int) {
  is >> sw_ >> cw_ >> Nmix_ >> iunit(_mw,GeV) >> sb_ >> cb_
     >> stop_ >> sbot_ >> stau_ >> Umix_ >> Vmix_ >> _includeOnShell;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<RPVNNPVertex,Helicity::GeneralFFVVertex>
describeRPVNNPVertex("Herwig::RPVNNPVertex", "HwSusy.so");

void RPVNNPVertex::Init() {

  static ClassDocumentation<RPVNNPVertex> documentation
    ("The loop-mediated coupling of the photon to a pair of neutralinos");

  static Switch<RPVNNPVertex,bool> interfaceIncludeOnShellIntermediates
    ("IncludeOnShellIntermediates",
     "Whether or not to include on-shell intermediate states",
     &RPVNNPVertex::_includeOnShell, false, false, false);
  static SwitchOption interfaceIncludeOnShellIntermediatesYes
    (interfaceIncludeOnShellIntermediates,
     "Yes",
     "Include them",
     true);
  static SwitchOption interfaceIncludeOnShellIntermediatesNo
    (interfaceIncludeOnShellIntermediates,
     "No",
     "Don't incldue them",
     false);

}

void RPVNNPVertex::setCoupling(Energy2 q2, tcPDPtr part1,
#ifndef NDEBUG
			      tcPDPtr part2,tcPDPtr part3) {
#else
			      tcPDPtr part2,tcPDPtr) {
#endif
  static long chargino[5]={ParticleID::SUSY_chi_1plus,ParticleID::SUSY_chi_2plus,11,13,15};
  int o[2]={1,0};
  long in1 = part1->id();
  long in2 = part2->id();
  Energy Mj = part1->mass();
  Energy Mi = part2->mass();
  // checks of the particle ids
  assert(part3->id()==ParticleID::gamma);
  assert(in1 == ParticleID::SUSY_chi_10 || in1 == ParticleID::SUSY_chi_20 ||
	 in1 == ParticleID::SUSY_chi_30 || in1 == ParticleID::SUSY_chi_40 ||
	 in1==12 || in1==14 || in1==16 || in1==17 || in1==18 || in1==19 );
  assert(in2 == ParticleID::SUSY_chi_10 || in2 == ParticleID::SUSY_chi_20 ||
	 in2 == ParticleID::SUSY_chi_30 || in2 == ParticleID::SUSY_chi_40 ||
	 in2==12 || in2==14 || in2==16 || in2==17 || in2==18 || in2==19 );
  // normal couplings are zero
  setLeft (0.);
  setRight(0.);
  if(in1==in2) {
    leftLast_  = ZERO;
    rightLast_ = ZERO;
    setLeftSigma (leftLast_ );
    setRightSigma(rightLast_);
    return;
  }
  if(q2 != q2Last_ || in1 != id1Last_ || in2 != id2Last_) {
    Looptools::clearcache();
    leftLast_  = ZERO;
    rightLast_ = ZERO;
    id1Last_ = in1;
    id2Last_ = in2;
    unsigned int neu1 = neutralinoIndex(in1);
    unsigned int neu2 = neutralinoIndex(in2);
    Complex n1prime[2] = { (*Nmix_)(neu2,0)*cw_ + (*Nmix_)(neu2,1)*sw_ ,
    			   (*Nmix_)(neu1,0)*cw_ + (*Nmix_)(neu1,1)*sw_ };
    Complex n2prime[2] = { (*Nmix_)(neu2,1)*cw_ - (*Nmix_)(neu2,0)*sw_ ,
    			   (*Nmix_)(neu1,1)*cw_ - (*Nmix_)(neu1,0)*sw_ };
    // sfermions
    unsigned int imax = stau_ ? 16 : 7;
    // sfermion/fermion loops
    for(long iferm=1;iferm<imax;++iferm) {
      // cerr << "testing in fermion loop " << iferm << "\n";
      if(iferm==7) iferm=11;
      if(iferm%2==0&&iferm>11) ++iferm;
      tcPDPtr smf = getParticleData(iferm);
      Energy mf = smf->mass();
      double qf = smf->charge()/eplus;
      double y = 0.5*dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel())->mass(q2,smf)/_mw;
      Complex bracketl[2] = { qf*sw_*( conj(n1prime[0]) - sw_*conj(n2prime[0])/cw_ ) ,
    			      qf*sw_*( conj(n1prime[1]) - sw_*conj(n2prime[1])/cw_ ) };
      double lambda(0.);
      //neutralino mixing element
      Complex nlf[2]={0.,0.};
      if( iferm % 2 == 0 ) {
    	y /= sb_;
    	lambda = -0.5 + qf*sqr(sw_);
    	nlf[0] = (*Nmix_)(neu2,3);
    	nlf[1] = (*Nmix_)(neu1,3);
      }
      else { 
    	y /= cb_;
    	lambda = 0.5 + qf*sqr(sw_);
    	nlf[0] = (*Nmix_)(neu2,2);
    	nlf[1] = (*Nmix_)(neu1,2);
      }
      Complex bracketr[2] = { sw_*qf*n1prime[0] - n2prime[0]*lambda/cw_ ,
    			      sw_*qf*n1prime[1] - n2prime[1]*lambda/cw_ };
      for(long iy=0;iy<2;++iy) {
    	long isf = 1000000*(1+iy)+iferm;
    	Energy msf = getParticleData(isf)->mass();
    	if(!_includeOnShell&&(mf+msf<Mj||mf+msf<Mi)) continue;
    	Complex g[2][2];
    	Complex ma1(0.), ma2(0.);
    	// heavy fermions
    	if( iferm == 5 || iferm == 6 || iferm == 15 ) {
    	  if( iferm == 5 ) {
    	    ma1 = (*sbot_)(iy,0);
    	    ma2 = (*sbot_)(iy,1); 
    	  } 
    	  else if( iferm == 6 ) {
    	    ma1 = (*stop_)(iy,0);
    	    ma2 = (*stop_)(iy,1);
    	  } 
    	  else {
    	    ma1 = (*stau_)(iy,0);
    	    ma2 = (*stau_)(iy,1);
    	  }
    	}
    	else if(iy==0) {
    	  ma1 = 1.;
    	}
    	else {
    	  ma2 = 1.;
    	}
    	for(unsigned int ix=0;ix<2;++ix) {
    	  g[ix][0] = y*conj(nlf[ix])*ma1 - ma2*bracketl[ix];
    	  g[ix][1] = y*nlf[ix]*ma2 + ma1*bracketr[ix];
    	  // cerr << "testing couplings " 
    	  //      << sqrt(2.)*weakCoupling(q2)*g[ix][0] << " " << sqrt(2.)*weakCoupling(q2)*g[ix][1] << "\n";
	  // cerr << "testing left/right " << ix<< " " << g[ix][0] << " " << g[ix][1] << "\n";
    	}
    	swap(g[0][0],g[0][1]);
    	complex<InvEnergy2> I,J,K,I2;
    	loopIntegrals(Mi,Mj,msf,mf,I,J,K,I2);
    	// cerr << "testing sfermion " << iy << " " << msf/GeV << " " << mf/GeV << "\n";
    	// cerr << "testing loop " << I*GeV2 << " " << J*GeV2 << " " << K*GeV2 << " " << I2*GeV2 << "\n";
    	complex<InvEnergy> coup[2];
    	for(unsigned int ix=0;ix<2;++ix) {
    	  coup[ix] =
	    Mj*(I2-K)*(g[0][ix]*g[1][o[ix]]-conj(g[0][o[ix]]*g[1][ix]))
    	    +Mi*K    *(g[0][o[ix]]*g[1][ix]-conj(g[0][ix]*g[1][o[ix]]))
    	    +mf*I    *(g[0][ix]*g[1][ix]-conj(g[0][o[ix]]*g[1][o[ix]]));
    	  // cerr << "testing coupling !!! " << ix << " " 
    	  //      << 4.*sqr(weakCoupling(q2))*coup[ix]*GeV << "\n";
    	}
    	double fact = 4.*qf;
    	if(iferm<=6) fact *=3.;

	// don't know about this
	fact *= 0.5;

    	leftLast_  += fact*coup[0];
	rightLast_ += fact*coup[1];
    	// cerr << "testing final " 
    	//      << 0.25*fact*sqr(weakCoupling(q2))*coup[0]*GeV  << " "
    	//      << 0.25*fact*sqr(weakCoupling(q2))*coup[1]*GeV   << "\n"
	  ;
      }
    }
    // cerr << "testing sfermion " 
    // 	 << 0.25*sqr(weakCoupling(q2))* leftLast_*GeV  << " "
    // 	 << 0.25*sqr(weakCoupling(q2))*rightLast_*GeV   << "\n";
    //leftLast_ = rightLast_ = ZERO;

    // the chargino W contribution
    for(unsigned int ic=0;ic<Vmix_->size().first;++ic) {
      long id = chargino[ic];
      Energy Mk = getParticleData(id)->mass();
      // cerr << "testing chargino " << Mk/GeV << "\n";
      if(!_includeOnShell&&(Mk+_mw<Mj||Mk+_mw<Mi)) continue;
      complex<InvEnergy2> I,J,K,I2;
      loopIntegrals(Mi,Mj,_mw,Mk,I,J,K,I2);
      Complex g[2][2];
      for(unsigned int ix=0;ix<2;++ix) {
     	unsigned int in = ix==0 ? neu2 : neu1;
 	g[ix][0] = 
 	  conj((*Nmix_)(in, 1))*(*Vmix_)(ic, 0) - 
 	  conj((*Nmix_)(in, 3))*(*Vmix_)(ic, 1)/sqrt(2);
 	g[ix][1] = 
 	  (*Nmix_)(in, 1)*conj((*Umix_)(ic, 0)) +
 	  (*Nmix_)(in, 2)*conj((*Umix_)(ic, 1))/sqrt(2);
	if(Vmix_->size().first==5) {
	  for(unsigned int k=0;k<3;++k)
	    g[ix][1] += ( (*Nmix_)(in, 4+k)*conj((*Umix_)(ic, 2+k))/sqrt(2));
	}
	// cerr << "testing left/right " << ix<< " " << g[ix][0] << " " << g[ix][1] << "\n";
      }
      complex<InvEnergy> coup[2];
      for(unsigned int ix=0;ix<2;++ix) {
    	coup[ix] = 
    	  Mj*(I2-J-K)*(g[0][o[ix]]*g[1][o[ix]]-conj(g[0][ix]*g[1][ix]))-
    	  Mi*(J-K)*(g[0][ix]*g[1][ix]-conj(g[0][o[ix]]*g[1][o[ix]]))+
    	  2.*Mk*J*(g[0][o[ix]]*g[1][ix]-conj(g[0][ix]*g[1][o[ix]]));
	// cerr << "testing coupling !!! pieces " << ix << " " 
	//      << sqr(weakCoupling(q2))*(g[0][o[ix]]*g[1][o[ix]]-conj(g[0][ix]*g[1][ix])) << " " 
	//      << sqr(weakCoupling(q2))*(g[0][ix]*g[1][ix]-conj(g[0][o[ix]]*g[1][o[ix]])) << " " 
	//      << sqr(weakCoupling(q2))*(g[0][o[ix]]*g[1][ix]-conj(g[0][ix]*g[1][o[ix]])) << "\n" ;
      }
      leftLast_  += 4.*coup[0];
      rightLast_ += 4.*coup[1];
    }
    // cerr << "testing total chargino " 
    // 	 << 0.25*sqr(weakCoupling(q2))* leftLast_*GeV  << " "
    // 	 << 0.25*sqr(weakCoupling(q2))*rightLast_*GeV   << "\n";



    // // the chargino charged higgs contribution
    // Energy mh = getParticleData(ParticleID::Hplus)->mass();
    // for(unsigned int ic=0;ic<2;++ic) {
    //   long id = ic==0 ? 
    // 	ParticleID::SUSY_chi_1plus : ParticleID::SUSY_chi_2plus;
    //   Energy Mk = getParticleData(id)->mass();
    //   if(!_includeOnShell&&(Mk+mh<Mj||Mk+mh<Mi)) continue;
    //   complex<InvEnergy2> I,J,K,I2;
    //   loopIntegrals(Mi,Mj,mh,Mk,I,J,K,I2);
    //   Complex g[2][2];
    //   for(unsigned int ix=0;ix<2;++ix) {
    // 	unsigned int in = ix==0 ? neu2 : neu1;
    // 	g[ix][0] =  (*Nmix_)(in, 3)*(*Vmix_)(ic,0) 
    // 	  +               ((*Nmix_)(in,1) + (*Nmix_)(in,0)*sw_/cw_)*
    // 	  (*Vmix_)(ic,1)/sqrt(2);
    // 	g[ix][0] *= cb_;
    // 	g[ix][1] = conj((*Nmix_)(in, 2)*(*Umix_)(ic,0) 
    // 			- ((*Nmix_)(in,1) + (*Nmix_)(in,0)*sw_/cw_)*
    // 			(*Umix_)(ic,1)/sqrt(2));
    // 	g[ix][1] *= sb_;
    //   }
    //   swap(g[1][0],g[1][1]);
    //   complex<InvEnergy> coup[2];
    //   for(unsigned int ix=0;ix<2;++ix) {
    // 	coup[ix] = Mj*(I2-K)*(g[0][ix]*g[1][o[ix]]-conj(g[0][o[ix]]*g[1][ix]))
    // 	  +Mi*K*(g[0][o[ix]]*g[1][ix]-conj(g[0][ix]*g[1][o[ix]]))
    // 	  +Mk*I*(g[0][ix]*g[1][ix]-conj(g[0][o[ix]]*g[1][o[ix]]));
    //   }
    //   leftLast_  += 2.*coup[0];
    //   rightLast_ += 2.*coup[1];
    // }
    // // the chargino goldstone contribution
    // for(unsigned int ic=0;ic<2;++ic) {
    //   long id = ic==0 ? 
    // 	ParticleID::SUSY_chi_1plus : ParticleID::SUSY_chi_2plus;
    //   Energy Mk = getParticleData(id)->mass();
    //   if(!_includeOnShell&&(Mk+_mw<Mj||Mk+_mw<Mi)) continue;
    //   complex<InvEnergy2> I,J,K,I2;
    //   loopIntegrals(Mi,Mj,_mw,Mk,I,J,K,I2);
    //   Complex g[2][2];
    //   for(unsigned int ix=0;ix<2;++ix) {
    // 	unsigned int in = ix==0 ? neu2 : neu1;
    // 	g[ix][0] = (*Nmix_)(in, 3)*(*Vmix_)(ic,0) 
    // 	  + ((*Nmix_)(in,1) + (*Nmix_)(in,0)*sw_/cw_)*
    // 	  (*Vmix_)(ic,1)/sqrt(2);
    // 	g[ix][0] *= sb_;
    // 	g[ix][1] = conj((*Nmix_)(in, 2)*(*Umix_)(ic,0) 
    // 			- ((*Nmix_)(in,1) + (*Nmix_)(in,0)*sw_/cw_)*
    // 			(*Umix_)(ic,1)/sqrt(2));
    // 	g[ix][1] *= cb_;
    //   }
    //   swap(g[1][0],g[1][1]);
    //   complex<InvEnergy> coup[2];
    //   for(unsigned int ix=0;ix<2;++ix) {
    // 	coup[ix] = Mj*(I2-K)*(g[0][ix]*g[1][o[ix]]-conj(g[0][o[ix]]*g[1][ix]))
    // 	  +Mi*K*(g[0][o[ix]]*g[1][ix]-conj(g[0][ix]*g[1][o[ix]]))
    // 	  +Mk*I*(g[0][ix]*g[1][ix]-conj(g[0][o[ix]]*g[1][o[ix]]));
    //   }
    //   leftLast_  += 2.*coup[0];
    //   rightLast_ += 2.*coup[1];
    // }
  }
  if(q2 != q2Last_ || coupLast_==0.) {
    q2Last_ = q2;
    coupLast_ = sqr(weakCoupling(q2))*
      electroMagneticCoupling(q2)/32./sqr(Constants::pi);
  }
  norm(coupLast_);
  setLeftSigma ( leftLast_);
  setRightSigma(rightLast_);  
}

void RPVNNPVertex::loopIntegrals(Energy Mi, Energy Mj, Energy M, Energy m,
				 complex<InvEnergy2> & I, complex<InvEnergy2> & J,
				 complex<InvEnergy2> & K, complex<InvEnergy2> & I2) {
  Energy2 m2(sqr(m)),M2(sqr(M)),Mi2(sqr(Mi)),Mj2(sqr(Mj));
  double min2  = Mj2*UnitRemoval::InvE2;
  double mout2 = Mi2*UnitRemoval::InvE2;
  double mf2   = m2 *UnitRemoval::InvE2;
  double ms2   = M2 *UnitRemoval::InvE2;
  I  = Looptools::C0i(Looptools::cc0,min2,mout2,0.,mf2,ms2,mf2)*UnitRemoval::InvE2;
  J  = Looptools::C0i(Looptools::cc0,min2,mout2,0.,ms2,mf2,ms2)*UnitRemoval::InvE2;
  I2 =-Looptools::C0i(Looptools::cc1,min2,mout2,0.,mf2,ms2,mf2)*UnitRemoval::InvE2;
  K  = (1.+Complex(m2*I+M2*J-Mj2*I2))/(Mi2-Mj2);
}
