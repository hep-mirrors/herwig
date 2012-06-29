// -*- C++ -*-
//
// SSNNPVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSNNPVertex class.
//

#include "SSNNPVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Models/Susy/MixingMatrix.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig++/Utilities/Maths.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSNNPVertex::SSNNPVertex() : _sw(0.), _cw(0.), _id1last(0), 
			     _id2last(0), _q2last(ZERO), _couplast(0.),
			     _leftlast(ZERO), _rightlast(ZERO) {
  orderInGem(3);
  orderInGs(0);
}

void SSNNPVertex::doinit() {
  int ineu[5] = {1000022,1000023,1000025,1000035,1000045};
  for(unsigned int i = 0; i < 5; ++i) {
    for(unsigned int j = 0; j < 5; ++j) {
      addToList(ineu[i], ineu[j], 22);
    }
  }
  GeneralFFVVertex::doinit();
  tMSSMPtr theSS = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
  if(!theSS)
    throw InitException() << "SSNNPVertex::doinit() - "
			  << "The model pointer is null."
			  << Exception::abortnow;
  
  _theN  = theSS->neutralinoMix();
  _theU = theSS->charginoUMix();
  _theV = theSS->charginoVMix();
  if(!_theN || !_theU || ! _theV)
    throw InitException() << "SSNNPVertex::doinit - The neutralino "
			  << "mixing matrix pointer is null." 
			  << Exception::abortnow;
  _sw = sqrt(sin2ThetaW());
  _cw = sqrt(1 - _sw*_sw);
  _mw = getParticleData(ParticleID::Wplus)->mass();
  double tb = theSS->tanBeta();
  _sb = tb/sqrt(1 + sqr(tb));
  _cb = sqrt(1 - sqr(_sb));
  _stop = theSS->stopMix();
  _sbot = theSS->sbottomMix();
  _stau = theSS->stauMix();
}

void SSNNPVertex::persistentOutput(PersistentOStream & os) const {
  os << _sw << _cw << _theN << ounit(_mw,GeV) << _sb << _cb
     << _stop << _sbot << _stau << _theU << _theV;
}

void SSNNPVertex::persistentInput(PersistentIStream & is, int) {
  is >> _sw >> _cw >> _theN >> iunit(_mw,GeV) >> _sb >> _cb
     >> _stop >> _sbot >> _stau >> _theU >> _theV;
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<SSNNPVertex,Helicity::GeneralFFVVertex>
describeSSNNPVertex("Herwig::SSNNPVertex", "HwSusy.so");

void SSNNPVertex::Init() {

  static ClassDocumentation<SSNNPVertex> documentation
    ("The coupling of a Z-boson to a pair of neutralinos");

}

void SSNNPVertex::setCoupling(Energy2 q2, tcPDPtr part1,
#ifndef NDEBUG
			      tcPDPtr part2,tcPDPtr part3) {
#else
			      tcPDPtr part2,tcPDPtr) {
#endif
  int o[2]={1,0};
  long in1 = part1->id();
  long in2 = part2->id();
  Energy Mj = part1->mass();
  Energy Mi = part2->mass();
  // checks of the particle ids
  assert(part3->id()==ParticleID::gamma);
  assert(in1 == ParticleID::SUSY_chi_10 || in1 == ParticleID::SUSY_chi_20 ||
	 in1 == ParticleID::SUSY_chi_30 || in1 == ParticleID::SUSY_chi_40 || 
	 in1 == 1000045                 );
  assert(in2 == ParticleID::SUSY_chi_10 || in2 == ParticleID::SUSY_chi_20 ||
	 in2 == ParticleID::SUSY_chi_30 || in2 == ParticleID::SUSY_chi_40 ||
	 in2 == 1000045              );
  // normal couplings are zero
  setLeft (0.);
  setRight(0.);
  if(in1==in2) {
    _leftlast  = ZERO;
    _rightlast = ZERO;
    setLeftSigma (_leftlast );
    setRightSigma(_rightlast);
    return;
  }
  if(q2 != _q2last || _couplast==0.) {
    _q2last = q2;
    _couplast = sqr(weakCoupling(q2))*
      electroMagneticCoupling(q2)/32./sqr(Constants::pi);
  }
  if(in1 != _id1last || in2 != _id2last) {
    _leftlast  = ZERO;
    _rightlast = ZERO;
    _id1last = in1;
    _id2last = in2;
    unsigned int neu1(in1 - 1000022), neu2(in2 - 1000022);
    if(neu1 > 1) neu1 = (in1-1000005)/10;
    if(neu2 > 1) neu2 = (in2-1000005)/10;
    Complex n1prime[2] = { (*_theN)(neu2,0)*_cw + (*_theN)(neu2,1)*_sw ,
			   (*_theN)(neu1,0)*_cw + (*_theN)(neu1,1)*_sw };
    Complex n2prime[2] = { (*_theN)(neu2,1)*_cw - (*_theN)(neu2,0)*_sw ,
			   (*_theN)(neu1,1)*_cw - (*_theN)(neu1,0)*_sw };
    // sfermion/fermion loops
    for(long iferm=1;iferm<16;++iferm) {
      if(iferm==7) iferm=11;
      if(iferm%2==0&&iferm>11) ++iferm;
      tcPDPtr smf = getParticleData(iferm);
      Energy mf = smf->mass();
      double qf = smf->charge()/eplus;
      double y = 0.5*mf/_mw;
      Complex bracketl[2] = { qf*_sw*( conj(n1prime[0]) - _sw*conj(n2prime[0])/_cw ) ,
			      qf*_sw*( conj(n1prime[1]) - _sw*conj(n2prime[1])/_cw ) };
      double lambda(0.);
      //neutralino mixing element
      Complex nlf[2]={0.,0.};
      if( iferm % 2 == 0 ) {
 	y /= _sb;
 	lambda = -0.5 + qf*sqr(_sw);
 	nlf[0] = (*_theN)(neu2,3);
 	nlf[1] = (*_theN)(neu1,3);
      }
      else { 
	y /= _cb;
	lambda = 0.5 + qf*sqr(_sw);
	nlf[0] = (*_theN)(neu2,2);
	nlf[1] = (*_theN)(neu1,2);
      }
      Complex bracketr[2] = { _sw*qf*n1prime[0] - n2prime[0]*lambda/_cw ,
			      _sw*qf*n1prime[1] - n2prime[1]*lambda/_cw };
      for(long iy=0;iy<2;++iy) {
	long isf = 1000000*(1+iy)+iferm;
	Energy msf = getParticleData(isf)->mass();
	if(mf+msf<Mj||mf+msf<Mi) continue;
	Complex g[2][2];
	Complex ma1(0.), ma2(0.);
	// heavy fermions
	if( iferm == 5 || iferm == 6 || iferm == 15 ) {
	  if( iferm == 5 ) {
	    ma1 = (*_sbot)(iy,0);
	    ma2 = (*_sbot)(iy,1); 
	  } 
	  else if( iferm == 6 ) {
	    ma1 = (*_stop)(iy,0);
	    ma2 = (*_stop)(iy,1);
	  } 
	  else {
	    ma1 = (*_stau)(iy,0);
	    ma2 = (*_stau)(iy,1);
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
	}
	swap(g[0][0],g[0][1]);
	complex<InvEnergy2> I,J,K,I2;
	loopIntegrals(Mi,Mj,msf,mf,I,J,K,I2);
	complex<InvEnergy> coup[2];
	for(unsigned int ix=0;ix<2;++ix) {
	  coup[ix] = Mj*(I2-K)*(g[0][ix]*g[1][o[ix]]-conj(g[0][o[ix]]*g[1][ix]))
	    +Mi*K*(g[0][o[ix]]*g[1][ix]-conj(g[0][ix]*g[1][o[ix]]))
	    +mf*I*(g[0][ix]*g[1][ix]-conj(g[0][o[ix]]*g[1][o[ix]]));
	}
	double fact = 4.*qf;
	if(iferm<=6) fact *=3.;
	_leftlast  += fact*coup[0];
	_rightlast += fact*coup[1];
      }
    }
    // the chargino W contribution
    for(unsigned int ic=0;ic<2;++ic) {
      long id = ic==0 ? 
	ParticleID::SUSY_chi_1plus : ParticleID::SUSY_chi_2plus;
      Energy Mk = getParticleData(id)->mass();
      if(Mk+_mw<Mj||Mk+_mw<Mi) continue;
      complex<InvEnergy2> I,J,K,I2;
      loopIntegrals(Mi,Mj,_mw,Mk,I,J,K,I2);
      Complex g[2][2];
      for(unsigned int ix=0;ix<2;++ix) {
	unsigned int in = ix==0 ? neu2 : neu1;
	g[ix][0] = 
	  conj((*_theN)(in, 1))*(*_theV)(ic, 0) - 
	  conj((*_theN)(in, 3))*(*_theV)(ic, 1)/sqrt(2);
	g[ix][1] = 
	  (*_theN)(in, 1)*conj((*_theU)(ic, 0)) +
	  (*_theN)(in, 2)*conj((*_theU)(ic, 1))/sqrt(2);
      }
      complex<InvEnergy> coup[2];
      for(unsigned int ix=0;ix<2;++ix) {
	coup[ix] = 
	  Mj*(I2-J-K)*(g[0][o[ix]]*g[1][o[ix]]-conj(g[0][ix]*g[1][ix]))-
	  Mi*(J-K)*(g[0][ix]*g[1][ix]-conj(g[0][o[ix]]*g[1][o[ix]]))+
	  2.*Mk*J*(g[0][o[ix]]*g[1][ix]-conj(g[0][ix]*g[1][o[ix]]));
      }
      _leftlast  += 4.*coup[0];
      _rightlast += 4.*coup[1];
    }
    // the chargino charged higgs contribution
    Energy mh = getParticleData(ParticleID::Hplus)->mass();
    for(unsigned int ic=0;ic<2;++ic) {
      long id = ic==0 ? 
	ParticleID::SUSY_chi_1plus : ParticleID::SUSY_chi_2plus;
      Energy Mk = getParticleData(id)->mass();
      if(Mk+mh<Mj||Mk+mh<Mi) continue;
      complex<InvEnergy2> I,J,K,I2;
      loopIntegrals(Mi,Mj,mh,Mk,I,J,K,I2);
      Complex g[2][2];
      for(unsigned int ix=0;ix<2;++ix) {
	unsigned int in = ix==0 ? neu2 : neu1;
	g[ix][0] =  (*_theN)(in, 3)*(*_theV)(ic,0) 
	  +               ((*_theN)(in,1) + (*_theN)(in,0)*_sw/_cw)*
	  (*_theV)(ic,1)/sqrt(2);
	g[ix][0] *= _cb;
	g[ix][1] = conj((*_theN)(in, 2)*(*_theU)(ic,0) 
			- ((*_theN)(in,1) + (*_theN)(in,0)*_sw/_cw)*
			(*_theU)(ic,1)/sqrt(2));
	g[ix][1] *= _sb;
      }
      swap(g[1][0],g[1][1]);
      complex<InvEnergy> coup[2];
      for(unsigned int ix=0;ix<2;++ix) {
	coup[ix] = Mj*(I2-K)*(g[0][ix]*g[1][o[ix]]-conj(g[0][o[ix]]*g[1][ix]))
	  +Mi*K*(g[0][o[ix]]*g[1][ix]-conj(g[0][ix]*g[1][o[ix]]))
	  +Mk*I*(g[0][ix]*g[1][ix]-conj(g[0][o[ix]]*g[1][o[ix]]));
      }
      _leftlast  += 2.*coup[0];
      _rightlast += 2.*coup[1];
    }
    // the chargino goldstone contribution
    for(unsigned int ic=0;ic<2;++ic) {
      long id = ic==0 ? 
	ParticleID::SUSY_chi_1plus : ParticleID::SUSY_chi_2plus;
      Energy Mk = getParticleData(id)->mass();
      if(Mk+_mw<Mj||Mk+_mw<Mi) continue;
      complex<InvEnergy2> I,J,K,I2;
      loopIntegrals(Mi,Mj,_mw,Mk,I,J,K,I2);
      Complex g[2][2];
      for(unsigned int ix=0;ix<2;++ix) {
	unsigned int in = ix==0 ? neu2 : neu1;
	g[ix][0] = (*_theN)(in, 3)*(*_theV)(ic,0) 
	  + ((*_theN)(in,1) + (*_theN)(in,0)*_sw/_cw)*
	  (*_theV)(ic,1)/sqrt(2);
	g[ix][0] *= _sb;
	g[ix][1] = conj((*_theN)(in, 2)*(*_theU)(ic,0) 
			- ((*_theN)(in,1) + (*_theN)(in,0)*_sw/_cw)*
			(*_theU)(ic,1)/sqrt(2));
	g[ix][1] *= _cb;
      }
      swap(g[1][0],g[1][1]);
      complex<InvEnergy> coup[2];
      for(unsigned int ix=0;ix<2;++ix) {
	coup[ix] = Mj*(I2-K)*(g[0][ix]*g[1][o[ix]]-conj(g[0][o[ix]]*g[1][ix]))
	  +Mi*K*(g[0][o[ix]]*g[1][ix]-conj(g[0][ix]*g[1][o[ix]]))
	  +Mk*I*(g[0][ix]*g[1][ix]-conj(g[0][o[ix]]*g[1][o[ix]]));
      }
      _leftlast  += 2.*coup[0];
      _rightlast += 2.*coup[1];
    }
  }
  norm(_couplast);
  setLeftSigma ( _leftlast);
  setRightSigma(_rightlast);  
}

void SSNNPVertex::loopIntegrals(Energy Mi, Energy Mj, Energy M, Energy m,
				complex<InvEnergy2> & I, complex<InvEnergy2> & J,
				complex<InvEnergy2> & K, complex<InvEnergy2> & I2) {
  static const Complex ii(0.,1.);
  static const Energy eps(100.*MeV);
  Energy2 m2(sqr(m)),M2(sqr(M)),Mi2(sqr(Mi)),Mj2(sqr(Mj));
  using Math::Li2;
  // general form
  if(m>eps) {
    Energy4 li = sqr(m2+M2-Mi2)-4.*sqr(m*M);
    complex<Energy2> rli = li<ZERO ? ii*sqrt(-li) : sqrt(li); 
    Energy4 lj = sqr(m2+M2-Mj2)-4.*sqr(m*M);
    complex<Energy2> rlj = lj<ZERO ? ii*sqrt(-lj) : sqrt(lj);
    Complex arg[6]={0.5/m2*(Mj2+m2-M2+rlj) ,0.5/m2*(Mj2+m2-M2-rlj),
		    0.5/m2*(Mi2+m2-M2+rli) ,0.5/m2*(Mi2+m2-M2-rli),
		    0.5/m/M*(m2+M2-Mj2+rlj),0.5/m/M*(m2+M2-Mi2+rli)};
    I  = 1./(Mi2-Mj2)*(Li2(arg[0])+Li2(arg[1])-
		       Li2(arg[2])-Li2(arg[3]));
    J  = 1./(Mj2-Mi2)*(sqr(log(arg[4]))-sqr(log(arg[5])))-I;
    Complex Itest[2];
    for(unsigned int ix=0;ix<2;++ix) {
      Complex a,b;
      if(ix==0) {
	a = 0.5*(M2+Mj2-m2)/Mj2;
	b = 0.5*rlj/Mj2;
      }
      else {
	a = 0.5*(M2+Mi2-m2)/Mi2; 
	b = 0.5*rli/Mi2;
      }
      Itest[ix] = -b*log(b-a)+a*log(b-a)+log(-a-b)*a
	+log(-a-b)*b+log(1.+b-a)+log(1.+b-a)*b-log(1.+b-a)*a+log(1.-a-b)
	-log(1.-a-b)*a-log(1.-a-b)*b-2.;
    }
    I2 = (Itest[0]-Itest[1]+log(Mj2/Mi2))/(Mj2-Mi2);
//     I2 = (M2-m2)/Mi2/Mj2*log(m/M)
//       +1./(Mj2-Mi2)*(0.5*rlj/Mj2*log((m2+M2-Mj2-rlj)/
// 				     (m2+M2-Mj2+rlj))-
// 		     0.5*rli/Mi2*log((m2+M2-Mi2-rli)/
// 				     (m2+M2-Mi2+rli)));
    K = 1./(Mi2-Mj2)*(1.+Complex(m2*I+M2*J-Mj2*I2));
  }
  // leading term for small m
  else {
    I  = 1./(Mj2-Mi2)*(-Li2(double(Mj2/(Mj2-M2)))+Li2(double(Mi2/(Mi2-M2)))
		       -2.*log(m/M)*log((M2-Mj2)/(M2-Mi2))
		       +0.5*sqr(log((M2-Mj2)/M2))-0.5*sqr(log((M2-Mi2)/M2)));
    J  = 1./(Mj2-Mi2)*(Li2(double(Mj2/(Mj2-M2)))-Li2(double(Mi2/(Mi2-M2)))
		       -0.5*sqr(log((M2-Mi2)/M2))+0.5*sqr(log((M2-Mj2)/M2)));
    I2 = 1./(Mj2-Mi2)*log((M2-Mj2)/(M2-Mi2))
      +M2/(Mj2-Mi2)/sqr(Mi*Mj)*(Mj2*log((M2-Mi2)/M2)-Mi2*log((M2-Mj2)/M2));
    K  = 1./(Mi2-Mj2)*(1.+Complex(M2*J)-Complex(Mj2*I2));
  }
}
