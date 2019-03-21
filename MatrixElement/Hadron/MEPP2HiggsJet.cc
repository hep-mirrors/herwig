// -*- C++ -*-
//
// MEPP2HiggsJet.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2HiggsJet class.
//

#include "MEPP2HiggsJet.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig/MatrixElement/HardVertex.h"

using namespace Herwig;

const Complex MEPP2HiggsJet::_epsi = Complex(0.,-1.e-20);

IBPtr MEPP2HiggsJet::clone() const {
  return new_ptr(*this);
}

IBPtr MEPP2HiggsJet::fullclone() const {
  return new_ptr(*this);
}

unsigned int MEPP2HiggsJet::orderInAlphaS() const {
  return 3;
}

unsigned int MEPP2HiggsJet::orderInAlphaEW() const {
  return 1;
}

void MEPP2HiggsJet::persistentOutput(PersistentOStream & os) const {
  os << _shapeopt << _maxflavour << _process << _minloop << _maxloop << _massopt
     << ounit(_mh,GeV) << ounit(_wh,GeV) << _hmass;
}

void MEPP2HiggsJet::persistentInput(PersistentIStream & is, int) {
  is >> _shapeopt >> _maxflavour >> _process >> _minloop >> _maxloop >> _massopt
     >> iunit(_mh,GeV) >> iunit(_wh,GeV) >> _hmass;
}

ClassDescription<MEPP2HiggsJet> MEPP2HiggsJet::initMEPP2HiggsJet;
// Definition of the static class description member.

void MEPP2HiggsJet::Init() {

  static ClassDocumentation<MEPP2HiggsJet> documentation
    ("The MEPP2HiggsJet class implements the matrix elements for"
     " Higgs+Jet production in hadron-hadron collisions.",
     "The theoretical calculations of \\cite{Baur:1989cm} and \\cite{Ellis:1987xu}"
     " were used for the Higgs+jet matrix element in hadron-hadron collisions.",
     "\\bibitem{Baur:1989cm} U.~Baur and E.~W.~N.~Glover,"
     "Nucl.\\ Phys.\\ B {\\bf 339} (1990) 38.\n"
     "\\bibitem{Ellis:1987xu} R.~K.~Ellis, I.~Hinchliffe, M.~Soldate and "
     "J.~J.~van der Bij, Nucl.\\ Phys.\\ B {\\bf 297} (1988) 221.");

  static Parameter<MEPP2HiggsJet,unsigned int> interfaceMaximumFlavour
    ("MaximumFlavour",
     "The maximum flavour of the quarks in the process",
     &MEPP2HiggsJet::_maxflavour, 5, 1, 5,
     false, false, Interface::limited);

  static Switch<MEPP2HiggsJet,unsigned int> interfaceShapeOption
    ("ShapeScheme",
     "Option for the treatment of the Higgs resonance shape",
     &MEPP2HiggsJet::_shapeopt, 1, false, false);
  static SwitchOption interfaceStandardShapeFixed
    (interfaceShapeOption,
     "FixedBreitWigner",
     "Breit-Wigner s-channel resonanse",
     1);
  static SwitchOption interfaceStandardShapeRunning
    (interfaceShapeOption,
     "MassGenerator",
     "Use the mass generator to give the shape",
     2);

  static Switch<MEPP2HiggsJet,unsigned int> interfaceProcess
    ("Process",
     "Which subprocesses to include",
     &MEPP2HiggsJet::_process, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Include all subprocesses",
     0);
  static SwitchOption interfaceProcess1
    (interfaceProcess,
     "qqbar",
     "Only include the incoming q qbar subprocess",
     1);
  static SwitchOption interfaceProcessqg
    (interfaceProcess,
     "qg",
     "Only include the incoming qg subprocess",
     2);
  static SwitchOption interfaceProcessqbarg
    (interfaceProcess,
     "qbarg",
     "Only include the incoming qbar g subprocess",
     3);
  static SwitchOption interfaceProcessgg
    (interfaceProcess,
     "gg",
     "Only include the incoming gg subprocess",
     4);

  static Parameter<MEPP2HiggsJet,int> interfaceMinimumInLoop
    ("MinimumInLoop",
     "The minimum flavour of the quarks to include in the loops",
     &MEPP2HiggsJet::_minloop, 6, 4, 6,
     false, false, Interface::limited);

  static Parameter<MEPP2HiggsJet,int> interfaceMaximumInLoop
    ("MaximumInLoop",
     "The maximum flavour of the quarks to include in the loops",
     &MEPP2HiggsJet::_maxloop, 6, 4, 6,
     false, false, Interface::limited);

  static Switch<MEPP2HiggsJet,unsigned int> interfaceMassOption
    ("MassOption",
     "Option for the treatment of the masses in the loop diagrams",
     &MEPP2HiggsJet::_massopt, 0, false, false);
  static SwitchOption interfaceMassOptionFull
    (interfaceMassOption,
     "Full",
     "Include the full mass dependence",
     0);
  static SwitchOption interfaceMassOptionLarge
    (interfaceMassOption,
     "Large",
     "Use the heavy mass limit",
     1);

}

bool MEPP2HiggsJet::generateKinematics(const double * r) {
  Energy ptmin = max(lastCuts().minKT(mePartonData()[2]),
   		     lastCuts().minKT(mePartonData()[3]));
  Energy e = sqrt(sHat())/2.0;
  // generate the mass of the higgs boson
  Energy2 mhmax2 = sHat()-4.*ptmin*e;
  Energy2 mhmin2 =ZERO;
  if(mhmax2<=mhmin2) return false;
  double rhomin = atan2((mhmin2-sqr(_mh)), _mh*_wh);
  double rhomax = atan2((mhmax2-sqr(_mh)), _mh*_wh);
  Energy mh = sqrt(_mh*_wh*tan(rhomin+r[1]*(rhomax-rhomin))+sqr(_mh));
  // assign masses
  if(mePartonData()[2]->id()!=ParticleID::h0) {
    meMomenta()[2].setMass(ZERO);
    meMomenta()[3].setMass(mh);
  }
  else {
    meMomenta()[3].setMass(ZERO);
    meMomenta()[2].setMass(mh);
  }

  Energy q = ZERO;
  try {
    q = SimplePhaseSpace::
      getMagnitude(sHat(), meMomenta()[2].mass(), meMomenta()[3].mass());
  } 
  catch ( ImpossibleKinematics & e ) {
    return false;
  }
		    
  Energy2 m22  = meMomenta()[2].mass2();
  Energy2 m32  = meMomenta()[3].mass2();
  Energy2 e0e2 = 2.0*e*sqrt(sqr(q) + m22);
  Energy2 e1e2 = 2.0*e*sqrt(sqr(q) + m22);
  Energy2 e0e3 = 2.0*e*sqrt(sqr(q) + m32);
  Energy2 e1e3 = 2.0*e*sqrt(sqr(q) + m32);
  Energy2 pq   = 2.0*e*q;
  double ctmin = -1.0,ctmax = 1.0;
  Energy2 thmin = lastCuts().minTij(mePartonData()[0], mePartonData()[2]);
  if ( thmin > ZERO ) ctmax = min(ctmax, (e0e2 - m22 - thmin)/pq);

  thmin = lastCuts().minTij(mePartonData()[1], mePartonData()[2]);
  if ( thmin > ZERO ) ctmin = max(ctmin, (thmin + m22 - e1e2)/pq);

  thmin = lastCuts().minTij(mePartonData()[1], mePartonData()[3]);
  if ( thmin > ZERO ) ctmax = min(ctmax, (e1e3 - m32 - thmin)/pq);

  thmin = lastCuts().minTij(mePartonData()[0], mePartonData()[3]);
  if ( thmin > ZERO ) ctmin = max(ctmin, (thmin + m32 - e0e3)/pq);
  if ( ptmin > ZERO ) {
    double ctm = 1.0 - sqr(ptmin/q);
    if ( ctm <= 0.0 ) return false;
    ctmin = max(ctmin, -sqrt(ctm));
    ctmax = min(ctmax, sqrt(ctm));
  }

  if ( ctmin >= ctmax ) return false;
    
  double cth = getCosTheta(ctmin, ctmax, r);

  Energy pt = q*sqrt(1.0-sqr(cth));
  phi(rnd(2.0*Constants::pi));
  meMomenta()[2].setX(pt*sin(phi()));
  meMomenta()[2].setY(pt*cos(phi()));
  meMomenta()[2].setZ(q*cth);

  meMomenta()[3].setX(-pt*sin(phi()));
  meMomenta()[3].setY(-pt*cos(phi()));
  meMomenta()[3].setZ(-q*cth);

  meMomenta()[2].rescaleEnergy();
  meMomenta()[3].rescaleEnergy();

  vector<LorentzMomentum> out(2);
  out[0] = meMomenta()[2];
  out[1] = meMomenta()[3];
  tcPDVector tout(2);
  tout[0] = mePartonData()[2];
  tout[1] = mePartonData()[3];
  if ( !lastCuts().passCuts(tout, out, mePartonData()[0], mePartonData()[1]) )
    return false;

  tHat(pq*cth + m22 - e0e2);
  uHat(m22 + m32 - sHat() - tHat());
  // main piece
  jacobian((pq/sHat())*Constants::pi*jacobian());
  // mass piece
  jacobian((rhomax-rhomin)*jacobian());
  return true;
}

void MEPP2HiggsJet::getDiagrams() const {
  tcPDPtr h0=getParticleData(ParticleID::h0);
  tcPDPtr g =getParticleData(ParticleID::g);
  tcPDPtr q[6],qb[6];
  for(int ix=0;ix<int(_maxflavour);++ix) {
    q [ix]=getParticleData( ix+1);
    qb[ix]=getParticleData(-ix-1);
  }
  // q qbar -> H g
  if(_process==0||_process==1)
    {for(unsigned int ix=0;ix<_maxflavour;++ix)
	{add(new_ptr((Tree2toNDiagram(2), q[ix], qb[ix], 1, g , 3, h0, 3, g, -1)));}}

  // q g ->  H g
  if(_process==0||_process==2)
    {for(unsigned int ix=0;ix<_maxflavour;++ix)
	{add(new_ptr((Tree2toNDiagram(3), q[ix], g, g, 2, h0, 1, q[ix], -2)));}}

  // qbar g -> H qbar 
  if(_process==0||_process==3)
    {for(unsigned int ix=0;ix<_maxflavour;++ix)
	{add(new_ptr((Tree2toNDiagram(3), qb[ix], g, g, 2, h0, 1, qb[ix], -3)));}}

  // g g -> H g
  if(_process==0||_process==4)
    {
      // t channel
      add(new_ptr((Tree2toNDiagram(3), g, g, g, 1, h0, 2, g, -4)));
      // u channel
      add(new_ptr((Tree2toNDiagram(3), g, g, g, 2, h0, 1, g, -5)));
      // s channel
      add(new_ptr((Tree2toNDiagram(2), g, g, 1, g , 3, h0, 3, g, -6)));
    }
}

Energy2 MEPP2HiggsJet::scale() const {
  return meMomenta()[2].perp2()+ meMomenta()[2].m2();
}

double MEPP2HiggsJet::me2() const {
  useMe();
  double output(0.);
  // g g to H g
  if(mePartonData()[0]->id()==ParticleID::g&&mePartonData()[1]->id()==ParticleID::g) {
    // order of the particles
    unsigned int ih(2),ig(3);
    if(mePartonData()[3]->id()==ParticleID::h0){ig=2;ih=3;}
    VectorWaveFunction glin1(meMomenta()[ 0],mePartonData()[ 0],incoming);
    VectorWaveFunction glin2(meMomenta()[ 1],mePartonData()[ 1],incoming);
    ScalarWaveFunction  hout(meMomenta()[ih],mePartonData()[ih],outgoing);
    VectorWaveFunction glout(meMomenta()[ig],mePartonData()[ig],outgoing);
    vector<VectorWaveFunction> g1,g2,g4;
    for(unsigned int ix=0;ix<2;++ix) {
      glin1.reset(2*ix);g1.push_back(glin1);
      glin2.reset(2*ix);g2.push_back(glin2);
      glout.reset(2*ix);g4.push_back(glout);
    }
    // calculate the matrix element
    output = ggME(g1,g2,hout,g4,false); 
  }
  // qg -> H q
  else if(mePartonData()[0]->id()>0&&mePartonData()[1]->id()==ParticleID::g) {
    // order of the particles
    unsigned int iq(0),iqb(3),ih(2),ig(1);
    if(mePartonData()[0]->id()==ParticleID::g){iq=1;ig=0;}
    if(mePartonData()[3]->id()==ParticleID::h0){iqb=2;ih=3;}
    // calculate the spinors and polarization vectors
    vector<SpinorWaveFunction> fin;
    vector<SpinorBarWaveFunction>  fout;
    vector<VectorWaveFunction> gin;
    SpinorWaveFunction    qin (meMomenta()[iq ],mePartonData()[iq ],incoming);
    VectorWaveFunction    glin(meMomenta()[ig ],mePartonData()[ig ],incoming);
    ScalarWaveFunction    hout(meMomenta()[ih ],mePartonData()[ih ],outgoing);
    SpinorBarWaveFunction qout(meMomenta()[iqb],mePartonData()[iqb],outgoing);
    for(unsigned int ix=0;ix<2;++ix) {
      qin.reset(ix)   ; fin.push_back( qin);
      qout.reset(ix)  ;fout.push_back(qout);
      glin.reset(2*ix); gin.push_back(glin);
    }
    // calculate the matrix element
    output = qgME(fin,gin,hout,fout,false); 
  }
  // qbar g -> H q
  else if(mePartonData()[0]->id()<0&&mePartonData()[1]->id()==ParticleID::g) {
    // order of the particles
    unsigned int iq(0),iqb(3),ih(2),ig(1);
    if(mePartonData()[0]->id()==ParticleID::g){iq=1;ig=0;}
    if(mePartonData()[3]->id()==ParticleID::h0){iqb=2;ih=3;}
    // calculate the spinors and polarization vectors
    vector<SpinorBarWaveFunction> fin;
    vector<SpinorWaveFunction>  fout;
    vector<VectorWaveFunction> gin;
    SpinorBarWaveFunction qin (meMomenta()[iq ],mePartonData()[iq ],incoming);
    VectorWaveFunction    glin(meMomenta()[ig ],mePartonData()[ig ],incoming);
    ScalarWaveFunction    hout(meMomenta()[ih ],mePartonData()[ih ],outgoing);
    SpinorWaveFunction    qout(meMomenta()[iqb],mePartonData()[iqb],outgoing);
    for(unsigned int ix=0;ix<2;++ix) {
      qin.reset(ix)   ; fin.push_back( qin);
      qout.reset(ix)  ;fout.push_back(qout);
      glin.reset(2*ix); gin.push_back(glin);
    }
    // calculate the matrix element
    output = qbargME(fin,gin,hout,fout,false); 
  }
  // q qbar to H g
  else if(mePartonData()[0]->id()==-mePartonData()[1]->id()) {
    // order of the particles
    unsigned int iq(0),iqb(1),ih(2),ig(3);
    if(mePartonData()[0]->id()<0){iq=1;iqb=0;}
    if(mePartonData()[2]->id()==ParticleID::g){ig=2;ih=3;}
    // calculate the spinors and polarization vectors
    vector<SpinorWaveFunction> fin;
    vector<SpinorBarWaveFunction>  ain;
    vector<VectorWaveFunction> gout;
    SpinorWaveFunction    qin (meMomenta()[iq ],mePartonData()[iq ],incoming);
    SpinorBarWaveFunction qbin(meMomenta()[iqb],mePartonData()[iqb],incoming);
    ScalarWaveFunction    hout(meMomenta()[ih ],mePartonData()[ih ],outgoing);
    VectorWaveFunction   glout(meMomenta()[ig ],mePartonData()[ig ],outgoing);
    for(unsigned int ix=0;ix<2;++ix) {
      qin.reset(ix)    ; fin.push_back(  qin);
      qbin.reset(ix)   ; ain.push_back( qbin);
      glout.reset(2*ix);gout.push_back(glout);
    }
    // calculate the matrix element
    output = qqbarME(fin,ain,hout,gout,false); 
  }
  else
    throw Exception() << "Unknown subprocess in MEPP2HiggsJet::me2()" 
		      << Exception::runerror;
  // return the answer
  return output;
}

double MEPP2HiggsJet::qqbarME(vector<SpinorWaveFunction>    & fin,
			      vector<SpinorBarWaveFunction> & ain,
			      ScalarWaveFunction & hout,
			      vector<VectorWaveFunction>    & gout,
			      bool calc) const {
  // the particles should be in the order
  // for the incoming 
  // 0 incoming fermion     (u    spinor)
  // 1 incoming antifermion (vbar spinor)
  // for the outgoing
  // 0 outgoing higgs       
  // 1 outgoing gluon
  // me to be returned
  ProductionMatrixElement newme(PDT::Spin1Half,PDT::Spin1Half,
				PDT::Spin0,PDT::Spin1);
  // get the kinematic invariants
  Energy2 s(sHat()),u(uHat()),t(tHat()),mh2(hout.m2()),et(scale());
  // calculate the loop function
  complex<Energy2> A5 = Energy2();
  for ( int ix=_minloop; ix<=_maxloop; ++ix ) {
    // full mass dependance
    if(_massopt==0) {
      Energy2 mf2=sqr(getParticleData(ix)->mass());
      A5+= mf2*(4.+4.*double(s/(u+t))*(W1(s,mf2)-W1(mh2,mf2))
		+(1.-4.*double(mf2/(u+t)))*(W2(s,mf2)-W2(mh2,mf2)));	
    }
    // infinite mass limit
    else {
      A5+=2.*(s-mh2)/3.;
    }
  }
  // multiply by the rest of the form factors
  using Constants::pi;
  double g(sqrt(4.*pi*SM().alphaEM(mh2)/SM().sin2ThetaW()));
  double gs(sqrt(4.*pi*SM().alphaS(et)));
  Energy mw(getParticleData(ParticleID::Wplus)->mass());
  complex<InvEnergy> A5c = A5 * Complex(0.,1.)*g*sqr(gs)*gs/(32.*s*sqr(pi)*mw);
  // compute the matrix element
  LorentzPolarizationVectorE fcurrent;
  complex<Energy2> fdotp;
  complex<Energy> epsdot[2];
  Complex diag;
  Lorentz5Momentum ps(fin[0].momentum()+ain[0].momentum());
  ps.rescaleMass();
  for(unsigned int ix=0;ix<2;++ix){epsdot[ix]=gout[ix].wave().dot(ps);}
  Energy2 denom(-ps*gout[0].momentum());
  LorentzSpinorBar<double> atemp;
  double output(0.);
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      // compute the fermion current
      atemp=ain[ihel2].wave();
      fcurrent=UnitRemoval::E*fin[ihel1].wave().vectorCurrent(atemp);
      fdotp = -(fcurrent.dot(gout[0].momentum()));
      for(unsigned int ghel=0;ghel<2;++ghel) {
	// calculate the matrix element
	diag = Complex(A5c*(fcurrent.dot(gout[ghel].wave())
			    -fdotp*epsdot[ghel]/denom));
	// calculate the matrix element
	output+=real(diag*conj(diag));
	if(calc) newme(ihel1,ihel2,0,2*ghel)=diag;
      }
    }
  }
  // test with glover form 
  // final colour/spin factors
  if(calc) _me.reset(newme);
  return output/9.;
}

double MEPP2HiggsJet::qgME(vector<SpinorWaveFunction> & fin,
			   vector<VectorWaveFunction> & gin,
			   ScalarWaveFunction & hout, 
			   vector<SpinorBarWaveFunction> & fout,bool calc) const {
  // the particles should be in the order
  // for the incoming 
  // 0 incoming fermion     (u    spinor)
  // 1 incoming gluon
  // for the outgoing
  // 0 outgoing higgs       
  // 1 outgoing fermion     (ubar spinor)
  // me to be returned
  ProductionMatrixElement newme(PDT::Spin1Half,PDT::Spin1,
				PDT::Spin0,PDT::Spin1Half);
  // get the kinematic invariants
  Energy2 s(sHat()),u(uHat()),t(tHat()),mh2(hout.m2()),et(scale());
  // calculate the loop function
  complex<Energy2> A5 = Energy2();
  for(int ix=_minloop;ix<=_maxloop;++ix) {
      if(_massopt==0) {
	Energy2 mf2=sqr(getParticleData(ix)->mass());
	A5+= mf2*(4.+4.*double(u/(s+t))*(W1(u,mf2)-W1(mh2,mf2))
		  +(1.-4.*double(mf2/(s+t)))*(W2(u,mf2)-W2(mh2,mf2)));
      }
      else {
	A5+=2.*(u-mh2)/3.;
      }
  }
  // multiply by the rest of the form factors
  using Constants::pi;
  double g(sqrt(4.*pi*SM().alphaEM(mh2)/SM().sin2ThetaW()));
  double gs(sqrt(4.*pi*SM().alphaS(et)));
  Energy mw(getParticleData(ParticleID::Wplus)->mass());
  complex<InvEnergy> A5c =A5*Complex(0.,1.)*g*sqr(gs)*gs/(32.*u*sqr(pi)*mw);
  // compute the matrix element
  LorentzPolarizationVectorE fcurrent;
  complex<Energy2> fdotp;
  complex<Energy> epsdot[2];
  Complex diag;
  Lorentz5Momentum pu(fin[0].momentum()+fout[0].momentum());
  pu.rescaleMass();
  for(unsigned int ix=0;ix<2;++ix){epsdot[ix]=gin[ix].wave().dot(pu);}
  Energy2 denom(pu*gin[0].momentum());
  LorentzSpinorBar<double> atemp;
  double output(0.);
  for(unsigned int ihel=0;ihel<2;++ihel) {
    for(unsigned int ohel=0;ohel<2;++ohel) {
      // compute the fermion current
      atemp=fout[ohel].wave();
      fcurrent=UnitRemoval::E*fin[ihel].wave().vectorCurrent(atemp);
      fdotp=fcurrent.dot(gin[0].momentum());
      for(unsigned int ghel=0;ghel<2;++ghel) {
	// calculate the matrix element
	diag = Complex(A5c*(fcurrent.dot(gin[ghel].wave())-fdotp*epsdot[ghel]/denom));
	// calculate the matrix element
	output+=real(diag*conj(diag));
	if(calc) newme(ihel,2*ghel,0,ohel)=diag;
      }
    }
  }
  // final colour/spin factors
  if(calc) _me.reset(newme);
  return output/24.;
}

double MEPP2HiggsJet::qbargME(vector<SpinorBarWaveFunction> & fin,
			      vector<VectorWaveFunction> & gin,
			      ScalarWaveFunction & hout,
			      vector<SpinorWaveFunction> & fout,bool calc) const {
  // the particles should be in the order
  // for the incoming 
  // 0 incoming antifermion (vbar spinor)       
  // 1 incoming gluon
  // for the outgoing
  // 0 outgoing higgs
  // 1 outgoing antifermion (v    spinor)
  // me to be returned
  ProductionMatrixElement newme(PDT::Spin1Half,PDT::Spin1,
				PDT::Spin0,PDT::Spin1Half);
  // get the kinematic invariants
  Energy2 s(sHat()),u(uHat()),t(tHat()),mh2(hout.m2()),et(scale());
  // calculate the loop function
  complex<Energy2> A5 = Energy2();
  for(int ix=_minloop;ix<=_maxloop;++ix) {
    if(_massopt==0) {
      Energy2 mf2=sqr(getParticleData(ix)->mass());
      A5+= mf2*(4.+4.*double(u/(s+t))*(W1(u,mf2)-W1(mh2,mf2))
		+(1.-4.*double(mf2/(s+t)))*(W2(u,mf2)-W2(mh2,mf2)));
    }
    else { 
      A5+=2.*(u-mh2)/3.;
    }
  }
  // multiply by the rest of the form factors
  using Constants::pi;
  double g(sqrt(4.*pi*SM().alphaEM(mh2)/SM().sin2ThetaW()));
  double gs(sqrt(4.*pi*SM().alphaS(et)));
  Energy mw(getParticleData(ParticleID::Wplus)->mass());
  complex<InvEnergy> A5c = A5*Complex(0.,1.)*g*sqr(gs)*gs/(32.*u*sqr(pi)*mw);
  // compute the matrix element
  LorentzPolarizationVectorE fcurrent;
  complex<Energy2> fdotp;
  complex<Energy> epsdot[2];
  Complex diag;
  Lorentz5Momentum pu(fin[0].momentum()+fout[0].momentum());
  pu.rescaleMass();
  for(unsigned int ix=0;ix<2;++ix){epsdot[ix]=gin[ix].wave().dot(pu);}
  Energy2 denom(pu*gin[0].momentum());
  LorentzSpinorBar<double> atemp;
  double output(0.);
  for(unsigned int ihel=0;ihel<2;++ihel) {
    for(unsigned int ohel=0;ohel<2;++ohel) {
      // compute the fermion current
      atemp=fin[ihel].wave();
      fcurrent=UnitRemoval::E*fout[ohel].wave().vectorCurrent(atemp);
      fdotp=fcurrent.dot(gin[0].momentum());
      for(unsigned int ghel=0;ghel<2;++ghel) {
	// calculate the matrix element
	diag = Complex(A5c*(fcurrent.dot(gin[ghel].wave())-fdotp*epsdot[ghel]/denom));
	// calculate the matrix element
	output+=real(diag*conj(diag));
	if(calc) newme(ihel,2*ghel,0,ohel)=diag;
      }
    }
  }
  // final colour/spin factors
  if(calc) _me.reset(newme);
  return output/24.;
}

Selector<MEBase::DiagramIndex>
MEPP2HiggsJet::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i )
    {
      if(abs(diags[i]->id())<4) sel.insert(1.0, i);
      else sel.insert(_diagwgt[abs(diags[i]->id())-4], i);
    }
  return sel;
}

Selector<const ColourLines *>
MEPP2HiggsJet::colourGeometries(tcDiagPtr diag) const {
  // colour lines for q qbar -> h0 g
  static const ColourLines cqqbar("1 3 5,-2 -3 -5");
  // colour lines for q g -> h0 q
  static const ColourLines cqg("1 2 -3, 3 -2 5");
  // colour lines for qbar q -> h0 qbar
  static const ColourLines cqbarg("-1 -2 3, -3 2 -5");
  // colour lines for g g -> h0 g
  static const ColourLines cgg[6]={ColourLines("1 2 5, -3 -5, 3 -2 -1"),
			     ColourLines("-1 -2 -5, 3 5, -3 2 1"),
			     ColourLines("1 5, -1 -2 3, -3 2 -5"),
			     ColourLines("-1 -5, 1 2 -3, 3 -2 5"),
			     ColourLines("1 3 5, -5 -3 -2, 2 -1"),
			     ColourLines("-1 -3 -5, 5 3 2 ,-2 1")};
  // select the colour flow
  Selector<const ColourLines *> sel;
  if ( diag->id() == -1)      sel.insert(1.0, &cqqbar);
  else if ( diag->id() == -2) sel.insert(1.0, &cqg);
  else if ( diag->id() == -3) sel.insert(1.0, &cqbarg);
  else
    {
      sel.insert(0.5, &cgg[2*(abs(diag->id())-4)  ]);
      sel.insert(0.5, &cgg[2*(abs(diag->id())-4)+1]);
    }
  // return the answer
  return sel;
}

double MEPP2HiggsJet::ggME(vector<VectorWaveFunction> g1, vector<VectorWaveFunction> g2,
			   ScalarWaveFunction & hout,     vector<VectorWaveFunction> g4,
			   bool calc) const {
  // the particles should be in the order
  // for the incoming 
  // 0 first  incoming gluon
  // 1 second incoming gluon
  // for the outgoing
  // 0 outgoing higgs
  // 1 outgoing gluon
  // me to be returned
  ProductionMatrixElement newme(PDT::Spin1,PDT::Spin1,
				PDT::Spin0,PDT::Spin1);
   // get the kinematic invariants
   Energy2 s(sHat()),u(uHat()),t(tHat()),mh2(hout.m2()),et(scale());
   // calculate the loop functions
   Complex A4stu(0.),A2stu(0.),A2tsu(0.),A2ust(0.);
   Complex A5s(0.),A5t(0.),A5u(0.);
   for(int ix=_minloop;ix<=_maxloop;++ix) {
     Energy2 mf2=sqr(getParticleData(ix)->mass());
     // loop functions
     if(_massopt==0) {
       A4stu+=A4(s,t,u,mf2);
       A2stu+=A2(s,t,u,mf2);
       A2tsu+=A2(u,s,t,mf2);
       A2ust+=A2(t,s,u,mf2);
       A5s+= double(mf2/s)*(4.+4.*double(s/(u+t))*(W1(s,mf2)-W1(mh2,mf2))
			    +(1.-4.*double(mf2/(u+t)))*(W2(s,mf2)-W2(mh2,mf2)));
       A5t+= double(mf2/t)*(4.+4.*double(t/(s+u))*(W1(t,mf2)-W1(mh2,mf2))
			    +(1.-4.*double(mf2/(s+u)))*(W2(t,mf2)-W2(mh2,mf2)));
       A5u+= double(mf2/u)*(4.+4.*double(u/(s+t))*(W1(u,mf2)-W1(mh2,mf2))
			    +(1.-4.*double(mf2/(s+t)))*(W2(u,mf2)-W2(mh2,mf2)));
     }
     else {
       A4stu=-1./3.;
       A2stu=-sqr(s/mh2)/3.;
       A2tsu=-sqr(t/mh2)/3.;
       A2ust=-sqr(u/mh2)/3.;
       A5s+=2.*(s-mh2)/3./s;
       A5t+=2.*(t-mh2)/3./t;
       A5u+=2.*(u-mh2)/3./u;
     }
   }
   Complex A3stu=0.5*(A2stu+A2ust+A2tsu-A4stu);
   // compute the dot products for the matrix element
   complex<InvEnergy> eps[3][4][2];
   Energy2 pdot[4][4];
   pdot[0][0]=ZERO;
   pdot[0][1]= g1[0].momentum()*g2[0].momentum();
   pdot[0][2]=-1.*g1[0].momentum()*g4[0].momentum();
   pdot[0][3]=-1.*g1[0].momentum()*hout.momentum();
   pdot[1][0]= pdot[0][1];
   pdot[1][1]= ZERO;
   pdot[1][2]=-1.*g2[0].momentum()*g4[0].momentum();
   pdot[1][3]=-1.*g2[0].momentum()*hout.momentum();
   pdot[2][0]= pdot[0][2];
   pdot[2][1]= pdot[1][2];
   pdot[2][2]= ZERO;
   pdot[2][3]= g4[0].momentum()*hout.momentum();
   pdot[3][0]=pdot[0][3];
   pdot[3][1]=pdot[1][3];
   pdot[3][2]=pdot[2][3];
   pdot[3][3]=mh2;
   for(unsigned int ix=0;ix<2;++ix)
     {
       eps[0][0][ix]=InvEnergy();
       eps[0][1][ix]=g1[ix].wave().dot(g2[0].momentum())/pdot[0][1];
       eps[0][2][ix]=-1.*g1[ix].wave().dot(g4[0].momentum())/pdot[0][2];
       eps[0][3][ix]=-1.*g1[ix].wave().dot(hout.momentum())/ pdot[0][3];
       eps[1][0][ix]=g2[ix].wave().dot(g1[0].momentum())/    pdot[1][0];
       eps[1][1][ix]=InvEnergy();
       eps[1][2][ix]=-1.*g2[ix].wave().dot(g4[0].momentum())/pdot[1][2];
       eps[1][3][ix]=-1.*g2[ix].wave().dot(hout.momentum())/ pdot[1][3];
       eps[2][0][ix]=g4[ix].wave().dot(g1[0].momentum())/    pdot[2][0];
       eps[2][1][ix]=g4[ix].wave().dot(g2[0].momentum())/    pdot[2][1];
       eps[2][2][ix]=InvEnergy();
       eps[2][3][ix]=-1.*g4[ix].wave().dot(hout.momentum())/     pdot[2][3];
     }
   // prefactors
   using Constants::pi;
   double g(sqrt(4.*pi*SM().alphaEM(mh2)/SM().sin2ThetaW()));
   double gs(sqrt(4.*pi*SM().alphaS(et)));
   Energy mw(getParticleData(ParticleID::Wplus)->mass());
   Energy3 pre=g*sqr(mh2)*gs*sqr(gs)/(32.*sqr(pi)*mw);
   // compute the matrix element
   double output(0.);
   Complex diag[4],wdot[3][3];
   _diagwgt[0]=0.;
   _diagwgt[1]=0.;
   _diagwgt[2]=0.;
   for(unsigned int ihel1=0;ihel1<2;++ihel1) {
     for(unsigned int ihel2=0;ihel2<2;++ihel2) {
       for(unsigned int ohel=0;ohel<2;++ohel) {
	 wdot[0][1]=g1[ihel1].wave().dot(g2[ihel2].wave());
	 wdot[0][2]=g1[ihel1].wave().dot(g4[ohel ].wave());
	 wdot[1][0]=wdot[0][1];
	 wdot[1][2]=g2[ihel2].wave().dot(g4[ohel ].wave());
	 wdot[2][0]=wdot[0][2];
	 wdot[2][1]=wdot[1][2];
	 // last piece
	 diag[3]= Complex(pre*A3stu*(eps[0][2][ihel1]*eps[1][0][ihel2]*eps[2][1][ohel]-
				     eps[0][1][ihel1]*eps[1][2][ihel2]*eps[2][0][ohel]+
				     (eps[2][0][ohel ]-eps[2][1][ohel ])*wdot[0][1]/pdot[0][1]+
				     (eps[1][2][ihel2]-eps[1][0][ihel2])*wdot[0][2]/pdot[0][2]+
				     (eps[0][1][ihel1]-eps[0][2][ihel1])*wdot[1][2]/pdot[1][2]));
	 // first piece
	 diag[3] += Complex(pre*(+A2stu*(eps[0][1][ihel1]*eps[1][0][ihel2]-wdot[0][1]/pdot[0][1])*
				 (eps[2][0][ohel ]-eps[2][1][ohel ])
				 +A2ust*(eps[0][2][ihel1]*eps[2][0][ohel ]-wdot[0][2]/pdot[0][2])*
				 (eps[1][2][ihel2]-eps[1][0][ihel2])
				 +A2tsu*(eps[1][2][ihel2]*eps[2][1][ohel ]-wdot[1][2]/pdot[1][2])*
				 (eps[0][1][ihel1]-eps[0][2][ihel1])
				 ));
	 output+=real(diag[3]*conj(diag[3]));
	 // matrix element if needed
	 if(calc) newme(2*ihel1,2*ihel2,0,2*ohel)=diag[3];
	 // different diagrams 
	 diag[0] = Complex(A5t*UnitRemoval::InvE*(-eps[0][3][ihel1]*
						  (-2.*eps[2][1][ohel ]*eps[1][0][ihel2]*pdot[2][1]*pdot[1][0]
						   -2.*eps[1][2][ihel2]*eps[2][0][ohel ]*pdot[1][2]*pdot[2][0]
						   +wdot[1][2]*(pdot[0][1]+pdot[0][2]))
						  -2.*eps[2][1][ohel ]*pdot[2][1]*wdot[0][1]
						  -2.*eps[1][2][ihel2]*pdot[1][2]*wdot[0][2]
						  +wdot[1][2]*(eps[0][1][ihel1]*pdot[0][1]+
							       eps[0][2][ihel1]*pdot[0][2])));
	 diag[1] = Complex(A5u*UnitRemoval::InvE*(-eps[1][3][ihel2]*
						  (+2.*eps[0][1][ihel1]*eps[2][0][ohel ]*pdot[0][1]*pdot[2][0]
						   +2.*eps[0][2][ihel1]*eps[2][1][ohel ]*pdot[0][2]*pdot[2][1]
						   -wdot[0][2]*(pdot[1][0]+pdot[1][2]))
						  +2.*eps[2][0][ohel ]*pdot[2][0]*wdot[0][1]
						  +2.*eps[0][2][ihel1]*pdot[0][2]*wdot[2][1]
						  -wdot[0][2]*(eps[1][0][ihel2]*pdot[1][0]+
							       eps[1][2][ihel2]*pdot[1][2])));
	 diag[2] = Complex(A5s*UnitRemoval::InvE*(-eps[2][3][ohel ]*
						  (+2.*eps[0][1][ihel1]*eps[1][2][ihel2]*pdot[0][1]*pdot[1][2]
						   -2.*eps[1][0][ihel2]*eps[0][2][ihel1]*pdot[1][0]*pdot[1][3]
						   +wdot[0][1]*(pdot[2][0]-pdot[2][1]))
						  +2.*eps[0][1][ihel1]*pdot[0][1]*wdot[1][2]
						  -2.*eps[1][0][ihel2]*pdot[1][0]*wdot[0][2]
						  +wdot[0][1]*(eps[2][0][ohel]*pdot[2][0]-
							       eps[2][1][ohel]*pdot[2][1])));
	 _diagwgt[0]+=real(diag[0]*conj(diag[0]));
	 _diagwgt[1]+=real(diag[1]*conj(diag[1]));
	 _diagwgt[2]+=real(diag[2]*conj(diag[2]));
       }
     }
   }
   // final colour and spin factors
   if(calc){_me.reset(newme);}
   return 3.*output/32.;
}

void MEPP2HiggsJet::constructVertex(tSubProPtr sub)
{
  // extract the particles in the hard process
  ParticleVector hard;
  hard.push_back(sub->incoming().first);hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);hard.push_back(sub->outgoing()[1]);
  // ensure correct order or particles
  if((hard[0]->id()==ParticleID::g&&hard[1]->id()!=ParticleID::g)|| 
     (hard[0]->id()<0&&hard[1]->id()<6)) swap(hard[0],hard[1]);
  if(hard[2]->id()!=ParticleID::h0) swap(hard[2],hard[3]);
  // different processes
  // g g to H g
  if(hard[0]->id()==ParticleID::g) {
    vector<VectorWaveFunction> g1,g2,g4;
    VectorWaveFunction(g1,hard[0],incoming,false,true,true);
    VectorWaveFunction(g2,hard[1],incoming,false,true,true);
    VectorWaveFunction(g4,hard[3],outgoing,true ,true,true);
    ScalarWaveFunction hout(hard[2],outgoing,true);
    g1[1]=g1[2];g2[1]=g2[2];g4[1]=g4[2];
    ggME(g1,g2,hout,g4,true);
  }
  // qg -> H q
  else if(hard[0]->id()>0&&hard[1]->id()==ParticleID::g) {
    vector<VectorWaveFunction> g2;
    vector<SpinorWaveFunction> qin;
    vector<SpinorBarWaveFunction> qout;
    SpinorWaveFunction(    qin,hard[0],incoming,false,true);
    VectorWaveFunction(     g2,hard[1],incoming,false,true,true);
    SpinorBarWaveFunction(qout,hard[3],outgoing,true ,true);
    ScalarWaveFunction hout(hard[2],outgoing,true);
    g2[1]=g2[2];
    qgME(qin,g2,hout,qout,true); 
  }
  // qbar g -> H q
  else if(hard[0]->id()<0&&hard[1]->id()==ParticleID::g) {
    vector<VectorWaveFunction> g2;
    vector<SpinorBarWaveFunction> qin;
    vector<SpinorWaveFunction> qout;
    SpinorBarWaveFunction( qin,hard[0],incoming,false,true);
    VectorWaveFunction(     g2,hard[1],incoming,false,true,true);
    SpinorWaveFunction(   qout,hard[3],outgoing,true ,true);
    ScalarWaveFunction hout(hard[2],outgoing,true);
    g2[1]=g2[2];
    qbargME(qin,g2,hout,qout,true); 
  }
  // q qbar to H g
  else if(hard[0]->id()==-hard[1]->id()) {
    vector<SpinorBarWaveFunction> qbar;
    vector<SpinorWaveFunction> q;
    vector<VectorWaveFunction> g4;
    SpinorWaveFunction(    q  ,hard[0],incoming,false,true);
    SpinorBarWaveFunction(qbar,hard[1],incoming,false,true);
    VectorWaveFunction(     g4,hard[3],outgoing,true ,true,true);
    ScalarWaveFunction hout(hard[2],outgoing,true);
    g4[1]=g4[2];
    qqbarME(q,qbar,hout,g4,true); 
  }
  else throw Exception() << "Unknown subprocess in MEPP2HiggsJet::constructVertex()" 
			 << Exception::runerror;
  // construct the vertex
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(_me);
  // set the pointers and to and from the vertex
  for(unsigned int ix=0;ix<4;++ix) {
    (hard[ix]->spinInfo())->
      productionVertex(hardvertex);
  }
}

int MEPP2HiggsJet::nDim() const {
  return 2;
}

void MEPP2HiggsJet::doinit() {
  ME2to2Base::doinit();
  tcPDPtr h0=getParticleData(ParticleID::h0);
  _mh = h0->mass();
  _wh = h0->generateWidth(_mh);
  if(h0->massGenerator()) {
    _hmass=dynamic_ptr_cast<GenericMassGeneratorPtr>(h0->massGenerator());
  }
  if(_shapeopt==2&&!_hmass) throw InitException()
    << "If using the mass generator for the line shape in MEPP2HiggsJet::doinit()"
    << "the mass generator must be an instance of the GenericMassGenerator class"
    << Exception::runerror;
}

CrossSection MEPP2HiggsJet::dSigHatDR() const {
  using Constants::pi;
  InvEnergy2 bwfact;
  Energy moff = mePartonData()[2]->id()==ParticleID::h0 ?
    meMomenta()[2].mass() : meMomenta()[3].mass();
  if(_shapeopt==1) {
    tcPDPtr h0 = mePartonData()[2]->id()==ParticleID::h0 ?
      mePartonData()[2] : mePartonData()[3];
    bwfact = h0->generateWidth(moff)*moff/pi/
      (sqr(sqr(moff)-sqr(_mh))+sqr(_mh*_wh));
  }
  else {
    bwfact = _hmass->BreitWignerWeight(moff);
  }
  return me2()*jacobian()/(16.0*sqr(Constants::pi)*sHat())*sqr(hbarc)*
    (sqr(sqr(moff)-sqr(_mh))+sqr(_mh*_wh))/(_mh*_wh)*bwfact;
}
