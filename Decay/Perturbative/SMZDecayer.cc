// -*- C++ -*-
//
// SMZDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMZDecayer class.
//

#include "SMZDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig/Decay/DecayVertex.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "Herwig/Shower/Base/ShowerTree.h"
#include "Herwig/Shower/Base/ShowerProgenitor.h"
#include "Herwig/Shower/Base/ShowerParticle.h"
#include "Herwig/Shower/Base/Branching.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

const double SMZDecayer::EPS_=0.00000001;

SMZDecayer::SMZDecayer() 
  : quarkWeight_(5,0.), leptonWeight_(6,0.) {
   quarkWeight_[0]  = 0.488029;
   quarkWeight_[1]  = 0.378461;
   quarkWeight_[2]  = 0.488019;
   quarkWeight_[3]  = 0.378027;
   quarkWeight_[4]  = 0.483207;
   leptonWeight_[0] = 0.110709;
   leptonWeight_[1] = 0.220276;
   leptonWeight_[2] = 0.110708;
   leptonWeight_[3] = 0.220276;
   leptonWeight_[4] = 0.110458;
   leptonWeight_[5] = 0.220276;
   // intermediates
   generateIntermediates(false);
   // QED corrections
  hasRealEmissionME(true);
  hasOneLoopME(true);
}

void SMZDecayer::doinit() {
  DecayIntegrator::doinit();
  // get the vertices from the Standard Model object
  tcHwSMPtr hwsm=dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if(!hwsm) throw InitException() << "Must have Herwig StandardModel object in"
				  << "SMZDecayer::doinit()"
				  << Exception::runerror;
  FFZvertex_ = dynamic_ptr_cast<FFVVertexPtr>(hwsm->vertexFFZ());
  FFPvertex_ = hwsm->vertexFFP();
  // make sure they are initialized
  FFZvertex_->init();
  FFPvertex_->init();
  // now set up the decay modes
  DecayPhaseSpaceModePtr mode;
  tPDVector extpart(3);
  vector<double> wgt(0);
  // the Z decay modes
  extpart[0]=getParticleData(ParticleID::Z0);
  // loop over the  quarks and the leptons
  for(int istep=0;istep<11;istep+=10) {
    for(int ix=1;ix<7;++ix) {
      int iy=istep+ix;
      if(iy==6) continue;
      extpart[1] = getParticleData(-iy);
      extpart[2] = getParticleData( iy);
      mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
      if(iy<=6) addMode(mode,  quarkWeight_.at(ix-1),wgt);
      else      addMode(mode,leptonWeight_.at(iy-11),wgt);
    }
  }
}

int SMZDecayer::modeNumber(bool & cc,tcPDPtr parent, 
			    const tPDVector & children) const {
  int imode(-1);
  if(children.size()!=2) return imode;
  int id0=parent->id();
  tPDVector::const_iterator pit = children.begin();
  int id1=(**pit).id();
  ++pit;
  int id2=(**pit).id();
  // Z to quarks or leptons
  cc =false;
  if(id0!=ParticleID::Z0) return imode;
  if(abs(id1)<6&&id1==-id2) {
    imode=abs(id1)-1;
  }
  else if(abs(id1)>=11&&abs(id1)<=16&&id1==-id2) {
    imode=abs(id1)-6;
  }
  cc = false;
  return imode;
}

void SMZDecayer::persistentOutput(PersistentOStream & os) const {
  os << FFZvertex_ << FFPvertex_ << quarkWeight_ << leptonWeight_ << alpha_;
}

void SMZDecayer::persistentInput(PersistentIStream & is, int) {
  is >> FFZvertex_ >> FFPvertex_ >> quarkWeight_ >> leptonWeight_ >> alpha_;
}

ClassDescription<SMZDecayer> SMZDecayer::initSMZDecayer;
// Definition of the static class description member.

void SMZDecayer::Init() {

  static ClassDocumentation<SMZDecayer> documentation
    ("The SMZDecayer class is the implementation of the decay"
     " Z boson to the Standard Model fermions.");

  static ParVector<SMZDecayer,double> interfaceZquarkMax
    ("QuarkMax",
     "The maximum weight for the decay of the Z to quarks",
     &SMZDecayer::quarkWeight_,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<SMZDecayer,double> interfaceZleptonMax
    ("LeptonMax",
     "The maximum weight for the decay of the Z to leptons",
     &SMZDecayer::leptonWeight_,
     0, 0, 0, -10000, 10000, false, false, true);

  static Reference<SMZDecayer,ShowerAlpha> interfaceCoupling
    ("Coupling",
     "Pointer to the object to calculate the coupling for the correction",
     &SMZDecayer::alpha_, false, false, true, false, false);
}


// return the matrix element squared
double SMZDecayer::me2(const int, const Particle & inpart,
			const ParticleVector & decay,
			MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half)));
  int iferm(1),ianti(0);
  if(decay[0]->id()>0) swap(iferm,ianti);
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(_vectors,_rho,
					       const_ptr_cast<tPPtr>(&inpart),
					       incoming,false);
  }
  if(meopt==Terminate) {
    VectorWaveFunction::constructSpinInfo(_vectors,const_ptr_cast<tPPtr>(&inpart),
					  incoming,true,false);
    SpinorBarWaveFunction::
      constructSpinInfo(_wavebar,decay[iferm],outgoing,true);
    SpinorWaveFunction::
      constructSpinInfo(_wave   ,decay[ianti],outgoing,true);
    return 0.;
  }
  SpinorBarWaveFunction::
    calculateWaveFunctions(_wavebar,decay[iferm],outgoing);
  SpinorWaveFunction::
    calculateWaveFunctions(_wave   ,decay[ianti],outgoing);
  // compute the matrix element
  Energy2 scale(sqr(inpart.mass()));
  unsigned int ifm,ia,vhel;
  for(ifm=0;ifm<2;++ifm) {
    for(ia=0;ia<2;++ia) {
      for(vhel=0;vhel<3;++vhel) {
	if(iferm>ianti) (*ME())(vhel,ia,ifm)=
	  FFZvertex_->evaluate(scale,_wave[ia],_wavebar[ifm],_vectors[vhel]);
	else            (*ME())(vhel,ifm,ia)=
	  FFZvertex_->evaluate(scale,_wave[ia],_wavebar[ifm],_vectors[vhel]);
      }
    }
  }
  double output=(ME()->contract(_rho)).real()*UnitRemoval::E2/scale;
  if(abs(decay[0]->id())<=6) output*=3.;
  if(decay[0]->hasColour())      decay[0]->antiColourNeighbour(decay[1]);
  else if(decay[1]->hasColour()) decay[1]->antiColourNeighbour(decay[0]);
  return output;
}

void SMZDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<numberModes();++ix) {
      if(ix<5)       quarkWeight_ [ix   ]=mode(ix)->maxWeight();
      else if(ix<11) leptonWeight_[ix-5 ]=mode(ix)->maxWeight();
    }
  }
}

void SMZDecayer::dataBaseOutput(ofstream & output,
				 bool header) const {
  if(header) output << "update decayers set parameters=\"";
  for(unsigned int ix=0;ix<quarkWeight_.size();++ix) {
    output << "newdef " << name() << ":QuarkMax " << ix << " "
	   << quarkWeight_[ix] << "\n";
  }
  for(unsigned int ix=0;ix<leptonWeight_.size();++ix) {
    output << "newdef " << name() << ":LeptonMax " << ix << " "
	   << leptonWeight_[ix] << "\n";
  }
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

InvEnergy2 SMZDecayer::
realEmissionME(unsigned int,const Particle &parent, 
	       ParticleVector &children,
	       unsigned int iemitter,
	       double ctheta, double stheta,
	       const LorentzRotation & rot1,
	       const LorentzRotation & rot2) {
  // check the number of products and parent
  assert(children.size()==3 && parent.id()==ParticleID::Z0);
  // the electric charge
  double e = sqrt(SM().alphaEM()*4.*Constants::pi);
  // azimuth of the photon
  double phi = children[2]->momentum().phi();
  // wavefunctions for the decaying particle in the rotated dipole frame
  vector<VectorWaveFunction> vec1 = _vectors;
  for(unsigned int ix=0;ix<vec1.size();++ix) {
    vec1[ix].transform(rot1);
    vec1[ix].transform(rot2);
  }
  // wavefunctions for the decaying particle in the rotated rest frame
  vector<VectorWaveFunction> vec2 = _vectors;
  for(unsigned int ix=0;ix<vec1.size();++ix) {
    vec2[ix].transform(rot2);
  }
  // find the outgoing particle and antiparticle
  unsigned int iferm(0),ianti(1);
  if(children[iferm]->id()<0) swap(iferm,ianti);
  // wavefunctions for the particles before the radiation
  // wavefunctions for the outgoing fermion
  SpinorBarWaveFunction wavebartemp;
  Lorentz5Momentum ptemp =  - _wavebar[0].momentum();
  ptemp *= rot2;
  if(ptemp.perp()/ptemp.e()<1e-10) {
    ptemp.setX(ZERO);
    ptemp.setY(ZERO);
  }
  wavebartemp = SpinorBarWaveFunction(ptemp,_wavebar[0].particle(),outgoing);
  // wavefunctions for the outgoing antifermion
  SpinorWaveFunction wavetemp;
  ptemp =  - _wave[0].momentum();
  ptemp *= rot2;
  if(ptemp.perp()/ptemp.e()<1e-10) {
    ptemp.setX(ZERO);
    ptemp.setY(ZERO);
  }
  wavetemp = SpinorWaveFunction(ptemp,_wave[0].particle(),outgoing);
  // loop over helicities
  vector<SpinorWaveFunction> wf_old;
  vector<SpinorBarWaveFunction> wfb_old;
  for(unsigned int ihel=0;ihel<2;++ihel) {
    wavetemp.reset(ihel);
    wf_old.push_back(wavetemp);
    wavebartemp.reset(ihel);
    wfb_old.push_back(wavebartemp);
  }
  // calculate the wave functions for the new fermions
  // ensure the momenta have pT=0
  for(unsigned int ix=0;ix<2;++ix) {
    Lorentz5Momentum ptemp = children[ix]->momentum();
    if(ptemp.perp()/ptemp.e()<1e-10) {
      ptemp.setX(ZERO);
      ptemp.setY(ZERO);
      children[ix]->set5Momentum(ptemp);
    }
  }
  // calculate the wavefunctions
  vector<SpinorBarWaveFunction> wfb;
  SpinorBarWaveFunction::calculateWaveFunctions(wfb,children[iferm],outgoing);
  vector<SpinorWaveFunction> wf;
  SpinorWaveFunction::calculateWaveFunctions   (wf ,children[ianti],outgoing);
  // wave functions for the photons
  vector<VectorWaveFunction> photon;
  VectorWaveFunction::calculateWaveFunctions(photon,children[2],outgoing,true);
  // loop to calculate the matrix elements
  Complex lome[3][2][2],diffme[3][2][2][2],summe[3][2][2][2];
  Energy2 scale(sqr(parent.mass()));
  Complex diff[2]={0.,0.};
  Complex sum [2]={0.,0.};
  for(unsigned int ifm=0;ifm<2;++ifm) {
    for(unsigned int ia=0;ia<2;++ia) {
      for(unsigned int vhel=0;vhel<3;++vhel) {
	// calculation of the leading-order matrix element
	Complex loamp  = FFZvertex_->evaluate(scale,wf_old[ia],
					      wfb_old[ifm],vec2[vhel]);
	Complex lotemp = FFZvertex_->evaluate(scale,wf[ia],
					      wfb[ifm],vec1[vhel]);
	if(iferm>ianti) lome[vhel][ia][ifm] = loamp;
	else            lome[vhel][ifm][ia] = loamp;
	// photon loop for the real emmision terms
	for(unsigned int phel=0;phel<2;++phel) {
	  // radiation from the antifermion
	  // normal case with small angle treatment
	  if(children[2    ]->momentum().z()/
	     children[iferm]->momentum().z()>=ZERO && iemitter == iferm ) {
	    Complex dipole = e*double(children[iferm]->dataPtr()->iCharge())/3.*
	      UnitRemoval::E*loamp*
	      (children[iferm]->momentum()*photon[2*phel].wave())/
	      (children[iferm]->momentum()*children[2]->momentum());
	    // sum and difference
	    SpinorBarWaveFunction foff =
	      FFPvertex_->evaluateSmall(ZERO,3,children[iferm]->dataPtr()->CC(),
					wfb[ifm],photon[2*phel],
					ifm,2*phel,ctheta,phi,stheta,false);
	    diff[0] = FFZvertex_->evaluate(scale,wf[ia],foff,vec1[vhel]) +
	      e*double(children[iferm]->dataPtr()->iCharge())/3.*
	      UnitRemoval::E*(lotemp-loamp)*
	      (children[iferm]->momentum()*photon[2*phel].wave())/
	      (children[iferm]->momentum()*children[2]->momentum());
	    sum [0] = diff[0]+2.*dipole;
	  }
	  // special if fermion backwards
	  else {
	    SpinorBarWaveFunction foff = 
	      FFPvertex_->evaluate(ZERO,3,children[iferm]->dataPtr()->CC(),
				   wfb[ifm],photon[2*phel]);
	    Complex diag = 
	      FFZvertex_->evaluate(scale,wf[ia],foff,vec1[vhel]);
	    Complex dipole = e*double(children[iferm]->dataPtr()->iCharge())/3.*
	      UnitRemoval::E*loamp*
	      (children[iferm]->momentum()*photon[2*phel].wave())/
	      (children[iferm]->momentum()*children[2]->momentum());
	    diff[0] = diag-dipole;
	    sum [0] = diag+dipole;
	  }
	  // radiation from the anti fermion 
	  // small angle case in general
	  if(children[2    ]->momentum().z()/
	     children[ianti]->momentum().z()>=ZERO && iemitter == ianti ) {
	    Complex dipole = e*double(children[ianti]->dataPtr()->iCharge())/3.*
	      UnitRemoval::E*loamp*
	      (children[ianti]->momentum()*photon[2*phel].wave())/
	      (children[ianti]->momentum()*children[2]->momentum());
	    // sum and difference
	    SpinorWaveFunction foff =
	      FFPvertex_->evaluateSmall(ZERO,3,children[ianti]->dataPtr()->CC(),
					wf[ia],photon[2*phel],
					ia,2*phel,ctheta,phi,stheta,false);
	    diff[1] = FFZvertex_->evaluate(scale,foff ,wfb[ifm],vec1[vhel]) +
	      e*double(children[ianti]->dataPtr()->iCharge())/3.*
	      UnitRemoval::E*(lotemp-loamp)*
	      (children[ianti]->momentum()*photon[2*phel].wave())/
	      (children[ianti]->momentum()*children[2]->momentum());
	    sum [1] = diff[1]+2.*dipole;
	  }	    
	  // special if fermion backwards after radiation
	  else {
	    SpinorWaveFunction foff = 
	      FFPvertex_->evaluate(ZERO,3,children[ianti]->dataPtr()->CC(),
				   wf[ia],photon[2*phel]);
	    Complex diag = 
	      FFZvertex_->evaluate(scale,foff ,wfb[ifm],vec1[vhel]);
	    Complex dipole = e*double(children[ianti]->dataPtr()->iCharge())/3.*
	      UnitRemoval::E*loamp*
	      (children[ianti]->momentum()*photon[2*phel].wave())/
	      (children[ianti]->momentum()*children[2]->momentum());
	    // sum and difference
	    diff[1] = diag - dipole;
	    sum [1] = diag + dipole;
	  }
	  // add to me
	  if(iferm>ianti) {
	    diffme[vhel][ia][ifm][phel] = diff[0] + diff[1];
	    summe [vhel][ia][ifm][phel] = sum[0]  + sum[1] ;
	  }
	  else {
	    diffme  [vhel][ifm][ia][phel] = diff[0] + diff[1];
	    summe   [vhel][ifm][ia][phel] = sum[0]  + sum[1] ;
	  }
	}
      }
    }
  }
//   cerr << parent << "\n";
//   for(unsigned int ix=0;ix<children.size();++ix) {
//     cerr << *children[ix] << "\n";
//   }
//   _rho = RhoDMatrix(PDT::Spin1);
  Complex lo(0.),difference(0.);
  for(unsigned int vhel1=0;vhel1<3;++vhel1) {
    for(unsigned int vhel2=0;vhel2<3;++vhel2) {
      for(unsigned int ifm=0;ifm<2;++ifm) {
	for(unsigned int ia=0;ia<2;++ia) {
	  lo += _rho(vhel1,vhel2)*lome[vhel1][ifm][ia]*conj(lome[vhel2][ifm][ia]);
	  for(unsigned int phel=0;phel<2;++phel) {
	    difference += 
	      _rho(vhel1,vhel2)*diffme[vhel1][ifm][ia][phel]*conj(summe[vhel2][ifm][ia][phel]);
	  }
	}
      }
    }
  }
//   // analytic result
//   double iCharge = children[0]->dataPtr()->iCharge()*
//     children[1]->dataPtr()->iCharge()/9.;
//   Energy2 ubar = 2.*children[0]->momentum()*children[2]->momentum();
//   Energy2 tbar = 2.*children[1]->momentum()*children[2]->momentum();
//   double mu2 = sqr(children[1]->mass()/parent.mass());
//   double gL = (FFZvertex_->left() *FFZvertex_->norm()).real();
//   double gR = (FFZvertex_->right()*FFZvertex_->norm()).real();
//   Energy2 den = sqr(parent.mass())*(((sqr(gL)+sqr(gR))*(1-mu2)+6.*mu2*gL*gR));

//   InvEnergy2 anal =  -iCharge*( 2.*(ubar/tbar+tbar/ubar)/sqr(parent.mass())+
// 				4.*mu2/den*((sqr(gL)+sqr(gR))*(1+ubar/tbar+tbar/ubar)
// 					    -2.*gL*gR*(1.+2.*(ubar/tbar+tbar/ubar))));
//   cerr << "testing ratio " << parent.PDGName() 
//        << " " << difference.real()/sqr(e)/lo.real()*UnitRemoval::InvE2/(anal) << "\n"
//        << stheta << " " << ctheta << "\n";
  return difference.real()/sqr(e)/lo.real()*UnitRemoval::InvE2;
}

double SMZDecayer::oneLoopVirtualME(unsigned int,
				    const Particle & parent, 
				    const ParticleVector & children) {
  assert(children.size()==2);
  // velocities of the particles
  double beta = sqrt(1.-4.*sqr(children[0]->mass()/parent.mass()));
  double opb = 1.+beta;
  double omb = 4.*sqr(children[0]->mass()/parent.mass())/opb;
  // couplings
  double gL = (FFZvertex_->left() *FFZvertex_->norm()).real();
  double gR = (FFZvertex_->right()*FFZvertex_->norm()).real();
  double gA = 0.5*(gL-gR);
  double gV = 0.5*(gL+gR);
  // correction terms
  double ln = log(omb/opb);
  double f1 = 1. + ln*beta;
  double fA = 1. + ln/beta;
  InvEnergy f2 = 0.5*sqrt(omb*opb)/parent.mass()/beta*ln;
  // momentum difference for the loop
  Lorentz5Momentum q = children[0]->momentum()-children[1]->momentum();
  if(children[0]->id()<0) q *= -1.;
  // spinors
  vector<LorentzSpinor   <SqrtEnergy> > sp;
  vector<LorentzSpinorBar<SqrtEnergy> > sbar;
  for(unsigned int ix=0;ix<2;++ix) {
    sp  .push_back(   _wave[ix].dimensionedWave());
    sbar.push_back(_wavebar[ix].dimensionedWave());
  }
  // polarization vectors
  vector<LorentzPolarizationVector> pol;
  for(unsigned int ix=0;ix<3;++ix)
    pol.push_back(_vectors[ix].wave());
  // matrix elements
  complex<Energy> lome[3][2][2],loopme[3][2][2];
  for(unsigned int vhel=0;vhel<3;++vhel) {
    for(unsigned int ihel1=0;ihel1<2;++ihel1) {
      for(unsigned int ihel2=0;ihel2<2;++ihel2) {
	complex<Energy> vector = 
	  sp[ihel1].generalCurrent(sbar[ihel2], 1.,1.).dot(pol[vhel]);
	complex<Energy>  axial = 
	  sp[ihel1].generalCurrent(sbar[ihel2],-1.,1.).dot(pol[vhel]);
	complex<Energy2> scalar =
	  sp[ihel1].scalar(sbar[ihel2])*(q*pol[vhel]);
	lome  [vhel][ihel1][ihel2] = gV*   vector-gA*   axial;
	loopme[vhel][ihel1][ihel2] = gV*f1*vector-gA*fA*axial+scalar*f2*gV;
      }
    }
  }
  // sum sums
  complex<Energy2> den(ZERO),num(ZERO);
  for(unsigned int vhel1=0;vhel1<3;++vhel1) {
    for(unsigned int vhel2=0;vhel2<3;++vhel2) {
      for(unsigned int ihel1=0;ihel1<2;++ihel1) {
	for(unsigned int ihel2=0;ihel2<2;++ihel2) {
	  num += _rho(vhel1,vhel2)*
	    (  lome[vhel1][ihel1][ihel2]*conj(loopme[vhel2][ihel1][ihel2])+
	     loopme[vhel1][ihel1][ihel2]*conj(  lome[vhel2][ihel1][ihel2]));
	  den += _rho(vhel1,vhel2)*
	    lome[vhel1][ihel1][ihel2]*conj(lome[vhel2][ihel1][ihel2]);
	}
      }
    }
  }
  // prefactor
  double iCharge = children[0]->dataPtr()->iCharge()*
                   children[1]->dataPtr()->iCharge()/9.;
  double pre = 0.5*SM().alphaEM()*iCharge/Constants::pi;
  // output
  return pre*num.real()/den.real();
}


void SMZDecayer::
initializeMECorrection(ShowerTreePtr tree, double & initial,
		       double & final) {
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
  // get the quark and antiquark
  ParticleVector qq; 
  for(cjt=tree->outgoingLines().begin();cjt!=tree->outgoingLines().end();++cjt)
    qq.push_back(cjt->first->copy());
  // ensure quark first
  if(qq[0]->id()<0) swap(qq[0],qq[1]);
  // centre of mass energy
  d_Q_ = (qq[0]->momentum() + qq[1]->momentum()).m();
  // quark mass
  d_m_ = 0.5*(qq[0]->momentum().m()+qq[1]->momentum().m());
  // set the other parameters
  setRho(sqr(d_m_/d_Q_));
  setKtildeSymm();
  // otherwise can do it
  initial=1.;
  final  =1.;
}

void SMZDecayer::
applyHardMatrixElementCorrection(ShowerTreePtr tree) {
  // get the quark and antiquark
  ParticleVector qq; 
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cit;
  for(cit=tree->outgoingLines().begin();cit!=tree->outgoingLines().end();++cit)
    qq.push_back(cit->first->copy());
  if(!qq[0]->dataPtr()->coloured()) return;
  // ensure quark first
  if(qq[0]->id()<0) swap(qq[0],qq[1]);
  // get the momenta
  vector<Lorentz5Momentum> newfs = applyHard(qq);
  // return if no emission
  if(newfs.size()!=3) return;
  // perform final check to ensure energy greater than constituent mass
  for (int i=0; i<2; i++) {
    if (newfs[i].e() < qq[i]->data().constituentMass()) return;
  }
  if (newfs[2].e() < getParticleData(ParticleID::g)->constituentMass())
    return;
  // set masses
  for (int i=0; i<2; i++) newfs[i].setMass(qq[i]->mass());
  newfs[2].setMass(ZERO);
  // decide which particle emits
  bool firstEmits=
    newfs[2].vect().perp2(newfs[0].vect())<
    newfs[2].vect().perp2(newfs[1].vect());
  // create the new quark, antiquark and gluon
  PPtr newg = getParticleData(ParticleID::g)->produceParticle(newfs[2]);
  PPtr newq,newa;
  if(firstEmits) {
    newq = getParticleData(abs(qq[0]->id()))->produceParticle(newfs[0]);
    newa = new_ptr(Particle(*qq[1]));
    qq[1]->antiColourLine()->removeAntiColoured(newa);
    newa->set5Momentum(newfs[1]);
  }
  else {
    newq = new_ptr(Particle(*qq[0]));
    qq[0]->colourLine()->removeColoured(newq);
    newq->set5Momentum(newfs[0]);
    newa = getParticleData(-abs(qq[0]->id()))->produceParticle(newfs[1]);
  }
  // get the original colour line
  ColinePtr col;
  if(qq[0]->id()>0) col=qq[0]->colourLine();
  else              col=qq[0]->antiColourLine();
  // set the colour lines
  if(firstEmits) {
    col->addColoured(newq);
    col->addAntiColoured(newg);
    newa->colourNeighbour(newg);
  }
  else {
    col->addAntiColoured(newa);
    col->addColoured(newg);
    newq->antiColourNeighbour(newg);
  }
  // change the existing quark and antiquark
  PPtr orig;
  for(cit=tree->outgoingLines().begin();cit!=tree->outgoingLines().end();++cit) {
    if(cit->first->progenitor()->id()==newq->id()) {
      // remove old particles from colour line
      col->removeColoured(cit->first->copy());
      col->removeColoured(cit->first->progenitor());
      // insert new particles
      cit->first->copy(newq);
      ShowerParticlePtr sp(new_ptr(ShowerParticle(*newq,1,true)));
      cit->first->progenitor(sp);
      tree->outgoingLines()[cit->first]=sp;
      cit->first->perturbative(!firstEmits);
      if(firstEmits) orig=cit->first->original();
    }
    else {
      // remove old particles from colour line
      col->removeAntiColoured(cit->first->copy());
      col->removeColoured(cit->first->progenitor());
      // insert new particles
      cit->first->copy(newa);
      ShowerParticlePtr sp(new_ptr(ShowerParticle(*newa,1,true)));
      cit->first->progenitor(sp);
      tree->outgoingLines()[cit->first]=sp;
      cit->first->perturbative(firstEmits);
      if(!firstEmits) orig=cit->first->original();
    }
  }
  // add the gluon
  ShowerParticlePtr sg=new_ptr(ShowerParticle(*newg,1,true));
  ShowerProgenitorPtr gluon=new_ptr(ShowerProgenitor(orig,newg,sg));
  gluon->perturbative(false);
  tree->outgoingLines().insert(make_pair(gluon,sg));
  tree->hardMatrixElementCorrection(true);
}

vector<Lorentz5Momentum> SMZDecayer::
applyHard(const ParticleVector &p) {
  double x, xbar;
  vector<Lorentz5Momentum> fs; 
  // return if no emission
  if (getHard(x, xbar) < UseRandom::rnd() || p.size() != 2) return fs; 
  // centre of mass energy
  Lorentz5Momentum pcm = p[0]->momentum() + p[1]->momentum(); 
  // momenta of quark,antiquark and gluon
  Lorentz5Momentum pq, pa, pg;
  if (p[0]->id() > 0) {
    pq = p[0]->momentum(); 
    pa = p[1]->momentum(); 
  } else {
    pa = p[0]->momentum(); 
    pq = p[1]->momentum(); 
  }
  // boost to boson rest frame
  Boost beta = (pcm.findBoostToCM()); 
  pq.boost(beta);    
  pa.boost(beta);
  // return if fails ?????
  double xg = 2.-x-xbar; 
  if((1.-x)*(1.-xbar)*(1.-xg) < d_rho_*xg*xg) return fs;
  Axis u1, u2, u3;
  // moduli of momenta in units of Q and cos theta
  // stick to q direction?
  // p1 is the one that is kept, p2 is the other fermion, p3 the gluon.
  Energy e1, e2, e3; 
  Energy pp1, pp2, pp3;
  bool keepq = true; 
  if (UseRandom::rnd() > sqr(x)/(sqr(x)+sqr(xbar))) 
    keepq = false; 
  if (keepq) {
    pp1 = d_Q_*sqrt(sqr(x)-4.*d_rho_)/2.;
    pp2 = d_Q_*sqrt(sqr(xbar)-4.*d_rho_)/2.;
    e1 = d_Q_*x/2.; 
    e2 = d_Q_*xbar/2.; 
    u1 = pq.vect().unit();
  } else {
    pp2 = d_Q_*sqrt(sqr(x)-4.*d_rho_)/2.;
    pp1 = d_Q_*sqrt(sqr(xbar)-4.*d_rho_)/2.;
    e2 = d_Q_*x/2.; 
    e1 = d_Q_*xbar/2.; 
    u1 = pa.vect().unit();
  }
  pp3 = d_Q_*xg/2.;       
  e3 = pp3; 
  u2 = u1.orthogonal();
  u2 /= u2.mag();
  u3 = u1.cross(u2);
  u3 /= u3.mag();
  double ct2=-2., ct3=-2.;
  if (pp1 == ZERO || pp2 == ZERO || pp3 == ZERO) {
    bool touched = false;
    if (pp1 == ZERO) {
      ct2 = 1; 
      ct3 = -1; 
      touched = true;
    } 
    if (pp2 == ZERO || pp3 == ZERO) {
      ct2 = 1; 
      ct3 = 1; 
      touched = true;
    }
    if (!touched) 
      throw Exception() << "SMZDecayer::applyHard()"
			<< " did not set ct2/3" 
			<< Exception::abortnow;
  } else {
    ct3 = (sqr(pp1)+sqr(pp3)-sqr(pp2))/(2.*pp1*pp3);
    ct2 = (sqr(pp1)+sqr(pp2)-sqr(pp3))/(2.*pp1*pp2);
  }
  double phi = Constants::twopi*UseRandom::rnd();
  double cphi = cos(phi);
  double sphi = sin(phi); 
  double st2 = sqrt(1.-sqr(ct2));
  double st3 = sqrt(1.-sqr(ct3));
  ThreeVector<Energy> pv1, pv2, pv3; 
  pv1 = pp1*u1;
  pv2 = -ct2*pp2*u1 + st2*cphi*pp2*u2 + st2*sphi*pp2*u3;
  pv3 = -ct3*pp3*u1 - st3*cphi*pp3*u2 - st3*sphi*pp3*u3;
  if (keepq) {
    pq = Lorentz5Momentum(pv1, e1);
    pa = Lorentz5Momentum(pv2, e2);
  } else {
    pa = Lorentz5Momentum(pv1, e1);
    pq = Lorentz5Momentum(pv2, e2);
  }
  pg = Lorentz5Momentum(pv3, e3);
  pq.boost(-beta);
  pa.boost(-beta);
  pg.boost(-beta);
  fs.push_back(pq); 
  fs.push_back(pa); 
  fs.push_back(pg); 
  return fs;
}

double SMZDecayer::getHard(double &x1, double &x2) {
  double w = 0.0;
  double y1 = UseRandom::rnd(),y2 = UseRandom::rnd(); 
  // simply double MC efficiency 
  // -> weight has to be divided by two (Jacobian)
  if (y1 + y2 > 1) {
    y1 = 1.-y1; 
    y2 = 1.-y2;
  }
  bool inSoft = false; 
  if (y1 < 0.25) { 
    if (y2 < 0.25) {
      inSoft = true; 
      if (y1 < y2) {
	y1 = 0.25-y1;
	y2 = y1*(1.5 - 2.*y2);
      }	else {
	y2 = 0.25 - y2;
	y1 = y2*(1.5 - 2.*y1);
      }
    } else {
      if (y2 < y1 + 2.*sqr(y1)) return w;
    }
  } else {
    if (y2 < 0.25) {
      if (y1 < y2 + 2.*sqr(y2)) return w;
    }
  } 
  // inside PS?
  x1 = 1.-y1;
  x2 = 1.-y2;
  if(y1*y2*(1.-y1-y2) < d_rho_*sqr(y1+y2)) return w;
  double k1 = getKfromX(x1, x2);
  double k2 = getKfromX(x2, x1);
  // Is it in the quark emission zone?
  if (k1 < d_kt1_) return 0.0;
  // No...is it in the anti-quark emission zone?
  if (k2 < d_kt2_) return 0.0;  
  // Point is in dead zone: compute q qbar g weight
  w = MEV(x1, x2); 
  // for axial: 
  //  w = MEA(x1, x2); 
  // Reweight soft region
  if (inSoft) { 
    if (y1 < y2) w *= 2.*y1;
    else w *= 2.*y2;
  }
  // alpha and colour factors
  Energy2 pt2 = sqr(d_Q_)*(1.-x1)*(1.-x2);
  w *= 1./3./Constants::pi*alpha_->value(pt2); 
  return w; 
}

bool SMZDecayer::
softMatrixElementVeto(ShowerProgenitorPtr initial,ShowerParticlePtr parent,Branching br) {
  // check we should be applying the veto
  if(parent->id()!=initial->progenitor()->id()||
     br.ids[0]!=br.ids[1]||
     br.ids[2]!=ParticleID::g) return false;
  // calculate pt
  double d_z = br.kinematics->z();
  Energy d_qt = br.kinematics->scale();
  Energy2 d_m2 = parent->momentum().m2();
  Energy pPerp = (1.-d_z)*sqrt( sqr(d_z*d_qt) - d_m2);
  // if not hardest so far don't apply veto
  if(pPerp<initial->highestpT()) return false;
  // calculate the weight
  double weight = 0.;
  if(parent->id()>0) weight = qWeightX(d_qt, d_z);
  else weight = qbarWeightX(d_qt, d_z);
  // compute veto from weight
  bool veto = !UseRandom::rndbool(weight);
  // if vetoing reset the scale
  if(veto) parent->vetoEmission(br.type,br.kinematics->scale());
  // return the veto
  return veto;
}


void SMZDecayer::setRho(double r) 
{ 
  d_rho_ = r;
  d_v_ = sqrt(1.-4.*d_rho_);
}

void SMZDecayer::setKtildeSymm() { 
  d_kt1_ = (1. + sqrt(1. - 4.*d_rho_))/2.;
  setKtilde2();
}

void SMZDecayer::setKtilde2() { 
   double num = d_rho_ * d_kt1_ + 0.25 * d_v_ *(1.+d_v_)*(1.+d_v_);
   double den = d_kt1_ - d_rho_;
   d_kt2_ = num/den;
}

double SMZDecayer::getZfromX(double x1, double x2) {
  double uval = u(x2);
  double num = x1 - (2. - x2)*uval;
  double den = sqrt(x2*x2 - 4.*d_rho_);
  return uval + num/den;
}

double SMZDecayer::getKfromX(double x1, double x2) {
   double zval = getZfromX(x1, x2);
   return (1.-x2)/(zval*(1.-zval));
}

double SMZDecayer::MEV(double x1, double x2) {
  // Vector part
  double num = (x1+2.*d_rho_)*(x1+2.*d_rho_) + (x2+2.*d_rho_)*(x2+2.*d_rho_) 
    - 8.*d_rho_*(1.+2.*d_rho_);
  double den = (1.+2.*d_rho_)*(1.-x1)*(1.-x2);
  return (num/den - 2.*d_rho_/((1.-x1)*(1.-x1)) 
	  - 2*d_rho_/((1.-x2)*(1.-x2)))/d_v_;
}

double SMZDecayer::MEA(double x1, double x2) {
  // Axial part
  double num = (x1+2.*d_rho_)*(x1+2.*d_rho_) + (x2+2.*d_rho_)*(x2+2.*d_rho_) 
    + 2.*d_rho_*((5.-x1-x2)*(5.-x1-x2) - 19.0 + 4*d_rho_);
  double den = d_v_*d_v_*(1.-x1)*(1.-x2);
  return (num/den - 2.*d_rho_/((1.-x1)*(1.-x1)) 
	  - 2*d_rho_/((1.-x2)*(1.-x2)))/d_v_;
}

double SMZDecayer::u(double x2) {
  return 0.5*(1. + d_rho_/(1.-x2+d_rho_));
}

void SMZDecayer::
getXXbar(double kti, double z, double &x, double &xbar) {
  double w = sqr(d_v_) + kti*(-1. + z)*z*(2. + kti*(-1. + z)*z);
  if (w < 0) {
    x = -1.; 
    xbar = -1;
  } else {
    x = (1. + sqr(d_v_)*(-1. + z) + sqr(kti*(-1. + z))*z*z*z 
	 + z*sqrt(w)
	 - kti*(-1. + z)*z*(2. + z*(-2 + sqrt(w))))/
      (1. - kti*(-1. + z)*z + sqrt(w));
    xbar = 1. + kti*(-1. + z)*z;
  }
}

double SMZDecayer::qWeight(double x, double xbar) {
  double rval; 
  double xg = 2. - xbar - x;
  // always return one in the soft gluon region
  if(xg < EPS_) return 1.0;
  // check it is in the phase space
  if((1.-x)*(1.-xbar)*(1.-xg) < d_rho_*xg*xg) return 0.0;
  double k1 = getKfromX(x, xbar);
  double k2 = getKfromX(xbar, x);
  // Is it in the quark emission zone?
  if(k1 < d_kt1_) {
    rval = MEV(x, xbar)/PS(x, xbar);
    // is it also in the anti-quark emission zone?
    if(k2 < d_kt2_) rval *= 0.5;
    return rval;
  }
  return 1.0;
}

double SMZDecayer::qbarWeight(double x, double xbar) {
  double rval; 
  double xg = 2. - xbar - x;
  // always return one in the soft gluon region
  if(xg < EPS_) return 1.0;
  // check it is in the phase space
  if((1.-x)*(1.-xbar)*(1.-xg) < d_rho_*xg*xg) return 0.0;
  double k1 = getKfromX(x, xbar);
  double k2 = getKfromX(xbar, x);
  // Is it in the antiquark emission zone?
  if(k2 < d_kt2_) {
    rval = MEV(x, xbar)/PS(xbar, x);
    // is it also in the quark emission zone?
    if(k1 < d_kt1_) rval *= 0.5;
    return rval;
  }
  return 1.0;
}

double SMZDecayer::qWeightX(Energy qtilde, double z) {
  double x, xb;
  getXXbar(sqr(qtilde/d_Q_), z, x, xb);
  // if exceptionally out of phase space, leave this emission, as there 
  // is no good interpretation for the soft ME correction. 
  if (x < 0 || xb < 0) return 1.0; 
  return qWeight(x, xb); 
}

double SMZDecayer::qbarWeightX(Energy qtilde, double z) {
  double x, xb;
  getXXbar(sqr(qtilde/d_Q_), z, xb, x);
  // see above in qWeightX. 
  if (x < 0 || xb < 0) return 1.0; 
  return qbarWeight(x, xb); 
}

double SMZDecayer::PS(double x, double xbar) {
  double u = 0.5*(1. + d_rho_ / (1.-xbar+d_rho_));
  double z = u + (x - (2.-xbar)*u)/sqrt(xbar*xbar - 4.*d_rho_);
  double brack = (1.+z*z)/(1.-z)- 2.*d_rho_/(1-xbar);
  // interesting: the splitting function without the subtraction
  // term. Actually gives a much worse approximation in the collinear
  // limit.  double brack = (1.+z*z)/(1.-z);
  double den = (1.-xbar)*sqrt(xbar*xbar - 4.*d_rho_);
  return brack/den;
}
