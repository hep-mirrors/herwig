// -*- C++ -*-
//
// SMZDecayer.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
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
#include "Herwig++/Decay/DecayVertex.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "Herwig++/Shower/Base/ShowerTree.h"
#include "Herwig++/Shower/Base/ShowerProgenitor.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"
#include "Herwig++/Shower/Base/Branching.h"
#include "Herwig++/Shower/Base/HardTree.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

const double SMZDecayer::EPS_=0.00000001;

SMZDecayer::SMZDecayer() 
  : quarkWeight_(5,0.), leptonWeight_(6,0.), pTmin_(GeV),
    preFactor_(6.) {
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
  if(!hwsm) throw InitException() << "Must have Herwig++ StandardModel object in"
				  << "SMZDecayer::doinit()"
				  << Exception::runerror;
  FFZvertex_ = dynamic_ptr_cast<FFVVertexPtr>(hwsm->vertexFFZ());
  FFPvertex_ = hwsm->vertexFFP();
  FFGVertex_ = hwsm->vertexFFG();
  // make sure they are initialized
  FFZvertex_->init();
  FFPvertex_->init();
  FFGVertex_->init();
  // extract ParticleData boejct for the gluon
  gluon_ = getParticleData(ParticleID::g);
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
  os << FFZvertex_ << FFPvertex_  << FFGVertex_ << gluon_ 
     << quarkWeight_ << leptonWeight_ << alpha_ << ounit(pTmin_,GeV)
     << preFactor_;
}

void SMZDecayer::persistentInput(PersistentIStream & is, int) {
  is >> FFZvertex_ >> FFPvertex_ >> FFGVertex_ >> gluon_  
     >> quarkWeight_ >> leptonWeight_ >> alpha_ >> iunit(pTmin_,GeV)
     >> preFactor_;
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
  int iferm(1),ianti(0);
  if(decay[0]->id()>0) swap(iferm,ianti);
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(_vectors,_rho,
					       const_ptr_cast<tPPtr>(&inpart),
					       incoming,false);
    ME(DecayMatrixElement(PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half));
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
	if(iferm>ianti) ME()(vhel,ia,ifm)=
	  FFZvertex_->evaluate(scale,_wave[ia],_wavebar[ifm],_vectors[vhel]);
	else            ME()(vhel,ifm,ia)=
	  FFZvertex_->evaluate(scale,_wave[ia],_wavebar[ifm],_vectors[vhel]);
      }
    }
  }
  double output=(ME().contract(_rho)).real()*UnitRemoval::E2/scale;
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
softMatrixElementVeto(ShowerProgenitorPtr initial,ShowerParticlePtr parent,Branching br)
{
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
  // if not vetoed reset max
  if(!veto) initial->highestpT(pPerp);
  // if vetoing reset the scale
  if(veto) parent->evolutionScale(br.kinematics->scale());
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

HardTreePtr SMZDecayer::generateHardest(ShowerTreePtr tree) {
  ShowerProgenitorPtr ZProgenitor = tree->incomingLines().begin() ->first; 
  // outgoing
  ShowerProgenitorPtr 
    qkProgenitor = tree->outgoingLines().begin() ->first,
    qbProgenitor = tree->outgoingLines().rbegin()->first;
  if(!qkProgenitor->progenitor()->dataPtr()->coloured())
    return HardTreePtr();
  // get the order right 
  if(qkProgenitor->id()<0) swap(qkProgenitor,qbProgenitor);
  // extract the momenta 
  loMomenta_.resize(0);
  loMomenta_.push_back(ZProgenitor->progenitor()->momentum());
  loMomenta_.push_back(qkProgenitor->progenitor()->momentum());
  loMomenta_.push_back(qbProgenitor->progenitor()->momentum());
  // and ParticleData objects
  partons_.resize(0);
  partons_.push_back(ZProgenitor->progenitor()->dataPtr());
  partons_.push_back(qkProgenitor->progenitor()->dataPtr());
  partons_.push_back(qbProgenitor->progenitor()->dataPtr());
  partons_.push_back(gluon_);
  // boost from lab to CMS frame with outgoing particles
  // along the z axis
  LorentzRotation eventFrame( ( loMomenta_[1] + loMomenta_[2] ).findBoostToCM() );
  Lorentz5Momentum spectator = eventFrame*loMomenta_[2];
  eventFrame.rotateZ( -spectator.phi() );
  eventFrame.rotateY( -spectator.theta()  );
  eventFrame.invert();
  // mass of the final-state system
  Energy2 M2 = (loMomenta_[1]+loMomenta_[2]).m2();
  Energy  M  = sqrt(M2);
  double mu1 = loMomenta_[1].mass()/M;
  double mu2 = loMomenta_[2].mass()/M;
  double mu12 = sqr(mu1), mu22 = sqr(mu2);
  double lambda = sqrt(1.+sqr(mu12)+sqr(mu22)-2.*mu12-2.*mu22-2.*mu12*mu22);
  // max pT
  Energy pTmax = 0.5*sqrt(M2)*
    (1.-sqr(loMomenta_[1].mass()+loMomenta_[2].mass())/M2);
  // max y
  double ymax = acosh(pTmax/pTmin_);
  // prefactor for the overestimate of the Sudakov
  double a = 4./3.*alpha_->overestimateValue()/Constants::twopi*
    2.*ymax*preFactor_;
  // variables for the emission
  Energy pT[2];
  double  y[2],phi[2],x3[2],x1[2][2],x2[2][2];
  double contrib[2][2];
  // storage of the real emmision momenta
  vector<Lorentz5Momentum> realMomenta[2][2]=
    {{vector<Lorentz5Momentum>(4),vector<Lorentz5Momentum>(4)},
     {vector<Lorentz5Momentum>(4),vector<Lorentz5Momentum>(4)}};
  for(unsigned int ix=0;ix<2;++ix)
    for(unsigned int iy=0;iy<2;++iy)
      realMomenta[ix][iy][0] = loMomenta_[0];
  // generate the emission
  for(unsigned int ix=0;ix<2;++ix) {
    if(ix==1) {
      swap(mu1 ,mu2 );
      swap(mu12,mu22);
    }
    pT[ix] = pTmax;
    y [ix] = 0.;
    bool reject = true;
    do {
      // generate pT
      pT[ix] *= pow(UseRandom::rnd(),1./a);
      if(pT[ix]<pTmin_) {
        pT[ix] = -GeV;
        break;
      }
      // generate y
      y[ix] = -ymax+2.*UseRandom::rnd()*ymax;
      // generate phi
      phi[ix] = UseRandom::rnd()*Constants::twopi;
      // calculate x3 and check in allowed region
      x3[ix] = 2.*pT[ix]*cosh(y[ix])/M;
      if(x3[ix] < 0. || x3[ix] > 1. -sqr( mu1 + mu2 ) ) continue;
      // find the possible solutions for x1
      double xT2 = sqr(2./M*pT[ix]);
      double root = (-sqr(x3[ix])+xT2)*
	(xT2*mu22+2.*x3[ix]-sqr(mu12)+2.*mu22+2.*mu12-sqr(x3[ix])-1.
	 +2.*mu12*mu22-sqr(mu22)-2.*mu22*x3[ix]-2.*mu12*x3[ix]);
      double c1=2.*sqr(x3[ix])-4.*mu22-6.*x3[ix]+4.*mu12-xT2*x3[ix]
	+2.*xT2-2.*mu12*x3[ix]+2.*mu22*x3[ix]+4.;
      if(root<0.) continue;
      x1[ix][0] = 1./(4.-4.*x3[ix]+xT2)*(c1-2.*sqrt(root));
      x1[ix][1] = 1./(4.-4.*x3[ix]+xT2)*(c1+2.*sqrt(root));
      // change sign of y if 2nd particle emits
      if(ix==1) y[ix] *=-1.;
      // loop over the solutions
      for(unsigned int iy=0;iy<2;++iy) {
	contrib[ix][iy]=0.;
	// check x1 value allowed
	if(x1[ix][iy]<2.*mu1||x1[ix][iy]>1.+mu12-mu22) continue;
	// calculate x2 value and check allowed
	x2[ix][iy] = 2.-x3[ix]-x1[ix][iy];
	double root = max(0.,sqr(x1[ix][iy])-4.*mu12);
	root = sqrt(root);
	double x2min = 1.+mu22-mu12
	  -0.5*(1.-x1[ix][iy]+mu12-mu22)/(1.-x1[ix][iy]+mu12)*(x1[ix][iy]-2.*mu12+root);
	double x2max = 1.+mu22-mu12
	  -0.5*(1.-x1[ix][iy]+mu12-mu22)/(1.-x1[ix][iy]+mu12)*(x1[ix][iy]-2.*mu12-root);
	if(x2[ix][iy]<x2min||x2[ix][iy]>x2max) continue;
	// check the z components
	double z1 =  sqrt(sqr(x1[ix][iy])-4.*mu12-xT2);
	double z2 = -sqrt(sqr(x2[ix][iy])-4.*mu22);
	double z3 =  pT[ix]*sinh(y[ix])*2./M;
	if(ix==1) z3 *=-1.;
	if(abs(-z1+z2+z3)<1e-9) z1 *= -1.;
	if(abs(z1+z2+z3)>1e-3) continue;
	// construct the momenta
	realMomenta[ix][iy][3] =
	  Lorentz5Momentum(pT[ix]*cos(phi[ix]),pT[ix]*sin(phi[ix]),
			   pT[ix]*sinh(y[ix]) ,pT[ix]*cosh(y[ix]),ZERO);
	if(ix==0) {
	  realMomenta[ix][iy][1] =
	    Lorentz5Momentum(-pT[ix]*cos(phi[ix]),-pT[ix]*sin(phi[ix]),
			     z1*0.5*M,x1[ix][iy]*0.5*M,M*mu1);
	  realMomenta[ix][iy][2] =
	    Lorentz5Momentum(ZERO,ZERO, z2*0.5*M,x2[ix][iy]*0.5*M,M*mu2);
	}
	else {
	  realMomenta[ix][iy][1] =
	    Lorentz5Momentum(ZERO,ZERO,-z2*0.5*M,x2[ix][iy]*0.5*M,M*mu2);
	  realMomenta[ix][iy][2] =
	    Lorentz5Momentum(-pT[ix]*cos(phi[ix]),-pT[ix]*sin(phi[ix]),
			     -z1*0.5*M,x1[ix][iy]*0.5*M,M*mu1);
	}
	// boost the momenta back to the lab
	for(unsigned int iz=1;iz<4;++iz)
	  realMomenta[ix][iy][iz] *= eventFrame;
	// jacobian and prefactors for the weight
	Energy J = M/sqrt(xT2)*abs(-x1[ix][iy]*x2[ix][iy]+2.*mu22*x1[ix][iy]
				   +x2[ix][iy]+x2[ix][iy]*mu12+mu22*x2[ix][iy]
				   -sqr(x2[ix][iy]))
	  /pow(sqr(x2[ix][iy])-4.*mu22,1.5);
	// prefactors etc
	contrib[ix][iy] = 0.5*pT[ix]/J/preFactor_/lambda;
	// matrix element piece
	contrib[ix][iy] *= alpha_->ratio(sqr(pT[ix]))*
	  meRatio(partons_,realMomenta[ix][iy],ix,false);
      }
      if(contrib[ix][0]+contrib[ix][1]>1.)
	cerr << "testing weight greater than one " 
	     << contrib[ix][0]+contrib[ix][1] << "\n";
      reject =  UseRandom::rnd() > contrib[ix][0] + contrib[ix][1];
    }
    while (reject);
    if(pT[ix]<pTmin_) pT[ix] = -GeV;
  }
  if(pT[0]<ZERO && pT[1]<ZERO) {
    qkProgenitor->maximumpT(pTmin_,ShowerInteraction::QCD);
    qbProgenitor->maximumpT(pTmin_,ShowerInteraction::QCD);
    return HardTreePtr();
  }
  // now pick the emission with highest pT
  vector<Lorentz5Momentum> emmision;
  unsigned int iemit=0,ispect=0;
  Energy pTveto;
  if(pT[0]>pT[1]) {
    iemit  = 1;
    ispect = 2;
    pTveto = pT[0];
    if(UseRandom::rnd()<contrib[0][0]/(contrib[0][0]+contrib[0][1]))
      emmision = realMomenta[0][0];
    else
      emmision = realMomenta[0][1];
  }
  else {
    iemit  = 2;
    ispect = 1;
    pTveto = pT[1];
    if(UseRandom::rnd()<contrib[1][0]/(contrib[1][0]+contrib[1][1]))
      emmision = realMomenta[1][0];
    else
      emmision = realMomenta[1][1];
  }
  // Make the particles for the hard tree
  ShowerParticleVector hardParticles;
  for(unsigned int ix=0;ix<partons_.size();++ix) {
    hardParticles.push_back(new_ptr(ShowerParticle(partons_[ix],ix>=1)));
    hardParticles.back()->set5Momentum(emmision[ix]);
  }
  ShowerParticlePtr parent(new_ptr(ShowerParticle(partons_[iemit],true)));
  Lorentz5Momentum parentMomentum(emmision[iemit]+emmision[3]);
  parentMomentum.setMass(partons_[iemit]->mass());
  parent->set5Momentum(parentMomentum);
  // Create the vectors of HardBranchings to create the HardTree:
  vector<HardBranchingPtr> spaceBranchings,allBranchings;
  // Incoming boson:
  spaceBranchings.push_back(new_ptr(HardBranching(hardParticles[0],SudakovPtr(),
						  HardBranchingPtr(),
						  HardBranching::Incoming)));
  allBranchings.push_back(spaceBranchings.back());
  // Outgoing particles from hard emission:
  HardBranchingPtr spectatorBranch(new_ptr(HardBranching(hardParticles[ispect],
 							 SudakovPtr(),HardBranchingPtr(),
 							 HardBranching::Outgoing)));
  HardBranchingPtr emitterBranch(new_ptr(HardBranching(parent,SudakovPtr(),
						       HardBranchingPtr(),
						       HardBranching::Outgoing)));
  emitterBranch->addChild(new_ptr(HardBranching(hardParticles[iemit], 
 						SudakovPtr(),HardBranchingPtr(),
 						HardBranching::Outgoing)));
  emitterBranch->addChild(new_ptr(HardBranching(hardParticles[3],
 						SudakovPtr(),HardBranchingPtr(),
 						HardBranching::Outgoing)));
  if(iemit==0) {
    allBranchings.push_back(emitterBranch);
    allBranchings.push_back(spectatorBranch);
  } 
  else {
    allBranchings.push_back( spectatorBranch );
    allBranchings.push_back( emitterBranch );
  }
  // Make the HardTree from the HardBranching vectors.
  HardTreePtr hardtree = new_ptr(HardTree(allBranchings,spaceBranchings,
 					   ShowerInteraction::QCD));
  // Set the maximum pt for all other emissions
  qkProgenitor->maximumpT(pTveto,ShowerInteraction::QCD);
  qbProgenitor->maximumpT(pTveto,ShowerInteraction::QCD);
  // Connect the particles with the branchings in the HardTree
  hardtree->connect( qkProgenitor->progenitor(), allBranchings[1] );
  hardtree->connect( qbProgenitor->progenitor(), allBranchings[2] );
  // colour flow
  ColinePtr newline=new_ptr(ColourLine());
  for(set<HardBranchingPtr>::const_iterator cit=hardtree->branchings().begin();
      cit!=hardtree->branchings().end();++cit) {
    if((**cit).branchingParticle()->dataPtr()->iColour()==PDT::Colour3)
      newline->addColoured((**cit).branchingParticle());
    else if((**cit).branchingParticle()->dataPtr()->iColour()==PDT::Colour3bar)
      newline->addAntiColoured((**cit).branchingParticle());
  }
  // Return the HardTree
  return hardtree;
}

double SMZDecayer::meRatio(vector<cPDPtr> partons, 
			   vector<Lorentz5Momentum> momenta,
			   unsigned int iemitter,bool subtract) const {
  Lorentz5Momentum q = momenta[1]+momenta[2]+momenta[3];
  Energy2 Q2=q.m2();
  Energy2 lambda = sqrt((Q2-sqr(momenta[1].mass()+momenta[2].mass()))*
			(Q2-sqr(momenta[1].mass()-momenta[2].mass())));
  InvEnergy2 D[2];
  Energy2 lome[2];
  for(unsigned int iemit=0;iemit<2;++iemit) {
    unsigned int ispect = iemit==0 ? 1 : 0;
    Energy2 pipj = momenta[3      ] * momenta[1+iemit ];
    Energy2 pipk = momenta[3      ] * momenta[1+ispect];
    Energy2 pjpk = momenta[1+iemit] * momenta[1+ispect];
    double y = pipj/(pipj+pipk+pjpk);
    double z = pipk/(     pipk+pjpk);
    Energy mij = sqrt(2.*pipj+sqr(momenta[1+iemit].mass()));
    Energy2 lamB = sqrt((Q2-sqr(mij+momenta[1+ispect].mass()))*
			(Q2-sqr(mij-momenta[1+ispect].mass())));
    Energy2 Qpk = q*momenta[1+ispect];
    Lorentz5Momentum pkt = 
      lambda/lamB*(momenta[1+ispect]-Qpk/Q2*q)
      +0.5/Q2*(Q2+sqr(momenta[1+ispect].mass())-sqr(momenta[1+ispect].mass()))*q;
    Lorentz5Momentum pijt = q-pkt;
    double muj = momenta[1+iemit ].mass()/sqrt(Q2);
    double muk = momenta[1+ispect].mass()/sqrt(Q2);
    double vt = sqrt((1.-sqr(muj+muk))*(1.-sqr(muj-muk)))/(1.-sqr(muj)-sqr(muk));
    double v  = sqrt(sqr(2.*sqr(muk)+(1.-sqr(muj)-sqr(muk))*(1.-y))-4.*sqr(muk))
      /(1.-y)/(1.-sqr(muj)-sqr(muk));
    // dipole term
    D[iemit] = 0.5/pipj*(2./(1.-(1.-z)*(1.-y))
			 -vt/v*(2.-z+sqr(momenta[1+iemit].mass())/pipj));
    // matrix element
    vector<Lorentz5Momentum> lomom(4);
    lomom[0] = momenta[0];
    if(iemit==0) {
      lomom[1] = pijt;
      lomom[2] = pkt ;
    }
    else {
      lomom[2] = pijt;
      lomom[1] = pkt ;
    }
    lome[iemit]  = loME(partons,lomom,false)/3.;
  }
  if(isnan(D[0]*GeV2) || isnan(D[1]*GeV2)) {
    for(unsigned int ix=0;ix<partons.size();++ix)
      cerr << partons[ix]->PDGName() << " " << momenta[ix]/GeV << "\n";
  }
  assert(!isnan(D[0]*GeV2) && !isnan(D[1]*GeV2));
  assert(!isinf(D[0]*GeV2) && !isinf(D[1]*GeV2));
  InvEnergy2 ratio = realME(partons,momenta)*abs(D[iemitter])/(abs(D[0]*lome[0])+abs(D[1]*lome[1]));
  assert(!isnan(ratio*GeV2));;
  assert(!isinf(ratio*GeV2));
  if(subtract) {
    return Q2*(ratio-2.*D[iemitter]);
  }
  else
    return Q2*ratio;
}

Energy2 SMZDecayer::loME(const vector<cPDPtr> & partons, 
			const vector<Lorentz5Momentum> & momenta,
			bool first) const {
  // compute the spinors
  vector<SpinorWaveFunction> aout;
  vector<SpinorBarWaveFunction>  fout;
  SpinorBarWaveFunction qkout(momenta[1],partons[1],outgoing);
  SpinorWaveFunction    qbout(momenta[2],partons[2],outgoing);
  for(unsigned int ix=0;ix<2;++ix) {
    qkout.reset(ix);
    fout.push_back(qkout);
    qbout.reset(ix);
    aout.push_back(qbout);
  }
  Energy2 scale(sqr(momenta[0].mass()));
  DecayMatrixElement me(PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half);
  for(unsigned int ifm=0;ifm<2;++ifm) {
    for(unsigned int ia=0;ia<2;++ia) {
      for(unsigned int vhel=0;vhel<3;++vhel) {
	me(vhel,ifm,ia)=
 	  FFZvertex_->evaluate(scale,aout[ia],fout[ifm],_vectors[vhel]);
      }
    }
  }
  Energy2 output=(me.contract(_rho)).real()*UnitRemoval::E2;
  if(abs(partons[1]->id())<=6) output*=3.;
  return output;
}

double SMZDecayer::realME(const vector<cPDPtr> & partons, 
			  const vector<Lorentz5Momentum> & momenta) const {
  // compute the spinors
  vector<SpinorWaveFunction>    aout;
  vector<SpinorBarWaveFunction> fout;
  vector<VectorWaveFunction> gout;
  SpinorBarWaveFunction qkout(momenta[1],partons[1],outgoing);
  SpinorWaveFunction    qbout(momenta[2],partons[2],outgoing);
  VectorWaveFunction    gluon(momenta[3],partons[3],outgoing);
  for(unsigned int ix=0;ix<2;++ix) {
    qkout.reset(ix);
    fout.push_back(qkout);
    qbout.reset(ix);
    aout.push_back(qbout);
    gluon.reset(2*ix);
    gout.push_back(gluon);
  }
  vector<Complex> diag(2,0.);
  Energy2 scale = sqr(momenta[0].mass());
  DecayMatrixElement output(PDT::Spin1,PDT::Spin1Half,
			    PDT::Spin1Half,PDT::Spin1);
  for(unsigned int outhel1=0;outhel1<2;++outhel1) {
    for(unsigned int outhel2=0;outhel2<2;++outhel2) {
      for(unsigned int outhel3=0;outhel3<2;++outhel3) {
	SpinorBarWaveFunction off1 =
	  FFGVertex_->evaluate(scale,3,partons[1],fout[outhel1],gout[outhel3]);
	SpinorWaveFunction    off2 = 
	  FFGVertex_->evaluate(scale,3,partons[2],aout[outhel2],gout[outhel3]);
	
	
	for(unsigned int vhel=0;vhel<3;++vhel) {
	  diag[0] = FFZvertex_->evaluate(scale,aout[outhel2],off1,_vectors[vhel]);
	  diag[1] = FFZvertex_->evaluate(scale,off2,fout[outhel1],_vectors[vhel]);
	  // sum of diagrams
	  output(vhel,outhel1,outhel2,outhel3) = diag[0] + diag[1];
	}
      }		
    }
  }
  // spin average
  double total = (output.contract(_rho)).real();
  // divide out the coupling
  total /= norm(FFGVertex_->norm());
  // return the total
  return total;
}
