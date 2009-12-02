// -*- C++ -*-
//
// SMZDecayer.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
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
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig++/Decay/DecayVertex.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

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
  if(!hwsm) throw InitException() << "Must have Herwig++ StandardModel object in"
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
      // check that the combination of particles is allowed
      if(!FFZvertex_->allowed(-iy,iy,ParticleID::Z0))
	throw InitException() << "SMZDecayer::doinit() the Z vertex " 
			      << "cannot handle all the modes" 
			      << Exception::abortnow;
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
  os << FFZvertex_ << FFPvertex_ << quarkWeight_ << leptonWeight_;
}

void SMZDecayer::persistentInput(PersistentIStream & is, int) {
  is >> FFZvertex_ >> FFPvertex_ >> quarkWeight_ >> leptonWeight_;
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
    output << "set " << name() << ":QuarkMax " << ix << " "
	   << quarkWeight_[ix] << "\n";
  }
  for(unsigned int ix=0;ix<leptonWeight_.size();++ix) {
    output << "set " << name() << ":LeptonMax " << ix << " "
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
	      FFPvertex_->evaluateSmall(ZERO,3,children[iferm]->dataPtr(),
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
	      FFPvertex_->evaluate(ZERO,3,children[iferm]->dataPtr(),
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
	      FFPvertex_->evaluateSmall(ZERO,3,children[ianti]->dataPtr(),
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
	      FFPvertex_->evaluate(ZERO,3,children[ianti]->dataPtr(),
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
