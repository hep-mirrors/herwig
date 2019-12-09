// -*- C++ -*-
//
// SMZDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMZDecayer class.
//

#include "SMZDecayer.h"
#include "Herwig/Utilities/Maths.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig/Decay/DecayVertex.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "Herwig/Shower/RealEmissionProcess.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
#include <numeric>

using namespace Herwig;
using namespace ThePEG::Helicity;

const double SMZDecayer::EPS_=0.00000001;

SMZDecayer::SMZDecayer() 
  : quarkWeight_(5,0.), leptonWeight_(6,0.), CF_(4./3.),
    NLO_(false) {
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
  PerturbativeDecayer::doinit();
  // get the vertices from the Standard Model object
  tcHwSMPtr hwsm=dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if(!hwsm) throw InitException() << "Must have Herwig StandardModel object in"
				  << "SMZDecayer::doinit()"
				  << Exception::runerror;
  // cast the vertices
  FFZVertex_ = dynamic_ptr_cast<FFVVertexPtr>(hwsm->vertexFFZ());
  FFZVertex_->init();
  FFGVertex_ = hwsm->vertexFFG();
  FFGVertex_->init();
  FFPVertex_ = hwsm->vertexFFP();
  FFPVertex_->init();
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
  os << FFZVertex_ << FFPVertex_ << FFGVertex_
     << quarkWeight_ << leptonWeight_ << NLO_
     << gluon_;
}

void SMZDecayer::persistentInput(PersistentIStream & is, int) {
  is >> FFZVertex_ >> FFPVertex_ >> FFGVertex_
     >> quarkWeight_ >> leptonWeight_ >> NLO_
     >> gluon_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SMZDecayer,PerturbativeDecayer>
describeHerwigSMZDecayer("Herwig::SMZDecayer", "HwPerturbativeDecay.so");

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

  static Switch<SMZDecayer,bool> interfaceNLO
    ("NLO",
     "Whether to return the LO or NLO result",
     &SMZDecayer::NLO_, false, false, false);
  static SwitchOption interfaceNLOLO
    (interfaceNLO,
     "No",
     "Leading-order result",
     false);
  static SwitchOption interfaceNLONLO
    (interfaceNLO,
     "Yes",
     "NLO result",
     true);

}

// return the matrix element squared
double SMZDecayer::me2(const int, const Particle & part,
			const ParticleVector & decay,
			MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half)));
  int iferm(1),ianti(0);
  if(decay[0]->id()>0) swap(iferm,ianti);
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(_vectors,_rho,
					       const_ptr_cast<tPPtr>(&part),
					       incoming,false);
    // fix rho if no correlations
    fixRho(_rho);
  }
  if(meopt==Terminate) {
    VectorWaveFunction::constructSpinInfo(_vectors,const_ptr_cast<tPPtr>(&part),
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
  Energy2 scale(sqr(part.mass()));
  unsigned int ifm,ia,vhel;
  for(ifm=0;ifm<2;++ifm) {
    for(ia=0;ia<2;++ia) {
      for(vhel=0;vhel<3;++vhel) {
	if(iferm>ianti) (*ME())(vhel,ia,ifm)=
	  FFZVertex_->evaluate(scale,_wave[ia],_wavebar[ifm],_vectors[vhel]);
	else            (*ME())(vhel,ifm,ia)=
	  FFZVertex_->evaluate(scale,_wave[ia],_wavebar[ifm],_vectors[vhel]);
      }
    }
  }
  double output=(ME()->contract(_rho)).real()*UnitRemoval::E2/scale;
  if(abs(decay[0]->id())<=6) output*=3.;
  if(decay[0]->hasColour())      decay[0]->antiColourNeighbour(decay[1]);
  else if(decay[1]->hasColour()) decay[1]->antiColourNeighbour(decay[0]);
  // if LO return
  if(!NLO_) return output;  // check decay products coloured, otherwise return
  if(!decay[0]->dataPtr()->coloured()) return output;
  // inital masses, couplings  etc
  // fermion mass
  Energy particleMass = decay[0]->dataPtr()->mass();
  // Z mass
  mZ_ = part.mass();
  // strong coupling
  aS_ = SM().alphaS(sqr(mZ_));
  // reduced mass
  mu_  = particleMass/mZ_;
  mu2_ = sqr(mu_);
  // scale
  scale_ = sqr(mZ_);
  // compute the spinors
  vector<SpinorWaveFunction>    aout;
  vector<SpinorBarWaveFunction> fout;
  vector<VectorWaveFunction>    vin;
  SpinorBarWaveFunction qkout(decay[0]->momentum(),decay[0]->dataPtr(),outgoing);
  SpinorWaveFunction    qbout(decay[1]->momentum(),decay[1]->dataPtr(),outgoing);
  VectorWaveFunction    zin  (part.momentum()     ,part.dataPtr()     ,incoming);
  for(unsigned int ix=0;ix<2;++ix){
    qkout.reset(ix);
    fout.push_back(qkout);
    qbout.reset(ix);
    aout.push_back(qbout);
  }
  for(unsigned int ix=0;ix<3;++ix){
    zin.reset(ix);
    vin.push_back(zin);
  }
  // temporary storage of the different diagrams
  // sum over helicities to get the matrix element
  double total=0.;
  if(mu_!=0.) {
    LorentzPolarizationVector momDiff = 
      (decay[0]->momentum()-decay[1]->momentum())/2./
      (decay[0]->momentum().mass()+decay[1]->momentum().mass());
    // scalars
    Complex scalar1 = zin.wave().dot(momDiff);
    for(unsigned int outhel1=0;outhel1<2;++outhel1) {
      for(unsigned int outhel2=0;outhel2<2;++outhel2) {		
	for(unsigned int inhel=0;inhel<3;++inhel) {
	  // first the LO bit
	  Complex diag1 = FFZVertex_->evaluate(scale_,aout[outhel2],
					       fout[outhel1],vin[inhel]);
	  // extra stuff for NLO 
	  LorentzPolarizationVector left  = 
	    aout[outhel2].wave().leftCurrent(fout[outhel1].wave());
	  LorentzPolarizationVector right = 
	    aout[outhel2].wave().rightCurrent(fout[outhel1].wave());
	  Complex scalar = 
	    aout[outhel2].wave().scalar(fout[outhel1].wave());
	  // nlo specific pieces
	  Complex diag3 =
	    Complex(0.,1.)*FFZVertex_->norm()*
	    (FFZVertex_->right()*( left.dot(zin.wave())) +
	     FFZVertex_-> left()*(right.dot(zin.wave())) -
	     ( FFZVertex_-> left()+FFZVertex_->right())*scalar1*scalar);
	  // nlo piece
	  total += real(diag1*conj(diag3) + diag3*conj(diag1));
	}
      }
    }
    // rescale
    total *= UnitRemoval::E2/scale_;
  }
  else {
    total = ZERO;
  }
  // now for the NLO bit
  double mu4 = sqr(mu2_);
  double lmu = mu_!=0. ? log(mu_) : 0.;
  double v = sqrt(1.-4.*mu2_),v2(sqr(v));
  double omv = 4.*mu2_/(1.+v);
  double f1,f2,fNS,VNS;
  double r = omv/(1.+v);
  double lr = mu_!=0. ? log(r) : 0.;
  // normal form
  if(mu_>1e-4) {
    f1 = CF_*aS_/Constants::pi*
      ( +1. + 3.*log(0.5*(1.+v)) - 1.5*log(0.5*(1.+v2)) + sqr(Constants::pi)/6.
	- 0.5*sqr(lr) - (1.+v2)/v*(lr*log(1.+v2) + sqr(Constants::pi)/12. 
				       -0.5*log(4.*mu2_)*lr + 0.25*sqr(lr)));
    fNS = -0.5*(1.+2.*v2)*lr/v + 1.5*lr - 2./3.*sqr(Constants::pi) + 0.5*sqr(lr)
      + (1.+v2)/v*(Herwig::Math::ReLi2(r) + sqr(Constants::pi)/3. - 0.25*sqr(lr) 
		   + lr*log((2.*v/ (1.+v))));
    VNS = 1.5*log(0.5*(1.+v2)) 
      + 0.5*(1.+v2)/v*( 2.*lr*log(2.*(1.+v2)/sqr(1.+v))  
			+ 2.*Herwig::Math::ReLi2(sqr(r)) 
			- 2.*Herwig::Math::ReLi2(2.*v/(1.+v)) - sqr(Constants::pi)/6.)
      + log(1.-mu_) - 2.*log(1.-2.*mu_) - 4.*mu2_/(1.+v2)*log(mu_/(1.-mu_)) 
      - mu_/(1.-mu_)
      + 4.*(2.*mu2_-mu_)/(1.+v2) + 0.5*sqr(Constants::pi); 
    f2 = CF_*aS_/Constants::pi*mu2_*lr/v;
  }
  // small mass limit
  else {
    f1 = -CF_*aS_/Constants::pi/6.*
      ( - 6. - 24.*lmu*mu2_ - 15.*mu4 - 12.*mu4*lmu - 24.*mu4*sqr(lmu) 
	+ 2.*mu4*sqr(Constants::pi) - 12.*mu2_*mu4 - 96.*mu2_*mu4*sqr(lmu) 
	+ 8.*mu2_*mu4*sqr(Constants::pi) - 80.*mu2_*mu4*lmu);
    fNS = - mu2_/18.*( + 36.*lmu - 36. - 45.*mu2_ + 216.*lmu*mu2_ - 24.*mu2_*sqr(Constants::pi) 
		      + 72.*mu2_*sqr(lmu) - 22.*mu4 + 1032.*mu4 * lmu
		      - 96.*mu4*sqr(Constants::pi) + 288.*mu4*sqr(lmu));
    VNS = - mu2_/1260.*(-6930. + 7560.*lmu + 2520.*mu_ - 16695.*mu2_ 
		       + 1260.*mu2_*sqr(Constants::pi) 
		       + 12600.*lmu*mu2_ + 1344.*mu_*mu2_ - 52780.*mu4 + 36960.*mu4*lmu 
		       + 5040.*mu4*sqr(Constants::pi) - 12216.*mu_*mu4);
    f2 = CF_*aS_*mu2_/Constants::pi*( 2.*lmu + 4.*mu2_*lmu + 2.*mu2_ + 12.*mu4*lmu + 7.*mu4);
  }
  // add up bits for f1
  f1 += CF_*aS_/Constants::pi*(fNS+VNS);
  double realFact(0.); 
  for(int iemit=0;iemit<2;++iemit) {
    // now for the real correction
    double phi  = UseRandom::rnd()*Constants::twopi;
    // calculate y
    double yminus = 0.; 
    double yplus  = 1.-2.*mu_*(1.-mu_)/(1.-2*mu2_);
    double y = yminus + UseRandom::rnd()*(yplus-yminus);
    // calculate z
    double v1  = sqrt(sqr(2.*mu2_+(1.-2.*mu2_)*(1.-y))-4.*mu2_)/(1.-2.*mu2_)/(1.-y);
    double zplus  = (1.+v1)*(1.-2.*mu2_)*y/2./(mu2_ +(1.-2.*mu2_)*y);
    double zminus = (1.-v1)*(1.-2.*mu2_)*y/2./(mu2_ +(1.-2.*mu2_)*y);
    double z = zminus + UseRandom::rnd()*(zplus-zminus);
    double jac = (1.-y)*(yplus-yminus)*(zplus-zminus);
    // calculate x1,x2
    double x2 = 1. - y*(1.-2.*mu2_);
    double x1 = 1. - z*(x2-2.*mu2_);
    // copy the particle objects over for calculateRealEmission
    vector<PPtr> hardProcess(3);
    hardProcess[0] = const_ptr_cast<PPtr>(&part);
    hardProcess[1] = decay[0];
    hardProcess[2] = decay[1];
    // total real emission contribution
    realFact += 0.25*jac*sqr(1.-2.*mu2_)/
    sqrt(1.-4.*mu2_)/Constants::twopi
      *2.*CF_*aS_*calculateRealEmission(x1, x2, hardProcess, phi,
					iemit, true);
  }
  // the born + virtual + real
  output = output*(1. + f1 + realFact) + f2*total;
  return output;
}

void SMZDecayer::doinitrun() {
  PerturbativeDecayer::doinitrun();
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
  // parameters for the PerturbativeDecayer base class
  PerturbativeDecayer::dataBaseOutput(output,false);
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
	Complex loamp  = FFZVertex_->evaluate(scale,wf_old[ia],
					      wfb_old[ifm],vec2[vhel]);
	Complex lotemp = FFZVertex_->evaluate(scale,wf[ia],
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
	      FFPVertex_->evaluateSmall(ZERO,3,children[iferm]->dataPtr()->CC(),
					wfb[ifm],photon[2*phel],
					ifm,2*phel,ctheta,phi,stheta,false);
	    diff[0] = FFZVertex_->evaluate(scale,wf[ia],foff,vec1[vhel]) +
	      e*double(children[iferm]->dataPtr()->iCharge())/3.*
	      UnitRemoval::E*(lotemp-loamp)*
	      (children[iferm]->momentum()*photon[2*phel].wave())/
	      (children[iferm]->momentum()*children[2]->momentum());
	    sum [0] = diff[0]+2.*dipole;
	  }
	  // special if fermion backwards
	  else {
	    SpinorBarWaveFunction foff = 
	      FFPVertex_->evaluate(ZERO,3,children[iferm]->dataPtr()->CC(),
				   wfb[ifm],photon[2*phel]);
	    Complex diag = 
	      FFZVertex_->evaluate(scale,wf[ia],foff,vec1[vhel]);
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
	      FFPVertex_->evaluateSmall(ZERO,3,children[ianti]->dataPtr()->CC(),
					wf[ia],photon[2*phel],
					ia,2*phel,ctheta,phi,stheta,false);
	    diff[1] = FFZVertex_->evaluate(scale,foff ,wfb[ifm],vec1[vhel]) +
	      e*double(children[ianti]->dataPtr()->iCharge())/3.*
	      UnitRemoval::E*(lotemp-loamp)*
	      (children[ianti]->momentum()*photon[2*phel].wave())/
	      (children[ianti]->momentum()*children[2]->momentum());
	    sum [1] = diff[1]+2.*dipole;
	  }	    
	  // special if fermion backwards after radiation
	  else {
	    SpinorWaveFunction foff = 
	      FFPVertex_->evaluate(ZERO,3,children[ianti]->dataPtr()->CC(),
				   wf[ia],photon[2*phel]);
	    Complex diag = 
	      FFZVertex_->evaluate(scale,foff ,wfb[ifm],vec1[vhel]);
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
//   double gL = (FFZVertex_->left() *FFZVertex_->norm()).real();
//   double gR = (FFZVertex_->right()*FFZVertex_->norm()).real();
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
  double gL = (FFZVertex_->left() *FFZVertex_->norm()).real();
  double gR = (FFZVertex_->right()*FFZVertex_->norm()).real();
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
initializeMECorrection(RealEmissionProcessPtr born, double & initial,
		       double & final) {
  // get the quark and antiquark
  ParticleVector qq; 
  for(unsigned int ix=0;ix<born->bornOutgoing().size();++ix)
    qq.push_back(born->bornOutgoing()[ix]);
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

bool SMZDecayer::softMatrixElementVeto(PPtr parent,
				       PPtr progenitor,
				       const bool & ,
				       const Energy & highestpT,
				       const vector<tcPDPtr> & ids,
				       const double & d_z,
				       const Energy & d_qt,
				       const Energy &) {
  // check we should be applying the veto
  if(parent->id()!=progenitor->id()||
     ids[0]->id()!=ids[1]->id()||
     ids[2]->id()!=ParticleID::g) return false;
  // calculate pt
  Energy2 d_m2 = parent->momentum().m2();
  Energy pPerp = (1.-d_z)*sqrt( sqr(d_z*d_qt) - d_m2);
  // if not hardest so far don't apply veto
  if(pPerp<highestpT) return false;
  // calculate the weight
  double weight = 0.;
  if(parent->id()>0) weight = qWeightX(d_qt, d_z);
  else weight = qbarWeightX(d_qt, d_z);
  // compute veto from weight and return 
  return !UseRandom::rndbool(weight);
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

double SMZDecayer::matrixElementRatio(const Particle & inpart, const ParticleVector & decay2,
				      const ParticleVector & decay3, MEOption,
				      ShowerInteraction inter) {
  // extract partons and LO momentas
  vector<cPDPtr> partons(1,inpart.dataPtr());
  vector<Lorentz5Momentum> lomom(1,inpart.momentum());
  for(unsigned int ix=0;ix<2;++ix) {
    partons.push_back(decay2[ix]->dataPtr());
    lomom.push_back(decay2[ix]->momentum());
  }
  vector<Lorentz5Momentum> realmom(1,inpart.momentum());
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==2) partons.push_back(decay3[ix]->dataPtr());
    realmom.push_back(decay3[ix]->momentum());
  }
  if(partons[0]->id()<0) {
    swap(partons[1],partons[2]);
    swap(lomom[1],lomom[2]);
    swap(realmom[1],realmom[2]);
  }
  scale_ = sqr(inpart.mass());
  double     lome = loME(partons,lomom);
  InvEnergy2 reme = realME(partons,realmom,inter);
  double ratio = reme/lome*sqr(inpart.mass())*4.*Constants::pi;
  if(inter==ShowerInteraction::QCD) ratio *= CF_;
  return ratio;
}

double SMZDecayer::meRatio(vector<cPDPtr> partons, 
					 vector<Lorentz5Momentum> momenta,
					 unsigned int iemitter, bool subtract) const {
  Lorentz5Momentum q = momenta[1]+momenta[2]+momenta[3];
  Energy2 Q2=q.m2();
  Energy2 lambda = sqrt((Q2-sqr(momenta[1].mass()+momenta[2].mass()))*
			(Q2-sqr(momenta[1].mass()-momenta[2].mass())));
  InvEnergy2 D[2];
  double lome[2];
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
    Lorentz5Momentum pijt = 
      q-pkt;
    double muj = momenta[1+iemit ].mass()/sqrt(Q2);
    double muk = momenta[1+ispect].mass()/sqrt(Q2);
    double vt = sqrt((1.-sqr(muj+muk))*(1.-sqr(muj-muk)))/(1.-sqr(muj)-sqr(muk));
    double v  = sqrt(sqr(2.*sqr(muk)+(1.-sqr(muj)-sqr(muk))*(1.-y))-4.*sqr(muk))
      /(1.-y)/(1.-sqr(muj)-sqr(muk));
    // dipole term
    D[iemit] = 0.5/pipj*(2./(1.-(1.-z)*(1.-y))
			 -vt/v*(2.-z+sqr(momenta[1+iemit].mass())/pipj));
    // matrix element
    vector<Lorentz5Momentum> lomom(3);
    lomom[0] = momenta[0];
    if(iemit==0) {
      lomom[1] = pijt;
      lomom[2] = pkt ;
    }
    else {
      lomom[2] = pijt;
      lomom[1] = pkt ;
    }
    lome[iemit]  = loME(partons,lomom);
  }
  InvEnergy2 ratio = realME(partons,momenta,ShowerInteraction::QCD)*abs(D[iemitter])
    /(abs(D[0]*lome[0])+abs(D[1]*lome[1]));
  if(subtract)
    return Q2*(ratio-2.*D[iemitter]);
  else
    return Q2*ratio;
}

double SMZDecayer::loME(const vector<cPDPtr> & partons, 
			const vector<Lorentz5Momentum> & momenta) const {
  // compute the spinors
  vector<VectorWaveFunction>    vin;
  vector<SpinorWaveFunction>    aout;
  vector<SpinorBarWaveFunction> fout;
  VectorWaveFunction    zin  (momenta[0],partons[0],incoming);
  SpinorBarWaveFunction qkout(momenta[1],partons[1],outgoing);
  SpinorWaveFunction    qbout(momenta[2],partons[2],outgoing);
  for(unsigned int ix=0;ix<2;++ix){
    qkout.reset(ix);
    fout.push_back(qkout);
    qbout.reset(ix);
    aout.push_back(qbout);
  }
  for(unsigned int ix=0;ix<3;++ix){
    zin.reset(ix);
    vin.push_back(zin);
  }
  // temporary storage of the different diagrams
  // sum over helicities to get the matrix element
  double total(0.);
  for(unsigned int inhel=0;inhel<3;++inhel) {
    for(unsigned int outhel1=0;outhel1<2;++outhel1) {
      for(unsigned int outhel2=0;outhel2<2;++outhel2) {
	Complex diag1 = FFZVertex_->evaluate(scale_,aout[outhel2],fout[outhel1],vin[inhel]);
	total += norm(diag1);
      }
    }
  }
  // return the answer
  return total;
}
 
InvEnergy2 SMZDecayer::realME(const vector<cPDPtr> & partons, 
			      const vector<Lorentz5Momentum> & momenta,
			      ShowerInteraction inter) const {
  AbstractFFVVertexPtr vertex = inter==ShowerInteraction::QCD ?
    FFGVertex_ : FFPVertex_;
  // compute the spinors
  vector<VectorWaveFunction>     vin;
  vector<SpinorWaveFunction>     aout;
  vector<SpinorBarWaveFunction>  fout;
  vector<VectorWaveFunction>     gout;
  VectorWaveFunction    zin  (momenta[0],partons[0],incoming);
  SpinorBarWaveFunction qkout(momenta[1],partons[1],outgoing);
  SpinorWaveFunction    qbout(momenta[2],partons[2],outgoing);
  VectorWaveFunction    gluon(momenta[3],partons[3],outgoing);
  for(unsigned int ix=0;ix<2;++ix){
    qkout.reset(ix);
    fout.push_back(qkout);
    qbout.reset(ix);
    aout.push_back(qbout);
    gluon.reset(2*ix);
    gout.push_back(gluon);
  }
  for(unsigned int ix=0;ix<3;++ix){
    zin.reset(ix);
    vin.push_back(zin);
  }
  vector<Complex> diag(2,0.);

  double total(0.);
  for(unsigned int inhel1=0;inhel1<3;++inhel1) {
    for(unsigned int outhel1=0;outhel1<2;++outhel1) {
      for(unsigned int outhel2=0;outhel2<2;++outhel2) {
	for(unsigned int outhel3=0;outhel3<2;++outhel3) {
	  SpinorBarWaveFunction off1 =
	    vertex->evaluate(scale_,3,partons[1]->CC(),fout[outhel1],gout[outhel3]);
	  diag[0] = FFZVertex_->evaluate(scale_,aout[outhel2],off1,vin[inhel1]);

	  SpinorWaveFunction off2 = 
	    vertex->evaluate(scale_,3,partons[2]->CC(),aout[outhel2],gout[outhel3]);
	  diag[1] = FFZVertex_->evaluate(scale_,off2,fout[outhel1],vin[inhel1]);

	  // sum of diagrams
	  Complex sum = std::accumulate(diag.begin(),diag.end(),Complex(0.));
	  // me2
	  total += norm(sum);
	}
      }
    }
  }

  // divide out the coupling
  total /= norm(vertex->norm());
  // return the total
  return total*UnitRemoval::InvE2;
}

double SMZDecayer::calculateRealEmission(double x1, double x2, 
					 vector<PPtr> hardProcess,
					 double phi,
					 bool subtract) const {
  // make partons data object for meRatio
  vector<cPDPtr> partons (3);
  for(int ix=0; ix<3; ++ix)
    partons[ix] = hardProcess[ix]->dataPtr();
  partons.push_back(gluon_);
  // calculate x3
  double x3 = 2.-x1-x2;
  double xT = sqrt(max(0.,sqr(x3) -0.25*sqr(sqr(x2)+sqr(x3)-sqr(x1))/(sqr(x2)-4.*mu2_)));
  // calculate the momenta
  Energy M = mZ_;
  Lorentz5Momentum pspect(ZERO,ZERO,-0.5*M*sqrt(max(sqr(x2)-4.*mu2_,0.)),0.5*M*x2,M*mu_); 
  Lorentz5Momentum pemit (-0.5*M*xT*cos(phi),-0.5*M*xT*sin(phi),
  			  0.5*M*sqrt(max(sqr(x1)-sqr(xT)-4.*mu2_,0.)),0.5*M*x1,M*mu_);
  Lorentz5Momentum pgluon(0.5*M*xT*cos(phi), 0.5*M*xT*sin(phi),
			  0.5*M*sqrt(max(sqr(x3)-sqr(xT),0.)),0.5*M*x3,ZERO);
  if(abs(pspect.z()+pemit.z()-pgluon.z())/M<1e-6) 
    pgluon.setZ(-pgluon.z());
  else if(abs(pspect.z()-pemit.z()+pgluon.z())/M<1e-6) 
    pemit .setZ(- pemit.z());
  // loop over the possible emitting partons
  double realwgt(0.);
  for(unsigned int iemit=0;iemit<2;++iemit) {
    // boost and rotate momenta
    LorentzRotation eventFrame( ( hardProcess[1]->momentum() +
				  hardProcess[2]->momentum() ).findBoostToCM() );
    Lorentz5Momentum spectator = eventFrame*hardProcess[iemit+1]->momentum();
    eventFrame.rotateZ( -spectator.phi()    );
    eventFrame.rotateY( -spectator.theta()  );
    eventFrame.invert();
    vector<Lorentz5Momentum> momenta(3);
    momenta[0]   = hardProcess[0]->momentum();
    if(iemit==0) {
      momenta[2] = eventFrame*pspect;
      momenta[1] = eventFrame*pemit ;
    }
    else {
      momenta[1] = eventFrame*pspect;
      momenta[2] = eventFrame*pemit ;
    }
    momenta.push_back(eventFrame*pgluon);
    // calculate the weight
    if(1.-x1>1e-5 && 1.-x2>1e-5) 
      realwgt += meRatio(partons,momenta,iemit,subtract);
  }
  
  // total real emission contribution
  return realwgt;
}

double SMZDecayer::calculateRealEmission(double x1, double x2, 
					 vector<PPtr> hardProcess,
					 double phi,
					 bool subtract,
					 int emitter) const {
  // make partons data object for meRatio
  vector<cPDPtr> partons (3);
  for(int ix=0; ix<3; ++ix)
    partons[ix] = hardProcess[ix]->dataPtr();
  partons.push_back(gluon_);
  // calculate x3
  double x3 = 2.-x1-x2;
  double xT = sqrt(max(0.,sqr(x3) -0.25*sqr(sqr(x2)+sqr(x3)-sqr(x1))/(sqr(x2)-4.*mu2_)));
  // calculate the momenta
  Energy M = mZ_;
  Lorentz5Momentum pspect(ZERO,ZERO,-0.5*M*sqrt(max(sqr(x2)-4.*mu2_,0.)),0.5*M*x2,M*mu_); 
  Lorentz5Momentum pemit (-0.5*M*xT*cos(phi),-0.5*M*xT*sin(phi),
  			  0.5*M*sqrt(max(sqr(x1)-sqr(xT)-4.*mu2_,0.)),0.5*M*x1,M*mu_);
  Lorentz5Momentum pgluon( 0.5*M*xT*cos(phi), 0.5*M*xT*sin(phi),
  			   0.5*M*sqrt(max(sqr(x3)-sqr(xT),0.)),0.5*M*x3,ZERO);
  if(abs(pspect.z()+pemit.z()-pgluon.z())/M<1e-6) 
    pgluon.setZ(-pgluon.z());
  else if(abs(pspect.z()-pemit.z()+pgluon.z())/M<1e-6) 
    pemit .setZ(- pemit.z());
  // boost and rotate momenta
  LorentzRotation eventFrame( ( hardProcess[1]->momentum() +
				hardProcess[2]->momentum() ).findBoostToCM() );
  Lorentz5Momentum spectator = eventFrame*hardProcess[emitter+1]->momentum();
  eventFrame.rotateZ( -spectator.phi()    );
  eventFrame.rotateY( -spectator.theta()  );
  eventFrame.invert();
  vector<Lorentz5Momentum> momenta(3);
  momenta[0]   = hardProcess[0]->momentum();
  if(emitter==0) {
    momenta[2] = eventFrame*pspect;
    momenta[1] = eventFrame*pemit ;
  }
  else {
    momenta[1] = eventFrame*pspect;
    momenta[2] = eventFrame*pemit ;
  }
  momenta.push_back(eventFrame*pgluon);
  // calculate the weight
  double realwgt(0.);
  if(1.-x1>1e-5 && 1.-x2>1e-5) 
    realwgt = meRatio(partons,momenta,emitter,subtract);  
  // total real emission contribution
  return realwgt;
}
