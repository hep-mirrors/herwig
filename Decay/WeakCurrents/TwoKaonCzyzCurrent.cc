// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TwoKaonCzyzCurrent class.
//

#include "TwoKaonCzyzCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Decay/ResonanceHelpers.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include <cmath>

using namespace Herwig;

HERWIG_INTERPOLATOR_CLASSDESC(TwoKaonCzyzCurrent,double,Energy2)

TwoKaonCzyzCurrent::TwoKaonCzyzCurrent()
// various parameters from 1002.0279, Fit 2 commented out                                                                                         
// substituted by own fit like Fit 2 but with new data + phi mass and width kept fixed
//: betaRho_(2.21), betaOmega_(2.75), betaPhi_(1.91),
//  nMax_(1000), etaPhi_(1.055), gammaOmega_(0.5), gammaPhi_(0.2), mpi_(140.*MeV) {
: betaRho_(2.20), betaOmega_(2.52), betaPhi_(1.96),
  nMax_(200), etaPhi_(1.03), gammaOmega_(0.5), gammaPhi_(0.2), mpi_(140.*MeV) {
  using Constants::pi;
  // rho parameter
  //rhoMag_    = {1.120, 0.107, 0.028,0.032};
  //rhoPhase_  = {0    ,    pi,    pi,    0};
  rhoMag_    = {1.120, 0.101, 0.0347,0.0803};//,0.0692,0.15488};
  rhoPhase_  = {0    ,    pi,    pi,    0.};//, pi,0.};
  rhoMasses_ = {775.49*MeV,1465.*MeV,1720.*MeV};//,2150*MeV};
  rhoWidths_ = {149.4 *MeV,400. *MeV, 250.*MeV};//,150*MeV};
  // omega parameters
  //omegaMag_    = {1.37, 0.173, 0.621,0.43};
  //omegaPhase_  = {0   ,    pi,    pi,   0};
  omegaMag_    = {1.28, 0.0267, 0.742,0.908};//,1.4831};
  omegaPhase_  = {0   ,    0,    pi,   0};//,pi};
  omegaMasses_ = {782.65*MeV,1425.*MeV,1670.*MeV};
  omegaWidths_ = {8.49  *MeV, 215.*MeV, 315.*MeV};
  // phi parameters
  //phiMag_    = {0.947,0.0136,0.0214};
  //phiPhase_  = {0.   ,0.    ,0.    };
  phiMag_    = {0.976,0.0138,0.00223};//,0.0842,0.079138};
  phiPhase_  = {0.   ,0.    ,0.};//   ,0. ,0.};
  phiMasses_ = {1019.415*MeV,1680.*MeV};//,2183*MeV};
  phiWidths_ = {4.22    *MeV, 150.*MeV};//,88*MeV};
  // set up for the modes in the base class
  addDecayMode(2,-1);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  setInitialModes(5);
}

IBPtr TwoKaonCzyzCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr TwoKaonCzyzCurrent::fullclone() const {
  return new_ptr(*this);
}

void TwoKaonCzyzCurrent::persistentOutput(PersistentOStream & os) const {
  os << betaRho_ << betaOmega_ << betaPhi_
     << rhoWgt_ << rhoMag_ << rhoPhase_
     << ounit(rhoMasses_,GeV) << ounit(rhoWidths_,GeV)
     << phiWgt_ << phiMag_ << phiPhase_
     << ounit(phiMasses_,GeV) << ounit(phiWidths_,GeV)
     << ounit(mass_,GeV) << ounit(width_,GeV) << coup_
     << dh_ << ounit(hres_,GeV2) << ounit(h0_,GeV2)
     << nMax_ << etaPhi_ << gammaOmega_ << gammaPhi_ << ounit(mpi_,GeV)
     << ounit(eMax_,GeV) << fKI0Re_ << fKI0Im_ << fKI1Re_ << fKI1Im_;
}

void TwoKaonCzyzCurrent::persistentInput(PersistentIStream & is, int) {
  is >> betaRho_ >> betaOmega_ >> betaPhi_
     >> rhoWgt_ >> rhoMag_ >> rhoPhase_
     >> iunit(rhoMasses_,GeV) >> iunit(rhoWidths_,GeV)
     >> phiWgt_ >> phiMag_ >> phiPhase_
     >> iunit(phiMasses_,GeV) >> iunit(phiWidths_,GeV)
     >> iunit(mass_,GeV) >> iunit(width_,GeV) >> coup_
     >> dh_ >> iunit(hres_,GeV2) >> iunit(h0_,GeV2)
     >> nMax_ >> etaPhi_ >> gammaOmega_ >> gammaPhi_ >> iunit(mpi_,GeV)
     >> iunit(eMax_,GeV) >> fKI0Re_ >> fKI0Im_ >> fKI1Re_ >> fKI1Im_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<TwoKaonCzyzCurrent,WeakCurrent>
describeHerwigTwoKaonCzyzCurrent("Herwig::TwoKaonCzyzCurrent", "HwWeakCurrents.so");

void TwoKaonCzyzCurrent::Init() {

  static ClassDocumentation<TwoKaonCzyzCurrent> documentation
    ("The TwoKaonCzyzCurrent class uses the currents from "
     "PRD 81 094014 for the weak current with two kaons",
     "The current for two kaons from \\cite{Czyz:2010hj} was used.",
     "%\\cite{Czyz:2010hj}\n"
     "\\bibitem{Czyz:2010hj}\n"
     "H.~Czyz, A.~Grzelinska and J.~H.~Kuhn,\n"
     "%``Narrow resonances studies with the radiative return method,''\n"
     "Phys.\\ Rev.\\ D {\\bf 81} (2010) 094014\n"
     "doi:10.1103/PhysRevD.81.094014\n"
     "[arXiv:1002.0279 [hep-ph]].\n"
     "%%CITATION = doi:10.1103/PhysRevD.81.094014;%%\n"
     "%28 citations counted in INSPIRE as of 30 Jul 2018\n");

  static ParVector<TwoKaonCzyzCurrent,Energy> interfaceRhoMasses
    ("RhoMasses",
     "The masses of the different rho resonances for the pi pi channel",
     &TwoKaonCzyzCurrent::rhoMasses_, MeV, -1, 775.8*MeV, ZERO, 10000.*MeV,
     false, false, true);

  static ParVector<TwoKaonCzyzCurrent,Energy> interfaceRhoWidths
    ("RhoWidths",
     "The widths of the different rho resonances for the pi pi channel",
     &TwoKaonCzyzCurrent::rhoWidths_, MeV, -1, 150.3*MeV, ZERO, 1000.*MeV,
     false, false, true);
  
  static ParVector<TwoKaonCzyzCurrent,double> interfaceRhoMagnitude
    ("RhoMagnitude",
     "Magnitude of the weight of the different resonances for the pi pi channel",
     &TwoKaonCzyzCurrent::rhoMag_, -1, 0., 0, 0,
     false, false, Interface::nolimits);

  static ParVector<TwoKaonCzyzCurrent,double> interfaceRhoPhase
    ("RhoPhase",
     "Phase of the weight of the different resonances for the pi pi channel",
     &TwoKaonCzyzCurrent::rhoPhase_, -1, 0., 0, 0,
     false, false, Interface::nolimits);

  static ParVector<TwoKaonCzyzCurrent,Energy> interfaceOmegaMasses
    ("OmegaMasses",
     "The masses of the different omega resonances for the pi pi channel",
     &TwoKaonCzyzCurrent::omegaMasses_, MeV, -1, 775.8*MeV, ZERO, 10000.*MeV,
     false, false, true);

  static ParVector<TwoKaonCzyzCurrent,Energy> interfaceOmegaWidths
    ("OmegaWidths",
     "The widths of the different omega resonances for the pi pi channel",
     &TwoKaonCzyzCurrent::omegaWidths_, MeV, -1, 150.3*MeV, ZERO, 1000.*MeV,
     false, false, true);
  
  static ParVector<TwoKaonCzyzCurrent,double> interfaceOmegaMagnitude
    ("OmegaMagnitude",
     "Magnitude of the weight of the different resonances for the pi pi channel",
     &TwoKaonCzyzCurrent::omegaMag_, -1, 0., 0, 0,
     false, false, Interface::nolimits);

  static ParVector<TwoKaonCzyzCurrent,double> interfaceOmegaPhase
    ("OmegaPhase",
     "Phase of the weight of the different resonances for the pi pi channel",
     &TwoKaonCzyzCurrent::omegaPhase_, -1, 0., 0, 0,
     false, false, Interface::nolimits);

  static ParVector<TwoKaonCzyzCurrent,Energy> interfacePhiMasses
    ("PhiMasses",
     "The masses of the different phi resonances for the pi pi channel",
     &TwoKaonCzyzCurrent::phiMasses_, MeV, -1, 775.8*MeV, ZERO, 10000.*MeV,
     false, false, true);

  static ParVector<TwoKaonCzyzCurrent,Energy> interfacePhiWidths
    ("PhiWidths",
     "The widths of the different phi resonances for the pi pi channel",
     &TwoKaonCzyzCurrent::phiWidths_, MeV, -1, 150.3*MeV, ZERO, 1000.*MeV,
     false, false, true);
  
  static ParVector<TwoKaonCzyzCurrent,double> interfacePhiMagnitude
    ("PhiMagnitude",
     "Magnitude of the weight of the different resonances for the pi pi channel",
     &TwoKaonCzyzCurrent::phiMag_, -1, 0., 0, 0,
     false, false, Interface::nolimits);

  static ParVector<TwoKaonCzyzCurrent,double> interfacePhiPhase
    ("PhiPhase",
     "Phase of the weight of the different resonances for the pi pi channel",
     &TwoKaonCzyzCurrent::phiPhase_, -1, 0., 0, 0,
     false, false, Interface::nolimits);

  static Parameter<TwoKaonCzyzCurrent,unsigned int> interfacenMax
    ("nMax",
     "The maximum number of resonances to include in the sum,"
     " should be approx infinity",
     &TwoKaonCzyzCurrent::nMax_, 1000, 10, 10000,
     false, false, Interface::limited);

  static Parameter<TwoKaonCzyzCurrent,double> interfacebetaRho
    ("betaRho",
     "The beta parameter for the rho couplings",
     &TwoKaonCzyzCurrent::betaRho_, 2.23, 0.0, 100.,
     false, false, Interface::limited);

  static Parameter<TwoKaonCzyzCurrent,double> interfacebetaOmega
    ("betaOmega",
     "The beta parameter for the rho couplings",
     &TwoKaonCzyzCurrent::betaOmega_, 2.23, 0.0, 100.,
     false, false, Interface::limited);

  static Parameter<TwoKaonCzyzCurrent,double> interfacebetaPhi
    ("betaPhi",
     "The beta parameter for the phi couplings",
     &TwoKaonCzyzCurrent::betaPhi_, 1.97, 0.0, 100.,
     false, false, Interface::limited);

  static Parameter<TwoKaonCzyzCurrent,double> interfaceEtaPhi
    ("EtaPhi",
     "The eta_phi mixing parameter",
     &TwoKaonCzyzCurrent::etaPhi_, 1.04, 0.0, 10.0,
     false, false, Interface::limited);
  
  static Parameter<TwoKaonCzyzCurrent,double> interfacegammaOmega
    ("gammaOmega",
     "The gamma parameter for the widths of omega resonances",
     &TwoKaonCzyzCurrent::gammaOmega_, 0.5, 0.0, 1.0,
     false, false, Interface::limited);
  
  static Parameter<TwoKaonCzyzCurrent,double> interfacegammaPhi
    ("gammaPhi",
     "The gamma parameter for the widths of phi resonances",
     &TwoKaonCzyzCurrent::gammaPhi_, 0.2, 0.0, 1.0,
     false, false, Interface::limited);

}

void TwoKaonCzyzCurrent::doinit() {
  WeakCurrent::doinit();
  // check consistency of parametrers
  if(rhoMasses_.size()   != rhoWidths_.size() ||
     omegaMasses_.size() != omegaWidths_.size() ||
     phiMasses_.size()   != phiWidths_.size() )
    throw InitException() << "Inconsistent parameters in TwoKaonCzyzCurrent"
			  << "::doinit()" << Exception::abortnow;
  // weights for the rho channels
  if(rhoMag_.size()!=rhoPhase_.size()) 
    throw InitException() << "The vectors containing the weights and phase for the"
  			  << " rho channel must be the same size in "
  			  << "TwoKaonCzyzCurrent::doinit()" << Exception::runerror;
  // combine mags and phase
  for(unsigned int ix=0;ix<rhoMag_.size();++ix) {
    rhoWgt_.push_back(rhoMag_[ix]*(cos(rhoPhase_[ix])+Complex(0.,1.)*sin(rhoPhase_[ix])));
  }
  for(unsigned int ix=0;ix<omegaMag_.size();++ix) {
    omegaWgt_.push_back(omegaMag_[ix]*(cos(omegaPhase_[ix])+Complex(0.,1.)*sin(omegaPhase_[ix])));
  }
  for(unsigned int ix=0;ix<phiMag_.size();++ix) {
    phiWgt_.push_back(phiMag_[ix]*(cos(phiPhase_[ix])+Complex(0.,1.)*sin(phiPhase_[ix])));
  }
  // pion mass
  mpi_ = getParticleData(211)->mass();
  // rho masses and couplings
  double gamB(std::tgamma(2.-betaRho_));
  mass_.push_back(vector<Energy>());
  width_.push_back(vector<Energy>());
  coup_.push_back(vector<Complex>());
  for(unsigned int ix=0;ix<nMax_;++ix) {
    // this is gam(2-beta+n)/gam(n+1)
    if(ix>0) {
      gamB *= ((1.-betaRho_+double(ix)))/double(ix);
    }
    Complex c_n = std::tgamma(betaRho_-0.5) /(0.5+double(ix)) / sqrt(Constants::pi) *
      sin(Constants::pi*(betaRho_-1.-double(ix)))/Constants::pi*gamB;
    if(ix%2!=0) c_n *= -1.;
    // couplings
    if(ix>=rhoWgt_.size()) {
      coup_[0].push_back(c_n);
    }
    else {
      coup_[0].push_back(rhoWgt_[ix]);
    }
    // set the masses and widths
    // calc for higher resonances
    if(ix>=rhoMasses_.size()) {
      mass_ [0].push_back(rhoMasses_[0]*sqrt(1.+2.*double(ix)));
      width_[0].push_back(rhoWidths_[0]/rhoMasses_[0]*mass_[0].back());
    }
    // input for lower ones
    else {
      mass_ [0].push_back(rhoMasses_[ix]);
      width_[0].push_back(rhoWidths_[ix]);
    }
    // parameters for the gs propagators
    hres_.push_back(Resonance::Hhat(sqr(mass_[0].back()),
	 			       mass_[0].back(),width_[0].back(),mpi_,mpi_));
    dh_  .push_back(Resonance::dHhatds(mass_[0].back(),width_[0].back(),mpi_,mpi_));
    h0_  .push_back(Resonance::H(ZERO,mass_[0].back(),width_[0].back(),
				 mpi_,mpi_,dh_.back(),hres_.back()));
  }
  // omega masses and couplings
  gamB = std::tgamma(2.-betaOmega_);
  mass_.push_back(vector<Energy>());
  width_.push_back(vector<Energy>());
  coup_.push_back(vector<Complex>());
  for(unsigned int ix=0;ix<nMax_;++ix) {
    // this is gam(2-beta+n)/gam(n+1)
    if(ix>0) {
      gamB *= ((1.-betaOmega_+double(ix)))/double(ix);
    }
    Complex c_n = std::tgamma(betaOmega_-0.5) /(0.5+double(ix)) / sqrt(Constants::pi) *
      sin(Constants::pi*(betaOmega_-1.-double(ix)))/Constants::pi*gamB;
    if(ix%2!=0) c_n *= -1.;
    // couplings
    if(ix>=omegaWgt_.size()) {
      coup_[1].push_back(c_n);
    }
    else {
      coup_[1].push_back(omegaWgt_[ix]);
    }
    // set the masses and widths
    // calc for higher resonances
    if(ix>=omegaMasses_.size()) {
      mass_ [1].push_back(omegaMasses_[0]*sqrt(1.+2.*double(ix)));
      width_[1].push_back(gammaOmega_*mass_[1].back());
    }
    // input for lower ones
    else {
      mass_ [1].push_back(omegaMasses_[ix]);
      width_[1].push_back(omegaWidths_[ix]);
    }
  }
  // phi masses and couplings
  gamB = std::tgamma(2.-betaPhi_);
  mass_.push_back(vector<Energy>());
  width_.push_back(vector<Energy>());
  coup_.push_back(vector<Complex>());
  for(unsigned int ix=0;ix<nMax_;++ix) {
    // this is gam(2-beta+n)/gam(n+1)
    if(ix>0) {
      gamB *= ((1.-betaPhi_+double(ix)))/double(ix);
    }
    Complex c_n = std::tgamma(betaPhi_-0.5) /(0.5+double(ix)) / sqrt(Constants::pi) *
      sin(Constants::pi*(betaPhi_-1.-double(ix)))/Constants::pi*gamB;
    if(ix%2!=0) c_n *= -1.;
    // couplings
    if(ix>=phiWgt_.size()) {
      coup_[2].push_back(c_n);
    }
    else {
      coup_[2].push_back(phiWgt_[ix]);
    }
    // set the masses and widths
    // calc for higher resonances
    if(ix>=phiMasses_.size()) {
      mass_ [2].push_back(phiMasses_[0]*sqrt(1.+2.*double(ix)));
      width_[2].push_back(gammaPhi_*mass_[2].back());
    }
    // input for lower ones
    else {
      mass_ [2].push_back(phiMasses_[ix]);
      width_[2].push_back(phiWidths_[ix]);
    }
  }
}

void TwoKaonCzyzCurrent::constructInterpolators() const {
  // construct the interpolators
  vector<Energy2> en;
  vector<double> re0,im0;
  vector<double> re1,im1;
  Energy mK = getParticleData(ParticleID::Kplus)->mass();
  Energy2 step = (sqr(eMax_)-sqr(2.*mK))/nMax_;
  Energy2 Q2 = sqr(2.*mK);
  for(unsigned int ix=0;ix<nMax_+1;++ix) {
    Complex value = FkaonRemainderI1(Q2);
    re1.push_back(value.real());
    im1.push_back(value.imag());
    value = FkaonRemainderI0(Q2,mK,mK);
    re0.push_back(value.real());
    im0.push_back(value.imag());
    en.push_back(Q2);
    Q2+=step;
  }
  fKI0Re_ = make_InterpolatorPtr(re0,en,3);
  fKI0Im_ = make_InterpolatorPtr(im0,en,3);
  fKI1Re_ = make_InterpolatorPtr(re1,en,3);
  fKI1Im_ = make_InterpolatorPtr(im1,en,3);
}

// complete the construction of the decay mode for integration
bool TwoKaonCzyzCurrent::createMode(int icharge, tcPDPtr resonance,
				    IsoSpin::IsoSpin Itotal, IsoSpin::I3 i3, Strangeness::Strange S,
				    unsigned int imode,PhaseSpaceModePtr mode,
				    unsigned int iloc,int ires,
				    PhaseSpaceChannel phase, Energy upp ) {
  // check the charge
  if((imode==0 && abs(icharge)!=3) ||
     (imode>0  && icharge !=0)) return false; 
  // check the total isospin
  if(Itotal!=IsoSpin::IUnknown) {
    if(Itotal!=IsoSpin::IOne && Itotal!=IsoSpin::IZero ) return false;
  }
  // check I_3
  if(i3!=IsoSpin::I3Unknown&&Itotal==IsoSpin::IOne) {
    switch(i3) {
    case IsoSpin::I3Zero:
      if(imode==0) return false;
      break;
    case IsoSpin::I3One:
      if(imode!=0 || icharge ==-3) return false;
      break;
    case IsoSpin::I3MinusOne:
      if(imode!=0 || icharge ==3) return false;
      break;
    default:
      return false;
    }
  }
  // make sure that the decays are kinematically allowed
  tPDPtr part[2];
  if(imode==0) {
    part[0]=getParticleData(ParticleID::Kplus);
    part[1]=getParticleData(ParticleID::Kbar0);
  }
  else if(imode==1|| imode==2)  {
    part[0]=getParticleData(ParticleID::K_S0);
    part[1]=getParticleData(ParticleID::K_L0);
  }
  else {
    part[0]=getParticleData(ParticleID::Kplus);
    part[1]=getParticleData(ParticleID::Kminus);
  }
  Energy min(part[0]->massMin()+part[1]->massMin());
  if(min>upp) return false;
  eMax_=upp;
  // set the resonances
  vector<tPDPtr> res;
  if(icharge==0) {
    res.push_back(getParticleData(113   ));
    res.push_back(getParticleData(100113));
    res.push_back(getParticleData(30113 ));
    res.push_back(getParticleData(   223));
    res.push_back(getParticleData(   333));
  }
  else {
    res.push_back(getParticleData(213   ));
    res.push_back(getParticleData(100213));
    res.push_back(getParticleData(30213 ));
    if(icharge==-3) {
      for(unsigned int ix=0;ix<3;++ix) {
  	if(res[ix]&&res[ix]->CC()) res[ix]=res[ix]->CC();
      }
    }
  }
  // create the channels
  for(unsigned int ix=0;ix<res.size();++ix) {
    if(!res[ix]) continue;
    if(resonance && resonance != res[ix]) continue;
    if(Itotal!=IsoSpin::IUnknown && Itotal!=IsoSpin::IOne  && ix < 3) continue;
    if(Itotal!=IsoSpin::IUnknown && Itotal!=IsoSpin::IZero && ix >=3) continue;
    PhaseSpaceChannel newChannel((PhaseSpaceChannel(phase),ires,res[ix],ires+1,iloc+1,ires+1,iloc+2));
    mode->addChannel(newChannel);
  }
  // reset the masses in the intergrators
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix<rhoMasses_.size()&&res[ix]) {
      mode->resetIntermediate(res[ix],rhoMasses_[ix],rhoWidths_[ix]);
    }
  }
  if(res.size()>3) {
    mode->resetIntermediate(res[3],omegaMasses_[0],omegaWidths_[0]);
    mode->resetIntermediate(res[4],phiMasses_  [0],  phiWidths_[0]);
  }
  // return if successful
  return true;
}

// the particles produced by the current
tPDVector TwoKaonCzyzCurrent::particles(int icharge, unsigned int imode,
					     int,int) {
  tPDVector output(2);
  if(imode==0) {
    output[0]=getParticleData(ParticleID::Kplus);
    output[1]=getParticleData(ParticleID::K0);
  }
  else if(imode==1||imode==2) {
    output[0]=getParticleData(ParticleID::K_S0);
    output[1]=getParticleData(ParticleID::K_L0);
  }
  else {
    output[0]=getParticleData(ParticleID::Kplus );
    output[1]=getParticleData(ParticleID::Kminus);
  }
  if(icharge==-3) {
    for(unsigned int ix=0;ix<output.size();++ix) {
      if(output[ix]->CC()) output[ix]=output[ix]->CC();
    }
  }
  return output;
}

// hadronic current   
vector<LorentzPolarizationVectorE> 
TwoKaonCzyzCurrent::current(tcPDPtr resonance,
			    IsoSpin::IsoSpin Itotal, IsoSpin::I3 i3, Strangeness::Strange S,
			    const int imode, const int ichan,Energy & scale, 
			    const tPDVector & outgoing,
			    const vector<Lorentz5Momentum> & momenta,
			    DecayIntegrator::MEOption) const {
  useMe();
  // check the total isospin
  if(Itotal!=IsoSpin::IUnknown) {
    if(Itotal!=IsoSpin::IOne && Itotal!=IsoSpin::IZero )
      return vector<LorentzPolarizationVectorE>();
  }
  // check I_3
  int icharge = outgoing[0]->iCharge()+outgoing[1]->iCharge();
  if(i3!=IsoSpin::I3Unknown&&Itotal==IsoSpin::IOne) {
    switch(i3) {
    case IsoSpin::I3Zero:
      if(imode==0) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3One:
      if(imode!=0 || icharge ==-3) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3MinusOne:
      if(imode!=0 || icharge ==3) return vector<LorentzPolarizationVectorE>();
      break;
    default:
      return vector<LorentzPolarizationVectorE>();
    }
  }
  // momentum difference and sum of the mesons
  Lorentz5Momentum pdiff(momenta[0]-momenta[1]);
  Lorentz5Momentum psum (momenta[0]+momenta[1]);
  psum.rescaleMass();
  scale=psum.mass();
  // mass2 of vector intermediate state
  Energy2 q2(psum.m2());
  double dot(psum*pdiff/q2);
  psum *=dot;
  // calculate the current
  Complex FK = Fkaon(q2,imode,ichan,Itotal,resonance,
		     momenta[0].mass(),momenta[1].mass());
  // compute the current
  pdiff -= psum;
  return vector<LorentzPolarizationVectorE>(1,FK*pdiff);
}
   
bool TwoKaonCzyzCurrent::accept(vector<int> id) {
  // check there are only two particles
  if(id.size()!=2) return false;
  // pion modes
  if((id[0]==ParticleID::Kminus && id[1]==ParticleID::K0)     ||
     (id[0]==ParticleID::K0     && id[1]==ParticleID::Kminus) ||
     (id[0]==ParticleID::Kplus  && id[1]==ParticleID::Kbar0)  ||
     (id[0]==ParticleID::Kbar0  && id[1]==ParticleID::Kplus))
    return true;
  else if((id[0]==ParticleID::Kminus && id[1]==ParticleID::Kplus) ||
	  (id[0]==ParticleID::Kplus  && id[1]==ParticleID::Kminus))
    return true;
  else if((id[0]==ParticleID::K_S0 && id[1]==ParticleID::K_L0) ||
	  (id[0]==ParticleID::K_L0 && id[1]==ParticleID::K_S0))
    return true;
  else
    return false;
}

// the decay mode
unsigned int TwoKaonCzyzCurrent::decayMode(vector<int> idout) {
  unsigned int nk0(0),nkp(0);
  for(unsigned int ix=0;ix<idout.size();++ix) {
    if(abs(idout[ix])==ParticleID::Kplus) ++nkp;
    else if(abs(idout[ix])==ParticleID::K0 ||
	    idout[ix]==ParticleID::K_L0 ||idout[ix]==ParticleID::K_S0 ) ++nk0;
  }
  if(nkp==1&&nk0==1) return 0;
  else if(nkp==2)    return 3;
  else if(nk0==2)    return 1;
  else return false;
}

// output the information for the database
void TwoKaonCzyzCurrent::dataBaseOutput(ofstream & output,bool header,
					bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::TwoKaonCzyzCurrent " 
		    << name() << " HwWeakCurrents.so\n";
  unsigned int ix;
  for(ix=0;ix<rhoMasses_.size();++ix) {
    if(ix<3)  output << "newdef ";
    else      output << "insert ";
    output << name() << ":RhoMasses " << ix << " " << rhoMasses_[ix]/MeV << "\n";
  }
  for(ix=0;ix<rhoWidths_.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoWidths " << ix << " " << rhoWidths_[ix]/MeV << "\n";
  }
  for(ix=0;ix<rhoWgt_.size();++ix) {
    if(ix<5) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoMagnitude " << ix << " " << rhoMag_[ix]   << "\n";
    if(ix<5) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoPhase "     << ix << " " << rhoPhase_[ix] << "\n";
  }
  for(ix=0;ix<omegaMasses_.size();++ix) {
    if(ix<3)  output << "newdef ";
    else      output << "insert ";
    output << name() << ":OmegaMasses " << ix << " " << omegaMasses_[ix]/MeV << "\n";
  }
  for(ix=0;ix<omegaWidths_.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":OmegaWidths " << ix << " " << omegaWidths_[ix]/MeV << "\n";
  }
  for(ix=0;ix<omegaWgt_.size();++ix) {
    if(ix<5) output << "newdef ";
    else     output << "insert ";
    output << name() << ":OmegaMagnitude " << ix << " " << omegaMag_[ix]   << "\n";
    if(ix<5) output << "newdef ";
    else     output << "insert ";
    output << name() << ":OmegaPhase "     << ix << " " << omegaPhase_[ix] << "\n";
  }
  for(ix=0;ix<phiMasses_.size();++ix) {
    if(ix<2)  output << "newdef ";
    else      output << "insert ";
    output << name() << ":PhiMasses " << ix << " " << phiMasses_[ix]/MeV << "\n";
  }
  for(ix=0;ix<phiWidths_.size();++ix) {
    if(ix<2) output << "newdef ";
    else     output << "insert ";
    output << name() << ":PhiWidths " << ix << " " << phiWidths_[ix]/MeV << "\n";
  }
  for(ix=0;ix<phiWgt_.size();++ix) {
    if(ix<4) output << "newdef ";
    else     output << "insert ";
    output << name() << ":PhiMagnitude " << ix << " " << phiMag_[ix]   << "\n";
    if(ix<4) output << "newdef ";
    else     output << "insert ";
    output << name() << ":PhiPhase "     << ix << " " << phiPhase_[ix] << "\n";
  }
  output << "newdef " << name() << ":betaRho " << betaRho_ << "\n";
  output << "newdef " << name() << ":betaOmega " << betaOmega_ << "\n";
  output << "newdef " << name() << ":betaPhi " << betaPhi_ << "\n";
  output << "newdef " << name() << ":gammaOmega " << gammaOmega_ << "\n";
  output << "newdef " << name() << ":gammaPhi " << gammaPhi_ << "\n";
  output << "newdef " << name() << ":etaPhi " << etaPhi_ << "\n";
  output << "newdef " << name() << ":nMax " << nMax_ << "\n";
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

Complex TwoKaonCzyzCurrent::Fkaon(Energy2 q2,const int imode, const int ichan,
				  IsoSpin::IsoSpin Itotal, tcPDPtr resonance,
				  Energy ma, Energy mb) const {
  unsigned int imin=0, imax = 4;
  bool on[3] = {(Itotal==IsoSpin::IUnknown || Itotal==IsoSpin::IOne),
		(Itotal==IsoSpin::IUnknown || Itotal==IsoSpin::IZero) && imode!=0,
		(Itotal==IsoSpin::IUnknown || Itotal==IsoSpin::IZero) && imode!=0};
  if(ichan>=0) {
    if(ichan<3) {
      on[1]=on[2]=false;
      imin = ichan;
      imax = ichan+1;
    }
    else if(ichan==3) {
      on[0]=on[2]=false;
      imin=0;
      imax=1;
    }
    else if(ichan==4) {
      on[0]=on[1]=false;
      imin=0;
      imax=1;
    }
    else
      assert(false);
  }
  if(resonance) {
    switch(resonance->id()%1000) {
    case 223:
      imin=0;
      on[0]=on[2]=false;
      break;
    case 333:
      imin=0;
      on[0]=on[1]=false;
      break;
    case 113:
      switch(resonance->id()/1000) {
      case 0:
	imin=0;
	break;
      case 100:
	imin = 1;
	break;
      case 30 :
	imin = 2;
	break;
      default :
	assert(false);
      }
      on[1]=on[2]=false;
      break;
    default:
      assert(false);
    } 
    imax = imin+1;
  }
  // calculate the form factor
  Complex FK(0.);
  for(unsigned int ix=imin;ix<imax;++ix) {
    // rho exchange
    if(on[0]) {
      Complex term = coup_[0][ix]*Resonance::BreitWignerGS(q2,mass_[0][ix],width_[0][ix],
							   mpi_,mpi_,h0_[ix],dh_[ix],hres_[ix]);
      FK += imode!=1 ? 0.5*term : -0.5*term;
    }
    // omega exchange
    if(on[1]) {
      Complex term = coup_[1][ix]*Resonance::BreitWignerFW(q2,mass_[1][ix],width_[1][ix]);
      FK += 1./6.*term;
    }
    // phi exchange
    if(on[2]) {
      Complex term = coup_[2][ix]*Resonance::BreitWignerPWave(q2,mass_[2][ix],width_[2][ix],ma,mb);
      if(ix==0 && imode==1 ) term *=etaPhi_;
      FK += term/3.;
    }
  }
  // remainder pieces
  if(imax==4) {
    if(!fKI1Re_) constructInterpolators();
    Complex i1((*fKI1Re_)(q2),(*fKI1Im_)(q2));
    FK += imode!=1 ? i1 : -i1;
    FK += Complex((*fKI0Re_)(q2),(*fKI0Im_)(q2));
  }
  // factor for cc mode
  if(imode==0) FK *= sqrt(2.0);
  return FK;
}

Complex  TwoKaonCzyzCurrent::FkaonRemainderI1(Energy2 q2) const {
  Complex output(0.);
  for(unsigned int ix=4;ix<coup_[0].size();++ix) {
    output += 0.5*coup_[0][ix]*Resonance::BreitWignerGS(q2,mass_[0][ix],width_[0][ix],
						    mpi_,mpi_,h0_[ix],dh_[ix],hres_[ix]);
  }
  return output;
}

Complex  TwoKaonCzyzCurrent::FkaonRemainderI0(Energy2 q2,Energy ma, Energy mb) const {
  Complex output(0.);
  // omega exchange
  for(unsigned int ix=4;ix<coup_[1].size();++ix) {
    output += 1./6.*coup_[1][ix]*Resonance::BreitWignerFW(q2,mass_[1][ix],width_[1][ix]);
  }
  // phi exchange
  for(unsigned int ix=4;ix<coup_[2].size();++ix) {
    output += 1./3.*coup_[2][ix]*Resonance::BreitWignerPWave(q2,mass_[2][ix],width_[2][ix],ma,mb);
  }
  return output;
}
