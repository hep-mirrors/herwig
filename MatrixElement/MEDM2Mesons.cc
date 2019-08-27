// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEDM2Mesons class.
//

#include "MEDM2Mesons.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"

using namespace Herwig;
typedef LorentzVector<complex<InvEnergy> > LorentzPolarizationVectorInvE;

MEDM2Mesons::MEDM2Mesons() {
  cSMmed_ = {1.0,1.0};
}

Energy2 MEDM2Mesons::scale() const {
  return sHat();
}

unsigned int MEDM2Mesons::orderInAlphaS() const {
  return 0;
}

unsigned int MEDM2Mesons::orderInAlphaEW() const {
  return 0;
}

IBPtr MEDM2Mesons::clone() const {
  return new_ptr(*this);
}

IBPtr MEDM2Mesons::fullclone() const {
  return new_ptr(*this);
}

void MEDM2Mesons::doinit() {
  // make sure the current got initialised
  current_->init();
  // max energy
  Energy Emax = generator()->maximumCMEnergy();
  // loop over the modes
  int nmode=0;
  for(unsigned int imode=0;imode<current_->numberOfModes();++imode) {
     // get the external particles for this mode
     int iq(0),ia(0);
     tPDVector out = current_->particles(0,imode,iq,ia);
    current_->decayModeInfo(imode,iq,ia);
     if(iq==2&&ia==-2) continue;
     PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(incomingA_,out,1.,
						     incomingB_,Emax));
     PhaseSpaceChannel channel(mode);
     if(!current_->createMode(0,tcPDPtr(), IsoSpin::IUnknown, IsoSpin::I3Unknown,
			      imode,mode,0,-1,channel,Emax)) continue;
     modeMap_[imode] = nmode;
     addMode(mode);
     ++nmode;
  }
  MEMultiChannel::doinit();
}

void MEDM2Mesons::persistentOutput(PersistentOStream & os) const {
  os << current_ << modeMap_ << incomingA_ << incomingB_ << Mediator_ << cDMmed_ << cSMmed_ << ounit(MMed_,GeV);
}

void MEDM2Mesons::persistentInput(PersistentIStream & is, int) {
  is >> current_ >> modeMap_ >> incomingA_ >> incomingB_ >> Mediator_ >> cDMmed_ >> cSMmed_ >> iunit(MMed_,GeV);
}

//The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEDM2Mesons,MEMultiChannel>
  describeHerwigMEDM2Mesons("Herwig::MEDM2Mesons", "Herwig.so");

void MEDM2Mesons::Init() {

  static ClassDocumentation<MEDM2Mesons> documentation
    ("The MEDM2Mesons class simulates the annhilation of"
     " DM particles to mesons at low energy");

  static Reference<MEDM2Mesons,WeakCurrent> interfaceWeakCurrent
    ("WeakCurrent",
     "The reference for the decay current to be used.",
     &MEDM2Mesons::current_, false, false, true, false, false);
  
  static Reference<MEDM2Mesons,ParticleData> interfaceIncomingA
    ("IncomingA",
     "First incoming particle",
     &MEDM2Mesons::incomingA_, false, false, true, false, false);
  
  static Reference<MEDM2Mesons,ParticleData> interfaceIncomingB
    ("IncomingB",
     "Second incoming particle",
     &MEDM2Mesons::incomingB_, false, false, true, false, false);

  static Reference<MEDM2Mesons,ParticleData> interfaceMediator_
    ("Mediator",
     "DM mediator",
     &MEDM2Mesons::Mediator_, false, false, true, false, false);

  static Parameter<MEDM2Mesons,double> interfacecDMmed
    ("cDMmed",
     "coupling of DM to dark mediator",
     &MEDM2Mesons::cDMmed_, 1.0, 0., 10., false, false, Interface::limited);

  static ParVector<MEDM2Mesons,double> interfacecSMmed
    ("cSMmed",
     "coupling of SM to dark mediator",
     &MEDM2Mesons::cSMmed_, -1 , 1.0 , 0. , 10. , false, false, Interface::limited);

}


double MEDM2Mesons::me2(const int ichan) const {
  // compute the incoming current
  LorentzPolarizationVectorInvE lepton[2][2];
  if(incomingA_->iSpin()==PDT::Spin1Half && incomingB_->iSpin()==PDT::Spin1Half) {
    SpinorWaveFunction    em_in( meMomenta()[0],mePartonData()[0],incoming);
    SpinorBarWaveFunction ep_in( meMomenta()[1],mePartonData()[1],incoming);
    vector<SpinorWaveFunction> f1;
    vector<SpinorBarWaveFunction> a1;
    for(unsigned int ix=0;ix<2;++ix) {
      em_in.reset(ix);
      f1.push_back(em_in);
      ep_in.reset(ix);
      a1.push_back(ep_in);
    }
    // this should be coupling of DM to mediator/ mediator propagator
    complex<Energy> mmed = Mediator_->mass();
    complex<Energy2> mmed2 = sqr(mmed);
    complex<Energy> mwid = Mediator_->width();
    complex<Energy2> prop = sHat()-mmed2+Complex(0.,1.)*mmed*mwid;
    complex<InvEnergy2> pre = cDMmed_/prop;
    for(unsigned ix=0;ix<2;++ix) {
      for(unsigned iy=0;iy<2;++iy) {
	lepton[ix][iy]= pre*f1[ix].dimensionedWave().vectorCurrent(a1[iy].dimensionedWave());
      }
    }
  }
  // TODO think about other spins for the DM
  else
    assert(false);
  // work out the mapping for the hadron vector
  int nOut = int(meMomenta().size())-2;
  vector<unsigned int> constants(nOut+1);
  vector<PDT::Spin   > iSpin(nOut);
  vector<int> hadrons(nOut);
  int itemp(1);
  int ix(nOut);
  do {
    --ix;
    iSpin[ix]      = mePartonData()[ix+2]->iSpin();
    itemp         *= iSpin[ix];
    constants[ix]  = itemp;
    hadrons[ix]   = mePartonData()[ix+2]->id();
  }
  while(ix>0);
  constants[nOut] = 1;
  // calculate the matrix element
  me_.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1Half,iSpin));
  // calculate the hadron current
  unsigned int imode = current_->decayMode(hadrons);
  Energy q = sqrt(sHat());
  vector<Lorentz5Momentum> momenta(meMomenta()   .begin()+2,   meMomenta().end());
  tPDVector out = mode(modeMap_.at(imode))->outgoing();
  if(ichan<0) iMode(modeMap_.at(imode));
  // get the hadronic currents for the I=1 and I=0 components
  vector<LorentzPolarizationVectorE> 
    hadronI0(current_->current(tcPDPtr(), IsoSpin::IZero, IsoSpin::I3Zero, imode,ichan,
			       q,out,momenta,DecayIntegrator::Calculate));
  vector<LorentzPolarizationVectorE> 
    hadronI1(current_->current(tcPDPtr(), IsoSpin::IOne, IsoSpin::I3Zero, imode,ichan,
			       q,out,momenta,DecayIntegrator::Calculate));
  // compute the matrix element
  vector<unsigned int> ihel(meMomenta().size());
  double output(0.);
  int hI0_size = hadronI0.size();
  int hI1_size = hadronI1.size();
  int maxsize = max(hadronI0.size(),hadronI1.size());
  for(unsigned int hhel=0;hhel<maxsize;++hhel) {
    // map the index for the hadrons to a helicity state
    for(int ix=nOut;ix>0;--ix) {
      ihel[ix+1]=(hhel%constants[ix-1])/constants[ix];
    }
    // loop over the helicities of the incoming particles
    for(ihel[1]=0;ihel[1]<2;++ihel[1]){
      for(ihel[0]=0;ihel[0]<2;++ihel[0]) {
	Complex amp;
	// work on coefficients for the I1 and I0 bits
  	if(hI0_size != 0 && hI1_size !=0){
	  amp = lepton[ihel[0]][ihel[1]].dot(cSMmed_[0]*hadronI0[hhel]+cSMmed_[1]*hadronI1[hhel]);
	}
	else if(hI0_size != 0 && hI1_size == 0){
	  amp = lepton[ihel[0]][ihel[1]].dot(cSMmed_[0]*hadronI0[hhel]);
	}
	else {
	  amp = lepton[ihel[0]][ihel[1]].dot(cSMmed_[1]*hadronI1[hhel]);
	}
	me_(ihel)= amp;
  	output += std::norm(amp);
      }
    }
  }
  // symmetry factors
  map<long,int> ncount;
  double symmetry(1.);
  for(tPDPtr o : out) ncount[o->id()]+=1;
  for(map<long,int>::const_iterator it=ncount.begin();it!=ncount.end();++it) {
    symmetry *= it->second;
  }
  // prefactors
  output *= 0.25*sqr(pow(sqrt(sHat())/q,int(momenta.size()-2)));
  return output/symmetry;
}


void MEDM2Mesons::constructVertex(tSubProPtr) {
}
