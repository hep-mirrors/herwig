// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEee2Mesons class.
//

#include "MEee2Mesons.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Decay/DecayIntegrator2.fh"
#include "Herwig/Decay/PhaseSpaceMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/PDF/PolarizedBeamParticleData.h"

using namespace Herwig;

typedef LorentzVector<complex<InvEnergy> > LorentzPolarizationVectorInvE;

MEee2Mesons::MEee2Mesons() {}

Energy2 MEee2Mesons::scale() const {
  return sHat();
}

unsigned int MEee2Mesons::orderInAlphaS() const {
  return 0;
}

unsigned int MEee2Mesons::orderInAlphaEW() const {
  return 0;
}

IBPtr MEee2Mesons::clone() const {
  return new_ptr(*this);
}

IBPtr MEee2Mesons::fullclone() const {
  return new_ptr(*this);
}

void MEee2Mesons::doinit() {
  // make sure the current got initialised
  current_->init();
  // max energy
  Energy Emax = generator()->maximumCMEnergy();
  // incoming particles
  tPDPtr em = getParticleData(ParticleID::eminus);
  tPDPtr ep = getParticleData(ParticleID::eplus);
  // loop over the modes
  int nmode=0;
  for(unsigned int imode=0;imode<current_->numberOfModes();++imode) {
    // get the external particles for this mode
    int iq(0),ia(0);
    tPDVector out = current_->particles(0,imode,iq,ia);
    current_->decayModeInfo(imode,iq,ia);
    if(iq==2&&ia==-2) continue;
    PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(em,out,1.,ep,Emax));
    PhaseSpaceChannel channel(mode);
    if(!current_->createMode(0,tcPDPtr(), IsoSpin::IUnknown, IsoSpin::I3Unknown,
			     imode,mode,0,-1,channel,Emax)) continue;
    modeMap_[imode] = nmode;
    addMode(mode);
    ++nmode;
  }
  MEMultiChannel::doinit();
}

void MEee2Mesons::persistentOutput(PersistentOStream & os) const {
os << current_ << modeMap_;
}

void MEee2Mesons::persistentInput(PersistentIStream & is, int) {
is >> current_ >> modeMap_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEee2Mesons,MEMultiChannel>
describeHerwigMEee2Mesons("Herwig::MEee2Mesons", "HwMELeptonLowEnergy.so");

void MEee2Mesons::Init() {

  static ClassDocumentation<MEee2Mesons> documentation
    ("The MEee2Mesons class simluation the production of low multiplicity"
     " events via the weak current");

  static Reference<MEee2Mesons,WeakCurrent> interfaceWeakCurrent
    ("WeakCurrent",
     "The reference for the decay current to be used.",
     &MEee2Mesons::current_, false, false, true, false, false);
  
}

double MEee2Mesons::me2(const int ichan) const {
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
  // compute the leptonic current
  LorentzPolarizationVectorInvE lepton[2][2];
  InvEnergy2 pre = SM().alphaEM(sHat())*4.*Constants::pi/sHat();
  for(unsigned ix=0;ix<2;++ix) {
    for(unsigned iy=0;iy<2;++iy) {
      lepton[ix][iy]= pre*f1[ix].dimensionedWave().vectorCurrent(a1[iy].dimensionedWave());
    }
  }
  // work out the mapping for the hadron vector
  vector<unsigned int> constants(meMomenta().size()+1);
  vector<PDT::Spin   > ispin(meMomenta().size()-2);
  vector<int> hadrons(meMomenta().size()-2);
  int itemp(1);
  unsigned int ix(meMomenta().size());
  do {
    --ix;
    ispin[ix-2]     = mePartonData()[ix]->iSpin();
    hadrons[ix-2]   = mePartonData()[ix]->id();
    itemp         *= ispin[ix-2];
    constants[ix] = itemp;
  }
  while(ix>2);
  constants[meMomenta().size()] = 1;
  constants[0] = constants[1] = constants[2];
  // calculate the matrix element
  me_.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1Half,ispin));
  // calculate the hadron current
  unsigned int imode = current_->decayMode(hadrons);
  Energy q = sqrt(sHat());
  vector<Lorentz5Momentum> momenta(meMomenta()   .begin()+2,   meMomenta().end());
  tPDVector out = mode(modeMap_.at(imode))->outgoing();
  if(ichan<0) iMode(modeMap_.at(imode));
  vector<LorentzPolarizationVectorE> 
    hadron(current_->current(tcPDPtr(), IsoSpin::IUnknown, IsoSpin::I3Unknown, imode,ichan,
			     q,out,momenta,DecayIntegrator2::Calculate));
  // compute the matrix element
  vector<unsigned int> ihel(meMomenta().size());
  double output(0.);
  for(unsigned int hhel=0;hhel<hadron.size();++hhel) {
    // map the index for the hadrons to a helicity state
    for(unsigned int ix=meMomenta().size()-1;ix>1;--ix) {
      ihel[ix]=(hhel%constants[ix-1])/constants[ix];
    }
    // loop over the helicities of the incoming leptons
    for(ihel[1]=0;ihel[1]<2;++ihel[1]){
      for(ihel[0]=0;ihel[0]<2;++ihel[0]) {
	Complex amp = lepton[ihel[0]][ihel[1]].dot(hadron[hhel]);
   	me_(ihel)= amp;
	output += std::norm(amp);
      }
    }
  }
  // prefactors
  output *= 0.25*sqr(pow(sqrt(sHat())/q,int(momenta.size()-2)));
  // polarization stuff
  tcPolarizedBeamPDPtr beam[2] = 
    {dynamic_ptr_cast<tcPolarizedBeamPDPtr>(mePartonData()[0]),
     dynamic_ptr_cast<tcPolarizedBeamPDPtr>(mePartonData()[1])};
  if( beam[0] || beam[1] ) {
    RhoDMatrix rho[2] = {beam[0] ? beam[0]->rhoMatrix() : RhoDMatrix(mePartonData()[0]->iSpin()),
			 beam[1] ? beam[1]->rhoMatrix() : RhoDMatrix(mePartonData()[1]->iSpin())};
    output = me_.average(rho[0],rho[1]);
  }
  return output;
}

void MEee2Mesons::constructVertex(tSubProPtr) {
}
