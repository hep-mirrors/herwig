// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEee2Mesons class.
//

#include "MEee2Mesons.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Decay/DecayIntegrator.fh"
#include "Herwig/Decay/PhaseSpaceMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/PDF/PolarizedBeamParticleData.h"
#include "Herwig/MatrixElement/HardVertex.h"

using namespace Herwig;

typedef LorentzVector<complex<InvEnergy> > LorentzPolarizationVectorInvE;

MEee2Mesons::MEee2Mesons() : flavOpt_(0)
{}

Energy2 MEee2Mesons::scale() const {
  return sHat();
}

unsigned int MEee2Mesons::orderInAlphaS() const {
  return 0;
}

unsigned int MEee2Mesons::orderInAlphaEW() const {
  return 2;
}

IBPtr MEee2Mesons::clone() const {
  return new_ptr(*this);
}

IBPtr MEee2Mesons::fullclone() const {
  return new_ptr(*this);
}

void MEee2Mesons::setFlavour() {
  flavour_ = FlavourInfo();
  if(flavOpt_==1) {
    flavour_.I  = IsoSpin::IZero;
    flavour_.I3 = IsoSpin::I3Zero;
    flavour_.strange = Strangeness::Zero;
    flavour_.charm   = Charm::Zero;
    flavour_.bottom  = Beauty::Zero;
  }
  else if(flavOpt_==2) {
    flavour_.I  = IsoSpin::IOne;
    flavour_.I3 = IsoSpin::I3Zero;
    flavour_.strange = Strangeness::Zero;
    flavour_.charm   = Charm::Zero;
    flavour_.bottom  = Beauty::Zero;
  }
  else if(flavOpt_==3) {
    flavour_.I  = IsoSpin::IZero;
    flavour_.I3 = IsoSpin::I3Zero;
    flavour_.strange = Strangeness::ssbar;
    flavour_.charm   = Charm::Zero;
    flavour_.bottom  = Beauty::Zero;
  }
  else if(flavOpt_==4) {
    flavour_.I  = IsoSpin::IZero;
    flavour_.I3 = IsoSpin::I3Zero;
    flavour_.strange = Strangeness::Zero;
    flavour_.charm   = Charm::ccbar;
    flavour_.bottom  = Beauty::Zero;
  }
  else if(flavOpt_==5) {
    flavour_.I  = IsoSpin::IZero;
    flavour_.I3 = IsoSpin::I3Zero;
    flavour_.strange = Strangeness::Zero;
    flavour_.charm   = Charm::Zero;
    flavour_.bottom  = Beauty::bbbar;
  }
}

void MEee2Mesons::doinitrun() {
  // calculate the flavour
  setFlavour();
  MEMultiChannel::doinitrun();
}

void MEee2Mesons::doinit() {
  // make sure the current got initialised
  current_->init();
  // calculate the flavour
  setFlavour();
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
    if(!current_->createMode(0,tcPDPtr(), flavour_, imode, mode,0,-1,channel,Emax)) continue;
    modeMap_[imode] = nmode;
    addMode(mode);
    ++nmode;
  }
  MEMultiChannel::doinit();
}

void MEee2Mesons::persistentOutput(PersistentOStream & os) const {
  os << current_ << modeMap_ << flavOpt_;
}

void MEee2Mesons::persistentInput(PersistentIStream & is, int) {
  is >> current_ >> modeMap_ >> flavOpt_;
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
  
  static Switch<MEee2Mesons,unsigned int> interfaceFlavour
    ("Flavour",
     "Control of the flavours of the particles in the hadronic currents",
     &MEee2Mesons::flavOpt_, 0, false, false);
  static SwitchOption interfaceFlavourAll
    (interfaceFlavour,
     "All",
     "All flavours",
     0);
  static SwitchOption interfaceFlavourI0
    (interfaceFlavour,
     "I0",
     "Only the I=0 non-strange",
     1);
  static SwitchOption interfaceFlavourI1
    (interfaceFlavour,
     "I1",
     "Only include the I=1 component",
     2);
  static SwitchOption interfaceFlavourStrange
    (interfaceFlavour,
     "Strange",
     "Only include the s sbar component",
     3);
  static SwitchOption interfaceFlavourCharm
    (interfaceFlavour,
     "Charm",
     "Only inlude the c cbar component",
     4);
  static SwitchOption interfaceFlavourBottom
    (interfaceFlavour,
     "Bottom",
     "Only include the b bbar component",
     5);

}

double MEee2Mesons::helicityME(const int ichan, const cPDVector & particles,
			       const vector<Lorentz5Momentum> & momenta) const {
  SpinorWaveFunction    em_in( momenta[0],particles[0],incoming);
  SpinorBarWaveFunction ep_in( momenta[1],particles[1],incoming);
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
  int nOut = int(momenta.size())-2;
  vector<unsigned int> constants(nOut+1);
  vector<PDT::Spin   > iSpin(nOut);
  vector<int> hadrons(nOut);
  int itemp(1);
  int ix(nOut);
  do {
    --ix;
    iSpin[ix]      = particles[ix+2]->iSpin();
    itemp         *= iSpin[ix];
    constants[ix]  = itemp;
    hadrons[ix]   = particles[ix+2]->id();
  }
  while(ix>0);
  constants[nOut] = 1;
  // calculate the matrix element
  me_.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1Half,iSpin));
  // calculate the hadron current
  unsigned int imode = current_->decayMode(hadrons);
  Energy q = sqrt(sHat());
  vector<Lorentz5Momentum> momenta2(momenta   .begin()+2,   momenta.end());
  tPDVector out = mode(modeMap_.at(imode))->outgoing();
  if(ichan<0) iMode(modeMap_.at(imode));
  vector<LorentzPolarizationVectorE> 
    hadron(current_->current(tcPDPtr(), flavour_, imode,ichan,
			     q,out,momenta2,DecayIntegrator::Calculate));
  // compute the matrix element
  vector<unsigned int> ihel(momenta.size());
  double output(0.);
  for(unsigned int hhel=0;hhel<hadron.size();++hhel) {
    // map the index for the hadrons to a helicity state
    for(int ix=nOut;ix>0;--ix) {
      ihel[ix+1]=(hhel%constants[ix-1])/constants[ix];
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
  // symmetry factors
  map<long,int> ncount;
  double symmetry(1.);
  for(tPDPtr o : out) ncount[o->id()]+=1;
  for(map<long,int>::const_iterator it=ncount.begin();it!=ncount.end();++it) {
    symmetry *= it->second;
  }
  // prefactors
  output *= 0.25*sqr(pow(sqrt(sHat())/q,int(momenta.size()-2)));
  // polarization stuff
  tcPolarizedBeamPDPtr beam[2] = 
    {dynamic_ptr_cast<tcPolarizedBeamPDPtr>(particles[0]),
     dynamic_ptr_cast<tcPolarizedBeamPDPtr>(particles[1])};
  if( beam[0] || beam[1] ) {
    RhoDMatrix rho[2] = {beam[0] ? beam[0]->rhoMatrix() : RhoDMatrix(particles[0]->iSpin()),
			 beam[1] ? beam[1]->rhoMatrix() : RhoDMatrix(particles[1]->iSpin())};
    output = me_.average(rho[0],rho[1]);
  }
  return output/symmetry;
  return output/symmetry;
}

double MEee2Mesons::me2(const int ichan) const {
  return helicityME(ichan,mePartonData(),meMomenta());
}

void MEee2Mesons::constructVertex(tSubProPtr sub) {
  // extract the particles in the hard process
  ParticleVector hard;
  hard.push_back(sub->incoming().first);
  hard.push_back(sub->incoming().second);
  for(unsigned int ix=0;ix<sub->outgoing().size();++ix)
    hard.push_back(sub->outgoing()[ix]);
  if(hard[0]->id()<hard[1]->id()) swap(hard[0],hard[1]);
  cPDVector particles;
  vector<Lorentz5Momentum> momenta;
  for(unsigned int ix=0;ix<hard.size();++ix) {
    particles.push_back(hard[ix]-> dataPtr());
    momenta  .push_back(hard[ix]->momentum());
  }
  helicityME(-1,particles,momenta);
  // construct the vertex
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(me_);  
  // wavefunctions for the incoming particles
  vector<SpinorWaveFunction>    fin;
  vector<SpinorBarWaveFunction> ain;
  SpinorWaveFunction(   fin ,hard[0],incoming,false,true);
  SpinorBarWaveFunction(ain ,hard[1],incoming,false,true);
  // and the outgoing particles
  current_->constructSpinInfo(ParticleVector(hard.begin()+2,hard.end()));
  // set the pointers and to and from the vertex
  for(unsigned int ix=0;ix<hard.size();++ix) {
    tSpinPtr spin = hard[ix]->spinInfo();
    if(ix<2) {
      tcPolarizedBeamPDPtr beam = 
	dynamic_ptr_cast<tcPolarizedBeamPDPtr>(hard[ix]->dataPtr());
      if(beam) spin->rhoMatrix() = beam->rhoMatrix();
    }
    spin->productionVertex(hardvertex);
  }
}
