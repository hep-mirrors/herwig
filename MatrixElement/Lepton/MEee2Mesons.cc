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
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig/Decay/DecayIntegrator.h"
#include "ThePEG/PDF/PolarizedBeamParticleData.h"

using namespace Herwig;

typedef LorentzVector<complex<InvEnergy> > LorentzPolarizationVectorInvE;

void MEee2Mesons::getDiagrams() const {
  // cerr << "statr of get diagrams\n";
  // make sure the current got initialised
  current_->init();
  tPDPtr gamma = getParticleData(ParticleID::gamma);
  Energy Emax = generator()->maximumCMEnergy();
  tcPDPtr em = getParticleData(ParticleID::eminus);
  tcPDPtr ep = getParticleData(ParticleID::eplus);
  // loop at the modes
  for(unsigned int imode=0;imode<current_->numberOfModes();++imode) {
    // get the external particles for this mode
    int iq(0),ia(0);
    tPDVector extpart = {gamma};
    tPDVector ptemp = current_->particles(0,imode,iq,ia);
    extpart.insert(std::end(extpart), std::begin(ptemp), std::end(ptemp));
    // create the mode
    DecayPhaseSpaceModePtr mode = new_ptr(DecayPhaseSpaceMode(extpart,DecayIntegratorPtr()));
    // create the first piece of the channel
    DecayPhaseSpaceChannelPtr channel = new_ptr(DecayPhaseSpaceChannel(mode));
    channel->addIntermediate(extpart[0],0,0.0,-1,1);
    if(!current_->createMode(0,imode,mode,1,1,channel,Emax)) continue;
    // for(unsigned int ix=0;ix<extpart.size();++ix)
    //   cerr << extpart[ix]->PDGName() << " ";
    // cerr << "\n";
    nmult_ = int(extpart.size())-1;
    int ndiag=0;
    for(unsigned int imode=0;imode<mode->numberChannels();++imode) {
      tcDecayPhaseSpaceChannelPtr iChannel = mode->channel(imode);
      // extract the intermediates
      vector<pair<unsigned int,unsigned int> > children;
      vector<tcPDPtr> intermediate;
      vector<unsigned int> jacType;
      vector<Energy> intMass;
      vector<Energy> intWidth;
      vector<double> intPower;
      for(unsigned int inter=1;inter<iChannel->numberOfIntermediates();++inter) {
	children.push_back(make_pair(0,0));
	intermediate.push_back(tcPDPtr());
	jacType.push_back(0);
	intMass.push_back(ZERO);
	intWidth.push_back(ZERO);
	intPower.push_back(0.);
	iChannel->intermediateInfo(inter,intermediate.back(),jacType.back(),
				   intMass.back(),intWidth.back(),intPower.back(),
				   children.back().first,children.back().second);
      }
      unsigned int isize=2;
      map<unsigned,unsigned int> ires;
      // create the diagram
      ThePEG::Ptr<ThePEG::Tree2toNDiagram>::pointer diag;
      for(unsigned int ix=0;ix<intermediate.size();++ix) {
	if(ix==0) {
	  diag = new_ptr((Tree2toNDiagram(2), em, ep,1,intermediate[ix]));
	  ndiag+=1;
	  // cerr << "creating " << ix << " " << diag->partons().size() << "\n";
	  isize+=1;
	  ires[ix] = isize;
	}
	if(children[ix].first>0) {
	  diag = new_ptr((*diag,ires[ix],extpart[children[ix].first]));
	  isize+=1;
	}
	else {
	  int iloc = -children[ix].first-1;
	  isize+=1;
	  ires[iloc] = isize;
	  diag = new_ptr((*diag,ires[ix],intermediate[iloc]));
	}
	if(children[ix].second>0) {
	  diag = new_ptr((*diag,ires[ix],extpart[children[ix].second]));
	  isize+=1;
	}
	else {
	  int iloc = -children[ix].second-1;
	  isize+=1;
	  ires[iloc] = isize;
	  diag = new_ptr((*diag,ires[ix],intermediate[iloc]));
	}
	// cerr << "testing in the loop " << ix << " " << diag->partons().size() << "\n";
      }
      diag = new_ptr((*diag,-ndiag));
      add(diag);
    }
  }
}

Energy2 MEee2Mesons::scale() const {
  return sHat();
}

int MEee2Mesons::nDim() const {
  return 3*nmult_-4;
}

void MEee2Mesons::setKinematics() {
  HwMEBase::setKinematics();
}

bool MEee2Mesons::generateKinematics(const double * r) {
  // cerr << "testing in generate\n";
  // Save the jacobian dPS/dr for later use.
  jacobian(1.0);
  if(nDim()==2) {
    // set up the momenta
    for ( int i = 2, N = meMomenta().size(); i < N; ++i ) {
      meMomenta()[i] = Lorentz5Momentum(mePartonData()[i]->mass());
    }
    double ctmin = -1.0, ctmax = 1.0;
    Energy q = ZERO;
    try {
      q = SimplePhaseSpace::
	getMagnitude(sHat(), meMomenta()[2].mass(), meMomenta()[3].mass());
    } 
    catch ( ImpossibleKinematics ) {
      return false;
    }
    Energy e = sqrt(sHat())/2.0;
    
    Energy2 m22 = meMomenta()[2].mass2();
    Energy2 m32 = meMomenta()[3].mass2();
    Energy2 e0e2 = 2.0*e*sqrt(sqr(q) + m22);
    Energy2 e1e2 = 2.0*e*sqrt(sqr(q) + m22);
    Energy2 e0e3 = 2.0*e*sqrt(sqr(q) + m32);
    Energy2 e1e3 = 2.0*e*sqrt(sqr(q) + m32);
    Energy2 pq = 2.0*e*q;
    
    Energy2 thmin = lastCuts().minTij(mePartonData()[0], mePartonData()[2]);
    if ( thmin > ZERO ) ctmax = min(ctmax, (e0e2 - m22 - thmin)/pq);
    
    thmin = lastCuts().minTij(mePartonData()[1], mePartonData()[2]);
    if ( thmin > ZERO ) ctmin = max(ctmin, (thmin + m22 - e1e2)/pq);
    
    thmin = lastCuts().minTij(mePartonData()[1], mePartonData()[3]);
    if ( thmin > ZERO ) ctmax = min(ctmax, (e1e3 - m32 - thmin)/pq);
    
    thmin = lastCuts().minTij(mePartonData()[0], mePartonData()[3]);
    if ( thmin > ZERO ) ctmin = max(ctmin, (thmin + m32 - e0e3)/pq);
    
    Energy ptmin = max(lastCuts().minKT(mePartonData()[2]),
		       lastCuts().minKT(mePartonData()[3]));
    if ( ptmin > ZERO ) {
      double ctm = 1.0 - sqr(ptmin/q);
      if ( ctm <= 0.0 ) return false;
      ctmin = max(ctmin, -sqrt(ctm));
      ctmax = min(ctmax, sqrt(ctm));
    }
    
    double ymin2 = lastCuts().minYStar(mePartonData()[2]);
    double ymax2 = lastCuts().maxYStar(mePartonData()[2]);
    double ymin3 = lastCuts().minYStar(mePartonData()[3]);
    double ymax3 = lastCuts().maxYStar(mePartonData()[3]);
    double ytot = lastCuts().Y() + lastCuts().currentYHat();
    if ( ymin2 + ytot > -0.9*Constants::MaxRapidity )
      ctmin = max(ctmin, sqrt(sqr(q) +  m22)*tanh(ymin2)/q);
    if ( ymax2 + ytot < 0.9*Constants::MaxRapidity )
      ctmax = min(ctmax, sqrt(sqr(q) +  m22)*tanh(ymax2)/q);
    if ( ymin3 + ytot > -0.9*Constants::MaxRapidity )
      ctmax = min(ctmax, sqrt(sqr(q) +  m32)*tanh(-ymin3)/q);
    if ( ymax3 + ytot < 0.9*Constants::MaxRapidity )
      ctmin = max(ctmin, sqrt(sqr(q) +  m32)*tanh(-ymax3)/q);
    
    if ( ctmin >= ctmax ) return false;
    
    double cth = getCosTheta(ctmin, ctmax, r[0]);
    Energy pt = q*sqrt(1.0-sqr(cth));
    phi(rnd(2.0*Constants::pi));
    meMomenta()[2].setVect(Momentum3( pt*sin(phi()),  pt*cos(phi()),  q*cth));
    meMomenta()[3].setVect(Momentum3(-pt*sin(phi()), -pt*cos(phi()), -q*cth));
    
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
    jacobian((pq/sHat())*Constants::pi*jacobian());
    return true;
  }
  else {
    assert(false);
  }
}

double MEee2Mesons::me2() const {
  // cerr << "testing in me2\n";
  // for(unsigned int ix=0;ix<mePartonData().size();++ix)
  //   cerr << mePartonData()[ix]->PDGName() << " " << meMomenta()[ix]/GeV << "\n";
  // wavefunctions for the incoming leptons
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
  // cerr << "testing constants\n";
  // for(unsigned int ix=0;ix<constants.size();++ix)
  //   cerr << constants[ix] << " ";
  // cerr << "\n";
  // cerr << "testing ispn\n";
  // for(unsigned int ix=0;ix<ispin.size();++ix)
  //   cerr << ispin[ix] << " ";
  // cerr << "\n";
  // calculate the matrix element
  me_.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1Half,ispin));
  // calculate the hadron current
  unsigned int imode = current_->decayMode(hadrons);
  Energy q = sqrt(sHat());
  // cerr << "testing mode " << imode << " " << q/GeV << "\n";
  ParticleVector hadpart;
  for(unsigned int ix=2;ix<mePartonData().size();++ix) {
    hadpart.push_back(mePartonData()[ix]->produceParticle(meMomenta()[ix]));
  }
  vector<LorentzPolarizationVectorE> 
    hadron(current_->current(imode,0,q,hadpart,DecayIntegrator::Calculate));
  // for(unsigned int ix=0;ix<hadron.size();++ix)
  //   cerr << hadron[ix].x()/GeV << " " << hadron[ix].y()/GeV << " " << hadron[ix].z()/GeV << " " << hadron[ix].t()/GeV << "\n"; 
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
	// cerr << "testing ME  ";
	// for(unsigned int ix=0;ix<ihel.size();++ix) cerr << ihel[ix] << " ";
	// cerr << me_(ihel) << "\n";
      }
    }
  }
  output *= 0.25;
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

CrossSection MEee2Mesons::dSigHatDR() const {
  return me2()*jacobian()/(16.0*sqr(Constants::pi)*sHat())*sqr(hbarc);
}

unsigned int MEee2Mesons::orderInAlphaS() const {
  return 0;
}

unsigned int MEee2Mesons::orderInAlphaEW() const {
  return 0;
}

Selector<MEBase::DiagramIndex>
MEee2Mesons::diagrams(const DiagramVector & diags) const {
  // This example corresponds to the diagrams specified in the example
  // in the getDiagrams() function.

  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) 
    if ( diags[i]->id() == -1 ) sel.insert(1.0, i);
    else if ( diags[i]->id() == -2 )  sel.insert(1.0, i);
    else if ( diags[i]->id() == -3 )  sel.insert(1.0, i);
  // You probably do not want equal weights here...
  return sel;

  // If there is only one possible diagram you can override the
  // MEBase::diagram function instead.

}

Selector<const ColourLines *>
MEee2Mesons::colourGeometries(tcDiagPtr) const {
  static ColourLines none("");
  Selector<const ColourLines *> sel;
  sel.insert(1.0, &none);
  return sel;
}

IBPtr MEee2Mesons::clone() const {
  return new_ptr(*this);
}

IBPtr MEee2Mesons::fullclone() const {
  return new_ptr(*this);
}

void MEee2Mesons::doinit() {
  HwMEBase::doinit();
}

void MEee2Mesons::persistentOutput(PersistentOStream & os) const {
  os << current_ << nmult_;
}

void MEee2Mesons::persistentInput(PersistentIStream & is, int) {
  is >> current_ >> nmult_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEee2Mesons,HwMEBase>
  describeHerwigMEee2Mesons("Herwig::MEee2Mesons", "HwMELeptonLowEnergy.so");

void MEee2Mesons::Init() {

  static ClassDocumentation<MEee2Mesons> documentation
    ("The MEee2Mesons class simluation the production of low multiplicity"
     " events via the weak current");

  static Reference<MEee2Mesons,WeakDecayCurrent> interfaceWeakCurrent
    ("WeakCurrent",
     "The reference for the decay current to be used.",
     &MEee2Mesons::current_, false, false, true, false, false);

}


void MEee2Mesons::constructVertex(tSubProPtr sub) {
}
