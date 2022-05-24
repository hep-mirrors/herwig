// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ScalarTo3ScalarDalitz class.
//

#include "ScalarTo3ScalarDalitz.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

ScalarTo3ScalarDalitz::ScalarTo3ScalarDalitz(InvEnergy rP, bool useResonanceMass)
  : rParent_(rP), useResonanceMass_(useResonanceMass), maxWgt_(1.),
    channel1_(-1), channel2_(-1), incoming_(0), outgoing_({0,0,0}) {
  // intermediates
  generateIntermediates(true);
}

IBPtr ScalarTo3ScalarDalitz::clone() const {
  return new_ptr(*this);
}

IBPtr ScalarTo3ScalarDalitz::fullclone() const {
  return new_ptr(*this);
}

void ScalarTo3ScalarDalitz::persistentOutput(PersistentOStream & os) const {
  os << resonances_ << maxWgt_ << weights_ << ounit(rParent_,1./GeV)
     << channel1_ << channel2_ << incoming_ << outgoing_ << useResonanceMass_;
}

void ScalarTo3ScalarDalitz::persistentInput(PersistentIStream & is, int) {
  is >> resonances_ >> maxWgt_ >> weights_ >> iunit(rParent_,1./GeV)
     >> channel1_ >> channel2_ >> incoming_ >> outgoing_ >> useResonanceMass_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ScalarTo3ScalarDalitz,DecayIntegrator>
describeHerwigScalarTo3ScalarDalitz("Herwig::ScalarTo3ScalarDalitz", "HwDalitzDecay.so");

void ScalarTo3ScalarDalitz::Init() {

  static ClassDocumentation<ScalarTo3ScalarDalitz> documentation
    ("The ScalarTo3ScalarDalitz class provides a base class for "
     "weak three-body decays of bottom and charm mesons");

  static Command<ScalarTo3ScalarDalitz> interfaceSetExternal
    ("SetExternal",
     "Set the external particles for the decay mode",
     &ScalarTo3ScalarDalitz::setExternal, false);
  
  static Command<ScalarTo3ScalarDalitz> interfaceAddChannel
    ("AddChannel",
     "Add a channel for the description of the matrix element",
     &ScalarTo3ScalarDalitz::addChannel, false);

  static Parameter<ScalarTo3ScalarDalitz,InvEnergy> interfaceParentRadius
    ("ParentRadius",
     "The radius parameter for the Blatt-Weisskopf form-factor for the D",
     &ScalarTo3ScalarDalitz::rParent_, 1./GeV, 5./GeV, ZERO, 10./GeV,
     false, false, Interface::limited);

  static Parameter<ScalarTo3ScalarDalitz,double> interfaceMaximumWeight
    ("MaximumWeight",
     "The maximum weight for the phase-space sampling",
     &ScalarTo3ScalarDalitz::maxWgt_, 1.0, 0.0, 1e20,
     false, false, Interface::limited);

  static ParVector<ScalarTo3ScalarDalitz,double> interfaceWeights
    ("Weights",
     "The weights for the different channels for the phase-space integration",
     &ScalarTo3ScalarDalitz::weights_, -1, 1.0, 0.0, 1.0,
     false, false, Interface::limited);
  
  static Parameter<ScalarTo3ScalarDalitz,int> interfaceChannel1
    ("Channel1",
     "The first allowed channel, for debugging/calculation of fit fractions only",
     &ScalarTo3ScalarDalitz::channel1_, -1, -1, 100,
     false, false, Interface::limited);
  
  static Parameter<ScalarTo3ScalarDalitz,int> interfaceChannel2
    ("Channel2",
     "The first allowed channel, for debugging/calculation of fit fractions only",
     &ScalarTo3ScalarDalitz::channel2_, -1, -1, 100,
     false, false, Interface::limited);

  static Switch<ScalarTo3ScalarDalitz,bool> interfaceResonanceMass
    ("ResonanceMass",
     "Whether to use the kinematic mass or the resonance pole mass for the denominator in kinematic expressions",
     &ScalarTo3ScalarDalitz::useResonanceMass_, false, false, false);
  static SwitchOption interfaceResonanceMassYes
    (interfaceResonanceMass,
     "Yes",
     "Use the resonance mass, to be avoided only use if do in experimental fit",
     true);
  static SwitchOption interfaceResonanceMassNo
    (interfaceResonanceMass,
     "No",
     "Use the correct kinematic mass",
     false);

}

void ScalarTo3ScalarDalitz::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // set up the spin information for the decay products
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					incoming,true);
  for(unsigned int ix=0;ix<3;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
}

void ScalarTo3ScalarDalitz::doinit() {
  if(incoming_!=0) {
    tPDPtr in = getParticleData(incoming_);
    vector<tPDPtr> out = {getParticleData(outgoing_[0]),
			  getParticleData(outgoing_[1]),
			  getParticleData(outgoing_[2])};
    createMode(in,out);
  }
  DecayIntegrator::doinit();
}

void ScalarTo3ScalarDalitz::doinitrun() {
  DecayIntegrator::doinitrun();
  weights_.resize(mode(0)->channels().size());
  maxWgt_ = mode(0)->maxWeight();
  for(unsigned int iz=0;iz<mode(0)->channels().size();++iz) {
    weights_[iz]=mode(0)->channels()[iz].weight();
  }
}

double ScalarTo3ScalarDalitz::me2(const int ichan, const Particle & part,
			    const tPDVector & ,
			    const vector<Lorentz5Momentum> & momenta,
			    MEOption meopt) const {
  static const Complex ii(0.,1.);
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin0,PDT::Spin0,PDT::Spin0)));
  useMe();
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&part),incoming);
  }
  // set the kinematics
  mD_ = part.mass();
  for(unsigned int ix=0;ix<momenta.size();++ix) {
    mOut_[ix]=momenta[ix].mass();
    for(unsigned int iy=ix+1;iy<momenta.size();++iy) {
      m2_[ix][iy]=(momenta[ix]+momenta[iy]).m();
      m2_[iy][ix]=m2_[ix][iy];
    }
  }
  // compute the amplitude
  Complex amp(0.);
  int iloc=-1;
  for(int ix=0;ix<int(resonances().size());++ix) {
    ++iloc;
    if(channel1_>=0) {
      if(ix!=channel1_ && ix!=channel2_) continue;
    }
    if(ichan>=0&&ichan!=iloc) continue;
    amp += resAmp(ix);
  }
  // now compute the matrix element
  (*ME())(0,0,0,0) = amp;
  return norm(amp);
}

void ScalarTo3ScalarDalitz::createMode(tPDPtr in, tPDVector out) {
  PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(in,out,maxWgt_));
  if(weights_.size()!=resonances_.size()) {
    weights_=vector<double>(resonances_.size(),1./double(resonances_.size()));
  }
  unsigned int ix=0;
  for(DalitzResonance res : resonances_) {
    tPDPtr resonance = getParticleData(res.id);
    if(resonance) {
      mode->addChannel((PhaseSpaceChannel(mode),0,resonance,0,res.spectator+1,1,res.daughter1+1,1,res.daughter2+1));
      resetIntermediate(resonance,res.mass,res.width);
      ++ix;
    }
  }
  addMode(mode);
}

Complex ScalarTo3ScalarDalitz::resAmp(unsigned int i) const {
  static const Complex ii = Complex(0.,1.);
  Complex output = resonances_[i].amp;
  if (resonances_[i].type==ResonanceType::NonResonant) return output;
  // locations of the outgoing particles
  const unsigned int &d1 = resonances_[i].daughter1;
  const unsigned int &d2 = resonances_[i].daughter2;
  const unsigned int &sp = resonances_[i].spectator;
  // mass and width of the resonance
  const Energy & mR = resonances_[i].mass ;
  const Energy & wR = resonances_[i].width;
  // momenta for the resonance decay
  // off-shell
  Energy pAB=sqrt(0.25*sqr(sqr(m2_[d1][d2]) -sqr(mOut_[d1])-sqr(mOut_[d2])) - sqr(mOut_[d1]*mOut_[d2]))/m2_[d1][d2];
  if(resonances_[i].type==ResonanceType::BABARf0) {
    double rho = 2.*pAB/m2_[d1][d2];
    output *= GeV2/(sqr(resonances_[i].mass)-sqr(m2_[d1][d2])-
		    ii*resonances_[i].mass*resonances_[i].width*rho);
    return output;
  }
  //  on-shell
  Energy  pR=sqrt(0.25*sqr(    mR*mR        -sqr(mOut_[d1])-sqr(mOut_[d2])) - sqr(mOut_[d1]*mOut_[d2]))/mR;
  // Blatt-Weisskopf factors
  double fR=1, fD=1;
  unsigned int power(1);
  if(resonances_[i].type!=ResonanceType::Spin0 &&
     resonances_[i].type!=ResonanceType::Spin0E691) {
    // for the D decay
    Energy pD  = sqrt(max(ZERO,(0.25*sqr(sqr(mD_)-sqr(mR)-sqr(mOut_[sp])) - sqr(mR*mOut_[sp]))/sqr(mD_)));
    Energy pDAB= sqrt( 0.25*sqr(sqr(mD_)-sqr(m2_[d1][d2])-sqr(mOut_[sp])) - sqr(m2_[d1][d2]*mOut_[sp]))/mD_;
    double r1A(resonances_[i].R*pR),r1B(resonances_[i].R*pAB );
    double r2A(rParent_   *pD),r2B(rParent_   *pDAB);
    // mass for thre denominator
    Energy mDen = useResonanceMass_ ? mR : m2_[d1][d2];
    // denominator for the older form of the amplitude
    Energy2 denom = GeV2;
    if (resonances_[i].type/10 == 1 ) { 
      Energy2 pa2 = 0.25*(sqr(m2_[d1][d2])-2.*(sqr(mOut_[d1])+sqr(mOut_[d2])) + sqr(sqr(mOut_[d1])-sqr(mOut_[d2]))/sqr(m2_[d1][d2]));
      Energy2 pc2 = 0.25*(sqr(m2_[d1][d2])-2.*(sqr(mD_      )+sqr(mOut_[sp])) + sqr(sqr(mD_      )-sqr(mOut_[sp]))/sqr(m2_[d1][d2]));
      denom = 4.*sqrt(pa2*pc2);
    }
    // Blatt-Weisskopf factors and spin piece
    switch (resonances_[i].type) {
    case ResonanceType::Spin0Gauss:
      fR = exp(-(r1B-r1A)/12.);
      fD = exp(-(r2B-r2A)/12.);
      break;
    case ResonanceType::Spin1: case ResonanceType::Spin1E691 :
      fR=sqrt( (1. + sqr(r1A)) / (1. + sqr(r1B)) );
      fD=sqrt( (1. + sqr(r2A)) / (1. + sqr(r2B)) );
      power=3;
      output *= fR*fD*(sqr(m2_[d2][sp])-sqr(m2_[d1][sp])
		       + (  sqr(mD_)-sqr(mOut_[sp]))*(sqr(mOut_[d1])-sqr(mOut_[d2]))/sqr(mDen) )/denom;
      break;
    case ResonanceType::Spin2: case ResonanceType::Spin2E691:
      fR = sqrt( (9. + sqr(r1A)*(3.+sqr(r1A))) / (9. + sqr(r1B)*(3.+sqr(r1B))));
      fD = sqrt( (9. + sqr(r2A)*(3.+sqr(r2A))) / (9. + sqr(r2B)*(3.+sqr(r2B))));
      power=5;
      output *= fR*fD/sqr(denom)*( sqr( sqr(m2_[d2][sp]) - sqr(m2_[d1][sp]) +(sqr(mD_)-sqr(mOut_[sp]))*(sqr(mOut_[d1])-sqr(mOut_[d2]))/(sqr(mDen))) -
				   (sqr(m2_[d1][d2])-2*      sqr(mD_)-2*sqr(mOut_[sp]) + sqr((sqr(      mD_) - sqr(mOut_[sp]))/mDen))*
				   (sqr(m2_[d1][d2])-2*sqr(mOut_[d1])-2*sqr(mOut_[d2]) + sqr((sqr(mOut_[d1]) - sqr(mOut_[d2]))/mDen))/3.);
      break;
    default :
      assert(false);
    }
  }
  // multiply by Breit-Wigner piece and return
  if (resonances_[i].type/10 == 1 ) {
    return output*sqrt(0.5*wR/GeV/Constants::pi)*GeV/(m2_[d1][d2]-mR-complex<Energy>(ZERO,0.5*wR));
  }
  else {
    Energy gam = wR*pow(pAB/pR,power)*(mR/m2_[d1][d2])*fR*fR;
    return output*GeV2/(sqr(mR)-sqr(m2_[d1][d2])-mR*gam*ii);
  }
}

void ScalarTo3ScalarDalitz::dataBaseOutput(ofstream & output, bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  output << "newdef " << name() << ":ParentRadius " << rParent_*GeV << "\n";
  output << "newdef " << name() << ":ResonanceMass " << useResonanceMass_ << "\n";
  output << "newdef " << name() << ":MaximumWeight " << maxWgt_ << "\n";
  for(unsigned int ix=0;ix<weights_.size();++ix) {
    output << "insert " << name() << ":Weights "
	   << ix << " " << weights_[ix] << "\n";
  }
  for(unsigned int ix=0;ix<resonances_.size();++ix) {
    output << "do " << name() << ":AddChannel "
	   << resonances_[ix].id << " " << oenum(resonances_[ix].type) << " "
	   << resonances_[ix].mass/GeV << " " << resonances_[ix].width/GeV << " "
	   << resonances_[ix].daughter1 << " " << resonances_[ix].daughter2 << " "
	   << resonances_[ix].spectator << " " 
	   << abs(resonances_[ix].amp) << " " << arg(resonances_[ix].amp) << " "
	   << resonances_[ix].R*GeV << "\n"; 
  }
  output << "do " << name() << ":SetExternal " << incoming_;
  for(unsigned int ix=0;ix<3;++ix) output << " " << outgoing_[ix];
  output << "\n";
  if(header) {
    output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
  }
}

string ScalarTo3ScalarDalitz::addChannel(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long id = stoi(stype);
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  ResonanceType::Type type = static_cast<ResonanceType::Type>(stoi(stype));
  // mass and width
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  Energy mass = stof(stype)*GeV;
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  Energy width = stof(stype)*GeV;
  // children
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  int d1 = stoi(stype);
  // children
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  int d2 = stoi(stype);
  // children
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  int sp = stoi(stype);
  if (sp==d1 || sp ==d2 || d1 == d2)
    return "Daughters and spectator must all be different not " + std::to_string(d1) + ", " + std::to_string(d2) + ", " + std::to_string(sp);
  // magnitude and phase
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double mag = stof(stype);
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double phi = stof(stype);
  // radius
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  InvEnergy r = stof(stype)/GeV;
  // add to list
  resonances_.push_back(DalitzResonance(id,type,mass,width,d1,d2,sp,mag,phi,r));
  // success
  return "";
}

string ScalarTo3ScalarDalitz::setExternal(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long id = stoi(stype);
  tPDPtr in = getParticleData(id);
  if(!in)
    return "Incoming particle with id " + std::to_string(id) + "does not exist";
  if(in->iSpin()!=PDT::Spin0)
    return "Incoming particle with id " + std::to_string(id) + "does not have spin 0";
  tPDVector out;
  for(unsigned int ix=0;ix<3;++ix) {
    string stype = StringUtils::car(arg);
    arg          = StringUtils::cdr(arg);
    long in = stoi(stype);
    tPDPtr pData = getParticleData(in);
    if(!pData)
      return "Outgoing particle with id " + std::to_string(in) + "does not exist";
    if(pData->iSpin()!=PDT::Spin0)
      return "Outgoing particle with id " + std::to_string(in) + "does not have spin 0";
    out.push_back(pData);
  }
  incoming_ = in->id();
  outgoing_ = {out[0]->id(),out[1]->id(),out[2]->id()};
  // success
  return "";
}

int ScalarTo3ScalarDalitz::modeNumber(bool & cc,tcPDPtr parent,
				const tPDVector & children) const {
  // must be three decay products
  if(children.size()!=3) return -1;
  // ids of the outgoing particles
  map<long,unsigned int> ids;
  for(unsigned int iy=0;iy<3;++iy) {
    if(ids.find(outgoing_[iy])!=ids.end())
      ids[outgoing_[iy]]+=1;
    else
      ids[outgoing_[iy]]+=1;
  }
  if(incoming_==parent->id()) {
    bool found=true;
    map<long,unsigned int> ids2;
    for(unsigned int iy=0;iy<children.size();++iy) {
      if(ids2.find(children[iy]->id())!=ids2.end())
	ids2[children[iy]->id()]+=1;
      else
	ids2[children[iy]->id()] =1;
    }
    for (const auto& kv : ids2) {
      if(ids.find(kv.first)==ids.end() ||
	 ids[kv.first]!=kv.second) {
	found=false;
	break;
      }
    }
    if (found) {
      cc = false;
      return 0;
    }
  }
  if( (!parent->CC() && incoming_==parent->id()) ||
      ( parent->CC() && incoming_==parent->CC()->id()) ) {
    map<long,unsigned int> ids2;
    for(unsigned int iy=0;iy<3;++iy) {
      tPDPtr part = children[iy]->CC() ? children[iy]->CC() : children[iy];
      if(ids2.find(part->id())!=ids2.end())
	ids2[part->id()]+=1;
      else
	ids2[part->id()]+=1;
    }
    bool found=true;
    for (const auto& kv : ids2) {
      if(ids.find(kv.first)==ids.end() ||
	 ids[kv.first]!=kv.second) {
	found=false;
	break;
      }
    }
    if (found) {
      cc = true;
      return 0;
    }
  }
  return -1;
}
