// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DalitzBase class.
//

#include "DalitzBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Decay/PhaseSpaceMode.h"

using namespace Herwig;

void DalitzBase::persistentOutput(PersistentOStream & os) const {
  os << ounit(rParent_,1./GeV) << resonances_ << maxWgt_ << weights_ 
     << channel1_ << channel2_ << incoming_ << outgoing_ << useAllK0_
     << kMatrix_;
}

void DalitzBase::persistentInput(PersistentIStream & is, int) {
  is >> iunit(rParent_,1./GeV) >> resonances_ >> maxWgt_ >> weights_ 
     >> channel1_ >> channel2_ >> incoming_ >> outgoing_ >> useAllK0_
     >> kMatrix_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<DalitzBase,DecayIntegrator>
describeHerwigDalitzBase("Herwig::DalitzBase", "HwDalitzDecay.so");

void DalitzBase::Init() {

  static ClassDocumentation<DalitzBase> documentation
    ("The DalitzBase class provides a base class for the implementation of three-body Dalitz decays.");

  static Command<DalitzBase> interfaceSetExternal
    ("SetExternal",
     "Set the external particles for the decay mode",
     &DalitzBase::setExternal, false);
  
  static Command<DalitzBase> interfaceAddChannel
    ("AddChannel",
     "Add a channel for the description of the matrix element",
     &DalitzBase::addChannel, false);

  static Parameter<DalitzBase,InvEnergy> interfaceParentRadius
    ("ParentRadius",
     "The radius parameter for the Blatt-Weisskopf form-factor for the D",
     &DalitzBase::rParent_, 1./GeV, 5./GeV, ZERO, 10./GeV,
     false, false, Interface::limited);

  static Parameter<DalitzBase,double> interfaceMaximumWeight
    ("MaximumWeight",
     "The maximum weight for the phase-space sampling",
     &DalitzBase::maxWgt_, 1.0, 0.0, 1e20,
     false, false, Interface::limited);

  static ParVector<DalitzBase,double> interfaceWeights
    ("Weights",
     "The weights for the different channels for the phase-space integration",
     &DalitzBase::weights_, -1, 1.0, 0.0, 1.0,
     false, false, Interface::limited);
  
  static Parameter<DalitzBase,int> interfaceChannel1
    ("Channel1",
     "The first allowed channel, for debugging/calculation of fit fractions only",
     &DalitzBase::channel1_, -1, -1, 100,
     false, false, Interface::limited);
  
  static Parameter<DalitzBase,int> interfaceChannel2
    ("Channel2",
     "The first allowed channel, for debugging/calculation of fit fractions only",
     &DalitzBase::channel2_, -1, -1, 100,
     false, false, Interface::limited);
  
  static Switch<DalitzBase,bool> interfaceUseAllK0
    ("UseAllK0",
     "Use all K0 mesons when matching the mode",
     &DalitzBase::useAllK0_, false, false, false);
  static SwitchOption interfaceUseAllK0No
    (interfaceUseAllK0,
     "No",
     "Just use the identified state",
     false);
  static SwitchOption interfaceUseAllK0Yes
    (interfaceUseAllK0,
     "Yes",
     "Use all the states",
     true);

  static RefVector<DalitzBase,KMatrix> interfaceKMatrices
    ("KMatrices",
     "Any K-matrices needed to simulate the decay",
     &DalitzBase::kMatrix_, -1, false, false, true, false, false);

}

void DalitzBase::doinit() {
  if(incoming_!=0) {
    tPDPtr in = getParticleData(incoming_);
    vector<tPDPtr> out = {getParticleData(outgoing_[0]),
			  getParticleData(outgoing_[1]),
			  getParticleData(outgoing_[2])};
    createMode(in,out);
  }
  DecayIntegrator::doinit();
}

void DalitzBase::doinitrun() {
  if(!kMatrix_.empty()) {
    for(unsigned int ix=0;ix<resonances().size();++ix) {
      Ptr<Herwig::DalitzKMatrix>::transient_pointer mat =
	dynamic_ptr_cast<Ptr<Herwig::DalitzKMatrix>::transient_pointer>(resonances()[ix]);
      if(mat) {
	mat->setKMatrix(kMatrix_[mat->imatrix()]);
	// Energy2 s=3.*GeV2;
	// Complex amp = resonances()[ix]->BreitWigner(sqrt(s),0.139*GeV,0.139*GeV);
	// Energy2 s=0.5*GeV2;
	// while(s<3.5*GeV2) {
	//   Complex amp = resonances()[ix]->BreitWigner(sqrt(s),0.139*GeV,0.139*GeV);
	//   cerr << s/GeV2 << " " << abs(amp) << "\n";
	//   s+=0.01*GeV2;
	// }
      }
    }
  }
  DecayIntegrator::doinitrun();
  weights_.resize(mode(0)->channels().size());
  maxWgt_ = mode(0)->maxWeight();
  for(unsigned int iz=0;iz<mode(0)->channels().size();++iz) {
    weights_[iz]=mode(0)->channels()[iz].weight();
  }
}

void DalitzBase::createMode(tPDPtr in, tPDVector out) {
  PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(in,out,maxWgt_));
  if(weights_.size()!=resonances().size()) {
    weights_=vector<double>(resonances_.size(),1./double(resonances_.size()));
  }
  unsigned int ix=0;
  for(DalitzResonancePtr res : resonances_) {
    tPDPtr resonance = getParticleData(res->id);
    if(resonance) {
      mode->addChannel((PhaseSpaceChannel(mode),0,resonance,0,res->spectator+1,1,res->daughter1+1,1,res->daughter2+1));
      resetIntermediate(resonance,res->mass,abs(res->width));
      ++ix;
    }
  }
  addMode(mode);
}

void DalitzBase::dataBaseOutput(ofstream & output, bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DalitzBase base class
  DecayIntegrator::dataBaseOutput(output,false);
  output << "newdef " << name() << ":ParentRadius " << rParent_*GeV << "\n";
  output << "newdef " << name() << ":UseAllK0 " << useAllK0_ << "\n";
  output << "newdef " << name() << ":MaximumWeight " << maxWgt_ << "\n";
  for(unsigned int ix=0;ix<kMatrix_.size();++ix) {
    output << "insert " << name() << ":KMatrices " << ix << " " << kMatrix_[ix]->fullName() << "\n";
  }
  for(unsigned int ix=0;ix<weights_.size();++ix) {
    output << "insert " << name() << ":Weights "
	   << ix << " " << weights_[ix] << "\n";
  }
  for(unsigned int ix=0;ix<resonances_.size();++ix) {
    output << "do " << name() << ":AddChannel ";
    resonances_[ix]->dataBaseOutput(output);
    output << "\n";
  }
  output << "do " << name() << ":SetExternal " << incoming_;
  for(unsigned int ix=0;ix<3;++ix) output << " " << outgoing_[ix];
  output << "\n";
  if(header) {
    output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
  }
}

string DalitzBase::addChannel(string arg) {
  string error;
  DalitzResonancePtr res = DalitzResonance::readResonance(arg,error);
  if (res)
    resonances_.push_back(res);
  else
    return error;
  // success
  return "";
}

string DalitzBase::setExternal(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long id = stoi(stype);
  tPDPtr in = getParticleData(id);
  if(!in)
    return "Incoming particle with id " + std::to_string(id) + "does not exist";
  tPDVector out;
  for(unsigned int ix=0;ix<3;++ix) {
    string stype = StringUtils::car(arg);
    arg          = StringUtils::cdr(arg);
    long in = stoi(stype);
    tPDPtr pData = getParticleData(in);
    if(!pData)
      return "Outgoing particle with id " + std::to_string(in) + "does not exist";
    out.push_back(pData);
  }
  incoming_ = in->id();
  outgoing_ = {out[0]->id(),out[1]->id(),out[2]->id()};
  // success
  return "";
}

int DalitzBase::modeNumber(bool & cc,tcPDPtr parent,
				const tPDVector & children) const {
  // must be three decay products
  if(children.size()!=3) return -1;
  // ids of the outgoing particles
  map<long,unsigned int> ids;
  for(unsigned int iy=0;iy<3;++iy) {
    long id = outgoing_[iy];
    if(useAllK0_ && (id==130 || id==310 || id==311 || id==-311)) id=130; 
    if(ids.find(id)!=ids.end())
      ids[id]+=1;
    else
      ids[id]+=1;
  }
  if(incoming_==parent->id()) {
    bool found=true;
    map<long,unsigned int> ids2;
    for(unsigned int iy=0;iy<children.size();++iy) {
    long id = children[iy]->id();
    if(useAllK0_ && (id==130 || id==310 || id==311 || id==-311)) id=130; 
      if(ids2.find(id)!=ids2.end())
	ids2[id]+=1;
      else
	ids2[id] =1;
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
      long id = part->id();
      if(useAllK0_ && (id==130 || id==310 || id==311 || id==-311)) id=130; 
      if(ids2.find(id)!=ids2.end())
	ids2[id]+=1;
      else
	ids2[id]+=1;
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
