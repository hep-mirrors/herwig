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
#include "FlatteResonance.h"
#include "MIPWA.h"
#include "PiPiI2.h"
#include "DalitzKMatrix.h"
#include "DalitzLASS.h"
#include "DalitzGS.h"

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
  // special for flate
  if (type==ResonanceType::Flattef0 ||
      type==ResonanceType::Flattea0 ||
      type==ResonanceType::FlatteKstar0) {
    // Flatte parameters
    // magnitude and phase
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    double fpi = stof(stype);
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    double fK  = stof(stype);
    // add to list
    resonances_.push_back(new_ptr(FlatteResonance(id,type,mass,width,d1,d2,sp,mag,phi,r,fpi,fK)));
  }
  // MIPWA
  else if(type==ResonanceType::Spin0MIPWA) {
    // no of entries in table
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    int nn = stoi(stype);
    vector<Energy> en; en.reserve(nn);
    vector<double> mag2,phase2; mag2.reserve(nn); phase2.reserve(nn);
    for(int ix=0;ix<nn;++ix) {
      stype = StringUtils::car(arg);
      arg   = StringUtils::cdr(arg);
      Energy ee = stof(stype)*GeV;
      en.push_back(ee);
      stype = StringUtils::car(arg);
      arg   = StringUtils::cdr(arg);
      double mm=stof(stype);
      mag2.push_back(mm);
      stype = StringUtils::car(arg);
      arg   = StringUtils::cdr(arg);
      double pp=stof(stype);
      phase2.push_back(pp);
    }
    resonances_.push_back(new_ptr(MIPWA(id,type,mass,width,d1,d2,sp,mag,phi,r,en,mag2,phase2)));
  }
  // I=2 pipi
  else if(type==ResonanceType::PiPiI2) {
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    InvEnergy a = stof(stype)/GeV;
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    InvEnergy2 b = stof(stype)/GeV2;
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    InvEnergy4 c = stof(stype)/GeV2/GeV2;
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    InvEnergy6 d = stof(stype)/GeV2/GeV2/GeV2;
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    Energy mmin = stof(stype)*GeV;
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    Energy mmax = stof(stype)*GeV;
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    double dEta=stof(stype);
    resonances_.push_back(new_ptr(PiPiI2(id,type,mass,width,d1,d2,sp,mag,phi,r,
					 a,b,c,d,mmin,mmax,dEta)));
  }
  // K-matrix
  else if(type==ResonanceType::KMatrix) {
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    unsigned int imat = stoi(stype);
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    unsigned int chan = stoi(stype);
    assert(imat<kMatrix_.size());
    // expansion point for the constants terms
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    Energy2 sc = GeV2*stof(stype);
    // type of expansion
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    unsigned int itype= stoi(stype);
    vector<pair<double,double> > beta;
    // first loop over the coefficients of the poles
    for(unsigned int ix=0;ix<kMatrix_[imat]->poles().size();++ix) {
      stype = StringUtils::car(arg);
      arg   = StringUtils::cdr(arg);
      double b = stof(stype);
      stype = StringUtils::car(arg);
      arg   = StringUtils::cdr(arg);
      beta.push_back(make_pair(b,stof(stype)));
    }
    // now over the power series for the different channels
    vector<pair<double,vector<double > > > coeffs(kMatrix_[imat]->numberOfChannels());
    for(unsigned int ix=0;ix<kMatrix_[imat]->numberOfChannels();++ix) {
      stype = StringUtils::car(arg);
      arg   = StringUtils::cdr(arg);
      unsigned int nterms = stoi(stype);
      for(unsigned int iy=0;iy<nterms;++iy) {
	stype = StringUtils::car(arg);
	arg   = StringUtils::cdr(arg);
	coeffs[ix].second.push_back(stof(stype));
      }
      stype = StringUtils::car(arg);
      arg   = StringUtils::cdr(arg);
      coeffs[ix].first = stof(stype);
    }
    // finally make the channel
    resonances_.push_back(new_ptr(DalitzKMatrix(id,type,mass,width,d1,d2,sp,mag,phi,r,imat,chan,sc,itype,beta,coeffs)));
  }
  // LASS
  else if(type==ResonanceType::LASS) {
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    unsigned int iopt = stoi(stype);
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    double FNR = stof(stype);
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    double phiNR = stof(stype);
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    double FRes = stof(stype);
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    double phiRes = stof(stype);
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    InvEnergy ascat = stof(stype)/GeV;
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    InvEnergy reff = stof(stype)/GeV;
    // finally make the channel
    resonances_.push_back(new_ptr(DalitzLASS(id,type,mass,width,d1,d2,sp,mag,phi,r,iopt,
					     FNR,phiNR,FRes,phiRes,ascat,reff)));
  }
  // GS form
  else if(type==ResonanceType::Spin1GS) {
    resonances_.push_back(new_ptr(DalitzGS(id,type,mass,width,d1,d2,sp,mag,phi,r)));
  }
  // otherwise add to list
  else {
    resonances_.push_back(new_ptr(DalitzResonance(id,type,mass,width,d1,d2,sp,mag,phi,r)));
  }
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
