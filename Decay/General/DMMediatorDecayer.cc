// -*- C++ -*-
//
// DMMediatorDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DMMediatorDecayer class.
//

#include "DMMediatorDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;

IBPtr DMMediatorDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr DMMediatorDecayer::fullclone() const {
  return new_ptr(*this);
}

void DMMediatorDecayer::persistentOutput(PersistentOStream & os) const {
  os << inpart_ << currentOut_ << current_ << mode_ << wgtloc_ << wgtmax_ << weights_ << cDMmed_ << cSMmed_;
}

void DMMediatorDecayer::persistentInput(PersistentIStream & is, int) {
  is >> inpart_ >> currentOut_ >> current_ >> mode_ >> wgtloc_ >> wgtmax_ >> weights_ >> cDMmed_ >> cSMmed_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<DMMediatorDecayer,DecayIntegrator>
describeHerwigDMMediatorDecayer("Herwig::DMMediatorDecayer", "Herwig.so");

void DMMediatorDecayer::Init() {

  static ClassDocumentation<DMMediatorDecayer> documentation
    ("The DMMediatorDecayer class is designed for the decays of low mass vector bosons");
  
}

void DMMediatorDecayer::setDecayInfo(PDPtr in, const vector<tPDPtr> & outCurrent,
				     WeakCurrentPtr current, double normin, vector<double> sm) {
  inpart_ = in;
  currentOut_ = outCurrent;
  current_ = current;
  cDMmed_=normin;
  cSMmed_=sm;
}



int DMMediatorDecayer::modeNumber(bool & cc, tcPDPtr parent, 
				  const tPDVector & children) const {
  vector<long> id;
  id.push_back(parent->id());
  for(unsigned int ix=0;ix<children.size();++ix) id.push_back(children[ix]->id());
  return modeNumber(cc,id);
}

void DMMediatorDecayer::doinitrun() {
  current_->initrun();
  DecayIntegrator::doinitrun();
}

void DMMediatorDecayer::doinit() {
  DecayIntegrator::doinit();
  // make sure the current got initialised
  current_->init();
  // find the mode
  for(unsigned int ix=0;ix<current_->numberOfModes();++ix) {
    // get the external particles for this mode
     int iq(0),ia(0);
     tPDVector ptemp  = current_->particles(inpart_->iCharge(),ix,iq,ia);
     // check this is the right mode
     if(ptemp.size()!=currentOut_.size()) continue;
     vector<bool> matched(ptemp.size(),false);
     bool match = true;
     for(unsigned int iy=0;iy<currentOut_.size();++iy) {
       bool found = false;
       for(unsigned int iz=0;iz<ptemp.size();++iz) {
	 if(!matched[iz]&&ptemp[iz]==currentOut_[iy]) {
	   found = true;
	   matched[iz] = true;
	   break;
	 }
       }
       if(!found) {
	 match = false;
	 break;
       }
     }
     if(!match) continue;
     tPDVector out = {};
     out.insert(std::end(out), std::begin(ptemp), std::end(ptemp));
     // create the mode
     PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(inpart_,out,1.));
     PhaseSpaceChannel channel(mode,true);
     bool done=current_->createMode(inpart_->iCharge(),tcPDPtr(),FlavourInfo(),
				    ix,mode,0,-1,channel,inpart_->mass());
     if(done) {
       // the maximum weight and the channel weights
       // the weights for the channel
       if(weights_.empty()) {
	 weights_.resize(mode->channels().size(),1./(mode->channels().size()));
       }
       mode_ = ix;
       // special for the two body modes
       if(out.size()==2) {
	 weights_.clear();
	 mode=new_ptr(PhaseSpaceMode(inpart_,out,1.));
       }
       mode->maxWeight(wgtmax_);
       mode->setWeights(weights_);
       addMode(mode);
     }
     break;
  }
}

int DMMediatorDecayer::modeNumber(bool & cc, vector<long> id) const {
  // incoming particle
  long idtemp;
  tPDPtr p0=getParticleData(id[0]);
  idtemp = p0->CC() ? -id[0] : id[0];
  if(     id[0] ==inpart_->id()) cc=false;
  else if(idtemp==inpart_->id()) cc=true ;
  else return -1;
  vector<int> idout;
  for(vector<long>::iterator it=++id.begin();it!=id.end();++it) {
    idout.push_back(*it);
  }
  unsigned int icurr=current_->decayMode(idout);
  if(mode_==icurr) return  0;
  else             return -1;
}

void DMMediatorDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  VectorWaveFunction::constructSpinInfo(vectors_,const_ptr_cast<tPPtr>(&part),
					Helicity::incoming,true,false);
  weakCurrent()->constructSpinInfo(ParticleVector(decay.begin(),decay.end()));
}

double DMMediatorDecayer::me2(const int ichan, const Particle & part,
			      const tPDVector & outgoing,
			      const vector<Lorentz5Momentum> & momenta,
			      MEOption meopt) const {
  using namespace ThePEG::Helicity;
  // polarization vectors for the incoming particle
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_,rho_,
  					       const_ptr_cast<tPPtr>(&part),
  					       incoming,false);
    // fix rho if no correlations
    fixRho(rho_);
  }
  // work out the mapping for the hadron vector
  int nOut = momenta.size();
  vector<unsigned int> constants(nOut+1);
  vector<PDT::Spin   > iSpin(nOut);
  vector<int> hadrons(nOut);
  int itemp(1);
  int ix(nOut);
  do {
    --ix;
    iSpin[ix]      = outgoing[ix]->iSpin();
    itemp         *= iSpin[ix];
    constants[ix]  = itemp;
    hadrons[ix]    = outgoing[ix]->id();
  }
  while(ix>0);
  constants[nOut] = 1;
  Energy2 scale(sqr(part.mass()));
  // calculate the hadron current
  Energy q = part.mass();
  // currents for the different flavour components
  vector<LorentzPolarizationVectorE> 
    hadronI0(current_->current(tcPDPtr(), FlavourInfo(IsoSpin::IZero, IsoSpin::I3Zero,Strangeness::Zero),
			       mode(),ichan,q,outgoing,momenta,DecayIntegrator::Calculate));
  vector<LorentzPolarizationVectorE> 
    hadronI1(current_->current(tcPDPtr(), FlavourInfo(IsoSpin::IOne, IsoSpin::I3Zero,Strangeness::Zero),
			       mode(),ichan,q,outgoing,momenta,DecayIntegrator::Calculate));
  vector<LorentzPolarizationVectorE> 
    hadronssbar(current_->current(tcPDPtr(), FlavourInfo(IsoSpin::IZero, IsoSpin::I3Zero,Strangeness::ssbar),
				  mode(),ichan,q,outgoing,momenta,DecayIntegrator::Calculate));
  // compute the matrix element
  GeneralDecayMEPtr newME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,iSpin)));
  vector<unsigned int> ihel(momenta.size()+1);
  unsigned int hI0_size = hadronI0.size();
  unsigned int hI1_size = hadronI1.size();
  unsigned int hss_size = hadronssbar.size();
  unsigned int maxsize = max(max(hadronI0.size(),hadronI1.size()),hss_size);
  for(unsigned int hhel=0;hhel<maxsize;++hhel) {
    // map the index for the hadrons to a helicity state
    for(int ix=nOut;ix>0;--ix) {
      ihel[ix+1]=(hhel%constants[ix-1])/constants[ix];
    }
    for(ihel[0]=0;ihel[0]<3;++ihel[0]) {
      Complex amp = 0.;
      // work on coefficients for the I1 and I0 bits
      if(hI0_size != 0 )
	amp += (cSMmed_[0]+cSMmed_[1])/sqrt(2.)/q*(vectors_[ihel[0]].wave().dot(hadronI0[hhel]));
      if(hI1_size !=0)
	amp += (cSMmed_[0]-cSMmed_[1])/sqrt(2.)/q*(vectors_[ihel[0]].wave().dot(hadronI1[hhel]));
      if(hss_size !=0)
	amp += cSMmed_[2]                      /q*(vectors_[ihel[0]].wave().dot(hadronssbar[hhel]));
      (*newME)(ihel) = cDMmed_*amp;
    }
  }
  // store the matrix element
  ME(newME);
  // return the answer
  double output = (ME()->contract(rho_)).real();
  return output;
}
 
Energy DMMediatorDecayer::partialWidth(tPDPtr part, vector<tPDPtr> out) {
  vector<long> id;
  id.push_back(part->id());
  for(unsigned int ix=0;ix<out.size();++ix) id.push_back(out[ix]->id());
  bool cc;
  int mode=modeNumber(cc,id);
  imode(mode);
  return initializePhaseSpaceMode(mode,true,true);
}
