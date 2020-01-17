// -*- C++ -*-
//
// SemiLeptonicScalarDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SemiLeptonicScalarDecayer class.
//

#include "SemiLeptonicScalarDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"
#include "ThePEG/Helicity/epsilon.h"
#include "ThePEG/Helicity/LorentzTensor.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

SemiLeptonicScalarDecayer::SemiLeptonicScalarDecayer() {
  // intermediates
  generateIntermediates(true);
}

void SemiLeptonicScalarDecayer::doinitrun() {
  _current->initrun();
  _form->initrun();
  DecayIntegrator::doinitrun();
  if(initialize()) {
    _maxwgt.clear();
    for(unsigned int ix=0;ix<numberModes();++ix) {
      _maxwgt.push_back(mode(ix)->maxWeight());
    }
  }
}

void SemiLeptonicScalarDecayer::doinit() {
  DecayIntegrator::doinit();
  // make sure the current got initialised
  _current->init();
  // and the form factors
  _form->init();
  // the channels
  tPDVector extpart,ptemp;
  _modemap.clear();
  double maxweight;
  vector<double> channelwgts(1,1.);
  int id0(0),id1(0),Wcharge(0);
  Energy min;
  DecayPhaseSpaceModePtr mode;
  DecayPhaseSpaceChannelPtr channel;
  int iq,ia; unsigned int ix,iy,iz;
  bool done;
  for(ix=0;ix<_form->numberOfFactors();++ix) {
    // get the external particles for this mode
    extpart.resize(2);
    _form->particleID(ix,id0,id1);
    extpart[0]=getParticleData(id0);
    extpart[1]=getParticleData(id1);
    _modemap.push_back(numberModes());
    if(!extpart[0]||!extpart[1]) continue;
    Wcharge =(extpart[0]->iCharge()-extpart[1]->iCharge());
    min = extpart[0]->mass()+extpart[0]->widthUpCut()
      -extpart[1]->mass()+extpart[1]->widthLoCut();
    for(iy=0;iy<_current->numberOfModes();++iy) {
      extpart.resize(2); 	
      _current->decayModeInfo(iy,iq,ia);
      ptemp=_current->particles(Wcharge,iy,iq,ia);
      for(iz=0;iz<ptemp.size();++iz) {
	extpart.push_back(ptemp[iz]);
      }
      // create the mode
      mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
      // create the first piece of the channel
      channel = new_ptr(DecayPhaseSpaceChannel(mode));
      channel->addIntermediate(extpart[0],0,0.0,1,-1);
      done=_current->createMode(Wcharge,iy,mode,2,1,channel,min);
      if(done) {
	// the maximum weight
	if(_maxwgt.size()>numberModes()) maxweight=_maxwgt[numberModes()];
	else                             maxweight=2.;
	channelwgts.resize(mode->numberChannels(),
			   1./(mode->numberChannels()));
	addMode(mode,maxweight,channelwgts);
      }
    }
  }
}

bool SemiLeptonicScalarDecayer::accept(tcPDPtr parent, 
				       const tPDVector & children) const {
  // find the non-lepton
  int imes(0),idtemp,idin(parent->id());
  vector<int> idother; bool dummy;
  tPDVector::const_iterator pit  = children.begin();
  tPDVector::const_iterator pend = children.end();
  for( ; pit!=pend;++pit) {
    idtemp=(**pit).id();
    if(abs(idtemp)>16) imes=idtemp;
    else               idother.push_back(idtemp);
  }
  // check that the form factor exists
  if(_form->formFactorNumber(idin,imes,dummy)<0) return false;
  // and the current
  return _current->accept(idother);
}

int  SemiLeptonicScalarDecayer::modeNumber(bool & cc,tcPDPtr parent,
					   const tPDVector & children) const {
  // find the ids of the particles for the decay current
  tPDVector::const_iterator pit  = children.begin();
  tPDVector::const_iterator pend = children.end();
  int idtemp,imes(0),idin(parent->id());
  vector<int> idother;
  cc=false;
  for( ; pit!=pend;++pit) {
    idtemp=(**pit).id();
    if(abs(idtemp)>16) imes=idtemp;
    else               idother.push_back(idtemp);
  }
  return _modemap[_form->formFactorNumber(idin,imes,cc)]
    +_current->decayMode(idother);  
}


void SemiLeptonicScalarDecayer::persistentOutput(PersistentOStream & os) const {
  os << _current << _form << _maxwgt << _modemap;
}

void SemiLeptonicScalarDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _current >> _form >> _maxwgt >> _modemap;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SemiLeptonicScalarDecayer,DecayIntegrator>
describeHerwigSemiLeptonicScalarDecayer("Herwig::SemiLeptonicScalarDecayer", "HwSMDecay.so");

void SemiLeptonicScalarDecayer::Init() {

  static ClassDocumentation<SemiLeptonicScalarDecayer> documentation
    ("The SemiLeptonicScalarDecayer class is designed for the"
    "semi-leptonic decay of a (pseudo)-scalar meson.");

  static Reference<SemiLeptonicScalarDecayer,LeptonNeutrinoCurrent> interfaceCurrent
    ("Current",
     "The current for the leptons produced in the decay.",
     &SemiLeptonicScalarDecayer::_current, true, true, true, false, false);

  static Reference<SemiLeptonicScalarDecayer,ScalarFormFactor> interfaceFormFactor
    ("FormFactor",
     "The form factor",
     &SemiLeptonicScalarDecayer::_form, true, true, true, false, false);

  static ParVector<SemiLeptonicScalarDecayer,double> interfaceMaximumWeight
    ("MaximumWeight",
     "The maximum weights for the decays",
     &SemiLeptonicScalarDecayer::_maxwgt,
     0, 0, 0, 0, 100., false, false, true);

}

// combine the currents and form-factors to give the matrix element
double SemiLeptonicScalarDecayer::me2(const int ichan,
				      const Particle & inpart,
				      const ParticleVector & decay,
				      MEOption meopt) const {
  // get the information on the form-factor
  int jspin(0),id0(inpart.id()),id1(decay[0]->id());
  bool cc(false);
  unsigned int iloc(_form->formFactorNumber(id0,id1,cc));
  int spect,iq,ia;
  _form->formFactorInfo(iloc,jspin,spect,iq,ia);
  // extract leptons for the lepton current
  Energy scale;
  ParticleVector leptons;
  leptons.push_back(decay[decay.size()-2]);
  leptons.push_back(decay[decay.size()-1]);
  int mode=(abs(decay[1]->id())-11)/2;
  if(!ME()) {
    if(jspin==0)
      ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin0,PDT::Spin1Half,PDT::Spin1Half)));
    else if(jspin==1)       
      ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half)));
    else if(jspin==2)       
      ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin2,PDT::Spin1Half,PDT::Spin1Half)));
  }
  // initialisation
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&inpart),incoming);
    // work out the mapping for the lepton vector
    _constants.resize(decay.size()+1);
    _ispin.resize(decay.size());
    _imes=0;
    unsigned int itemp(1);
    for(int ix=int(decay.size()-1);ix>=0;--ix) {
      _ispin[ix]=decay[ix]->data().iSpin();
      if(abs(decay[ix]->id())<=16) {
	itemp*=_ispin[ix];
	_constants[ix]=itemp;
      }
      else _imes=ix;
    }
    _constants[decay.size()]=1;
    _constants[_imes]=_constants[_imes+1];
  }
  if(meopt==Terminate) {
    // set up the spin information for the decay products
    ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&inpart),
					  incoming,true);
    if(jspin==0)
      ScalarWaveFunction::
	constructSpinInfo(decay[0],outgoing,true);
    else if(jspin==1)
      VectorWaveFunction::
	constructSpinInfo(_vectors,decay[0],outgoing,true,false);
    else if(jspin==2)
      TensorWaveFunction::
	constructSpinInfo(_tensors,decay[0],outgoing,true,false);
    _current->current(mode,ichan,scale,leptons,meopt);
    return 0.;
  }
  // get the wavefunctions of the decay products
  switch(decay[0]->dataPtr()->iSpin()) {
  case 1:
    break;
  case 3:
    VectorWaveFunction::
      calculateWaveFunctions(_vectors,decay[0],outgoing,false);
    break;
  case 5:
    TensorWaveFunction::
      calculateWaveFunctions(_tensors,decay[0],outgoing,false);
    break;
  default:
    assert(false);
  }
  // work out the value of q and calculate the form factors
  Lorentz5Momentum q(inpart.momentum()-decay[0]->momentum());
  q.rescaleMass();
  Energy2 q2(q.mass2());
  Lorentz5Momentum sum(inpart.momentum()+decay[0]->momentum());
  // calculate the hadronic current for the decay
  Complex ii(0.,1.);
  vector<LorentzPolarizationVectorE> hadron;
  if(jspin==0) {
    Complex fp,f0;
    _form->ScalarScalarFormFactor(q2,iloc,id0,id1,inpart.mass(),decay[0]->mass(),
				  f0,fp);
    Complex pre((sqr(inpart.mass())-sqr(decay[0]->mass()))/q2*(f0-fp));
    hadron.push_back(fp*sum+(pre*q));
  }
  else if(jspin==1) {
    Complex A0,A1,A2,A3,V;
    complex<Energy> dot;
    Energy MP(inpart.mass()),MV(decay[0]->mass()),msum(MP+MV),mdiff(MP-MV);
    _form->ScalarVectorFormFactor(q2,iloc,id0,id1,MP,MV,A0,A1,A2,V);
    A3 = Complex(0.5/MV*(msum*A1-mdiff*A2));
    if(cc) V*=-1.;
    // compute the hadron currents
    for(unsigned int ix=0;ix<3;++ix) {
      // dot product
      dot = _vectors[ix]*inpart.momentum();
      // current
      hadron.push_back(-ii*msum*A1*_vectors[ix]
		       +ii*A2/msum*dot*sum
		       +2.*ii*MV/q2*(A3-A0)*dot*q
		       +2.*V/msum*Helicity::epsilon(_vectors[ix],inpart.momentum(),
						    decay[0]->momentum()));
    }
  }
  else if(jspin==2) {
    complex<InvEnergy2> h,bp,bm;
    complex<double> k;
    complex<Energy2> dot;
    _form->ScalarTensorFormFactor(q2,iloc,id0,id1,inpart.mass(),decay[0]->mass(),
				  h,k,bp,bm);
    if(!cc) h*=-1.;
    LorentzPolarizationVectorE dotv;
    // compute the hadron currents
    for(unsigned int ix=0;ix<5;++ix) {
      dotv = _tensors[ix]*inpart.momentum();
      dot = dotv*inpart.momentum();
      hadron.push_back(ii*h*Helicity::epsilon(dotv,sum,q)
		       -k*dotv-bp*dot*sum-bm*dot*q);
    }
  }
  // construct the lepton current
  vector<LorentzPolarizationVectorE> 
    lepton(_current->current(mode,ichan,scale,leptons,meopt));
  // compute the matrix element
  vector<unsigned int> ihel(decay.size()+1);
  for(unsigned int mhel=0;mhel<hadron.size();++mhel) {
    for(unsigned int lhel=0;lhel<lepton.size();++lhel) {
      // map the index for the leptons to a helicity state
      for(unsigned int ix=decay.size();ix>0;--ix) {
	if(ix-1!=_imes) ihel[ix]=(lhel%_constants[ix-1])/_constants[ix];
      }
      // helicities of mesons
      ihel[0]=0;
      ihel[_imes+1]=mhel;
      (*ME())(ihel) = Complex(lepton[lhel].dot(hadron[mhel])*SM().fermiConstant());
    }
  }
  // store the matrix element
  double ckm(1.);
  if(iq<=6) {
    if(iq%2==0) ckm = SM().CKM(abs(iq)/2-1,(abs(ia)-1)/2);
    else        ckm = SM().CKM(abs(ia)/2-1,(abs(iq)-1)/2);
  }
  // return the answer
  return 0.5*(ME()->contract(_rho)).real()*ckm; 
}
 
// output the setup information for the particle database
void SemiLeptonicScalarDecayer::dataBaseOutput(ofstream & output,
					       bool header) const {
  if(header) output << "update decayers set parameters=\"";
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<_maxwgt.size();++ix) {
    output << "insert " << name() << ":MaximumWeight " << ix << " " 
	   << _maxwgt[ix] << "\n";
  }
  _current->dataBaseOutput(output,false,true);
  output << "newdef " << name() << ":Current " << _current->name() << " \n";
  _form->dataBaseOutput(output,false,true);
  output << "newdef " << name() << ":FormFactor " << _form->name() << " \n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}
