// -*- C++ -*-
//
// SemiLeptonicScalarDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
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
#include "ThePEG/Helicity/HelicityFunctions.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

SemiLeptonicScalarDecayer::SemiLeptonicScalarDecayer() {
  // intermediates
  generateIntermediates(true);
}

void SemiLeptonicScalarDecayer::doinitrun() {
  _current->initrun();
  _form->initrun();
  DecayIntegrator2::doinitrun();
  if(initialize()) {
    _maxwgt.clear();
    for(unsigned int ix=0;ix<numberModes();++ix) {
      _maxwgt.push_back(mode(ix)->maxWeight());
    }
  }
}

void SemiLeptonicScalarDecayer::doinit() {
  DecayIntegrator2::doinit();
  // make sure the current got initialised
  _current->init();
  // and the form factors
  _form->init();
  _modemap.clear();
  for(unsigned int ix=0;ix<_form->numberOfFactors();++ix) {
    // get the external particles for this mode
    int id0(0),id1(0);
    _form->particleID(ix,id0,id1);
    tPDPtr  in = getParticleData(id0);
    tPDPtr out = getParticleData(id1);
    _modemap.push_back(numberModes());
    if(!in || !out) continue;
    int Wcharge =(in->iCharge()-out->iCharge());
    Energy min = in->mass()+in->widthUpCut()
      -out->mass()+out->widthLoCut();
    for(unsigned int iy=0;iy<_current->numberOfModes();++iy) {
      int iq(0),ia(0);
      _current->decayModeInfo(iy,iq,ia);
      tPDVector outV = {out};
      tPDVector ptemp=_current->particles(Wcharge,iy,iq,ia);
      outV.insert(std::end(outV), std::begin(ptemp), std::end(ptemp));
      // create the mode
      PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(in,outV,1.));
      // create the first piece of the channel
      PhaseSpaceChannel channel((PhaseSpaceChannel(mode),0,1));
      // and the rest
      bool done = _current->createMode(Wcharge,tcPDPtr(),IsoSpin::IUnknown,IsoSpin::I3Unknown,
				       iy,mode,1,0,channel,min);
      if(done) {
	// the maximum weight
	double maxweight = _maxwgt.size()>numberModes() ? _maxwgt[numberModes()] : 2.;
	mode->maxWeight(maxweight);
	addMode(mode);
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
DescribeClass<SemiLeptonicScalarDecayer,DecayIntegrator2>
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

void SemiLeptonicScalarDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // set up the spin information for the decay products
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					incoming,true);
  if(decay[0]->dataPtr()->iSpin()==PDT::Spin0)
    ScalarWaveFunction::
      constructSpinInfo(decay[0],outgoing,true);
  else if(decay[0]->dataPtr()->iSpin()==PDT::Spin1)
    VectorWaveFunction::
      constructSpinInfo(_vectors,decay[0],outgoing,true,false);
  else if(decay[0]->dataPtr()->iSpin()==PDT::Spin2)
    TensorWaveFunction::
      constructSpinInfo(_tensors,decay[0],outgoing,true,false);
  // and the stuff from the current
  _current->constructSpinInfo(ParticleVector(decay.begin()+1,decay.end()));
}

// combine the currents and form-factors to give the matrix element
double SemiLeptonicScalarDecayer::me2(const int , const Particle & part,
				      const tPDVector & outgoing,
				      const vector<Lorentz5Momentum> & momenta,
				      MEOption meopt) const {
  // get the information on the form-factor
  int jspin(0),id0(part.id()),id1(outgoing[0]->id());
  bool cc(false);
  unsigned int iloc(_form->formFactorNumber(id0,id1,cc));
  int spect,iq,ia;
  _form->formFactorInfo(iloc,jspin,spect,iq,ia);
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
      calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&part),incoming);
    // work out the mapping for the lepton vector
    _constants.resize(outgoing.size()+1);
    _ispin.resize(outgoing.size());
    _imes=0;
    unsigned int itemp(1);
    for(int ix=int(outgoing.size()-1);ix>=0;--ix) {
      _ispin[ix]=outgoing[ix]->iSpin();
      if(abs(outgoing[ix]->id())<=16) {
  	itemp*=_ispin[ix];
  	_constants[ix]=itemp;
      }
      else _imes=ix;
    }
    _constants[outgoing.size()]=1;
    _constants[_imes]=_constants[_imes+1];
  }
  // get the wavefunctions of the decay products
  switch(outgoing[0]->iSpin()) {
  case PDT::Spin0:
    break;
  case PDT::Spin1:
    _vectors.resize(3);
    for(unsigned int ihel=0;ihel<3;++ihel)
      _vectors[ihel] = HelicityFunctions::polarizationVector(-momenta[0],ihel,
							     Helicity::outgoing);
    break;
  case PDT::Spin2:
    {
      TensorWaveFunction twave(momenta[0],outgoing[0],Helicity::outgoing);
      _tensors.resize(5);
      for(unsigned int ihel=0;ihel<5;++ihel) {
	twave.reset(ihel);
	_tensors[ihel] = twave.wave();
      }
    }
    break;
  default:
    assert(false);
  }
  // work out the value of q and calculate the form factors
  Lorentz5Momentum q(part.momentum()-momenta[0]);
  q.rescaleMass();
  Energy2 q2(q.mass2());
  Lorentz5Momentum sum(part.momentum()+momenta[0]);
  // calculate the hadronic current for the decay
  Complex ii(0.,1.);
  vector<LorentzPolarizationVectorE> hadron;
  if(jspin==0) {
    Complex fp,f0;
    _form->ScalarScalarFormFactor(q2,iloc,id0,id1,part.mass(),momenta[0].mass(),
  				  f0,fp);
    Complex pre((sqr(part.mass())-sqr(momenta[0].mass()))/q2*(f0-fp));
    hadron.push_back(fp*sum+(pre*q));
  }
  else if(jspin==1) {
    Complex A0,A1,A2,A3,V;
    complex<Energy> dot;
    Energy MP(part.mass()),MV(momenta[0].mass()),msum(MP+MV),mdiff(MP-MV);
    _form->ScalarVectorFormFactor(q2,iloc,id0,id1,MP,MV,A0,A1,A2,V);
    A3 = 0.5/MV*(msum*A1-mdiff*A2);
    if(cc) V*=-1.;
    // compute the hadron currents
    for(unsigned int ix=0;ix<3;++ix) {
      // dot product
      dot = _vectors[ix]*part.momentum();
      // current
      hadron.push_back(-ii*msum*A1*_vectors[ix]
  		       +ii*A2/msum*dot*sum
  		       +2.*ii*MV/q2*(A3-A0)*dot*q
  		       +2.*V/msum*Helicity::epsilon(_vectors[ix],part.momentum(),
  						    momenta[0]));
    }
  }
  else if(jspin==2) {
    complex<InvEnergy2> h,bp,bm;
    complex<double> k;
    complex<Energy2> dot;
    _form->ScalarTensorFormFactor(q2,iloc,id0,id1,part.mass(),momenta[0].mass(),
  				  h,k,bp,bm);
    if(!cc) h*=-1.;
    LorentzPolarizationVectorE dotv;
    // compute the hadron currents
    for(unsigned int ix=0;ix<5;++ix) {
      dotv = _tensors[ix]*part.momentum();
      dot = dotv*part.momentum();
      hadron.push_back(ii*h*Helicity::epsilon(dotv,sum,q)
  		       -k*dotv-bp*dot*sum-bm*dot*q);
    }
  }
  Energy scale;
  int mode=(abs(outgoing[1]->id())-11)/2;
  vector<LorentzPolarizationVectorE> 
    lepton(_current->current(tcPDPtr(),IsoSpin::IUnknown,IsoSpin::I3Unknown,
			     mode,-1,scale,tPDVector(outgoing.begin()+1,outgoing.end()),
			     vector<Lorentz5Momentum>(momenta.begin()+1,momenta.end()),
			     meopt));
  // compute the matrix element
  vector<unsigned int> ihel(outgoing.size()+1);
  for(unsigned int mhel=0;mhel<hadron.size();++mhel) {
    for(unsigned int lhel=0;lhel<lepton.size();++lhel) {
      // map the index for the leptons to a helicity state
      for(unsigned int ix=outgoing.size();ix>0;--ix) {
  	if(ix-1!=_imes) ihel[ix]=(lhel%_constants[ix-1])/_constants[ix];
      }
      // helicities of mesons
      ihel[0]=0;
      ihel[_imes+1]=mhel;
      (*ME())(ihel)= lepton[lhel].dot(hadron[mhel])*SM().fermiConstant();
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
  DecayIntegrator2::dataBaseOutput(output,false);
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
