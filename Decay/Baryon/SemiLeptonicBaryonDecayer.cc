// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SemiLeptonicBaryonDecayer class.
//

#include "SemiLeptonicBaryonDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorBarWaveFunction.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/Helicity/RSFermionSpinInfo.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
#include "ThePEG/Helicity/HelicityFunctions.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

SemiLeptonicBaryonDecayer::SemiLeptonicBaryonDecayer() {
  // intermediates
  generateIntermediates(true);
}

void SemiLeptonicBaryonDecayer::doinitrun() {
  _current->initrun();
  _form->initrun();
  DecayIntegrator::doinitrun();
  if(initialize()) {
    _maxwgt.clear();
    for(unsigned int ix=0;ix<numberModes();++ix)
      _maxwgt.push_back(mode(ix)->maxWeight());
  }
}

void SemiLeptonicBaryonDecayer::doinit() {
  DecayIntegrator::doinit();
  // make sure the current got initialised
  _current->init();
  // and the form factors
  _form->init();
  // the channels
  _modemap.clear();
  for(unsigned int ix=0;ix<_form->numberOfFactors();++ix) {
    int id0(0),id1(0);
    // get the external particles for this mode
    _form->particleID(ix,id0,id1);
    int inspin,spect1,spect2,inquark,outquark,outspin;
    _form->formFactorInfo(ix,inspin,outspin,spect1,spect2,inquark,outquark);
    // incoming  and outgoing particles
    tPDPtr in  = getParticleData(id0);
    tPDPtr out = getParticleData(id1);
    // charge of the W
    int Wcharge =(in->iCharge()-out->iCharge());
    Energy min = in->mass()+in->widthUpCut()
      -out->mass()+out->widthLoCut();
    _modemap.push_back(numberModes());
    for(unsigned int iy=0;iy<_current->numberOfModes();++iy) {
      int iq(0),ia(0);
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
      // check the result
      if(done&&abs(Wcharge)==3&&inspin==2&&(outspin==2||outspin==4)) {
	// the maximum weight
	double maxweight = _maxwgt.size()>numberModes() ? _maxwgt[numberModes()] : 2.;
	mode->maxWeight(maxweight);
	addMode(mode);
      }
    }
  }
}

bool SemiLeptonicBaryonDecayer::accept(tcPDPtr parent,
				       const tPDVector & children) const {
  // find the non-lepton
  int ibar(0),idtemp,idin(parent->id());
  vector<int> idother; bool dummy;
  tPDVector::const_iterator pit  = children.begin();
  tPDVector::const_iterator pend = children.end();
  for( ; pit!=pend;++pit) {
    idtemp=(**pit).id();
    if(abs(idtemp)>16) ibar=idtemp;
    else idother.push_back(idtemp);
  }
  // check that the form factor exists
  if(_form->formFactorNumber(idin,ibar,dummy)<0) return false;
  // and the current
  return _current->accept(idother);
}

int SemiLeptonicBaryonDecayer::modeNumber(bool & cc,tcPDPtr parent,
					  const tPDVector & children) const {
  // find the ids of the particles for the decay current
  tPDVector::const_iterator pit  = children.begin();
  tPDVector::const_iterator pend = children.end();
  int idtemp,ibar(0),idin(parent->id());
  vector<int> idother;
  cc=false;
  for( ; pit!=pend;++pit) {
    idtemp=(**pit).id();
    if(abs(idtemp)>16) ibar=idtemp;
    else idother.push_back(idtemp);
  }
  return _modemap[_form->formFactorNumber(idin,ibar,cc)]
    +_current->decayMode(idother);
}


void SemiLeptonicBaryonDecayer::persistentOutput(PersistentOStream & os) const {
  os << _current << _form << _maxwgt << _modemap;
}

void SemiLeptonicBaryonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _current >> _form >> _maxwgt >> _modemap;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SemiLeptonicBaryonDecayer,DecayIntegrator>
describeHerwigSemiLeptonicBaryonDecayer("Herwig::SemiLeptonicBaryonDecayer", "HwBaryonDecay.so");

void SemiLeptonicBaryonDecayer::Init() {

  static ClassDocumentation<SemiLeptonicBaryonDecayer> documentation
    ("The SemiLeptonicBaryonDecayer class is designed for"
     " the semi-leptonic decay of the baryons.");

  static Reference<SemiLeptonicBaryonDecayer,LeptonNeutrinoCurrent> interfaceCurrent
    ("Current",
     "The current for the leptons produced in the decay.",
     &SemiLeptonicBaryonDecayer::_current, true, true, true, false, false);


  static Reference<SemiLeptonicBaryonDecayer,BaryonFormFactor> interfaceFormFactor
    ("FormFactor",
     "The form factor",
     &SemiLeptonicBaryonDecayer::_form, true, true, true, false, false);

  static ParVector<SemiLeptonicBaryonDecayer,double> interfaceMaximumWeight
    ("MaximumWeight",
     "The maximum weights for the decays",
     &SemiLeptonicBaryonDecayer::_maxwgt,
     0, 0, 0, 0, 10000, false, false, true);

}

void SemiLeptonicBaryonDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // for the decaying particle
  if(part.id()>0) {
    SpinorWaveFunction::
      constructSpinInfo(_inHalf,const_ptr_cast<tPPtr>(&part),incoming,true);
  }
  else {
    SpinorBarWaveFunction::
      constructSpinInfo(_inHalfBar,const_ptr_cast<tPPtr>(&part),incoming,true);
  }
  // decay product
  // spin 1/2
  if(decay[0]->dataPtr()->iSpin()==PDT::Spin1Half) {
    if(part.id()>0) {
      SpinorBarWaveFunction::constructSpinInfo(_inHalfBar,decay[0],outgoing,true);
    }
    else {
      SpinorWaveFunction::constructSpinInfo(_inHalf,decay[0],outgoing,true);
    }
  }
  // spin 3/2
  else if(decay[0]->dataPtr()->iSpin()==PDT::Spin3Half) {
    if(part.id()>0) {
      RSSpinorBarWaveFunction::constructSpinInfo(_inThreeHalfBar,
  						 decay[0],outgoing,true);
      
    }
    else {
      RSSpinorWaveFunction::constructSpinInfo(_inThreeHalf,
  					      decay[0],outgoing,true);
    }
  }
  else
    assert(false);
  // and the stuff from the current
  _current->constructSpinInfo(ParticleVector(decay.begin()+1,decay.end()));
}

// combine the currents and form-factors to give the matrix element
double SemiLeptonicBaryonDecayer::me2(const int , const Particle & part,
				      const tPDVector & outgoing,
				      const vector<Lorentz5Momentum> & momenta,
				      MEOption meopt) const {
  assert(part.dataPtr()->iSpin()==PDT::Spin1Half);
  double me(0.);
  if(outgoing[0]->iSpin()==PDT::Spin1Half)
    me =      halfHalf(part,outgoing,momenta,meopt);
  else if(outgoing[0]->iSpin()==PDT::Spin3Half)
    me = halfThreeHalf(part,outgoing,momenta,meopt);
  else
    assert(false);
  return me;
}

// matrix element for a 1/2 -> 1/2 semi-leptonic decay
double SemiLeptonicBaryonDecayer::halfHalf(const Particle & part,
					   const tPDVector & outgoing,
					   const vector<Lorentz5Momentum> & momenta,
					   MEOption meopt) const {
  // spinors etc for the decaying particle
  if(meopt==Initialize) {
    // spinors and rho
    if(part.id()>0)
      SpinorWaveFunction   ::calculateWaveFunctions(_inHalf,_rho,
  						    const_ptr_cast<tPPtr>(&part),
  						    incoming);
    else
      SpinorBarWaveFunction::calculateWaveFunctions(_inHalfBar,_rho,
  						    const_ptr_cast<tPPtr>(&part),
  						    incoming);
    // work out the mapping for the lepton vector
    _constants.resize(outgoing.size()+1);
    _ispin.resize(outgoing.size());
    int itemp(1);
    _ibar=0;
    for(int ix=int(outgoing.size()-1);ix>=0;--ix) {
      _ispin[ix]=outgoing[ix]->iSpin();
      if(abs(outgoing[ix]->id())<=16) {
  	itemp*=_ispin[ix];
  	_constants[ix]=itemp;
      }
      else _ibar=ix;
    }
    _constants[outgoing.size()]=1;
    _constants[_ibar]=_constants[_ibar+1];
  }
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half,_ispin)));
  // spinors for the decay product
  if(part.id()>0) {
    _inHalfBar.resize(2);
    for(unsigned int ihel=0;ihel<2;++ihel)
      _inHalfBar[ihel] = HelicityFunctions::dimensionedSpinorBar(-momenta[0],ihel,Helicity::outgoing);
  }
  else {
    _inHalf.resize(2);
    for(unsigned int ihel=0;ihel<2;++ihel)
      _inHalf[ihel] = HelicityFunctions::dimensionedSpinor   (-momenta[0],ihel,Helicity::outgoing);
  }
  // get the information on the form-factor
  int spinin(0),spinout(0),spect1,spect2,inquark,outquark;
  int id0(part.id()),id1(outgoing[0]->id());
  bool cc; int iloc(_form->formFactorNumber(id0,id1,cc));
  _form->formFactorInfo(iloc,spinin,spinout,spect1,spect2,inquark,outquark);
  // work out the value of q and calculate the form factors
  Lorentz5Momentum q(part.momentum()-momenta[0]);
  q.rescaleMass();
  Energy m0(part.mass()),m1(momenta[0].mass());
  Energy2 q2(q.mass2());
  Lorentz5Momentum sum(part.momentum()+momenta[0]);
  // calculate the form factors
  Complex f1v,f2v,f3v,f1a,f2a,f3a;
  _form->SpinHalfSpinHalfFormFactor(q2,iloc,id0,id1,m0,m1,
  				    f1v,f2v,f3v,f1a,f2a,f3a);
  // calculate the hadronic current for the decay
  vector<LorentzPolarizationVectorE> hadron(4);
  Complex left  =f1v-f1a-f2v-double((m0-m1)/(m0+m1))*f2a;
  Complex right =f1v+f1a-f2v+double((m0-m1)/(m0+m1))*f2a;
  LorentzPolarizationVectorE vtemp;
  for(unsigned int ix=0;ix<2;++ix) {
    for(unsigned int iy=0;iy<2;++iy) {
      vtemp = _inHalf[ix].generalCurrent(_inHalfBar[iy],left,right);
      complex<Energy> vspin = _inHalf[ix].scalar(_inHalfBar[iy]);
      complex<Energy> aspin = _inHalf[ix].pseudoScalar(_inHalfBar[iy]);
      // the momentum like pieces
      if(part.id()>0) {
  	vtemp+= (f2v*vspin+f2a*aspin)/(m0+m1)*sum;
  	vtemp+= (f3v*vspin+f3a*aspin)/(m0+m1)*q;
      }
      else {
  	vtemp-= (f2v*vspin-f2a*aspin)/(m0+m1)*sum;
  	vtemp+= (f3v*vspin-f3a*aspin)/(m0+m1)*q;
      }
      if(part.id()>0) hadron[2*ix+iy]=vtemp;
      else            hadron[2*iy+ix]=vtemp;
    }
  }
  // construct the lepton current
  Energy scale;
  int mode((abs(outgoing[1]->id())-11)/2);
  vector<LorentzPolarizationVectorE> 
    lepton(_current->current(tcPDPtr(),IsoSpin::IUnknown,IsoSpin::I3Unknown,
			     mode,-1,scale,tPDVector(outgoing.begin()+1,outgoing.end()),
			     vector<Lorentz5Momentum>(momenta.begin()+1,momenta.end()),
			     meopt));
  // matrix element
  vector<unsigned int> ihel(outgoing.size()+1);
  unsigned int mhel,ix,lhel;
  for(mhel=0;mhel<hadron.size();++mhel) {
    ihel[0      ]=mhel/spinout;
    ihel[_ibar+1]=mhel%spinout;
    for(lhel=0;lhel<lepton.size();++lhel) {
      // map the index for the leptons to a helicity state
      for(ix=outgoing.size();ix>0;--ix) {
  	if(ix-1!=_ibar) ihel[ix]=(lhel%_constants[ix-1])/_constants[ix];
      }
      (*ME())(ihel)= lepton[lhel].dot(hadron[mhel])*SM().fermiConstant();
    }
  }
  // ckm factor
  double ckm(1.);
  if(inquark<=6) {
    if(inquark%2==0) ckm = SM().CKM(inquark/2-1,(abs(outquark)-1)/2);
    else             ckm = SM().CKM(abs(outquark)/2-1,(inquark-1)/2);
  }
  // return the answer
  return 0.5*(ME()->contract(_rho)).real()*ckm; 
}



// matrix element for a 1/2 -> 3/2 semi-leptonic decay
double SemiLeptonicBaryonDecayer::halfThreeHalf(const Particle & part,
						const tPDVector & outgoing,
						const vector<Lorentz5Momentum> & momenta,
						MEOption meopt) const {
  // spinors etc for the decaying particle
  if(meopt==Initialize) {
    // spinors and rho
    if(part.id()>0)
      SpinorWaveFunction   ::calculateWaveFunctions(_inHalf,_rho,
  						    const_ptr_cast<tPPtr>(&part),
  						    incoming);
    else
      SpinorBarWaveFunction::calculateWaveFunctions(_inHalfBar,_rho,
  						    const_ptr_cast<tPPtr>(&part),
  						    incoming);
    // work out the mapping for the lepton vector
    _constants.resize(outgoing.size()+1);
    _ispin.resize(outgoing.size());
    int itemp(1);
    _ibar=0;
    for(int ix=int(outgoing.size()-1);ix>=0;--ix) {
      _ispin[ix]=outgoing[ix]->iSpin();
      if(abs(outgoing[ix]->id())<=16) {
  	itemp*=_ispin[ix];
  	_constants[ix]=itemp;
      }
      else _ibar=ix;
    }
    _constants[outgoing.size()]=1;
    _constants[_ibar]=_constants[_ibar+1];
  }
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half,_ispin)));
  // spinors for the decay product
  LorentzPolarizationVector in=UnitRemoval::InvE*part.momentum();
  if(part.id()>0) {
    RSSpinorBarWaveFunction swave(momenta[0],outgoing[0],Helicity::outgoing);
    _inThreeHalfBar.resize(4);
    _inHalfBar.resize(4);
    for(unsigned int ihel=0;ihel<4;++ihel) {
      swave.reset(ihel);
      _inThreeHalfBar[ihel]=swave.dimensionedWf();
      _inHalfBar[ihel] = _inThreeHalfBar[ihel].dot(in);
    }
  }
  else {
    RSSpinorWaveFunction swave(momenta[0],outgoing[0],Helicity::outgoing);
    _inThreeHalf.resize(4);
    _inHalf.resize(4);
    for(unsigned int ihel=0;ihel<4;++ihel) {
      swave.reset(ihel);
      _inThreeHalf[ihel]=swave.dimensionedWf();
      _inHalf[ihel] = _inThreeHalf[ihel].dot(in);
    }
  }
  // get the information on the form-factor
  int spinin(0),spinout(0),inquark,outquark,spect1,spect2;
  int id0(part.id()),id1(outgoing[0]->id());
  bool cc; int iloc(_form->formFactorNumber(id0,id1,cc));
  _form->formFactorInfo(iloc,spinin,spinout,spect1,spect2,inquark,outquark);
  // work out the value of q and calculate the form factors
  Lorentz5Momentum q(part.momentum()-momenta[0]);
  q.rescaleMass();
  Energy m0(part.mass()),m1(momenta[0].mass());
  Energy2 q2(q.mass2());
  Lorentz5Momentum sum(part.momentum()+momenta[0]);
  // calculate the form factors
  Complex f1v,f2v,f3v,f4v,f1a,f2a,f3a,f4a;
  _form->SpinHalfSpinThreeHalfFormFactor(q2,iloc,id0,id1,m0,m1,
  					 f1v,f2v,f3v,f4v,f1a,f2a,f3a,f4a);
  LorentzPolarizationVector vtemp;
  complex<InvEnergy2> lS1,lS2,rS1,rS2;
  complex<InvEnergy> lV,rV;
  Complex left,right;
  InvEnergy ms(1./(m0+m1));
  InvEnergy2 ms2(ms*ms);
  if(part.id()>0) {
    left  = f1a-f1v;
    right = f1a+f1v; 
    lS1   = ms2*(f3a-f4a-f3v+f4v);
    rS1   = ms2*(f3a-f4a+f3v-f4v);
    lS2   = ms2*(f4a-f4v);
    rS2   = ms2*(f4a+f4v);
    lV    = ms*(f2a-f2v);
    rV    = ms*(f2a+f2v);
  }
  else {
    left  = conj(f1a+f1v);
    right = conj(f1a-f1v); 
    lS1   = ms2*conj(f3a-f4a+f3v-f4v);
    rS1   = ms2*conj(f3a-f4a-f3v+f4v);
    lS2   = ms2*conj(f4a-f4v);
    rS2   = ms2*conj(f4a+f4v);
    lV    = ms *conj(f2a-f2v);
    rV    = ms *conj(f2a+f2v);
  }
  // calculate the hadronic current for the decay
  LorentzPolarizationVectorE hadron[4][2];
  // construct the vectors for the decay
  Complex scalar1,scalar2;
  complex<Energy> lfact,rfact;
  LorentzPolarizationVectorE tvec;
  LorentzPolarizationVector svec;
  for(unsigned int iya=0;iya<4;++iya) {
    for(unsigned int ixa=0;ixa<2;++ixa) {
      unsigned int ix=iya,iy=ixa;
      if(outgoing[0]->id()<0) swap(ix,iy);
      // scalar like terms
      lfact = _inHalf[iy].leftScalar(_inHalfBar[ix]);
      rfact = _inHalf[iy].rightScalar(_inHalfBar[ix]);
      scalar1 = (lS1*lfact+rS1*rfact)*UnitRemoval::E;
      scalar2 = (lS2*lfact+rS2*rfact)*UnitRemoval::E;
      svec = _inHalf[iy].generalCurrent(_inHalfBar[ix],lV/ms,rV/ms)*ms;
      if(part.id()>0) tvec=_inThreeHalfBar[ix].generalCurrent(_inHalf[iy],left,right);
      else              tvec=_inThreeHalf[iy].generalCurrent(_inHalfBar[ix],left,right);
      hadron[iya][ixa] = tvec+svec*UnitRemoval::E+scalar1*momenta[0]
  	+scalar2*part.momentum();
    }
  }
  // construct the lepton current
  Energy scale;
  int mode((abs(outgoing[1]->id())-11)/12);
  vector<LorentzPolarizationVectorE> 
    lepton(_current->current(tcPDPtr(),IsoSpin::IUnknown,IsoSpin::I3Unknown,
			     mode,-1,scale,tPDVector(outgoing.begin()+1,outgoing.end()),
			     vector<Lorentz5Momentum>(momenta.begin()+1,momenta.end()),meopt));
  vector<unsigned int> ihel(outgoing.size()+1);
  for(unsigned int iya=0;iya<4;++iya) {
    ihel[1]=iya;
    for(unsigned int ixa=0;ixa<2;++ixa) {
      ihel[0]=ixa;
      for(unsigned int lhel=0;lhel<lepton.size();++lhel) {
  	ihel[2] = lhel/2;
  	ihel[3] = lhel%2;
  	(*ME())(ihel) = lepton[lhel].dot(hadron[iya][ixa])*SM().fermiConstant();
      }
    }  
  }
  // ckm factor
  double ckm(1.);
  if(inquark<=6) {
    if(inquark%2==0){ckm = SM().CKM(inquark/2-1,(abs(outquark)-1)/2);}
    else{ckm = SM().CKM(abs(outquark)/2-1,(inquark-1)/2);}
  }
  // return the answer
  return 0.5*(ME()->contract(_rho)).real()*ckm;
}

// output the setup information for the particle database
void SemiLeptonicBaryonDecayer::dataBaseOutput(ofstream & output,bool header) const {
  if(header) output << "update decayers set parameters=\"";
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<_maxwgt.size();++ix) {
    output << "insert " << name() << ":MaximumWeight " << ix << " " 
	   << _maxwgt[ix] << " \n";
  }
  _current->dataBaseOutput(output,false,true);
  output << "newdef " << name() << ":Current " << _current->name() << " \n";
  _form->dataBaseOutput(output,false,true);
  output << "newdef " << name() << ":FormFactor " << _form->name() << " \n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
