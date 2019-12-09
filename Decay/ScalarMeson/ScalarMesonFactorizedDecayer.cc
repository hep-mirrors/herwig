// -*- C++ -*-
//
// ScalarMesonFactorizedDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ScalarMesonFactorizedDecayer class.
//

#include "ScalarMesonFactorizedDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/epsilon.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;


ScalarMesonFactorizedDecayer::ScalarMesonFactorizedDecayer() 
// default values of the couplings (taken from ZPC34, 103)
  : _a1b(1.10), _a2b(-0.24), _a1c(1.30), _a2c(-0.55) { 
  // intermediates
  generateIntermediates(true);
}

void ScalarMesonFactorizedDecayer::rebind(const TranslationMap & trans)
  {
  _ckm = trans.translate(_ckm);
  DecayIntegrator::rebind(trans);
}

IVector ScalarMesonFactorizedDecayer::getReferences() {
  IVector ret = DecayIntegrator::getReferences();
  ret.push_back(_ckm);
  return ret;
}

void ScalarMesonFactorizedDecayer::doinit() {
  DecayIntegrator::doinit();
  // get the ckm object
  _ckm=dynamic_ptr_cast<Ptr<StandardCKM>::pointer>(SM().CKM());
  if(!_ckm) throw InitException() << "ScalarMesonFactorizedDecayer::doinit() "
				  << "the CKM object must be the Herwig one"
				  << Exception::runerror;
  unsigned int ix,iy,iz,icurr,iform;
  // get the CKM matrix (unsquared for interference)
  Complex ckmmat[3][3];
  vector< vector<Complex > > CKM(_ckm->getUnsquaredMatrix(SM().families()));
  for(ix=0;ix<3;++ix){for(iy=0;iy<3;++iy){ckmmat[ix][iy]=CKM[ix][iy];}}
  int id0,id1,Wcharge,iq,ia,jspin,spect,inq,outq;
  // make sure the currents and form factors got initialised
  for(ix=0;ix<_current.size();++ix) _current[ix]->init();
  for(ix=0;ix<_form.size();++ix)    _form[ix]->init();
  // find all the possible modes
  vector<unsigned int> tformmap[2],tcurrmap[2];
  vector<int> inquark,outquark,currq,curra;
  vector<tPDVector> particles;
  tPDVector extpart,ptemp;
  Energy min,minb;
  // loop over the modes in the form factors and currents
  for(iform=0;iform<_form.size();++iform) {
    for(ix=0;ix<_form[iform]->numberOfFactors();++ix) {
      // particles from the form-factor
      extpart.resize(2);
      _form[iform]->particleID(ix,id0,id1);
      _form[iform]->formFactorInfo(ix,jspin,spect,inq,outq);
      // particles from the form factor
      extpart[0]=getParticleData(id0);
      extpart[1]=getParticleData(id1);
      // charge of the decay products
      Wcharge =extpart[0]->iCharge()-extpart[1]->iCharge();
      // max mass for the particles in the current
      min = extpart[0]->massMax()-extpart[1]->massMin();
      for(icurr=0;icurr<_current.size();++icurr) {
	for(iy=0;iy<_current[icurr]->numberOfModes();++iy) {
	  extpart.resize(2);
	  // get the particles from the current
	  _current[icurr]->decayModeInfo(iy,iq,ia);
	  ptemp=_current[icurr]->particles(Wcharge,iy,iq,ia);
	  minb=ZERO;
	  for(iz=0;iz<ptemp.size();++iz) {
	    extpart.push_back(ptemp[iz]);
	    minb+=ptemp[iz]->massMin();
	  }
	  // add this mode to the list
	  if(extpart.size()>2&&minb<min&&
	     (Wcharge!=0||(Wcharge==0&&((inq>0&&inq%2!=iq%2)||
					(inq<0&&abs(inq)%2!=abs(ia)%2))))) {
	    tformmap[0].push_back(iform);tformmap[1].push_back(ix);
	    tcurrmap[0].push_back(icurr);tcurrmap[1].push_back(iy);
	    particles.push_back(extpart);
	    inquark.push_back(inq);outquark.push_back(outq);
	    currq.push_back(iq);curra.push_back(ia);
	  }
	  // if the meson in the current is neutral try the CC mode
	  if(Wcharge==0&&iq!=-ia&&((inq>0&&inq%2!=iq%2)||
				   (inq<0&&abs(inq)%2!=abs(ia)%2))) {
	    extpart.resize(2);
	    // get the particles from the current
	    ptemp=_current[icurr]->particles(Wcharge,iy,-ia,-iq);
	    minb=ZERO;
	    for(iz=0;iz<ptemp.size();++iz) {
	      extpart.push_back(ptemp[iz]);
	      minb+=ptemp[iz]->massMin();
	    }
	    if(extpart.size()>2&&minb<min) {
	      tformmap[0].push_back(iform);tformmap[1].push_back(ix);
	      tcurrmap[0].push_back(icurr);tcurrmap[1].push_back(iy);
	      particles.push_back(extpart);
	      inquark.push_back(inq);outquark.push_back(outq);
	      currq.push_back(-ia);curra.push_back(-iq);
	    }
	  }
	}
      }
    }
  }
  // loop over the modes and find the dupliciates
  vector<bool> modecc; vector<unsigned int> modeloc,tformpart,ttform[2],ttcurr[2];
  vector<Complex> tCKM; Complex ckm;
  DecayPhaseSpaceModePtr mode;
  DecayPhaseSpaceChannelPtr channel;
  bool done;
  int id,idbar;
  double maxweight;
  vector<double> channelwgts;
  unsigned int isize;double ort(sqrt(0.5));
  vector<double>::iterator start,end;
  for(ix=0;ix<particles.size();++ix) {
    while (true) {
      if(particles[ix].empty()) break;
      findModes(ix,particles,modeloc,modecc);
      // if more than three particles only allow one diagram
      if ( particles[ix].size()>3 && !modeloc.empty() ) break;
      // create the mode and set the particles as for the first instance
      mode=new_ptr(DecayPhaseSpaceMode(particles[ix],this));
      channel = new_ptr(DecayPhaseSpaceChannel(mode));
      channel->addIntermediate(particles[ix][0],0,0.0,1,-1);
      min = particles[ix][0]->massMax()-particles[ix][1]->massMin();
      Wcharge = particles[ix][0]->iCharge()-particles[ix][1]->iCharge();
      done=_current[tcurrmap[0][ix]]->createMode(Wcharge,tcurrmap[1][ix],
						 mode,2,1,channel,min);
      if(!done) throw InitException() << "Failed to construct mode in "
				      << "ScalarMesonFactorizedDecayer::doinit()." 
				      << Exception::abortnow;
      // set the parameters for the additional modes
      ttform[0].clear();ttform[1].clear();
      ttcurr[0].clear();ttcurr[1].clear();
      ttform[0].push_back(tformmap[0][ix]);ttform[1].push_back(tformmap[1][ix]);
      ttcurr[0].push_back(tcurrmap[0][ix]);ttcurr[1].push_back(tcurrmap[1][ix]);
      tformpart.clear();tformpart.push_back(0);
      id=particles[ix][1]->id();
      idbar = particles[ix][1]->CC() ? particles[ix][1]->CC()->id() : id;
      for(iy=0;iy<modeloc.size();++iy) {
	ttform[0].push_back(tformmap[0][modeloc[iy]]);
	ttform[1].push_back(tformmap[1][modeloc[iy]]);
	ttcurr[0].push_back(tcurrmap[0][modeloc[iy]]);
	ttcurr[1].push_back(tcurrmap[1][modeloc[iy]]);
	iz=1;
	do {
	  if(( modecc[iy]&&particles[modeloc[iy]][iz]->id()==idbar)||
	     (!modecc[iy]&&particles[modeloc[iy]][iz]->id()==id))
	    tformpart.push_back(iz-1);
	  ++iz;
	}
	while(tformpart.size()!=iy+2&&iz<3);
      }
      // calculate ckm factors
      tCKM.clear();
      for(iy=0;iy<ttcurr[0].size();++iy) {
	// get the quarks involved in the process
	if(iy==0) {
	  iq=currq[ix];ia=curra[ix];
	  inq=inquark[ix];outq=outquark[ix];
	}
	else {
	  if(!modecc[iy-1]) {
	    iq=currq[modeloc[iy-1]];ia=curra[modeloc[iy-1]];
	    inq=inquark[modeloc[iy-1]];outq=outquark[modeloc[iy-1]];
	  }
	  else {
	    ia=-currq[modeloc[iy-1]];iq=-curra[modeloc[iy-1]];
	    inq=-inquark[modeloc[iy-1]];outq=-outquark[modeloc[iy-1]];
	  }
	}
	_form[ttform[0][iy]]->particleID(ttform[1][iy],id0,id1);
	Wcharge = getParticleData(id0)->iCharge()-getParticleData(id1)->iCharge();
	ckm=1.;
	if(Wcharge!=0) {
	  ckm=1.;
	  if(abs(iq)%2==0)  ckm *= conj(ckmmat[abs(iq)/2-1][(abs(ia)-1)/2]);
	  else              ckm *= conj(ckmmat[abs(ia)/2-1][(abs(iq)-1)/2]);
	  if(abs(inq)%2==0) ckm *= ckmmat[abs(inq)/2-1][(abs(outq)-1)/2];
	  else              ckm *= ckmmat[abs(outq)/2-1][(abs(inq)-1)/2];
	  if(abs(inq)==5)   ckm*=_a1b;
	  else              ckm*=_a1c;
	}
	else {
	  ckm=1.;
	  if(inq>0) {
	    if(abs(inq)%2==0)  ckm *= ckmmat[abs(inq)/2-1][(abs(iq)-1)/2];
	    else               ckm *= ckmmat[abs(iq)/2-1][(abs(inq)-1)/2];
	    if(abs(outq)%2==0) ckm *= conj(ckmmat[abs(outq)/2-1][(abs(ia)-1)/2]);
	    else               ckm *= conj(ckmmat[abs(ia)/2-1][(abs(outq)-1)/2]);
	  }
	  else {
	    if(abs(inq)%2==0)  ckm *= ckmmat[abs(inq)/2-1][(abs(ia)-1)/2];
	    else               ckm *= ckmmat[abs(ia)/2-1][(abs(inq)-1)/2];
	    if(abs(outq)%2==0) ckm *= conj(ckmmat[abs(outq)/2-1][(abs(iq)-1)/2]);
	    else               ckm *= conj(ckmmat[abs(iq)/2-1][(abs(outq)-1)/2]);
	  }
	  if(abs(inq)==5) ckm*=_a2b;
	  else            ckm*=_a2c;
	}
	if((abs(inq)%2==0&&inq<0)||(abs(inq)%2!=0&&inq>0)){ckm=conj(ckm);}
	tCKM.push_back(ckm);
      }
      // special if the particles are idential add additional modes and 
      // identical particle factors
      if(particles[ix][1]->id()==particles[ix][2]->id()&&particles[ix].size()==3) {
	isize=ttcurr[0].size();
	for(iy=0;iy<isize;++iy) {
	  ttcurr[0].push_back(ttcurr[0][iy]);ttcurr[1].push_back(ttcurr[1][iy]);
	  ttform[0].push_back(ttform[0][iy]);ttform[1].push_back(ttform[1][iy]);
	  if(tformpart[iy]==0){tformpart.push_back(1);}
	  else{tformpart.push_back(0);}
	  tCKM[iy]*=ort;tCKM.push_back(tCKM[iy]);
	}
      }
      // add the parameters for the mode to the list
      _currentmapA.push_back(ttcurr[0]);_currentmapB.push_back(ttcurr[1]);
      _formmapA.push_back(ttform[0]);_formmapB.push_back(ttform[1]);
      _formpart.push_back(tformpart);
      _CKMfact.push_back(tCKM);
      // add the mode to the list
      if(_wgtmax.size()>numberModes()){maxweight=_wgtmax[numberModes()];}
      else{maxweight=0.;}
      // the weights for the channel
      if(_wgtloc.size()>numberModes()&&
	 _wgtloc[numberModes()]+mode->numberChannels()<=_weights.size()) {
	start=_weights.begin()+_wgtloc[numberModes()];
	end  = start+mode->numberChannels();
	channelwgts=vector<double>(start,end);
      }
      else {
	channelwgts.resize(mode->numberChannels(),1./(mode->numberChannels()));
      }
      // don't need channels for two body decays
      if(particles[ix].size()==3) {
	channelwgts.clear();
	mode=new_ptr(DecayPhaseSpaceMode(particles[ix],this));
      }
      addMode(mode,maxweight,channelwgts);
      // resize the duplicate modes to remove them
      for(iy=0;iy<modeloc.size();++iy) particles[modeloc[iy]] = tPDVector();
      break;
    }
  }
}

void ScalarMesonFactorizedDecayer::doinitrun() {
  unsigned int ix,iy;
  for(ix=0;ix<_current.size();++ix) _current[ix]->initrun();
  for(ix=0;ix<_form.size();++ix)    _form[ix]->initrun();
  DecayIntegrator::doinitrun();
  if(initialize()) {
    _weights.clear();
    _wgtloc.clear();
    _wgtmax.clear();
    for(ix=0;ix<numberModes();++ix) {
      _wgtmax.push_back(mode(ix)->maxWeight());
      _wgtloc.push_back(_weights.size());
      for(iy=0;iy<mode(ix)->numberChannels();++iy) {
	_weights.push_back(mode(ix)->channelWeight(iy));
      }
    }
  }
}

bool ScalarMesonFactorizedDecayer::accept(tcPDPtr parent,
					  const tPDVector & children) const {
  // N.B. this is a necessary but not sufficient test
  bool allowed(false),dummy;
  // find the ids of the particles
  tPDVector::const_iterator pit  = children.begin();
  tPDVector::const_iterator pend = children.end();
  vector<int> ids,idcurr;
  int id(parent->id());
  for( ; pit!=pend;++pit) ids.push_back((**pit).id());
  // loop over the possible particles in the formfactor
  unsigned int ipart(0),iform,icurr,ix;
  do {
    idcurr.clear();
    for(ix=0;ix<ids.size();++ix){if(ix!=ipart){idcurr.push_back(ids[ix]);}}
    iform=0;
    do {
      // check if possible from the form factor
      if(_form[iform]->formFactorNumber(id,ids[ipart],dummy)>=0) {
	// check if possible from the current
	icurr=0;
	do {
	  allowed=_current[icurr]->accept(idcurr);
	  ++icurr;
	}
	while(!allowed&&icurr<_current.size());
      }
      ++iform;
    }
    while(!allowed&&iform<_form.size());
    ++ipart;
  }
  while(!allowed&&ipart<ids.size());
  return allowed;
}

int ScalarMesonFactorizedDecayer::modeNumber(bool & cc,tcPDPtr parent,
					     const tPDVector & children) const {
  int imode(-1);
  // id's of the particles and CC
  // of the parent
  int id0(parent->id()),id0bar(id0);
  if(parent->CC()){id0bar=parent->CC()->id();}
  // of the products
  vector<int> ids,idbars;
  tPDVector::const_iterator pit  = children.begin();
  tPDVector::const_iterator pend = children.end();
  for( ;pit!=pend;++pit) {
    ids.push_back((**pit).id());
    if((**pit).CC()) idbars.push_back((**pit).CC()->id());
    else             idbars.push_back(ids.back());
  }
  // loop over the modes
  vector<bool> done(ids.size(),false);
  unsigned int nfound,ix,iy,iz;
  int idtemp;
  bool found;
  cc=false;
  ix=0;
  do {
    // particle mode
    if(id0==mode(ix)->externalParticles(0)->id()&&
       ids.size()+1==mode(ix)->numberofParticles()) {
      nfound=0;
      for(iy=0;iy<ids.size();++iy){done[iy]=false;}
      for(iy=0;iy<ids.size();++iy) {
	idtemp=mode(ix)->externalParticles(iy+1)->id();
	iz=0;found=false;
	do{if(idtemp==ids[iz]&&!done[iz]){done[iz]=true;found=true;}++iz;}
	while(iz<ids.size()&&!found);
	if(found){++nfound;}
      }
      if(nfound==ids.size()){cc=false;imode=ix;}
    }
    // CC mode
    if(id0bar==mode(ix)->externalParticles(0)->id()&&
       ids.size()+1==mode(ix)->numberofParticles()) {
      nfound=0;
      for(iy=0;iy<idbars.size();++iy) done[iy]=false;
      for(iy=0;iy<idbars.size();++iy) {
	idtemp=mode(ix)->externalParticles(iy+1)->id();
	iz=0;found=false;
	do {
	  if(idtemp==idbars[iz]&&!done[iz]) {
	    done[iz]=true;
	    found=true;
	  }
	  ++iz;
	}
	while(iz<idbars.size()&&!found);
	if(found) ++nfound;
      }
      if(nfound==idbars.size()) {
	cc=true;
	imode=ix;
      }
    }
    ++ix;
  }
  while(imode<0&&ix<numberModes());
  if(imode<0) {
    string mode = parent->PDGName() + "->";
    for(unsigned int ix=0;ix<children.size();++ix) 
      mode += children[ix]->PDGName() +",";
    throw DecayIntegratorError() << "Unable to find the mode " << mode << " in " 
				 << name() 
				 << " ScalarMesonFactorizedDecayer::decay()" 
				 << Exception::abortnow;
  }
  return imode;
}


void ScalarMesonFactorizedDecayer::persistentOutput(PersistentOStream & os) const {
  os << _current << _form << _ckm 
     << _a1b << _a2b << _a1c << _a2c 
     << _currentmapA << _currentmapB 
     << _formmapA << _formmapB << _formpart << _wgtloc 
     << _wgtmax << _weights << _CKMfact ;
}

void ScalarMesonFactorizedDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _current >> _form >> _ckm 
     >> _a1b >> _a2b >> _a1c >> _a2c 
     >> _currentmapA >> _currentmapB 
     >> _formmapA >> _formmapB >> _formpart >> _wgtloc
     >> _wgtmax >> _weights >> _CKMfact;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ScalarMesonFactorizedDecayer,DecayIntegrator>
describeHerwigScalarMesonFactorizedDecayer("Herwig::ScalarMesonFactorizedDecayer", "HwSMDecay.so");

void ScalarMesonFactorizedDecayer::Init() {

  static ClassDocumentation<ScalarMesonFactorizedDecayer> documentation
    ("The ScalarMesonFactorizedDecayer class is designed for the weak decay of"
     " scalar mesons using the factorization approximation.");

  static RefVector<ScalarMesonFactorizedDecayer,WeakDecayCurrent> interfaceCurrents
    ("Currents",
     "A vector of references to the currents",
     &ScalarMesonFactorizedDecayer::_current, -1, false, false, true, false, false);

  static RefVector<ScalarMesonFactorizedDecayer,ScalarFormFactor> interfaceFormFactors
    ("FormFactors",
     "A vector of references to the form-factors",
     &ScalarMesonFactorizedDecayer::_form, -1, false, false, true, false, false);

  static Parameter<ScalarMesonFactorizedDecayer,double> interfacea1Bottom
    ("a1Bottom",
     "The factorization paramter a_1 for decays of bottom baryons",
     &ScalarMesonFactorizedDecayer::_a1b, 1.1, -10.0, 10.0,
     false, false, true);

  static Parameter<ScalarMesonFactorizedDecayer,double> interfacea2Bottom
    ("a2Bottom",
     "The factorization paramter a_2 for decays of bottom baryons",
     &ScalarMesonFactorizedDecayer::_a2b, -0.24, -10.0, 10.0,
     false, false, true);

  static Parameter<ScalarMesonFactorizedDecayer,double> interfacea1Charm
    ("a1Charm",
     "The factorization paramter a_1 for decays of charm baryons",
     &ScalarMesonFactorizedDecayer::_a1c, 1.3, -10.0, 10.0,
     false, false, true);

  static Parameter<ScalarMesonFactorizedDecayer,double> interfacea2Charm
    ("a2Charm",
     "The factorization paramter a_2 for decays of charm baryons",
     &ScalarMesonFactorizedDecayer::_a2c, -0.55, -10.0, 10.0,
     false, false, true);

  static ParVector<ScalarMesonFactorizedDecayer,int> interfaceWeightLocation
    ("WeightLocation",
     "The locations of the weights for a given channel in the vector",
     &ScalarMesonFactorizedDecayer::_wgtloc,
     0, 0, 0, 0, 10000, false, false, true);

  static ParVector<ScalarMesonFactorizedDecayer,double> interfaceWeightMax
    ("MaximumWeight",
     "The maximum weight for a given channel.",
     &ScalarMesonFactorizedDecayer::_wgtmax,
     0, 0, 0, 0., 100., false, false, true);

  static ParVector<ScalarMesonFactorizedDecayer,double> interfaceWeights
    ("Weights",
     "The weights for the integration.",
     &ScalarMesonFactorizedDecayer::_weights,
     0, 0, 0, 0., 1., false, false, true);
}
 
double ScalarMesonFactorizedDecayer::me2(const int ichan,
					 const Particle & part,
					 const ParticleVector & decay,
					 MEOption meopt) const {
  if(!ME()) {
    // create the matrix element
    vector<PDT::Spin> spin;
    for(unsigned int ix=0;ix<decay.size();++ix)
      spin.push_back(decay[ix]->dataPtr()->iSpin());
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,spin)));
  }
  // initialisation
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&part),incoming);
    _vectors.resize(decay.size());
    _tensors.resize(decay.size());
  }
  if(meopt==Terminate) {
    // set up the spin information for the decay products
    ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					  incoming,true);
    // get the wavefunctions of the decay products
    for(unsigned int ix=0;ix<decay.size();++ix) {
      switch(decay[ix]->dataPtr()->iSpin()) {
      case 1:
	ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
	break;
      case 3:
	VectorWaveFunction::constructSpinInfo(_vectors[ix],decay[ix],outgoing,
					      true,false);
	break;
      case 5:
	TensorWaveFunction::constructSpinInfo(_tensors[ix],decay[ix],outgoing,
					      true,false);
	break;
      default:
	assert(false);
      }
    }
    return 0.;
  }
  // get the wavefunctions of the decay products
  for(unsigned int ix=0;ix<decay.size();++ix) {
    switch(decay[ix]->dataPtr()->iSpin()) {
    case 1:
      break;
    case 3:
      VectorWaveFunction::
	calculateWaveFunctions(_vectors[ix],decay[ix],outgoing,false);
      break;
    case 5:
      TensorWaveFunction::
	calculateWaveFunctions(_tensors[ix],decay[ix],outgoing,false);
      break;
    default:
      assert(false);
    }
  }
  ME()->zero();
  // find the mode
  unsigned int mode(imode()),chel,fhel;
  int id0(part.id()),id1;
  Complex ii(0.,1.);
  vector<unsigned int> ihel(decay.size());
  // loop over the different diagrams
  vector<LorentzPolarizationVectorE> form;
  Complex fp,f0,A0,A1,A2,A3,V,k;
  complex<InvEnergy2> h,bp,bm;
  // complex<Energy2> dot;
  Lorentz5Momentum q,sum; 
  Energy2 q2;
  Energy MP(part.mass()),MV,msum,mdiff,scale;
  LorentzPolarizationVectorE dotv;
  double pre;
  ParticleVector cpart;
  for(unsigned int iy=0;iy<_CKMfact[mode].size();++iy) {
    MV=decay[_formpart[mode][iy]]->mass();
    id1=decay[_formpart[mode][iy]]->id();
    int id0t,id1t;
    _form[_formmapA[mode][iy]]->particleID(_formmapB[mode][iy],id0t,id1t);
    bool cc(id0t!=id0);
    // calculate the form-factor part
    form.clear();
    q   = part.momentum()-decay[_formpart[mode][iy]]->momentum();  q.rescaleMass();
    sum = part.momentum()+decay[_formpart[mode][iy]]->momentum();sum.rescaleMass();
    q2=q.mass2();
    if(decay[_formpart[mode][iy]]->dataPtr()->iSpin()==1) {
      _form[_formmapA[mode][iy]]->ScalarScalarFormFactor(q2,_formmapB[mode][iy],
							 id0,id1,MP,MV,f0,fp);
      pre=(MP*MP-MV*MV)/q2;
      form.push_back(fp*sum+pre*(f0-fp)*q);
    }
    else if(decay[_formpart[mode][iy]]->dataPtr()->iSpin()==3) {
      msum=MP+MV;mdiff=MP-MV;
      _form[_formmapA[mode][iy]]->ScalarVectorFormFactor(q2,_formmapB[mode][iy],id0,
							 id1,MP,MV,A0,A1,A2,V);
      if(cc){V=-V;}
      A3 = Complex(0.5/MV*(msum*A1-mdiff*A2));
      // compute the hadron currents
      for(unsigned int ix=0;ix<3;++ix) {
	// dot product
	complex<Energy> dot = _vectors[_formpart[mode][iy]][ix]*part.momentum();
	// current
	form.push_back(-ii*msum*A1*_vectors[_formpart[mode][iy]][ix]
		       +ii*A2/msum*dot*sum
		       +2.*ii*MV/q2*(A3-A0)*dot*q
		       +2.*V/msum*epsilon(_vectors[_formpart[mode][iy]][ix],
					  part.momentum(),
					  decay[_formpart[mode][iy]]->momentum())); 
      }
    }
    else if(decay[_formpart[mode][iy]]->dataPtr()->iSpin()==5) {
      _form[_formmapA[mode][iy]]->ScalarTensorFormFactor(q2,_formmapB[mode][iy],
							 id0,id1,MP,MV,h,k,bp,bm);
      if(cc){h=-h;}
      // compute the hadron currents
      for(unsigned int ix=0;ix<5;++ix) {
	dotv = _tensors[_formpart[mode][iy]][ix]*part.momentum();
	complex<Energy2> dot = dotv*part.momentum();
	form.push_back(ii*h*epsilon(dotv,sum,q)-k*dotv
		       -bp*dot*sum-bm*dot*q);
      }
    }
    // find the particles for the current
    cpart.clear();
    for(unsigned int ix=0;ix<decay.size();++ix)
      {if(ix!=_formpart[mode][iy]){cpart.push_back(decay[ix]);}}
    unsigned int ix=decay.size();
    vector<unsigned int> constants(decay.size()+1),ihel(decay.size()+1);
    int itemp(1);
    do {
      --ix;
      if(ix!=_formpart[mode][iy]) {
	itemp*=decay[ix]->data().iSpin();
	constants[ix]=itemp;
      }
    }
    while(ix!=0);
    constants[decay.size()]=1;
    if(_formpart[mode][iy]!=decay.size())
      constants[_formpart[mode][iy]]=constants[_formpart[mode][iy]+1];
    // calculate the current
    vector<LorentzPolarizationVectorE>
      curr=_current[_currentmapA[mode][iy]]->
      current(_currentmapB[mode][iy],ichan,scale,cpart,meopt);
    pre = (pow(part.mass()/scale,int(cpart.size()-2)));
    // loop over the helicities to calculate the matrix element
    ihel[0]=0;
    for(chel=0;chel<curr.size();++chel) {
      for(ix=decay.size();ix>0;--ix) {
	if(ix!=_formpart[mode][iy]+1)
	  ihel[ix]=(chel%constants[ix-1])/constants[ix];
      }
      for(fhel=0;fhel<form.size();++fhel) {
	ihel[_formpart[mode][iy]+1]=fhel;
	(*ME())(ihel) += Complex(pre*_CKMfact[mode][iy]*
				 form[fhel].dot(curr[chel])*SM().fermiConstant());
      }
    }
  }
  // perform the contraction
  return 0.5*(ME()->contract(_rho)).real();
}
  
void ScalarMesonFactorizedDecayer::findModes(unsigned int imode,
					     vector<tPDVector> & particles,
					     vector<unsigned int> & loc,
					     vector<bool> & cc) {
  unsigned int ix,iy,nfound,iz;
  // resize the vectors
  loc.clear();cc.clear();
  // get the id's for the mode
  vector<int> id,idbar;
  int idtemp; bool found;
  for(ix=0;ix<particles[imode].size();++ix) {
    id.push_back(particles[imode][ix]->id());
    if(particles[imode][ix]->CC()) idbar.push_back(particles[imode][ix]->CC()->id());
    else                           idbar.push_back(id[ix]);
  }
  vector<bool> done(id.size(),false);
  // loop over the modes
  for(ix=0;ix<particles.size();++ix) {
    if(ix==imode||particles[ix].empty()) continue;
    assert(!particles[ix].empty());
    assert(particles[ix][0]);
    // the particle mode
    if(particles[ix][0]->id()==id[0]&&particles[ix].size()==id.size()) {
      nfound=1;
      for(iy=0;iy<id.size();++iy){done[iy]=false;}
      for(iy=1;iy<id.size();++iy) {
	idtemp=particles[ix][iy]->id();
	iz=1;
	found=false;
	do {
	  if(idtemp==id[iz]&&!done[iz]) {
	    done[iz]=true;
	    found=true;
	  }
	  ++iz;
	}
	while(iz<id.size()&&!found);
	if(found) ++nfound;
      }
      if(nfound==id.size()) {
	cc.push_back(false);
	loc.push_back(ix);
      }
    }
    // the charge conjugate mode
    if(particles[ix][0]->id()==idbar[0]&&particles[ix].size()==idbar.size()) {
      nfound=1;
      for(iy=0;iy<idbar.size();++iy) done[iy]=false;
      for(iy=1;iy<idbar.size();++iy) {
	idtemp=particles[ix][iy]->id();
	iz=1;
	found=false;
	do {
	  if(idtemp==idbar[iz]&&!done[iz]) {
	    done[iz]=true;
	    found=true;
	  }
	  ++iz;
	}
	while(iz<idbar.size()&&!found);
	if(found){++nfound;}
      }
      if(nfound==idbar.size()){cc.push_back(false);loc.push_back(ix);}
    }
  }
}

void ScalarMesonFactorizedDecayer::dataBaseOutput(ofstream & output,
						  bool header) const {
  unsigned int ix;
  if(header) output << "update decayers set parameters=\"";
  DecayIntegrator::dataBaseOutput(output,false);
  output << "newdef " << name() << ":a1Bottom "  << _a1b << "\n";
  output << "newdef " << name() << ":a2Bottom "  << _a2b << "\n";
  output << "newdef " << name() << ":a1Charm "   << _a1c << "\n";
  output << "newdef " << name() << ":a2Charm "   << _a2c << "\n";
  for(ix=0;ix<_current.size();++ix) {
    _current[ix]->dataBaseOutput(output,false,true);
    output << "insert " << name() << ":Currents " << ix << " " 
	   << _current[ix]->name() << " \n";
  }
  for(ix=0;ix<_form.size();++ix) {
    _form[ix]->dataBaseOutput(output,false,true);
    output << "insert " << name() << ":FormFactors " << ix << " " 
	   << _form[ix]->name() << " \n";
  }
  for(ix=0;ix<_wgtloc.size();++ix) {
    output << "insert " << name() << ":WeightLocation " << ix << " " 
	   << _wgtloc[ix] << "\n";
  }
  for(ix=0;ix<_wgtmax.size();++ix) {
    output << "insert " << name() << ":MaximumWeight "  << ix << " " 
	   << _wgtmax[ix] << "\n";
  }
  for(ix=0;ix<_weights.size();++ix) {
    output << "insert " << name() << ":Weights "        << ix << " " 
	   << _weights[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
