// -*- C++ -*-
//
// ScalarMesonFactorizedDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
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
  DecayIntegrator2::rebind(trans);
}

IVector ScalarMesonFactorizedDecayer::getReferences() {
  IVector ret = DecayIntegrator2::getReferences();
  ret.push_back(_ckm);
  return ret;
}

void ScalarMesonFactorizedDecayer::doinit() {
  DecayIntegrator2::doinit();
  // get the ckm object
  _ckm=dynamic_ptr_cast<Ptr<StandardCKM>::pointer>(SM().CKM());
  if(!_ckm) throw InitException() << "ScalarMesonFactorizedDecayer::doinit() "
  				  << "the CKM object must be the Herwig one"
  				  << Exception::runerror;
  // get the CKM matrix (unsquared for interference)
  Complex ckmmat[3][3];
  vector< vector<Complex > > CKM(_ckm->getUnsquaredMatrix(SM().families()));
  for(unsigned int ix=0;ix<3;++ix) {
    for(unsigned int iy=0;iy<3;++iy){
      ckmmat[ix][iy]=CKM[ix][iy];
    }
  }
  // make sure the currents and form factors got initialised
  for(unsigned int ix=0;ix<_current.size();++ix)
    _current[ix]->init();
  for(unsigned int ix=0;ix<_form.size();++ix)
    _form[ix]->init();
  // find all the possible modes
  vector<unsigned int> tformmap[2],tcurrmap[2];
  vector<int> inquark,outquark,currq,curra;
  tPDVector incoming;
  vector<tPDVector> outgoing;
  // loop over the modes in the form factors and currents
  for(unsigned int iform=0;iform<_form.size();++iform) {
    for(unsigned int ix=0;ix<_form[iform]->numberOfFactors();++ix) {
      // information from the form-factor
      int id0,id1,jspin,spect,inq,outq;
      _form[iform]->particleID(ix,id0,id1);
      _form[iform]->formFactorInfo(ix,jspin,spect,inq,outq);
      // particles from the form factor
      tPDPtr in  = getParticleData(id0);
      tPDPtr out = getParticleData(id1);
      // charge of the decay products
      int Wcharge = in->iCharge()-out->iCharge();
      // max mass for the particles in the current
      Energy min = in->massMax()-out->massMin();
      for(unsigned int icurr=0;icurr<_current.size();++icurr) {
  	for(unsigned int iy=0;iy<_current[icurr]->numberOfModes();++iy) {
	  // get the particles from the current
	  int iq,ia;
	  _current[icurr]->decayModeInfo(iy,iq,ia);
	  tPDVector ptemp=_current[icurr]->particles(Wcharge,iy,iq,ia);
	  tPDVector outV = {out};
	  outV.insert(std::end(outV), std::begin(ptemp), std::end(ptemp));
	  // create the mode
	  PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(in,outV,1.));
	  // create the first piece of the channel
	  PhaseSpaceChannel channel((PhaseSpaceChannel(mode),0,1)); 	  
   	  Energy minb=ZERO;
   	  for(unsigned int iz=0;iz<ptemp.size();++iz)
	    minb += ptemp[iz]->massMin();
  	  // add this mode to the list
   	  if(outV.size()>1&&minb<min&&
  	     (Wcharge!=0||(Wcharge==0&&((inq>0&&inq%2!=iq%2)||
   					(inq<0&&abs(inq)%2!=abs(ia)%2))))) {
   	    tformmap[0].push_back(iform);tformmap[1].push_back(ix);
   	    tcurrmap[0].push_back(icurr);tcurrmap[1].push_back(iy);
 	    incoming.push_back( in );
 	    outgoing.push_back(outV);
 	    inquark.push_back(inq);outquark.push_back(outq);
 	    currq.push_back(iq);curra.push_back(ia);
 	  }
   	  // if the meson in the current is neutral try the CC mode
   	  if(Wcharge==0&&iq!=-ia&&((inq>0&&inq%2!=iq%2)||
   				   (inq<0&&abs(inq)%2!=abs(ia)%2))) {
   	    // get the particles from the current
   	    tPDVector ptemp=_current[icurr]->particles(Wcharge,iy,-ia,-iq);
	    outV = {out};
	    outV.insert(std::end(outV), std::begin(ptemp), std::end(ptemp));
   	    minb=ZERO;
   	    for(unsigned int iz=0;iz<ptemp.size();++iz)
	      minb+=ptemp[iz]->massMin();
   	    if(outV.size()>1&&minb<min) {
   	      tformmap[0].push_back(iform);tformmap[1].push_back(ix);
   	      tcurrmap[0].push_back(icurr);tcurrmap[1].push_back(iy);
	      incoming.push_back(in);
   	      outgoing.push_back(outV);
   	      inquark.push_back(inq);outquark.push_back(outq);
   	      currq.push_back(-ia);curra.push_back(-iq);
   	    }
   	  }
  	}
      }
    }
  }
  // loop over the modes and find the dupliciates
  static const double ort(sqrt(0.5));
  for(unsigned int ix=0;ix<outgoing.size();++ix) {
    while (true) {
      if(outgoing[ix].empty()) break;
      vector<bool> modecc;
      vector<unsigned int> modeloc;
      findModes(ix,incoming,outgoing,modeloc,modecc);
      // if more than two outgoing only allow one diagram
      if ( outgoing[ix].size()>2 && !modeloc.empty() ) break;
      // create the mode and set the particles as for the first instance
      PhaseSpaceModePtr mode=new_ptr(PhaseSpaceMode(incoming[ix],outgoing[ix],1.));
      PhaseSpaceChannel channel((PhaseSpaceChannel(mode),0,1));
      Energy min = incoming[ix]->massMax()-outgoing[ix][0]->massMin();
      int Wcharge = incoming[ix]->iCharge()-outgoing[ix][0]->iCharge();
      bool done = _current[tcurrmap[0][ix]]->
	createMode(Wcharge,tcPDPtr(),IsoSpin::IUnknown,IsoSpin::I3Unknown,
		   tcurrmap[1][ix],mode,1,0,channel,min);
      if(!done) throw InitException() << "Failed to construct mode in "
   				      << "ScalarMesonFactorizedDecayer::doinit()." 
   				      << Exception::abortnow;
      // set the parameters for the additional modes
      vector<unsigned int> tformpart(1,0),ttform[2],ttcurr[2];
      ttform[0].push_back(tformmap[0][ix]);ttform[1].push_back(tformmap[1][ix]);
      ttcurr[0].push_back(tcurrmap[0][ix]);ttcurr[1].push_back(tcurrmap[1][ix]);
      int id    = outgoing[ix][0]->id();
      int idbar = outgoing[ix][0]->CC() ? outgoing[ix][0]->CC()->id() : id;
      for(unsigned int iy=0;iy<modeloc.size();++iy) {
   	ttform[0].push_back(tformmap[0][modeloc[iy]]);
   	ttform[1].push_back(tformmap[1][modeloc[iy]]);
   	ttcurr[0].push_back(tcurrmap[0][modeloc[iy]]);
   	ttcurr[1].push_back(tcurrmap[1][modeloc[iy]]);
   	unsigned int iz=0;
   	do {
   	  if(( modecc[iy]&&outgoing[modeloc[iy]][iz]->id()==idbar)||
   	     (!modecc[iy]&&outgoing[modeloc[iy]][iz]->id()==id))
   	    tformpart.push_back(iz);
   	  ++iz;
   	}
   	while(tformpart.size()!=iy+2&&iz<2);
      }
      // calculate ckm factors
      vector<Complex> tCKM;
      for(unsigned int iy=0;iy<ttcurr[0].size();++iy) {
  	// get the quarks involved in the process
	int iq,ia,inq,outq;
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
	int id0,id1;
	_form[ttform[0][iy]]->particleID(ttform[1][iy],id0,id1);
  	int Wcharge = getParticleData(id0)->iCharge()-getParticleData(id1)->iCharge();
   	Complex ckm=1.;
   	if(Wcharge!=0) {
   	  if(abs(iq)%2==0)  ckm *= conj(ckmmat[abs(iq)/2-1][(abs(ia)-1)/2]);
   	  else              ckm *= conj(ckmmat[abs(ia)/2-1][(abs(iq)-1)/2]);
   	  if(abs(inq)%2==0) ckm *= ckmmat[abs(inq)/2-1][(abs(outq)-1)/2];
   	  else              ckm *= ckmmat[abs(outq)/2-1][(abs(inq)-1)/2];
   	  if(abs(inq)==5)   ckm*=_a1b;
   	  else              ckm*=_a1c;
   	}
   	else {
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
      if(outgoing[ix][0]->id()==outgoing[ix][1]->id()&&outgoing[ix].size()==2) {
  	unsigned int isize=ttcurr[0].size();
   	for(unsigned int iy=0;iy<isize;++iy) {
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
      double maxweight(0.);
      if(_wgtmax.size()>numberModes()) maxweight=_wgtmax[numberModes()];
      // the weights for the channels
      vector<double> channelwgts;
     if(_wgtloc.size()>numberModes()&&
 	 _wgtloc[numberModes()]+mode->channels().size()<=_weights.size()) {
       vector<double>::iterator start=_weights.begin()+_wgtloc[numberModes()];
       vector<double>::iterator end  = start+mode->channels().size();
 	channelwgts=vector<double>(start,end);
     }
     else {
 	channelwgts.resize(mode->channels().size(),1./(mode->channels().size()));
     }
     // don't need channels for two body decays
     if(outgoing[ix].size()==2) {
       channelwgts.clear();
       mode = new_ptr(PhaseSpaceMode(incoming[ix],outgoing[ix],maxweight));
     }
     else {
       mode->maxWeight(maxweight);
       mode->setWeights(channelwgts);
     }
     addMode(mode);
     // resize the duplicate modes to remove them
     for(unsigned int iy=0;iy<modeloc.size();++iy) outgoing[modeloc[iy]] = tPDVector();
     break;
    }
  }
}

void ScalarMesonFactorizedDecayer::doinitrun() {
  unsigned int ix,iy;
  for(ix=0;ix<_current.size();++ix) _current[ix]->initrun();
  for(ix=0;ix<_form.size();++ix)    _form[ix]->initrun();
  DecayIntegrator2::doinitrun();
  if(initialize()) {
    _weights.clear();
    _wgtloc.clear();
    _wgtmax.clear();
    for(ix=0;ix<numberModes();++ix) {
      _wgtmax.push_back(mode(ix)->maxWeight());
      _wgtloc.push_back(_weights.size());
      for(iy=0;iy<mode(ix)->channels().size();++iy) {
	_weights.push_back(mode(ix)->channels()[iy].weight());
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
  if(parent->CC())  id0bar = parent->CC()->id();
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
  cc=false;
  unsigned int ix=0;
  do {
    // particle mode
    if(id0==mode(ix)->incoming().first->id()&&
       ids.size()==mode(ix)->outgoing().size()) {
      unsigned int nfound=0;
      vector<bool> done(ids.size(),false);
      for(unsigned int iy=0;iy<ids.size();++iy) {
   	int idtemp = mode(ix)->outgoing()[iy]->id();
	unsigned int iz=0;
	bool found=false;
	do{
	  if(idtemp==ids[iz]&&!done[iz]) {
	    done[iz]=true;
	    found=true;
	  }
	  ++iz;
	}
	while(iz<ids.size()&&!found);
	if(found) ++nfound;
      }
      if(nfound==ids.size()) {
	cc=false;
	imode=ix;
      }
    }
    // CC mode
    if(id0bar==mode(ix)->incoming().first->id()&&
       ids.size()==mode(ix)->outgoing().size()) {
      unsigned int nfound=0;
      vector<bool> done(ids.size(),false);
      for(unsigned int iy=0;iy<idbars.size();++iy) {
  	int idtemp=mode(ix)->outgoing()[iy]->id();
	unsigned int iz=0;
	bool found=false;
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
    throw DecayIntegrator2Error() << "Unable to find the mode " << mode << " in " 
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
DescribeClass<ScalarMesonFactorizedDecayer,DecayIntegrator2>
describeHerwigScalarMesonFactorizedDecayer("Herwig::ScalarMesonFactorizedDecayer", "HwSMDecay.so");

void ScalarMesonFactorizedDecayer::Init() {

  static ClassDocumentation<ScalarMesonFactorizedDecayer> documentation
    ("The ScalarMesonFactorizedDecayer class is designed for the weak decay of"
     " scalar mesons using the factorization approximation.");

  static RefVector<ScalarMesonFactorizedDecayer,WeakCurrent> interfaceCurrents
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
 
void ScalarMesonFactorizedDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // set up the spin information for the decay products
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
  					incoming,true);
  // get the wavefunctions of the decay products
  for(unsigned int ix=0;ix<decay.size();++ix) {
    switch(decay[ix]->dataPtr()->iSpin()) {
    case PDT::Spin0:
      ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
      break;
    case PDT::Spin1:
      VectorWaveFunction::constructSpinInfo(_vectors[ix],decay[ix],outgoing,
					    true,false);
      break;
    case PDT::Spin2:
      TensorWaveFunction::constructSpinInfo(_tensors[ix],decay[ix],outgoing,
					    true,false);
      break;
    default:
      assert(false);
    }
  }
}

double ScalarMesonFactorizedDecayer::me2(const int ichan, const Particle & part,
					 const tPDVector & outgoing,
					 const vector<Lorentz5Momentum> & momenta,
					 MEOption meopt) const {
  if(!ME()) {
    // create the matrix element
    vector<PDT::Spin> spin;
    for(unsigned int ix=0;ix<outgoing.size();++ix)
      spin.push_back(outgoing[ix]->iSpin());
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,spin)));
  }
  // initialisation
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&part),incoming);
    _vectors.resize(outgoing.size());
    _tensors.resize(outgoing.size());
  }
  // get the wavefunctions of the decay products
  for(unsigned int ix=0;ix<outgoing.size();++ix) {
    switch(outgoing[ix]->iSpin()) {
    case PDT::Spin0:
      break;
    case PDT::Spin1:
      _vectors[ix].resize(3);
      for(unsigned int ihel=0;ihel<3;++ihel)
	_vectors[ix][ihel] = HelicityFunctions::polarizationVector(-momenta[ix],ihel,
								   Helicity::outgoing);
      break;
    case PDT::Spin2:
      {
	TensorWaveFunction twave(momenta[ix],outgoing[ix],Helicity::outgoing);
	_tensors[ix].resize(5);
	for(unsigned int ihel=0;ihel<5;++ihel) {
	  twave.reset(ihel);
	  _tensors[ix][ihel] = twave.wave();
	}
      }
    break;
    default:
      assert(false);
    }
  }
  ME()->zero();
  // find the mode
  unsigned int mode(imode());
  int id0(part.id());
  Complex ii(0.,1.);
  vector<unsigned int> ihel(outgoing.size());
  // loop over the different diagrams
  Energy MP(part.mass()),scale;
  double pre;
  for(unsigned int iy=0;iy<_CKMfact[mode].size();++iy) {
    Energy MV = momenta[_formpart[mode][iy]].mass();
    int id1   = outgoing[_formpart[mode][iy]]->id();
    int id0t,id1t;
    _form[_formmapA[mode][iy]]->particleID(_formmapB[mode][iy],id0t,id1t);
    bool cc(id0t!=id0);
    // calculate the form-factor part
    vector<LorentzPolarizationVectorE> form;
    Lorentz5Momentum q   = part.momentum()-momenta[_formpart[mode][iy]];
    q.rescaleMass();
    Lorentz5Momentum sum = part.momentum()+momenta[_formpart[mode][iy]];
    sum.rescaleMass();
    Energy2 q2=q.mass2();
    if(outgoing[_formpart[mode][iy]]->iSpin()==1) {
      Complex fp,f0;
      _form[_formmapA[mode][iy]]->ScalarScalarFormFactor(q2,_formmapB[mode][iy],
   							 id0,id1,MP,MV,f0,fp);
      pre=(MP*MP-MV*MV)/q2;
      form.push_back(fp*sum+pre*(f0-fp)*q);
    }
    else if(outgoing[_formpart[mode][iy]]->iSpin()==3) {
      Energy msum  = MP+MV;
      Energy mdiff = MP-MV;
      Complex A0,A1,A2,V;
       _form[_formmapA[mode][iy]]->ScalarVectorFormFactor(q2,_formmapB[mode][iy],id0,
							  id1,MP,MV,A0,A1,A2,V);
       if(cc) V=-V;
       Complex A3 = 0.5/MV*(msum*A1-mdiff*A2);
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
					   momenta[_formpart[mode][iy]])); 
       }
    }
    else if(outgoing[_formpart[mode][iy]]->iSpin()==5) {
      Complex k;
      complex<InvEnergy2> h,bp,bm;
      _form[_formmapA[mode][iy]]->ScalarTensorFormFactor(q2,_formmapB[mode][iy],
  							 id0,id1,MP,MV,h,k,bp,bm);
      if(cc) h=-h;
      // compute the hadron currents
      for(unsigned int ix=0;ix<5;++ix) {
	LorentzPolarizationVectorE dotv =
	  _tensors[_formpart[mode][iy]][ix]*part.momentum();
   	complex<Energy2> dot = dotv*part.momentum();
   	form.push_back(ii*h*epsilon(dotv,sum,q)-k*dotv
   		       -bp*dot*sum-bm*dot*q);
      }
    }
    // find the particles for the current
    tPDVector cpart;
    vector<Lorentz5Momentum> cmom;
    for(unsigned int ix=0;ix<outgoing.size();++ix) {
      if(ix!=_formpart[mode][iy]) {
	cpart.push_back(outgoing[ix]);
	cmom .push_back(momenta[ix]);
      }
    }
    unsigned int ix=outgoing.size();
    vector<unsigned int> constants(outgoing.size()+1),ihel(outgoing.size()+1);
    int itemp(1);
    do {
      --ix;
      if(ix!=_formpart[mode][iy]) {
	itemp*=outgoing[ix]->iSpin();
	constants[ix]=itemp;
      }
    }
    while(ix!=0);
    constants[outgoing.size()]=1;
    if(_formpart[mode][iy]!=outgoing.size())
      constants[_formpart[mode][iy]]=constants[_formpart[mode][iy]+1];
    // calculate the current
    vector<LorentzPolarizationVectorE>
      curr=_current[_currentmapA[mode][iy]]->
      current(tcPDPtr(),IsoSpin::IUnknown,IsoSpin::I3Unknown,
	      _currentmapB[mode][iy],ichan,scale,cpart,cmom,meopt);
    pre = (pow(part.mass()/scale,int(cpart.size()-2)));
    // loop over the helicities to calculate the matrix element
    ihel[0]=0;
    for(unsigned int chel=0;chel<curr.size();++chel) {
      for(ix=outgoing.size();ix>0;--ix) {
	if(ix!=_formpart[mode][iy]+1)
	  ihel[ix]=(chel%constants[ix-1])/constants[ix];
      }
      for(unsigned int fhel=0;fhel<form.size();++fhel) {
	ihel[_formpart[mode][iy]+1]=fhel;
	(*ME())(ihel) +=pre*_CKMfact[mode][iy]*
	  form[fhel].dot(curr[chel])*SM().fermiConstant();
      }
    }
  }
  // perform the contraction
  return 0.5*(ME()->contract(_rho)).real();
}
  
void ScalarMesonFactorizedDecayer::findModes(unsigned int imode,
					     tPDVector & incoming,
					     vector<tPDVector> & outgoing,
					     vector<unsigned int> & loc,
					     vector<bool> & cc) {
  // get the id's for the mode
  // incoming
  int id_in    = incoming[imode]->id();
  int idbar_in = incoming[imode]->CC() ?
    incoming[imode]->CC()->id() : incoming[imode]->id();
  // outgoing
  vector<int> id_out,idbar_out;
  for(unsigned int ix=0;ix<outgoing[imode].size();++ix) {
    id_out.push_back(outgoing[imode][ix]->id());
    if(outgoing[imode][ix]->CC())
      idbar_out.push_back(outgoing[imode][ix]->CC()->id());
    else
      idbar_out.push_back(id_out[ix]);
  }
  // loop over the modes
  for(unsigned int ix=0;ix<outgoing.size();++ix) {
    if(ix==imode||outgoing[ix].empty()) continue;
    assert(!outgoing[ix].empty());
    assert(incoming[ix]);
    // the particle mode
    if(incoming[ix]->id()==id_in&&outgoing[ix].size()==id_out.size()) {
      vector<bool> done(id_out.size(),false);
      unsigned int nfound = 0;
      for(unsigned int iy=0;iy<id_out.size();++iy) {
	int idtemp=outgoing[ix][iy]->id();
     	unsigned int iz(0);
     	bool found=false;
    	do {
	  if(idtemp==id_out[iz]&&!done[iz]) {
     	    done[iz]=true;
     	    found=true;
    	  }
    	  ++iz;
    	}
     	while(iz<id_out.size()&&!found);
     	if(found) ++nfound;
	if(nfound==id_out.size()) {
	  cc.push_back(false);
	  loc.push_back(ix);
        }
      }
    }
    // the charge conjugate mode
    if(incoming[ix]->id()==idbar_in&&outgoing[ix].size()==idbar_out.size()) {
      vector<bool> done(id_out.size(),false);
      unsigned int nfound = 0;
      for(unsigned int iy=0;iy<idbar_out.size();++iy) {
    	int idtemp=outgoing[ix][iy]->id();
	unsigned int iz(0);
   	bool found=false;
   	do {
     	  if(idtemp==idbar_out[iz]&&!done[iz]) {
     	    done[iz]=true;
     	    found=true;
     	  }
     	  ++iz;
     	}
     	while(iz<idbar_out.size()&&!found);
    	if(found) ++nfound;
      }
      if(nfound==idbar_out.size()) {
	cc.push_back(false);
	loc.push_back(ix);
      }
    }
  }
}

void ScalarMesonFactorizedDecayer::dataBaseOutput(ofstream & output,
						  bool header) const {
  unsigned int ix;
  if(header) output << "update decayers set parameters=\"";
  DecayIntegrator2::dataBaseOutput(output,false);
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
