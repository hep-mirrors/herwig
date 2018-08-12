// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BaryonFactorizedDecayer class.
//

#include "BaryonFactorizedDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorBarWaveFunction.h"
#include "Herwig/Decay/DecayVertex.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/Helicity/RSFermionSpinInfo.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

BaryonFactorizedDecayer::BaryonFactorizedDecayer() {
  // default values taken from PRD56, 2799
  _a1c= 1.1;
  _a2c=-0.5;
  _a1b= 1.0;
  _a2b= 0.28;
  // intermediates
  generateIntermediates(true);
}

void BaryonFactorizedDecayer::doinitrun() {
  _current->initrun();
  _form->initrun();
  DecayIntegrator2::doinitrun();
  _weights.clear();_wgtloc.clear();_wgtmax.clear();
  for(unsigned int ix=0;ix<numberModes();++ix) {
    _wgtmax.push_back(mode(ix)->maxWeight());
    _wgtloc.push_back(_weights.size());
    for(unsigned int iy=0;iy<mode(ix)->channels().size();++iy)
      _weights.push_back(mode(ix)->channels()[iy].weight());
  }
}

void BaryonFactorizedDecayer::doinit() {
  DecayIntegrator2::doinit();
  // get the CKM matrix (unsquared for interference)
  Complex ckmmat[3][3];
  vector< vector<Complex > > CKM(_theCKM->getUnsquaredMatrix(SM().families()));
  for(unsigned int ix=0;ix<3;++ix) {
    for(unsigned int iy=0;iy<3;++iy) {
      ckmmat[ix][iy]=CKM[ix][iy];
    }
  }
  // make sure the current and form factor got initialised
  _current->init();
  _form->init();
  // find all the possible modes
  vector<unsigned int> tformmap,tcurrmap;
  vector<int> inquark,outquark,currq,curra;
  tPDVector incoming;
  vector<tPDVector> outgoing;
  for(unsigned int iform=0;iform<_form->numberOfFactors();++iform) {
    // particles from the form factor
    int id0,id1;
    _form->particleID (iform,id0,id1);
    int spect1,spect2,inq,outq,ispin,ospin;
    _form->formFactorInfo(iform,ispin,ospin,spect1,spect2,inq,outq);
    // particles from the form factor
    tPDPtr in  = getParticleData(id0);
    tPDPtr out = getParticleData(id1);
    // the charge of the decay products
    int Wcharge = in->iCharge()-out->iCharge();
    // max mass for the particles in the current
    Energy min = in->massMax()-out->massMin();
    for(unsigned int icurr=0;icurr<_current->numberOfModes();++icurr) {
      // get the particles from the current
      int iq,ia;
      _current->decayModeInfo(icurr,iq,ia);
      tPDVector ptemp=_current->particles(Wcharge,icurr,iq,ia);
      tPDVector outV = {out};
      outV.insert(std::end(outV), std::begin(ptemp), std::end(ptemp));
      Energy minb=ZERO;
      for(unsigned int iz=0;iz<ptemp.size();++iz) minb+=ptemp[iz]->massMin();
      // valid mode
      if(outV.size()>1&&minb<min&&
	 (Wcharge!=0||(Wcharge==0&&((inq>0&&inq%2!=iq%2)||
				    (inq<0&&abs(inq)%2!=abs(ia)%2))))) {
	tformmap.push_back(iform);tcurrmap.push_back(icurr);
	incoming.push_back(in);
	outgoing.push_back(outV);
	inquark.push_back(inq);outquark.push_back(outq);
	currq.push_back(iq);curra.push_back(ia);
      }
      // if the meson is neutral try the CC mode
      if(Wcharge==0&&iq!=-ia&&((inq>0&&inq%2!=iq%2)||
			       (inq<0&&abs(inq)%2!=abs(ia)%2))) {
	ptemp=_current->particles(Wcharge,icurr,-ia,-iq);
	tPDVector outV = {out};
	outV.insert(std::end(outV), std::begin(ptemp), std::end(ptemp));
	Energy minb=ZERO;
	for(unsigned int iz=0;iz<ptemp.size();++iz) minb+=ptemp[iz]->massMin();
	if(outV.size()>1&&minb<min) {
  	  tformmap.push_back(iform);tcurrmap.push_back(icurr);
	  incoming.push_back(in);
  	  outgoing.push_back(outV);
  	  inquark.push_back(inq);outquark.push_back(outq);
  	  currq.push_back(-ia);curra.push_back(-iq);
  	}
      }
    }
  }
  _formmap.clear();
  _currentmap.clear();
  // loop over the modes and find the dupliciates
  for(unsigned int ix=0;ix<outgoing.size();++ix) {
    while(true) {
      if ( outgoing[ix].empty() ) break;
      vector<unsigned int> modeloc;
      vector<bool> modecc;
      findModes(ix,incoming,outgoing,modeloc,modecc);
      // if more than two outgoing particles only allow one diagram
      if ( outgoing[ix].size() > 2 && !modeloc.empty() ) {break;}
      // create the mode and set the particles as for the first instance
      PhaseSpaceModePtr mode=new_ptr(PhaseSpaceMode(incoming[ix],outgoing[ix],1.));
      PhaseSpaceChannel channel((PhaseSpaceChannel(mode),0,1));
      Energy min = incoming[ix]->massMax()-outgoing[ix][0]->massMin();
      int Wcharge = incoming[ix]->iCharge()-outgoing[ix][0]->iCharge();
      bool done = _current->createMode(Wcharge,tcPDPtr(),IsoSpin::IUnknown,IsoSpin::I3Unknown,
				       tcurrmap[ix],mode,1,0,channel,min);
      if(!done){throw InitException() << "Failed to construct mode in "
				      << "BaryonFactorizedDecayer::doinit()." 
				      << Exception::abortnow;}
      // set the parameters for the additional modes
      vector<unsigned int>ttform,ttcurr;
      ttform.push_back(tformmap[ix]);ttcurr.push_back(tcurrmap[ix]);
      for(unsigned int iy=0;iy<modeloc.size();++iy) {
	ttform.push_back(tformmap[modeloc[iy]]);
	ttcurr.push_back(tcurrmap[modeloc[iy]]);
      }
      vector<Complex> tCKM; Complex ckm;
      for(unsigned int iy=0;iy<ttcurr.size();++iy) {
	// get the quarks involved in the process
	int iq,ia,inq,outq;
	if(iy==0) {
	  iq=currq[ix];
	  ia=curra[ix];
	  inq=inquark[ix];
	  outq=outquark[ix];
	}
	else {
	  if(!modecc[iy-1]) {
	    iq=currq[modeloc[iy-1]];
	    ia=curra[modeloc[iy-1]];
	    inq=inquark[modeloc[iy-1]];
	    outq=outquark[modeloc[iy-1]];
	  }
	  else {
	    ia=-currq[modeloc[iy-1]];
	    iq=-curra[modeloc[iy-1]];
	    inq=-inquark[modeloc[iy-1]];
	    outq=-outquark[modeloc[iy-1]];
	  }
 	}
	int id0,id1;
	_form->particleID(ttform[iy],id0,id1);
	int Wcharge = getParticleData(id0)->iCharge()-getParticleData(id1)->iCharge();
	Complex ckm=1.;
	if(Wcharge!=0) {
	  if(abs(iq)%2==0){ckm *= conj(ckmmat[abs(iq)/2-1][(abs(ia)-1)/2]);}
	  else{ckm *= conj(ckmmat[abs(ia)/2-1][(abs(iq)-1)/2]);}
	  if(abs(inq)%2==0){ckm *= ckmmat[abs(inq)/2-1][(abs(outq)-1)/2];}
	  else{ckm *= ckmmat[abs(outq)/2-1][(abs(inq)-1)/2];}
	  if(abs(inq)==5){ckm*=_a1b;}
	  else{ckm*=_a1c;}
	}
	else {
	  if(inq>0) {
	    if(abs(inq)%2==0){ckm *= ckmmat[abs(inq)/2-1][(abs(iq)-1)/2];}
	    else{ckm *= ckmmat[abs(iq)/2-1][(abs(inq)-1)/2];}
	    if(abs(outq)%2==0)
	      {ckm *= conj(ckmmat[abs(outq)/2-1][(abs(ia)-1)/2]);}
	    else{ckm *= conj(ckmmat[abs(ia)/2-1][(abs(outq)-1)/2]);}
	  }
	  else {
	    if(abs(inq)%2==0){ckm *= ckmmat[abs(inq)/2-1][(abs(ia)-1)/2];}
	    else{ckm *= ckmmat[abs(ia)/2-1][(abs(inq)-1)/2];}
	    if(abs(outq)%2==0)
	      {ckm *= conj(ckmmat[abs(outq)/2-1][(abs(iq)-1)/2]);}
	    else{ckm *= conj(ckmmat[abs(iq)/2-1][(abs(outq)-1)/2]);}
	  }
	  if(abs(inq)==5){ckm*=_a2b;}
	  else{ckm*=_a2c;}
	}
	if((abs(inq)%2==0&&inq<0)||(abs(inq)%2!=0&&inq>0)){ckm=conj(ckm);}
	tCKM.push_back(ckm);
      }
      // add the parameters for the mode to the list
      _currentmap.push_back(ttcurr);
      _formmap.push_back(ttform);
      _factCKM.push_back(tCKM);
      double maxweight(0.);
      // add the mode to the list
      if(_wgtmax.size()>numberModes()) maxweight=_wgtmax[numberModes()];
      // the weights for the channel
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
      for(unsigned int iy=0;iy<modeloc.size();++iy) outgoing[modeloc[iy]]=tPDVector();
      break;
    }
  }
}

bool BaryonFactorizedDecayer::accept(tcPDPtr parent, const tPDVector & children) const {
  bool allowed=false;
  unsigned int iform(0),ix;
  int idin(parent->id()),ibaryon,foundb,id0,id1;
  vector<int> idall,idother;
  tPDVector::const_iterator pit  = children.begin();
  tPDVector::const_iterator pend = children.end();
  for( ; pit!=pend;++pit){idall.push_back((**pit).id());}
  // loop over the particles in the form factor
  do {
    _form->particleID(iform,id0,id1);
    ibaryon=0;
    if(id0==idin){ibaryon=id1;}
    else if(id0==-idin){ibaryon=-id1;}
    if(ibaryon!=0) {
      foundb=false;
      idother.clear();
      for(ix=0;ix<idall.size();++ix) {
	if(idall[ix]==ibaryon){foundb=true;}
	else{idother.push_back(idall[ix]);}
      }
      if(foundb){allowed=_current->accept(idother);}
    }
    ++iform;
  }
  while(!allowed&&iform<_form->numberOfFactors());
  return allowed;
}

int BaryonFactorizedDecayer::modeNumber(bool & cc,tcPDPtr parent,
					const tPDVector & children) const {
  unsigned int ix,iy;
  int idin(parent->id()),ibaryon,foundb,id0,id1,icurr(-1),iform(0);
  vector<int> idall,idother;
  tPDVector::const_iterator pit  = children.begin();
  tPDVector::const_iterator pend = children.end();
  for( ; pit!=pend;++pit){idall.push_back((**pit).id());}
  // loop over the particles in the form factor
  do
    {
      _form->particleID(iform,id0,id1);
      ibaryon=0;
      if(id0==idin){ibaryon=id1;}
      else if(id0==-idin){ibaryon=-id1;}
      ++iform;
      foundb=false;
      idother.clear();
      for(ix=0;ix<idall.size();++ix)
	{
	  if(idall[ix]==ibaryon){foundb=true;}
	  else{idother.push_back(idall[ix]);}
	}
      if(foundb){icurr=_current->decayMode(idother);}
    }
  while(icurr<0&&iform<int(_form->numberOfFactors()));
  // now find the mode
  int imode=-1;
  ix=0;
  --iform;
  do
    {
      for(iy=0;iy<_currentmap[ix].size();++iy)
	{if(int(_currentmap[ix][iy])==icurr&&int(_formmap[ix][iy])==iform){imode=ix;}}
      ++ix;
    }
  while(imode<0&&ix<numberModes());
  if(imode<0){throw DecayIntegrator2Error() << "Unable to find the mode in " 
					   << "BaryonFactorizedDecayer::decay()" 
					   << Exception::abortnow;}
  // generate the mode
  cc=id0!=idin;
  return imode;
}


void BaryonFactorizedDecayer::persistentOutput(PersistentOStream & os) const {
  os << _current << _form << _a1b << _a2b <<_a1c <<_a2c 
     << _currentmap << _formmap << _factCKM << _wgtloc << _wgtmax << _weights 
     << _theCKM;
}

void BaryonFactorizedDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _current >> _form >> _a1b >> _a2b >>_a1c >>_a2c 
     >> _currentmap >> _formmap >> _factCKM >> _wgtloc >> _wgtmax >> _weights 
     >> _theCKM;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<BaryonFactorizedDecayer,DecayIntegrator2>
describeHerwigBaryonFactorizedDecayer("Herwig::BaryonFactorizedDecayer", "HwBaryonDecay.so");

void BaryonFactorizedDecayer::Init() {

  static ClassDocumentation<BaryonFactorizedDecayer> documentation
    ("The BaryonFactorizedDecayer class combines the baryon form factor and a"
     " weak current to perform a decay in the naive factorization approximation.");

  static Reference<BaryonFactorizedDecayer,WeakCurrent> interfaceWeakCurrent
    ("Current",
     "The reference for the decay current to be used.",
     &BaryonFactorizedDecayer::_current, false, false, true, false, false);

  static ParVector<BaryonFactorizedDecayer,int> interfaceWeightLocation
    ("WeightLocation",
     "The locations of the weights for a given channel in the vector",
     &BaryonFactorizedDecayer::_wgtloc,
     0, 0, 0, 0, 10000, false, false, true);

  static ParVector<BaryonFactorizedDecayer,double> interfaceWeightMax
    ("MaximumWeight",
     "The maximum weight for a given channel.",
     &BaryonFactorizedDecayer::_wgtmax,
     0, 0, 0, 0., 100., false, false, true);

  static ParVector<BaryonFactorizedDecayer,double> interfaceWeights
    ("Weights",
     "The weights for the integration.",
     &BaryonFactorizedDecayer::_weights,
     0, 0, 0, 0., 1., false, false, true);

  static Reference<BaryonFactorizedDecayer,BaryonFormFactor> interfaceFormFactor
    ("FormFactor",
     "The form-factor",
     &BaryonFactorizedDecayer::_form, true, true, true, false, false);

  static Parameter<BaryonFactorizedDecayer,double> interfacea1Bottom
    ("a1Bottom",
     "The factorization paramter a_1 for decays of bottom baryons",
     &BaryonFactorizedDecayer::_a1b, 1., -10.0, 10.0,
     false, false, true);

  static Parameter<BaryonFactorizedDecayer,double> interfacea2Bottom
    ("a2Bottom",
     "The factorization paramter a_2 for decays of bottom baryons",
     &BaryonFactorizedDecayer::_a2b, 0.28, -10.0, 10.0,
     false, false, true);

  static Parameter<BaryonFactorizedDecayer,double> interfacea1Charm
    ("a1Charm",
     "The factorization paramter a_1 for decays of charm baryons",
     &BaryonFactorizedDecayer::_a1c, 1.1, -10.0, 10.0,
     false, false, true);

  static Parameter<BaryonFactorizedDecayer,double> interfacea2Charm
    ("a2Charm",
     "The factorization paramter a_2 for decays of charm baryons",
     &BaryonFactorizedDecayer::_a2c, -0.5, -10.0, 10.0,
     false, false, true);

  static Reference<BaryonFactorizedDecayer,StandardCKM> interfaceCKM
    ("CKM",
     "Reference to the Standard Model object",
     &BaryonFactorizedDecayer::_theCKM, false, false, true, false);
}

void BaryonFactorizedDecayer::
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

double BaryonFactorizedDecayer::me2(const int ichan, const Particle & part,
				    const tPDVector & outgoing,
				    const vector<Lorentz5Momentum> & momenta,
				    MEOption meopt) const {
  double me(0.);
  assert(part.dataPtr()->iSpin()==PDT::Spin1Half);
  if(outgoing[0]->iSpin()==PDT::Spin1Half)
    me=halfHalf(ichan,part,outgoing,momenta,meopt);
  else if(outgoing[0]->iSpin()==PDT::Spin3Half)
    me=halfThreeHalf(ichan,part,outgoing,momenta,meopt);
  else
    assert(false);
  return me;
}

// matrix element for a 1/2 -> 1/2 decay
double BaryonFactorizedDecayer::halfHalf(const int ichan, const Particle & part,
					 const tPDVector & outgoing,
					 const vector<Lorentz5Momentum> & momenta,
					 MEOption meopt) const {
  Energy scale;
  // extract the spins of the particles
  vector<PDT::Spin> spin;
  for(unsigned ix=0;ix<outgoing.size();++ix) 
    spin.push_back(outgoing[ix]->iSpin());
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half,spin)));
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
  }
  ME()->zero();
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
  int id0(part.id()),id1(outgoing[0]->id());
  // work out the value of q and calculate the form factors
  Lorentz5Momentum q(part.momentum()-momenta[0]);
  q.rescaleMass();
  Energy m0(part.mass()),m1(momenta[0].mass());
  Energy2 q2(q.mass2());
  Lorentz5Momentum sum(part.momentum()+momenta[0]);
  // calculate the baryonic part of the current for the decay
  double pre(0.);
  for(unsigned int mode=0;mode<_formmap[imode()].size();++mode) {
    Complex f1v,f2v,f3v,f1a,f2a,f3a;
    // calculate the form factor piece
    _form->SpinHalfSpinHalfFormFactor(q2,_formmap[imode()][mode],id0,id1,m0,m1,
				      f1v,f2v,f3v,f1a,f2a,f3a);
    Complex  left = f1v-f1a-f2v-double((m0-m1)/(m0+m1))*f2a;
    Complex right = f1v+f1a-f2v+double((m0-m1)/(m0+m1))*f2a;
    vector<LorentzPolarizationVectorE> baryon(4);
    for(unsigned int ix=0;ix<2;++ix) {
      for(unsigned int iy=0;iy<2;++iy) {
	LorentzPolarizationVectorE 
	  vtemp = _inHalf[ix].generalCurrent(_inHalfBar[iy],left,right);
	complex<Energy> vspin=_inHalf[ix].scalar(_inHalfBar[iy]);      
	complex<Energy> aspin=_inHalf[ix].pseudoScalar(_inHalfBar[iy]);
	// the momentum like pieces
 	if(part.id()>0) {
 	  vtemp+= (f2v*vspin+f2a*aspin)/(m0+m1)*sum;
 	  vtemp+= (f3v*vspin+f3a*aspin)/(m0+m1)*q;
 	}
 	else {
 	  vtemp+= (f2v*vspin-f2a*aspin)/(m0+m1)*sum;
 	  vtemp+= (f3v*vspin-f3a*aspin)/(m0+m1)*q;
 	}
 	if(part.id()>0) baryon[2*ix+iy]=vtemp;
 	else            baryon[2*iy+ix]=vtemp;
      }
    }
    // construct the weak current
    vector<LorentzPolarizationVectorE> hadron =
      _current->current(tcPDPtr(),IsoSpin::IUnknown,IsoSpin::I3Unknown,
			_currentmap[imode()][mode],ichan,scale,
			tPDVector(outgoing.begin()+1,outgoing.end()),
			vector<Lorentz5Momentum>(momenta.begin()+1,momenta.end()),meopt);
    pre=pow(part.mass()/scale,int(outgoing.size()-3));pre*=pre;
    vector<unsigned int> constants(outgoing.size()+1),ihel(outgoing.size()+1);
    int itemp=1;
    unsigned int ibar=0;
    for(int iz=int(outgoing.size()-1);iz>=0;--iz) {
      if(abs(outgoing[iz]->id())!=id1) {
 	itemp *= outgoing[iz]->iSpin();
 	constants[iz]=itemp;
      }
      else ibar=iz;
      constants[outgoing.size()]=1;
      constants[ibar]=constants[ibar+1];
    }
    for(unsigned int mhel=0;mhel<baryon.size();++mhel) {
      ihel[0     ]=mhel/2;
      ihel[ibar+1]=mhel%2;
      for(unsigned int lhel=0;lhel<hadron.size();++lhel) {
 	// map the index for the hadrons to a helicity state
 	for(unsigned int ix=outgoing.size();ix>0;--ix) {
 	  if(ix-1!=ibar){ihel[ix]=(lhel%constants[ix-1])/constants[ix];}}
 	(*ME())(ihel) += hadron[lhel].dot(baryon[mhel])*
 	  _factCKM[imode()][mode]*SM().fermiConstant();
      }
    }
  }
  // return the answer
  return 0.5*pre*(ME()->contract(_rho)).real();
}

// matrix element for a 1/2 -> 3/2 decay
double BaryonFactorizedDecayer::halfThreeHalf(const int ichan, const Particle & part,
					      const tPDVector & outgoing,
					      const vector<Lorentz5Momentum> & momenta,
					      MEOption meopt) const {
  // spins
  Energy scale;
  vector<PDT::Spin> spin(outgoing.size());
  for(unsigned int ix=0;ix<outgoing.size();++ix)
    spin[ix]=outgoing[ix]->iSpin();
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half,spin)));
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
  }
  ME()->zero();
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
  int id0(part.id()),id1(outgoing[0]->id());
  // work out the value of q and calculate the form factors
  Lorentz5Momentum q(part.momentum()-momenta[0]);
  q.rescaleMass();
  Energy m0(part.mass()),m1(momenta[0].mass());
  Energy2 q2(q.mass2());
  Lorentz5Momentum sum(part.momentum()+momenta[0]);
  InvEnergy ms(1./(m0+m1));
  InvEnergy2 ms2(ms*ms);
  double pre(0.);
  for(unsigned int mode=0;mode<_formmap[imode()].size();++mode) {
    // calculate the form factors
    Complex f1v,f2v,f3v,f4v,f1a,f2a,f3a,f4a;
    _form->SpinHalfSpinThreeHalfFormFactor(q2,_formmap[imode()][mode],id0,id1,m0,m1,
  					   f1v,f2v,f3v,f4v,f1a,f2a,f3a,f4a);
    complex<InvEnergy2> lS1,lS2,rS1,rS2;
    Complex left,right;
    complex<InvEnergy> lV,rV;
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
    // construct the vectors for the decay
    LorentzPolarizationVectorE baryon[4][2];
    for(unsigned int iya=0;iya<4;++iya) {
      for(unsigned int ixa=0;ixa<2;++ixa) {
	unsigned int ix,iy;
  	if(outgoing[0]->id()>0) {
	  ix=iya;
	  iy=ixa;
	}
	else {
	  ix=ixa;
	  iy=iya;
	}
  	// scalar like terms
	complex<Energy> lfact = _inHalf[iy].leftScalar( _inHalfBar[ix]);
	complex<Energy> rfact = _inHalf[iy].rightScalar(_inHalfBar[ix]);
	Complex       scalar1 = (lS1*lfact+rS1*rfact)*UnitRemoval::E;
	Complex       scalar2 = (lS2*lfact+rS2*rfact)*UnitRemoval::E;
	LorentzPolarizationVector  svec = _inHalf[iy].generalCurrent(_inHalfBar[ix],lV/ms,rV/ms)*ms;
	LorentzPolarizationVectorE tvec;
   	if(part.id()>0) {
   	  tvec=_inThreeHalfBar[ix].generalCurrent(_inHalf[iy],left,right);
   	}
   	else {
   	  tvec=_inThreeHalf[iy].generalCurrent(_inHalfBar[ix],left,right);
   	}
  	baryon[iya][ixa] = tvec+svec*UnitRemoval::E
  	  +scalar1*momenta[0]+scalar2*part.momentum();
      }
    }
    vector<LorentzPolarizationVectorE> hadron =
      _current->current(tcPDPtr(),IsoSpin::IUnknown,IsoSpin::I3Unknown,
			_currentmap[imode()][mode],ichan,scale,
			tPDVector(outgoing.begin()+1,outgoing.end()),
			vector<Lorentz5Momentum>(momenta.begin()+1,momenta.end()),meopt);
    // prefactor
    pre  = pow(part.mass()/scale,int(outgoing.size()-3));
    pre *= pre;
    // work out the mapping for the hadron vector
    vector<unsigned int> constants(outgoing.size()+1),ihel(outgoing.size()+1);
    int itemp = 1;
    int ibar  = 0;
    for(int ix=int(outgoing.size()-1);ix>=0;--ix) {
      if(abs(outgoing[ix]->id())!=id1) {
   	itemp*=outgoing[ix]->iSpin();
   	constants[ix]=itemp;
      }
      else{ibar=ix;}
    }
    constants[outgoing.size()]=1;
    constants[ibar]=constants[ibar+1];
    for(unsigned int iya=0;iya<4;++iya) {
      ihel[1]=iya;
      for(unsigned int ixa=0;ixa<2;++ixa) {
   	ihel[0]=ixa;
	for(unsigned int lhel=0;lhel<hadron.size();++lhel) {
	  // map the index for the hadrons to a helicity state
	  for(int ix=int(outgoing.size());ix>0;--ix)
	    {if(ix-1!=ibar){ihel[ix]=(lhel%constants[ix-1])/constants[ix];}}
	  (*ME())(ihel) += hadron[lhel].dot(baryon[iya][ixa])*
	    _factCKM[imode()][mode]*SM().fermiConstant();
	}
      }
    }
  }
  // return the answer
  return 0.5*pre*(ME()->contract(_rho)).real();
}

void BaryonFactorizedDecayer::findModes(unsigned int imode,
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
	int idtemp = outgoing[ix][iy]->id();
	unsigned int iz = 0;
	bool found = false;
	do {
	  if(idtemp==id_out[iz]&&!done[iz]) {
	    done[iz]=true;
	    found=true;
	  }
	  ++iz;
	}
	while(iz<id_out.size()&&!found);
	if(found) ++nfound;
      }
      if(nfound==id_out.size()) {
	cc.push_back(false);
	loc.push_back(ix);
      }
    }
    // the charge conjugate mode
    if(incoming[ix]->id()==idbar_in&&outgoing[ix].size()==idbar_out.size()) {
      vector<bool> done(id_out.size(),false);
      unsigned int nfound=0;
      for(unsigned int iy=0;iy<idbar_out.size();++iy) {
	int idtemp=outgoing[ix][iy]->id();
	unsigned int iz=0;
	bool found = false;
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

// output the setup information for the particle database
void BaryonFactorizedDecayer::dataBaseOutput(ofstream & output, bool header) const {
  unsigned int ix;
  if(header){output << "update decayers set parameters=\"";}
  DecayIntegrator2::dataBaseOutput(output,false);
  output << "newdef " << name() << ":a1Bottom "  << _a1b << "\n";
  output << "newdef " << name() << ":a2Bottom "  << _a2b << "\n";
  output << "newdef " << name() << ":a1Charm "   << _a1c << "\n";
  output << "newdef " << name() << ":a2Charm "   << _a2c << "\n";
  output << "newdef " << name() << ":CKM "       << _theCKM->name() << " \n";
  for(ix=0;ix<_wgtloc.size();++ix)
    {output << "insert " << name() << ":WeightLocation " << ix << " " 
	    << _wgtloc[ix] << "\n";}
  for(ix=0;ix<_wgtmax.size();++ix)
    {output << "insert " << name() << ":MaximumWeight "  << ix << " " 
	    << _wgtmax[ix] << "\n";}
  for(ix=0;ix<_weights.size();++ix)
    {output << "insert " << name() << ":Weights "        << ix << " " 
	    << _weights[ix] << "\n";}
  _current->dataBaseOutput(output,false,true);
  output << "newdef " << name() << ":Current " << _current->name() << " \n";
  _form->dataBaseOutput(output,false,true);
  output << "newdef " << name() << ":FormFactor " << _form->name() << " \n";
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}
