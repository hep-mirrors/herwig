// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BaryonFactorizedDecayer class.
//

#include "BaryonFactorizedDecayer.h"
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
  DecayIntegrator::doinitrun();
  _weights.clear();_wgtloc.clear();_wgtmax.clear();
  unsigned int ix,iy;
  for(ix=0;ix<numberModes();++ix) {
    _wgtmax.push_back(mode(ix)->maxWeight());
    _wgtloc.push_back(_weights.size());
    for(iy=0;iy<mode(ix)->numberChannels();++iy)
      _weights.push_back(mode(ix)->channelWeight(iy));
  }
}

void BaryonFactorizedDecayer::doinit() {
  DecayIntegrator::doinit();
  // get the CKM matrix elements
  unsigned int ix,iy,iz,iform,icurr;
  // get the CKM matrix (unsquared for interference)
  Complex ckmmat[3][3];
  vector< vector<Complex > > CKM(_theCKM->getUnsquaredMatrix(SM().families()));
  for(ix=0;ix<3;++ix){for(iy=0;iy<3;++iy){ckmmat[ix][iy]=CKM[ix][iy];}}
  // make sure the current and form factor got initialised
  _current->init();
  _form->init();
  // find all the possible modes
  vector<unsigned int> tformmap,tcurrmap;
  vector<int> inquark,outquark,currq,curra;
  vector<tPDVector> particles;
  tPDVector extpart,extpartb,ptemp;
  Energy min,minb;
  int iq,ia,spect1,spect2,inq,outq,id0,id1,Wcharge,ispin,ospin;
  for(iform=0;iform<_form->numberOfFactors();++iform)
    {
      // particles from the form factor
      extpart.resize(2);
      _form->particleID (iform,id0,id1);
      _form->formFactorInfo(iform,ispin,ospin,spect1,spect2,inq,outq);
      // particles from the form factor
      extpart[0]=getParticleData(id0);
      extpart[1]=getParticleData(id1);
      // the charge of the decay products
      Wcharge = extpart[0]->iCharge()-extpart[1]->iCharge();
      // max mass for the particles in the current
      min = extpart[0]->massMax()-extpart[1]->massMin();
      for(icurr=0;icurr<_current->numberOfModes();++icurr)
	{
	  extpart.resize(2);
	  // get the particles from the current
	  _current->decayModeInfo(icurr,iq,ia);
	  ptemp=_current->particles(Wcharge,icurr,iq,ia);
	  minb=ZERO;
	  for(iz=0;iz<ptemp.size();++iz)
	    {extpart.push_back(ptemp[iz]);minb+=ptemp[iz]->massMin();}
	  // valid mode
	  if(extpart.size()>2&&minb<min&&
	     (Wcharge!=0||(Wcharge==0&&((inq>0&&inq%2!=iq%2)||
					(inq<0&&abs(inq)%2!=abs(ia)%2)))))
	    {
	      tformmap.push_back(iform);tcurrmap.push_back(icurr);
	      particles.push_back(extpart);
	      inquark.push_back(inq);outquark.push_back(outq);
	      currq.push_back(iq);curra.push_back(ia);
	    }
	  // if the meson is neutral try the CC mode
	  if(Wcharge==0&&iq!=-ia&&((inq>0&&inq%2!=iq%2)||
				   (inq<0&&abs(inq)%2!=abs(ia)%2)))
	    {
	      extpart.resize(2);
	      ptemp=_current->particles(Wcharge,icurr,-ia,-iq);
	      minb=ZERO;
	      for(iz=0;iz<ptemp.size();++iz)
		{extpart.push_back(ptemp[iz]);minb+=ptemp[iz]->massMin();}
	      if(extpart.size()>2&&minb<min)
		{
		  tformmap.push_back(iform);tcurrmap.push_back(icurr);
		  particles.push_back(extpart);
		  inquark.push_back(inq);outquark.push_back(outq);
		  currq.push_back(-ia);curra.push_back(-iq);
		}
	    }
	}
    } 
  vector<bool> modecc;
  vector<unsigned int> modeloc,ttform,ttcurr;
  vector<Complex> tCKM; Complex ckm;
  bool done;
  DecayPhaseSpaceModePtr mode;
  DecayPhaseSpaceChannelPtr channel;
  vector<double>::iterator start,end;
  double maxweight;
  vector<double> channelwgts;
  _formmap.clear();
  _currentmap.clear();
  // loop over the modes and find the dupliciates
  for(ix=0;ix<particles.size();++ix)
    {
      while(true)
	{
	  if ( particles[ix].empty() ) {break;}
	  findModes(ix,particles,modeloc,modecc);
	  // if more than three particles only allow one diagram
	  if ( particles[ix].size() > 3 && !modeloc.empty() ) {break;}
	  // create the mode and set the particles as for the first instance
	  mode=new_ptr(DecayPhaseSpaceMode(particles[ix],this));
	  channel = new_ptr(DecayPhaseSpaceChannel(mode));
	  channel->addIntermediate(particles[ix][0],0,0.0,1,-1);
	  min = particles[ix][0]->massMax()-particles[ix][1]->massMin();
	  Wcharge = particles[ix][0]->iCharge()-particles[ix][1]->iCharge();
	  done=_current->createMode(Wcharge,tcurrmap[ix],mode,2,1,channel,min);
	  if(!done){throw InitException() << "Failed to construct mode in "
					  << "BaryonFactorizedDecayer::doinit()." 
					  << Exception::abortnow;}
	  // set the parameters for the additional modes
	  ttform.clear();ttcurr.clear();
	  ttform.push_back(tformmap[ix]);ttcurr.push_back(tcurrmap[ix]);
	  for(iy=0;iy<modeloc.size();++iy)
	    {
	      ttform.push_back(tformmap[modeloc[iy]]);
	      ttcurr.push_back(tcurrmap[modeloc[iy]]);
	    }
	  tCKM.clear();
	  for(iy=0;iy<ttcurr.size();++iy)
	    {
	      // get the quarks involved in the process
	      if(iy==0)
		{iq=currq[ix];ia=curra[ix];inq=inquark[ix];outq=outquark[ix];}
	      else
		{
		  if(!modecc[iy-1])
		    {
		      iq=currq[modeloc[iy-1]];ia=curra[modeloc[iy-1]];
		      inq=inquark[modeloc[iy-1]];outq=outquark[modeloc[iy-1]];
		    }
		  else
		    {
		      ia=-currq[modeloc[iy-1]];iq=-curra[modeloc[iy-1]];
		      inq=-inquark[modeloc[iy-1]];outq=-outquark[modeloc[iy-1]];
		    }
		}
	      _form->particleID(ttform[iy],id0,id1);
	      Wcharge = getParticleData(id0)->iCharge()-getParticleData(id1)->iCharge();
	      ckm=1.;
	      if(Wcharge!=0)
		{
		  ckm=1.;
		  if(abs(iq)%2==0){ckm *= conj(ckmmat[abs(iq)/2-1][(abs(ia)-1)/2]);}
		  else{ckm *= conj(ckmmat[abs(ia)/2-1][(abs(iq)-1)/2]);}
		  if(abs(inq)%2==0){ckm *= ckmmat[abs(inq)/2-1][(abs(outq)-1)/2];}
		  else{ckm *= ckmmat[abs(outq)/2-1][(abs(inq)-1)/2];}
		  if(abs(inq)==5){ckm*=_a1b;}
		  else{ckm*=_a1c;}
		}
	      else
		{
		  ckm=1.;
		  if(inq>0)
		    {
		      if(abs(inq)%2==0){ckm *= ckmmat[abs(inq)/2-1][(abs(iq)-1)/2];}
		      else{ckm *= ckmmat[abs(iq)/2-1][(abs(inq)-1)/2];}
		      if(abs(outq)%2==0)
			{ckm *= conj(ckmmat[abs(outq)/2-1][(abs(ia)-1)/2]);}
		      else{ckm *= conj(ckmmat[abs(ia)/2-1][(abs(outq)-1)/2]);}
		    }
		  else
		    {
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
	  _currentmap.push_back(ttcurr);_formmap.push_back(ttform);
	  _factCKM.push_back(tCKM);
	  // add the mode to the list
	  if(_wgtmax.size()>numberModes()){maxweight=_wgtmax[numberModes()];}
	  else{maxweight=0.;}
	  // the weights for the channel
	  if(_wgtloc.size()>numberModes()&&
	     _wgtloc[numberModes()]+mode->numberChannels()<=_weights.size())
	    {
	      start=_weights.begin()+_wgtloc[numberModes()];
	      end  = start+mode->numberChannels();
	      channelwgts=vector<double>(start,end);
	    }
	  else
	    {channelwgts.resize(mode->numberChannels(),1./(mode->numberChannels()));}
	  // don't need channels for two body decays
	  if(particles[ix].size()==3)
	    {
	      channelwgts.clear();
	      mode=new_ptr(DecayPhaseSpaceMode(particles[ix],this));
	    }
	  addMode(mode,maxweight,channelwgts);
	  // resize the duplicate modes to remove them
	  for(iy=0;iy<modeloc.size();++iy){particles[modeloc[iy]]=tPDVector();}
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
  do
    {
      _form->particleID(iform,id0,id1);
      ibaryon=0;
      if(id0==idin){ibaryon=id1;}
      else if(id0==-idin){ibaryon=-id1;}
      if(ibaryon!=0)
	{
	  foundb=false;
	  idother.clear();
	  for(ix=0;ix<idall.size();++ix)
	    {
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
  if(imode<0){throw DecayIntegratorError() << "Unable to find the mode in " 
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

ClassDescription<BaryonFactorizedDecayer> BaryonFactorizedDecayer::initBaryonFactorizedDecayer;
// Definition of the static class description member.

void BaryonFactorizedDecayer::Init() {

  static ClassDocumentation<BaryonFactorizedDecayer> documentation
    ("The BaryonFactorizedDecayer class combines the baryon form factor and a"
     " weak current to perform a decay in the naive factorization approximation.");

  static Reference<BaryonFactorizedDecayer,WeakDecayCurrent> interfaceWeakCurrent
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

double BaryonFactorizedDecayer::me2(const int ichan,
				    const Particle & inpart,
				    const ParticleVector & decay,
				    MEOption meopt) const {
  double me(0.);
  assert(inpart.dataPtr()->iSpin()==2);
  if(decay[0]->dataPtr()->iSpin()==2)
    me=halfHalf(ichan,inpart,decay,meopt);
  else if(decay[0]->dataPtr()->iSpin()==4)
    me=halfThreeHalf(ichan,inpart,decay,meopt);
  else
    assert(false);
  return me;
}

// matrix element for a 1/2 -> 1/2 decay
double BaryonFactorizedDecayer::halfHalf(const int ichan,
					 const Particle & inpart,
					 const ParticleVector & decay,
					 MEOption meopt) const {
  Energy scale;
  // extract the spins of the particles
  vector<PDT::Spin> spin;
  for(unsigned ix=0;ix<decay.size();++ix) 
    spin.push_back(decay[ix]->dataPtr()->iSpin());
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half,spin)));
  if(meopt==Initialize) {
    // spinors and rho
    if(inpart.id()>0)
      SpinorWaveFunction   ::calculateWaveFunctions(_inHalf,_rho,
						    const_ptr_cast<tPPtr>(&inpart),
						    incoming);
    else
      SpinorBarWaveFunction::calculateWaveFunctions(_inHalfBar,_rho,
						    const_ptr_cast<tPPtr>(&inpart),
						    incoming);
  }
  // setup spin info when needed
  if(meopt==Terminate) {
    // for the decaying particle
    if(inpart.id()>0) {
      SpinorWaveFunction::
	constructSpinInfo(_inHalf,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      SpinorBarWaveFunction::constructSpinInfo(_inHalfBar,decay[0],outgoing,true);
    }
    else {
      SpinorBarWaveFunction::
	constructSpinInfo(_inHalfBar,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      SpinorWaveFunction::constructSpinInfo(_inHalf,decay[0],outgoing,true);
    }
    ParticleVector::const_iterator start = decay.begin()+1,
                                   end   = decay.end();
    ParticleVector hadpart(start,end);
    _current->current(_currentmap[imode()][0],
		      ichan,scale,hadpart,meopt);
    return 0.;
  }
  ME()->zero();
  // spinors for the decay product
  if(inpart.id()>0) {
    SpinorBarWaveFunction::calculateWaveFunctions(_inHalfBar,decay[0],outgoing);
  }
  else {
    SpinorWaveFunction   ::calculateWaveFunctions(_inHalf,decay[0],outgoing);
  }
  // get the information on the form-factor
  int id0(inpart.id()),id1(decay[0]->id());
  // work out the value of q and calculate the form factors
  Lorentz5Momentum q(inpart.momentum()-decay[0]->momentum());q.rescaleMass();
  Energy m0(inpart.mass()),m1(decay[0]->mass());
  Energy2 q2(q.mass2());
  Lorentz5Momentum sum(inpart.momentum()+decay[0]->momentum());
  // calculate the baryon part of the current for the decay
  vector<LorentzPolarizationVectorE> baryon;
  Complex f1v,f2v,f3v,f1a,f2a,f3a;
  baryon.resize(4);
  Complex left,right;
  ParticleVector::const_iterator start,end;
  ParticleVector hadpart;
  vector<LorentzPolarizationVectorE> hadron;
  double pre(0.);
  unsigned int mhel,ix,iy,lhel;
  vector<unsigned int> constants,ihel;
  int itemp; unsigned int ibar;
  for(unsigned int mode=0;mode<_formmap[imode()].size();++mode) {
    // calculate the form factor piece
    _form->SpinHalfSpinHalfFormFactor(q2,_formmap[imode()][mode],id0,id1,m0,m1,
				      f1v,f2v,f3v,f1a,f2a,f3a);
    left =  f1v-f1a-f2v-double((m0-m1)/(m0+m1))*f2a;
    right = f1v+f1a-f2v+double((m0-m1)/(m0+m1))*f2a;
    for(ix=0;ix<2;++ix) {
      for(iy=0;iy<2;++iy) {
	LorentzPolarizationVectorE 
	  vtemp = _inHalf[ix].generalCurrent(_inHalfBar[iy],left,right);
	complex<Energy> vspin=_inHalf[ix].scalar(_inHalfBar[iy]);      
	complex<Energy> aspin=_inHalf[ix].pseudoScalar(_inHalfBar[iy]);
	// the momentum like pieces
	if(inpart.id()>0) {
	  vtemp+= (f2v*vspin+f2a*aspin)/(m0+m1)*sum;
	  vtemp+= (f3v*vspin+f3a*aspin)/(m0+m1)*q;
	}
	else {
	  vtemp+= (f2v*vspin-f2a*aspin)/(m0+m1)*sum;
	  vtemp+= (f3v*vspin-f3a*aspin)/(m0+m1)*q;
	}
	if(inpart.id()>0){baryon[2*ix+iy]=vtemp;}
	else{baryon[2*iy+ix]=vtemp;}
      }
    }
    // construct the weak current
    start=decay.begin()+1;
    end  =decay.end();
    hadpart = ParticleVector(start,end);
    hadron=_current->current(_currentmap[imode()][mode],
			     ichan,scale,hadpart,meopt);
    pre=pow(inpart.mass()/scale,int(hadpart.size()-2));pre*=pre;
    constants.resize(decay.size()+1);ihel.resize(decay.size()+1);
    itemp=1;ibar=0;
    for(int iz=int(decay.size()-1);iz>=0;--iz) {
      if(abs(decay[iz]->id())!=id1) {
	itemp *= decay[iz]->data().iSpin();
	constants[iz]=itemp;
      }
      else ibar=iz;
      constants[decay.size()]=1;
      constants[ibar]=constants[ibar+1];
    }
    for(mhel=0;mhel<baryon.size();++mhel) {
      ihel[0     ]=mhel/2;
      ihel[ibar+1]=mhel%2;
      for(lhel=0;lhel<hadron.size();++lhel) {
	// map the index for the hadrons to a helicity state
	for(ix=decay.size();ix>0;--ix) {
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
double BaryonFactorizedDecayer::halfThreeHalf(const int ichan,
					      const Particle & inpart,
					      const ParticleVector & decay,
					      MEOption meopt) const {
  // spins
  Energy scale;
  vector<PDT::Spin> spin(decay.size());
  for(unsigned int ix=0;ix<decay.size();++ix)
    spin[ix]=decay[ix]->data().iSpin();
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half,spin)));
  // spinors etc for the decaying particle
  if(meopt==Initialize) {
    // spinors and rho
    if(inpart.id()>0)
      SpinorWaveFunction   ::calculateWaveFunctions(_inHalf,_rho,
						    const_ptr_cast<tPPtr>(&inpart),
						    incoming);
    else
      SpinorBarWaveFunction::calculateWaveFunctions(_inHalfBar,_rho,
						    const_ptr_cast<tPPtr>(&inpart),
						    incoming);
  }
  // setup spin info when needed
  if(meopt==Terminate) {
    // for the decaying particle
    if(inpart.id()>0) {
      SpinorWaveFunction::
	constructSpinInfo(_inHalf,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      RSSpinorBarWaveFunction::constructSpinInfo(_inThreeHalfBar,
						 decay[0],outgoing,true);

    }
    else {
      SpinorBarWaveFunction::
	constructSpinInfo(_inHalfBar,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      RSSpinorWaveFunction::constructSpinInfo(_inThreeHalf,
					      decay[0],outgoing,true);
    }
    ParticleVector::const_iterator start = decay.begin()+1,
                                   end   = decay.end();
    ParticleVector hadpart(start,end);
    _current->current(_currentmap[imode()][0],
		      ichan,scale,hadpart,meopt);
    return 0.;
  }
  ME()->zero();
  // spinors for the decay product
  LorentzPolarizationVector in=UnitRemoval::InvE*inpart.momentum();
  if(inpart.id()>0) {
    RSSpinorBarWaveFunction::
      calculateWaveFunctions(_inThreeHalfBar,decay[0],outgoing);
    _inHalfBar.resize(_inThreeHalfBar.size());
    for(unsigned int ix=0;ix<_inThreeHalfBar.size();++ix)
      _inHalfBar[ix] = _inThreeHalfBar[ix].dot(in);
    }
  else {
    RSSpinorWaveFunction::
      calculateWaveFunctions(_inThreeHalf,decay[0],outgoing);
    _inHalf.resize(_inThreeHalf.size());
    for(unsigned int ix=0;ix<_inThreeHalf.size();++ix)
      _inHalf[ix] = _inThreeHalf[ix].dot(in);
  }
  // get the information on the form-factor
  int id0(inpart.id()),id1(decay[0]->id());
  // work out the value of q and calculate the form factors
  Lorentz5Momentum q(inpart.momentum()-decay[0]->momentum());q.rescaleMass();
  Energy m0(inpart.mass()),m1(decay[0]->mass());
  Energy2 q2(q.mass2());
  Lorentz5Momentum sum(inpart.momentum()+decay[0]->momentum());
  LorentzPolarizationVectorE baryon[4][2];
  Complex f1v,f2v,f3v,f4v,f1a,f2a,f3a,f4a;
  LorentzPolarizationVector vtemp;
  complex<InvEnergy2> lS1,lS2,rS1,rS2;
  Complex left,right;
  complex<InvEnergy> lV,rV;
  InvEnergy ms(1./(m0+m1));
  InvEnergy2 ms2(ms*ms);
  unsigned int ix,iy,ixa,iya;
  ParticleVector::const_iterator start,end;
  ParticleVector hadpart;
  vector<LorentzPolarizationVectorE> hadron;
  double pre(0.);
  vector<unsigned int> constants,ihel;
  int itemp; unsigned int ibar;
  for(unsigned int mode=0;mode<_formmap[imode()].size();++mode) {
    // calculate the form factors
    _form->SpinHalfSpinThreeHalfFormFactor(q2,_formmap[imode()][mode],id0,id1,m0,m1,
					   f1v,f2v,f3v,f4v,f1a,f2a,f3a,f4a);
    if(inpart.id()>0) {
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
    Complex scalar1,scalar2;
    complex<Energy> lfact,rfact;
    LorentzPolarizationVectorE tvec;
    LorentzPolarizationVector  svec;
    for(iya=0;iya<4;++iya) {
      for(ixa=0;ixa<2;++ixa) {
	if(decay[0]->id()>0){ix=iya;iy=ixa;}
	else{ix=ixa;iy=iya;}
	// scalar like terms
	lfact = _inHalf[iy].leftScalar( _inHalfBar[ix]);
	rfact = _inHalf[iy].rightScalar(_inHalfBar[ix]);
	scalar1 = (lS1*lfact+rS1*rfact)*UnitRemoval::E;
	scalar2 = (lS2*lfact+rS2*rfact)*UnitRemoval::E;
	svec = _inHalf[iy].generalCurrent(_inHalfBar[ix],lV/ms,rV/ms)*ms;
	if(inpart.id()>0) {
	  tvec=_inThreeHalfBar[ix].generalCurrent(_inHalf[iy],left,right);
	}
	else {
	  tvec=_inThreeHalf[iy].generalCurrent(_inHalfBar[ix],left,right);
	}
	baryon[iya][ixa] = tvec+svec*UnitRemoval::E
	  +scalar1*decay[0]->momentum()+scalar2*inpart.momentum();
      }
    }
    start=decay.begin()+1;
    end  =decay.end();
    hadpart=ParticleVector(start,end);
    hadron=_current->current(_currentmap[imode()][mode],
			     ichan,scale,hadpart,meopt);
    // prefactor
    pre = pow(inpart.mass()/scale,int(hadpart.size()-2));pre*=pre;
    // work out the mapping for the hadron vector
    constants.resize(decay.size()+1);ihel.resize(decay.size()+1);
    itemp=1;ibar=0;
    for(int ix=int(decay.size()-1);ix>=0;--ix) {
      if(abs(decay[ix]->id())!=id1) {
	itemp*=decay[ix]->data().iSpin();
	constants[ix]=itemp;
      }
      else{ibar=ix;}
    }
    constants[decay.size()]=1;
    constants[ibar]=constants[ibar+1];
    for(iya=0;iya<4;++iya) {
      ihel[1]=iya;
      for(ixa=0;ixa<2;++ixa) {
	ihel[0]=ixa;
	for(unsigned int lhel=0;lhel<hadron.size();++lhel) {
	  // map the index for the hadrons to a helicity state
	  for(unsigned int ix=decay.size();ix>0;--ix)
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
					vector<tPDVector> & particles,
					vector<unsigned int> & loc,
					vector<bool> & cc) {
  unsigned int ix,iy,nfound,iz;
  // resize the vectors
  loc.clear();cc.clear();
  // get the id's for the mode
  vector<int> id,idbar;
  int idtemp; bool found;
  for(ix=0;ix<particles[imode].size();++ix)
    {
      id.push_back(particles[imode][ix]->id());
      if(particles[imode][ix]->CC()){idbar.push_back(particles[imode][ix]->CC()->id());}
      else{idbar.push_back(id[ix]);}
    }
  vector<bool> done(id.size(),false);
  // loop over the modes
  for(ix=0;ix<particles.size();++ix)
    {
      if(ix==imode||particles[ix].empty()) continue;
	  assert(!particles[ix].empty());
	  assert(particles[ix][0]);
	  // the particle mode
	  if(particles[ix][0]->id()==id[0]&&particles[ix].size()==id.size())
	    {
	      nfound=1;
	      for(iy=0;iy<id.size();++iy){done[iy]=false;}
	      for(iy=1;iy<id.size();++iy)
		{
		  idtemp=particles[ix][iy]->id();
		  iz=1;found=false;
		  do{if(idtemp==id[iz]&&!done[iz]){done[iz]=true;found=true;}++iz;}
		  while(iz<id.size()&&!found);
		  if(found){++nfound;}
		}
	      if(nfound==id.size()){cc.push_back(false);loc.push_back(ix);}
	    }
	  // the charge conjugate mode
	  if(particles[ix][0]->id()==idbar[0]&&particles[ix].size()==idbar.size())
	    {
	      nfound=1;
	      for(iy=0;iy<idbar.size();++iy){done[iy]=false;}
	      for(iy=1;iy<idbar.size();++iy)
		{
		  idtemp=particles[ix][iy]->id();
		  iz=1;found=false;
		  do{if(idtemp==idbar[iz]&&!done[iz]){done[iz]=true;found=true;}++iz;}
		  while(iz<idbar.size()&&!found);
		  if(found){++nfound;}
		}
	      if(nfound==idbar.size()){cc.push_back(false);loc.push_back(ix);}
	    }
    }
}

// output the setup information for the particle database
void BaryonFactorizedDecayer::dataBaseOutput(ofstream & output, bool header) const
{
  unsigned int ix;
  if(header){output << "update decayers set parameters=\"";}
  DecayIntegrator::dataBaseOutput(output,false);
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
