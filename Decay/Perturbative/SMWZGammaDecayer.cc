// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMWZGammaDecayer class.
//

#include "SMWZGammaDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/Interface/ParVector.h"
#include "Herwig++/Utilities/Kinematics.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"

using namespace Herwig;
using namespace ThePEG::Helicity;
using namespace ThePEG::Helicity;


void SMWZGammaDecayer::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  // get the vertices from the Standard Model object
  tcHwSMPtr hwsm=dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if(hwsm) {
    _wvertex = hwsm->vertexFFW();
    _zvertex = hwsm->vertexFFZ();
    _pvertex = hwsm->vertexFFP();
    _wwwvertex = hwsm->vertexWWW();
    // make sure they are initialized
    _wvertex->init();
    _zvertex->init();
    _pvertex->init();
    _wwwvertex->init();
  }
  else {
    throw InitException();
  }
  // now set up the decay modes
  DecayPhaseSpaceModePtr mode;
  PDVector extpart(4);
  // the Z decay modes
  extpart[0]=getParticleData(ParticleID::Z0);
  extpart[3]=getParticleData(ParticleID::gamma);
  // loop over the  quarks and the leptons
  DecayPhaseSpaceChannelPtr newchannel;
  unsigned int ix,istep=0,iy,iloc;
  vector<double> wgt(2);
  for( ;istep<11;istep+=10)  {
    iloc=0;
    for(ix=1;ix<7;++ix) {
      iy=istep+ix;
      if(iy!=6&&!(istep==10&&ix%2==0)) {
	// check that the combination of particles is allowed
	if(_zvertex->allowed(-iy,iy,ParticleID::Z0)) {
	  extpart[1] = getParticleData(-iy);
	  extpart[2] = getParticleData( iy);
	  mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
	  newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
	  newchannel->addIntermediate(extpart[0],0, 0.0,-1,1);
	  newchannel->addIntermediate(extpart[2],1,-0.9, 2,3);
	  mode->addChannel(newchannel);
	  newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
	  newchannel->addIntermediate(extpart[0],0, 0.0,-1,2);
	  newchannel->addIntermediate(extpart[1],1,-0.9, 1,3);
	  mode->addChannel(newchannel);
	  if(iy<=6) {
	    if(iloc+2<=_Zquarkchannels.size()) {
	      wgt[0]=_Zquarkchannels[iloc];
	      wgt[1]=_Zquarkchannels[iloc+1];
	      iloc+=2;
	    }
	    else {
	      wgt[0]=0.5;
	      wgt[1]=0.5;
	      iloc+=2;
	    }
	    addMode(mode,_Zquarkwgt[ix-1],wgt);
	  }
	  else {
	    if(iloc+2<=_Zleptonchannels.size()) {
	      wgt[0]=_Zleptonchannels[iloc];
	      wgt[1]=_Zleptonchannels[iloc+1];
	      iloc+=2;
	    }
	    else {
	      wgt[0]=0.5;
	      wgt[1]=0.5;
	      iloc+=2;
	    }
	    addMode(mode,_Zleptonwgt[ix-11],wgt);
	  }
	}
	else {
	  throw InitException() << "SMWZGammaDecayer::doinit() the Z vertex" 
				<< "cannot handle all the modes" 
				<< Exception::abortnow;
	}
      }
    }
  }
  // and the W modes
  extpart[0]=getParticleData(ParticleID::Wplus);
  // loop for the quarks
  unsigned int iz=0;
  for(ix=1;ix<6;ix+=2) {
    for(iy=2;iy<6;iy+=2) {
      // check that the combination of particles is allowed
      if(_wvertex->allowed(-ix,iy,ParticleID::Wminus)) {
	extpart[1] = getParticleData(-ix);
	extpart[2] = getParticleData( iy);
	mode = new DecayPhaseSpaceMode(extpart,this);
	addMode(mode,_Wquarkwgt[iz],wgt);
	++iz;
      }
      else {
	throw InitException() << "SMWZDecayer::doinit() the W vertex" 
			      << "cannot handle all the quark modes" 
			      << Exception::abortnow;
      }
    }
  }
  for(ix=11;ix<17;ix+=2) {
    // check that the combination of particles is allowed
    if(_wvertex->allowed(-ix,ix+1,ParticleID::Wminus)) {
      extpart[1] = getParticleData(-ix);
      extpart[2] = getParticleData(ix+1);
      mode = new DecayPhaseSpaceMode(extpart,this);
      addMode(mode,_Wleptonwgt[(ix-11)/2],wgt);
    }
    else {
      throw InitException() << "SMWZDecayer::doinit() the W vertex" 
			    << "cannot handle all the lepton modes" 
			    << Exception::abortnow;
    }
  }
}

void SMWZGammaDecayer::persistentOutput(PersistentOStream & os) const {
  os << _wvertex << _zvertex << _pvertex << _wwwvertex  
     << _Zquarkwgt << _Wquarkwgt << _Zleptonwgt << _Wleptonwgt
     << _Zquarkchannels << _Wquarkchannels << _Zleptonchannels << _Wleptonchannels;
}

void SMWZGammaDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _wvertex >> _zvertex >> _pvertex >> _wwwvertex
     >> _Zquarkwgt >> _Wquarkwgt >> _Zleptonwgt >> _Wleptonwgt
     >> _Zquarkchannels >> _Wquarkchannels >> _Zleptonchannels >> _Wleptonchannels;
}

ClassDescription<SMWZGammaDecayer> SMWZGammaDecayer::initSMWZGammaDecayer;
// Definition of the static class description member.

void SMWZGammaDecayer::Init() {

  static ClassDocumentation<SMWZGammaDecayer> documentation
    ("There is no documentation for the SMWZGammaDecayer class");

  static ParVector<SMWZGammaDecayer,double> interfaceZquarkMax
    ("ZquarkMax",
     "The maximum weight for the decay of the Z to quarks",
     &SMWZGammaDecayer::_Zquarkwgt,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<SMWZGammaDecayer,double> interfaceWquarkMax
    ("WquarkMax",
     "The maximum weight for the decay of the W to quarks",
     &SMWZGammaDecayer::_Wquarkwgt,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<SMWZGammaDecayer,double> interfaceZleptonMax
    ("ZleptonMax",
     "The maximum weight for the decay of the Z to leptons",
     &SMWZGammaDecayer::_Zleptonwgt,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<SMWZGammaDecayer,double> interfaceWleptonMax
    ("WleptonMax",
     "The maximum weight for the decay of the W to leptons",
     &SMWZGammaDecayer::_Wleptonwgt,
     0, 0, 0, -10000, 10000, false, false, true);

}

int SMWZGammaDecayer::modeNumber(bool & cc,tcPDPtr parent, 
				 const PDVector & children) const {
  cc = false;
  int imode(-1);
  int id0(parent->id());
  if(children.size()!=3){return imode;}
  PDVector::const_iterator pit = children.begin();
  int id1=(**pit).id();
  ++pit;
  int id2=(**pit).id();
  ++pit;
  // check last particle a photon
  if((**pit).id()!=ParticleID::gamma){return imode;}
  // Z to quarks or leptons
  if(id0==ParticleID::Z0) {
    if(id1==-id2&&(abs(id1)<=5||(abs(id1)>=11&&abs(id1)<=16))) {
      if(abs(id1)<6) {
	imode=abs(id1)-1;
      }
      else if(abs(id1)>=11&&abs(id1)<=16) {
	imode=(abs(id1)-10)/2+5;
      }
    }
  }
  return imode;
}

double SMWZGammaDecayer::me2(bool vertex, const int ichan, const Particle & inpart,
			     const ParticleVector & decay) const {
  if(decay[2]->momentum().e()<_emin) return 0.;
  // check if the incoming particle has a spin info 
  tcVectorSpinPtr inspin;
  if(inpart.spinInfo()) {
    inspin = dynamic_ptr_cast<tcVectorSpinPtr>(inpart.spinInfo());
  }
  RhoDMatrix rhoin(PDT::Spin1);rhoin.average();
  VectorWaveFunction inwave[3];
  // if the spin info object exists use it
  if(inspin&&inpart.spinInfo()) {
    for(unsigned int ix=0;ix<3;++ix) {
      inwave[ix]=VectorWaveFunction(inpart.momentum(),inpart.dataPtr(),
				    inspin->getDecayBasisState(ix),incoming);
    }
    inspin->decay();
    rhoin = inspin->rhoMatrix();
  }
  else {
    // if has spin info but not the right type issue warning and throw away
    if(inpart.spinInfo()) {
      throw DecayIntegratorError() << "Wrong type of spin info for the "
				   << "incoming particle in SMWZDecayer::me2()" 
				   << Exception::warning;}
    SpinPtr newspin=new_ptr(VectorSpinInfo(inpart.momentum(),true));
    inspin = dynamic_ptr_cast<tcVectorSpinPtr>(newspin);
    inspin->decayed(true);
    VectorWaveFunction temp=VectorWaveFunction(inpart.momentum(),inpart.dataPtr(),
					       incoming);
    for(int ix=0;ix<3;++ix) {
      temp.reset(ix);
      inwave[ix]=temp;
      inspin->setDecayState(ix,temp.wave());
    }
    const_ptr_cast<tPPtr>(&inpart)->spinInfo(newspin);
  }
  // locations of the fermions
  int iferm,ianti;
  if(decay[0]->id()<0) {
    iferm=1;
    ianti=0;
  }
  else {
    iferm=0;
    ianti=1;
  }
  SpinorWaveFunction awave;
  SpinorBarWaveFunction fwave;
  VectorWaveFunction vwave;  // compute the matrix element
  DecayMatrixElement newme(PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1);
  Energy2 scale(inpart.mass()*inpart.mass());
  Complex me[2];
  // the full matrix element
  if(_iopt==0) {
    // wavefunctions
    awave=SpinorWaveFunction(   decay[ianti]->momentum(),decay[ianti]->dataPtr(),
				outgoing);
    fwave=SpinorBarWaveFunction(decay[iferm]->momentum(),decay[iferm]->dataPtr(),
				outgoing);
    vwave=VectorWaveFunction(   decay[2]->momentum()    ,decay[2]->dataPtr(),
				outgoing);
    unsigned int ifm,ia,vhel,phel;
    SpinorWaveFunction atemp;
    SpinorBarWaveFunction ftemp;
    for(ifm=0;ifm<2;++ifm) {
      fwave.reset(ifm);
      for(ia=0;ia<2;++ia) {
	awave.reset(ia);
	for(phel=0;phel<3;phel+=2) {
	  vwave.reset(phel);
	  // calculate the combinations of photon and (anti)fermion
	  atemp=_pvertex->evaluate(0.*GeV2,3,decay[ianti]->dataPtr(),awave,vwave);
	  ftemp=_pvertex->evaluate(0.*GeV2,3,decay[iferm]->dataPtr(),fwave,vwave);
	  // calculate the matrix element
	  for(vhel=0;vhel<3;++vhel) {
	    if(inpart.id()==ParticleID::Z0) {
	      // the matrix elements
	      me[0]=_zvertex->evaluate(scale,awave,ftemp,inwave[vhel]);
	      me[1]=_zvertex->evaluate(scale,atemp,fwave,inwave[vhel]);
	      if(ichan<0){me[0]+=me[1];}
	      else{me[0]=me[ichan];}
	      if(iferm>ianti){newme(vhel,ia,ifm,phel)=me[0];}
	      else           {newme(vhel,ifm,ia,phel)=me[0];}
	    }
	    else {
	      throw Exception() << "W not yet implemented in SMWZGammaDecayer"
				<< "::me2() " 
				<< Exception::abortnow;
	    }
	  }
	}
      }
    }
  }
  // dipole approximation
  else {
    Energy mout=decay[0]->mass();
    Energy pcm(Kinematics::pstarTwoBodyDecay(inpart.mass(),mout,mout));
    double costh,phi;
    Kinematics::generateAngles(costh,phi);
    double sinth(sqrt(1.-costh*costh));
    Lorentz5Momentum p1(+sinth*cos(phi)*pcm,+sinth*sin(phi)*pcm,+costh*pcm,
			0.5*inpart.mass(),mout);
    Lorentz5Momentum p2(-sinth*cos(phi)*pcm,-sinth*sin(phi)*pcm,-costh*pcm,
			0.5*inpart.mass(),mout);
    // wavefunctions
    awave=SpinorWaveFunction(   p1,decay[ianti]->dataPtr(),outgoing);
    fwave=SpinorBarWaveFunction(p2,decay[iferm]->dataPtr(),outgoing);
    vwave=VectorWaveFunction(decay[2]->momentum(),decay[2]->dataPtr(),outgoing);
    complex<Energy> me0[3][2][2];
    complex<InvEnergy> dipole[3][2];
    // leading order matrix element
    unsigned int ifm,ia,vhel,phel;
    for(ifm=0;ifm<2;++ifm) {
      fwave.reset(ifm);
      for(ia=0;ia<2;++ia) {
	awave.reset(ia);
	for(vhel=0;vhel<3;++vhel) {
	  me0[ifm][ia][vhel]=_zvertex->evaluate(scale,awave,fwave,inwave[vhel])
	    * UnitRemoval::E;}
      }
    }
    // dipole term
    Energy2 
      p1dot(decay[2]->momentum()*decay[0]->momentum()),
      p2dot(decay[2]->momentum()*decay[1]->momentum());
    complex<Energy> eps1,eps2;
    double fact(sqrt(4.*Constants::pi*standardModel()->alphaEM())*
		decay[0]->dataPtr()->iCharge()/3.);
    for(phel=0;phel<3;phel+=2) {
      vwave.reset(phel);
      eps1=vwave.wave()*decay[0]->momentum();
      eps2=vwave.wave()*decay[1]->momentum();
      dipole[phel][0]= eps1/p1dot;
      dipole[phel][1]=-eps2/p2dot;
    }
    // compute the matrix element
    for(ifm=0;ifm<2;++ifm){
      for(ia=0;ia<2;++ia) {
	for(phel=0;phel<3;phel+=2) {
	  for(vhel=0;vhel<3;++vhel) {
	    if(inpart.id()==ParticleID::Z0) {
	      // the matrix elements
	      me[0]=fact*dipole[phel][0]*me0[ifm][ia][vhel];
	      me[1]=fact*dipole[phel][1]*me0[ifm][ia][vhel];
	      if(ichan<0){me[0]+=me[1];}
	      else{me[0]=me[ichan];}
	      if(iferm>ianti){newme(vhel,ia,ifm,phel)=me[0];}
	      else           {newme(vhel,ifm,ia,phel)=me[0];}
	    }
	    else {
	      throw Exception() << "W not yet implemented in SMWZGammaDecayer"
				<< "::me2() " 
				<< Exception::abortnow;
	    }
	  }
	}
      } 
    }
  }
  ME(newme);
  double output((newme.contract(rhoin)).real());
  // spin info for the outgoing particles
  if(vertex) {
    FermionSpinPtr ferm,anti;
    VectorSpinPtr vect;
    ferm = new_ptr(FermionSpinInfo(decay[iferm]->momentum(),true));
    anti = new_ptr(FermionSpinInfo(decay[ianti]->momentum(),true));
    vect = new_ptr( VectorSpinInfo(decay[2]->momentum()    ,true));
    vect->setBasisState(1,LorentzPolarizationVector());
    decay[iferm]->spinInfo(ferm);
    decay[ianti]->spinInfo(anti);
    decay[2]->spinInfo(vect);
    awave.reset(decay[ianti]->momentum());
    fwave.reset(decay[iferm]->momentum());
    for(unsigned int ix=0;ix<2;++ix) {
      awave.reset(ix);
      fwave.reset(ix);
      vwave.reset(2*ix);
      ferm->setBasisState(ix  ,fwave.dimensionedWave().bar());
      anti->setBasisState(ix  ,awave.dimensionedWave());
      vect->setBasisState(2*ix,vwave.wave());
    }
  }
  // final colour factors and connections if needed
  if(abs(decay[0]->id())<=6) output*=3.;
  if(decay[0]->hasColour()) {
    decay[0]->antiColourNeighbour(decay[1]);
  }
  else if(decay[1]->hasColour()) {
    decay[1]->antiColourNeighbour(decay[0]);
  }
  return output;
}
