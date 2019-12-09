// -*- C++ -*-
//
// WeakPartonicDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the WeakPartonicDecayer class.
//

#include "WeakPartonicDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/ConstituentParticleData.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

WeakPartonicDecayer::WeakPartonicDecayer() : MECode(0), _radprob(0.0), _maxtry(300), 
					     _threemax(3.), _fourmax(3.) 
{}

IBPtr WeakPartonicDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr WeakPartonicDecayer::fullclone() const {
  return new_ptr(*this);
}

bool WeakPartonicDecayer::accept(tcPDPtr parent, const tPDVector & prod) const {
  // check we can find the flavours of the quarks in the decaying meson
  long id = parent->id();
  int flav1, flav2;
  if((id / 1000)%10) {
    flav1 = (id/1000)%10;
    flav2 = (id/10)%100;
  } 
  else {
    flav1 = id/100;  
    flav2 = (id/10)%10;
  }
  if(!flav1 || !flav2) return false;
  // if two decay products one must be in triplet and one antitriplet
  if(prod.size()==2) {
    if((prod[0]->iColour()==PDT::Colour3&&prod[1]->iColour()==PDT::Colour3bar)||
       (prod[0]->iColour()==PDT::Colour3bar&&prod[1]->iColour()==PDT::Colour3))
      return true;
  }
  else if(prod.size()==3) {
    if(((prod[0]->iColour()==PDT::Colour3   &&prod[2]->iColour()==PDT::Colour3bar)||
	(prod[0]->iColour()==PDT::Colour3bar&&prod[2]->iColour()==PDT::Colour3))
       &&prod[1]->iColour()==PDT::Colour8) return true;
  }
  else if(prod.size()==4) {
    // first two particles should be leptons or q qbar
    if((prod[0]->id()>=11&&prod[0]->id()<=16&&prod[1]->id()<=-11&&prod[1]->id()>=-16)||
       (prod[1]->id()>=11&&prod[1]->id()<=16&&prod[0]->id()<=-11&&prod[0]->id()>=-16)||
       (prod[0]->iColour()==PDT::Colour3    &&prod[1]->iColour()==PDT::Colour3bar   )||
       (prod[1]->iColour()==PDT::Colour3    &&prod[0]->iColour()==PDT::Colour3bar   )) {
      // third particle quark and fourth colour anti-triplet or
      // thrid particle antiquark and fourth colour triplet
      if((prod[2]->iColour()==PDT::Colour3bar&&prod[3]->iColour()==PDT::Colour3   )||
	 (prod[2]->iColour()==PDT::Colour3   &&prod[3]->iColour()==PDT::Colour3bar))
	return true;
    }
    else return false;
  }
  return false;
}

ParticleVector WeakPartonicDecayer::decay(const Particle & parent,
					  const tPDVector & children) const {
  static tcPDPtr gluon=getParticleData(ParticleID::g);
  // make the particles
  ParticleVector partons;
  for(unsigned int ix=0;ix<children.size();++ix) {
    partons.push_back(children[ix]->produceParticle());
    // these products have the mass but should have constituent mass
    partons[ix]->set5Momentum(Lorentz5Momentum(children[ix]->constituentMass()));
  }
  // 2-body decays
  if(partons.size()==2) {
    // no gluon if not select based on probability or if three body not allowed
    if(UseRandom::rnd()>_radprob||
       parent.mass()<gluon->constituentMass()+partons[0]->mass()+partons[1]->mass()) {
      double ctheta,phi;
      Lorentz5Momentum pout[2];
      for(unsigned int ix=0;ix<2;++ix) pout[ix].setMass(partons[ix]->mass());
      Kinematics::generateAngles(ctheta,phi);
      Kinematics::twoBodyDecay(parent.momentum(),pout[0].mass(),pout[1].mass(),
			       ctheta,phi,pout[0],pout[1]);
      for(unsigned int ix=0; ix<2;++ix) partons[ix]->setMomentum(pout[ix]);
      if(partons[0]->dataPtr()->iColour()==PDT::Colour3) {
	partons[0]->antiColourNeighbour(partons[1]);
      }
      else { 
	partons[0]->    colourNeighbour(partons[1]);
      }
    }
    else {
      Lorentz5Momentum pout[3];
      for(unsigned int ix=0;ix<2;++ix) pout[ix].setMass(partons[ix]->mass());
      // add gluon
      partons.push_back(gluon->produceParticle());
      partons.back()->set5Momentum(gluon->constituentMass());
      // momentum of gluon
      pout[2] = Lorentz5Momentum(gluon->constituentMass());
      Kinematics::threeBodyDecay(parent.momentum(),pout[1],pout[0],pout[2]);
      for(unsigned int ix=0; ix<3;++ix) partons[ix]->setMomentum(pout[ix]);
      if(partons[0]->dataPtr()->iColour()==PDT::Colour3) {
	partons[0]->antiColourNeighbour(partons[2]);
	partons[1]->    colourNeighbour(partons[2]);
      }
      else { 
	partons[0]->    colourNeighbour(partons[2]);
	partons[1]->antiColourNeighbour(partons[2]);
      }
    }
  }
  // 3-body decays
  else if(partons.size()==3) {
    // set masses of products
    Lorentz5Momentum pout[3],pin(parent.momentum());
    for(unsigned int ix=0;ix<3;++ix) pout[ix].setMass(partons[ix]->mass());
    double xs(partons[2]->mass()/pin.e()),xb(1.-xs);
    pout[2]=xs*pin;
    // Get the particle quark that is decaying
    long idQ, idSpec;
    idSpec = partons[2]->id();
    idQ = (parent.id()/1000)%10;
    if(!idQ) idQ = (parent.id()/100)%10;
    // Now the odd case of a B_c where the c decays, not the b
    if(idSpec == idQ) idQ = (parent.id()/10)%10;
    // momentum of the decaying quark
    PPtr inter = getParticleData(idQ)->produceParticle(parent.momentum()*xb);
    // two body decay of heavy quark
    double ctheta,phi;
    Kinematics::generateAngles(ctheta,phi);
    Kinematics::twoBodyDecay(inter->momentum(),pout[0].mass(),pout[1].mass(),
			     ctheta,phi,pout[0],pout[1]);
    // set the momenta of the decay products
    for(unsigned int ix=0; ix<3;++ix) partons[ix]->setMomentum(pout[ix]);
    // make the colour connections
    // quark first
    if(partons[0]->data().iColour()==PDT::Colour3) {
      partons[0]->antiColourNeighbour(partons[1]);
      partons[1]->colourNeighbour(partons[0]);
      partons[1]->antiColourNeighbour(partons[2]);
      partons[2]->colourNeighbour(partons[1]);
    }
    // antiquark first
    else {
      partons[0]->colourNeighbour(partons[1]);
      partons[1]->antiColourNeighbour(partons[0]);
      partons[1]->colourNeighbour(partons[2]);
      partons[2]->antiColourNeighbour(partons[1]);
    }
  }
  // 4-body decays
  else if(partons.size()==4) {  
    // swap 0 and 1 if needed
    if(partons[1]->dataPtr()->iColour()!=PDT::Colour0&&
       partons[1]->dataPtr()->iColour()!=partons[2]->dataPtr()->iColour())
      swap(partons[0],partons[1]);
    // get the momenta of the decaying quark and the spectator
    Lorentz5Momentum pin(parent.momentum());
    double xs(partons[3]->mass()/pin.e()),xb(1.-xs);
    Lorentz5Momentum pspect(xs*pin),pdec(xb*pin);
    pspect.setMass(partons[3]->mass());
    pdec.rescaleMass();
    // Get the particle quark that is decaying
    long idQ, idSpec;
    idSpec = partons[3]->id();
    idQ = (abs(parent.id())/1000)%10;
    if(!idQ) idQ = (abs(parent.id())/100)%10;
    // Now the odd case of a B_c where the c decays, not the b
    if(abs(idSpec) == idQ) idQ = (abs(parent.id())/10)%10;
    // change sign if spectator quark or antidiquark
    if((idSpec>0&&idSpec<6)||idSpec<-6) idQ = -idQ;
    // check if W products coloured
    bool Wcol = partons[0]->coloured();
    // particle data object
    tcPDPtr dec = getParticleData(idQ);
    // momenta of the decay products
    vector<Lorentz5Momentum> pout(3,Lorentz5Momentum());
    for(unsigned int ix=0;ix<3;++ix) pout[ix].setMass(partons[ix]->mass());
    // charges of the exchanged boson and check if colour rearranged
    int c1 = dec                  ->iCharge()-partons[2]->dataPtr()->iCharge();
    int c2 = partons[0]->dataPtr()->iCharge()+partons[1]->dataPtr()->iCharge();
    bool rearranged = !(c1==c2&&abs(c1)==3);
    if(MECode==0) rearranged=false;
    if(rearranged) {
      int c3 = dec                  ->iCharge()-partons[1]->dataPtr()->iCharge();
      int c4 = partons[0]->dataPtr()->iCharge()+partons[2]->dataPtr()->iCharge();
      if(!(c3==c4&&abs(c3)==3)) {
	generator()->log() << "Unknown order for colour rearranged decay"
			   << " in WeakPartonicDecayer::decay()\n";
	generator()->log() << c1 << " " << c2 << " " << c3 << " " << c4 << "\n";
	generator()->log() << parent << "\n" << dec->PDGName() << "\n";
	for(unsigned int ix=0;ix<4;++ix) generator()->log() << *partons[ix] << "\n";
	throw Exception()  << "Unknown order for colour rearranged decay"
			   << " in WeakPartonicDecayer::decay() "
			   << Exception::runerror;
      }
      swap(pout[1]   ,pout[2]   );
      swap(partons[1],partons[2]);
    }
    // decide if three or four body using prob
    bool threeBody = UseRandom::rnd() > _radprob;
    // if four body not kinematically possible must be four body
    if(pdec.mass()<gluon->constituentMass()+pout[0].mass()+
       pout[1].mass()+pout[2].mass()) threeBody=true;
    // if code ==0 always three body
    if(MECode==0) threeBody=true;
    // three body decay
    if( threeBody ) {
      if(MECode==0) {
	Kinematics::threeBodyDecay(pdec,pout[1],pout[0],pout[2]);
      }
      else {
	// generate the kinematics
	double wgt(0.);
	Energy2 mb2max = sqr(pdec.mass()    - pout[2].mass());
	Energy2 mb2min = sqr(pout[0].mass() + pout[1].mass());
	unsigned int ntry = 0;
	do {
	  ++ntry;
	  Energy2 mb2 = (mb2max-mb2min)*UseRandom::rnd()+mb2min;
	  double CosAngle, AzmAngle;
	  // perform first decay
	  Lorentz5Momentum p01;
	  p01.setMass(sqrt(mb2));
	  Kinematics::generateAngles(CosAngle,AzmAngle);
	  Kinematics::twoBodyDecay(pdec,pout[2].mass(),p01.mass(),
				   CosAngle,AzmAngle,pout[2],p01);
	  // perform second decay
	  Kinematics::generateAngles(CosAngle,AzmAngle);
	  Kinematics::twoBodyDecay(p01,pout[0].mass(),pout[1].mass(),
				   CosAngle,AzmAngle,pout[0],pout[1]);
	  // kinematic piece of the weight
	  wgt = 
	    Kinematics::pstarTwoBodyDecay(pdec.mass(),p01    .mass(),pout[2].mass())/pdec.mass()*
	    Kinematics::pstarTwoBodyDecay(p01 .mass(),pout[0].mass(),pout[1].mass())/p01.mass();
	  // piece to improve weight variation
	  wgt *= pdec.mass()/Kinematics::pstarTwoBodyDecay(pdec.mass(),sqrt(mb2min),pout[2].mass());
	  // matrix element piece
	  wgt *= 16.*(pdec*pout[1])*(pout[0]*pout[2])/sqr(mb2max-mb2min);
	  // check doesn't violate max
	  if(wgt>_threemax) {
	    ostringstream message;
	    message << "Maximum weight for three-body decay "
		    << "violated in WeakPartonicDecayer2::decay()"
		    << "Maximum = " << _threemax << " weight = " << wgt;
	    generator()->logWarning( Exception(message.str(),Exception::warning) );
	  }
	}
	while( wgt < _threemax*UseRandom::rnd() && ntry < _maxtry );
	if(ntry==_maxtry) throw Exception() 
	  << "Too many attempts to generate three body kinematics in "
	  << "WeakPartonicDecayer2::decay()" << Exception::eventerror;
      }
      // set momenta of particles
      for(unsigned int ix=0;ix<pout.size();++ix) partons[ix]->setMomentum(pout[ix]);
      partons[3]->setMomentum(pspect);
      // special for tau leptons to get correlations
      if(abs(partons[0]->id())==ParticleID::tauminus||
	 abs(partons[1]->id())==ParticleID::tauminus)
	threeBodyMatrixElement(dec,pdec,partons);
      // set up the colour connections
      if(rearranged) swap(partons[1],partons[2]);
      if(Wcol) {
	if(partons[0]->data().iColour()==PDT::Colour3)
	  partons[0]->antiColourNeighbour(partons[1]);
	else
	  partons[0]->    colourNeighbour(partons[1]);
      }
      if(partons[2]->data().iColour()==PDT::Colour3) {
	partons[2]->antiColourNeighbour(partons[3]);
      }
      else {
	partons[2]->    colourNeighbour(partons[3]);
      }
    } 
    // four body decay
    else {
      // generate the extra gluon
      partons.push_back(gluon->produceParticle());
      partons.back()->set5Momentum(gluon->constituentMass());
      // momentum of gluon
      pout.push_back(Lorentz5Momentum(gluon->constituentMass()));
      // generate the kinematics
      Energy2 ms2min(sqr(pout[0].mass()+pout[1].mass()+pout[2].mass()));
      Energy2 ms2max(sqr(pdec.mass()-pout[3].mass()));
      double wgt(0.);
      unsigned int ntry=0;
      bool initial = true;
      do {
	++ntry;
	Energy2 ms2 = ms2min+UseRandom::rnd()*(ms2max-ms2min);
	Energy ms  = sqrt(ms2);
	// and the W
	Energy2 mb2max = sqr(ms            -pout[2].mass());
	Energy2 mb2min = sqr(pout[0].mass()+pout[1].mass());
	Energy2 mb2 = (mb2max-mb2min)*UseRandom::rnd()+mb2min;
	wgt = (mb2max-mb2min)/(ms2max-mb2min);
	// perform first decay
	Lorentz5Momentum ps;
	double CosAngle,AzmAngle;
	Kinematics::generateAngles(CosAngle,AzmAngle);
	Kinematics::twoBodyDecay(pdec,pout[3].mass(),ms,CosAngle,AzmAngle,pout[3],ps);
	// generate the kinematics
	// perform second decay
	Kinematics::generateAngles(CosAngle,AzmAngle);
	Lorentz5Momentum p01;
	p01.setMass(sqrt(mb2));
	Kinematics::twoBodyDecay(ps,pout[2].mass(),p01.mass(),
				 CosAngle,AzmAngle,pout[2],p01);
	// perform third decay
	Kinematics::generateAngles(CosAngle,AzmAngle);
	Kinematics::twoBodyDecay(p01,pout[0].mass(),pout[1].mass(),
				 CosAngle,AzmAngle,pout[0],pout[1]);
	// kinematic piece of the weight
	wgt *= 16.*
	  Kinematics::pstarTwoBodyDecay(pdec.mass(),pout[3].mass(),ms            )/pdec.mass()*
	  Kinematics::pstarTwoBodyDecay(ms         ,p01    .mass(),pout[2].mass())/ms*
	  Kinematics::pstarTwoBodyDecay(p01 .mass(),pout[0].mass(),pout[1].mass())/p01.mass();
	wgt *= fourBodyMatrixElement(pdec,pout[2],pout[0],pout[1],pout[3],Wcol,initial);
	// check doesn't violate max
	if(wgt>_threemax) {
	  ostringstream message;
	  message << "Maximum weight for four-body decay "
		  << "violated in WeakPartonicDecayer::decay()"
		  << "Maximum = " << _fourmax << " weight = " << wgt;
	  generator()->logWarning( Exception(message.str(),Exception::warning) );
	}
      }
      while ( wgt < _fourmax*UseRandom::rnd() && ntry < _maxtry );
      if(ntry==_maxtry) throw Exception() 
	<< "Too many attempts to generate four body kinematics in "
	<< "WeakPartonicDecayer::decay()" << Exception::eventerror;
      // set momenta of particles
      for(unsigned int ix=0;ix<3;++ix) partons[ix]->setMomentum(pout[ix]);
      partons[3]->setMomentum(pspect);
      partons[4]->setMomentum(pout[3]);
      // special for tau leptons to get correlations
      if(abs(partons[0]->id())==ParticleID::tauminus||
	 abs(partons[1]->id())==ParticleID::tauminus)
	threeBodyMatrixElement(dec,pdec,partons);
      // set up the colour connections
      if(rearranged) swap(partons[1],partons[2]);
      // radiation from initial-state
      if(initial) {
	if(Wcol) {
	  if(partons[0]->data().iColour()==PDT::Colour3)
	    partons[0]->antiColourNeighbour(partons[1]);
	  else
	    partons[0]->    colourNeighbour(partons[1]);
	}
	if(partons[2]->data().iColour()==PDT::Colour3) {
	  partons[2]->antiColourNeighbour(partons[4]);
	  partons[3]->    colourNeighbour(partons[4]);
	}
	else {
	  partons[2]->    colourNeighbour(partons[4]);
	  partons[3]->antiColourNeighbour(partons[4]);
	}
      }
      // radiation from final-state
      else {
	if(partons[0]->data().iColour()==PDT::Colour3) {
	  partons[0]->antiColourNeighbour(partons[4]);
	  partons[1]->    colourNeighbour(partons[4]);
	}
	else {
	  partons[0]->    colourNeighbour(partons[4]);
	  partons[1]->antiColourNeighbour(partons[4]);
	}
	if(partons[2]->data().iColour()==PDT::Colour3) {
	  partons[2]->antiColourNeighbour(partons[3]);
	}
	else {
	  partons[2]->    colourNeighbour(partons[3]);
	}
      }
    }
  }
  return partons;
}
  
void WeakPartonicDecayer::persistentOutput(PersistentOStream & os) const {
  os << MECode << _radprob << _maxtry << _threemax << _fourmax;
}

void WeakPartonicDecayer::persistentInput(PersistentIStream & is, int) {
  is >> MECode >> _radprob >> _maxtry >> _threemax >> _fourmax;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<WeakPartonicDecayer,PartonicDecayerBase>
describeHerwigWeakPartonicDecayer("Herwig::WeakPartonicDecayer", "HwPartonicDecay.so");

void WeakPartonicDecayer::Init() {

  static ClassDocumentation<WeakPartonicDecayer> documentation
    ("The WeakPartonicDecayer class performs partonic decays of hadrons containing a "
     "heavy quark.");

  static Switch<WeakPartonicDecayer,int> interfaceMECode
    ("MECode",
     "The code for the type of matrix element to be used.",
     &WeakPartonicDecayer::MECode, 0, false, false);
  static SwitchOption interfaceMECodePhaseSpace
    (interfaceMECode,
     "PhaseSpace",
     "Phase space decays",
     0);
  static SwitchOption interfaceMECodeWeak
    (interfaceMECode,
     "Weak",
     "Weak matrix element",
     100);

  static Parameter<WeakPartonicDecayer,double> interfaceRadiationProbability
    ("RadiationProbability",
     "The probability that QCD radiation produces an extra q qbar pair",
     &WeakPartonicDecayer::_radprob, 0., 0.0, 1.,
     false, false, Interface::limited);

  static Parameter<WeakPartonicDecayer,unsigned int> interfaceMaxTry
    ("MaxTry",
     "The maximum number of attempts to generate the kinematics",
     &WeakPartonicDecayer::_maxtry, 300, 10, 1000,
     false, false, Interface::limited);

  static Parameter<WeakPartonicDecayer,double> interfaceThreeMax
    ("ThreeMax",
     "Maximum weight for sampling of three-body decays",
     &WeakPartonicDecayer::_threemax, 3.0, 1.0, 1000.0,
     false, false, Interface::limited);

  static Parameter<WeakPartonicDecayer,double> interfaceFourMax
    ("FourMax",
     "Maximum weight for sampling of four-body decays",
     &WeakPartonicDecayer::_fourmax, 3.0, 1.0, 1000.0,
     false, false, Interface::limited);

}

double WeakPartonicDecayer::VAWt(Energy2 t0, Energy2 t1, Energy2 t2, InvEnergy4 t3) {
  return (t1-t0)*(t0-t2)*t3; 
}

void WeakPartonicDecayer::dataBaseOutput(ofstream & output,
					 bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the PartonicDecayerBase base class
  PartonicDecayerBase::dataBaseOutput(output,false);
  output << "newdef " << name() << ":MECode " << MECode << " \n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}


void WeakPartonicDecayer::
threeBodyMatrixElement(tcPDPtr dec,Lorentz5Momentum & pdec,
		       ParticleVector & partons) const {
  // spinors
  LorentzSpinor   <SqrtEnergy> w0[2],w2[2];
  LorentzSpinorBar<SqrtEnergy> w1[2],w3[2];
  // spinors for the decaying particle and first product
  if(dec->id()>0) {
    SpinorWaveFunction    win(pdec,dec,0,incoming);
    w0[0] = win.dimensionedWave();
    win.reset(1);
    w0[1] = win.dimensionedWave();
    SpinorBarWaveFunction wout(partons[2]->momentum(),
			       partons[2]->dataPtr(),0,outgoing);
    w1[0] = wout.dimensionedWave();
    wout.reset(1);
    w1[1] = wout.dimensionedWave();
  }
  else {
    SpinorBarWaveFunction win(pdec,dec,0,incoming);
    w1[0] = win.dimensionedWave();
    win.reset(1);
    w1[1] = win.dimensionedWave();
    SpinorWaveFunction wout(partons[2]->momentum(),
			    partons[2]->dataPtr(),0,outgoing);
    w0[0] = wout.dimensionedWave();
    wout.reset(1);
    w0[1] = wout.dimensionedWave();
  }
  // spinors for the W decay products
  bool taufirst = true;
  if(partons[0]->id()<0) {
    SpinorWaveFunction    wout2(partons[0]->momentum(),
				partons[0]->dataPtr(),0,outgoing);
    SpinorBarWaveFunction wout3(partons[1]->momentum(),
				partons[1]->dataPtr(),0,outgoing);
    w2[0] = wout2.dimensionedWave();
    w3[0] = wout3.dimensionedWave();
    wout2.reset(1);
    wout3.reset(1);
    w2[1] = wout2.dimensionedWave();
    w3[1] = wout3.dimensionedWave();
    if(abs(partons[0]->id())!=ParticleID::tauminus) taufirst=false;
  }
  else {
    SpinorWaveFunction    wout2(partons[1]->momentum(),
				partons[1]->dataPtr(),0,outgoing);
    SpinorBarWaveFunction wout3(partons[0]->momentum(),
				partons[0]->dataPtr(),0,outgoing);
    w2[0] = wout2.dimensionedWave();
    w3[0] = wout3.dimensionedWave();
    wout2.reset(1);
    wout3.reset(1);
    w2[1] = wout2.dimensionedWave();
    w3[1] = wout3.dimensionedWave();
    if(abs(partons[1]->id())!=ParticleID::tauminus) taufirst=false;
  }
  // calculate the currents
  LorentzPolarizationVectorE Jbc[2][2],Jdec[2][2];
  for(unsigned int ix=0;ix<2;++ix) {
    for(unsigned int iy=0;iy<2;++iy) {
      Jbc [ix][iy] = w0[ix].leftCurrent(w1[iy]);
      Jdec[ix][iy] = w2[ix].leftCurrent(w3[iy]);
    }
  }
  // compute the matrix element
  Complex me[2][2][2][2];
  for(unsigned int i0=0;i0<2;++i0) {
    for(unsigned int i1=0;i1<2;++i1) {
      for(unsigned int i2=0;i2<2;++i2) {
	for(unsigned int i3=0;i3<2;++i3) {
	  me[i0][i1][i2][i3] = Jbc[i0][i1].dot(Jdec[i2][i3])/sqr(pdec.mass());
	}
      }
    }
  }
  RhoDMatrix rho(PDT::Spin1Half);
  for(unsigned int it1=0;it1<2;++it1) {
    for(unsigned int it2=0;it2<2;++it2) {
      for(unsigned int i0=0;i0<2;++i0) {
	for(unsigned int i1=0;i1<2;++i1) {
	  for(unsigned int i2=0;i2<2;++i2) {
	    rho(it1,it2) += taufirst ? 
	      me[i0][i1][it1][i2 ]*conj(me[i0][i1][it2][i2 ]) :
	      me[i0][i1][i2 ][it1]*conj(me[i0][i1][i2 ][it2]);
	  }
	}
      }
    }
  }
  // normalize matrix to unit trace
  rho.normalize();
  for(unsigned int ix=0;ix<2;++ix) {
    if(abs(partons[ix]->id())!=ParticleID::tauminus) continue;
    bool loc = partons[ix]->id() < 0;
    // create the spin info object
    FermionSpinPtr spin = new_ptr(FermionSpinInfo(partons[ix]->momentum(),true));
    // assign spinors
    for(unsigned int iy=0;iy<2;++iy) {
      spin->setBasisState(iy, loc ? w2[iy] : w3[iy].bar());
    }
    // assign rho
    spin->rhoMatrix() = rho;
    // assign spin info
    partons[ix]->spinInfo(spin);
  }
}

double WeakPartonicDecayer::
fourBodyMatrixElement(Lorentz5Momentum & p0,Lorentz5Momentum & p1,
		      Lorentz5Momentum & p2,Lorentz5Momentum & p3,
		      Lorentz5Momentum & pg, bool Wcol, bool & initial) const {
  Energy2 d01(p0*p1),d02(p0*p2),d03(p0*p3),d0g(p0*pg);
  Energy2 d12(p1*p2),d13(p1*p3),d1g(p1*pg);
  Energy2 d23(p2*p3),d2g(p2*pg),d3g(p3*pg);
  Energy2 m02(sqr(p0.mass())),m12(sqr(p1.mass())),m22(sqr(p2.mass())),
    m32(sqr(p3.mass()));
  Energy2 mei = 
    +1./d0g/d1g  *( -d01*d12*d3g+d01*d03*d2g+2*d01*d03*d12 )
    +1./d0g      *( d12*d3g-d03*d12-d02*d03 )
    +1./d1g      *( d12*d13+d03*d2g+d03*d12 )
    +m12/sqr(d1g)*( -d03*d2g-d03*d12 )
    +m02/sqr(d0g)*(  d12*d3g-d03*d12 );
  Energy2 mef = !Wcol ? ZERO : 
    +1./d2g/d3g  *( d0g*d12*d23+d03*d1g*d23+2*d03*d12*d23 )
    +1./d2g      *( d03*d1g+d03*d12-d02*d12 )
    +1./d3g      *( d0g*d12-d03*d13+d03*d12 )
    +m32/sqr(d3g)*( -d0g*d12-d03*d12 )
    +m22/sqr(d2g)*( -d03*d1g-d03*d12 );
  initial = mef/(mei+mef)<UseRandom::rnd();
  return 0.5*(mei+mef)/sqr(p0.mass()-p1.mass()-p2.mass()-p3.mass()-pg.mass());
}
