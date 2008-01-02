// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GeneralThreeBodyDecayer class.
//

#include "GeneralThreeBodyDecayer.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/PDT/ThreeBodyAllOnCalculator.h"

using namespace Herwig;

void GeneralThreeBodyDecayer::persistentOutput(PersistentOStream & os) const {
  os << _incoming << _outgoing << _diagrams << _colour << _nflow;
}

void GeneralThreeBodyDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _incoming >> _outgoing >> _diagrams >> _colour >> _nflow;
}

AbstractClassDescription<GeneralThreeBodyDecayer> GeneralThreeBodyDecayer::initGeneralThreeBodyDecayer;
// Definition of the static class description member.

void GeneralThreeBodyDecayer::Init() {

  static ClassDocumentation<GeneralThreeBodyDecayer> documentation
    ("There is no documentation for the GeneralThreeBodyDecayer class");

}

int  GeneralThreeBodyDecayer::
modeNumber(bool & cc, tcPDPtr in, const PDVector & out) const {
  // check number of outgoing particles
  if(out.size()!=3) return -1;
  // check incoming particle
  if(abs(in->id())!=_incoming->id()) return false;
  // check outgoing particles
  pair<bool,bool> allowed=make_pair(true,true);
  for(unsigned int ix=0;ix<3;++ix) {
    if(out[ix]!=_outgoing[ix]) 
      allowed.first  = false;
    if( ( out[ix]->CC() && out[ix]->CC() == _outgoing[ix] ) ||
	(!out[ix]->CC() && out[ix]       == _outgoing[ix] ) )
      allowed.second = false;
  }
  if(allowed.first  && in->id() == _incoming->id()) {
    cc = false;
    return 0;
  }
  if(allowed.second && ( (!in->CC() && in->id() ==  _incoming->id()) ||
			 ( in->CC() && in->id() == -_incoming->id()) )) {
    cc = true;
    return 0.;
  }
  return -1;
}

void GeneralThreeBodyDecayer::setDecayInfo(PDPtr incoming,
					   vector<PDPtr> outgoing,
					   const vector<TBDiagram> & process,
					   const vector<DVector> & factors,
					   const unsigned int ncf) {
  // set the member variables from the info supplied
  _incoming = incoming;
  _outgoing = outgoing;
  _diagrams = process;
  _colour   = factors;
  _nflow    = ncf;
  cerr << "testing in the setdecayinfo " << _diagrams.size() << "\n";
  for(unsigned int ix=0;ix<_outgoing.size();++ix) {
    cerr << "testing outgoing " << _outgoing[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_diagrams.size();++ix) {
    cerr << "testing in the loop " << _diagrams[ix].outgoing 
	 << " " << _diagrams[ix].channelType << "\n";
    unsigned int iy=0;
    for(;iy<3;++iy) 
      if(_diagrams[ix].outgoing == _outgoing[iy]->id()) break;
    if(_diagrams[ix].channelType == TBDiagram::UNDEFINED) {
      _diagrams[ix].channelType = TBDiagram::Channel(iy);
      if( ( iy == 0 && outgoing[1]->id() != _diagrams[ix].outgoingPair.first)||
	  ( iy == 1 && outgoing[0]->id() != _diagrams[ix].outgoingPair.first)|| 
	  ( iy == 2 && outgoing[0]->id() != _diagrams[ix].outgoingPair.first) ) 
	swap(_diagrams[ix].outgoingPair.first, _diagrams[ix].outgoingPair.second);
    }
  }
}

void GeneralThreeBodyDecayer::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  // create the phase space integrator
  PDVector extpart(1,_incoming);
  extpart.insert(extpart.end(),_outgoing.begin(),_outgoing.end());
  // create the integration channels for the decay
  DecayPhaseSpaceModePtr mode(new_ptr(DecayPhaseSpaceMode(extpart,this)));
  DecayPhaseSpaceChannelPtr newchannel;
  // create the phase-space channels for the integration
  unsigned int nmode(0);
  for(unsigned int ix=0;ix<_diagrams.size();++ix) {
    if(_diagrams[ix].channelType==TBDiagram::fourPoint||
       _diagrams[ix].channelType==TBDiagram::UNDEFINED) continue;
    // create the new channel
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    if(_diagrams[ix].channelType==TBDiagram::channel23) {
      newchannel->addIntermediate(extpart[0],0,0.0,-1,1);
      newchannel->addIntermediate(_diagrams[ix].intermediate,0,0.0, 2,3);
    }
    else if(_diagrams[ix].channelType==TBDiagram::channel13) {
      newchannel->addIntermediate(extpart[0],0,0.0,-1,2);
      newchannel->addIntermediate(_diagrams[ix].intermediate,0,0.0, 1,3);
    }
    else if(_diagrams[ix].channelType==TBDiagram::channel12) {
      newchannel->addIntermediate(extpart[0],0,0.0,-1,3);
      newchannel->addIntermediate(_diagrams[ix].intermediate,0,0.0, 1,2);
    }
    mode->addChannel(newchannel);
    ++nmode;
  }
  // add the mode
  vector<double> wgt(nmode,1./double(nmode));
  addMode(mode,1.,wgt);
  cerr << *this;
}

double GeneralThreeBodyDecayer::
threeBodyMatrixElement(const int imode,  const Energy2 q2,
		       const Energy2 s3, const Energy2 s2, 
		       const Energy2 s1, const Energy  m1, 
		       const Energy  m2, const Energy  m3) const {
  // calculate the momenta of the outgoing particles
  Energy m0=sqrt(q2);
  // energies
  Energy eout[3] = {0.5*(q2+sqr(m1)-s1)/m0,
		    0.5*(q2+sqr(m2)-s2)/m0,
		    0.5*(q2+sqr(m3)-s3)/m0};
  // magnitudes of the momenta
  Energy pout[3] = {sqrt(sqr(eout[0])-sqr(m1)),
		    sqrt(sqr(eout[1])-sqr(m2)),
		    sqrt(sqr(eout[2])-sqr(m3))};
  double cos2 = 0.5*(sqr(pout[0])+sqr(pout[1])-sqr(pout[2]))/pout[0]/pout[1];
  double cos3 = 0.5*(sqr(pout[0])-sqr(pout[1])+sqr(pout[2]))/pout[0]/pout[2];
  double sin2 = sqrt(1.-sqr(cos2)), sin3 = sqrt(1.-sqr(cos3));
  Lorentz5Momentum out[3]=
    {Lorentz5Momentum(      0.*GeV   , 0.*GeV ,  pout[0]      , eout[0] , m1),
     Lorentz5Momentum(  pout[1]*sin2 , 0.*GeV , -pout[1]*cos2 , eout[1] , m2),
     Lorentz5Momentum( -pout[2]*sin3 , 0.*GeV , -pout[2]*cos3 , eout[2] , m3)};
  // create the incoming
  PPtr inpart=mode(imode)->externalParticles(0)->
    produceParticle(Lorentz5Momentum(sqrt(q2)));
  // and outgoing particles
  ParticleVector decay;
  for(unsigned int ix=1;ix<4;++ix)
    decay.push_back(mode(imode)->externalParticles(ix)->produceParticle(out[ix-1]));
  // return the matrix element
  return me2(false,-1,*inpart,decay);
}

double GeneralThreeBodyDecayer::brat(const DecayMode &, const Particle & p,
				     double oldbrat) const {
  ParticleVector children = p.children();
  if( children.size() != 3 || !p.data().widthGenerator() ) 
    return oldbrat;
  // partial width for this mode
  Energy scale = p.mass();
  Energy pwidth = 
    partialWidth( make_pair(p.dataPtr(), scale),
		  make_pair(children[0]->dataPtr(), children[0]->mass()),
		  make_pair(children[1]->dataPtr(), children[1]->mass()),
		  make_pair(children[2]->dataPtr(), children[2]->mass()) );
  Energy width = p.data().widthGenerator()->width(p.data(), scale);
  return pwidth/width;
}

Energy GeneralThreeBodyDecayer::partialWidth(PMPair inpart, PMPair outa, 
					     PMPair outb, PMPair outc) const {
  if(inpart.second<outa.second+outb.second+outc.second) return 0.*GeV;
  // create the object to calculate the width if it doesn't all ready exist
  if(!_widthcalc) {
    string tag = _incoming->PDGName() + "->";
    tag += _outgoing[0]->PDGName() + "," + _outgoing[1]->PDGName() + ","
      + _outgoing[2]->PDGName() + ";";
    DMPtr dm = generator()->findDecayMode(tag);
    _widthcalc = threeBodyMEIntegrator(*dm);
  }
  return _widthcalc->partialWidth(sqr(inpart.second));
}

