
// -*- C++ -*-
//
// GeneralTwoBodyDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GeneralTwoBodyDecayer class.
//

#include "GeneralTwoBodyDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Utilities/Exception.h"
#include "Herwig/Shower/RealEmissionProcess.h"

using namespace Herwig;

ParticleVector GeneralTwoBodyDecayer::decay(const Particle & parent,
					    const tPDVector & children) const {

  // return empty vector if products heavier than parent
  Energy mout(ZERO);
  for(tPDVector::const_iterator it=children.begin();
      it!=children.end();++it) mout+=(**it).massMin();
  if(mout>parent.mass()) return ParticleVector();
  // generate the decay
  bool cc;
  int imode=modeNumber(cc,parent.dataPtr(),children);
  // generate the kinematics
  ParticleVector decay=generate(generateIntermediates(),cc,imode,parent);
  // make the colour connections
  colourConnections(parent, decay);
  // return the answer
  return decay;
}

void GeneralTwoBodyDecayer::doinit() {
  DecayIntegrator::doinit();
  assert( vertex_ );
  assert( _incoming && _outgoing.size()==2);
  vertex_->init();
  //create phase space mode
  tPDVector extpart(3);
  extpart[0] = _incoming;
  extpart[1] = _outgoing[0];
  extpart[2] = _outgoing[1];
  addMode(new_ptr(DecayPhaseSpaceMode(extpart, this)), _maxweight, vector<double>());
}

int GeneralTwoBodyDecayer::modeNumber(bool & cc, tcPDPtr parent, 
				      const tPDVector & children) const {
  long parentID = parent->id();
  long id1 = children[0]->id();
  long id2 = children[1]->id();
  cc = false;
  long out1 = _outgoing[0]->id();
  long out2 = _outgoing[1]->id();
  if( parentID == _incoming->id() && 
      ((id1 == out1 && id2 == out2) || 
       (id1 == out2 && id2 == out1)) ) {
    return 0;
  }
  else if(_incoming->CC() && parentID == _incoming->CC()->id()) {
    cc = true;
    if( _outgoing[0]->CC()) out1 = _outgoing[0]->CC()->id();
    if( _outgoing[1]->CC()) out2 = _outgoing[1]->CC()->id();
    if((id1 == out1 && id2 == out2) || 
       (id1 == out2 && id2 == out1)) return 0;
  }
  return -1;
}

void GeneralTwoBodyDecayer::
colourConnections(const Particle & parent,
		  const ParticleVector & out) const {
  PDT::Colour incColour(parent.data().iColour());
  PDT::Colour outaColour(out[0]->data().iColour());
  PDT::Colour outbColour(out[1]->data().iColour());
  //incoming colour singlet
  if(incColour == PDT::Colour0) {
    // colour triplet-colourantitriplet
    if((outaColour == PDT::Colour3 && outbColour == PDT::Colour3bar) ||
       (outaColour == PDT::Colour3bar && outbColour == PDT::Colour3)) {
      bool ac(out[0]->id() < 0);
      out[0]->colourNeighbour(out[1],!ac);
    }
    //colour octet
    else if(outaColour == PDT::Colour8 && outbColour == PDT::Colour8) {
      out[0]->colourNeighbour(out[1]);
      out[0]->antiColourNeighbour(out[1]);
    }
    // colour singlets
    else if(outaColour == PDT::Colour0 && outbColour == PDT::Colour0) {
    }
    // unknown
    else
      throw Exception() << "Unknown outgoing colours for decaying "
			<< "colour singlet in "
			<< "GeneralTwoBodyDecayer::colourConnections "
			<< outaColour << " " << outbColour
			<< Exception::runerror; 
  }
  //incoming colour triplet
  else if(incColour == PDT::Colour3) {
    // colour triplet + singlet
    if(outaColour == PDT::Colour3 && outbColour == PDT::Colour0) {
      out[0]->incomingColour(const_ptr_cast<tPPtr>(&parent));
    }
    //opposite order
    else if(outaColour == PDT::Colour0 && outbColour == PDT::Colour3) {
      out[1]->incomingColour(const_ptr_cast<tPPtr>(&parent));
    }
    // octet + triplet
    else if(outaColour == PDT::Colour8 && outbColour == PDT::Colour3) {
      out[0]->incomingColour(const_ptr_cast<tPPtr>(&parent));
      out[1]->antiColourNeighbour(out[0]);
    }
    //opposite order
    else if(outaColour == PDT::Colour3 && outbColour == PDT::Colour8) {
      out[1]->incomingColour(const_ptr_cast<tPPtr>(&parent));
      out[0]->antiColourNeighbour(out[1]);
    }
    else if(outaColour == PDT::Colour3bar && outaColour == PDT::Colour3bar) {
      tColinePtr col[2] = {ColourLine::create(out[0],true),
			   ColourLine::create(out[1],true)};
      parent.colourLine()->setSinkNeighbours(col[0],col[1]);
    }
    else
      throw Exception() << "Unknown outgoing colours for decaying "
			<< "colour triplet in "
			<< "GeneralTwoBodyDecayer::colourConnections() "
			<< outaColour << " " << outbColour
			<< Exception::runerror; 
  }
  // incoming colour anti triplet
  else if(incColour == PDT::Colour3bar) {
    // colour antitriplet +singlet
    if(outaColour == PDT::Colour3bar && outbColour == PDT::Colour0) {
      out[0]->incomingAntiColour(const_ptr_cast<tPPtr>(&parent));
    }
    //opposite order
    else if(outaColour == PDT::Colour0 && outbColour == PDT::Colour3bar) {
      out[1]->incomingAntiColour(const_ptr_cast<tPPtr>(&parent));
    }
    //octet + antitriplet
    else if(outaColour == PDT::Colour3bar && outbColour == PDT::Colour8) {
      out[1]->incomingAntiColour(const_ptr_cast<tPPtr>(&parent));
      out[0]->colourNeighbour(out[1]);
    }
    //opposite order
    else if(outaColour == PDT::Colour8 && outbColour == PDT::Colour3bar) {
      out[0]->incomingAntiColour(const_ptr_cast<tPPtr>(&parent));
      out[1]->colourNeighbour(out[0]);
    }
    else if(outaColour == PDT::Colour3 && outbColour == PDT::Colour3) {
      tColinePtr col[2] = {ColourLine::create(out[0]),
			   ColourLine::create(out[1])};
      parent.antiColourLine()->setSourceNeighbours(col[0],col[1]);
    }
    else
      throw Exception() << "Unknown outgoing colours for decaying "
			<< "colour antitriplet "
			<< "in GeneralTwoBodyDecayer::colourConnections() "
			<< outaColour << " " << outbColour
			<< Exception::runerror; 
  }
  //incoming colour octet
  else if(incColour == PDT::Colour8) {
    // triplet-antitriplet
    if(outaColour == PDT::Colour3&&outbColour == PDT::Colour3bar) {
      out[0]->incomingColour(const_ptr_cast<tPPtr>(&parent));
      out[1]->incomingAntiColour(const_ptr_cast<tPPtr>(&parent));
    }
    // opposite order
    else if(outbColour == PDT::Colour3&&outaColour == PDT::Colour3bar) {
      out[0]->incomingAntiColour(const_ptr_cast<tPPtr>(&parent));
      out[1]->incomingColour(const_ptr_cast<tPPtr>(&parent));
    }
    // neutral octet
    else if(outaColour == PDT::Colour0&&outbColour == PDT::Colour8) {
      out[1]->incomingColour(const_ptr_cast<tPPtr>(&parent));
      out[1]->incomingAntiColour(const_ptr_cast<tPPtr>(&parent));
    }
    else if(outbColour == PDT::Colour0&&outaColour == PDT::Colour8) {
      out[0]->incomingColour(const_ptr_cast<tPPtr>(&parent));
      out[0]->incomingAntiColour(const_ptr_cast<tPPtr>(&parent));
    }
    else
      throw Exception() << "Unknown outgoing colours for decaying "
                        << "colour octet "
                        << "in GeneralTwoBodyDecayer::colourConnections() "
                        << outaColour << " " << outbColour
                        << Exception::runerror;
  }
  else if(incColour == PDT::Colour6) {
    if(outaColour == PDT::Colour3 && outbColour == PDT::Colour3) {
      tPPtr tempParent = const_ptr_cast<tPPtr>(&parent);
      Ptr<MultiColour>::pointer parentColour = 
      	dynamic_ptr_cast<Ptr<MultiColour>::pointer>
      	(tempParent->colourInfo());

      tColinePtr line1 = const_ptr_cast<tColinePtr>(parentColour->colourLines()[0]);
      line1->addColoured(dynamic_ptr_cast<tPPtr>(out[0]));

      tColinePtr line2 = const_ptr_cast<tColinePtr>(parentColour->colourLines()[1]);
      line2->addColoured(dynamic_ptr_cast<tPPtr>(out[1]));
    }
    else
      throw Exception() << "Unknown outgoing colours for decaying "
                        << "colour sextet "
                        << "in GeneralTwoBodyDecayer::colourConnections() "
                        << outaColour << " " << outbColour
                        << Exception::runerror;
  }
  else if(incColour == PDT::Colour6bar) {
    if(outaColour == PDT::Colour3bar && outbColour == PDT::Colour3bar) {
     tPPtr tempParent = const_ptr_cast<tPPtr>(&parent);
      Ptr<MultiColour>::pointer parentColour = 
      	dynamic_ptr_cast<Ptr<MultiColour>::pointer>
      	(tempParent->colourInfo());

      tColinePtr line1 = const_ptr_cast<tColinePtr>(parentColour->antiColourLines()[0]);
      line1->addAntiColoured(dynamic_ptr_cast<tPPtr>(out[0]));

      tColinePtr line2 = const_ptr_cast<tColinePtr>(parentColour->antiColourLines()[1]);
      line2->addAntiColoured(dynamic_ptr_cast<tPPtr>(out[1]));
    }
    else
      throw Exception() << "Unknown outgoing colours for decaying "
                        << "colour anti-sextet "
                        << "in GeneralTwoBodyDecayer::colourConnections() "
                        << outaColour << " " << outbColour
                        << Exception::runerror;
  }
  else
    throw Exception() << "Unknown incoming colour in "
		      << "GeneralTwoBodyDecayer::colourConnections() "
		      << incColour
		      << Exception::runerror; 
}

bool GeneralTwoBodyDecayer::twoBodyMEcode(const DecayMode & dm, int & mecode,
					  double & coupling) const {
  assert(dm.parent()->id() == _incoming->id());
  ParticleMSet::const_iterator pit = dm.products().begin();
  long id1 = (*pit)->id();
  ++pit;
  long id2 = (*pit)->id();
  long id1t(_outgoing[0]->id()), id2t(_outgoing[1]->id());
  mecode = -1;
  coupling = 1.;
  if( id1 == id1t && id2 == id2t ) {
    return true;
  }
  else if( id1 == id2t && id2 == id1t ) {
    return false;
  }
  else
    assert(false);
  return false;
}


void GeneralTwoBodyDecayer::persistentOutput(PersistentOStream & os) const {
  os << vertex_ << _incoming << _outgoing << _maxweight << ounit(pTmin_,GeV)
     << coupling_ << incomingVertex_ << outgoingVertices_ << fourPointVertex_;
}

void GeneralTwoBodyDecayer::persistentInput(PersistentIStream & is, int) {
  is >> vertex_ >> _incoming >> _outgoing >> _maxweight >> iunit(pTmin_,GeV) 
     >> coupling_ >> incomingVertex_ >> outgoingVertices_ >> fourPointVertex_;
}

AbstractClassDescription<GeneralTwoBodyDecayer> 
GeneralTwoBodyDecayer::initGeneralTwoBodyDecayer;
// Definition of the static class description member.

void GeneralTwoBodyDecayer::Init() {

  static ClassDocumentation<GeneralTwoBodyDecayer> documentation
    ("This class is designed to be a base class for all 2 body decays"
     "in a general model");

  static Parameter<GeneralTwoBodyDecayer,Energy> interfacepTmin
    ("pTmin",
     "Minimum transverse momentum from gluon radiation",
     &GeneralTwoBodyDecayer::pTmin_, GeV, 1.0*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);
 
  static Reference<GeneralTwoBodyDecayer,ShowerAlpha> interfaceCoupling
    ("Coupling",
     "Object for the coupling in the generation of hard radiation",
     &GeneralTwoBodyDecayer::coupling_, false, false, true, false, false);

}

double GeneralTwoBodyDecayer::brat(const DecayMode &, const Particle & p,
				   double oldbrat) const {
  ParticleVector children = p.children();
  if( children.size() != 2 || !p.data().widthGenerator() ) 
    return oldbrat;
  
  // partial width for this mode
  Energy scale = p.mass();
  Energy pwidth = 
    partialWidth( make_pair(p.dataPtr(), scale),
		  make_pair(children[0]->dataPtr(), children[0]->mass()),
		  make_pair(children[1]->dataPtr(), children[1]->mass()) );
  Energy width = p.data().widthGenerator()->width(p.data(), scale);
  return pwidth/width;
}

void GeneralTwoBodyDecayer::doinitrun() {
  vertex_->initrun();
  DecayIntegrator::doinitrun();
  for(unsigned int ix=0;ix<numberModes();++ix) {
    double fact = pow(1.5,int(mode(ix)->externalParticles(0)->iSpin())-1);
    mode(ix)->setMaxWeight(fact*mode(ix)->maxWeight());
  }
}

double GeneralTwoBodyDecayer::colourFactor(tcPDPtr in, tcPDPtr out1,
					   tcPDPtr out2) const {
  // identical particle symmetry factor
  double output = out1->id()==out2->id() ? 0.5 : 1.;
  // colour neutral incoming particle
  if(in->iColour()==PDT::Colour0) {
    // both colour neutral
    if(out1->iColour()==PDT::Colour0 && out2->iColour()==PDT::Colour0)
      output *= 1.;
    // colour triplet/ antitriplet
    else if((out1->iColour()==PDT::Colour3    && out2->iColour()==PDT::Colour3bar) ||
	    (out1->iColour()==PDT::Colour3bar && out2->iColour()==PDT::Colour3   ) ) {
      output *= 3.;
    }
    // colour octet colour octet
    else if(out1->iColour()==PDT::Colour8 && out2->iColour()==PDT::Colour8 ) {
      output *= 8.;
    }
    else 
      throw Exception() << "Unknown colour for the outgoing particles"
			<< " for decay colour neutral particle in "
			<< "GeneralTwoBodyDecayer::colourFactor() for "
			<< in->PDGName() << " -> "
			<< out1->PDGName() << " " << out2->PDGName() 
			<< Exception::runerror;
  }
  // triplet
  else if(in->iColour()==PDT::Colour3) {
    // colour triplet + neutral
    if((out1->iColour()==PDT::Colour0 && out2->iColour()==PDT::Colour3) ||
       (out1->iColour()==PDT::Colour3 && out2->iColour()==PDT::Colour0) ) {
      output *= 1.;
    }
    // colour triplet + octet
    else if((out1->iColour()==PDT::Colour8 && out2->iColour()==PDT::Colour3) ||
	    (out1->iColour()==PDT::Colour3 && out2->iColour()==PDT::Colour8) ) {
      output *= 4./3.;
    }
    // colour anti triplet anti triplet
    else if(out1->iColour()==PDT::Colour3bar && 
	    out2->iColour()==PDT::Colour3bar) {
      output *= 2.;
    }
    else
      throw Exception() << "Unknown colour for the outgoing particles"
			<< " for decay colour triplet particle in "
			<< "GeneralTwoBodyDecayer::colourFactor() for "
			<< in->PDGName() << " -> "
			<< out1->PDGName() << " " << out2->PDGName() 
			<< Exception::runerror;
  }
  // anti triplet
  else if(in->iColour()==PDT::Colour3bar) {
    // colour anti triplet + neutral
    if((out1->iColour()==PDT::Colour0    && out2->iColour()==PDT::Colour3bar ) ||
       (out1->iColour()==PDT::Colour3bar && out2->iColour()==PDT::Colour0    ) ) {
      output *= 1.;
    }
    // colour anti triplet + octet
    else if((out1->iColour()==PDT::Colour8    && out2->iColour()==PDT::Colour3bar ) ||
	    (out1->iColour()==PDT::Colour3bar && out2->iColour()==PDT::Colour8    ) ) {
      output *= 4./3.;
    }
    // colour triplet triplet
    else if(out1->iColour()==PDT::Colour3 && 
	    out2->iColour()==PDT::Colour3) {
      output *= 2.;
    }
    else
      throw Exception() << "Unknown colour for the outgoing particles"
			<< " for decay colour anti triplet particle in "
			<< "GeneralTwoBodyDecayer::colourFactor() for "
			<< in->PDGName() << " -> "
			<< out1->PDGName() << " " << out2->PDGName() 
			<< Exception::runerror;
  }
  else if(in->iColour()==PDT::Colour8) {
    // colour octet + neutral
    if((out1->iColour()==PDT::Colour0 && out2->iColour()==PDT::Colour8 ) ||
       (out1->iColour()==PDT::Colour8 && out2->iColour()==PDT::Colour0 ) ) {
      output *= 1.;
    }
    // colour triplet/antitriplet
    else if((out1->iColour()==PDT::Colour3    && out2->iColour()==PDT::Colour3bar) ||
	    (out1->iColour()==PDT::Colour3bar && out2->iColour()==PDT::Colour3   ) ) {
      output *= 0.5;
    }
    else
      throw Exception() << "Unknown colour for the outgoing particles"
			<< " for decay colour octet particle in "
			<< "GeneralTwoBodyDecayer::colourFactor() for "
			<< in->PDGName() << " -> "
			<< out1->PDGName() << " " << out2->PDGName() 
			<< Exception::runerror;
  }
  else if(in->iColour()==PDT::Colour6) {
    // colour sextet -> triplet triplet
    if( out1->iColour()==PDT::Colour3 && out2->iColour()==PDT::Colour3 ) {
      output *= 1.;
    }
    else
      throw Exception() << "Unknown colour for the outgoing particles"
			<< " for decay colour sextet particle in "
			<< "GeneralTwoBodyDecayer::colourFactor() for "
			<< in->PDGName() << " -> "
			<< out1->PDGName() << " " << out2->PDGName() 
			<< Exception::runerror;
  }
  else if(in->iColour()==PDT::Colour6bar) {
    // colour sextet -> triplet triplet
    if( out1->iColour()==PDT::Colour3bar && out2->iColour()==PDT::Colour3bar ) {
      output *= 1.;
    }
    else
      throw Exception() << "Unknown colour for the outgoing particles"
			<< " for decay colour anti-sextet particle in "
			<< "GeneralTwoBodyDecayer::colourFactor() for "
			<< in->PDGName() << " -> "
			<< out1->PDGName() << " " << out2->PDGName() 
			<< Exception::runerror;
  }
  else
    throw Exception() << "Unknown colour "
		      << in->iColour() << " for the decaying particle in "
		      << "GeneralTwoBodyDecayer::colourFactor() for "
		      << in->PDGName() << " -> "
		      << out1->PDGName() << " " << out2->PDGName() 
		      << Exception::runerror;
  return output;
}

Energy GeneralTwoBodyDecayer::partialWidth(PMPair inpart, PMPair outa, 
					   PMPair outb) const {
  // select the number of the mode
  tPDVector children;
  children.push_back(const_ptr_cast<PDPtr>(outa.first));
  children.push_back(const_ptr_cast<PDPtr>(outb.first));
  bool cc;
  int nmode=modeNumber(cc,inpart.first,children);
  tcPDPtr newchild[2] = {mode(nmode)->externalParticles(1),
			 mode(nmode)->externalParticles(2)};
  // make the particles
  Lorentz5Momentum pparent = Lorentz5Momentum(inpart.second);
  PPtr parent = inpart.first->produceParticle(pparent);
  Lorentz5Momentum pout[2];
  double ctheta,phi;
  Kinematics::generateAngles(ctheta,phi);
  Kinematics::twoBodyDecay(pparent, outa.second, outb.second,
 			   ctheta, phi,pout[0],pout[1]);
  if( ( !cc && outa.first!=newchild[0]) ||
      (  cc && !((  outa.first->CC() && outa.first->CC() == newchild[0])||
		 ( !outa.first->CC() && outa.first       == newchild[0]) )))
    swap(pout[0],pout[1]);
  ParticleVector decay;
  decay.push_back(newchild[0]->produceParticle(pout[0]));
  decay.push_back(newchild[1]->produceParticle(pout[1]));
  double me =  me2(-1,*parent,decay,Initialize);
  Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second,
					     outa.second, outb.second);
  return me/(8.*Constants::pi)*pcm;
}

void GeneralTwoBodyDecayer::setDecayInfo(PDPtr incoming,PDPair outgoing,
					 VertexBasePtr vertex, VertexBasePtr inV,
					 const vector<VertexBasePtr> & outV,
					 VertexBasePtr fourV) {
  _incoming=incoming;
  _outgoing.clear();
  _outgoing.push_back(outgoing.first );
  _outgoing.push_back(outgoing.second);
  vertex_ = vertex;
  incomingVertex_ = inV;
  outgoingVertices_ = outV;
  fourPointVertex_ = fourV;

}

RealEmissionProcessPtr GeneralTwoBodyDecayer::generateHardest(RealEmissionProcessPtr born) {
  // check one incoming
  assert(born->bornIncoming().size()==1);
  // check exactly two outgoing particles
  assert(born->bornOutgoing().size()==2);
  // ignore effective vertices
  if (vertex_ && (vertex_->orderInGem()+vertex_->orderInGs())>1) 
    return RealEmissionProcessPtr();
  // search for coloured particles
  bool colouredParticles=born->bornIncoming()[0]->dataPtr()->coloured();
  if(!colouredParticles) {
    for(unsigned int ix=0;ix<born->bornOutgoing().size();++ix) {
      if(born->bornOutgoing()[ix]->dataPtr()->coloured()) {
	colouredParticles=true;
	break;
      }
    }
  }
  // if no coloured particles return
  if ( !colouredParticles ) return RealEmissionProcessPtr();
  // for decay b -> a c 
  // set progenitors
  PPtr cProgenitor = born->bornOutgoing()[0];
  PPtr aProgenitor = born->bornOutgoing()[1];
  // get the decaying particle
  PPtr bProgenitor = born->bornIncoming()[0];

  // identify which dipoles are required
  vector<dipoleType> dipoles;
  if(!identifyDipoles(dipoles,aProgenitor,bProgenitor,cProgenitor)) {
    return RealEmissionProcessPtr();
  }

  Energy trialpT = pTmin_;
  LorentzRotation eventFrame;
  vector<Lorentz5Momentum> momenta;
  vector<Lorentz5Momentum> trialMomenta(4);
  PPtr finalEmitter, finalSpectator;
  PPtr trialEmitter, trialSpectator;
  dipoleType finalType(FFa);

  for (int i=0; i<int(dipoles.size()); ++i){

    // assign emitter and spectator based on current dipole
    if (dipoles[i]==FFc || dipoles[i]==IFc || dipoles[i]==IFbc){
      trialEmitter   = cProgenitor;
      trialSpectator = aProgenitor;
    }
    else if (dipoles[i]==FFa || dipoles[i]==IFa || dipoles[i]==IFba){
      trialEmitter   = aProgenitor;
      trialSpectator = cProgenitor;
    }

    // find rotation from lab to frame with the spectator along -z
    LorentzRotation trialEventFrame(bProgenitor->momentum().findBoostToCM());
    Lorentz5Momentum pspectator = (trialEventFrame*trialSpectator->momentum());
    trialEventFrame.rotateZ( -pspectator.phi() );
    trialEventFrame.rotateY( -pspectator.theta() - Constants::pi );
    // invert it
    trialEventFrame.invert();

    // try to generate an emission
    pT_ = pTmin_;
    vector<Lorentz5Momentum> trialMomenta 
      = hardMomenta(bProgenitor, trialEmitter, trialSpectator, dipoles, i);
  
    // select dipole which gives highest pT emission
    if(pT_>trialpT) {
      trialpT        = pT_;
      momenta        = trialMomenta;
      eventFrame     = trialEventFrame;
      finalEmitter   = trialEmitter;
      finalSpectator = trialSpectator;
      finalType      = dipoles[i];

      if (dipoles[i]==FFc || dipoles[i]==FFa ) {
      	if((momenta[3]+momenta[1]).m2()-momenta[1].m2()>(momenta[3]+momenta[2]).m2()-momenta[2].m2()) {
      	  swap(finalEmitter,finalSpectator);
      	  swap(momenta[1],momenta[2]);
      	}
      }
    }
  }
  pT_ = trialpT;
  // if no emission return
  if(momenta.empty()) {
    born->pT()[ShowerInteraction::QCD] = pTmin_;
    return born;
  }

  // rotate momenta back to the lab
  for(unsigned int ix=0;ix<momenta.size();++ix) {
    momenta[ix] *= eventFrame;
  }
 
  // set maximum pT for subsequent branchings
  born->pT()[ShowerInteraction::QCD] = pT_;

  // get ParticleData objects
  tcPDPtr b = bProgenitor   ->dataPtr();
  tcPDPtr e = finalEmitter  ->dataPtr();
  tcPDPtr s = finalSpectator->dataPtr();
  tcPDPtr gluon  = getParticleData(ParticleID::g);

  // create new ShowerParticles
  PPtr emitter   = e    ->produceParticle(momenta[1]);
  PPtr spectator = s    ->produceParticle(momenta[2]);
  PPtr gauge     = gluon->produceParticle(momenta[3]);
  PPtr incoming  = b    ->produceParticle(bProgenitor->momentum());

  // insert the particles
  born->incoming().push_back(incoming);
  unsigned int iemit(0),ispect(0);
  for(unsigned int ix=0;ix<born->bornOutgoing().size();++ix) {
    if(born->bornOutgoing()[ix]==finalEmitter) {
      born->outgoing().push_back(emitter);
      iemit = born->outgoing().size();
    }
    else if(born->bornOutgoing()[ix]==finalSpectator) {
      born->outgoing().push_back(spectator);
      ispect = born->outgoing().size();
    }
  }
  born->outgoing().push_back(gauge);
  if(!spectator->dataPtr()->coloured() ||
     (finalType != FFa && finalType!=FFc) ) ispect = 0;
  born->emitter(iemit);
  born->spectator(ispect);
  born->emitted(3);
  // set up colour lines
  getColourLines(born);
  // return the tree
  born->interaction(ShowerInteraction::QCD);
  return born;
}

double GeneralTwoBodyDecayer::threeBodyME(const int , const Particle &,
					  const ParticleVector &,MEOption) {
  throw Exception() << "Base class  GeneralTwoBodyDecayer::threeBodyME() "
		    << "called, should have an implementation in the inheriting class"
		    << Exception::runerror;
  return 0.;
}

vector<Lorentz5Momentum>  GeneralTwoBodyDecayer::hardMomenta(const PPtr &in, 
							     const PPtr &emitter, 
							     const PPtr &spectator, 
							     const vector<dipoleType>  &dipoles, 
							     int i) {
  double C    = 6.3;
  double ymax = 10.;
  double ymin = -ymax;

  // get masses of the particles
  mb_  = in       ->momentum().mass();
  e_   = emitter  ->momentum().mass()/mb_;
  s_   = spectator->momentum().mass()/mb_;
  e2_  = sqr(e_);
  s2_  = sqr(s_);
  
  vector<Lorentz5Momentum> particleMomenta (4);
  Energy2 lambda = sqr(mb_)*sqrt(1.+sqr(s2_)+sqr(e2_)-2.*s2_-2.*e2_-2.*s2_*e2_);    

  // calculate A
  double A = (ymax-ymin)*C*(coupling()->overestimateValue()/(2.*Constants::pi)); 
  Energy pTmax = mb_*(sqr(1.-s_)-e2_)/(2.*(1.-s_));
  
  // if no possible branching return
  if ( pTmax < pTmin_ ) {
    particleMomenta.clear(); 
    return particleMomenta;
  }

  while (pTmax >= pTmin_) {  
    // generate pT, y and phi values
    Energy pT = pTmax*pow(UseRandom::rnd(),(1./A)); 
    if (pT < pTmin_) {
      particleMomenta.clear(); 
      break;
    }

    double phi = UseRandom::rnd()*Constants::twopi;
    double y   = ymin+UseRandom::rnd()*(ymax-ymin);

    double weight[2] = {0.,0.};
    double xs[2], xe[2], xe_z[2], xg;
    
    for (unsigned int j=0; j<2; j++) {

      // check if the momenta are physical
      if (!calcMomenta(j, pT, y, phi, xg, xs[j], xe[j], xe_z[j], particleMomenta)) 
	continue;
      
      // check if point lies within phase space
      if (!psCheck(xg, xs[j])) 
	continue;
   
      // decay products for 3 body decay
      PPtr inpart   = in        ->dataPtr()->produceParticle(particleMomenta[0]);     
      ParticleVector decay3;
      decay3.push_back(emitter  ->dataPtr()->produceParticle(particleMomenta[1]));
      decay3.push_back(spectator->dataPtr()->produceParticle(particleMomenta[2]));
      decay3.push_back(getParticleData(ParticleID::g    )->produceParticle(particleMomenta[3]));
      
      // decay products for 2 body decay
      Lorentz5Momentum p1(ZERO,ZERO, lambda/2./mb_,(mb_/2.)*(1.+e2_-s2_),mb_*e_);
      Lorentz5Momentum p2(ZERO,ZERO,-lambda/2./mb_,(mb_/2.)*(1.+s2_-e2_),mb_*s_);
      ParticleVector decay2;
      decay2.push_back(emitter  ->dataPtr()->produceParticle(p1));
      decay2.push_back(spectator->dataPtr()->produceParticle(p2));
      if (decay2[0]->dataPtr()->iSpin()!=PDT::Spin1Half &&
	  decay2[1]->dataPtr()->iSpin()==PDT::Spin1Half) swap(decay2[0], decay2[1]);
  
      // calculate matrix element ratio R/B
      double meRatio = matrixElementRatio(*inpart,decay2,decay3,Initialize);
   
      // calculate dipole factor
      InvEnergy2 dipoleSum = ZERO;
      InvEnergy2 numerator = ZERO;
      for (int k=0; k<int(dipoles.size()); ++k){
	InvEnergy2 dipole = abs(calculateDipole(dipoles[k],*inpart,decay3,dipoles[i]));
	dipoleSum += dipole;
	if (k==i) numerator = dipole;
      }
      meRatio *= numerator/dipoleSum;
   
      // calculate jacobian
      Energy2 denom = (mb_-particleMomenta[3].e())*particleMomenta[2].vect().mag() -
		       particleMomenta[2].e()*particleMomenta[3].z(); 
      InvEnergy2  J  = (particleMomenta[2].vect().mag2())/(lambda*denom);     
      // calculate weight
      weight[j] = meRatio*fabs(sqr(pT)*J)*coupling()->ratio(pT*pT)/C/Constants::twopi; 
    }
    // accept point if weight > R
    if (weight[0] + weight[1] > UseRandom::rnd()) {
      if (weight[0] > (weight[0] + weight[1])*UseRandom::rnd()) {
	particleMomenta[1].setE( (mb_/2.)*xe  [0]);
	particleMomenta[1].setZ( (mb_/2.)*xe_z[0]);
	particleMomenta[2].setE( (mb_/2.)*xs  [0]);
	particleMomenta[2].setZ(-(mb_/2.)*sqrt(sqr(xs[0])-4.*s2_));
      }
      else {
	particleMomenta[1].setE( (mb_/2.)*xe  [1]);
	particleMomenta[1].setZ( (mb_/2.)*xe_z[1]);
	particleMomenta[2].setE( (mb_/2.)*xs  [1]);
	particleMomenta[2].setZ(-(mb_/2.)*sqrt(sqr(xs[1])-4.*s2_));
      }
      pT_ = pT;
      break;   
    }
    // if there's no branching lower the pT
    pTmax = pT; 
    
  }
  return particleMomenta;
}

double GeneralTwoBodyDecayer::matrixElementRatio(const Particle & inpart, 
						 const ParticleVector & decay2,
						 const ParticleVector & decay3, 
						 MEOption meopt) {
  // calculate R/B
  double B = me2        (0, inpart, decay2, meopt);    
  double R = threeBodyME(0, inpart, decay3, meopt);
  return R/B;
  
}

bool GeneralTwoBodyDecayer::calcMomenta(int j, Energy pT, double y, double phi,
					double& xg, double& xs, double& xe, double& xe_z,
					vector<Lorentz5Momentum>& particleMomenta){
  
  // calculate xg
  xg = 2.*pT*cosh(y) / mb_;
  if (xg>(1. - sqr(e_ + s_)) || xg<0.) return false;

  // calculate the two values of xs
  double xT  = 2.*pT / mb_;
  double A   = 4.-4.*xg+sqr(xT);
  double B   = 4.*(3.*xg-2.+2.*e2_-2.*s2_-sqr(xg)-xg*e2_+xg*s2_);
  double L   = 1.+sqr(s2_)+sqr(e2_)-2.*s2_-2.*e2_-2.*s2_*e2_;
  double det = 16.*( -L*sqr(xT)+pow(xT,4)*s2_+2.*xg*sqr(xT)*(1.-s2_-e2_)+ 
		      L*sqr(xg)-sqr(xg*xT)*(1.+s2_)+pow(xg,4)+ 
		      2.*pow(xg,3)*(-1.+s2_+e2_) );

  if (det<0.) return false;
  if (j==0) xs = (-B+sqrt(det))/(2.*A);
  if (j==1) xs = (-B-sqrt(det))/(2.*A);  
  // check value of xs is physical
  if (xs>(1.+s2_-e2_) || xs<2.*s_) return false;

  // calculate xe
  xe = 2.-xs-xg;     
  // check value of xe is physical
  if (xe>(1.+e2_-s2_) || xe<2.*e_) return false;       

  // calculate xe_z
  double root1 = sqrt(max(0.,sqr(xs)-4.*s2_)), root2 = sqrt(max(0.,sqr(xe)-4.*e2_-sqr(xT)));
  double epsilon_p =  -root1+xT*sinh(y)+root2;
  double epsilon_m =  -root1+xT*sinh(y)-root2;

  // find direction of emitter
  if      (fabs(epsilon_p) < 1.e-10) xe_z =  sqrt(sqr(xe)-4.*e2_-sqr(xT));
  else if (fabs(epsilon_m) < 1.e-10) xe_z = -sqrt(sqr(xe)-4.*e2_-sqr(xT));
  else return false;

  // check the emitter is on shell
  if (fabs((sqr(xe)-sqr(xT)-sqr(xe_z)-4.*e2_))>1.e-10) return false;

  // calculate 4 momenta
  particleMomenta[0].setE   ( mb_);
  particleMomenta[0].setX   ( ZERO);
  particleMomenta[0].setY   ( ZERO);
  particleMomenta[0].setZ   ( ZERO);
  particleMomenta[0].setMass( mb_);

  particleMomenta[1].setE   ( mb_*xe/2.);
  particleMomenta[1].setX   (-pT*cos(phi));
  particleMomenta[1].setY   (-pT*sin(phi));
  particleMomenta[1].setZ   ( mb_*xe_z/2.);
  particleMomenta[1].setMass( mb_*e_);

  particleMomenta[2].setE   ( mb_*xs/2.);
  particleMomenta[2].setX   ( ZERO);
  particleMomenta[2].setY   ( ZERO);
  particleMomenta[2].setZ   (-mb_*sqrt(sqr(xs)-4.*s2_)/2.);
  particleMomenta[2].setMass( mb_*s_);

  particleMomenta[3].setE   ( pT*cosh(y));
  particleMomenta[3].setX   ( pT*cos(phi));
  particleMomenta[3].setY   ( pT*sin(phi));
  particleMomenta[3].setZ   ( pT*sinh(y));
  particleMomenta[3].setMass( ZERO);
 
  return true;
}


bool GeneralTwoBodyDecayer::psCheck(const double xg, const double xs) {
  
  // check is point is in allowed region of phase space
  double xe_star = (1.-s2_+e2_-xg)/sqrt(1.-xg);
  double xg_star = xg/sqrt(1.-xg);

  if ((sqr(xe_star)-4.*e2_) < 1e-10) return false;
  double xs_max = (4.+4.*s2_-sqr(xe_star+xg_star)+ 
		   sqr(sqrt(sqr(xe_star)-4.*e2_)+xg_star))/ 4.;
  double xs_min = (4.+4.*s2_-sqr(xe_star+xg_star)+ 
		   sqr(sqrt(sqr(xe_star)-4.*e2_)-xg_star))/ 4.;

  if (xs < xs_min || xs > xs_max) return false;

  return true;
}


double GeneralTwoBodyDecayer::colourCoeff(const PDT::Colour emitter,
					  const PDT::Colour spectator,
					  const PDT::Colour other){

  // calculate the colour factor of the dipole
  double numerator=1.;
  double denominator=1.;
  if (emitter!=PDT::Colour0 && spectator!=PDT::Colour0 && other!=PDT::Colour0){
    if      (emitter  ==PDT::Colour3 || emitter  ==PDT::Colour3bar) numerator=-4./3;
    else if (emitter  ==PDT::Colour8)                               numerator=-3. ;
    denominator=-1.*numerator;
    if      (spectator==PDT::Colour3 || spectator==PDT::Colour3bar) numerator-=4./3;
    else if (spectator==PDT::Colour8)                               numerator-=3. ;
    if      (other    ==PDT::Colour3 || other    ==PDT::Colour3bar) numerator+=4./3;
    else if (other    ==PDT::Colour8)                               numerator+=3. ;
    numerator*=(-1./2.);				  
  }

  if      (emitter==PDT::Colour3 || emitter==  PDT::Colour3bar) numerator*=4./3.;
  else if (emitter==PDT::Colour8 && spectator!=PDT::Colour8)    numerator*=3.;
  else if (emitter==PDT::Colour8 && spectator==PDT::Colour8)    numerator*=6.;
  
  return (numerator/denominator);
}


InvEnergy2 GeneralTwoBodyDecayer::calculateDipole(const dipoleType & dipoleId,  
						  const Particle & inpart,
						  const ParticleVector & decay3, 
						  const dipoleType & emittingDipole){
  // calculate dipole for decay b->ac
  InvEnergy2 dipole = ZERO;
  double xe = 2.*decay3[0]->momentum().e()/mb_;
  double xs = 2.*decay3[1]->momentum().e()/mb_;
  double xg = 2.*decay3[2]->momentum().e()/mb_;
  double coeff = 8.*Constants::pi*coupling()->value(mb_*mb_); 

  // radiation from b with initial-final connection 
  if (dipoleId==IFba || dipoleId==IFbc){
    dipole  = -2./sqr(mb_*xg);
    dipole *= colourCoeff(inpart.dataPtr()->iColour(), decay3[0]->dataPtr()->iColour(),
			  decay3[1]->dataPtr()->iColour());
  }

  // radiation from a/c with initial-final connection
  else if ((dipoleId==IFa && 
	    (emittingDipole==IFba || emittingDipole==IFa || emittingDipole==FFa)) || 
	   (dipoleId==IFc && 
	    (emittingDipole==IFbc || emittingDipole==IFc || emittingDipole==FFc))){
    double z  = 1. - xg/(1.-s2_+e2_);    
    dipole = (-2.*e2_/sqr(1.-xs+s2_-e2_)/sqr(mb_) + (1./(1.-xs+s2_-e2_)/sqr(mb_))*
	      (2./(1.-z)-dipoleSpinFactor(decay3[0],z)));

    dipole *= colourCoeff(decay3[0]->dataPtr()->iColour(),inpart.dataPtr()->iColour(), 
			  decay3[1]->dataPtr()->iColour());
  }
  else if (dipoleId==IFa || dipoleId==IFc){
    double z  = 1. - xg/(1.-e2_+s2_);
    dipole = (-2.*s2_/sqr(1.-xe+e2_-s2_)/sqr(mb_)+(1./(1.-xe+e2_-s2_)/sqr(mb_))*
	      (2./(1.-z)-dipoleSpinFactor(decay3[1],z)));
    dipole *= colourCoeff(decay3[1]->dataPtr()->iColour(),inpart.dataPtr()->iColour(), 
			  decay3[0]->dataPtr()->iColour());  
  }
  // radiation from a/c with final-final connection
  else if ((dipoleId==FFa && 
	    (emittingDipole==IFba || emittingDipole==IFa || emittingDipole==FFa)) || 
	   (dipoleId==FFc && 
	    (emittingDipole==IFbc || emittingDipole==IFc || emittingDipole==FFc))){
    double z = 1. + ((xe-1.+s2_-e2_)/(xs-2.*s2_));
    dipole = (1./(1.-xs+s2_-e2_)/sqr(mb_))*((2./(1.-z))-dipoleSpinFactor(decay3[0],z)-
					     (2.*e2_/(1.+s2_-e2_-xs)) );
    dipole *= colourCoeff(decay3[0]->dataPtr()->iColour(), 
			  decay3[1]->dataPtr()->iColour(),
			  inpart.dataPtr()->iColour());
  }
  else if (dipoleId==FFa || dipoleId==FFc) { 
    double z = 1. + ((xs-1.+e2_-s2_)/(xe-2.*e2_));
    dipole = (1./(1.-xe+e2_-s2_)/sqr(mb_))*((2./(1.-z))-dipoleSpinFactor(decay3[1],z)-
					     (2.*s2_/(1.+e2_-s2_-xe)) );
    dipole *= colourCoeff(decay3[1]->dataPtr()->iColour(), 
			  decay3[0]->dataPtr()->iColour(),
			  inpart.dataPtr()->iColour());
  }

  dipole *= coeff;
  return dipole;
}

double GeneralTwoBodyDecayer::dipoleSpinFactor(const PPtr & emitter, double z){
  // calculate the spin dependent component of the dipole  
  if      (emitter->dataPtr()->iSpin()==PDT::Spin0)
    return 2.;
  else if (emitter->dataPtr()->iSpin()==PDT::Spin1Half)
    return (1. + z);
  else if (emitter->dataPtr()->iSpin()==PDT::Spin1)
    return (2.*z*(1.-z) - 1./(1.-z) + 1./z -2.);
  return 0.;
}


const vector<DVector> & GeneralTwoBodyDecayer::getColourFactors(const Particle & inpart, 
								const ParticleVector & decay,
								unsigned int & nflow){  
  // calculate the colour factors for the three-body decay
  vector<int> sing,trip,atrip,oct;
  for(unsigned int it=0;it<decay.size();++it) {
    if     (decay[it]->dataPtr()->iColour() == PDT::Colour0    ) sing. push_back(it);
    else if(decay[it]->dataPtr()->iColour() == PDT::Colour3    ) trip. push_back(it);
    else if(decay[it]->dataPtr()->iColour() == PDT::Colour3bar ) atrip.push_back(it);
    else if(decay[it]->dataPtr()->iColour() == PDT::Colour8    ) oct.  push_back(it);
  }
  // require at least one gluon
  assert(oct.size()>=1);

  // identical particle symmetry factor
  double symFactor=1.;
  if (( sing.size()==2 && decay[ sing[0]]->id()==decay[ sing[1]]->id()) ||
      ( trip.size()==2 && decay[ trip[0]]->id()==decay[ trip[1]]->id()) ||
      (atrip.size()==2 && decay[atrip[0]]->id()==decay[atrip[1]]->id()) ||
      (  oct.size()==2 && decay[  oct[0]]->id()==decay[  oct[1]]->id()))
    symFactor/=2.;
  else if (oct.size()==3 && 
	   decay[oct[0]]->id()==decay[oct[1]]->id() &&
	   decay[oct[0]]->id()==decay[oct[2]]->id())
    symFactor/=6.;
  
  colour_ = vector<DVector>(1,DVector(1,symFactor*1.));
  
  // decaying colour singlet   
  if(inpart.dataPtr()->iColour() == PDT::Colour0) {
    if(trip.size()==1 && atrip.size()==1 && oct.size()==1) {
      nflow = 1;
      colour_ = vector<DVector>(1,DVector(1,symFactor*4.));
    }
    else if (oct.size()==3){
      nflow = 1.;
      colour_ = vector<DVector>(1,DVector(1,symFactor*24.));
    }
    else
      throw Exception() << "Unknown colour for the outgoing particles"
			<< " for decay colour scalar particle in "
			<< "GeneralTwoBodyDecayer::getColourFactors() for "
			<< inpart.   dataPtr()->PDGName() << " -> "
			<< decay[0]->dataPtr()->PDGName() << " " 
			<< decay[1]->dataPtr()->PDGName() << " "  
			<< decay[2]->dataPtr()->PDGName()
			<< Exception::runerror;
  }
  // decaying colour triplet
  else if(inpart.dataPtr()->iColour() == PDT::Colour3) {
    if(trip.size()==1 && sing.size()==1 && oct.size()==1) {
      nflow = 1;
      colour_ = vector<DVector>(1,DVector(1,symFactor*4./3.));
    }
    else if(trip.size()==1 && oct.size()==2) {
      nflow = 2;
      colour_.clear();
      colour_.resize(2,DVector(2,0.));
      colour_[0][0] =  symFactor*16./9.; colour_[0][1] = -symFactor*2./9.;
      colour_[1][0] = -symFactor*2./9.;  colour_[1][1] =  symFactor*16./9.;
    }
    else
      throw Exception() << "Unknown colour for the outgoing particles"
			<< " for decay colour triplet particle in "
			<< "GeneralTwoBodyDecayer::getColourFactors() for "
			<< inpart.   dataPtr()->PDGName() << " -> "
			<< decay[0]->dataPtr()->PDGName() << " " 
			<< decay[1]->dataPtr()->PDGName() << " "  
			<< decay[2]->dataPtr()->PDGName()
			<< Exception::runerror;
  }
  // decaying colour anti-triplet
  else if(inpart.dataPtr()->iColour() == PDT::Colour3bar) {
    if(atrip.size()==1 && sing.size()==1 && oct.size()==1) {
      nflow = 1;
      colour_ = vector<DVector>(1,DVector(1,symFactor*4./3.));
    }
    else if(atrip.size()==1 && oct.size()==2){
      nflow = 2;
      colour_.clear();
      colour_ .resize(2,DVector(2,0.));
      colour_[0][0] =  symFactor*16./9.; colour_[0][1] = -symFactor*2./9.;
      colour_[1][0] = -symFactor*2./9.;  colour_[1][1] =  symFactor*16./9.;
    }
    else
      throw Exception() << "Unknown colour for the outgoing particles"
			<< " for decay colour anti-triplet particle in "
			<< "GeneralTwoBodyDecayer::getColourFactors() for "
			<< inpart.   dataPtr()->PDGName() << " -> "
			<< decay[0]->dataPtr()->PDGName() << " " 
			<< decay[1]->dataPtr()->PDGName() << " "  
			<< decay[2]->dataPtr()->PDGName()
			<< Exception::runerror;
  }
  // decaying colour octet
  else if(inpart.dataPtr()->iColour() == PDT::Colour8) {
    if(oct.size()==1 && trip.size()==1 && atrip.size()==1) {
      nflow = 2;
      colour_.clear();
      colour_.resize(2,DVector(2,0.));
      colour_[0][0] =  symFactor*2./3. ; colour_[0][1] = -symFactor*1./12.;
      colour_[1][0] = -symFactor*1./12.; colour_[1][1] =  symFactor*2./3. ;
    }
    else if (oct.size()==2 && sing.size()==1){
      nflow = 1;
      colour_ = vector<DVector>(1,DVector(1,symFactor*3.));
    }
    else
      throw Exception() << "Unknown colour for the outgoing particles"
			<< " for decay colour octet particle in "
			<< "GeneralTwoBodyDecayer::getColourFactors() for "
			<< inpart.   dataPtr()->PDGName() << " -> "
			<< decay[0]->dataPtr()->PDGName() << " " 
			<< decay[1]->dataPtr()->PDGName() << " "  
			<< decay[2]->dataPtr()->PDGName()
			<< Exception::runerror;
  }
  else
    throw Exception() << "Unknown colour for the decaying particle in "
  		      << "GeneralTwoBodyDecayer::getColourFactors() for "
  		      << inpart.   dataPtr()->PDGName() << " -> "
		      << decay[0]->dataPtr()->PDGName() << " " 
		      << decay[1]->dataPtr()->PDGName() << " "  
		      << decay[2]->dataPtr()->PDGName() 
  		      << Exception::runerror;
  return colour_;
}


bool GeneralTwoBodyDecayer::identifyDipoles(vector<dipoleType>  & dipoles,
					    PPtr & aProgenitor,
					    PPtr & bProgenitor,
					    PPtr & cProgenitor) const {

  PDT::Colour bColour = bProgenitor->dataPtr()->iColour();
  PDT::Colour cColour = cProgenitor->dataPtr()->iColour();
  PDT::Colour aColour = aProgenitor->dataPtr()->iColour();

  // decaying colour singlet
  if    (bColour==PDT::Colour0 ) {
    if ((cColour==PDT::Colour3    && aColour==PDT::Colour3bar) ||
	(cColour==PDT::Colour3bar && aColour==PDT::Colour3)    ||
	(cColour==PDT::Colour8    && aColour==PDT::Colour8)){
      dipoles.push_back(FFa);
      dipoles.push_back(FFc);
    }
  }
  // decaying colour triplet
  else if (bColour==PDT::Colour3 ) {
    if (cColour==PDT::Colour3 && aColour==PDT::Colour0){
      dipoles.push_back(IFbc);
      dipoles.push_back(IFc );
    }
    else if (cColour==PDT::Colour0 && aColour==PDT::Colour3){
      dipoles.push_back(IFba);
      dipoles.push_back(IFa );
    }
    else if (cColour==PDT::Colour8 && aColour==PDT::Colour3){
      dipoles.push_back(IFbc);
      dipoles.push_back(IFc );
      dipoles.push_back(FFc );
      dipoles.push_back(FFa );
    }
    else if (cColour==PDT::Colour3 && aColour==PDT::Colour8){
      dipoles.push_back(IFba);
      dipoles.push_back(IFa );
      dipoles.push_back(FFc );
      dipoles.push_back(FFa );
    }
  }
  // decaying colour anti-triplet 
  else if (bColour==PDT::Colour3bar) {
    if ((cColour==PDT::Colour3bar && aColour==PDT::Colour0)){
      dipoles.push_back(IFbc);
      dipoles.push_back(IFc );
    }
    else if ((cColour==PDT::Colour0 && aColour==PDT::Colour3bar)){
      dipoles.push_back(IFba);
      dipoles.push_back(IFa );      
    }
    else if (cColour==PDT::Colour8 && aColour==PDT::Colour3bar){
      dipoles.push_back(IFbc);
      dipoles.push_back(IFc );
      dipoles.push_back(FFc );
      dipoles.push_back(FFa );
    }
    else if (cColour==PDT::Colour3bar && aColour==PDT::Colour8){
      dipoles.push_back(IFba);
      dipoles.push_back(IFa );
      dipoles.push_back(FFc );
      dipoles.push_back(FFa );
    }
  }
  // decaying colour octet
  else if (bColour==PDT::Colour8){
    if ((cColour==PDT::Colour3    && aColour==PDT::Colour3bar) ||
	(cColour==PDT::Colour3bar && aColour==PDT::Colour3)){
      dipoles.push_back(IFba);
      dipoles.push_back(IFbc);
      dipoles.push_back(IFa);
      dipoles.push_back(IFc);
    }
    else if (cColour==PDT::Colour8 && aColour==PDT::Colour0){
      dipoles.push_back(IFbc);
      dipoles.push_back(IFc);
    }
    else if (cColour==PDT::Colour0 && aColour==PDT::Colour8){
      dipoles.push_back(IFba);
      dipoles.push_back(IFa);
    }
  }
  // check colour structure is allowed
  return !dipoles.empty();
}

const GeneralTwoBodyDecayer::CFlow & 
GeneralTwoBodyDecayer::colourFlows(const Particle & inpart,
				   const ParticleVector & decay) {

  // static initialization of commonly used colour structures
  static const CFlow init = CFlow(3, CFlowPairVec(1, make_pair(0, 1.)));
  static CFlow tripflow = init;
  static CFlow atripflow = init;
  static CFlow octflow = init;
  static const CFlow fpflow = CFlow(4, CFlowPairVec(1, make_pair(0, 1.)));

  static bool initialized = false;

  if (! initialized) {
    tripflow[2].resize(2, make_pair(0,1.));
    tripflow[2][0] = make_pair(0, 1.);
    tripflow[2][1] = make_pair(1,-1.);
    tripflow[1][0] = make_pair(1, 1.);
    
    atripflow[1].resize(2, make_pair(0,1.));
    atripflow[1][0] = make_pair(0, 1.);
    atripflow[1][1] = make_pair(1,-1.);
    atripflow[2][0] = make_pair(1, 1.);
    
    octflow[0].resize(2, make_pair(0,1.));   
    octflow[0][0] = make_pair(0,-1.);
    octflow[0][1] = make_pair(1, 1.);
    octflow[2][0] = make_pair(1, 1.);

    initialized = true;
  }
  

  // main function body
  int sing=0,trip=0,atrip=0,oct=0;
  for (size_t it=0; it<decay.size(); ++it) {
    switch ( decay[it]->dataPtr()->iColour() ) {
        case PDT::Colour0:     ++sing; break;
        case PDT::Colour3:     ++trip; break;
        case PDT::Colour3bar: ++atrip; break;
        case PDT::Colour8:      ++oct; break;
        /// @todo: handle these better
        case PDT::ColourUndefined:     break;
        case PDT::Coloured:            break;
        case PDT::Colour6:             break;
        case PDT::Colour6bar:          break;
    }
  }

  // require a gluon
  assert(oct>=1);

  const CFlow * retval = 0;
  bool inconsistent4PV = true;
  // decaying colour triplet
  if(inpart.dataPtr()->iColour() == PDT::Colour3 &&
     trip==1 && oct==2) {
    retval = &tripflow;
  }
  // decaying colour anti-triplet
  else if(inpart.dataPtr()->iColour() == PDT::Colour3bar &&
	  atrip==1 && oct==2){
    retval = &atripflow;
  }
  // decaying colour octet
  else if(inpart.dataPtr()->iColour() == PDT::Colour8 &&
          oct==1 && trip==1 && atrip==1) {
    retval = &octflow;
  }  
  else {
    inconsistent4PV = false;
    retval = &init;  
  }
 
  // if a 4 point vertex exists, add a colour flow for it
  if ( fourPointVertex_ ) {
    if ( inconsistent4PV )   
      throw Exception() << "Unknown colour flows for 4 point vertex in "
			<< "GeneralTwoBodyDecayer::colourFlows()"
			<< Exception::runerror;
    else {
      retval = &fpflow;
    }
  }

  return *retval;
}

void GeneralTwoBodyDecayer::getColourLines(RealEmissionProcessPtr real) {
  // extract the particles
  vector<PPtr> branchingPart;
  branchingPart.push_back(real->incoming()[0]);
  for(unsigned int ix=0;ix<real->outgoing().size();++ix)
    branchingPart.push_back(real->outgoing()[ix]);

  vector<unsigned int> sing,trip,atrip,oct;
  for (size_t ib=0;ib<branchingPart.size()-1;++ib) {
    if     (branchingPart[ib]->dataPtr()->iColour()==PDT::Colour0   ) sing. push_back(ib);
    else if(branchingPart[ib]->dataPtr()->iColour()==PDT::Colour3   ) trip. push_back(ib);
    else if(branchingPart[ib]->dataPtr()->iColour()==PDT::Colour3bar) atrip.push_back(ib);
    else if(branchingPart[ib]->dataPtr()->iColour()==PDT::Colour8   ) oct.  push_back(ib);
  }
  // decaying colour singlet
  if (branchingPart[0]->dataPtr()->iColour()==PDT::Colour0) {
    // 0 -> 3 3bar
    if (trip.size()==1 && atrip.size()==1) {
      branchingPart[atrip[0]]->colourConnect(branchingPart[   3   ]);
      branchingPart[    3   ]->colourConnect(branchingPart[trip[0]]);
    }
    // 0 -> 8 8
    else if (oct.size()==2 ) {
      bool col = UseRandom::rndbool();
      branchingPart[oct[0]]->colourConnect(branchingPart[   3  ],col);
      branchingPart[   3  ]->colourConnect(branchingPart[oct[1]],col);
      branchingPart[oct[1]]->colourConnect(branchingPart[oct[0]],col);
    }
    else 
      assert(false);
  }
  // decaying colour triplet
  else if (branchingPart[0]->dataPtr()->iColour()==PDT::Colour3 ){
    // 3 -> 3 0
    if (trip.size()==2 && sing.size()==1) {
      branchingPart[3]->incomingColour(branchingPart[trip[0]]);
      branchingPart[3]-> colourConnect(branchingPart[trip[1]]);
    }
    // 3 -> 3 8
    else if (trip.size()==2 && oct.size()==1) {
      // 8 emit incoming partner
      if(real->emitter()==oct[0]&&real->spectator()==0) {
	branchingPart[  3   ]->incomingColour(branchingPart[trip[0]]);
	branchingPart[  3   ]-> colourConnect(branchingPart[oct[0] ]);
	branchingPart[oct[0]]-> colourConnect(branchingPart[trip[1]]);
      }
      // 8 emit final spectator or vice veras
      else {
	branchingPart[oct[0]]->incomingColour(branchingPart[trip[0]]);
	branchingPart[oct[0]]-> colourConnect(branchingPart[   3   ]);
	branchingPart[   3  ]-> colourConnect(branchingPart[trip[1]]);
      }
    }
    else
      assert(false);
  }
  // decaying colour anti-triplet
  else if (branchingPart[0]->dataPtr()->iColour()==PDT::Colour3bar) {
    // 3bar -> 3bar 0
    if (atrip.size()==2 && sing.size()==1) {      
      branchingPart[3]->incomingColour(branchingPart[atrip[0]],true);
      branchingPart[3]-> colourConnect(branchingPart[atrip[1]],true);
    }
    // 3 -> 3 8
    else if (atrip.size()==2 && oct.size()==1){
      // 8 emit incoming partner
      if(real->emitter()==oct[0]&&real->spectator()==0) {
	branchingPart[   3  ]->incomingColour(branchingPart[atrip[0]],true);
	branchingPart[   3  ]-> colourConnect(branchingPart[oct[0]  ],true);
	branchingPart[oct[0]]-> colourConnect(branchingPart[atrip[1]],true);
      }
      // 8 emit final spectator or vice veras
      else {
	branchingPart[oct[0]]->incomingColour(branchingPart[atrip[0]],true);
	branchingPart[oct[0]]-> colourConnect(branchingPart[   3    ],true);
	branchingPart[3]-> colourConnect(branchingPart[atrip[1]]     ,true);
      }
    }
    else
      assert(false);
  }
  // decaying colour octet
  else if(branchingPart[0]->dataPtr()->iColour()==PDT::Colour8 ) {
    // 8 -> 3 3bar
    if (trip.size()==1 && atrip.size()==1) {
      // 3 emits
      if(trip[0]==real->emitter()) {
	branchingPart[3]       ->incomingColour(branchingPart[oct[0]] );
	branchingPart[3]       -> colourConnect(branchingPart[trip[0]]);
	branchingPart[atrip[0]]->incomingColour(branchingPart[oct[0]],true);
      }
      // 3bar emits
      else {
	branchingPart[3]       ->incomingColour(branchingPart[oct[0]]  ,true);
	branchingPart[3]       -> colourConnect(branchingPart[atrip[0]],true);
	branchingPart[trip[0]]->incomingColour(branchingPart[oct[0]]  );
      }
    }
    // 8 -> 8 0 
    else if (sing.size()==1 && oct.size()==2) {
      bool col = UseRandom::rndbool();
      branchingPart[   3  ]->colourConnect (branchingPart[oct[1]], col);
      branchingPart[   3  ]->incomingColour(branchingPart[oct[0]], col);
      branchingPart[oct[1]]->incomingColour(branchingPart[oct[0]],!col);
    }
    else
      assert(false);
  }
}
