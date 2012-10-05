// -*- C++ -*-
//
// GeneralTwoBodyDecayer.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
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
     << coupling_ << incomingVertex_ << outgoingVertices_;
}

void GeneralTwoBodyDecayer::persistentInput(PersistentIStream & is, int) {
  is >> vertex_ >> _incoming >> _outgoing >> _maxweight >> iunit(pTmin_,GeV) 
     >> coupling_ >> incomingVertex_ >> outgoingVertices_;
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
    throw Exception() << "Unknown colour for the decaying particle in "
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
					 const vector<VertexBasePtr> & outV) {
  _incoming=incoming;
  _outgoing.clear();
  _outgoing.push_back(outgoing.first );
  _outgoing.push_back(outgoing.second);
  vertex_ = vertex;
  incomingVertex_ = inV;
  outgoingVertices_ = outV;
}

HardTreePtr GeneralTwoBodyDecayer::generateHardest(ShowerTreePtr tree) {
  // search for the coloured particles
  bool colouredParticles=false;
  vector<ShowerProgenitorPtr> Progenitors = tree->extractProgenitors();
  for (unsigned int it=0; it<Progenitors.size(); ++it){
    if (Progenitors[it]->progenitor()->dataPtr()->coloured()){
      colouredParticles=true;
      break;
    }
  }
  // if no coloured particles return
  if ( !colouredParticles ) {
    for (unsigned int it=0; it<Progenitors.size(); ++it){
      Progenitors[it]->maximumpT(pTmin_);
    }
    return HardTreePtr();
  }
  // check exactly two outgoing particles
  if (tree->outgoingLines().size()!=2) {
    throw Exception()
      << "Number of outgoing particles is not equal to 2 in "
      << "GeneralTwoBodyDecayer::generateHardest()" 
      << Exception::runerror;
  }
  // for decay b -> a c 
  // set progenitors
  // outgoing
  ShowerProgenitorPtr 
    cProgenitor = tree->outgoingLines(). begin()->first,
    aProgenitor = tree->outgoingLines().rbegin()->first;
  // Get the decaying particle
  ShowerProgenitorPtr bProgenitor = tree->incomingLines().begin()->first;
  // Identify the dipoles
  vector<dipoleType> dipoles;
  identifyDipoles(dipoles,aProgenitor,bProgenitor,cProgenitor);

  Energy trialpT = pTmin_;
  LorentzRotation eventFrame;
  vector<Lorentz5Momentum> momenta;
  vector<Lorentz5Momentum> trialMomenta(4);
  ShowerProgenitorPtr finalEmitter, finalSpectator;
  ShowerProgenitorPtr trialEmitter, trialSpectator;

  for (int i=0; i<dipoles.size(); ++i){
    //assign emitter and spectator depending on current dipole
    if (dipoles[i]==FFc || dipoles[i]==IFc || dipoles[i]==IFbc){
      trialEmitter   = cProgenitor;
      trialSpectator = aProgenitor;
    }
    else if (dipoles[i]==FFa || dipoles[i]==IFa || dipoles[i]==IFba){
      trialEmitter   = aProgenitor;
      trialSpectator = cProgenitor;
    }

    // find rotation from lab to frame with a along -z
    LorentzRotation trialEventFrame = ( bProgenitor->progenitor()->momentum().findBoostToCM() );
    Lorentz5Momentum pspectator = trialEventFrame*trialSpectator->progenitor()->momentum();
    trialEventFrame.rotateZ( -pspectator.phi() );
    trialEventFrame.rotateY( -pspectator.theta() - Constants::pi );

    //invert it
    trialEventFrame.invert();

    trialMomenta.clear();
    pT_ = pTmin_;
    trialMomenta = hardMomenta(bProgenitor, trialEmitter, trialSpectator, dipoles, i);

    if(pT_>trialpT){
      trialpT        = pT_;
      momenta        = trialMomenta;
      eventFrame     = trialEventFrame;
      finalEmitter   = trialEmitter;
      finalSpectator = trialSpectator;
    }
  }  
  pT_=trialpT;

  // if no emission return
  if(momenta.empty()) {
    bProgenitor   ->maximumpT(pTmin_);
    aProgenitor   ->maximumpT(pTmin_);
    cProgenitor   ->maximumpT(pTmin_);
    return HardTreePtr();
  }

  // rotate momenta back to the lab
  for(unsigned int ix=0;ix<momenta.size();++ix) {
    momenta[ix] *= eventFrame;
  }
  
  bProgenitor   ->maximumpT(pT_);
  finalEmitter  ->maximumpT(pT_);
  finalSpectator->maximumpT(pT_);

  // get ParticleData objects
  tcPDPtr b = bProgenitor   ->progenitor()->dataPtr();
  tcPDPtr e = finalEmitter  ->progenitor()->dataPtr();
  tcPDPtr s = finalSpectator->progenitor()->dataPtr();
  tcPDPtr gluon  = getParticleData(ParticleID::g);

  // create new ShowerParticles
  ShowerParticlePtr emitter  (new_ptr(ShowerParticle(e,     true )));
  ShowerParticlePtr spectator(new_ptr(ShowerParticle(s,     true )));
  ShowerParticlePtr gauge    (new_ptr(ShowerParticle(gluon, true )));
  ShowerParticlePtr incoming (new_ptr(ShowerParticle(b,     false)));
  ShowerParticlePtr parent   (new_ptr(ShowerParticle(e,     true )));

  // set momenta
  emitter  ->set5Momentum(momenta[1]); 
  spectator->set5Momentum(momenta[2]);  
  gauge    ->set5Momentum(momenta[3]); 
  incoming ->set5Momentum(bProgenitor->progenitor()->momentum());  
  Lorentz5Momentum parentMomentum(momenta[1]+momenta[3]);
  parentMomentum.rescaleMass();
  parent->set5Momentum(parentMomentum);

  // Create the vectors of HardBranchings to create the HardTree:
  vector<HardBranchingPtr> spaceBranchings,allBranchings;
  // Incoming particle b
  spaceBranchings.push_back(new_ptr(HardBranching(incoming,SudakovPtr(),
						  HardBranchingPtr(),
						  HardBranching::Incoming)));
  // Outgoing particles from hard emission:
  HardBranchingPtr spectatorBranch(new_ptr(HardBranching(spectator,SudakovPtr(),
							 HardBranchingPtr(),
							 HardBranching::Outgoing)));
  HardBranchingPtr emitterBranch(new_ptr(HardBranching(parent,SudakovPtr(),
						       HardBranchingPtr(),
						       HardBranching::Outgoing)));
  emitterBranch->addChild(new_ptr(HardBranching(emitter,SudakovPtr(),
						HardBranchingPtr(),
						HardBranching::Outgoing)));
  emitterBranch->addChild(new_ptr(HardBranching(gauge,SudakovPtr(),
						HardBranchingPtr(),
						HardBranching::Outgoing)));
  allBranchings.push_back(spaceBranchings[0]);
  allBranchings.push_back(emitterBranch);
  allBranchings.push_back(spectatorBranch);
  // Make the HardTree from the HardBranching vectors.
  HardTreePtr hardtree = new_ptr(HardTree(allBranchings,spaceBranchings,
					  ShowerInteraction::QCD));
  
  // Connect the particles with the branchings in the HardTree
  hardtree->connect( bProgenitor   ->progenitor(), spaceBranchings[0] );
  hardtree->connect( finalEmitter  ->progenitor(),   allBranchings[1] );
  hardtree->connect( finalSpectator->progenitor(),   allBranchings[2] );
  // colour flow
  ColinePtr newline=new_ptr(ColourLine());
  for(set<HardBranchingPtr>::const_iterator cit=hardtree->branchings().begin();
      cit!=hardtree->branchings().end();++cit) {

    if((**cit).branchingParticle()->dataPtr()->iColour()==PDT::Colour3)
      newline->addColoured((**cit).branchingParticle());
    else if((**cit).branchingParticle()->dataPtr()->iColour()==PDT::Colour3bar)
      newline->addAntiColoured((**cit).branchingParticle());
  }

  //return the tree
  return hardtree;
}

double GeneralTwoBodyDecayer::threeBodyME(const int , const Particle & inpart,
					  const ParticleVector & decay,MEOption meopt) {
  throw Exception() << "Base class  GeneralTwoBodyDecayer::threeBodyME() "
		    << "called, should have an implementation in the inheriting class"
		    << Exception::runerror;
  return 0.;
}

vector<Lorentz5Momentum>  GeneralTwoBodyDecayer::hardMomenta(const ShowerProgenitorPtr &in, 
							     const ShowerProgenitorPtr &emitter, 
							     const ShowerProgenitorPtr &spectator, 
							     vector<dipoleType> dipoles, int i) {
  double C    = 6.3;
  double ymax = 10.;
  double ymin = -ymax;

  // masses of the particles
  mb_  = in       ->progenitor()->momentum().mass();
  e_   = emitter  ->progenitor()->momentum().mass() / mb_;
  s_   = spectator->progenitor()->momentum().mass() / mb_;
  e2_  = sqr(e_);
  s2_  = sqr(s_);
  
  vector<Lorentz5Momentum> particleMomenta (4);
  Energy2 lambda = sqr(mb_)* sqrt( 1. + sqr(s2_) + sqr(e2_) - 2.*s2_ - 2.*e2_ - 2.*s2_*e2_);    

  //Calculate A
  double A = (ymax - ymin) * C * (coupling_->overestimateValue() / (2.*Constants::pi));
 
  Energy pTmax = mb_* (sqr(1.-s_) - e2_) / (2.*(1.-s_));
  if (pTmax < pTmin_) particleMomenta.clear(); 

  while (pTmax >= pTmin_) {  
    //Generate pT, y and phi values
    Energy pT = pTmax * pow(UseRandom::rnd() , (1./A));  
    if (pT < pTmin_) {particleMomenta.clear(); break;}

    double phi = UseRandom::rnd() * Constants::twopi;
    double y   = ymin + UseRandom::rnd() * (ymax-ymin);

    double weight[2] = {0.,0.};
    double xs[2], xe[2], xe_z[2], xg;
    
    for (unsigned int j=0; j<2; j++) {

      //Check if the momenta are physical
      bool physical = calcMomenta(j, pT, y, phi, xg, xs[j], xe[j], xe_z[j], particleMomenta);
      if (!physical) continue;
      
      //Check if point lies within phase space
      bool inPS = psCheck(xg, xs[j]);
      if (!inPS) continue;
      
      //Calculate the ratio R/B
      //decay products for 3 body decay
      Particle inpart   = Particle(in       ->progenitor()->dataPtr());      
      Particle outparte = Particle(emitter  ->progenitor()->dataPtr()); 
      Particle outparts = Particle(spectator->progenitor()->dataPtr()); 
      Particle g        = Particle(getParticleData(ParticleID::g));
      inpart.  set5Momentum(particleMomenta[0]);
      outparte.set5Momentum(particleMomenta[1]);
      outparts.set5Momentum(particleMomenta[2]);
      g.       set5Momentum(particleMomenta[3]);
      ParticleVector decay3;
      decay3.push_back(&outparte); decay3.push_back(&outparts); decay3.push_back(&g);
      
      //decay products for 2 body decay
      Lorentz5Momentum p1, p2; 
      p1.setE((mb_/2.)*(1.+e2_-s2_)); p1.setX(ZERO); p1.setY(ZERO); p1.setZ( lambda/2./mb_);
      p1.setMass(mb_*e_);    
      p2.setE((mb_/2.)*(1.+s2_-e2_)); p2.setX(ZERO); p2.setY(ZERO); p2.setZ(-lambda/2./mb_);
      p2.setMass(mb_*s_);
      Particle outpart1 = Particle(emitter  ->progenitor()->dataPtr()); outpart1.set5Momentum(p1);
      Particle outpart2 = Particle(spectator->progenitor()->dataPtr()); outpart2.set5Momentum(p2);
      ParticleVector decay2; decay2.push_back(&outpart1); decay2.push_back(&outpart2);

      //calculate matrix element ratio R/B
      double meRatio = matrixElementRatio(inpart,decay2,decay3,Initialize);

      //calculate dipole factor
      InvEnergy2 dipoleSum = 1./GeV/GeV;

      //S->FF
      for (int k=0; k<dipoles.size(); ++k){
	if (dipoles.size()==2){
	  if ((dipoles[k]==FFa && dipoles[i]==FFc) ||
	      (dipoles[k]==FFc && dipoles[i]==FFa)){
	    dipoleSum = abs(calculateDipole(dipoles[i], xe[j], xs[j], xe_z[j]))
	              + abs(calculateDipole(dipoles[k], xs[j], xe[j], xe_z[j]));
	  }
	  else if ((dipoles[i]==IFba && dipoles[k]==IFa ) ||
		   (dipoles[i]==IFbc && dipoles[k]==IFc ) ||
		   (dipoles[i]==IFa  && dipoles[k]==IFba) ||
		   (dipoles[i]==IFc  && dipoles[k]==IFbc)){
	    dipoleSum = abs(calculateDipole(dipoles[i], xe[j], xs[j], xe_z[j]))
	              + abs(calculateDipole(dipoles[k], xe[j], xs[j], xe_z[j]));
	  }	  
	}
	else if (dipoles.size()==4){
	  if (dipoles[k]==IFa){
	    dipoleSum = abs(calculateDipole(IFba, xe[j], xs[j], xe_z[j]))
	              + abs(calculateDipole(IFa , xe[j], xs[j], xe_z[j]))
	              + abs(calculateDipole(FFa , xe[j], xs[j], xe_z[j]))
	              + abs(calculateDipole(FFc , xs[j], xe[j], xe_z[j]));
	  }
	  else if (dipoles[k]==IFc){
	    dipoleSum = abs(calculateDipole(IFbc, xe[j], xs[j], xe_z[j]))
	              + abs(calculateDipole(IFc , xe[j], xs[j], xe_z[j]))
	              + abs(calculateDipole(FFc , xe[j], xs[j], xe_z[j]))
	              + abs(calculateDipole(FFa , xs[j], xe[j], xe_z[j]));
	  }
	}	
      }

      double dipoleFactor = abs(calculateDipole(dipoles[i], xe[j], xs[j], xe_z[j]))/dipoleSum;
      meRatio *= dipoleFactor;
      //Calculate jacobian
      Energy2 denom = (mb_ - particleMomenta[3].e()) * particleMomenta[2].vect().mag() -
		       particleMomenta[2].e() * particleMomenta[3].z(); 
      InvEnergy2  J  = (particleMomenta[2].vect().mag2()) / (lambda * denom);     
      //Calculate weight
      weight[j] = meRatio * fabs(sqr(pT)*J) * coupling_->ratio(pT*pT) / C / Constants::twopi; 
     }

    //ofstream weights;
    //if (weight[0] + weight[1] > 1.){
    //weights.open("weights.top", ios::app);
    //weights << weight[0]+weight[1] << endl;
    //}

    //Accept point if weight > R
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
    //If there's no splitting lower the pT
    pTmax = pT; 
    
  }
  return particleMomenta;
}

double GeneralTwoBodyDecayer::matrixElementRatio(const Particle & inpart, 
						 const ParticleVector & decay2,
						 const ParticleVector & decay3, 
						 MEOption meopt) {
  double B = me2(0, inpart, decay2, meopt);
  double R = threeBodyME(0, inpart, decay3, meopt);
  return R/B;
}

bool GeneralTwoBodyDecayer::calcMomenta(int j, Energy pT, double y, double phi,
   		        double& xg, double& xs, double& xe, double& xe_z,
			vector<Lorentz5Momentum>& particleMomenta){
  
  //Calculate xg
  xg = 2.*pT*cosh(y) / mb_;
  if (xg>(1. - sqr(e_ + s_)) || xg<0.) return false;

  //Calculate xs
  double xT  = 2.*pT / mb_;
  double A   = 4. - 4.*xg + sqr(xT);
  double B   = 4.*(3.*xg - 2. + 2.*e2_ - 2.*s2_ - sqr(xg) - xg*e2_ + xg*s2_);
  double L   = 1. + sqr(s2_) + sqr(e2_) - 2.*s2_ - 2.*e2_ - 2.*s2_*e2_;
  double det = 16.*( -L*sqr(xT) + pow(xT,4)*s2_ + 2.*xg*sqr(xT)*(1.-s2_-e2_) + 
		      L*sqr(xg) - sqr(xg*xT)*(1. + s2_) + pow(xg,4) + 
		      2.*pow(xg,3)*(- 1. + s2_ + e2_) );

  if (det<0.) return false;
  if (j==0) xs = (-B + sqrt(det))/(2.*A);
  if (j==1) xs = (-B - sqrt(det))/(2.*A);  
  if (xs>(1. + s2_ - e2_) || xs<2.*s_) return false;

  //Calculate xe
  xe = 2. - xs - xg;     
  if (xe>(1. + e2_ - s2_) || xe<2.*e_) return false;       

  //Calculate xe_z  
  double epsilon_p =  -sqrt(sqr(xs) - 4.*s2_) + xT*sinh(y) + sqrt(sqr(xe) - 4.*e2_ - sqr(xT));
  double epsilon_m =  -sqrt(sqr(xs) - 4.*s2_) + xT*sinh(y) - sqrt(sqr(xe) - 4.*e2_ - sqr(xT));

  if      (fabs(epsilon_p) < 1.e-10) xe_z =  sqrt(sqr(xe) - 4.*e2_ - sqr(xT));
  else if (fabs(epsilon_m) < 1.e-10) xe_z = -sqrt(sqr(xe) - 4.*e2_ - sqr(xT));
  
  else return false;

  //Check b is on shell
  if (fabs((sqr(xe) - sqr(xT) - sqr(xe_z) - 4.*e2_))>1.e-10) return false;

  //Calculate 4 momenta
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
  particleMomenta[2].setZ   (-mb_*sqrt(sqr(xs) - 4.*s2_)/2.);
  particleMomenta[2].setMass( mb_*s_);

  particleMomenta[3].setE   ( pT*cosh(y));
  particleMomenta[3].setX   ( pT*cos(phi));
  particleMomenta[3].setY   ( pT*sin(phi));
  particleMomenta[3].setZ   ( pT*sinh(y));
  particleMomenta[3].setMass( ZERO);
 
  return true;
}


bool GeneralTwoBodyDecayer::psCheck(double xg, double xs) {
  
  //Check is point is in allowed region of phase space
  double xe_star = (1. - s2_ + e2_ - xg) / sqrt(1. - xg);
  double xg_star = xg / sqrt(1. - xg);

  if ((sqr(xe_star) - 4.*e2_) < 1e-10) return false;
  double xs_max = (4. + 4.*s2_ - sqr(xe_star + xg_star) + 
		   sqr(sqrt(sqr(xe_star) - 4.*e2_) + xg_star)) / 4.;
  double xs_min = (4. + 4.*s2_ - sqr(xe_star + xg_star) + 
		   sqr(sqrt(sqr(xe_star) - 4.*e2_) - xg_star)) / 4.;

  if (xs < xs_min || xs > xs_max) return false;

  return true;
}

InvEnergy2 GeneralTwoBodyDecayer::calculateDipole(dipoleType dipoleId, double xe, double xs,
						  double xe_z){

  InvEnergy2 dipole = 1./GeV/GeV;
  double coeff = 4./3.*8.*Constants::pi*coupling_->value(mb_*mb_);

  if (dipoleId==FFa || dipoleId==FFc){
    double lambda2 = 1. + sqr(s2_) + sqr(e2_) - 2.*s2_ - 2.*e2_ - 2.*s2_*e2_; 
    dipole = (1./(1.-xs+s2_-e2_)/sqr(mb_))*( (2.*(1.-s2_-e2_)/(2.-xs-xe))-
					      sqrt((lambda2)/(sqr(xs)-4.*s2_))*
					    ((xs-2.*s2_)/(1-s2_-e2_))*
					     (2.+ ((xe-1.+s2_-e2_)/(xs-2.*s2_)) + 
					     (2.*e2_/(1.+s2_-e2_-xs))) );
  }
  else if (dipoleId==IFba || dipoleId==IFbc){
    double xg = 2. - xe - xs;
    dipole = -4./ sqr(mb_*xg);
  }
  else if (dipoleId==IFa || dipoleId==IFc){
    double xg = 2. - xe - xs;
    double xT = sqrt(sqr(xe) - sqr(xe_z) - 4.*e2_);
    double z  = 1. - xg/(1. - s2_ + e2_);
    Energy2 pepg = sqr(mb_/2.)*(xg*xe + sqr(xT) - xe_z*sqrt(sqr(xg)-sqr(xT)));
    dipole = - sqr(mb_)*e2_ / sqr(pepg) + (1./pepg)*(2./(1.-z) - 1. - z);
  }

  dipole *= coeff;
  return dipole;
}


const vector<DVector> & GeneralTwoBodyDecayer::getColourFactors(const Particle & inpart, 
								const ParticleVector & decay,
								unsigned int & nflow){  
  vector<int> sing,trip,atrip,oct;
  
  for(unsigned int it=0;it<decay.size();++it) {
    if     (decay[it]->dataPtr()->iColour() == PDT::Colour0    ) sing. push_back(it);
    else if(decay[it]->dataPtr()->iColour() == PDT::Colour3    ) trip. push_back(it);
    else if(decay[it]->dataPtr()->iColour() == PDT::Colour3bar ) atrip.push_back(it);
    else if(decay[it]->dataPtr()->iColour() == PDT::Colour8    ) oct.  push_back(it);
  }
  //require a gluon
  assert(oct.size()>=1);

  colour_ = vector<DVector>(1,DVector(1,1));
  
  //decaying colour singlet   
  if(inpart.dataPtr()->iColour() == PDT::Colour0) {
    if(trip.size()==1 && atrip.size()==1 && oct.size()==1) {
      nflow = 1;
      colour_ = vector<DVector>(1,DVector(1,4.));
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
  //decaying colour triplet
  else if(inpart.dataPtr()->iColour() == PDT::Colour3) {
    if(trip.size()==1 && sing.size()==1 && oct.size()==1) {
      nflow = 1;
      colour_ = vector<DVector>(1,DVector(1,4./3.));
    }
    else if(trip.size()==1 && oct.size()==2) {
      nflow = 2;
      colour_ .resize(2,DVector(2,0.));
      colour_[0][0] = 16./3.; colour_[0][1] = -2./9.;
      colour_[1][0] = -2./9.; colour_[1][1] = 16./3.;
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
  //decaying anti-triplet
  else if(inpart.dataPtr()->iColour() == PDT::Colour3bar) {
    if(atrip.size()==1 && sing.size()==1 && oct.size()==1) {
      nflow = 1;
      colour_ = vector<DVector>(1,DVector(1,4./3.));
    }
    else if(atrip.size()==1 && oct.size()==2){
      nflow = 2;
      colour_ .resize(2,DVector(2,0.));
      colour_[0][0] = 16./9.; colour_[0][1] = -2./9.;
      colour_[1][0] = -2./9.; colour_[1][1] = 16./9.;
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
  //decaying octet
  else if(inpart.dataPtr()->iColour() == PDT::Colour8) {
    if(oct.size()==2 && sing.size() == 1) {
      nflow = 1;
      colour_ = vector<DVector>(1,DVector(1,3.));
    }
    else if(oct.size()==1 && trip.size()==1 && atrip.size()==1) {
      nflow = 2;
      colour_ .resize(2,DVector(2,0.));
      colour_[0][0] =  2./3. ; colour_[0][1] = -1./12.;
      colour_[1][0] = -1./12.; colour_[1][1] =  2./3. ;
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
  //decaying sextet
  else if(inpart.dataPtr()->iColour() == PDT::Colour6) {
    if(oct.size()==1 && trip.size()==2) {
      nflow = 1;
      //colour_ = vector<DVector>(1,DVector(1,1));
    }
    else
      throw Exception() << "Unknown colour for the outgoing particles"
			<< " for decay colour sextet particle in "
			<< "GeneralTwoBodyDecayer::getColourFactors() for "
			<< inpart.   dataPtr()->PDGName() << " -> "
			<< decay[0]->dataPtr()->PDGName() << " " 
			<< decay[1]->dataPtr()->PDGName() << " "  
			<< decay[2]->dataPtr()->PDGName()
			<< Exception::runerror;
  }
  //decaying anti-sextet
  else if(inpart.dataPtr()->iColour() == PDT::Colour6bar) {
    if(oct.size()==1 && atrip.size()==2) {
      nflow = 1;
      //colour_ = vector<DVector>(1,DVector(1,1));
    }
    else
      throw Exception() << "Unknown colour for the outgoing particles"
			<< " for decay colour anti-sextet particle in "
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


void GeneralTwoBodyDecayer::identifyDipoles(vector<dipoleType> & dipoles,
					    ShowerProgenitorPtr & aProgenitor,
					    ShowerProgenitorPtr & bProgenitor,
					    ShowerProgenitorPtr & cProgenitor) const {

  // neutral decaying particle
  if (bProgenitor->progenitor()->dataPtr()->iColour()==PDT::Colour0 ) {
    if ((cProgenitor->progenitor()->dataPtr()->iColour()==PDT::Colour3     &&
	 aProgenitor->progenitor()->dataPtr()->iColour()==PDT::Colour3bar) ||
	(cProgenitor->progenitor()->dataPtr()->iColour()==PDT::Colour3bar  &&
	 aProgenitor->progenitor()->dataPtr()->iColour()==PDT::Colour3)){
      dipoles.push_back(FFa);
      dipoles.push_back(FFc);
    }
  }
  // triplet decaying particle
  else if (bProgenitor->progenitor()->dataPtr()->iColour()==PDT::Colour3 ) {
    if (cProgenitor->progenitor()->dataPtr()->iColour()==PDT::Colour3 &&
	aProgenitor->progenitor()->dataPtr()->iColour()==PDT::Colour0){
      dipoles.push_back(IFbc);
      dipoles.push_back(IFc );
    }
    else if (cProgenitor->progenitor()->dataPtr()->iColour()==PDT::Colour0 &&
	     aProgenitor->progenitor()->dataPtr()->iColour()==PDT::Colour3){
      dipoles.push_back(IFba);
      dipoles.push_back(IFa );
    }
    else if (cProgenitor->progenitor()->dataPtr()->iColour()==PDT::Colour8 &&
	     aProgenitor->progenitor()->dataPtr()->iColour()==PDT::Colour3){
      dipoles.push_back(IFbc);
      dipoles.push_back(IFc );
      dipoles.push_back(FFc );
      dipoles.push_back(FFa );
    }
    else if (cProgenitor->progenitor()->dataPtr()->iColour()==PDT::Colour3 &&
	     aProgenitor->progenitor()->dataPtr()->iColour()==PDT::Colour8){
      dipoles.push_back(IFba);
      dipoles.push_back(IFa );
      dipoles.push_back(FFc );
      dipoles.push_back(FFa );
    }
  }
  // antitriplet decaying particle 
  else if (bProgenitor->progenitor()->dataPtr()->iColour()==PDT::Colour3bar) {
    if ((cProgenitor->progenitor()->dataPtr()->iColour()==PDT::Colour3bar &&
	 aProgenitor->progenitor()->dataPtr()->iColour()==PDT::Colour0)){
      dipoles.push_back(IFbc);
      dipoles.push_back(IFc );
    }
    else if ((cProgenitor->progenitor()->dataPtr()->iColour()==PDT::Colour0 &&
	      aProgenitor->progenitor()->dataPtr()->iColour()==PDT::Colour3bar)){
      dipoles.push_back(IFba);
      dipoles.push_back(IFa );      
    }
    else if (cProgenitor->progenitor()->dataPtr()->iColour()==PDT::Colour8 &&
	     aProgenitor->progenitor()->dataPtr()->iColour()==PDT::Colour3bar){
      dipoles.push_back(IFbc);
      dipoles.push_back(IFc );
      dipoles.push_back(FFc );
      dipoles.push_back(FFa );
    }
    else if (cProgenitor->progenitor()->dataPtr()->iColour()==PDT::Colour3bar &&
	     aProgenitor->progenitor()->dataPtr()->iColour()==PDT::Colour8){
      dipoles.push_back(IFba);
      dipoles.push_back(IFa );
      dipoles.push_back(FFc );
      dipoles.push_back(FFa );
    }
  }
  // check this is allowed
  if (dipoles.size()==0) 
    throw Exception() << "Unknown colour structure in 3 boday decay in "
		      << "GeneralTwoBodyDecayer::generateHardest()"
		      << Exception::runerror;
}
