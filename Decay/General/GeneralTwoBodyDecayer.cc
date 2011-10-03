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
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
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
  assert( _theVertex );
  assert( _incoming && _outgoing.size()==2);
  _theVertex->init();
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
  else
    throw Exception() << "Unknown incoming colour in "
		      << "GeneralTwoBodyDecayer::colourConnections() "
		      << incColour
		      << Exception::runerror; 
}

bool GeneralTwoBodyDecayer::twoBodyMEcode(const DecayMode & dm, int & mecode,
					  double & coupling) const {
  long parent = dm.parent()->id();
  ParticleMSet::const_iterator pit = dm.products().begin();
  long id1 = (*pit)->id();
  ++pit;
  long id2 = (*pit)->id();
  assert(parent == _incoming->id());
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
  os << _theVertex << _incoming << _outgoing << _maxweight;
}

void GeneralTwoBodyDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _theVertex >> _incoming >> _outgoing >> _maxweight;
}

AbstractClassDescription<GeneralTwoBodyDecayer> 
GeneralTwoBodyDecayer::initGeneralTwoBodyDecayer;
// Definition of the static class description member.

void GeneralTwoBodyDecayer::Init() {

  static ClassDocumentation<GeneralTwoBodyDecayer> documentation
    ("This class is designed to be a base class for all 2 body decays"
     "in a general model");
 
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
					 VertexBasePtr vertex) {
  _incoming=incoming;
  _outgoing.clear();
  _outgoing.push_back(outgoing.first );
  _outgoing.push_back(outgoing.second);
  _theVertex = vertex;
}
