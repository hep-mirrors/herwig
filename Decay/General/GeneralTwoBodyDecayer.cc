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
     << coupling_;
}

void GeneralTwoBodyDecayer::persistentInput(PersistentIStream & is, int) {
  is >> vertex_ >> _incoming >> _outgoing >> _maxweight >> iunit(pTmin_,GeV) 
     >> coupling_;
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

  if (tree->outgoingLines().size()!=2){
  
    vector<ShowerProgenitorPtr> check = tree->extractProgenitors();
    cerr << check.size() << endl;
    ShowerProgenitorPtr gProg = check[0];
    ShowerProgenitorPtr hProg = check[1];
    ShowerProgenitorPtr iProg = check[2];
    ShowerProgenitorPtr jProg = check[3];
  
    cerr << gProg->progenitor()->dataPtr()->PDGName() << "\t->\t"
	 << hProg->progenitor()->dataPtr()->PDGName() << "\t"
	 << iProg->progenitor()->dataPtr()->PDGName() << "\t"
	 << jProg->progenitor()->dataPtr()->PDGName() << endl;

    gProg  ->maximumpT(pTmin_);
    hProg  ->maximumpT(pTmin_);
    iProg  ->maximumpT(pTmin_);
    jProg  ->maximumpT(pTmin_);
    return HardTreePtr();
  }

  //for decay b -> a c where a is the colour singlet
  assert(tree->outgoingLines().size()==2);
  
  ShowerProgenitorPtr 
    cProgenitor = tree->outgoingLines(). begin()->first,
    aProgenitor = tree->outgoingLines().rbegin()->first;
  
  if((aProgenitor->progenitor()->dataPtr()->iColour()==PDT::Colour3 || 
      aProgenitor->progenitor()->dataPtr()->iColour()==PDT::Colour3bar) &&
      cProgenitor->progenitor()->dataPtr()->iColour()==PDT::Colour0) 
    swap(cProgenitor,aProgenitor);
  // Get the decaying particle
  ShowerProgenitorPtr bProgenitor = tree->incomingLines().begin()->first;

  //  cerr << bProgenitor->progenitor()->dataPtr()->PDGName() << "\t->\t"
  //   << cProgenitor->progenitor()->dataPtr()->PDGName() << "\t"
  //   << aProgenitor->progenitor()->dataPtr()->PDGName() << endl;

  // masses of the particles
  mb_ = bProgenitor  ->progenitor()->momentum().mass();
  a_  = aProgenitor  ->progenitor()->momentum().mass() / mb_;
  c_  = cProgenitor  ->progenitor()->momentum().mass() / mb_; 
  a2_ = sqr(a_);
  c2_ = sqr(c_);

  // find rotation fgrom lab to frame with a along -z
  LorentzRotation eventFrame( bProgenitor->progenitor()->momentum().findBoostToCM() );
  Lorentz5Momentum pspectator = eventFrame*aProgenitor->progenitor()->momentum();
  eventFrame.rotateZ( -pspectator.phi() );
  eventFrame.rotateY( -pspectator.theta() - Constants::pi );

  //invert it
  eventFrame.invert();
  
  //generate the hard emission
  vector<Lorentz5Momentum> momenta = hardMomenta();

  // if no emission return
  if(momenta.empty()) {
    bProgenitor  ->maximumpT(pTmin_);
    cProgenitor  ->maximumpT(pTmin_);
    return HardTreePtr();
  }
  
  // rotate momenta back to the lab
  for(unsigned int ix=0;ix<momenta.size();++ix) {
    momenta[ix] *= eventFrame;
  
  }
  // get ParticleData objects
  tcPDPtr b = bProgenitor  ->progenitor()->dataPtr();
  tcPDPtr c = cProgenitor  ->progenitor()->dataPtr();
  tcPDPtr a = aProgenitor  ->progenitor()->dataPtr();
  tcPDPtr gluon  = getParticleData(ParticleID::g);
  // create new ShowerParticles
  ShowerParticlePtr emitter  (new_ptr(ShowerParticle(c,true )));
  ShowerParticlePtr spectator(new_ptr(ShowerParticle(a,true )));
  ShowerParticlePtr gauge    (new_ptr(ShowerParticle(gluon ,true )));
  ShowerParticlePtr incoming (new_ptr(ShowerParticle(b   ,false)));
  ShowerParticlePtr parent   (new_ptr(ShowerParticle(c,true )));
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
  // Set the maximum pt for all other emissions
  bProgenitor->maximumpT(pT_);
  cProgenitor  ->maximumpT(pT_);
  // Connect the particles with the branchings in the HardTree
  hardtree->connect( bProgenitor->progenitor(), spaceBranchings[0] );
  hardtree->connect( cProgenitor->progenitor(),   allBranchings[1] );
  hardtree->connect( aProgenitor->progenitor(),   allBranchings[2] );
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

double GeneralTwoBodyDecayer::threeBodyME(const ParticleVector & ) {
  throw Exception() << "Base class  GeneralTwoBodyDecayer::threeBodyME() "
		    << "called, should have an implementation in the inheriting class"
		    << Exception::runerror;
  return 0.;
}

vector<Lorentz5Momentum>  GeneralTwoBodyDecayer::hardMomenta() {

  double C    = 6.3;
  double ymax = 10.;
  double ymin = -ymax;
  
  vector<Lorentz5Momentum> particleMomenta (4);
  Energy2 lambda = sqr(mb_)* sqrt( 1. + sqr(a2_) + sqr(c2_) -
			           2.*a2_ - 2.*c2_ - 2.*a2_*c2_);    

  //Calculate A
  double A = (ymax - ymin) * C * (coupling_->overestimateValue() / 
				  (2.*Constants::pi));
 
  Energy pTmax = mb_* (sqr(1.-a_) - c2_) / (2.*(1.-a_));
  if (pTmax < pTmin_) particleMomenta.clear();

  while (pTmax >= pTmin_) {  
    //Generate pT, y and phi values
    Energy pT = pTmax * pow(UseRandom::rnd() , (1./A));  
    
    if (pT < pTmin_) {particleMomenta.clear(); break;}

    double phi = UseRandom::rnd() * Constants::twopi;
    double y   = ymin + UseRandom::rnd() * (ymax-ymin);

    double weight[2] = {0.,0.};
    double xa[2], xc[2], xc_z[2], xg;
    
    for (unsigned int j=0; j<2; j++) {
      //Check if the momenta are physical
      bool physical = calcMomenta(j, pT, y, phi, xg, xa[j], xc[j], xc_z[j], 
				  particleMomenta);
      if (not physical) continue;
      
      //Check if point lies within phase space
      bool inPS = psCheck(xg, xa[j]);
      if (not inPS) continue;
      
      //Calculate the ratio R/B
      double meRatio = matrixElementRatio();
      
      //Calculate jacobian
      Energy2 denom = (mb_ - particleMomenta[3].e()) * 
	               particleMomenta[2].vect().mag() -
		       particleMomenta[2].e() * particleMomenta[3].z(); 

      InvEnergy2 J  = (particleMomenta[2].vect().mag2()) / (2.* lambda * denom);
     
      //Calculate weight
      weight[j] = meRatio * fabs(sqr(pT)*J) * coupling_->ratio(pT*pT) / C;          
     }

    //    ofstream weights;
    //if (weight[0] + weight[1] > 1.){
    //weights.open("weights.top", ios::app);
    //weights << weight[0]+weight[1] << endl;
    //}

    //Accept point if weight > R
    if (weight[0] + weight[1] > UseRandom::rnd()) {
      if (weight[0] > (weight[0] + weight[1])*UseRandom::rnd()) {
	particleMomenta[1].setE( (mb_/2.)*xc  [0]);
	particleMomenta[1].setZ( (mb_/2.)*xc_z[0]);
	particleMomenta[2].setE( (mb_/2.)*xa  [0]);
	particleMomenta[2].setZ(-(mb_/2.)*sqrt(sqr(xa[0])-4.*a2_));
      }
      else {
	particleMomenta[1].setE( (mb_/2.)*xc  [1]);
	particleMomenta[1].setZ( (mb_/2.)*xc_z[1]);
	particleMomenta[2].setE( (mb_/2.)*xa  [1]);
	particleMomenta[2].setZ(-(mb_/2.)*sqrt(sqr(xa[1])-4.*a2_));
      }
      pT_ = pT;
      break;   
    }
    //If there's no splitting lower the pT
    pTmax = pT; 
  }
  return particleMomenta;
}

double GeneralTwoBodyDecayer::matrixElementRatio() {
  
  return 1.;
}

bool GeneralTwoBodyDecayer::calcMomenta(int j, Energy pT, double y, double phi,
   		        double& xg, double& xa, double& xc, double& xc_z,
			vector<Lorentz5Momentum>& particleMomenta){
  
  //Calculate xg
  xg = 2.*pT*cosh(y) / mb_;
  if (xg>(1. - sqr(c_ + a_)) || xg<0.) return false;

  //Calculate xa
  double xT  = 2.*pT / mb_;
  double A   = 4. - 4.*xg + sqr(xT);
  double B   = 4.*(3.*xg - 2. + 2.*c2_ - 2.*a2_ - sqr(xg) - xg*c2_ + xg*a2_);
  double L   = 1. + sqr(a2_) + sqr(c2_) - 2.*a2_ - 2.*c2_ - 2.*a2_*c2_;
  double det = 16.*( -L*sqr(xT) + pow(xT,4)*a2_ + 2.*xg*sqr(xT)*(1.-a2_-c2_) + 
		      L*sqr(xg) - sqr(xg*xT)*(1. + a2_) + pow(xg,4) + 
		      2.*pow(xg,3)*(- 1. + a2_ + c2_) );

  if (det<0.) return false;
  if (j==0) xa = (-B + sqrt(det))/(2.*A);
  if (j==1) xa = (-B - sqrt(det))/(2.*A);  
  if (xa>(1. + a2_ - c2_) || xa<2.*a_) return false;

  //Calculate xc
  xc = 2. - xa - xg;     
  if (xc>(1. + c2_ - a2_) || xc<2.*c_) return false;       

  //Calculate xc_z  
  double epsilon_p =  -sqrt(sqr(xa) - 4.*a2_) + xT*sinh(y) +
                       sqrt(sqr(xc) - 4.*c2_  - sqr(xT));
  double epsilon_m =  -sqrt(sqr(xa) - 4.*a2_) + xT*sinh(y) - 
                       sqrt(sqr(xc) - 4.*c2_  - sqr(xT));

  if (fabs(epsilon_p) < 1.e-10){
    xc_z =  sqrt(sqr(xc) - 4.*c2_ - sqr(xT));
  }
  else if (fabs(epsilon_m) < 1.e-10){
    xc_z = -sqrt(sqr(xc) - 4.*c2_ - sqr(xT));
  }
  else return false;

  //Check b is on shell
  if (fabs((sqr(xc) - sqr(xT) - sqr(xc_z) - 4.*c2_))>1.e-10) return false;

  //Calculate 4 momenta
  particleMomenta[0].setE   ( mb_);
  particleMomenta[0].setX   ( ZERO);
  particleMomenta[0].setY   ( ZERO);
  particleMomenta[0].setZ   ( ZERO);
  particleMomenta[0].setMass( mb_);

  particleMomenta[1].setE   ( mb_*xc/2.);
  particleMomenta[1].setX   (-pT*cos(phi));
  particleMomenta[1].setY   (-pT*sin(phi));
  particleMomenta[1].setZ   ( mb_*xc_z/2.);
  particleMomenta[1].setMass( mb_*c_);

  particleMomenta[2].setE   ( mb_*xa/2.);
  particleMomenta[2].setX   ( ZERO);
  particleMomenta[2].setY   ( ZERO);
  particleMomenta[2].setZ   (-mb_*sqrt(sqr(xa) - 4.*a2_)/2.);
  particleMomenta[2].setMass( mb_*a_);

  particleMomenta[3].setE   ( pT*cosh(y));
  particleMomenta[3].setX   ( pT*cos(phi));
  particleMomenta[3].setY   ( pT*sin(phi));
  particleMomenta[3].setZ   ( pT*sinh(y));
  particleMomenta[3].setMass( ZERO);
 
  return true;
}


bool GeneralTwoBodyDecayer::psCheck(double xg, double xa) {
  
  //Check is point is in allowed region of phase space
  double xc_star = (1. - a2_ + c2_ - xg) / sqrt(1. - xg);
  double xg_star = xg / sqrt(1. - xg);

  if ((sqr(xc_star) - 4.*c2_) < 1e-10) return false;
  double xa_max = (4. + 4.*a2_ - sqr(xc_star + xg_star) + 
		   sqr(sqrt(sqr(xc_star) - 4.*c2_) + xg_star)) / 4.;
  double xa_min = (4. + 4.*a2_ - sqr(xc_star + xg_star) + 
		   sqr(sqrt(sqr(xc_star) - 4.*c2_) - xg_star)) / 4.;

  if (xa < xa_min || xa > xa_max) return false;

  return true;
}
