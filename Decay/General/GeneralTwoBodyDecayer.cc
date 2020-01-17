
// -*- C++ -*-
//
// GeneralTwoBodyDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GeneralTwoBodyDecayer class.
//

#include "GeneralTwoBodyDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
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
  PerturbativeDecayer::doinit();
  assert( incoming_ && outgoing_.size()==2);
  //create phase space mode
  tPDVector extpart(3);
  extpart[0] = incoming_;
  extpart[1] = outgoing_[0];
  extpart[2] = outgoing_[1];
  addMode(new_ptr(DecayPhaseSpaceMode(extpart, this)), maxWeight_, vector<double>());
}

int GeneralTwoBodyDecayer::modeNumber(bool & cc, tcPDPtr parent, 
				      const tPDVector & children) const {
  long parentID = parent->id();
  long id1 = children[0]->id();
  long id2 = children[1]->id();
  cc = false;
  long out1 = outgoing_[0]->id();
  long out2 = outgoing_[1]->id();
  if( parentID == incoming_->id() && 
      ((id1 == out1 && id2 == out2) || 
       (id1 == out2 && id2 == out1)) ) {
    return 0;
  }
  else if(incoming_->CC() && parentID == incoming_->CC()->id()) {
    cc = true;
    if( outgoing_[0]->CC()) out1 = outgoing_[0]->CC()->id();
    if( outgoing_[1]->CC()) out2 = outgoing_[1]->CC()->id();
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
  assert(dm.parent()->id() == incoming_->id());
  ParticleMSet::const_iterator pit = dm.products().begin();
  long id1 = (*pit)->id();
  ++pit;
  long id2 = (*pit)->id();
  long id1t(outgoing_[0]->id()), id2t(outgoing_[1]->id());
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
  os << incoming_ << outgoing_ << maxWeight_;
}

void GeneralTwoBodyDecayer::persistentInput(PersistentIStream & is, int) {
  is >> incoming_ >> outgoing_ >> maxWeight_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<GeneralTwoBodyDecayer,PerturbativeDecayer>
describeHerwigGeneralTwoBodyDecayer("Herwig::GeneralTwoBodyDecayer", "Herwig.so");

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
  PerturbativeDecayer::doinitrun();
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

void GeneralTwoBodyDecayer::decayInfo(PDPtr incoming, PDPair outgoing) {
  incoming_=incoming;
  outgoing_.clear();
  outgoing_.push_back(outgoing.first );
  outgoing_.push_back(outgoing.second);
}

double GeneralTwoBodyDecayer::matrixElementRatio(const Particle & inpart, 
						 const ParticleVector & decay2,
						 const ParticleVector & decay3, 
						 MEOption meopt,
						 ShowerInteraction inter) {
  // calculate R/B
  double B = me2        (0, inpart, decay2, meopt);    
  double R = threeBodyME(0, inpart, decay3, inter, meopt);
  return R/B;
  
}

const vector<DVector> & GeneralTwoBodyDecayer::getColourFactors(const Particle & inpart, 
								const ParticleVector & decay,
								unsigned int & nflow) {
  // calculate the colour factors for the three-body decay
  vector<int> sing,trip,atrip,oct,sex,asex;
  for(unsigned int it=0;it<decay.size();++it) {
    if     (decay[it]->dataPtr()->iColour() == PDT::Colour0    ) sing. push_back(it);
    else if(decay[it]->dataPtr()->iColour() == PDT::Colour3    ) trip. push_back(it);
    else if(decay[it]->dataPtr()->iColour() == PDT::Colour3bar ) atrip.push_back(it);
    else if(decay[it]->dataPtr()->iColour() == PDT::Colour8    ) oct.  push_back(it);
    else if(decay[it]->dataPtr()->iColour() == PDT::Colour6    ) sex.  push_back(it);
    else if(decay[it]->dataPtr()->iColour() == PDT::Colour6bar ) asex. push_back(it);
  }

  // identical particle symmetry factor
  double symFactor=1.;
  if (( sing.size()==2 && decay[ sing[0]]->id()==decay[ sing[1]]->id()) ||
      ( trip.size()==2 && decay[ trip[0]]->id()==decay[ trip[1]]->id()) ||
      (atrip.size()==2 && decay[atrip[0]]->id()==decay[atrip[1]]->id()) ||
      (  oct.size()==2 && decay[  oct[0]]->id()==decay[  oct[1]]->id()) ||
      (  sex.size()==2 && decay[  sex[0]]->id()==decay[  sex[1]]->id()) ||
      ( asex.size()==2 && decay[ asex[0]]->id()==decay[ asex[1]]->id()))
    symFactor /= 2.;
  else if (oct.size()==3 && 
	   decay[oct[0]]->id()==decay[oct[1]]->id() &&
	   decay[oct[0]]->id()==decay[oct[2]]->id())
    symFactor /= 6.;
  
  colour_ = vector<DVector>(1,DVector(1,symFactor*1.));
  
  // decaying colour singlet   
  if(inpart.dataPtr()->iColour() == PDT::Colour0) {
    if(trip.size()==1 && atrip.size()==1 && oct.size()==1) {
      nflow = 1;
      colour_ = vector<DVector>(1,DVector(1,symFactor*4.));
    }
    else if(trip.size()==1 && atrip.size()==1 && sing.size()==1) {
      nflow = 1;
      colour_ = vector<DVector>(1,DVector(1,symFactor*3.));
    }
    else if (oct.size()==3){
      nflow = 1;
      colour_ = vector<DVector>(1,DVector(1,symFactor*24.));
    }
    else if(sing.size()==3) {
      nflow = 1;
      colour_ = vector<DVector>(1,DVector(1,symFactor));
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
    else if(trip.size()==1 && sing.size()==2) {
      nflow = 1;
      colour_ = vector<DVector>(1,DVector(1,symFactor));
    }
    else if(trip.size()==1 && oct.size()==2) {
      nflow = 2;
      colour_.clear();
      colour_.resize(2,DVector(2,0.));
      colour_[0][0] =  symFactor*16./9.; colour_[0][1] = -symFactor*2./9.;
      colour_[1][0] = -symFactor*2./9.;  colour_[1][1] =  symFactor*16./9.;
    }
    else if(atrip.size()==2 && sing.size()==1) {
      nflow = 1;
      colour_ = vector<DVector>(1,DVector(1,symFactor*2.));
    }
    else if(atrip.size()==2 && oct.size()==1) {
      nflow = 2;
      colour_.clear();
      colour_.resize(2,DVector(2,0.));
      colour_[0][0] =  8./3.*symFactor; colour_[0][1] =-4./3.*symFactor;
      colour_[1][0] = -4./3.*symFactor; colour_[1][1] = 8./3.*symFactor;
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
    else if(atrip.size()==1 && sing.size()==2) {
      nflow = 1;
      colour_ = vector<DVector>(1,DVector(1,symFactor));
    }
    else if(atrip.size()==1 && oct.size()==2){
      nflow = 2;
      colour_.clear();
      colour_ .resize(2,DVector(2,0.));
      colour_[0][0] =  symFactor*16./9.; colour_[0][1] = -symFactor*2./9.;
      colour_[1][0] = -symFactor*2./9.;  colour_[1][1] =  symFactor*16./9.;
    }
    else if(trip.size()==2 && sing.size()==1) {
      nflow = 1;
      colour_ = vector<DVector>(1,DVector(1,symFactor*2.));
    }
    else if(trip.size()==2 && oct.size()==1) {
      nflow = 2;
      colour_.clear();
      colour_.resize(2,DVector(2,0.));
      colour_[0][0] =  8./3.*symFactor; colour_[0][1] =-4./3.*symFactor;
      colour_[1][0] = -4./3.*symFactor; colour_[1][1] = 8./3.*symFactor;
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
    else if (sing.size()==1 && trip.size()==1 && atrip.size()==1) {
      nflow = 1;
      colour_ = vector<DVector>(1,DVector(1,symFactor*0.5));
    }
    else if (oct.size()==2 && sing.size()==1) {
      nflow = 1;
      colour_ = vector<DVector>(1,DVector(1,symFactor*3.));
    }
    else
      throw Exception() << "Unknown colour for the outgoing particles"
			<< " for a decaying colour octet particle in "
			<< "GeneralTwoBodyDecayer::getColourFactors() for "
			<< inpart.   dataPtr()->PDGName() << " -> "
			<< decay[0]->dataPtr()->PDGName() << " " 
			<< decay[1]->dataPtr()->PDGName() << " "  
			<< decay[2]->dataPtr()->PDGName()
			<< Exception::runerror;
  }
  // Sextet
  else if(inpart.dataPtr()->iColour() == PDT::Colour6) {
    if(trip.size()==2 && sing.size()==1) {
      nflow = 1;
      colour_ = vector<DVector>(1,DVector(1,symFactor));
    }
    else if(trip.size()==2 && oct.size()==1) {
      nflow = 2;
      colour_.clear();
      colour_.resize(2,DVector(2,0.));
      colour_[0][0] =  4./3.*symFactor; colour_[0][1] = 1./3.*symFactor;
      colour_[1][0] =  1./3.*symFactor; colour_[1][1] = 4./3.*symFactor;
    }
    else
      throw Exception() << "Unknown colour for the outgoing particles"
			<< " for a decaying colour sextet particle in "
			<< "GeneralTwoBodyDecayer::getColourFactors() for "
			<< inpart.   dataPtr()->PDGName() << " -> "
			<< decay[0]->dataPtr()->PDGName() << " " 
			<< decay[1]->dataPtr()->PDGName() << " "  
			<< decay[2]->dataPtr()->PDGName()
			<< Exception::runerror;
  }
  // anti Sextet
  else if(inpart.dataPtr()->iColour() == PDT::Colour6bar) {
    if(atrip.size()==2 && sing.size()==1) {
      nflow = 1;
      colour_ = vector<DVector>(1,DVector(1,symFactor));
    }
    else if(atrip.size()==2 && oct.size()==1) {
      nflow = 2;
      colour_.clear();
      colour_.resize(2,DVector(2,0.));
      colour_[0][0] =  4./3.*symFactor; colour_[0][1] = 1./3.*symFactor;
      colour_[1][0] =  1./3.*symFactor; colour_[1][1] = 4./3.*symFactor;
    }
    else
      throw Exception() << "Unknown colour for the outgoing particles"
			<< " for a decaying colour anti-sextet particle in "
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

const GeneralTwoBodyDecayer::CFlow & 
GeneralTwoBodyDecayer::colourFlows(const Particle & inpart,
				   const ParticleVector & decay) {
  // static initialization of commonly used colour structures
  static const CFlow init = CFlow(3, CFlowPairVec(1, make_pair(0, 1.)));
  static CFlow tripflow = init;
  static CFlow atripflow = init;
  static CFlow octflow = init;
  static CFlow sexflow = init;
  static CFlow epsflow = init;
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

    sexflow[0].resize(2, make_pair(0,1.));
    sexflow[0][1] = make_pair(1,1.);
    sexflow[2][0] = make_pair(1,1.);

    epsflow[0].resize(2, make_pair(0,-1.));
    epsflow[0][0] = make_pair(0,-1.);
    epsflow[0][1] = make_pair(1,-1.);
    epsflow[2][0] = make_pair(1,1.);
    initialized = true;
  }
  

  // main function body
  int sing=0,trip=0,atrip=0,oct=0,sex=0,asex=0;
  for (size_t it=0; it<decay.size(); ++it) {
    switch ( decay[it]->dataPtr()->iColour() ) {
    case PDT::Colour0:     ++sing; break;
    case PDT::Colour3:     ++trip; break;
    case PDT::Colour3bar: ++atrip; break;
    case PDT::Colour8:      ++oct; break;
    case PDT::Colour6:      ++sex; break;
    case PDT::Colour6bar:  ++asex; break;
      /// @todo: handle these better
    case PDT::ColourUndefined:     break;
    case PDT::Coloured:            break;
    }
  }

  const CFlow * retval = 0;
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
  // decaying colour sextet
  else if(inpart.dataPtr()->iColour() ==PDT::Colour6 &&
	  oct==1 && trip==2) {
    retval = &sexflow;
  }
  // decaying anti colour sextet
  else if(inpart.dataPtr()->iColour() ==PDT::Colour6bar &&
	  oct==1 && atrip==2) {
    retval = &sexflow;
  }
  // decaying colour triplet (eps)
  else if(inpart.dataPtr()->iColour() == PDT::Colour3 &&
     atrip==2 && oct==1) {
    retval = &epsflow;
  }
  // decaying colour anti-triplet (eps)
  else if(inpart.dataPtr()->iColour() == PDT::Colour3bar &&
	  trip==2 && oct==1){
    retval = &epsflow;
  }
  else {
    retval = &fpflow;  
  }

  return *retval;
}

double GeneralTwoBodyDecayer::threeBodyME(const int , const Particle &,
					  const ParticleVector &,
					  ShowerInteraction, MEOption) {
  throw Exception() << "Base class PerturbativeDecayer::threeBodyME() "
		    << "called, should have an implementation in the inheriting class"
		    << Exception::runerror;
  return 0.;
}
