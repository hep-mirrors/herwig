// -*- C++ -*-
//
// GeneralTwoBodyDecayer.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
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
  _theVertex->init();
  vector<double> wgt;  
  set<tPDPtr> parents = _theVertex->incoming();
  for( set<tPDPtr>::const_iterator it = parents.begin(); 
       it != parents.end(); ++it ) {
    long pid = (*it)->id();
    if( pid < 0 ) continue;
    tcPDPtr inpart = *it;
    Energy m1 = inpart->mass();
    tPDVector decaylist;
    for(unsigned int il = 0; il< _thelist.size(); ++il) {
       tPDVector temp = _theVertex->search(_thelist[il], inpart);
       decaylist.insert(decaylist.end(),temp.begin(),temp.end());
    }
    tPDVector::size_type ndec = decaylist.size();
    for( tPDVector::size_type j = 0; j < ndec; j +=3 ) {
      tPDPtr pa(decaylist[j]), pb(decaylist[j + 1]), pc(decaylist[j + 2]);
      if( pb->id() == pid ) swap(pa, pb);
      if( pc->id() == pid ) swap(pa, pc);
      //allowed on-shell decay?
      if( m1 <= pb->mass() + pc->mass() ) continue;
      //vertices are defined with all particles incoming
      if( pb->CC() ) pb = pb->CC();
      if( pc->CC() ) pc = pc->CC();
      //store ids so that the decayer knows what it is allowed to 
      //decay
      _inpart.push_back(pid); _outparta.push_back(pb->id());
      _outpartb.push_back(pc->id());
      //create phase space mode
      tPDVector extpart(3);
      extpart[0] = pa;
      extpart[1] = pb;
      extpart[2] = pc;
      addMode(new_ptr(DecayPhaseSpaceMode(extpart, this)), 
	      _maxweight[0], wgt);
      
    }
  }
  unsigned int isize(_inpart.size()), oasize(_outparta.size()),
    obsize(_outpartb.size());
  if(  isize == 0 ||  oasize == 0 || obsize == 0 )
    throw InitException()
      << "GeneralTwoBodyDecayer::doinit() - Atleast one of the particle "
      << "vectors has zero size, cannot continue." 
      << isize << " " << oasize << " " << obsize 
      << Exception::abortnow;
  
  if(  isize != oasize || isize != obsize )
    throw InitException()
      << "GeneralTwoBodyDecayer::doinit() - The particle vectors have "
      << "different sizes. " << isize << " " << oasize << " " << obsize
      << Exception::abortnow;

}

int GeneralTwoBodyDecayer::modeNumber(bool & cc, tcPDPtr parent, 
				      const tPDVector & children) const {
  long parentID = parent->id();
  long id1 = children[0]->id();
  long id2 = children[1]->id();
  int imode(-1);
  unsigned ii(0), nipart(_inpart.size());
  cc = false;
  do {
    long listpid(_inpart[ii]), listid1(_outparta[ii]),
      listid2(_outpartb[ii]);
    if( parentID == listpid && 
	((id1 == listid1 && id2 == listid2) || 
	 (id1 == listid2 && id2 == listid1)) )
      imode = ii;
    //cc-mode
    else if(parentID == -listpid) {
      cc = true;
      if((id1 == -listid1 && id2 == -listid2) || 
	 (id1 == -listid2 && id2 == -listid1) ||
	 (id1 == listid1 && id2 == -listid2)  || 
	 (id1 == -listid1 && id2 == listid2)  ||
	 (id1 == listid2 && id2 == -listid1)  || 
	 (id1 == -listid2 && id2 == listid1) )
	imode = ii;
      else ++ii;
    }
    else ++ii;	
  }
  while( imode < 0 && ii < nipart );
  return imode;
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
  bool order(false);
  vector<int>::size_type ix(0);
  do {
    if( parent == _inpart[ix] ) {
      long id1t(_outparta[ix]), id2t(_outpartb[ix]);
      if( id1 == id1t && id2 == id2t ) {
	order = true;
	break;
      }
      if( id1 == id2t && id2 == id1t ) {
	order = false;
	break;
      }
    }
    ++ix;
  }
  while( ix < _inpart.size() );
  mecode = -1;
  coupling = 1.;
  return order;
}


void GeneralTwoBodyDecayer::persistentOutput(PersistentOStream & os) const {
  os << _thelist << _theVertex << _inpart << _outparta << _outpartb
     << _maxweight;
}

void GeneralTwoBodyDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _thelist >> _theVertex >>_inpart >>_outparta >>_outpartb
     >> _maxweight;
}

AbstractClassDescription<GeneralTwoBodyDecayer> 
GeneralTwoBodyDecayer::initGeneralTwoBodyDecayer;
// Definition of the static class description member.

void GeneralTwoBodyDecayer::Init() {

  static ClassDocumentation<GeneralTwoBodyDecayer> documentation
    ("This class is designed to be a base class for all 2 body decays"
     "in a general model");

  static Reference<GeneralTwoBodyDecayer,Helicity::VertexBase> interfaceDecayVertex
    ("DecayVertex",
     "Pointer to decayer vertex",
     &GeneralTwoBodyDecayer::_theVertex, false, false, true, false);
  
  static ParVector<GeneralTwoBodyDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "Maximum weight for integration",
     &GeneralTwoBodyDecayer::_maxweight, 1.0, -1, 0, 0,
     false, false, false,&GeneralTwoBodyDecayer::setWeight, 0 ,0, 0, 0);

  static ParVector<GeneralTwoBodyDecayer,int> interfaceIncomingPart
    ("Incoming",
     "PDG Codes for incoming particles",
     &GeneralTwoBodyDecayer::_inpart, 0, -1, 0, 0,
     false, false, false);

  static ParVector<GeneralTwoBodyDecayer,int> interfaceOutgoingPartA
    ("OutgoingA",
     "PDG Codes for first set of outgoing particles",
     &GeneralTwoBodyDecayer::_outparta, 0, -1, 0, 0,
     false, false, false);

  static ParVector<GeneralTwoBodyDecayer,int> interfaceOutgoingPartB
    ("OutgoingB",
     "PDG Codes for second set of outgoing particles",
     &GeneralTwoBodyDecayer::_outpartb, 0, -1, 0, 0,
     false, false, false);
 
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
