// -*- C++ -*-
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

void GeneralTwoBodyDecayer::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  if(_theVertex) {
    _theVertex->init();
  }
  else 
    throw InitException() << "GeneralTwoBodyDecayer::doinit() - "
			  << "Null vertex pointer!\n";  

  vector<double> wgt(0);  
  PDVector parents = _theVertex->getIncoming();
  PDVector::size_type np = parents.size();
  PDVector extpart(3);
  for( tPDVector::size_type i = 0; i < np; ++i ) {
    tPDPtr inpart = parents[i];
    long pid = inpart->id();
    if( pid < 0 ) continue;
    Energy m1 = inpart->mass();
    PDVector decaylist;
    for(unsigned int il = 0; il< _thelist.size(); ++il) {
       PDVector temp = _theVertex->search(_thelist[il], pid);
       decaylist.insert(decaylist.end(),temp.begin(),temp.end());
    }
    PDVector::size_type ndec = decaylist.size();
    for( PDVector::size_type j = 0; j < ndec; j +=3 ) {
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
				      const PDVector & children) const {
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
			<< "in GeneralTwoBodyDecayer::decay() "
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
    else
      throw Exception() << "Unknown outgoing colours for decaying "
			<< "colour octet "
			<< "in GeneralTwoBodyDecayer::decay() "
			<< outaColour << " " << outbColour
			<< Exception::runerror;
  }
  else
    throw Exception() << "Unknown incoming colour in "
		      << "GeneralTwoBodyDecayer::decay() "
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
    ("IncomingPart",
     "PDG Codes for incoming particles",
     &GeneralTwoBodyDecayer::_inpart, 0, -1, 0, 0,
     false, false, false);

  static ParVector<GeneralTwoBodyDecayer,int> interfaceOutgoingPartA
    ("OutgoingPartA",
     "PDG Codes for first set of outgoing particles",
     &GeneralTwoBodyDecayer::_outparta, 0, -1, 0, 0,
     false, false, false);

  static ParVector<GeneralTwoBodyDecayer,int> interfaceOutgoingPartB
    ("OutgoingPartB",
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
