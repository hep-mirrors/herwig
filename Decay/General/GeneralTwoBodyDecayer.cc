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
  PDVector inc(_theVertex->getIncoming());
  for(unsigned int i = 0; i < inc.size(); ++i) {
    int id = inc[i]->id();
    if(id < 0)  continue;
    Energy m1 = inc[i]->mass();
    PDVector decaylist(0);
    for(unsigned int il=0;il<_thelist.size();++il) {
      PDVector dtemp = _theVertex->search(_thelist[il],id);
      decaylist.insert(decaylist.end(),dtemp.begin(),dtemp.end());
    }
    PDVector extpart(3);
    for(PDVector::iterator iter=decaylist.begin();iter!=decaylist.end();) {
      Energy m2(0.*MeV),m3(0.*MeV);
      bool cc1(false),cc2(false),cc3(false);
      if((*iter)->CC()) {cc1=true;}
      if((*(iter+1))->CC()) {cc2=true;}
      if((*(iter+2))->CC()) {cc3=true;}

      if((*iter)->id()==id) {
	m2 = (*(iter+1))->mass();
	m3 = (*(iter+2))->mass();
	extpart[0] = *iter;
	if(cc1) {
	  if(cc2) {extpart[1] = (*(iter+1))->CC();}
	  else {extpart[1] = *(iter+1);}
	  
	  if(cc3) {extpart[2] = (*(iter+2))->CC();}
	  else {extpart[2] = *(iter+2);}
	}
	else {
	  extpart[1] = *(iter+1);
	  extpart[2] = *(iter+2);	  
	}
      }
      else if((*(iter+1))->id()==id) {
	m2 = (*iter)->mass();
	m3 = (*(iter+2))->mass();
 	extpart[0] = *(iter+1);
	if(cc2) {
	  if(cc1) {extpart[1] = (*iter)->CC();}
	  else {extpart[1] = *iter;}
	  
	  if(cc3) {extpart[2] = (*(iter+2))->CC();}
	  else {extpart[2] = *(iter+2);}
	}
	else {
	  extpart[1] = *iter;
	  extpart[2] = *(iter+2);	  
	}
      }
      else {
	m2 = (*iter)->mass();
	m3 = (*(iter+1))->mass();
	extpart[0]=*(iter+2);

	if(cc3) {
	  if(cc1) {extpart[1] = (*iter)->CC();}
	  else {extpart[1] = *iter;}
	  
	  if(cc2) {extpart[2] = (*(iter+1))->CC();}
	  else {extpart[2] = *(iter+1);}
	}
	else {
	  extpart[1] = *iter;
	  extpart[2] = *(iter+1);	  
	}
	
      }
      if(m1 <= (m2 + m3)) { 
	decaylist.erase(iter,iter+3);
      }
      else {
	_inpart.push_back(extpart[0]->id());
	_outparta.push_back(extpart[1]->id());
	_outpartb.push_back(extpart[2]->id());
	DecayPhaseSpaceModePtr mode;
	mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
	addMode(mode,_maxweight[0],wgt);
	iter+=3;
      }
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


