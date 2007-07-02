// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HardProcessConstructor class.
//

#include "HardProcessConstructor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "Herwig++/MatrixElement/GeneralHardME.h"

using namespace Herwig;

HardProcessConstructor::HardProcessConstructor() : 
  theNout(0), theNv(0),theProcesses(0), theAllDiagrams(false), 
  the33bto33b(4, DVector(4, 0.)),  the33bpto33bp(3, DVector(3, 0.)),
  the33bto88(2, DVector(4, 0.)), the88to88(2, DVector(4, 0.)) {
  //set-up colour factor matrices
  //33b->33b and similar
  the33bto33b[0][0] = the33bto33b[1][1] = 2.; 
  the33bto33b[2][2] = the33bto33b[3][3] = 9.;
  the33bto33b[0][1] = the33bto33b[1][0] = -2./3.;
  the33bto33b[0][2] = the33bto33b[2][0] = 0.;
  the33bto33b[0][3] = the33bto33b[3][0] = 4.;
  the33bto33b[1][2] = the33bto33b[2][1] = 4.;
  the33bto33b[1][3] = the33bto33b[3][1] = 0.;
  the33bto33b[2][3] = the33bto33b[3][2] = 3.;
  
  //33b'->33b' and similar
  the33bpto33bp[0][0] = 2;
  the33bpto33bp[1][1] = the33bpto33bp[2][2] = 9.;
  the33bpto33bp[0][1] = the33bpto33bp[1][0] = 0.;
  the33bpto33bp[0][2] = the33bpto33bp[2][0] = 4.;
  the33bpto33bp[1][2] = the33bpto33bp[1][2] = 3.;

  //33b->88 and similar
  the33bto88[0][0] = the33bto88[1][1] = 16./3.;
  the33bto88[0][1] = the33bto88[1][0] = -2./3.;

  //88->88 and similar
  the88to88[0][0] = the88to88[1][1] = 72.;
  the88to88[0][1] = the88to88[1][0] = 36.;
  
}

void HardProcessConstructor::doinit() throw(InitException) {
  Interfaced::doinit();
  theNout = theOutgoing.size();
  PDVector::size_type ninc = theIncoming.size();
  // exit if nothing to do
  if(theNout==0||ninc==0) return;
  //create vector of initial-state pairs
  for(PDVector::size_type i = 0; i < ninc; ++i) {
    for(PDVector::size_type j = 0; j < ninc; ++j) {
      tcPDPair inc = make_pair(theIncoming[i], theIncoming[j]);
      
      if( (inc.first->iSpin() > inc.second->iSpin()) ||
	  (inc.first->iSpin() == inc.second->iSpin() &&
	   inc.first->id() < inc.second->id()) )
	swap(inc.first, inc.second);

      if( !duplicate(inc) ) {
	theIncPairs.push_back(inc);
      }
    }
  }
  //set up pointers
  EGPtr eg = generator();
  theModel = dynamic_ptr_cast<tHwSMPtr>(eg->standardModel());
  if(!theModel)
    throw InitException() << "HardProcessConstructor::doinit - "
			  << "The model pointer is null!"
			  << Exception::abortnow;
  if(eg->eventHandler()) {
    string subProcessName = 
      eg->preinitInterface(eg->eventHandler(), "SubProcessHandlers", "get", "");
    theSubProcess = eg->getObject<SubProcessHandler>(subProcessName);
    if(!theSubProcess)
      throw InitException() << "HardProcessConstructor::doinit - "
			    << "The SubProcessHandler pointer is null!"
			    << Exception::abortnow;
  }
  else
    throw InitException() 
      << "HardProcessConstructor::doinit - generator()->eventHandler()"
      << " has returned a null pointer. Cannot access SubProcessHandler."
      << Exception::abortnow;
}

bool HardProcessConstructor::duplicate(tcPDPair ppair) const {
  vector<tcPDPair>::size_type psize = theIncPairs.size();
  if(psize == 0) return false;
  bool found(false);
  for(vector<tcPDPair>::size_type i = 0; i < psize; ++i) {
    long id1 = theIncPairs[i].first->id();
    long id2 = theIncPairs[i].second->id();
    if( (ppair.first->id() == id1 && ppair.second->id() == id2) ||
	(ppair.first->id() == id2 && ppair.second->id() == id1) ) {
      found = true;
      break;
    }
  }
  return found;
}

void HardProcessConstructor::persistentOutput(PersistentOStream & os) const {
  os << theIncoming << theOutgoing << theModel << theAllDiagrams << theSubProcess
     << the33bto33b << the33bpto33bp << the33bto88 << the88to88;
}

void HardProcessConstructor::persistentInput(PersistentIStream & is, int) {
  is >> theIncoming >> theOutgoing  >> theModel >> theAllDiagrams >> theSubProcess
     >> the33bto33b >> the33bpto33bp >> the33bto88 >> the88to88;
  theNout = 0;
  theNv = 0;
}

ClassDescription<HardProcessConstructor> 
HardProcessConstructor::initHardProcessConstructor;
// Definition of the static class description member.

void HardProcessConstructor::Init() {

  static ClassDocumentation<HardProcessConstructor> documentation
    ("HardProcessConstructor constructs the possible diagrams for "
     "a process given the external particles");
 
  static RefVector<HardProcessConstructor,ThePEG::ParticleData> interfaceIn
    ("Incoming",
     "Pointers to incoming particles",
     &HardProcessConstructor::theIncoming, -1, false, false, true, false);

 static RefVector<HardProcessConstructor,ThePEG::ParticleData> interfaceOut
    ("Outgoing",
     "Pointers to incoming particles",
     &HardProcessConstructor::theOutgoing, -1, false, false, true, false);
 
 static Switch<HardProcessConstructor,bool> interfaceIncludeAllDiagrams
   ("QCDandEW",
    "Switch to decide which diagrams to include in ME calc.",
     &HardProcessConstructor::theAllDiagrams, 0, false, false);
 static SwitchOption interfaceIncludeAllDiagramsOff
   (interfaceIncludeAllDiagrams,
    "Off",
    "Do not include all diagrams in ME calc, only those with strong coupling in them",
    0);
 static SwitchOption interfaceIncludeAllDiagramsOn
   (interfaceIncludeAllDiagrams,
     "On",
    "Include all diagrams in the ME calculation",
    1);
}

void HardProcessConstructor::constructDiagrams() {
  if(theIncPairs.empty() || theOutgoing.empty()) return;
  // delete the matrix elements we already have
  int nme = theSubProcess->MEs().size();
  for(int ix = nme - 1; ix >= 0; --ix)
    generator()->preinitInterface(theSubProcess, "MatrixElements", 
				  ix, "erase", "");
  theModel->init();
  theNv = theModel->numberOfVertices();
  //make sure  vertices are initialised
  for(unsigned int ix = 0; ix < theNv; ++ix )
    theModel->vertex(ix)->init();
  //Create necessary diagrams
  vector<tcPDPair>::size_type is;
  PDVector::size_type os;
  for(is = 0; is < theIncPairs.size(); ++is) {
    tcPDPair ppi = theIncPairs[is]; 
    for(os = 0; os < theNout; ++os) { 
      long fs = theOutgoing[os]->id();
      for(size_t iv = 0; iv < theNv; ++iv) {
	tVertexBasePtr vertexA = theModel->vertex(iv);
	if( vertexA->getName() == GeneralSVV ||
	    (!theAllDiagrams && vertexA->orderInGs() == 0) ) 
	  continue;

	if(vertexA->getNpoint() == 3) {
	  //scattering diagrams
	  createTChannels(ppi, fs, vertexA);
	  
	  //resonance diagrams
	  if( vertexA->incoming(ppi.first->id()) &&  
	      vertexA->incoming(ppi.second->id()) )
	    createSChannels(ppi, fs, vertexA);
	}
	else 
	  makeFourPointDiagrams(ppi.first->id(), ppi.second->id(),
				fs, vertexA);
      }
    }
  }
  // We now have a vector of diagrams that need to be grouped together
  HPDVector::iterator ita, itb;
  for(ita = theProcesses.begin(); ita != theProcesses.end(); ) {
    HPDVector group;
    HPDiagram current = *ita;
    group.push_back(current);
    for(itb = ita + 1; itb != theProcesses.end();) {
      if( (*itb).sameProcess(current) ) {
	group.push_back(*itb);
	theProcesses.erase(itb);
      }
      else
	++itb;
    }
    theProcesses.erase(ita);
    if( !group.empty() ) {
      createMatrixElement(group);
      group.clear();
    }
  }
  
}

void HardProcessConstructor::
createSChannels(tcPDPair inpp, long fs, tVertexBasePtr vertex) {
  //Have 2 incoming particle and a vertex, find the possible offshell
  //particles
  pair<long,long> inc = make_pair(inpp.first->id(), inpp.second->id());
  PDSet offshells = search(vertex, inpp.first->id(), incoming,
			   inpp.second->id(), incoming, outgoing);
  PDSet::const_iterator it;
  for(it = offshells.begin(); it != offshells.end(); ++it) {
    for(size_t iv = 0; iv < theNv; ++iv) {
      tVertexBasePtr vertexB = theModel->vertex(iv);
      if( vertexB->getNpoint() != 3) continue;
      if( vertexB->getName() == GeneralSVV ||
	  (!theAllDiagrams && vertexB->orderInGs() == 0) ) 
	continue;
      
      PDSet final;
      if( vertexB->outgoing(fs) &&
	  vertexB->incoming((*it)->id()) )
	final = search(vertexB, (*it)->id(), incoming, fs,
		       outgoing, outgoing);
      //Now make diagrams
	if(!final.empty()) 
	  makeDiagrams(inc, fs, final, *it, HPDiagram::sChannel,
		       make_pair(vertex, vertexB), make_pair(true,true));
    }
  }

}

void HardProcessConstructor::
createTChannels(tcPDPair inpp, long fs, tVertexBasePtr vertex) {
  pair<long,long> inc = make_pair(inpp.first->id(), inpp.second->id());
  //first try a with c
  PDSet offshells = search(vertex, inpp.first->id(), incoming, fs,
			   outgoing, outgoing);
  PDSet::const_iterator it;
  for(it = offshells.begin(); it != offshells.end(); ++it) {
     for(size_t iv = 0; iv < theNv; ++iv) {
       tVertexBasePtr vertexB = theModel->vertex(iv);
       if( vertexB->getNpoint() != 3 ) continue;
       if( vertexB->getName() == GeneralSVV ||
	  (!theAllDiagrams && vertexB->orderInGs() == 0) ) 
	continue;
       PDSet final;
       if( vertexB->incoming(inc.second) )
	 final = search(vertexB, inc.second, incoming, (*it)->id(),
			incoming, outgoing);
       if( !final.empty() )
	 makeDiagrams(inc, fs, final, *it, HPDiagram::tChannel,
		      make_pair(vertex, vertexB), make_pair(true,true));
     }
  }
  //now try b with c
  offshells = search(vertex, inpp.second->id(), incoming, fs,
			   outgoing, incoming);
  for(it = offshells.begin(); it != offshells.end(); ++it) {
    for(size_t iv = 0; iv < theNv; ++iv) {
       tVertexBasePtr vertexB = theModel->vertex(iv);
       if( vertexB->getNpoint() != 3 ) continue;
       if( vertexB->getName() == GeneralSVV ||
	  (!theAllDiagrams && vertexB->orderInGs() == 0) ) 
	continue;

       PDSet final;
       if( vertexB->incoming(inc.first) )
	 final = search(vertexB, inc.first, incoming, (*it)->id(),
			outgoing, outgoing);
       if( !final.empty() )
	 makeDiagrams(inc, fs, final, *it, HPDiagram::tChannel,
		      make_pair(vertexB, vertex), make_pair(true, false));
    }
  }

}

void HardProcessConstructor::makeFourPointDiagrams(long parta, long partb,
						   long partc, VBPtr vert) {
  PDSet ext = search(vert, parta, incoming, partb,incoming, partc, outgoing);
  if(ext.size() > 0) {
    IDPair in(parta, partb);
    for(PDSet::const_iterator iter=ext.begin(); iter!=ext.end();
	++iter) {
      HPDiagram nhp;
      nhp.incoming = in; 
      nhp.outgoing = make_pair(partc, (*iter)->id());
      nhp.vertices = make_pair(vert, vert);
      nhp.channelType = HPDiagram::fourPoint;
      fixFSOrder(nhp);
      if( !duplicate(nhp, theProcesses) ) {
	assignToCF(nhp);
	theProcesses.push_back(nhp);
      }
    }
  }
}

void 
HardProcessConstructor::makeDiagrams(IDPair in, long out1, const PDSet & out2, 
				     PDPtr inter, HPDiagram::Channel chan, 
				     VBPair vertexpair, BPair cross) {
  for(PDSet::const_iterator it = out2.begin(); it != out2.end(); ++it) {
    HPDiagram nhp;
    nhp.incoming  = in;
    nhp.outgoing = make_pair(out1, (*it)->id());
    nhp.intermediate = inter;
    nhp.vertices = vertexpair;
    nhp.channelType = chan;
    nhp.ordered = cross;
    fixFSOrder(nhp);
    if( !duplicate(nhp, theProcesses) ) {
      assignToCF(nhp);
      theProcesses.push_back(nhp);
    }
  }

}

set<PDPtr> 
HardProcessConstructor::search(VBPtr vertex, long part1, direction d1, 
			       long part2, direction d2, direction d3) {
  if(vertex->getNpoint() != 3)
    return PDSet();
  if(d1 == incoming && getParticleData(part1)->CC()) part1 = -part1;
  if(d2 == incoming && getParticleData(part2)->CC()) part2 = -part2;
  PDVector ext;
  PDSet third;
  for(unsigned int ix = 0;ix < 3; ++ix) {
    PDVector pdlist = vertex->search(ix, part1);
    ext.insert(ext.end(), pdlist.begin(), pdlist.end());
  }
  for(unsigned int ix = 0; ix < ext.size(); ix += 3) {
    long id0 = ext.at(ix)->id();
    long id1 = ext.at(ix+1)->id();
    long id2 = ext.at(ix+2)->id();
    int pos;
    if((id0 == part1 && id1 == part2) ||
       (id0 == part2 && id1 == part1))
      pos = ix + 2;
    else if((id0 == part1 && id2 == part2) ||
	    (id0 == part2 && id2 == part1))
      pos = ix + 1;
    else if((id1 == part1 && id2 == part2) ||
	    (id1 == part2 && id2 == part1))
      pos = ix;
    else
      pos = -1;
    if(pos >= 0) {
      if(d3 == incoming && ext[pos]->CC()) ext[pos] = ext[pos]->CC();
      third.insert(ext[pos]);
    }
  }
  
  return third;
}

set<PDPtr>
HardProcessConstructor::search(VBPtr vertex, long part1, direction d1,
			       long part2, direction d2, long part3, direction d3,
			       direction d4) {
  if(vertex->getNpoint() != 4)
    return PDSet();
  if(d1 == incoming && getParticleData(part1)->CC()) part1 = -part1;
  if(d2 == incoming && getParticleData(part2)->CC()) part2 = -part2;
  if(d3 == incoming && getParticleData(part3)->CC()) part3 = -part3;
  PDVector ext;
  PDSet fourth;
  for(unsigned int ix = 0;ix < 4; ++ix) {
    PDVector pdlist = vertex->search(ix, part1);
    ext.insert(ext.end(), pdlist.begin(), pdlist.end());
  }
  for(unsigned int ix = 0;ix < ext.size(); ix += 4) {
    long id0 = ext.at(ix)->id(); long id1 = ext.at(ix + 1)->id();
    long id2 = ext.at(ix + 2)->id(); long id3 = ext.at(ix + 3)->id();
    int pos;
    if((id0 == part1 && id1 == part2 && id2 == part3) ||
       (id0 == part1 && id1 == part3 && id2 == part2) ||
       (id0 == part2 && id1 == part1 && id2 == part3) ||
       (id0 == part2 && id1 == part3 && id2 == part1) ||
       (id0 == part3 && id1 == part1 && id2 == part2) ||
       (id0 == part3 && id1 == part2 && id2 == part1))
      pos = ix + 3;
    else  if((id0 == part1 && id1 == part2 && id3 == part3) ||
	     (id0 == part1 && id1 == part3 && id3 == part2) ||
	     (id0 == part2 && id1 == part1 && id3 == part3) ||
	     (id0 == part2 && id1 == part3 && id3 == part1) ||
	     (id0 == part3 && id1 == part1 && id3 == part2) ||
	     (id0 == part3 && id1 == part2 && id3 == part1))
      pos = ix + 2;
    else if((id0 == part1 && id2 == part2 && id3 == part3) ||
	    (id0 == part1 && id2 == part3 && id3 == part2) ||
	    (id0 == part2 && id2 == part1 && id3 == part3) ||
	    (id0 == part2 && id2 == part3 && id3 == part1) ||
	    (id0 == part3 && id2 == part1 && id3 == part2) ||
	    (id0 == part3 && id2 == part2 && id3 == part1))
      pos = ix + 1;
    else if((id1 == part1 && id2 == part2 && id3 == part3) ||
	    (id1 == part1 && id2 == part3 && id3 == part2) ||
	    (id1 == part2 && id2 == part1 && id3 == part3) ||
	    (id1 == part2 && id2 == part3 && id3 == part1) ||
	    (id1 == part3 && id2 == part1 && id3 == part2) ||
	    (id1 == part3 && id2 == part2 && id3 == part1))
      pos = ix;
    else 
      pos = -1;
    
    if(pos >= 0) {
      if(d4 == incoming && ext[pos]->CC()) 
	ext[pos] = ext[pos]->CC();
      fourth.insert(ext[pos]);
    }
  } 
  return fourth;
}

void HardProcessConstructor::fixFSOrder(HPDiagram & diag) {
  tcPDPtr psa = getParticleData(diag.incoming.first);
  tcPDPtr psb = getParticleData(diag.incoming.second);
  tcPDPtr psc = getParticleData(diag.outgoing.first);
  tcPDPtr psd = getParticleData(diag.outgoing.second);

  //fix a spin order
  if( psc->iSpin() < psd->iSpin() ) {
    swap(diag.outgoing.first, diag.outgoing.second);
    if(diag.channelType == HPDiagram::tChannel) {
      diag.ordered.second = !diag.ordered.second;
    }
    return;
  }
  
  //for diagrams with different flavour incoming states
  if( psa->iSpin() == PDT::Spin1Half && psb->iSpin() == PDT::Spin1Half &&
      !sameQuarkFlavour(psa->id(), psb->id()) &&
      sameQuarkFlavour(psa->id(), psd->id()) ) {
    swap(diag.outgoing.first, diag.outgoing.second);
    if(diag.channelType == HPDiagram::tChannel) {
      diag.ordered.second = !diag.ordered.second;
    }
    return;
  }
  
  if(  psc->id() < 0 && psd->id() > 0 ) {
    swap(diag.outgoing.first, diag.outgoing.second);
    if(diag.channelType == HPDiagram::tChannel) {
      diag.ordered.second = !diag.ordered.second;
    }
    return;
  }

}

void HardProcessConstructor::assignToCF(HPDiagram & diag) {

  if(diag.channelType == HPDiagram::tChannel) {
    if(diag.ordered.second) tChannelCF(diag);
    else uChannelCF(diag);
  }
  else if(diag.channelType == HPDiagram::sChannel)
    sChannelCF(diag);
  else {
    vector<CFPair> cfv(2);
    cfv[0] = make_pair(1, 1);
    cfv[1] = make_pair(2, 1);
    diag.colourFlow = cfv;
  }
}

void HardProcessConstructor::tChannelCF(HPDiagram & diag) {
  vector<CFPair> cfv(1, make_pair(1, 1.));
  if(diag.intermediate->iColour() == PDT::Colour0) {
    long id1 = abs(diag.incoming.first);
    long id2 = abs(diag.incoming.second);
    long id3 = abs(diag.outgoing.first);
    long id4 = abs(diag.outgoing.second);
    if( getParticleData(id1)->iColour() == PDT::Colour3 && 
	getParticleData(id2)->iColour() == PDT::Colour3 &&
	getParticleData(id3)->iColour() == PDT::Colour3 && 
	getParticleData(id4)->iColour() == PDT::Colour3 ) {
      if( !sameQuarkFlavour(diag.incoming.first, diag.incoming.second) )
	cfv[0].first = 2;
      else
	cfv[0].first = 3;
    }
  }
  diag.colourFlow = cfv;
}
 
void HardProcessConstructor::uChannelCF(HPDiagram & diag) {
  PDT::Colour offshell = diag.intermediate->iColour();
  PDT::Colour outa = getParticleData(diag.outgoing.first)->iColour();
  PDT::Colour outb = getParticleData(diag.outgoing.second)->iColour();
  vector<CFPair> cfv(1, make_pair(2, 1.));
  if(offshell == PDT::Colour8 && (outa != outb) ) {
    cfv[0].first = 1;
    cfv.push_back(make_pair(2, -1.));
  }
  else {
    long id1 = abs(diag.incoming.first);
    long id2 = abs(diag.incoming.second);
    if( getParticleData(id1)->iColour() == PDT::Colour3 && 
	getParticleData(id2)->iColour() == PDT::Colour3 &&
	outa == PDT::Colour3 && outb == PDT::Colour3 ) {
      if( offshell == PDT::Colour0 ) {
	if( sameQuarkFlavour(id1, id2) ) 
	  cfv[0].first = 4;
	else
	  cfv[0].first = 3;
      }
    }
    if( outa == PDT::Colour0 || outb == PDT::Colour0 )
      cfv[0].first = 1;
  }
  diag.colourFlow = cfv;
}

void HardProcessConstructor::sChannelCF(HPDiagram & diag) {
  tcPDPtr pa = getParticleData(diag.incoming.first);
  tcPDPtr pb = getParticleData(diag.incoming.second);
  PDT::Colour ina = pa->iColour();
  PDT::Colour inb = pb->iColour();
  PDT::Colour offshell = diag.intermediate->iColour();
  tcPDPtr pc = getParticleData(diag.outgoing.first);
  tcPDPtr pd = getParticleData(diag.outgoing.second);
  PDT::Colour outa = pc->iColour();
  PDT::Colour outb = pd->iColour();

  vector<CFPair> cfv(1);
  if(offshell == PDT::Colour8) {
    if( (ina == PDT::Colour8 && inb == PDT::Colour8) ||
	(outa == PDT::Colour8 && outb == PDT::Colour8) ) { 
      //Require an additional minus sign for a fermion 33bar final state
      //due to the way the vertex rules are defined.
      int prefact(1);
      if( (pc->iSpin() == PDT::Spin1Half && pd->iSpin() == PDT::Spin1Half) &&
	  (pc->iColour() == PDT::Colour3 && pd->iColour() == PDT::Colour3bar) )
	prefact = -1;
      
      cfv[0].first = 1;
      cfv[0].second = -prefact;
      cfv.push_back(make_pair(2, prefact));
    }
    else if( !sameQuarkFlavour(diag.incoming.first, diag.outgoing.first) )
      cfv[0] = make_pair(1, 1);
    else
      cfv[0] = make_pair(2, 1);
  }
  else if(offshell == PDT::Colour0) {
    if( ina == PDT::Colour0 || inb == PDT::Colour0 ||
	outa == PDT::Colour0 || outb == PDT::Colour0 )
      cfv[0] = make_pair(1, 1);
    else {
      if( sameQuarkFlavour(diag.incoming.first, diag.incoming.second) ) {
	if( sameQuarkFlavour(diag.incoming.first, diag.outgoing.first) )
	  cfv[0] = make_pair(4, 1);
	else
	  cfv[0] = make_pair(2, 1);
      }
      else
	cfv[0] = make_pair(3, 1);
    }
  }
  else {
    if(outa == PDT::Colour0)
      cfv[0] = make_pair(1, 1);
    else
      cfv[0] = make_pair(2, 1);
  }
  
  diag.colourFlow = cfv; 
}

void 
HardProcessConstructor::createMatrixElement(const HPDVector & process) const {
  for(HPDVector::size_type d = 0; d < process.size(); ++d) {
    HPDiagram diag = process[d];
     cout << getParticleData(diag.incoming.first)->PDGName() << ","
	 << getParticleData(diag.incoming.second)->PDGName() << "->";
    if(diag.intermediate)
      cout << diag.intermediate->PDGName() << "->";
    
    cout << getParticleData(diag.outgoing.first)->PDGName() << ","
	 << getParticleData(diag.outgoing.second)->PDGName()
	 << "  channel " << diag.channelType;
    if(diag.channelType == HPDiagram::tChannel) {
      cout << "  ordering " << diag.ordered.first << " " 
	   << diag.ordered.second << "   ";
      for(unsigned int cf = 0; cf < diag.colourFlow.size(); ++cf) 
	cout << "(" << diag.colourFlow[cf].first << "," 
	     <<diag.colourFlow[cf].second << ")  ";
      cout << "\n\n";
    }
    else {
      cout << "   ";
      for(unsigned int cf = 0; cf < diag.colourFlow.size(); ++cf) 
	cout << "(" << diag.colourFlow[cf].first << "," 
	     <<diag.colourFlow[cf].second << ")  ";
      cout << "\n\n";
    }
  }
  cout << "---------------------------\n";
  cout << flush;
  tcPDVector extpart(4);
  extpart[0] = getParticleData(process[0].incoming.first);
  extpart[1] = getParticleData(process[0].incoming.second);
  extpart[2] = getParticleData(process[0].outgoing.first);
  extpart[3] = getParticleData(process[0].outgoing.second);
  string objectname ("/Defaults/MatrixElements/");
  string classname = MEClassname(extpart, objectname);
  GeneralHardMEPtr matrixElement = dynamic_ptr_cast<GeneralHardMEPtr>
      (generator()->preinitCreate(classname, objectname));
  if(matrixElement) {
    unsigned int ncf(0);
    vector<DVector> cfactors = getColourFactors(extpart, ncf);
    matrixElement->setProcessInfo(process, cfactors, ncf);
    generator()->preinitInterface(theSubProcess, "MatrixElements", 
				  theSubProcess->MEs().size(),
				  "insert", matrixElement->fullName()); 
  }
  else 
    throw HardProcessConstructorError() 
      << "createMatrixElement - No matrix element object could be created for "
      << "the process " 
      << extpart[0]->PDGName() << " " << extpart[0]->iSpin() << "," 
      << extpart[1]->PDGName() << " " << extpart[1]->iSpin() << "->" 
      << extpart[2]->PDGName() << " " << extpart[2]->iSpin() << "," 
      << extpart[3]->PDGName() << " " << extpart[3]->iSpin() 
      << ".  No class for this spin-structure exists! \n"
      << Exception::warning;
}

vector<DVector> HardProcessConstructor::
getColourFactors(const tcPDVector & extpart, unsigned int & ncf) const {
  PDT::Colour ina = extpart[0]->iColour();
  PDT::Colour inb = extpart[1]->iColour();
  PDT::Colour outa = extpart[2]->iColour();
  PDT::Colour outb = extpart[3]->iColour();

  vector<DVector> scf(1, DVector(1, 1.));
  if( ina == PDT::Colour0 || inb == PDT::Colour0 ||
      outa == PDT::Colour0 || outb == PDT::Colour0 ) {
    ncf = 1;
    if(outa == PDT::Colour8 || outb == PDT::Colour8)
      scf[0][0] = 4.;
    else
      scf[0][0] = 3.;
  }
  else if( ina == PDT::Colour8 || inb == PDT::Colour8 ||
	   outa == PDT::Colour8 || outb == PDT::Colour8 ) {
    ncf = 2;
    if( ina == PDT::Colour8 && inb == PDT::Colour8 &&
        outa == PDT::Colour8 && outb == PDT::Colour8 )
      return the88to88;
    else
      return the33bto88;
  }
  else {
    if( sameQuarkFlavour(extpart[0]->id(), extpart[1]->id()) &&
	sameQuarkFlavour(extpart[0]->id(), extpart[2]->id()) ) {
      if(theAllDiagrams) ncf = 4;
      else ncf = 2;
      return the33bto33b;
    }
    else {
      if(theAllDiagrams) ncf = 3;
      else ncf = 1;
      return the33bpto33bp;
    }
      
  }
  return scf;
}

bool HardProcessConstructor::duplicate(const HPDiagram & diag, 
				       const HPDVector & group) const {
  unsigned int Nd = group.size(), ix(0);
  if(Nd == 0) return false;
  bool copy(false);
  do {
    if( group[ix] == diag ) copy = true;
    ++ix;
  }
  while(copy == false && ix < Nd);
  return copy;
} 

string HardProcessConstructor::MEClassname(const vector<tcPDPtr> & extpart, 
					   string & objname) const {
 string classname("/Herwig++/ME");
  for(vector<tcPDPtr>::size_type ix = 0; ix < extpart.size(); ++ix) {
    if(ix == 2) classname += "2";
    if(extpart[ix]->iSpin() == PDT::Spin0) classname += "s";
    else if(extpart[ix]->iSpin() == PDT::Spin1) classname += "v";
    else if(extpart[ix]->iSpin() == PDT::Spin1Half) classname += "f";
    else if(extpart[ix]->iSpin() == PDT::Spin2) classname += "t";
    else
      throw HardProcessConstructorError() 
	<< "MEClassname() : Encountered an unknown spin for "
	<< extpart[ix]->PDGName() << " while constructing MatrixElement "
	<< "classname " << extpart[ix]->iSpin() << Exception::warning;
  }
  objname += "ME" + extpart[0]->PDGName() + extpart[1]->PDGName() + "2" 
    + extpart[2]->PDGName() + extpart[3]->PDGName();
  return classname;  
}

bool HardProcessConstructor::sameQuarkFlavour(long id1, long id2) const {
  //if ids are not quarks or quark partners return false
  if( getParticleData(abs(id1))->iColour() != PDT::Colour3 && 
      getParticleData(abs(id2))->iColour() != PDT::Colour3 ) return false;
  long diff = abs(abs(id2) - abs(id1));
  if(diff == 0 || diff % 10 == 0) return true;
  return false;
}
