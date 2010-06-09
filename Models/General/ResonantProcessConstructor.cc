// -*- C++ -*-
//
// ResonantProcessConstructor.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ResonantProcessConstructor class.
//

#include "ResonantProcessConstructor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Switch.h"

using namespace Herwig;

IBPtr ResonantProcessConstructor::clone() const {
  return new_ptr(*this);
}

IBPtr ResonantProcessConstructor::fullclone() const {
  return new_ptr(*this);
}

void ResonantProcessConstructor::persistentOutput(PersistentOStream & os) const {
  os << theIncoming << theIntermediates << theOutgoing;
}

void ResonantProcessConstructor::persistentInput(PersistentIStream & is, int) {
  is >> theIncoming >> theIntermediates >> theOutgoing;
}

ClassDescription<ResonantProcessConstructor> 
ResonantProcessConstructor::initResonantProcessConstructor;
// Definition of the static class description member.

void ResonantProcessConstructor::Init() {

  static ClassDocumentation<ResonantProcessConstructor> documentation
    ("This class is designed solely to contruct resonant processes using"
     "a provided set of intermediate particles");
  
  static RefVector<ResonantProcessConstructor, ParticleData> interfaceOffshell
    ("Intermediates",
     "A vector of offshell particles for resonant diagrams",
     &ResonantProcessConstructor::theIntermediates, -1, false, false, true, 
     false);

  static RefVector<ResonantProcessConstructor, ParticleData> interfaceIncoming
    ("Incoming",
     "A vector of incoming particles for resonant diagrams",
     &ResonantProcessConstructor::theIncoming, -1, false, false, true, 
     false);

  static RefVector<ResonantProcessConstructor, ParticleData> interfaceOutgoing
    ("Outgoing",
     "A vector of outgoin particles for resonant diagrams",
     &ResonantProcessConstructor::theOutgoing, -1, false, false, true, 
     false);
}

namespace {
  // Helper functor for find_if in duplicate function.
  class SameIncomingAs {
  public:
    SameIncomingAs(tPDPair in) : a(in.first->id()), b(in.second->id())  {}
    bool operator()(tPDPair ppair) const {
      long id1(ppair.first->id()), id2(ppair.second->id());
      return ( id1 == a && id2 == b ) || ( id1 == b && id2 == a );
    }
  private:
    long a, b;
  };

  bool duplicateIncoming(tPDPair ppair,const vector<tPDPair> &incPairs) {
    vector<tPDPair>::const_iterator it = 
      find_if( incPairs.begin(), incPairs.end(), SameIncomingAs(ppair) );
    return it != incPairs.end(); 
  }
}

void ResonantProcessConstructor::constructDiagrams() {
  size_t ninc = theIncoming.size() , ninter = theIntermediates.size();
  if(ninc == 0 || ninter == 0 || theOutgoing.size() == 0) return;
  // find the incoming particle pairs
  vector<tPDPair> incPairs;
  for(PDVector::size_type i = 0; i < ninc; ++i) {
    for(PDVector::size_type j = 0; j < ninc; ++j) {
      tPDPair inc = make_pair(theIncoming[i], theIncoming[j]);
      if( (inc.first->iSpin() > inc.second->iSpin()) ||
	  (inc.first->iSpin() == inc.second->iSpin() &&
	   inc.first->id() < inc.second->id()) )
	swap(inc.first, inc.second);
      if( !duplicateIncoming(inc,incPairs) ) {
	incPairs.push_back(inc);
      }
    }
  }
  //Remove matrix elements already present
  EGPtr eg = generator();
  int nme = subProcess()->MEs().size();
  for(int ix = nme - 1; ix >= 0; --ix)
    eg->preinitInterface(subProcess(), "MatrixElements", ix, "erase", "");
  // Get a pointer to the model in use
  model()->init();
  size_t nvertices = model()->numberOfVertices(); 
  //init the vertices
  for(size_t iv = 0; iv < nvertices; ++iv)
    model()->vertex(iv)->init();
  //To construct resonant diagrams loop over the incoming particles, intermediates
  //and vertices to find allowed diagrams. Need to exclude the diagrams that have 
  //the intermediate as an external particle aswell
  for(vector<tcPDPair>::size_type is = 0; is < incPairs.size(); ++is) {
    tPDPair ppi = incPairs[is]; 
    for(vector<PDPtr>::size_type ik = 0; ik < ninter ; ++ik) {
      long ipart = theIntermediates[ik]->id();
      for(size_t iv = 0; iv < nvertices; ++iv) {
	VBPtr vertex = model()->vertex(iv);
	if(vertex->getNpoint() > 3) continue;
	long part1 = ppi.first->CC() ? -ppi.first->id() : ppi.first->id();
	long part2 = ppi.second->CC() ? -ppi.second->id() : ppi.second->id();
	if(vertex->allowed(part1, part2, ipart) || 
	   vertex->allowed(part1, ipart, part2) ||
	   vertex->allowed(part2, part1, ipart) ||
	   vertex->allowed(part2, ipart, part1) ||
	   vertex->allowed(ipart, part1, part2) ||
	   vertex->allowed(ipart, part2, part1) ) {
	  constructVertex2(make_pair(ppi.first->id(), ppi.second->id()), vertex, 
			   theIntermediates[ik]);
	}
      }
    }
  }
  //Create matrix element for each process
  const HPDVector::size_type ndiags = theDiagrams.size();
  for(HPDVector::size_type ix = 0; ix < ndiags; ++ix)
    createMatrixElement(theDiagrams[ix]);
}

void ResonantProcessConstructor::
constructVertex2(IDPair in, VertexBasePtr vertex, 
		 PDPtr partc) {
  //We have the left hand part of the diagram, just need all the possibilities
  //for the RHS
  size_t nvertices = model()->numberOfVertices(); 
  for(size_t io = 0; io < theOutgoing.size(); ++io) {
    tcPDPtr outa = theOutgoing[io];
    for(size_t iv = 0; iv < nvertices; ++iv) {
      VBPtr vertex2 = model()->vertex(iv);
      if(vertex2->getNpoint() > 3) continue;
      tPDSet outb = search(vertex2, partc->id(), incoming, outa->id(), outgoing,
			  outgoing);
      makeResonantDiagrams(in, partc, outa->id(), outb, 
			   make_pair(vertex, vertex2));
    }
  }
}

void ResonantProcessConstructor::
makeResonantDiagrams(IDPair in, PDPtr offshell, long outa, const tPDSet & out, 
		     VBPair vertpair) {
  for(tPDSet::const_iterator ita = out.begin(); ita != out.end(); ++ita) {
    if( abs(outa) == abs(offshell->id()) || 
	abs((*ita)->id()) == abs(offshell->id())) continue;
    HPDiagram newdiag(in,make_pair(outa, (*ita)->id()) );
    newdiag.intermediate = offshell;
    newdiag.vertices = vertpair;
    newdiag.channelType = HPDiagram::sChannel;
    fixFSOrder(newdiag);
    assignToCF(newdiag);
    if(!duplicate(newdiag,theDiagrams))
      theDiagrams.push_back(newdiag);
  }
}
	
set<tPDPtr> 
ResonantProcessConstructor::search(VBPtr vertex, long part1, direction d1, 
				   long part2, direction d2, direction d3) {
  if(d1 == incoming && getParticleData(part1)->CC()) part1 = -part1;
  if(d2 == incoming && getParticleData(part2)->CC()) part2 = -part2;
  vector<long> ext;
  tPDSet third;
  for(unsigned int ix = 0;ix < 3; ++ix) {
    vector<long> pdlist = vertex->search(ix, part1);
    ext.insert(ext.end(), pdlist.begin(), pdlist.end());
  }
  for(unsigned int ix = 0; ix < ext.size(); ix += 3) {
    long id0 = ext.at(ix);
    long id1 = ext.at(ix+1);
    long id2 = ext.at(ix+2);
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
      tPDPtr p = getParticleData(ext[pos]);
      if(d3 == incoming && p->CC()) p = p->CC();
      third.insert(p);
    }
  }
  
  return third;
}
						
IDPair ResonantProcessConstructor::
find(long part, const vector<PDPtr> & out) const {
  vector<PDPtr>::size_type iloc(0);
  bool found(false);
  do {
    if(out[iloc]->id() == part) found = true;
    else ++iloc;
  }
  while(found == false && iloc < out.size());
  //found offshell
  IDPair outids;
  if(iloc == 0)
    outids = make_pair(out[1]->id(), out[2]->id());
  else if(iloc == 1)
    outids = make_pair(out[0]->id(), out[2]->id());
  else
    outids = make_pair(out[0]->id(), out[1]->id());
  return outids;
} 

void ResonantProcessConstructor::
createMatrixElement(const HPDiagram & diag) const {
  vector<tcPDPtr> extpart(4);
  extpart[0] = getParticleData(diag.incoming.first);
  extpart[1] = getParticleData(diag.incoming.second);
  extpart[2] = getParticleData(diag.outgoing.first);
  extpart[3] = getParticleData(diag.outgoing.second);
  string objectname ("/Herwig/MatrixElements/");
  string classname = MEClassname(extpart, diag.intermediate, objectname);
  GeneralHardMEPtr matrixElement = dynamic_ptr_cast<GeneralHardMEPtr>
    (generator()->preinitCreate(classname, objectname));
  if( !matrixElement ) {
    throw RPConstructorError() 
      << "ResonantProcessConstructor::createMatrixElement "
      << "- No matrix element object could be created for "
      << "the process " 
      << extpart[0]->PDGName() << " " << extpart[0]->iSpin() << "," 
      << extpart[1]->PDGName() << " " << extpart[1]->iSpin() << "->" 
      << extpart[2]->PDGName() << " " << extpart[2]->iSpin() << "," 
      << extpart[3]->PDGName() << " " << extpart[3]->iSpin() 
      << ".  Constructed class name: \"" << classname << "\"\n"
      << Exception::warning;
    return;
  }
  matrixElement->setProcessInfo(HPDVector(1, diag),
				colourFlow(extpart), debug(),1);
  generator()->preinitInterface(subProcess(), "MatrixElements", 
				subProcess()->MEs().size(),
				"insert", matrixElement->fullName()); 
}

string ResonantProcessConstructor::
MEClassname(const vector<tcPDPtr> & extpart, tcPDPtr inter,
	    string & objname) const {
 string classname("Herwig::ME");
  for(vector<tcPDPtr>::size_type ix = 0; ix < extpart.size(); ++ix) {
    if(ix == 2) classname += "2";
    if(extpart[ix]->iSpin() == PDT::Spin0) classname += "s";
    else if(extpart[ix]->iSpin() == PDT::Spin1) classname += "v";
    else if(extpart[ix]->iSpin() == PDT::Spin1Half) classname += "f";
    else if(extpart[ix]->iSpin() == PDT::Spin2) classname += "t";
    else
      throw RPConstructorError()
	<< "MEClassname() : Encountered an unknown spin for "
	<< extpart[ix]->PDGName() << " while constructing MatrixElement "
	<< "classname " << extpart[ix]->iSpin() << Exception::warning;
  }
  objname += "ME" + extpart[0]->PDGName() + extpart[1]->PDGName() + "2"
    + inter->PDGName() + "2"
    + extpart[2]->PDGName() + extpart[3]->PDGName();
  return classname;  
}
