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
#include "Herwig++/MatrixElement/General/GeneralHardME.h"

using namespace Herwig;

IBPtr ResonantProcessConstructor::clone() const {
  return new_ptr(*this);
}

IBPtr ResonantProcessConstructor::fullclone() const {
  return new_ptr(*this);
}

void ResonantProcessConstructor::doinit() {
  Interfaced::doinit();
  EGPtr eg = generator();
  theModel = dynamic_ptr_cast<HwSMPtr>(eg->standardModel());
  if(!theModel)
    throw InitException() << "ResonantProcessConstructor:: doinit() - "
			  << "The model pointer is null!"
			  << Exception::abortnow;
  if(eg->eventHandler()) {
    string subProcessName = 
      eg->preinitInterface(eg->eventHandler(), "SubProcessHandlers", "get","");
   theSubProcess = eg->getObject<SubProcessHandler>(subProcessName);
   if(!theSubProcess)
     throw InitException() << "ResonantProcessConstructor:: doinit() - "
			   << "There was an error getting the SubProcessHandler "
			   << "from the current event handler. "
			   << Exception::abortnow;
  }
  else
    throw
      InitException() << "ResonantProcessConstructor:: doinit() - "
		      << "The eventHandler pointer was null therefore "
		      << "could not get SubProcessHandler pointer " 
		      << Exception::abortnow;
}

void ResonantProcessConstructor::persistentOutput(PersistentOStream & os) const {
  os << theIncoming << theIntermediates << theOutgoing << theModel 
     << theSubProcess;
}

void ResonantProcessConstructor::persistentInput(PersistentIStream & is, int) {
  is >> theIncoming >> theIntermediates >> theOutgoing >> theModel 
     >> theSubProcess;
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

  static Switch<ResonantProcessConstructor,bool> interfaceDebugME
    ("DebugME",
     "Print comparison with analytical ME",
     &ResonantProcessConstructor::theDebug, false, false, false);
  static SwitchOption interfaceDebugMEYes
    (interfaceDebugME,
     "Yes",
     "Print the debug information",
     true);
  static SwitchOption interfaceDebugMENo
    (interfaceDebugME,
     "No",
     "Do not print the debug information",
     false);
}

void ResonantProcessConstructor::constructResonances() {
  size_t ninc = theIncoming.size() , ninter = theIntermediates.size();
  if(ninc == 0 || ninter == 0 || theOutgoing.size() == 0) return;
  //Remove matrix elements already present
  EGPtr eg = generator();
  int nme = theSubProcess->MEs().size();
  for(int ix = nme - 1; ix >= 0; --ix)
    eg->preinitInterface(theSubProcess, "MatrixElements", ix, "erase", "");
    // Get a pointer to the model in use
  theModel->init();
  size_t nvertices = theModel->numberOfVertices(); 
  //init the vertices
  for(size_t iv = 0; iv < nvertices; ++iv)
    theModel->vertex(iv)->init();
  //To construct resonant diagrams loop over the incoming particles, intermediates
  //and vertices to find allowed diagrams. Need to exclude the diagrams that have 
  //the intermediate as an external particle aswell
  for(vector<PDPtr>::size_type ii = 0; ii < ninc; ++ii) {
    tcPDPtr parta = theIncoming[ii];
    for(vector<PDPtr>::size_type ij = ii; ij < ninc; ++ij) {
      tcPDPtr partb = theIncoming[ij];
      for(vector<PDPtr>::size_type ik = 0; ik < ninter ; ++ik) {
	long ipart = theIntermediates[ik]->id();
	for(size_t iv = 0; iv < nvertices; ++iv) {
	  VBPtr vertex = theModel->vertex(iv);
	  if(vertex->getNpoint() > 3) continue;
	  long part1 = parta->CC() ? -parta->id() : parta->id();
	  long part2 = partb->CC() ? -partb->id() : partb->id();
	  if(vertex->allowed(part1, part2, ipart) || 
	     vertex->allowed(part1, ipart, part2) ||
	     vertex->allowed(part2, part1, ipart) ||
	     vertex->allowed(part2, ipart, part1) ||
	     vertex->allowed(ipart, part1, part2) ||
	     vertex->allowed(ipart, part2, part1) ) {
	    constructVertex2(make_pair(parta->id(), partb->id()), vertex, 
			     theIntermediates[ik]);
	  }
	}
      }
    }
  }
//   cout << "Created " << theDiagrams.size() << " resonant processes\n";
//   for(HPDVector::size_type ix = 0; ix < theDiagrams.size(); ++ix)
//     cout << getParticleData(theDiagrams[ix].incoming.first)->PDGName() << ","
// 	 << getParticleData(theDiagrams[ix].incoming.second)->PDGName()  << "->"
// 	 << theDiagrams[ix].intermediate->PDGName() << "->" 
// 	 << getParticleData(theDiagrams[ix].outgoing.first)->PDGName() << ","
// 	 << getParticleData(theDiagrams[ix].outgoing.second)->PDGName()  << "\n";
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
  size_t nvertices = theModel->numberOfVertices(); 
  for(size_t io = 0; io < theOutgoing.size(); ++io) {
    tcPDPtr outa = theOutgoing[io];
    for(size_t iv = 0; iv < nvertices; ++iv) {
      VBPtr vertex2 = theModel->vertex(iv);
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
    newdiag.colourFlow = vector<CFPair>(1, make_pair(1, 1.));
    theDiagrams.push_back(newdiag);
  }
}
	
set<tPDPtr> 
ResonantProcessConstructor::search(VBPtr vertex, long part1, direction d1, 
				   long part2, direction d2, direction d3) {
  if(d1 == incoming && getParticleData(part1)->CC()) part1 = -part1;
  if(d2 == incoming && getParticleData(part2)->CC()) part2 = -part2;
  tPDVector ext;
  tPDSet third;
  for(unsigned int ix = 0;ix < 3; ++ix) {
    tPDVector pdlist = vertex->search(ix, part1);
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
      if(d3 == incoming && ext.at(pos)->CC()) ext.at(pos) = ext.at(pos)->CC();
      third.insert(ext.at(pos));
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
  string classname = MEClassname(extpart, objectname);
  GeneralHardMEPtr matrixElement = dynamic_ptr_cast<GeneralHardMEPtr>
    (generator()->preinitCreate(classname, objectname));
  if( !matrixElement ) {
    throw RPConstructorError() 
      << "createMatrixElement - No matrix element object could be created for "
      << "the process " 
      << extpart[0]->PDGName() << " " << extpart[0]->iSpin() << "," 
      << extpart[1]->PDGName() << " " << extpart[1]->iSpin() << "->" 
      << extpart[2]->PDGName() << " " << extpart[2]->iSpin() << "," 
      << extpart[3]->PDGName() << " " << extpart[3]->iSpin() 
      << ".  Constructed class name: \"" << classname << "\""
      << Exception::warning;
    return;
  }
  matrixElement->setProcessInfo(HPDVector(1, diag),
				colourFactor(extpart), 1, theDebug,1);
  generator()->preinitInterface(theSubProcess, "MatrixElements", 
				theSubProcess->MEs().size(),
				"insert", matrixElement->fullName()); 
}

string ResonantProcessConstructor::MEClassname(const vector<tcPDPtr> & extpart, 
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
    + extpart[2]->PDGName() + extpart[3]->PDGName();
  return classname;  
}

vector<DVector> ResonantProcessConstructor::
colourFactor(const tcPDVector & extpart) const {
  vector<DVector > cfactor(1, DVector(1, 0.));

  if(extpart[0]->iColour() == PDT::Colour3 ||
     extpart[1]->iColour() == PDT::Colour3) {
    if(extpart[2]->iColour() == PDT::Colour3 || 
       extpart[3]->iColour() == PDT::Colour3)
      cfactor[0][0] = 9.;
    else if(extpart[2]->iColour() == PDT::Colour8 || 
	    extpart[3]->iColour() == PDT::Colour8)
      cfactor[0][0] = 24.;
    else
      cfactor[0][0] = 3.;
  }
  else if(extpart[0]->iColour() == PDT::Colour8) {
    if(extpart[2]->iColour() == PDT::Colour3 || 
       extpart[3]->iColour() == PDT::Colour3)
      cfactor[0][0] = 24.;
    else if(extpart[2]->iColour() == PDT::Colour8 || 
	    extpart[3]->iColour() == PDT::Colour8)
      cfactor[0][0] = 64.;
    else 
      cfactor[0][0] = 8.;
  }
  else if(extpart[0]->iColour() == PDT::Colour0) {
    if(extpart[2]->iColour() == PDT::Colour3 || 
       extpart[3]->iColour() == PDT::Colour3)
      cfactor[0][0] = 3.;
    else if(extpart[2]->iColour() == PDT::Colour8 || 
	    extpart[3]->iColour() == PDT::Colour8)
      cfactor[0][0] = 8.;
    else
      cfactor[0][0] = 1.;
  }
  return cfactor;
}
