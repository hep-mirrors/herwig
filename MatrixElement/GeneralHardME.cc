// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GeneralHardME class.
//

#include "GeneralHardME.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"

using namespace Herwig;

GeneralHardME::GeneralHardME() : theIncoming(0, 0), theOutgoing(0, 0),
				 theDiagrams(0), theNDiags(0), 
				 theColour(0), theNcf(0) , 
				 theDebug(false) {}
  
void GeneralHardME::setProcessInfo(const vector<HPDiagram> & alldiagrams,
				   const vector<vector<double> > & factors,
				   const unsigned int ncf,
				   bool debug) {
  theIncoming = alldiagrams.at(0).incoming;
  theOutgoing = alldiagrams.at(0).outgoing;
  theDiagrams = alldiagrams;
  theColour = factors;
  theNDiags = alldiagrams.size();
  theNcf = ncf;
  theDebug = debug;
}


void GeneralHardME::getDiagrams() const {
  //get ParticleData pointers for external particles
  tcPDPtr ina = getParticleData(getIncoming().first);
  tcPDPtr inb = getParticleData(getIncoming().second);
  tcPDPtr outa = getParticleData(getOutgoing().first);
  tcPDPtr outb = getParticleData(getOutgoing().second);
  
  for(HPCount idx = 0; idx < theNDiags; ++idx) {
    HPDiagram current = getProcessInfo()[idx];
    tcPDPtr offshell = current.intermediate;
    if(!offshell) continue;
    //t-channel
    if(current.channelType == HPDiagram::tChannel) {
      if(offshell->id() < 0) offshell = offshell->CC();
      if(current.ordered.second)
	add(new_ptr((Tree2toNDiagram(3), ina, offshell,
		     inb, 1, outa, 2, outb, -(idx+1))));
      else 
	add(new_ptr((Tree2toNDiagram(3), ina, offshell,
		     inb, 2, outa, 1, outb, -(idx+1))));
    }
    //s-channel
    else if(current.channelType == HPDiagram::sChannel) 
      add(new_ptr((Tree2toNDiagram(2), ina, inb, 1, offshell,
		   3, outa, 3, outb, -(idx+1))));
    else
      throw MEException() << "getDiagrams() - Unknown diagram in matrix element "
			  << fullName() << Exception::runerror;			  
  }
}

unsigned int GeneralHardME::orderInAlphaS() const {
  unsigned int order(0);
  for(HPCount idx = 0; idx < theNDiags; ++idx) {
    unsigned int tOrder = theDiagrams[idx].vertices.first->orderInGs() + 
      theDiagrams[idx].vertices.second->orderInGs();
    if(tOrder > order) order = tOrder;
  }
  return order;
}

unsigned int GeneralHardME::orderInAlphaEW() const {
  unsigned int order(0);
  for(HPCount idx = 0; idx < theNDiags; ++idx) {
    unsigned int tOrder = theDiagrams[idx].vertices.first->orderInGem() + 
      theDiagrams[idx].vertices.second->orderInGem();
    if(tOrder > order) order = tOrder;
  }
  return order;
}

Selector<MEBase::DiagramIndex>
GeneralHardME::diagrams(const DiagramVector & diags) const {
  vector<double> singleME;
  if(lastXCombPtr()) {
    singleME.assign(meInfo().begin(), meInfo().end());
  }
  else {
    throw MEException() << fullName() << "::diagrams - Null lastXCombPtr! "
			<< Exception::runerror;
  } 
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i )
    sel.insert(singleME.at(abs(diags[i]->id()) - 1), i);
  return sel;
}

void GeneralHardME::persistentOutput(PersistentOStream & os) const {
  os << theIncoming << theOutgoing << theDiagrams << theColour 
     << theNDiags << theNcf << theDebug;
}

void GeneralHardME::persistentInput(PersistentIStream & is, int) {
  is >> theIncoming >> theOutgoing >> theDiagrams >> theColour 
     >> theNDiags >> theNcf >> theDebug;
}

AbstractClassDescription<GeneralHardME> GeneralHardME::initGeneralHardME;
// Definition of the static class description member.

void GeneralHardME::Init() {

  static ClassDocumentation<GeneralHardME> documentation
    ("This class is designed to be a base class for a specific spin "
     "configuration where no matrix element exists, i.e. when processes "
     "are created automaticlly for a different model.");

}

void GeneralHardME::debug(double) const {
}
