// -*- C++ -*-
//
// ResonantProcessConstructor.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
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
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"

using namespace Herwig;

IBPtr ResonantProcessConstructor::clone() const {
  return new_ptr(*this);
}

IBPtr ResonantProcessConstructor::fullclone() const {
  return new_ptr(*this);
}

void ResonantProcessConstructor::persistentOutput(PersistentOStream & os) const {
  os << incoming_ << intermediates_ << outgoing_ 
     << processOption_ << scaleChoice_ << scaleFactor_;
}

void ResonantProcessConstructor::persistentInput(PersistentIStream & is, int) {
  is >> incoming_ >> intermediates_ >> outgoing_ 
     >> processOption_ >> scaleChoice_ >> scaleFactor_;
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
     &ResonantProcessConstructor::intermediates_, -1, false, false, true, 
     false);

  static RefVector<ResonantProcessConstructor, ParticleData> interfaceIncoming
    ("Incoming",
     "A vector of incoming particles for resonant diagrams",
     &ResonantProcessConstructor::incoming_, -1, false, false, true, 
     false);

  static RefVector<ResonantProcessConstructor, ParticleData> interfaceOutgoing
    ("Outgoing",
     "A vector of outgoing particles for resonant diagrams",
     &ResonantProcessConstructor::outgoing_, -1, false, false, true, 
     false);

  static Switch<ResonantProcessConstructor,unsigned int> interfaceProcesses
    ("Processes",
     "Whether to generate inclusive or exclusive processes",
     &ResonantProcessConstructor::processOption_, 0, false, false);
  static SwitchOption interfaceProcessesSingleParticleInclusive
    (interfaceProcesses,
     "SingleParticleInclusive",
     "Require at least one particle from the list of outgoing particles"
     " in the hard process",
     0);
  static SwitchOption interfaceProcessesTwoParticleInclusive
    (interfaceProcesses,
     "TwoParticleInclusive",
     "Require that both the particles in the hard processes are in the"
     " list of outgoing particles",
     1);
  static SwitchOption interfaceProcessesExclusive
    (interfaceProcesses,
     "Exclusive",
     "Require that both the particles in the hard processes are in the"
     " list of outgoing particles in every hard process",
     2);
  static SwitchOption interfaceProcessesInclusive
    (interfaceProcesses,
     "Inclusive",
     "Generate all modes which are allowed for the on-shell intermediate particle",
     3);

  static Switch<ResonantProcessConstructor,unsigned int> interfaceScaleChoice
    ("ScaleChoice",
     "&ResonantProcessConstructor::scaleChoice_",
     &ResonantProcessConstructor::scaleChoice_, 1, false, false);
  static SwitchOption interfaceScaleChoicesHat
    (interfaceScaleChoice,
     "sHat",
     "Always use sHat",
     1);
  static SwitchOption interfaceScaleChoiceTransverseMass
    (interfaceScaleChoice,
     "TransverseMass",
     "Always use the transverse mass",
     2);
  static SwitchOption interfaceScaleChoiceGeometicMean
    (interfaceScaleChoice,
     "MaxMT",
     "Use the maximum of m^2+p_T^2 for the two particles",
     3);

  static Parameter<ResonantProcessConstructor,double> interfaceScaleFactor
    ("ScaleFactor",
     "The prefactor used in the scale calculation. The scale used is"
     " sHat multiplied by this prefactor",
     &ResonantProcessConstructor::scaleFactor_, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

}

void ResonantProcessConstructor::doinit() {
  HardProcessConstructor::doinit();
  if(processOption_==2&&outgoing_.size()!=2)
    throw InitException() 
      << "Exclusive processes require exactly"
      << " two outgoing particles but " << outgoing_.size()
      << "have been inserted in ResonantProcessConstructor::doinit()."
      << Exception::runerror;
}

void ResonantProcessConstructor::constructDiagrams() {
  size_t ninc = incoming_.size() , ninter = intermediates_.size();
  if(ninc == 0 || ninter == 0  || !subProcess() ) return;
  // find the incoming particle pairs
  vector<tPDPair> incPairs;
  for(PDVector::size_type i = 0; i < ninc; ++i) {
    for(PDVector::size_type j = 0; j < ninc; ++j) {
      tPDPair inc = make_pair(incoming_[i], incoming_[j]);
      if( (inc.first->iSpin() > inc.second->iSpin()) ||
	  (inc.first->iSpin() == inc.second->iSpin() &&
	   inc.first->id() < inc.second->id()) )
	swap(inc.first, inc.second);
      if( !HPC_helper::duplicateIncoming(inc,incPairs) ) {
	incPairs.push_back(inc);
      }
    }
  }
  size_t nvertices = model()->numberOfVertices(); 
  //To construct resonant diagrams loop over the incoming particles, intermediates
  //and vertices to find allowed diagrams. Need to exclude the diagrams that have 
  //the intermediate as an external particle as well
  for(vector<tcPDPair>::size_type is = 0; is < incPairs.size(); ++is) {
    tPDPair ppi = incPairs[is]; 
    for(vector<PDPtr>::size_type ik = 0; ik < ninter ; ++ik) {
      long ipart = intermediates_[ik]->id();
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
			   intermediates_[ik]);
	}
      }
    }
  }
  //Create matrix element for each process
  const HPDVector::size_type ndiags = diagrams_.size();
  for(HPDVector::size_type ix = 0; ix < ndiags; ++ix)
    createMatrixElement(diagrams_[ix]);
}

void ResonantProcessConstructor::
constructVertex2(IDPair in, VertexBasePtr vertex, 
		 PDPtr partc) {
  //We have the left hand part of the diagram, just need all the possibilities
  //for the RHS
  size_t nvertices = model()->numberOfVertices(); 
  if(processOption_!=3) {
    for(size_t io = 0; io < outgoing_.size(); ++io) {
      tcPDPtr outa = outgoing_[io];
      for(size_t iv = 0; iv < nvertices; ++iv) {
	VBPtr vertex2 = model()->vertex(iv);
	if(vertex2->getNpoint() > 3) continue;
	tPDSet outb = search(vertex2, partc->id(), incoming, outa->id(), outgoing, 
			     outgoing);
	for(tPDSet::const_iterator ita = outb.begin(); ita != outb.end(); ++ita)
	  makeResonantDiagram(in, partc, outa->id(),(**ita).id(), 
			      make_pair(vertex, vertex2));
      }
    }
  }
  else {
    long idRes = !partc->CC() ? partc->id() : partc->CC()->id();
    for(size_t iv = 0; iv < nvertices; ++iv) {
      VBPtr vertex2 = model()->vertex(iv);
      if(vertex2->getNpoint() > 3) continue;
      for(unsigned int ix = 0;ix < 3; ++ix) {
	vector<long> pdlist = vertex2->search(ix, idRes);
	for(unsigned int iy=0;iy<pdlist.size();iy+=3) {
	  long out1 = ix==0 ? pdlist.at(iy+1) : pdlist.at(iy  );
	  long out2 = ix==2 ? pdlist.at(iy+1) : pdlist.at(iy+2);
	  if(partc->mass() < getParticleData(out1)->mass() + 
	     getParticleData(out2)->mass()) continue;
	  makeResonantDiagram(in, partc, out1, out2, 
			      make_pair(vertex, vertex2));
	}
      }
    }
  }
}

void ResonantProcessConstructor::
makeResonantDiagram(IDPair in, PDPtr offshell, long outa, long outb, 
		     VBPair vertpair) {
  assert(vertpair.first && vertpair.second);
  if( abs(outa) == abs(offshell->id()) || 
      abs(outb) == abs(offshell->id())) return;
  HPDiagram newdiag(in,make_pair(outa,outb));
  newdiag.intermediate = offshell;
  newdiag.vertices = vertpair;
  newdiag.channelType = HPDiagram::sChannel;
  fixFSOrder(newdiag);
  assignToCF(newdiag);
  if(duplicate(newdiag,diagrams_)) return;
  // if inclusive enforce the exclusivity
  if(processOption_==1) {
    PDVector::const_iterator loc = 
      std::find(outgoing_.begin(),outgoing_.end(),
		getParticleData(newdiag.outgoing. first));
    if(loc==outgoing_.end()) return;
    loc = 
      std::find(outgoing_.begin(),outgoing_.end(),
		getParticleData(newdiag.outgoing.second));
    if(loc==outgoing_.end()) return;
  }
  else if(processOption_==2) {
    if(!((newdiag.outgoing. first==outgoing_[0]->id()&&
	  newdiag.outgoing.second==outgoing_[1]->id())||
	 (newdiag.outgoing. first==outgoing_[1]->id()&&
	  newdiag.outgoing.second==outgoing_[0]->id())))
      return;
  }
  // add to the list
  diagrams_.push_back(newdiag);
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
				colourFlow(extpart), debug(),scaleChoice_-1,scaleFactor_);
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
