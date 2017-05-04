// -*- C++ -*-
//
// TwoToTwoProcessConstructor.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TwoToTwoProcessConstructor class.
//

#include "TwoToTwoProcessConstructor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include <sstream>

using std::stringstream;

using namespace Herwig;

TwoToTwoProcessConstructor::TwoToTwoProcessConstructor() : 
  Nout_(0), nv_(0), allDiagrams_(true),
  processOption_(0), scaleChoice_(0), scaleFactor_(1.) 
{}

IBPtr TwoToTwoProcessConstructor::clone() const {
  return new_ptr(*this);
}

IBPtr TwoToTwoProcessConstructor::fullclone() const {
  return new_ptr(*this);
}

void TwoToTwoProcessConstructor::doinit() {
  HardProcessConstructor::doinit();
  if(processOption_==2&&outgoing_.size()!=2)
    throw InitException() 
      << "Exclusive processes require exactly"
      << " two outgoing particles but " << outgoing_.size()
      << "have been inserted in TwoToTwoProcessConstructor::doinit()." 
      << Exception::runerror;
  Nout_ = outgoing_.size();
  PDVector::size_type ninc = incoming_.size();
  // exit if nothing to do
  if(Nout_==0||ninc==0) return;
  //create vector of initial-state pairs
  for(PDVector::size_type i = 0; i < ninc; ++i) {
    for(PDVector::size_type j = 0; j < ninc; ++j) {
      tPDPair inc = make_pair(incoming_[i], incoming_[j]);
      
      if( (inc.first->iSpin() > inc.second->iSpin()) ||
	  (inc.first->iSpin() == inc.second->iSpin() &&
	   inc.first->id() < inc.second->id()) )
	swap(inc.first, inc.second);

      if( !HPC_helper::duplicateIncoming(inc,incPairs_) ) {
	incPairs_.push_back(inc);
      }
    }
  }
  // excluded vertices
  excludedVertexSet_ = 
    set<VertexBasePtr>(excludedVertexVector_.begin(),
		       excludedVertexVector_.end());
}


void TwoToTwoProcessConstructor::persistentOutput(PersistentOStream & os) const {
  os << vertices_ << incoming_ << outgoing_
     << allDiagrams_ << processOption_
     << scaleChoice_ << scaleFactor_ << excluded_ << excludedExternal_
     << excludedVertexVector_ << excludedVertexSet_;
}

void TwoToTwoProcessConstructor::persistentInput(PersistentIStream & is, int) {
  is >> vertices_ >> incoming_ >> outgoing_
     >> allDiagrams_ >> processOption_
     >> scaleChoice_ >> scaleFactor_ >> excluded_ >> excludedExternal_
     >> excludedVertexVector_ >> excludedVertexSet_;
}

ClassDescription<TwoToTwoProcessConstructor> 
TwoToTwoProcessConstructor::initTwoToTwoProcessConstructor;
// Definition of the static class description member.

void TwoToTwoProcessConstructor::Init() {

  static ClassDocumentation<TwoToTwoProcessConstructor> documentation
    ("TwoToTwoProcessConstructor constructs the possible diagrams for "
     "a process given the external particles");
 
  static RefVector<TwoToTwoProcessConstructor,ThePEG::ParticleData> interfaceIn
    ("Incoming",
     "Pointers to incoming particles",
     &TwoToTwoProcessConstructor::incoming_, -1, false, false, true, false);

  static RefVector<TwoToTwoProcessConstructor,ThePEG::ParticleData> interfaceOut
    ("Outgoing",
     "Pointers to incoming particles",
     &TwoToTwoProcessConstructor::outgoing_, -1, false, false, true, false);
 
  static Switch<TwoToTwoProcessConstructor,bool> interfaceIncludeAllDiagrams
    ("IncludeEW",
     "Switch to decide which diagrams to include in ME calc.",
     &TwoToTwoProcessConstructor::allDiagrams_, true, false, false);
  static SwitchOption interfaceIncludeAllDiagramsNo
    (interfaceIncludeAllDiagrams,
     "No",
     "Only include QCD diagrams",
     false);
  static SwitchOption interfaceIncludeAllDiagramsYes
   (interfaceIncludeAllDiagrams,
     "Yes",
    "Include EW+QCD.",
    true);

  static Switch<TwoToTwoProcessConstructor,unsigned int> interfaceProcesses
    ("Processes",
     "Whether to generate inclusive or exclusive processes",
     &TwoToTwoProcessConstructor::processOption_, 0, false, false);
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

  static Switch<TwoToTwoProcessConstructor,unsigned int> interfaceScaleChoice
    ("ScaleChoice",
     "&TwoToTwoProcessConstructor::scaleChoice_",
     &TwoToTwoProcessConstructor::scaleChoice_, 0, false, false);
  static SwitchOption interfaceScaleChoiceDefault
    (interfaceScaleChoice,
     "Default",
     "Use if sHat if intermediates all colour neutral, otherwise the transverse mass",
     0);
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

  static Parameter<TwoToTwoProcessConstructor,double> interfaceScaleFactor
    ("ScaleFactor",
     "The prefactor used in the scale calculation. The scale used is"
     " that defined by scaleChoice multiplied by this prefactor",
     &TwoToTwoProcessConstructor::scaleFactor_, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static RefVector<TwoToTwoProcessConstructor,ThePEG::ParticleData> interfaceExcluded
    ("Excluded",
     "Particles which are not allowed as intermediates",
     &TwoToTwoProcessConstructor::excluded_, -1, false, false, true, false, false);

  static RefVector<TwoToTwoProcessConstructor,ParticleData> interfaceExcludedExternal
    ("ExcludedExternal",
     "Particles which are not allowed as outgoing particles",
     &TwoToTwoProcessConstructor::excludedExternal_, -1,
     false, false, true, false, false);

  static RefVector<TwoToTwoProcessConstructor,VertexBase> interfaceExcludedVertices
    ("ExcludedVertices",
     "Vertices which are not included in the 2 -> 2 scatterings",
     &TwoToTwoProcessConstructor::excludedVertexVector_, -1, false, false, true, true, false);

}

void TwoToTwoProcessConstructor::constructDiagrams() {
  if(incPairs_.empty() || outgoing_.empty() || !subProcess() ) return;
  nv_ = model()->numberOfVertices();
  //make sure  vertices are initialised
  for(unsigned int ix = 0; ix < nv_; ++ix ) {
    VertexBasePtr vertex = model()->vertex(ix);
    if(excludedVertexSet_.find(vertex) != 
       excludedVertexSet_.end()) continue;
    vertices_.push_back(vertex);
  }
  nv_ = vertices_.size();
  //Create necessary diagrams
  vector<tcPDPair>::size_type is;
  PDVector::size_type os;
  for(is = 0; is < incPairs_.size(); ++is) {
    tPDPair ppi = incPairs_[is]; 
    for(os = 0; os < Nout_; ++os) { 
      long fs = outgoing_[os]->id();
      for(size_t iv = 0; iv < nv_; ++iv) {
	tVertexBasePtr vertexA = vertices_[iv];

	//This skips an effective vertex and the EW ones if 
	// we only want the strong diagrams
	if( !allDiagrams_ && vertexA->orderInGs() == 0 ) 
	  continue;

	if(vertexA->getNpoint() == 3) {
	  //scattering diagrams
	  createTChannels(ppi, fs, vertexA);
	  
	  //resonance diagrams
	  if( vertexA->isIncoming(ppi.first) &&  
	      vertexA->isIncoming(ppi.second) )
	    createSChannels(ppi, fs, vertexA);
	}
	else 
	  makeFourPointDiagrams(ppi.first->id(), ppi.second->id(),
				fs, vertexA);
      }
    }
  }
  //need to find all of the diagrams that relate to the same process
  //first insert them into a map which uses the '<' operator 
  //to sort the diagrams 
  multiset<HPDiagram> grouped;
  HPDVector::iterator dit = processes_.begin();
  HPDVector::iterator dend = processes_.end();
  bool abort=false;
  for( ; dit != dend; ++dit) {
    // check for on-shell s-channel
    tPDPtr out1 = getParticleData(dit->outgoing.first );
    tPDPtr out2 = getParticleData(dit->outgoing.second);
    if(dit->channelType == HPDiagram::sChannel && 
       dit->intermediate->width()==ZERO &&
       dit->intermediate->mass() > out1->mass()+ out2->mass()) {
      tPDPtr in1 = getParticleData(dit->incoming.first );
      tPDPtr in2 = getParticleData(dit->incoming.second);
      generator()->log() << dit->intermediate->PDGName() 
			 << " can be on-shell in the process "
			 << in1 ->PDGName() << " " <<  in2->PDGName() << " -> "
			 << out1->PDGName() << " " << out2->PDGName() 
			 << " but has zero width.\nEither set the width, enable "
			 << "calculation of its decays, and hence the width,\n"
			 << "or disable it as a potential intermediate using\n"
			 << "insert " << fullName() << ":Excluded 0 "
			 << dit->intermediate->fullName() << "\n---\n";
      abort = true;
    }
    grouped.insert(*dit);
  }
  if(abort) throw Exception() << "One or more processes with zero width"
			      << " resonant intermediates"
			      << Exception::runerror;
  assert( processes_.size() == grouped.size() );
  processes_.clear();
  typedef multiset<HPDiagram>::const_iterator set_iter;
  set_iter it = grouped.begin(), iend = grouped.end();
  while( it != iend ) {
    pair<set_iter,set_iter> range = grouped.equal_range(*it);
    set_iter itb = range.first;
    HPDVector process;
    for( ; itb != range.second; ++itb ) {
      process.push_back(*itb);
    }
    // if inclusive enforce the exclusivity
    if(processOption_==2) {
      if(!((process[0].outgoing. first==outgoing_[0]->id()&&
	    process[0].outgoing.second==outgoing_[1]->id())||
	   (process[0].outgoing. first==outgoing_[1]->id()&&
	    process[0].outgoing.second==outgoing_[0]->id()))) {
	process.clear();
	it = range.second;
	continue;
      }
    }
    if(find(excludedExternal_.begin(),excludedExternal_.end(),
	    getParticleData(process[0].outgoing. first))!=excludedExternal_.end()) {
      process.clear();
      it = range.second;
      continue;
    }
    if(find(excludedExternal_.begin(),excludedExternal_.end(),
	    getParticleData(process[0].outgoing.second))!=excludedExternal_.end()) {
      process.clear();
      it = range.second;
      continue;
    }
    // finally if the process is allow assign the colour flows
    for(unsigned int ix=0;ix<process.size();++ix) assignToCF(process[ix]);
    // create the matrix element
    createMatrixElement(process);
    process.clear();
    it = range.second;
  }
}

void TwoToTwoProcessConstructor::
createSChannels(tcPDPair inpp, long fs, tVertexBasePtr vertex) {
  //Have 2 incoming particle and a vertex, find the possible offshell
  //particles
  pair<long,long> inc = make_pair(inpp.first->id(), inpp.second->id());
  tPDSet offshells = search(vertex, inpp.first->id(), incoming,
			   inpp.second->id(), incoming, outgoing);
  tPDSet::const_iterator it;
  for(it = offshells.begin(); it != offshells.end(); ++it) {
    if(find(excluded_.begin(),excluded_.end(),*it)!=excluded_.end()) continue;
    for(size_t iv = 0; iv < nv_; ++iv) {
      tVertexBasePtr vertexB = vertices_[iv];
      if( vertexB->getNpoint() != 3) continue;
      if( !allDiagrams_ && vertexB->orderInGs() == 0 ) continue;
      
      tPDSet final;
      if( vertexB->isOutgoing(getParticleData(fs)) &&
	  vertexB->isIncoming(*it) )
	final = search(vertexB, (*it)->id(), incoming, fs,
		       outgoing, outgoing);
      //Now make diagrams
      if(!final.empty()) 
	makeDiagrams(inc, fs, final, *it, HPDiagram::sChannel,
		     make_pair(vertex, vertexB), make_pair(true,true));
    }
  }

}

void TwoToTwoProcessConstructor::
createTChannels(tPDPair inpp, long fs, tVertexBasePtr vertex) {
  pair<long,long> inc = make_pair(inpp.first->id(), inpp.second->id());
  //first try a with c
  tPDSet offshells = search(vertex, inpp.first->id(), incoming, fs,
			   outgoing, outgoing);
  tPDSet::const_iterator it;
  for(it = offshells.begin(); it != offshells.end(); ++it) {
    if(find(excluded_.begin(),excluded_.end(),*it)!=excluded_.end()) continue;
     for(size_t iv = 0; iv < nv_; ++iv) {
       tVertexBasePtr vertexB = vertices_[iv];
       if( vertexB->getNpoint() != 3 ) continue;
       if( !allDiagrams_ && vertexB->orderInGs() == 0 ) continue;
       tPDSet final;
       if( vertexB->isIncoming(inpp.second) )
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
    if(find(excluded_.begin(),excluded_.end(),*it)!=excluded_.end()) continue;
    for(size_t iv = 0; iv < nv_; ++iv) {
       tVertexBasePtr vertexB = vertices_[iv];
       if( vertexB->getNpoint() != 3 ) continue;
       if( !allDiagrams_ && vertexB->orderInGs() == 0 ) continue;

       tPDSet final;
       if( vertexB->isIncoming(inpp.first) )
	 final = search(vertexB, inc.first, incoming, (*it)->id(),
			outgoing, outgoing);
       if( !final.empty() )
	 makeDiagrams(inc, fs, final, *it, HPDiagram::tChannel,
		      make_pair(vertexB, vertex), make_pair(true, false));
    }
  }

}

void TwoToTwoProcessConstructor::makeFourPointDiagrams(long parta, long partb,
						   long partc, VBPtr vert) {
  if(processOption_>=1) {
    PDVector::const_iterator loc = find(outgoing_.begin(),outgoing_.end(),
					getParticleData(partc));
    if(loc==outgoing_.end()) return;
  }
  tPDSet ext = search(vert, parta, incoming, partb,incoming, partc, outgoing);
  if( ext.empty() ) return;
  IDPair in(parta, partb);
  for(tPDSet::const_iterator iter=ext.begin(); iter!=ext.end();
      ++iter) {
    if(processOption_>=1) {
      PDVector::const_iterator loc = find(outgoing_.begin(),outgoing_.end(),
					  *iter);
      if(loc==outgoing_.end()) continue;
    }
    HPDiagram nhp(in,make_pair(partc, (*iter)->id()));
    nhp.vertices = make_pair(vert, vert);
    nhp.channelType = HPDiagram::fourPoint;
    fixFSOrder(nhp);
    if( !duplicate(nhp, processes_) ) processes_.push_back(nhp);
  }
}

void 
TwoToTwoProcessConstructor::makeDiagrams(IDPair in, long out1, const tPDSet & out2, 
				     PDPtr inter, HPDiagram::Channel chan, 
				     VBPair vertexpair, BPair cross) {
  if(processOption_>=1) {
    PDVector::const_iterator loc = find(outgoing_.begin(),outgoing_.end(),
					getParticleData(out1));
    if(loc==outgoing_.end()) return;
  }
  for(tPDSet::const_iterator it = out2.begin(); it != out2.end(); ++it) {
    if(processOption_>=1) {
      PDVector::const_iterator loc = find(outgoing_.begin(),outgoing_.end(),
					  *it);
      if(loc==outgoing_.end()) continue;
    }
    HPDiagram nhp( in,make_pair(out1, (*it)->id()) );
    nhp.intermediate = inter;
    nhp.vertices = vertexpair;
    nhp.channelType = chan;
    nhp.ordered = cross;
    fixFSOrder(nhp);
    if( !duplicate(nhp, processes_) ) processes_.push_back(nhp);
  }
}

set<tPDPtr> 
TwoToTwoProcessConstructor::search(VBPtr vertex, long part1, direction d1, 
			       long part2, direction d2, direction d3) {
  if(vertex->getNpoint() != 3) return tPDSet();
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

set<tPDPtr>
TwoToTwoProcessConstructor::search(VBPtr vertex,
				   long part1, direction d1,
				   long part2, direction d2,
				   long part3, direction d3,
				   direction d4) {
  if(vertex->getNpoint() != 4) return tPDSet();
  if(d1 == incoming && getParticleData(part1)->CC()) part1 = -part1;
  if(d2 == incoming && getParticleData(part2)->CC()) part2 = -part2;
  if(d3 == incoming && getParticleData(part3)->CC()) part3 = -part3;
  vector<long> ext;
  tPDSet fourth;
  for(unsigned int ix = 0;ix < 4; ++ix) {
    vector<long> pdlist = vertex->search(ix, part1);
    ext.insert(ext.end(), pdlist.begin(), pdlist.end());
  }
  for(unsigned int ix = 0;ix < ext.size(); ix += 4) {
    long id0 = ext.at(ix); long id1 = ext.at(ix + 1);
    long id2 = ext.at(ix + 2); long id3 = ext.at(ix + 3);
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
      tPDPtr p = getParticleData(ext[pos]);
      if(d4 == incoming && p->CC()) 
	p = p->CC();
      fourth.insert(p);
    }
  } 
  return fourth;
}

void 
TwoToTwoProcessConstructor::createMatrixElement(const HPDVector & process) const {
  if ( process.empty() ) return;
  // external particles
  tcPDVector extpart(4);
  extpart[0] = getParticleData(process[0].incoming.first);
  extpart[1] = getParticleData(process[0].incoming.second);
  extpart[2] = getParticleData(process[0].outgoing.first);
  extpart[3] = getParticleData(process[0].outgoing.second);
  // create the object
  string objectname ("/Herwig/MatrixElements/");
  string classname = MEClassname(extpart, objectname);
  GeneralHardMEPtr matrixElement = dynamic_ptr_cast<GeneralHardMEPtr>
      (generator()->preinitCreate(classname, objectname));
  if( !matrixElement ) {
    std::stringstream message;
    message << "TwoToTwoProcessConstructor::createMatrixElement "
	    << "- No matrix element object could be created for "
	    << "the process " 
	    << extpart[0]->PDGName() << " " << extpart[0]->iSpin() << "," 
	    << extpart[1]->PDGName() << " " << extpart[1]->iSpin() << "->" 
	    << extpart[2]->PDGName() << " " << extpart[2]->iSpin() << "," 
	    << extpart[3]->PDGName() << " " << extpart[3]->iSpin() 
	    << ".  Constructed class name: \"" << classname << "\"";
    generator()->logWarning(TwoToTwoProcessConstructorError(message.str(),Exception::warning));
    return;
  }
  // choice for the scale
  unsigned int scale;
  if(scaleChoice_==0) {
    // check coloured initial and final state
    bool inColour  = ( extpart[0]->coloured() ||
		       extpart[1]->coloured());
    bool outColour = ( extpart[2]->coloured() ||
		       extpart[3]->coloured());
    if(inColour&&outColour) {
      bool coloured = false;
      for(unsigned int ix=0;ix<process.size();++ix) {
	if(process[ix].intermediate&&
	   process[ix].intermediate->coloured()) {
	  coloured = true;
	  break;
	}
      }
      scale = coloured ? 1 : 0;
    }
    else {
      scale = 0;
    } 
  }
  else {
    scale = scaleChoice_-1;
  }
  // set the information
  matrixElement->setProcessInfo(process, colourFlow(extpart),
				debug(), scale, scaleFactor_);
  // insert it
  generator()->preinitInterface(subProcess(), "MatrixElements", 
				subProcess()->MEs().size(),
				"insert", matrixElement->fullName()); 
}

string TwoToTwoProcessConstructor::MEClassname(const vector<tcPDPtr> & extpart, 
					   string & objname) const {
  string classname("Herwig::ME");
  for(vector<tcPDPtr>::size_type ix = 0; ix < extpart.size(); ++ix) {
    if(ix == 2) classname += "2";
    if(extpart[ix]->iSpin() == PDT::Spin0) classname += "s";
    else if(extpart[ix]->iSpin() == PDT::Spin1) classname += "v";
    else if(extpart[ix]->iSpin() == PDT::Spin1Half) classname += "f";
    else if(extpart[ix]->iSpin() == PDT::Spin2) classname += "t";
    else {
      std::stringstream message;
      message << "MEClassname() : Encountered an unknown spin for "
	      << extpart[ix]->PDGName() << " while constructing MatrixElement "
	      << "classname " << extpart[ix]->iSpin();
      generator()->logWarning(TwoToTwoProcessConstructorError(message.str(),Exception::warning));
    }
  }
  objname += "ME" + extpart[0]->PDGName() + extpart[1]->PDGName() + "2" 
    + extpart[2]->PDGName() + extpart[3]->PDGName();
  return classname;  
}
