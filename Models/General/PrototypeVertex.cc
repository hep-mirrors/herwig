// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FourBodyDecayConstructor class.
//

#include "ThePEG/Helicity/Vertex/VertexBase.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Utilities/Exception.h"
#include "PrototypeVertex.h"
#include "NBodyDecayConstructorBase.h"

using namespace Herwig;

void PrototypeVertex::createPrototypes(tPDPtr inpart, VertexBasePtr vertex,
				       std::stack<PrototypeVertexPtr> & prototypes,
				       NBodyDecayConstructorBasePtr decayCon) {
  if(!vertex->isIncoming(inpart)) return;
  tPDPtr ccpart = inpart->CC() ? inpart->CC() : inpart;
  long id = ccpart->id();
  Energy2 scale = sqr(inpart->mass());
  for(unsigned int list=0;list<vertex->getNpoint();++list) {
    tPDVector decaylist = vertex->search(list, ccpart);
    tPDVector::size_type nd = decaylist.size();
    for( tPDVector::size_type i = 0; i < nd; i += vertex->getNpoint() ) {
      tPDVector pout(decaylist.begin()+i,
		     decaylist.begin()+i+vertex->getNpoint());
      OrderedVertices out;
      for(unsigned int ix=1;ix<pout.size();++ix) {
	if(pout[ix]->id() == id ) swap(pout[0], pout[ix]);
	out.insert(make_pair(pout[ix],PrototypeVertexPtr()));
      }
      // remove flavour changing vertices
      if(decayCon->removeFlavourChangingVertices()&&vertex->getNpoint()==3) {
	if((pout[1]->id()==ParticleID::gamma||
	    pout[1]->id()==ParticleID::g) &&
	   inpart!=pout[2]) continue;
	else if((pout[2]->id()==ParticleID::gamma||
		 pout[2]->id()==ParticleID::g) &&
		inpart!=pout[1]) continue;
      }
      // remove small vertices if needed
      if(decayCon->removeSmallVertices()) {
	if(vertex->getNpoint()==3)
	  vertex->setCoupling(scale,pout[0],pout[1],pout[2]);
	else if(vertex->getNpoint()==4)
	  vertex->setCoupling(scale,pout[0],pout[1],pout[2],pout[3]);
	if(abs(vertex->norm())<decayCon->minVertexNorm()) continue;
      }
      if(vertex->getNpoint()==3) {
	// remove radiation
	if((inpart==pout[1] && (pout[2]->id()==ParticleID::gamma||
				pout[2]->id()==ParticleID::g||
				pout[2]->id()==ParticleID::Z0)) 
	   ||
	   (inpart==pout[2] && (pout[1]->id()==ParticleID::gamma||
				pout[1]->id()==ParticleID::g||
				pout[1]->id()==ParticleID::Z0)))
	  continue;
      }
      prototypes.push(new_ptr(PrototypeVertex(inpart,out,vertex,
					      int(vertex->getNpoint())-1)));
    }
  }
}

PrototypeVertexPtr PrototypeVertex::replicateTree(PrototypeVertexPtr parent,
						  PrototypeVertexPtr oldChild,
						  PrototypeVertexPtr & newChild) {
  PrototypeVertexPtr newParent = 
    new_ptr(PrototypeVertex(parent->incoming,OrderedVertices(),
			    parent->vertex,parent->npart));
  for(OrderedVertices::const_iterator it = parent->outgoing.begin();
      it!=parent->outgoing.end();++it) {
    PrototypeVertexPtr child = it->second ? 
      replicateTree(it->second,oldChild,newChild) :
      PrototypeVertexPtr();
    newParent->outgoing.insert(make_pair(it->first,child));
    if(child) child->parent = newParent;
    if(it->second==oldChild) newChild=child;
  }
  if(oldChild==parent) newChild=newParent;
  return newParent;
}

void PrototypeVertex::expandPrototypes(PrototypeVertexPtr proto, VertexBasePtr vertex,
				       std::stack<PrototypeVertexPtr> & prototypes,
				       const set<PDPtr> & excluded,
				       NBodyDecayConstructorBasePtr decayCon) {
  for(OrderedVertices::const_iterator it = proto->outgoing.begin();
      it!=proto->outgoing.end();++it) {
    if(it->second) {
      expandPrototypes(it->second,vertex,prototypes,excluded,decayCon);
    }
    else {
      if(!vertex->isIncoming(it->first)) continue;
      if(excluded.find(it->first)!=excluded.end()) continue;
      if(it->first->CC() && 
	 excluded.find(it->first->CC())!=excluded.end()) continue;

      tcPDPtr ccpart = it->first->CC() ? it->first->CC() : it->first;
      long id = ccpart->id();

      PrototypeVertexPtr parent=proto;
      while(parent->parent) parent=parent->parent;
      Energy2 scale(sqr(parent->incoming->mass()));
      for(unsigned int il = 0; il < vertex->getNpoint(); ++il) {
	tPDVector decaylist = vertex->search(il,ccpart);
	tPDVector::size_type nd = decaylist.size();
	for( tPDVector::size_type i = 0; i < nd; i += vertex->getNpoint() ) {
	  tPDVector pout(decaylist.begin()+i,
			 decaylist.begin()+i+vertex->getNpoint());
	  OrderedVertices outgoing;
	  for(unsigned int iy=1;iy<pout.size();++iy) {
	    if(pout[iy]->id() == id ) swap(pout[0], pout[iy]);
	    outgoing.insert(make_pair(pout[iy],PrototypeVertexPtr()));
	  }
	  if(vertex->getNpoint()==3) {
	    if((it->first==pout[1] && (pout[2]->id()==ParticleID::gamma||
				       pout[2]->id()==ParticleID::g||
				       pout[2]->id()==ParticleID::Z0)) 
	       ||
	       (it->first==pout[2] && (pout[1]->id()==ParticleID::gamma||
				       pout[1]->id()==ParticleID::g||
				       pout[1]->id()==ParticleID::Z0)))
	      continue;
	    // remove weak decays of quarks other than top
	    if(StandardQCDPartonMatcher::Check(pout[0]->id()) &&
	       ((StandardQCDPartonMatcher::Check(pout[1]->id()) &&
		 abs(pout[2]->id())==ParticleID::Wplus)||
		(StandardQCDPartonMatcher::Check(pout[2]->id())&&
		 abs(pout[1]->id())==ParticleID::Wplus))) continue;
	    // remove weak decays of leptons
	    if((abs(pout[0]->id())>=11&&abs(pout[0]->id())<=16) &&
		(((abs(pout[1]->id())>=11&&abs(pout[1]->id())<=16) &&
		  abs(pout[2]->id())==ParticleID::Wplus)||
		((abs(pout[2]->id())>=11&&abs(pout[2]->id())<=16)&&
		 abs(pout[1]->id())==ParticleID::Wplus))) continue;
	    }
	  // remove flavour changing vertices
	  if(decayCon->removeFlavourChangingVertices()&&vertex->getNpoint()==3) {
	    if((pout[1]->id()==ParticleID::gamma||
		pout[1]->id()==ParticleID::g) &&
	       it->first!=pout[2]) continue;
	    else if((pout[2]->id()==ParticleID::gamma||
		     pout[2]->id()==ParticleID::g) &&
		    it->first!=pout[1]) continue;
	  }
	  // remove small vertices if needed
	  if(decayCon->removeSmallVertices()) {
	    if(vertex->getNpoint()==3)
	      vertex->setCoupling(scale,pout[0],pout[1],pout[2]);
	    else if(vertex->getNpoint()==4)
	      vertex->setCoupling(scale,pout[0],pout[1],pout[2],pout[3]);
	    if(abs(vertex->norm())<decayCon->minVertexNorm()) continue;
	  }
	  PrototypeVertexPtr newBranch = 
	    new_ptr(PrototypeVertex(it->first,
				    outgoing,vertex,int(vertex->getNpoint())-1));
	  PrototypeVertexPtr newChild;
	  PrototypeVertexPtr newVertex = replicateTree(parent,proto,newChild);
	  newBranch->parent = newChild;
	  OrderedVertices::iterator kt = newChild->outgoing.begin();
	  for(OrderedVertices::const_iterator jt = proto->outgoing.begin();
	      jt!=it;++jt,++kt) {;}
	  pair< tPDPtr, PrototypeVertexPtr > newPair = make_pair(kt->first,newBranch);
	  newChild->outgoing.erase(kt);
	  newChild->outgoing.insert(newPair);
	  newChild->incrementN(newBranch->npart-1);
	  prototypes.push(newVertex);
	}
      }
    }
  }
}

bool PrototypeVertex::canBeOnShell(unsigned int opt,Energy maxMass,bool first) {
  Energy outMass = outgoingMass();
  if(!first) { 
    bool in  = maxMass>incomingMass();
    bool out = incomingMass()>outMass;
    if(opt!=0) {
      if( in && ( out || opt==2 ) ) return true;
    }
    else if (incoming->width() == ZERO) {
      tPrototypeVertexPtr original = this;
      while(original->parent) {
	original=original->parent;
      };
      ostringstream name;
      name << original->incoming->PDGName() << " -> ";
      for(OrderedParticles::const_iterator it = original->outPart.begin();
	  it!= original->outPart.end();++it)
	name << (**it).PDGName() << " ";
      Throw<InitException>() 
	<< "Trying to include on-shell diagram in decay"
	<< name.str() << "including on-shell "
	<< incoming->PDGName() << " with zero width.\n"
	<< "You should make sure that the width for the intermediate is either"
	<< " read from an SLHA file or the intermediate is included in the "
	<< "DecayParticles list of the ModelGenerator.\n"
	<< Exception::runerror;
    }
  }
  else maxMass = incomingMass();
  // check the decay products
  for(OrderedVertices::const_iterator it = outgoing.begin();
      it!=outgoing.end();++it) {
    if(!it->second) continue;
    Energy newMax = maxMass-outMass+it->second->outgoingMass();
    if(it->second->canBeOnShell(opt,newMax,false)) {
      if(first) possibleOnShell = true;
      return true;
    }
  }
  return false;
}

bool PrototypeVertex::sameDecay(const PrototypeVertex & x) const {
  if(incoming != x.incoming) return false;
  if(outPart.empty())     setOutgoing();
  if(x.outPart.empty()) x.setOutgoing();
  OrderedParticles::const_iterator cit =   outPart.begin();
  OrderedParticles::const_iterator cjt = x.outPart.begin();
  if(x.npart!=npart) return false;
  for(;cit!=outPart.end();++cit,++cjt) {
    if(*cit!=*cjt) return false;
  }
  return true;
}

namespace {
  void constructTag(vector<unsigned int> & order,NBVertex & vertex,
		    const OrderedParticles & outgoing,
		    vector<bool> & matched) {
    for(list<pair<PDPtr,NBVertex> >::iterator it=vertex.vertices.begin();
	it!=vertex.vertices.end();++it) {
      // search down the tree
      if(it->second.vertex) {
	constructTag(order,it->second,outgoing,matched);
      }
      // identify this particle
      else {
	unsigned int ix=0;
	for(OrderedParticles::const_iterator pit=outgoing.begin();
	    pit!=outgoing.end();++pit) {
	  if(*pit==it->first&&!matched[ix]) {
	    matched[ix] = true;
	    order.push_back(ix+1);
	    break;
	  }
	  ++ix;
	}
      }
    }
  }
}

NBDiagram::NBDiagram(PrototypeVertexPtr proto) 
  : NBVertex(proto),
    colourFlow       (vector<CFPair>(1,make_pair(1,1.))),
    largeNcColourFlow(vector<CFPair>(1,make_pair(1,1.))) {
  if(!proto) return;
  // finally work out the channel and the ordering
  vector<bool> matched(outgoing.size(),false);
  for(list<pair<PDPtr,NBVertex> >::iterator it=vertices.begin();
      it!=vertices.end();++it) {
    // search down the tree
    if(it->second.vertex) {
      constructTag(channelType,it->second,outgoing,matched);
    }
    // identify this particle
    else {
      unsigned int ix=0;
      for(OrderedParticles::const_iterator pit=outgoing.begin();
	  pit!=outgoing.end();++pit) {
	if(*pit==it->first&&!matched[ix]) {
	  matched[ix] = true;
	  channelType.push_back(ix+1);
	  break;
	}
	++ix;
      }
    }
  }
}

NBVertex::NBVertex(PrototypeVertexPtr proto) { 
  if(!proto) return;
  incoming = proto->incoming;
  outgoing = proto->outPart;
  vertex   = proto->vertex;
  // create the vertices
  for(OrderedVertices::const_iterator it=proto->outgoing.begin();
      it!=proto->outgoing.end();++it) {
    vertices.push_back(make_pair(it->first,NBVertex(it->second)));
  }
  // now let's re-order so that branchings are at the end
  for(list<pair<PDPtr,NBVertex> >::iterator it=vertices.begin();
      it!=vertices.end();++it) {
    if(!it->second.incoming) continue;
    list<pair<PDPtr,NBVertex> >::iterator jt=it;
    for( ; jt!=vertices.end();++jt) {
      if(jt==it) continue;
      if(!jt->second.incoming) {
	break;
      }
    }
    if(jt!=vertices.end()) {
      list<pair<PDPtr,NBVertex> >::iterator kt = it;
      while(kt!=jt) {
	list<pair<PDPtr,NBVertex> >::iterator lt = kt;
	++lt;
	swap(*kt,*lt);
	kt=lt;
      }
    }
  }
}
