// -*- C++ -*-
//
// PrototypeVertex.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_PrototypeVertex_H
#define HERWIG_PrototypeVertex_H
#include <stack>
#include "ThePEG/Helicity/Vertex/VertexBase.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/EnumIO.h"
#include "NBodyDecayConstructorBase.fh"
//
// This is the declaration of the PrototypeVertex class.
//

namespace Herwig {
using namespace ThePEG;
using Helicity::VertexBasePtr;

class PrototypeVertex;
ThePEG_DECLARE_POINTERS(Herwig::PrototypeVertex,PrototypeVertexPtr);
  
/** Pair of int,double */
typedef pair<unsigned int, double> CFPair;


/**
 *  A struct to order the particles in the same way as in the DecayModes
 */
struct ParticleOrdering {
  /**
   *  Operator for the ordering
   * @param p1 The first ParticleData object
   * @param p2 The second ParticleData object
   */
  bool operator() (tcPDPtr p1, tcPDPtr p2) const {
    return abs(p1->id()) > abs(p2->id()) ||
      ( abs(p1->id()) == abs(p2->id()) && p1->id() > p2->id() ) ||
      ( p1->id() == p2->id() && p1->fullName() > p2->fullName() );
  }
};

/**
 *  A struct to order the particles in the same way as in the DecayMode's
 */
struct VertexOrdering {
  /**
   *  Operator for the ordering
   * @param p1 The first ParticleData object
   * @param p2 The second ParticleData object
   */
  bool operator() (const pair< tPDPtr, PrototypeVertexPtr > & p1,
		   const pair< tPDPtr, PrototypeVertexPtr > & p2) const  {
    return  abs(p1.first->id()) > abs(p2.first->id()) ||
      ( abs(p1.first->id()) == abs(p2.first->id()) && p1.first->id() > p2.first->id() ) ||
      ( p1.first->id() == p2.first->id() && p1.first->fullName() > p2.first->fullName() );
  }
};
  
typedef multiset<pair< tPDPtr, PrototypeVertexPtr >,VertexOrdering > OrderedVertices;

/**
 * A set of ParticleData objects ordered as for the DecayMode's
 */
typedef multiset<PDPtr,ParticleOrdering> OrderedParticles;

/**
 *  Storage of a potenital n-body decay
 */
class PrototypeVertex : public Base {
  
public:

  /**
   *  Default Constructor
   */
  PrototypeVertex() : npart(0), possibleOnShell(false) {}

  /**
   *  Constructor
   */
  PrototypeVertex(tPDPtr in, OrderedVertices out,
		  VertexBasePtr v, int n) :
    incoming(in), outgoing(out), vertex(v), npart(n),
    possibleOnShell(false) {}

  /**
   *  Incoming particle
   */
  tPDPtr incoming;

  /**
   *  Outgoing particles
   */
  OrderedVertices outgoing;

  /**
   *  The vertex for the interaction
   */
  VertexBasePtr vertex;

  /**
   *  The parent of the vertex
   */
  tPrototypeVertexPtr parent;
  
  /**
   *  Number of particles
   */
  unsigned int npart;

  /**
   *  Outgoing particles
   */
  mutable OrderedParticles outPart;

  /**
   *  Can have on-shell intermediates
   */
  bool possibleOnShell;

  /**
   *  Increment the number of particles
   */
  void incrementN(int in) {
    npart += in;
    if(parent) parent->incrementN(in);
  }

  /**
   *  Mass of the incoming particle
   */
  Energy incomingMass() {
    return incoming->mass();
  }

  /**
   *  Total mass of all the outgoing particles
   */
  Energy outgoingMass() {
    Energy mass(ZERO);
    for(OrderedVertices::const_iterator it = outgoing.begin();
	it!=outgoing.end();++it) {
      mass += it->second ? 
	it->second->outgoingMass() : it->first->mass();
    }
    return mass;
  }

  /**
   * Total constituent mass of all the outgoing particles
   */
  Energy outgoingConstituentMass() {
    Energy mass(ZERO);
    for(OrderedVertices::const_iterator it = outgoing.begin();
	it!=outgoing.end();++it) {
      mass += it->second ? 
	it->second->outgoingConstituentMass() : it->first->constituentMass();
    }
    return mass;
  }

  /**
   * Check the external particles
   */
  bool checkExternal(bool first=true) {
    if(outPart.empty())     setOutgoing();
    if(first&&outPart.find(incoming)!=outPart.end()) return false;
    bool output = true;
    for(OrderedVertices::const_iterator it = outgoing.begin();
	it!=outgoing.end();++it) {
      if(it->second&& !it->second->checkExternal(false)) output = false;
    }
    return output;
  }
  
  /**
   * Set the outgoing particles
   */
  void setOutgoing() const {
    assert(outPart.empty());
    for(OrderedVertices::const_iterator it = outgoing.begin();
	it!=outgoing.end();++it) {
      if(it->second) {
	it->second->setOutgoing();
	outPart.insert(it->second->outPart.begin(),
		       it->second->outPart.end());
      }
      else
	outPart.insert(it->first);
    }
  }

  /**
   *  Are there potential on-shell intermediates?
   */
  bool canBeOnShell(unsigned int opt,Energy maxMass,bool first);

  /**
   *  Check if same external particles
   */
  bool sameDecay(const PrototypeVertex & x) const;

  /**
   *  Create a \f$1\to2\f$ prototype
   */
  static  void createPrototypes(tPDPtr inpart, VertexBasePtr vertex,
				std::stack<PrototypeVertexPtr> & prototypes,
				NBodyDecayConstructorBasePtr decayCon);

  /**
   *  Expand the prototypes by adding more legs
   */
  static void expandPrototypes(PrototypeVertexPtr proto, VertexBasePtr vertex,
			       std::stack<PrototypeVertexPtr> & prototypes,
			       const set<PDPtr> & excluded,
			       NBodyDecayConstructorBasePtr decayCon);

  /**
   *  Copy the whole structure with a new branching
   */
  static PrototypeVertexPtr replicateTree(PrototypeVertexPtr parent,
					  PrototypeVertexPtr oldChild,
					  PrototypeVertexPtr & newChild);

};

/**
 * Output to a stream 
 */
inline ostream & operator<<(ostream & os, const PrototypeVertex & diag) {
  os << diag.incoming->PDGName() << " -> ";
  bool seq=false;
  for(OrderedVertices::const_iterator it = diag.outgoing.begin();
      it!=diag.outgoing.end();++it) {
    os << it->first->PDGName() << " ";
    if(it->second) seq = true;
  }
  os << " decays via "
     << diag.vertex->fullName() << " in a " 
     << diag.npart << "-body decay\n";
  if(!seq) return os;
  os << "Followed by\n";
  for(OrderedVertices::const_iterator it = diag.outgoing.begin();
      it!=diag.outgoing.end();++it) {
    if(it->second) os << *it->second;
  }
  return os;
}

/**
 * Test whether two diagrams are identical.
 */
inline bool operator==(const PrototypeVertex & x, const PrototypeVertex & y) {
  if(x.incoming != y.incoming) return false;
  if(x.vertex != y.vertex) return false;
  if(x.npart != y.npart) return false;
  if(x.outgoing.empty()&&y.outgoing.empty()) return true;
  if(x.outgoing.size() != y.outgoing.size()) return false;
  OrderedVertices::const_iterator xt = x.outgoing.begin();
  OrderedVertices::const_iterator yt = y.outgoing.begin();
  while(xt!=x.outgoing.end()) {
    if(xt->first != yt->first) return false;
    // special for identical particles
    OrderedVertices::const_iterator lxt = x.outgoing.lower_bound(*xt);
    OrderedVertices::const_iterator uxt = x.outgoing.upper_bound(*xt);
    --uxt;
    // just one particle
    if(lxt==uxt) {
      if(xt->second && yt->second) {
	if(*(xt->second)==*(yt->second)) {
	  ++xt;
	  ++yt;
	  continue;
	}
	else return false;
      }
      else if(xt->second || yt->second)
	return false;
      ++xt;
      ++yt;
    }
    // identical particles
    else {
      ++uxt;
      OrderedVertices::const_iterator lyt = y.outgoing.lower_bound(*xt);
      OrderedVertices::const_iterator uyt = y.outgoing.upper_bound(*xt);
      unsigned int nx=0;
      for(OrderedVertices::const_iterator ixt=lxt;ixt!=uxt;++ixt) {++nx;}
      unsigned int ny=0;
      for(OrderedVertices::const_iterator iyt=lyt;iyt!=uyt;++iyt) {++ny;}
      if(nx!=ny) return false;
      vector<bool> matched(ny,false);
      for(OrderedVertices::const_iterator ixt=lxt;ixt!=uxt;++ixt) {
	bool found = false;
	unsigned int iy=0;
	for(OrderedVertices::const_iterator iyt=lyt;iyt!=uyt;++iyt) {
	  if(matched[iy]) {
	    ++iy;
	    continue;
	  }
	  if( (!ixt->second &&!iyt->second)  ||
	      ( ixt->second&&iyt->second &&
		*(ixt->second)==*(iyt->second)) ) {
	    matched[iy] = true;
	    found = true;
	    break;
	  }
	  ++iy;
	}
	if(!found) return false;
      }
      xt=uxt;
      yt=uyt;
    }
  }
  return true;
}

/**
 *  A simple vertex for the N-body diagram
 */
struct NBVertex {

  /**
   * Constructor taking a prototype vertex as the arguments
   */
  NBVertex(PrototypeVertexPtr proto = PrototypeVertexPtr() );

  /**
   * Incoming particle
   */
  tPDPtr incoming;

  /**
   *  Outgoing particles
   */
  mutable OrderedParticles outgoing;

  /**
   *  The vertices
   */
  list<pair<PDPtr,NBVertex> > vertices;

  /**
   *  The vertex
   */
  VertexBasePtr vertex;
};

/**
 * The NBDiagram struct contains information about a \f$1\to n\f$ decay
 * that has been automatically generated.
 */
  struct NBDiagram : public NBVertex {

  /**
   * Constructor taking a prototype vertex as the arguments*/
  NBDiagram(PrototypeVertexPtr proto=PrototypeVertexPtr());

  /**
   *  The type of channel
   */
  vector<unsigned int> channelType;

  /** Store colour flow at \f$N_c=3\f$ information */
  mutable vector<CFPair> colourFlow;
  
  /** Store colour flow at \f$N_c=\infty\f$ information */
  mutable vector<CFPair> largeNcColourFlow;
};

/** 
 * Output operator to allow the structure to be persistently written
 * @param os The output stream
 * @param x The NBVertex 
 */
inline PersistentOStream & operator<<(PersistentOStream & os, 
				      const NBVertex  & x) {
  os << x.incoming << x.outgoing << x.vertices << x.vertex;
  return os;
}
  
/** 
 * Input operator to allow persistently written data to be read in
 * @param is The input stream
 * @param x The NBVertex 
 */
inline PersistentIStream & operator>>(PersistentIStream & is,
				      NBVertex & x) {
  is >> x.incoming >> x.outgoing >> x.vertices >> x.vertex;
  return is;
}

/** 
 * Output operator to allow the structure to be persistently written
 * @param os The output stream
 * @param x The NBDiagram 
 */
inline PersistentOStream & operator<<(PersistentOStream & os, 
				      const NBDiagram  & x) {
  os << x.incoming << x.channelType << x.outgoing << x.vertices << x.vertex
     << x.colourFlow << x.largeNcColourFlow;
  return os;
}
  
/** 
 * Input operator to allow persistently written data to be read in
 * @param is The input stream
 * @param x The NBDiagram 
 */
inline PersistentIStream & operator>>(PersistentIStream & is,
				      NBDiagram & x) {
  is >> x.incoming >> x.channelType >> x.outgoing >> x.vertices >> x.vertex
     >> x.colourFlow >> x.largeNcColourFlow;
  return is;
}

/**
 * Output a NBVertex to a stream 
 */
inline ostream & operator<<(ostream & os, const NBVertex & vertex) {
  os << vertex.incoming->PDGName() << " -> ";
  bool seq=false;
  for(list<pair<PDPtr,NBVertex> >::const_iterator it=vertex.vertices.begin();
      it!=vertex.vertices.end();++it) {
    os << it->first->PDGName() << " ";
    if(it->second.incoming) seq = true;
  }
  os << "via vertex " << vertex.vertex->fullName() << "\n";
  if(!seq) return os;
  os << "Followed by\n";
  for(list<pair<PDPtr,NBVertex> >::const_iterator it=vertex.vertices.begin();
      it!=vertex.vertices.end();++it) {
    if(it->second.incoming) os << it->second;
  }
  return os;
}

/**
 * Output a NBDiagram to a stream 
 */
inline ostream & operator<<(ostream & os, const NBDiagram & diag) {
  os << diag.incoming->PDGName() << " -> ";
  for(OrderedParticles::const_iterator it=diag.outgoing.begin();
      it!=diag.outgoing.end();++it) {
    os << (**it).PDGName() << " ";
  }
  os << " has order ";
  for(unsigned int ix=0;ix<diag.channelType.size();++ix)
    os << diag.channelType[ix] << " ";
  os << "\n";
  os << "First decay " << diag.incoming->PDGName() << " -> ";
  bool seq=false;
  for(list<pair<PDPtr,NBVertex> >::const_iterator it=diag.vertices.begin();
      it!=diag.vertices.end();++it) {
    os << it->first->PDGName() << " ";
    if(it->second.incoming) seq = true;
  }
  os << "via vertex " << diag.vertex->fullName() << "\n";
  if(!seq) return os;
  os << "Followed by\n";
  for(list<pair<PDPtr,NBVertex> >::const_iterator it=diag.vertices.begin();
      it!=diag.vertices.end();++it) {
    if(it->second.incoming) os << it->second;
  }
  return os;
}

}

#endif /* HERWIG_PrototypeVertex_H */
