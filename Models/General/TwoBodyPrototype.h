// -*- C++ -*-
//
// TwoBodyPrototype.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_TwoBodyPrototype_H
#define HERWIG_TwoBodyPrototype_H
#include <queue>
#include "ThePEG/Helicity/Vertex/VertexBase.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/EnumIO.h"
//
// This is the declaration of the TwoBodyPrototype struct.
//

namespace Herwig {
using namespace ThePEG;
using Helicity::VertexBasePtr;

class PrototypeVertex;
ThePEG_DECLARE_POINTERS(Herwig::PrototypeVertex,PrototypeVertexPtr);
  
/** Pair of int,double */
typedef pair<unsigned int, double> CFPair;


/**
 *  A struct to order the particles in the same way as in the DecayMode's
 */
struct ParticleOrdering {
  /**
   *  Operator for the ordering
   * @param p1 The first ParticleData object
   * @param p2 The second ParticleData object
   */
  bool operator()(PDPtr p1, PDPtr p2) {
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
  bool operator()(pair< tPDPtr, PrototypeVertexPtr > p1,
		  pair< tPDPtr, PrototypeVertexPtr > p2) {
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
 * A two body decay mode which is a prototype for the 
 * three body mode
 */
struct TwoBodyPrototype {

  /**
   *  Constructor
   */
  TwoBodyPrototype(tPDPtr in, tPDPair out, VertexBasePtr v) :
    incoming(in), outgoing(out), vertex(v) {}

  /**
   *  Incoming particle
   */
  tPDPtr incoming;

  /**
   *  Outgoing particles
   */
  tPDPair outgoing;

  /**
   *  The vertex for the interaction
   */
  VertexBasePtr vertex;

  static  vector<TwoBodyPrototype> 
  createPrototypes(tPDPtr inpart, VertexBasePtr vertex, unsigned int list,
		   Energy weakCut) {
    int id = inpart->id();
    if( id < 0 || !vertex->isIncoming(inpart) || vertex->getNpoint() != 3 )
      return vector<TwoBodyPrototype>();
    tPDVector decaylist = vertex->search(list, inpart);
    vector<TwoBodyPrototype> decays;
    tPDVector::size_type nd = decaylist.size();
    for( tPDVector::size_type i = 0; i < nd; i += 3 ) {
      tPDPtr pa(decaylist[i]), pb(decaylist[i + 1]), pc(decaylist[i + 2]);
      if( pb->id() == id ) swap(pa, pb);
      if( pc->id() == id ) swap(pa, pc);
      //vertices are defined with all particles incoming
      if( pb->CC() ) pb = pb->CC();
      if( pc->CC() ) pc = pc->CC();
      // remove weak processes simulated using the current
      if(weakCut>ZERO) {
	if(abs(pb->id())==ParticleID::Wplus && pc->mass() < pa->mass() &&
	   pa->mass()-pc->mass()<weakCut) continue;
	if(abs(pc->id())==ParticleID::Wplus && pb->mass() < pa->mass() &&
	   pa->mass()-pb->mass()<weakCut) continue;
      }
      decays.push_back(TwoBodyPrototype(inpart,make_pair(pb,pc),vertex));
    }
    return decays;
  }

};

class PrototypeVertex : public Base {
  
public:

  /**
   *  Default Constructor
   */
  PrototypeVertex() {}

  /**
   *  Constructor
   */
  PrototypeVertex(tPDPtr in, OrderedVertices out,
		  VertexBasePtr v, int n) :
  incoming(in), outgoing(out), vertex(v), npart(n) {}

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
  int npart;

  /**
   *  Outgoing particles
   */
  mutable OrderedParticles outPart;

  void incrementN(int in) {
    npart += in;
    if(parent) parent->incrementN(in);
  }

  Energy incomingMass() {
    return incoming->mass();
  }

  Energy outgoingMass() {
    Energy mass(ZERO);
    for(OrderedVertices::const_iterator it = outgoing.begin();
	it!=outgoing.end();++it) {
      mass += it->second ? 
	it->second->outgoingMass() : it->first->mass();
    }
    return mass;
  }

  bool checkExternal() {
    if(outPart.empty())     setOutgoing();
    if(outPart.find(incoming)!=outPart.end()) return false;
    bool output = true;
    for(OrderedVertices::const_iterator it = outgoing.begin();
	it!=outgoing.end();++it) {
      if(it->second&& !it->second->checkExternal()) output = false;
    }
    return output;
  }
  
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


  bool sameDecay(const PrototypeVertex & x) const {
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

  static  void createPrototypes(tPDPtr inpart, VertexBasePtr vertex,
				std::queue<PrototypeVertexPtr> & prototypes,
				Energy weakCut);

  static void expandPrototypes(PrototypeVertexPtr proto, VertexBasePtr vertex,
			       std::queue<PrototypeVertexPtr> & prototypes);

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
  os << "\n Followed by\n";
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
  if(x.outgoing.size() != y.outgoing.size()) return false;
  OrderedVertices::const_iterator xt = x.outgoing.begin();
  OrderedVertices::const_iterator yt = y.outgoing.begin();
  for(;xt!=x.outgoing.end();++xt,++yt) {
    if(xt->first != yt->first) return false;
    if(xt->second && yt->second) {
      if(*(xt->second)==*(yt->second)) continue;
      else return false;
    }
    else if(xt->second || yt->second)
      return false;
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
 * The NBDiagram struct contains information about a \f$1\ton\f$ decay
 * that has been automatically generated.
 */
struct NBDiagram {

  /**
   * Enumeration for the channel type
   */
  enum Channel {UNDEFINED = -1, 
		TBchannel23=0, TBchannel13=1, 
		TBchannel12=2, TBfourPoint=3};

  /**
   * Standard Constructor
   */
  NBDiagram() : channelType(UNDEFINED) {}

  /**
   * Constructor taking a prototype vertex as the arguments*/
  NBDiagram(PrototypeVertexPtr proto);
  
  /**
   * Incoming particle
   */
  tPDPtr incoming;

  /**
   *  The type of channel
   */
  Channel channelType;

  /**
   *  Outgoing particles
   */
  mutable OrderedParticles outgoing;

  /**
   *  The vertices
   */
  list<pair<PDPtr,NBVertex> > vertices;

  /**
   *  The vertex for the parent branching
   */
  VertexBasePtr vertex;

  /** Store colour flow at \f$N_c=3\f$ information */
  mutable vector<CFPair> colourFlow;
  
  /** Store colour flow at \f$N_c=\infty\f$ information */
  mutable vector<CFPair> largeNcColourFlow;
};
//   /**
//    * Test whether this and x are the same decay
//    * @param x The other process to check
//    */
//   bool sameDecay(const TBDiagram & x) const {
//     if(ids[0] != x.ids[0]) return false;
//     bool used[4]={false,false,false,false};
//     for(unsigned int ix=1;ix<4;++ix) {
//       bool found=false;
//       for(unsigned int iy=1;iy<4;++iy) {
// 	if(used[iy]) continue;
// 	if(ids[ix]==x.ids[iy]) {
// 	  used[iy]=true;
// 	  found=true;
// 	  break;
// 	}
//       }
//       if(!found) return false;
//     }
//     return true;
//   }

/** 
 * Output operator to allow the structure to be persistently written
 * @param os The output stream
 * @param x The TBVertex 
 */
inline PersistentOStream & operator<<(PersistentOStream & os, 
				      const NBVertex  & x) {
  os << x.incoming << x.outgoing << x.vertices << x.vertex;
  return os;
}
  
/** 
 * Input operator to allow persistently written data to be read in
 * @param is The input stream
 * @param x The TBVertex 
 */
inline PersistentIStream & operator>>(PersistentIStream & is,
				      NBVertex & x) {
  is >> x.incoming >> x.outgoing >> x.vertices >> x.vertex;
  return is;
}

/** 
 * Output operator to allow the structure to be persistently written
 * @param os The output stream
 * @param x The TBDiagram 
 */
inline PersistentOStream & operator<<(PersistentOStream & os, 
				      const NBDiagram  & x) {
  os << x.incoming << oenum(x.channelType) << x.outgoing << x.vertices << x.vertex
     << x.colourFlow << x.largeNcColourFlow;
  return os;
}
  
/** 
 * Input operator to allow persistently written data to be read in
 * @param is The input stream
 * @param x The TBDiagram 
 */
inline PersistentIStream & operator>>(PersistentIStream & is,
				      NBDiagram & x) {
  is >> x.incoming >> ienum(x.channelType) >> x.outgoing >> x.vertices >> x.vertex
     >> x.colourFlow >> x.largeNcColourFlow;
  return is;
}

}

#endif /* HERWIG_TwoBodyPrototype_H */
