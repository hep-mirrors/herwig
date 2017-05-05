// -*- C++ -*-
//
// Tree2toNGenerator.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Tree2toNGenerator class.
//

#include "Tree2toNGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Utilities/StringUtils.h"

using namespace Herwig;

Tree2toNGenerator::Tree2toNGenerator() 
  : maxOrderGs(0), maxOrderGem(0), prepared(false) {}

Tree2toNGenerator::~Tree2toNGenerator() {}

IBPtr Tree2toNGenerator::clone() const {
  return new_ptr(*this);
}

IBPtr Tree2toNGenerator::fullclone() const {
  return new_ptr(*this);
}

vector<Ptr<Tree2toNDiagram>::ptr> Tree2toNGenerator::
generate(const PDVector& legs,
	 unsigned int orderInGs,
	 unsigned int orderInGem) {

  vector<Ptr<Tree2toNDiagram>::ptr> res;

  list<vector<Vertex> > prog = clusterAll(legs,orderInGs,orderInGem);

  int count = 1;
  for ( auto  & d : prog ) {
    assert(d.size() == 1);
    Tree2toNDiagram diag = d.front().generate(count);


    if ( !spaceLikeAllowed.empty() ) {
      map<tcPDPtr,int> counts;
      for ( int k = 1; k < diag.nSpace()-1; ++k )
	counts[diag.allPartons()[k]] += 1;
      for ( auto & m : spaceLikeAllowed ) {
	    m.reset();
        for ( auto const & c : counts )
	      m.add(c.first,c.second);
      }
      bool failed = false;
      for ( auto const & m : spaceLikeAllowed) {
	     if ( !m.check() ) {
	       failed = true;
	       break;
	     }
      }
      if ( failed )
	continue;
    }

    if ( !timeLikeAllowed.empty() ) {
      map<tcPDPtr,int> counts;
      int all = diag.allPartons().size();
      for ( int k = diag.nSpace(); k < all; ++k ) {
	if ( diag.children(k).first < 0 )
	  continue;
	counts[diag.allPartons()[k]] += 1;
      }
      for ( auto  & m : timeLikeAllowed ) {
	    m.reset();
        for ( auto const & c : counts )
	       m.add(c.first,c.second);
      }
      bool failed = false;
      for ( auto const &  m : timeLikeAllowed ) {
	if ( !m.check() ) {
	  failed = true;
	  break;
	}
      }
      if ( failed )
	continue;
    }


    bool internalVeto = false;
    set<int> external;
    int nex = diag.partons().size();
    for ( int i = 0; i < nex; ++i ) {
      external.insert(diag.diagramId(i));
    }
    int n = diag.allPartons().size();
    for ( int i = 0; i < n; ++i ) {
      if ( external.find(i) != external.end() )
	continue;
      if ( find(excludeInternal().begin(), excludeInternal().end(), diag.allPartons()[i])
	   != excludeInternal().end() ) {
	internalVeto = true;
	break;
      }
    }
    if ( internalVeto )
      continue;
    bool gotit = false;
    for ( auto const & d : res) {
      map<int,int> checkPermutation;
      if ( diag.isSame(d,checkPermutation) ) {
	gotit = true;
        for ( auto const & p : checkPermutation )
	  if ( p.first != p.second )
	    gotit = false;
	if ( gotit )
	  break;
      }
    }
    if ( !gotit ) {
      res.push_back(new_ptr(diag));
      ++count;
    }
  }

  return res;

}

list<vector<Tree2toNGenerator::Vertex> > Tree2toNGenerator::
cluster(const vector<Tree2toNGenerator::Vertex>& children,
	unsigned int orderInGs,
	unsigned int orderInGem) const {

  list<vector<Vertex> > res;

  bool externalCluster = children[1].externalId != -1;

  if ( children.size() == 3 ) {
    for ( auto const & v : theVertices ) {
	if ( v->getNpoint() != 3 )
	  continue;
	if ( find(theExcludeVertices.begin(), theExcludeVertices.end(), v) !=
	     theExcludeVertices.end() )
	  continue;
	bool noMatch =
	  v->orderInGs() != orderInGs ||
	  v->orderInGem() != orderInGem ||
	  !v->isIncoming(children[0].parent);
	long idij = children[0].parent->id();
	long idi  = children[2].parent->id();
	long idj  = children[1].parent->id();
	if ( externalCluster && children[1].parent->CC() )
	  idj = -idj;
	if ( children[0].parent->CC() )
	  idij = -idij;
	if ( !externalCluster )
	  noMatch |=
	    !v->isOutgoing(children[1].parent) ||
	    !v->isOutgoing(children[2].parent);
	else
	  noMatch |=
	    !v->isIncoming(children[1].parent) ||
	    !v->isOutgoing(children[2].parent);
	noMatch |=
	  !( v->allowed(idij,idi,idj) ||
	     v->allowed(idj,idij,idi) ||
	     v->allowed(idi,idj,idij) ||
	     v->allowed(idij,idj,idi) ||
	     v->allowed(idi,idij,idj) ||
	     v->allowed(idj,idi,idij) );
	if ( noMatch )
	  continue;
	Vertex last;
	last.spacelike = true;
	last.parent = children[0].parent;
	last.externalId = 0;
	last.children.push_back(children[1]);
	last.children.push_back(children[2]);
	res.push_back(vector<Vertex>(1,last));
	// only one possible
	break;
      }
      return res;
  }

  // spacelike clusterings (cluster on second one)
  for ( size_t i = 2; i < children.size(); ++i ) {
    for ( auto const & v : theVertices ) {
      if ( v->getNpoint() != 3 )
	continue;
      if ( find(theExcludeVertices.begin(), theExcludeVertices.end(), v) !=
	   theExcludeVertices.end() )
	continue;
      bool noMatch = false;
      noMatch |=
	v->orderInGs() != orderInGs ||
	v->orderInGem() != orderInGem;
      if ( !externalCluster )
	noMatch |=
	  !v->isOutgoing(children[1].parent) ||
	  !v->isOutgoing(children[i].parent);
      else
	noMatch |=
	  !v->isIncoming(children[1].parent) ||
	  !v->isOutgoing(children[i].parent);
      if ( noMatch )
	continue;
      long idi = children[i].parent->id();
      long idj = children[1].parent->id();
      if ( externalCluster && children[1].parent->CC() )
	idj = -idj;
      for ( set<tPDPtr>::const_iterator pij =
	      v->outgoing().begin(); pij != v->outgoing().end() ; ++pij ) {
	long idij = (**pij).id();
	if ( v->allowed(idij,idi,idj) ||
	     v->allowed(idj,idij,idi) ||
	     v->allowed(idi,idj,idij) ||
	     v->allowed(idij,idj,idi) ||
	     v->allowed(idi,idij,idj) ||
	     v->allowed(idj,idi,idij) ) {
	  PDPtr dij = (**pij).CC() ? (**pij).CC() : *pij;      
	  vector<Vertex> cled;
	  for ( size_t k = 0; k < children.size(); ++k ) {
	    if ( k != 1 && k != i )
	      cled.push_back(children[k]);
	    if ( k == 1 ) {
	      Vertex merge;
	      merge.children.push_back(children[1]);
	      merge.children.push_back(children[i]);
	      merge.parent = dij;
	      merge.spacelike = true;
	      cled.push_back(merge);
	    }
	    if ( k == i )
	      continue;
	  }
	  res.push_back(cled);
	}
      }
    }
  }

  // timelike clusterings
  for ( size_t i = 2; i < children.size(); ++i ) {
    for ( size_t j = i+1; j < children.size(); ++j ) {
      for ( auto const & v : theVertices ) {
	if ( v->getNpoint() != 3 )
	  continue;
	if ( find(theExcludeVertices.begin(), theExcludeVertices.end(), v) !=
	     theExcludeVertices.end() )
	  continue;
	if ( v->orderInGs() != orderInGs ||
	     v->orderInGem() != orderInGem ||
	     !v->isOutgoing(children[i].parent) ||
	     !v->isOutgoing(children[j].parent) )
	  continue;
	long idi = children[i].parent->id();
	long idj = children[j].parent->id();
	for ( set<tPDPtr>::const_iterator pij =
		v->outgoing().begin(); pij != v->outgoing().end() ; ++pij ) {
	  long idij = (**pij).id();
	  if ( v->allowed(idij,idi,idj) ||
	       v->allowed(idj,idij,idi) ||
	       v->allowed(idi,idj,idij) ||
	       v->allowed(idij,idj,idi) ||
	       v->allowed(idi,idij,idj) ||
	       v->allowed(idj,idi,idij) ) {
	    PDPtr dij = (**pij).CC() ? (**pij).CC() : *pij;
	    vector<Vertex> cled;
	    for ( size_t k = 0; k < children.size(); ++k ) {
	      if ( k != i && k != j )
		cled.push_back(children[k]);
	      if ( k == i ) {
		Vertex merge;
		merge.children.push_back(children[i]);
		merge.children.push_back(children[j]);
		merge.parent = dij;
		merge.spacelike = false;
		cled.push_back(merge);
	      }
	      if ( k == j )
		continue;
	    }
	    res.push_back(cled);
	  }
	}
      }
    }
  }

  return res;

}

list<vector<Tree2toNGenerator::Vertex> > Tree2toNGenerator::
clusterAll(const list<vector<Tree2toNGenerator::Vertex> >& current,
	   unsigned int orderInGs,
	   unsigned int orderInGem) const {

  list<vector<Vertex> > res;
  for ( list<vector<Vertex> >::const_iterator c = current.begin();
	c != current.end(); ++c ) {
    if ( c->size() == 1 ) {
      if ( orderInGs == 0 && orderInGem == 0 )
	res.push_back(*c);
      continue;
    }
    for ( unsigned int gs = 0; gs <= maxOrderGs; ++gs )
      for ( unsigned int gem = 0; gem <= maxOrderGem; ++gem ) {
	if ( gs == 0 && gem == 0 )
	  continue;
	if ( gs > orderInGs || gem > orderInGem )
	  continue;
	list<vector<Vertex> > next = cluster(*c,gs,gem);
	if ( next.empty() )
	  continue;
	list<vector<Vertex> > cled = clusterAll(next,orderInGs-gs,orderInGem-gem);
	copy(cled.begin(),cled.end(),back_inserter(res));
      }
  }

  return res;

}

list<vector<Tree2toNGenerator::Vertex> > Tree2toNGenerator::
clusterAll(const PDVector& external,
	   unsigned int orderInGs,
	   unsigned int orderInGem) {

  if ( !prepared ) {
    for ( auto & v : theVertices ) {
      if ( find(theExcludeVertices.begin(), theExcludeVertices.end(), v) !=
	   theExcludeVertices.end() )
	continue;
      v->init();
      maxOrderGs = max(maxOrderGs,v->orderInGs());
      maxOrderGem = max(maxOrderGem,v->orderInGem());
    }
    for ( auto & m : spaceLikeAllowed)
      m.rebind(this);
    
    for ( auto & m : timeLikeAllowed )
      m.rebind(this);
    
    prepared = true;
  }

  vector<Vertex> legs;

  for ( unsigned int k = 0; k < external.size(); ++k ) {
    Vertex v;
    v.parent = external[k];
    v.externalId = k;
    v.spacelike = k < 2;
    legs.push_back(v);
  }

  list<vector<Vertex> > firstlegs;
  firstlegs.push_back(legs);

  return clusterAll(firstlegs,orderInGs,orderInGem);

}

string Tree2toNGenerator::doSpaceLikeRange(string range) {
  if ( theRestrictLines.empty() )
    return "No particle data specified to restrict internal lines.";
  vector<string> bounds = StringUtils::split(range);
  if ( bounds.empty() || bounds.size() > 2 )
    return "Need to specify a minimum, or a minimum and maximum number of internal lines.";
  pair<int,int> irange(0,-1);
  istringstream in1(bounds[0]);
  in1 >> irange.first;
  if ( bounds.size() == 2 ) {
    istringstream in2(bounds[1]);
    in2 >> irange.second;
  } else {
    irange.second = irange.first;
  }
  if ( irange.second >= 0 && irange.first > irange.second )
    return "invalid range specified";
  spaceLikeAllowed.push_back(LineMatcher(theRestrictLines,irange));
  return "";
}

string Tree2toNGenerator::doTimeLikeRange(string range) {
  if ( theRestrictLines.empty() )
    return "No particle data specified to restrict internal lines.";
  vector<string> bounds = StringUtils::split(range);
  if ( bounds.empty() || bounds.size() > 2 )
    return "Need to specify a minimum, or a minimum and maximum number of internal lines.";
  pair<int,int> irange(0,-1);
  istringstream in1(bounds[0]);
  in1 >> irange.first;
  if ( bounds.size() == 2 ) {
    istringstream in2(bounds[1]);
    in2 >> irange.second;
  } else {
    irange.second = irange.first;
  }
  if ( irange.second >= 0 && irange.first > irange.second )
    return "invalid range specified";
  timeLikeAllowed.push_back(LineMatcher(theRestrictLines,irange));
  return "";
}

string Tree2toNGenerator::doClearRestrictLines(string) {
  theRestrictLines.clear();
  return "";
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void Tree2toNGenerator::persistentOutput(PersistentOStream & os) const {
  os << theVertices << theExcludeInternal << maxOrderGs << maxOrderGem << prepared
     << theExcludeVertices << spaceLikeAllowed << timeLikeAllowed;
}

void Tree2toNGenerator::persistentInput(PersistentIStream & is, int) {
  is >> theVertices >> theExcludeInternal >> maxOrderGs >> maxOrderGem >> prepared
     >> theExcludeVertices >> spaceLikeAllowed >> timeLikeAllowed;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<Tree2toNGenerator,HandlerBase>
  describeHerwigTree2toNGenerator("Herwig::Tree2toNGenerator", "Herwig.so");

void Tree2toNGenerator::Init() {

  static ClassDocumentation<Tree2toNGenerator> documentation
    ("Generate Tree2toNDiagrams for a given process.");


  static RefVector<Tree2toNGenerator,Helicity::VertexBase> interfaceVertices
    ("Vertices",
     "All vertices to consider.",
     &Tree2toNGenerator::theVertices, -1, false, false, true, false, false);

  static RefVector<Tree2toNGenerator,Helicity::VertexBase> interfaceExcludeVertices
    ("ExcludeVertices",
     "The vertices to exclude.",
     &Tree2toNGenerator::theExcludeVertices, -1, false, false, true, false, false);

  static RefVector<Tree2toNGenerator,ParticleData> interfaceExcludeInternal
    ("ExcludeInternal",
     "Particles to be exluded from becoming internal lines.",
     &Tree2toNGenerator::theExcludeInternal, -1, false, false, true, false, false);

  static RefVector<Tree2toNGenerator,ParticleData> interfaceRestrictLines
    ("RestrictLines",
     "Particles to be exluded from becoming internal lines.",
     &Tree2toNGenerator::theRestrictLines, -1, false, false, true, false, false);

  static Command<Tree2toNGenerator> interfaceSpaceLikeRange
    ("SpaceLikeRange",
     "Limit the number of spacelike occurences of the specified particle.",
     &Tree2toNGenerator::doSpaceLikeRange, false);

  static Command<Tree2toNGenerator> interfaceTimeLikeRange
    ("TimeLikeRange",
     "Limit the number of timelike occurences of the specified particle.",
     &Tree2toNGenerator::doTimeLikeRange, false);

  static Command<Tree2toNGenerator> interfaceClearRestrictLines
    ("ClearRestrictLines",
     "Clear the container of lines to be considered for restrictions.",
     &Tree2toNGenerator::doClearRestrictLines, false);

}

