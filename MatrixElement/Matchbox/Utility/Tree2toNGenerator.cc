// -*- C++ -*-
//
// Tree2toNGenerator.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
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
#include "Herwig/MatrixElement/Matchbox/Utility/DiagramDrawer.h"

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

  list<vector<Vertex> > prog =
    clusterAll(legs,orderInGs,orderInGem);

  int count = 1;
  for ( list<vector<Vertex> >::iterator d =
	  prog.begin(); d != prog.end(); ++d ) {
    assert(d->size() == 1);
    Tree2toNDiagram diag = d->front().generate(count);


    if ( !spaceLikeAllowed.empty() ) {
      map<tcPDPtr,int> counts;
      for ( int k = 1; k < diag.nSpace()-1; ++k ) {
	if ( counts.find(diag.allPartons()[k]) != counts.end() ) {
	  counts[diag.allPartons()[k]] += 1;
	} else {
	  counts[diag.allPartons()[k]] += 1;
	}
      }
      for ( vector<LineMatcher>::iterator m = spaceLikeAllowed.begin();
	    m != spaceLikeAllowed.end(); ++m ) {
	m->reset();
	for ( map<tcPDPtr,int>::const_iterator c = counts.begin();
	      c != counts.end(); ++c )
	  m->add(c->first,c->second);
      }
      bool failed = false;
      for ( vector<LineMatcher>::iterator m = spaceLikeAllowed.begin();
	    m != spaceLikeAllowed.end(); ++m ) {
	if ( !m->check() ) {
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
	if ( diag.children(k).empty() )
	  continue;
	assert(diag.children(k)[0]>=0);
	if ( counts.find(diag.allPartons()[k]) != counts.end() ) {
	  counts[diag.allPartons()[k]] += 1;
	} else {
	  counts[diag.allPartons()[k]] += 1;
	}
      }
      for ( vector<LineMatcher>::iterator m = timeLikeAllowed.begin();
	    m != timeLikeAllowed.end(); ++m ) {
	m->reset();
	for ( map<tcPDPtr,int>::const_iterator c = counts.begin();
	      c != counts.end(); ++c )
	  m->add(c->first,c->second);
      }
      bool failed = false;
      for ( vector<LineMatcher>::iterator m = timeLikeAllowed.begin();
	    m != timeLikeAllowed.end(); ++m ) {
	if ( !m->check() ) {
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
    for ( vector<Ptr<Tree2toNDiagram>::ptr>::const_iterator
	    d = res.begin(); d != res.end(); ++d ) {
      map<int,int> checkPermutation;
      if ( diag.isSame(*d,checkPermutation) ) {
	gotit = true;
//        cout<<"\nis same:\n";
//         DiagramDrawer::drawDiag(cout,diag);
	for ( map<int,int>::const_iterator p = checkPermutation.begin();
	      p != checkPermutation.end(); ++p ){
	  if ( p->first != p->second )
	    gotit = false;
        }
	if ( gotit )
	  break;
      }
    }
    if ( !gotit && diag.external().size()==legs.size()) {
      res.push_back(new_ptr(diag));
      ++count;
    }else if(!gotit){
      cout<<"\n external "<<diag.external().size()<<" legs "<<legs.size();
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
  bool externalCluster2 = children[2].externalId != -1;
 

  if ( children.size() == 3 ) {
      for ( VertexVector::const_iterator v = theVertices.begin();
	    v != theVertices.end(); ++v ) {
	if ( (**v).getNpoint() != 3 )
	  continue;
	if ( find(theExcludeVertices.begin(), theExcludeVertices.end(), *v) !=
	     theExcludeVertices.end() )
	  continue;
	bool noMatch =
	  (**v).orderInGs() != orderInGs ||
	  (**v).orderInGem() != orderInGem ||
	  !(**v).isIncoming(children[0].parent);
	long idij = children[0].parent->id();
	long idi  = children[2].parent->id();
	long idj  = children[1].parent->id();
	if ( externalCluster && children[1].parent->CC() )
	  idj = -idj;
	if ( children[0].parent->CC() )
	  idij = -idij;
	if ( !externalCluster )
	  noMatch |=
	    !(**v).isOutgoing(children[1].parent) ||
	    !(**v).isOutgoing(children[2].parent);
	else
	  noMatch |=
	    !(**v).isIncoming(children[1].parent) ||
	    !(**v).isOutgoing(children[2].parent);
	noMatch |=
	  !( (**v).allowed(idij,idi,idj) ||
	     (**v).allowed(idj,idij,idi) ||
	     (**v).allowed(idi,idj,idij) ||
	     (**v).allowed(idij,idj,idi) ||
	     (**v).allowed(idi,idij,idj) ||
	     (**v).allowed(idj,idi,idij) );
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
  
  
  if ( children.size() == 4 ) {
    for ( VertexVector::const_iterator v = theVertices.begin();
         v != theVertices.end(); ++v ) {
      if ( (**v).getNpoint() != 4 )
        continue;
      assert(false);
      if ( find(theExcludeVertices.begin(), theExcludeVertices.end(), *v) !=
          theExcludeVertices.end() )
        continue;
      bool noMatch =
      (**v).orderInGs() != orderInGs ||
      (**v).orderInGem() != orderInGem ||
      !(**v).isIncoming(children[0].parent);
      cout<<"\nwhat1";
      long idijk = children[0].parent->id();
      long idi  = children[3].parent->id();
      long idj  = children[2].parent->id();
      long idk  = children[1].parent->id();
      if ( externalCluster && children[1].parent->CC() )
        idk = -idk;
      if ( externalCluster2 && children[2].parent->CC() )
        idj = -idj;
      if ( children[0].parent->CC() )
        idijk = -idijk;
      cout<<"\nwhat2";
      if ( !externalCluster )
        noMatch |=
        !(**v).isOutgoing(children[1].parent) ||
        !(**v).isOutgoing(children[2].parent) ||
        !(**v).isOutgoing(children[3].parent);
      else
        noMatch |=
        !((**v).allowed(idijk,idi,idj,idk) ||
          (**v).allowed(idijk,idi,idk,idj) ||
          (**v).allowed(idijk,idj,idi,idk) ||
          (**v).allowed(idijk,idj,idk,idi) ||
          (**v).allowed(idijk,idk,idi,idj) ||
          (**v).allowed(idijk,idk,idj,idi) ||
          (**v).allowed(idi,idijk,idj,idk) ||
          (**v).allowed(idi,idijk,idk,idj) ||
          (**v).allowed(idi,idj,idijk,idk) ||
          (**v).allowed(idi,idj,idk,idijk) ||
          (**v).allowed(idi,idk,idijk,idj) ||
          (**v).allowed(idi,idk,idj,idijk) ||
          (**v).allowed(idj,idijk,idi,idk) ||
          (**v).allowed(idj,idijk,idk,idi) ||
          (**v).allowed(idj,idi,idijk,idk) ||
          (**v).allowed(idj,idi,idk,idijk) ||
          (**v).allowed(idj,idk,idijk,idi) ||
          (**v).allowed(idj,idk,idi,idijk) ||
          (**v).allowed(idk,idijk,idi,idj) ||
          (**v).allowed(idk,idijk,idj,idi) ||
          (**v).allowed(idk,idi,idijk,idj) ||
          (**v).allowed(idk,idi,idj,idijk) ||
          (**v).allowed(idk,idj,idijk,idi) ||
          (**v).allowed(idk,idj,idi,idijk) );
      if ( noMatch )
        continue;
      cout<<"\nfound one";
      Vertex last;
      last.spacelike = true;
      last.parent = children[0].parent;
      last.externalId = 0;
      last.children.push_back(children[1]);
      last.children.push_back(children[2]);
      last.children.push_back(children[3]);
      res.push_back(vector<Vertex>(1,last));
        // only one possible
      break;
    }
  }
  
    // spacelike clusterings (cluster on second one)
  for ( size_t i = 2; i < children.size(); ++i ) {
    for ( VertexVector::const_iterator v = theVertices.begin();
	  v != theVertices.end(); ++v ) {
      if ( (**v).getNpoint() != 3 )
	continue;
      if ( find(theExcludeVertices.begin(), theExcludeVertices.end(), *v) !=
	   theExcludeVertices.end() )
	continue;
      bool noMatch = false;
      noMatch |=
	(**v).orderInGs() != orderInGs ||
	(**v).orderInGem() != orderInGem;
      if ( !externalCluster )
	noMatch |=
	  !(**v).isOutgoing(children[1].parent) ||
	  !(**v).isOutgoing(children[i].parent);
      else
	noMatch |=
	  !(**v).isIncoming(children[1].parent) ||
	  !(**v).isOutgoing(children[i].parent);
      if ( noMatch )
	continue;
      long idi = children[i].parent->id();
      long idj = children[1].parent->id();
      if ( externalCluster && children[1].parent->CC() )
	idj = -idj;
      for ( set<tPDPtr>::const_iterator pij =
	      (**v).outgoing().begin(); pij != (**v).outgoing().end() ; ++pij ) {
	long idij = (**pij).id();
	if ( (**v).allowed(idij,idi,idj) ||
	     (**v).allowed(idj,idij,idi) ||
	     (**v).allowed(idi,idj,idij) ||
	     (**v).allowed(idij,idj,idi) ||
	     (**v).allowed(idi,idij,idj) ||
	     (**v).allowed(idj,idi,idij) ) {
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
  
    // spacelike clusterings (cluster on second one)  for 4 point functions
  if (children.size() != 4)
    for ( size_t i = 2; i < children.size(); ++i ) {
      for ( size_t j = i+1; j < children.size(); ++j ) {
        for ( VertexVector::const_iterator v = theVertices.begin();
             v != theVertices.end(); ++v ) {
          if ( (**v).getNpoint() != 4 )
            continue;
          assert(false);
          if ( find(theExcludeVertices.begin(), theExcludeVertices.end(), *v) !=
              theExcludeVertices.end() )
            continue;
          bool noMatch = false;
          noMatch |= (**v).orderInGs() != orderInGs ||        (**v).orderInGem() != orderInGem;
          if ( !externalCluster )
            noMatch |=
            !(**v).isOutgoing(children[1].parent) ||          !(**v).isOutgoing(children[i].parent)||
            !(**v).isOutgoing(children[j].parent);
          else
            noMatch |=
            !(**v).isIncoming(children[1].parent) ||          !(**v).isOutgoing(children[i].parent) ||
            !(**v).isOutgoing(children[j].parent);
          if ( noMatch )
            continue;
          long idi = children[i].parent->id();
          long idj = children[j].parent->id();
          long idk = children[1].parent->id();
          if ( externalCluster && children[1].parent->CC() )
            idk = -idk;
          for ( set<tPDPtr>::const_iterator pijk =
               (**v).outgoing().begin(); pijk != (**v).outgoing().end() ; ++pijk ) {
            long idijk = (**pijk).id();
            if ( ((**v).allowed(idijk,idi,idj,idk) || (**v).allowed(idijk,idi,idk,idj) || (**v).allowed(idijk,idj,idi,idk) || (**v).allowed(idijk,idj,idk,idi) || (**v).allowed(idijk,idk,idi,idj) || (**v).allowed(idijk,idk,idj,idi) || (**v).allowed(idi,idijk,idj,idk) || (**v).allowed(idi,idijk,idk,idj) || (**v).allowed(idi,idj,idijk,idk) || (**v).allowed(idi,idj,idk,idijk) || (**v).allowed(idi,idk,idijk,idj) ||  (**v).allowed(idi,idk,idj,idijk) || (**v).allowed(idj,idijk,idi,idk) || (**v).allowed(idj,idijk,idk,idi) || (**v).allowed(idj,idi,idijk,idk) || (**v).allowed(idj,idi,idk,idijk) || (**v).allowed(idj,idk,idijk,idi) || (**v).allowed(idj,idk,idi,idijk) ||  (**v).allowed(idk,idijk,idi,idj) ||  (**v).allowed(idk,idijk,idj,idi) ||
                  (**v).allowed(idk,idi,idijk,idj) ||  (**v).allowed(idk,idi,idj,idijk) || (**v).allowed(idk,idj,idijk,idi) || (**v).allowed(idk,idj,idi,idijk) ) ) {
              PDPtr dijk = (**pijk).CC() ? (**pijk).CC() : *pijk;
              vector<Vertex> cled;
              for ( size_t l = 0; l < children.size(); ++l ) {
                if ( l != 1 && l != i && l != j )
                  cled.push_back(children[l]);
                if ( l == 1 ) {
                  Vertex merge;
                  merge.children.push_back(children[1]);
                  merge.children.push_back(children[i]);
                  merge.children.push_back(children[j]);
                  merge.parent = dijk;
                  merge.spacelike = true;
                  cled.push_back(merge);
                }
                if ( l == i ||l == j)
                  continue;
              }
              cout<<"\nwhat55 no";
              res.push_back(cled);
            }
          }
        }
      }
    }
  
  
    // timelike clusterings
  for ( size_t i = 2; i < children.size(); ++i ) {
    for ( size_t j = i+1; j < children.size(); ++j ) {
      for ( VertexVector::const_iterator v = theVertices.begin();
	    v != theVertices.end(); ++v ) {
	if ( (**v).getNpoint() != 3 )
	  continue;
	if ( find(theExcludeVertices.begin(), theExcludeVertices.end(), *v) !=
	     theExcludeVertices.end() )
	  continue;
	if ( (**v).orderInGs() != orderInGs ||
	     (**v).orderInGem() != orderInGem ||
	     !(**v).isOutgoing(children[i].parent) ||
	     !(**v).isOutgoing(children[j].parent) )
	  continue;
	long idi = children[i].parent->id();
	long idj = children[j].parent->id();
	for ( set<tPDPtr>::const_iterator pij =
		(**v).outgoing().begin(); pij != (**v).outgoing().end() ; ++pij ) {
	  long idij = (**pij).id();
	  if ( (**v).allowed(idij,idi,idj) ||
	       (**v).allowed(idj,idij,idi) ||
	       (**v).allowed(idi,idj,idij) ||
	       (**v).allowed(idij,idj,idi) ||
	       (**v).allowed(idi,idij,idj) ||
	       (**v).allowed(idj,idi,idij) ) {
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
  
  
    // timelike clusterings  4 points
  for ( size_t i = 2; i < children.size(); ++i ) {
    for ( size_t j = i+1; j < children.size(); ++j ) {
      for ( size_t k = j+1; k < children.size(); ++k ) {
        for ( VertexVector::const_iterator v = theVertices.begin();
             v != theVertices.end(); ++v ) {
          if ( (**v).getNpoint() != 4 )
            continue;
          assert(false);
          if ( find(theExcludeVertices.begin(), theExcludeVertices.end(), *v) !=
              theExcludeVertices.end() )
            continue;
          if ( (**v).orderInGs() != orderInGs ||
              (**v).orderInGem() != orderInGem ||
              !(**v).isOutgoing(children[i].parent) ||
              !(**v).isOutgoing(children[j].parent) ||
              !(**v).isOutgoing(children[k].parent))
            continue;
          long idi = children[i].parent->id();
          long idj = children[j].parent->id();
          long idk = children[k].parent->id();
          for ( set<tPDPtr>::const_iterator pijk =
               (**v).outgoing().begin(); pijk != (**v).outgoing().end() ; ++pijk ) {
            long idijk = (**pijk).id();
            if ( ((**v).allowed(idijk,idi,idj,idk) ||                  (**v).allowed(idijk,idi,idk,idj) ||                  (**v).allowed(idijk,idj,idi,idk) ||                  (**v).allowed(idijk,idj,idk,idi) ||
                  (**v).allowed(idijk,idk,idi,idj) ||                  (**v).allowed(idijk,idk,idj,idi) ||                  (**v).allowed(idi,idijk,idj,idk) ||                  (**v).allowed(idi,idijk,idk,idj) ||
                  (**v).allowed(idi,idj,idijk,idk) ||                  (**v).allowed(idi,idj,idk,idijk) ||                  (**v).allowed(idi,idk,idijk,idj) ||                  (**v).allowed(idi,idk,idj,idijk) ||
                  (**v).allowed(idj,idijk,idi,idk) ||                  (**v).allowed(idj,idijk,idk,idi) ||                  (**v).allowed(idj,idi,idijk,idk) ||                  (**v).allowed(idj,idi,idk,idijk) ||
                  (**v).allowed(idj,idk,idijk,idi) ||                  (**v).allowed(idj,idk,idi,idijk) ||                  (**v).allowed(idk,idijk,idi,idj) ||                  (**v).allowed(idk,idijk,idj,idi) ||
                  (**v).allowed(idk,idi,idijk,idj) ||                  (**v).allowed(idk,idi,idj,idijk) ||                  (**v).allowed(idk,idj,idijk,idi) ||                  (**v).allowed(idk,idj,idi,idijk) ) ) {
              PDPtr dijk = (**pijk).CC() ? (**pijk).CC() : *pijk;
              vector<Vertex> cled;
              for ( size_t l = 0; l < children.size(); ++l ) {
                if ( l != i && l != j && l != k )
                  cled.push_back(children[l]);
                if ( l == i ) {
                  Vertex merge;
                  merge.children.push_back(children[i]);
                  merge.children.push_back(children[j]);
                  merge.children.push_back(children[k]);
                  merge.parent = dijk;
                  merge.spacelike = false;
                  cled.push_back(merge);
                }
                if ( k == j ||k == i )
                  continue;
              }
              
              cout<<"\nwhat666 no";
              res.push_back(cled);
            }
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
    for ( VertexVector::iterator v = theVertices.begin();
	  v != theVertices.end(); ++v ) {
      if ( find(theExcludeVertices.begin(), theExcludeVertices.end(), *v) !=
	   theExcludeVertices.end() )
	continue;
      (**v).init();
      maxOrderGs = max(maxOrderGs,(**v).orderInGs());
      maxOrderGem = max(maxOrderGem,(**v).orderInGem());
    }
    for ( vector<LineMatcher>::iterator m = spaceLikeAllowed.begin();
	  m != spaceLikeAllowed.end(); ++m ) {
      m->rebind(this);
    }
    for ( vector<LineMatcher>::iterator m = timeLikeAllowed.begin();
	  m != timeLikeAllowed.end(); ++m ) {
      m->rebind(this);
    }
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



/**
 * Update diagram returning a map of external ids to diagram id
 * parents.
 */
void Tree2toNGenerator::Vertex::update(Tree2toNDiagram& diag,
                                       map<int,pair<int,PDPtr> >& outgoing,
                                       int& lastUsed) {
  if ( externalId == 0 ) {
    assert(lastUsed==0);
    ++lastUsed;
    diag.operator,(parent);
    children[0].parentId = lastUsed;
    children[1].parentId = lastUsed;
    if (children.size()==3)children[2].parentId = lastUsed;
    children[0].update(diag,outgoing,lastUsed);
    children[1].update(diag,outgoing,lastUsed);
    if (children.size()==3)children[2].update(diag,outgoing,lastUsed);
    for ( map<int,pair<int,PDPtr> >::iterator out =
         outgoing.begin(); out != outgoing.end(); ++out ) {
      diag.operator,(out->second.first);
      diag.operator,(out->second.second);
    }
    return;
  }








  if ( spacelike ) {
    ++lastUsed;
    diag.operator,(parent);
    if ( externalId == 1 )
      return;
    children[0].parentId = lastUsed;
    children[1].parentId = lastUsed;
    if (children.size()==3)children[2].parentId = lastUsed;
    children[0].update(diag,outgoing,lastUsed);
    children[1].update(diag,outgoing,lastUsed);
    if (children.size()==3)children[2].update(diag,outgoing,lastUsed);
    return;
  }
  if ( children.empty() ) {
    outgoing[externalId] =
    make_pair(parentId,parent);
    return;
  }
  diag.operator,(parentId);
  diag.operator,(parent);
  ++lastUsed; 
  children[0].parentId = lastUsed;
  children[1].parentId = lastUsed;
  if (children.size()==3)children[2].parentId = lastUsed;
  children[0].update(diag,outgoing,lastUsed);
  children[1].update(diag,outgoing,lastUsed);
  if (children.size()==3)children[2].update(diag,outgoing,lastUsed);
}

/**
 * Generate a diagram of given id.
 */
Tree2toNDiagram Tree2toNGenerator::Vertex::generate(int id) {
  int nsp = nspace();
  Tree2toNDiagram res(nsp);
  int diagid = 0;
  map<int,pair<int,PDPtr> > out;
  update(res,out,diagid);
  res.operator,(-id);
  return res;
}





void Tree2toNGenerator::Vertex::print(ostream& os, const string& prefix ) const {
  os << prefix << parent->PDGName()
	 << "[" << (spacelike ? "s" : "t") << "] (";
  if ( externalId < 0 )
    os << "x)\n";
  else
    os << externalId << ")\n";
  if ( !children.empty() ) {
    os << prefix << "|__\n";
    children[0].print(os,prefix + "|  ");
    os << prefix << "|__\n";
    children[1].print(os,prefix + "|  ");
  }
}

