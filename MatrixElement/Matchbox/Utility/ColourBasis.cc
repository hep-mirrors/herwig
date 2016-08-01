// -*- C++ -*-
//
// ColourBasis.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ColourBasis class.
//

#include "ColourBasis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Handlers/SamplerBase.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/MatrixElement/Matchbox/MatchboxFactory.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <iterator>
using std::ostream_iterator;

#include "DiagramDrawer.h"

using namespace Herwig;

using boost::numeric::ublas::trans;

// default gcc on SLC6 confuses this with std::conj, 
// use explicit namespacing in the code instead
//
// using boost::numeric::ublas::conj;

using boost::numeric::ublas::row;
using boost::numeric::ublas::column;
using boost::numeric::ublas::prod;

Ptr<MatchboxFactory>::tptr ColourBasis::factory() const {
  return theFactory;
}

void ColourBasis::factory(Ptr<MatchboxFactory>::tptr f) {
  theFactory = f;
}

ColourBasis::ColourBasis() 
  : theLargeN(false), didRead(false), didWrite(false), theSearchPath("") {}

ColourBasis::~ColourBasis() {
  for ( map<Ptr<Tree2toNDiagram>::tcptr,vector<ColourLines*> >::iterator cl =
	  theColourLineMap.begin(); cl != theColourLineMap.end(); ++cl ) {
    for ( vector<ColourLines*>::iterator c = cl->second.begin();
	  c != cl->second.end(); ++c ) {
      if ( *c )
	delete *c;
    }
  }
  theColourLineMap.clear();
}

void ColourBasis::clear() {
  theLargeN = false;
  theNormalOrderedLegs.clear();
  theIndexMap.clear();
  theScalarProducts.clear();
  theCharges.clear();
  theChargeNonZeros.clear();
  theCorrelators.clear();
  theFlowMap.clear();
  theColourLineMap.clear();
  theOrderingStringIdentifiers.clear();
  theOrderingIdentifiers.clear();
  didRead = false;
  didWrite = false;
  tmp.clear();
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

bool ColourBasis::colourConnected(const cPDVector& sub,
				  const vector<PDT::Colour>& basis,
				  const pair<int,bool>& i, 
				  const pair<int,bool>& j, 
				  size_t a) const {

  // translate process to basis ids
  map<cPDVector,map<size_t,size_t> >::const_iterator trans
    = indexMap().find(sub);
  assert(trans != indexMap().end());

  int idColoured = i.second ? j.first : i.first;
  idColoured = trans->second.find(idColoured)->second;
  int idAntiColoured = i.second ? i.first : j.first;
  idAntiColoured = trans->second.find(idAntiColoured)->second;

  return colourConnected(basis,idColoured,idAntiColoured,a);

}

const string& ColourBasis::orderingString(const cPDVector& sub, 
					  const map<size_t,size_t>& colourToAmplitude,
					  size_t tensorId) {

  map<size_t,string>& tensors = theOrderingStringIdentifiers[sub];
  if ( !tensors.empty() ) {
    assert(tensors.find(tensorId) != tensors.end());
    return tensors[tensorId];
  }

  const set<vector<size_t> >& xordering = ordering(sub,colourToAmplitude,tensorId);

  ostringstream os;
  os << "[";
  for ( set<vector<size_t> >::const_iterator t = xordering.begin();
	t != xordering.end(); ++t ) {
    os << "[";
    for ( vector<size_t>::const_iterator s = t->begin();
	  s != t->end(); ++s ) {
      os << *s << (s != --t->end() ? "," : "");
    }
    os << "]" << (t != --xordering.end() ? "," : "");
  }
  os << "]";

  tensors[tensorId] = os.str();
  return tensors[tensorId];

}

const set<vector<size_t> >& ColourBasis::ordering(const cPDVector& sub, 
						  const map<size_t,size_t>& colourToAmplitude,
						  size_t tensorId, size_t shift) {

  map<size_t,set<vector<size_t> > >& tensors = theOrderingIdentifiers[sub];
  if ( !tensors.empty() ) {
    assert(tensors.find(tensorId) != tensors.end());
    return tensors[tensorId];
  }

  const vector<PDT::Colour>& basisId = normalOrderedLegs(sub);

  map<size_t,vector<vector<size_t> > > labels = basisList(basisId);

  for ( map<size_t,vector<vector<size_t> > >::const_iterator t =
	  labels.begin(); t != labels.end(); ++t ) {
    set<vector<size_t> > xordering;
    for ( vector<vector<size_t> >::const_iterator s = t->second.begin();
	  s != t->second.end(); ++s ) {
      vector<size_t> crossed;
      for ( vector<size_t>::const_iterator l = s->begin();
	    l != s->end(); ++l ) {
	map<size_t,size_t>::const_iterator trans = 
	  colourToAmplitude.find(*l);
	assert(trans != colourToAmplitude.end());
	crossed.push_back(trans->second + shift);
      }
      xordering.insert(crossed);
    }
    tensors[t->first] = xordering;
  }

  assert(tensors.find(tensorId) != tensors.end());
  return tensors[tensorId];

}

vector<PDT::Colour> ColourBasis::normalOrderMap(const cPDVector& sub) {

  vector<PDT::Colour> allLegs = projectColour(sub);
  vector<PDT::Colour> legs = normalOrder(allLegs);

  if ( allLegs[0] == PDT::Colour3 )
    allLegs[0] = PDT::Colour3bar;
  else if ( allLegs[0] == PDT::Colour3bar )
    allLegs[0] = PDT::Colour3;
  if ( allLegs[1] == PDT::Colour3 )
    allLegs[1] = PDT::Colour3bar;
  else if ( allLegs[1] == PDT::Colour3bar )
    allLegs[1] = PDT::Colour3;

  if ( theIndexMap.find(sub) == theIndexMap.end() ) {
    map<size_t,size_t> trans;
    vector<PDT::Colour> checkLegs = legs;
    size_t n = checkLegs.size();
    for ( size_t i = 0; i < allLegs.size(); ++i ) {
      size_t j = 0;
      while ( checkLegs[j] != allLegs[i] ) {
	++j; if ( j == n ) break;
      }
      if ( j == n ) continue;
      trans[i] = j;
      checkLegs[j] = PDT::ColourUndefined;
    }  

    theIndexMap[sub] = trans;

  }

  return legs;

}

const vector<PDT::Colour>& ColourBasis::normalOrderedLegs(const cPDVector& sub) const {
  static vector<PDT::Colour> empty;
  map<cPDVector,vector<PDT::Colour> >::const_iterator n =
    theNormalOrderedLegs.find(sub);
  if ( n != theNormalOrderedLegs.end() )
    return n->second;
  return empty;
}

size_t ColourBasis::prepare(const cPDVector& sub,
			    bool noCorrelations) {

  vector<PDT::Colour> legs = normalOrderMap(sub);

  bool doPrepare = false;

  if ( theNormalOrderedLegs.find(sub) == theNormalOrderedLegs.end() )
    theNormalOrderedLegs[sub] = legs;
  
  if ( theScalarProducts.find(legs) == theScalarProducts.end() )
    doPrepare = true;

  if ( doPrepare )
    doPrepare = !readBasis(legs);

  size_t dim = doPrepare ? prepareBasis(legs) : theScalarProducts[legs].size1();

  if ( theCharges.find(legs) != theCharges.end() )
    return dim;

  if ( !doPrepare && noCorrelations )
    return dim;

  symmetric_matrix<double,upper>& sp = 
    theScalarProducts.insert(make_pair(legs,symmetric_matrix<double,upper>(dim,dim))).first->second;
  
  for ( size_t a = 0; a < dim; ++a )
    for ( size_t b = a; b < dim; ++b )
      sp(a,b) = scalarProduct(a,b,legs);

  if ( noCorrelations )
    return dim;

  vector<PDT::Colour> legsPlus = legs;
  legsPlus.push_back(PDT::Colour8);
  legsPlus = normalOrder(legsPlus);

  bool doPreparePlus = theScalarProducts.find(legsPlus) == theScalarProducts.end();

  size_t dimPlus = doPreparePlus ? prepareBasis(legsPlus) : theScalarProducts[legsPlus].size1();

  symmetric_matrix<double,upper>& spPlus = 
    doPreparePlus ?
    theScalarProducts.insert(make_pair(legsPlus,symmetric_matrix<double,upper>(dimPlus,dimPlus))).first->second :
    theScalarProducts[legsPlus];

  if ( doPreparePlus ) {
    for ( size_t a = 0; a < dimPlus; ++a )
      for ( size_t b = a; b < dimPlus; ++b )
	spPlus(a,b) = scalarProduct(a,b,legsPlus);
  }

  typedef map<size_t,compressed_matrix<double> > cMap;
  cMap& cm = theCharges.insert(make_pair(legs,cMap())).first->second;

  typedef map<size_t,vector<pair<size_t,size_t> > > ccMap;
  ccMap& ccm = theChargeNonZeros.insert(make_pair(legs,ccMap())).first->second;

  tmp.resize(dimPlus,dim);
  for ( size_t i = 0; i < legs.size(); ++i ) {
    size_t nonZero = 0;
    vector<pair<size_t,size_t> > nonZeros;
    for ( size_t a = 0; a < dimPlus; ++a )
      for ( size_t b = 0; b < dim; ++b ) {
	tmp(a,b) = tMatrixElement(i,a,b,legsPlus,legs);
	if ( tmp(a,b) != 0. ) {
	  ++nonZero;
	  nonZeros.push_back(make_pair(a,b));
	}
      }
    ccm.insert(make_pair(i,nonZeros));
    compressed_matrix<double>& tm = 
      cm.insert(make_pair(i,compressed_matrix<double>(dimPlus,dim,nonZero))).first->second;
    for ( size_t a = 0; a < dimPlus; ++a )
      for ( size_t b = 0; b < dim; ++b ) {
	if ( tmp(a,b) != 0. )
	  tm(a,b) = tmp(a,b);
      }
  }

  map<pair<size_t,size_t>,symmetric_matrix<double,upper> >& xm = theCorrelators[legs];
  for ( size_t i = 0; i < legs.size(); ++i )
    for ( size_t j = i+1; j < legs.size(); ++j ) {
      symmetric_matrix<double,upper>& mm =
	xm.insert(make_pair(make_pair(i,j),symmetric_matrix<double,upper>(dim,dim))).first->second;
      chargeProduct(cm[i],ccm[i],spPlus,cm[j],ccm[j],mm);
    }

  return dim;

}

void ColourBasis::chargeProduct(const compressed_matrix<double>& ti,
				const vector<pair<size_t,size_t> >& tiNonZero,
				const symmetric_matrix<double,upper>& X,
				const compressed_matrix<double>& tj,
				const vector<pair<size_t,size_t> >& tjNonZero,
				symmetric_matrix<double,upper>& result) const {
  for ( size_t i = 0; i < result.size1(); ++i )
    for ( size_t j = i; j < result.size1(); ++j )
      result(i,j) = 0.;
  for ( vector<pair<size_t,size_t> >::const_iterator i = tiNonZero.begin();
	i != tiNonZero.end(); ++i )
    for ( vector<pair<size_t,size_t> >::const_iterator j = tjNonZero.begin();
	  j != tjNonZero.end(); ++j ) {
      if ( j->second < i->second )
	continue;
      result(i->second,j->second) += 
	ti(i->first,i->second)*tj(j->first,j->second)*X(i->first,j->first);
    }
}

void ColourBasis::chargeProductAdd(const compressed_matrix<double>& ti,
				   const vector<pair<size_t,size_t> >& tiNonZero,
				   const matrix<Complex>& X,
				   const compressed_matrix<double>& tj,
				   const vector<pair<size_t,size_t> >& tjNonZero,
				   matrix<Complex>& result,
				   double factor) const {
  for ( vector<pair<size_t,size_t> >::const_iterator i = tiNonZero.begin();
	i != tiNonZero.end(); ++i )
    for ( vector<pair<size_t,size_t> >::const_iterator j = tjNonZero.begin();
	  j != tjNonZero.end(); ++j ) {
      result(i->first,j->first) += factor*
	ti(i->first,i->second)*tj(j->first,j->second)*X(i->second,j->second);
    }
}

string ColourBasis::cfstring(const list<list<pair<int,bool> > >& flow) {
  ostringstream out("");
  for ( list<list<pair<int,bool> > >::const_iterator line =
	  flow.begin(); line != flow.end(); ++line ) {
    for ( list<pair<int,bool> >::const_iterator node = 
	    line->begin(); node != line->end(); ++node ) {
      out << (node->second ? "-" : "") << (node->first+1) << " ";
    }
    if ( line != --(flow.end()) )
      out << ", ";
  }
  return out.str();
}

vector<string> ColourBasis::makeFlows(Ptr<Tree2toNDiagram>::tcptr diag,
				      size_t dim) const {

  vector<string> res(dim);

  list<list<list<pair<int,bool> > > > fdata =
    colourFlows(diag);

  cPDVector ext;
  tcPDVector dext = diag->external();
  copy(dext.begin(),dext.end(),back_inserter(ext));

  vector<PDT::Colour> colouredLegs =
    normalOrder(projectColour(ext));

  for ( list<list<list<pair<int,bool> > > >::const_iterator flow =
	  fdata.begin(); flow != fdata.end(); ++flow ) {
    for ( size_t i = 0; i < dim; ++i ) {
      bool matches = true;
      for ( list<list<pair<int,bool> > >::const_iterator line =
	      flow->begin(); line != flow->end(); ++line ) {
	pair<int,bool> front(diag->externalId(line->front().first),line->front().second);
	if ( front.first < 2  )
	  front.second = !front.second;
	pair<int,bool> back(diag->externalId(line->back().first),line->back().second);
	if ( back.first < 2 )
	  back.second = !back.second;
	if ( !colourConnected(ext,colouredLegs,front,back,i) ) {
	  matches = false;
	  break;
	}
      }
      if ( matches ) {
	assert(res[i] == "" && 
	       "only support colour bases with unique mapping to large-N colour flows");
	res[i] = cfstring(*flow);
      }
    }
  }

  bool gotone = false;
  for ( vector<string>::const_iterator f = res.begin();
	f != res.end(); ++f ) {
    if ( *f != "" ) {
      gotone = true;
      break;
    }
  }
  if ( !gotone ) {
    generator()->log() << "warning no color flow found for diagram\n";
    DiagramDrawer::drawDiag(generator()->log(),*diag);
  }

  return res;
  
}

size_t ColourBasis::prepare(const MEBase::DiagramVector& diags,
			    bool noCorrelations) {

  size_t dim = 0;

  for ( MEBase::DiagramVector::const_iterator d = diags.begin();
	d != diags.end(); ++d ) {
    Ptr<Tree2toNDiagram>::tcptr dd = dynamic_ptr_cast<Ptr<Tree2toNDiagram>::ptr>(*d);
    assert(dd);
    dim = prepare(dd->partons(),noCorrelations);
    if ( !haveColourFlows() || theFlowMap.find(dd) != theFlowMap.end() )
      continue;
    theFlowMap[dd] = makeFlows(dd,dim);
  }

  return dim;

}

bool matchEnd(int a, pair<int,bool> b,
	      Ptr<Tree2toNDiagram>::tcptr diag) {

  if ( a != b.first )
    return false;

  if ( b.first != diag->nSpace()-1 ) {
    return
      !b.second ? 
      diag->allPartons()[b.first]->hasColour() :
      diag->allPartons()[b.first]->hasAntiColour();
  } else {
    return
      !b.second ? 
      diag->allPartons()[b.first]->hasAntiColour() :
      diag->allPartons()[b.first]->hasColour();
  }

  return false;

}

bool findPath(pair<int,bool> a, pair<int,bool> b,
	      Ptr<Tree2toNDiagram>::tcptr diag,
	      list<pair<int,bool> >& path,
	      bool backward) {

  assert(a.first==0 ? !backward : true);

  if ( path.empty() )
    path.push_back(a);

  if ( !backward ) {

    if ( diag->children(a.first).first == -1 )
      return matchEnd(a.first,b,diag);

    pair<int,int> children = diag->children(a.first);

    bool cc = (children.first == diag->nSpace()-1);
    if ( diag->allPartons()[children.first]->coloured() )
      if ( !cc ? 
	   (!a.second ?
	    diag->allPartons()[children.first]->hasColour() :
	    diag->allPartons()[children.first]->hasAntiColour()) :
	   (!a.second ?
	    diag->allPartons()[children.first]->hasAntiColour() :
	    diag->allPartons()[children.first]->hasColour())  ) {
	pair<int,bool> next(children.first,a.second);
	path.push_back(next);
	if ( !findPath(next,b,diag,path,false) ) {
	  path.pop_back();
	} else return true;
      }

    cc = (children.second == diag->nSpace()-1);
    if ( diag->allPartons()[children.second]->coloured() )
      if ( !cc ? 
	   (!a.second ?
	    diag->allPartons()[children.second]->hasColour() :
	    diag->allPartons()[children.second]->hasAntiColour()) :
	   (!a.second ?
	    diag->allPartons()[children.second]->hasAntiColour() :
	    diag->allPartons()[children.second]->hasColour())  ) {
	pair<int,bool> next(children.second,a.second);
	path.push_back(next);
	if ( !findPath(next,b,diag,path,false) ) {
	  path.pop_back();
	} else return true;
      }

    if ( path.size() == 1 )
      path.pop_back();
    return false;

  } else {

    int parent = diag->parent(a.first);
    pair<int,int> neighbours = diag->children(parent);
    int neighbour = a.first == neighbours.first ? neighbours.second : neighbours.first;

    if ( matchEnd(parent,b,diag) ) {
      path.push_back(b);
      return true;
    }

    if ( matchEnd(neighbour,b,diag) ) {
      path.push_back(b);
      return true;
    }

    if ( diag->allPartons()[neighbour]->coloured() ) 
      if ( a.second ?
	   diag->allPartons()[neighbour]->hasColour() :
	   diag->allPartons()[neighbour]->hasAntiColour() ) {
	pair<int,bool> next(neighbour,!a.second);
	path.push_back(next);
	if ( !findPath(next,b,diag,path,false) ) {
	  path.pop_back();
	} else return true;
      }

    if ( parent == 0 ) {
      if ( path.size() == 1 )
	path.pop_back();
      return false;
    }

    if ( diag->allPartons()[parent]->coloured() ) 
      if ( !a.second ?
	   diag->allPartons()[parent]->hasColour() :
	   diag->allPartons()[parent]->hasAntiColour() ) {
	pair<int,bool> next(parent,a.second);
	path.push_back(next);
	if ( !findPath(next,b,diag,path,true) ) {
	  path.pop_back();
	} else return true;
      }

    if ( path.size() == 1 )
      path.pop_back();
    return false;

  }

  return false;

}


list<pair<int,bool> > ColourBasis::colouredPath(pair<int,bool> a, pair<int,bool> b,
						Ptr<Tree2toNDiagram>::tcptr diag) {

  list<pair<int,bool> > res;

  if ( a.first == b.first )
    return res;

  bool aIn = (a.first < 2);
  bool bIn = (b.first < 2);

  if ( (aIn && bIn) || (!aIn && !bIn) )
    if ( (a.second && b.second) ||
	 (!a.second && !b.second) )
      return res;

  if ( (aIn && !bIn) || (!aIn && bIn) )
    if ( (!a.second && b.second) ||
	 (a.second && !b.second) )
      return res;

  if ( a.first > b.first )
    swap(a,b);

  a.first = diag->diagramId(a.first);
  b.first = diag->diagramId(b.first);

  if ( a.first == diag->nSpace()-1 )
    a.second = !a.second;

  if ( b.first == diag->nSpace()-1 )
    b.second = !b.second;

  if ( !findPath(a,b,diag,res,a.first != 0) )
    return res;

  if ( b.first == diag->nSpace()-1 ) {
    res.back().second = !res.back().second;
  }

  if ( a.first == diag->nSpace()-1 ) {
    res.front().second = !res.front().second;
  }

  return res;

}

list<list<list<pair<int,bool> > > >
ColourBasis::colourFlows(Ptr<Tree2toNDiagram>::tcptr diag) {

  vector<pair<int,bool> > connectSource;
  vector<pair<int,bool> > connectSink;
  for ( size_t i = 0; i != diag->partons().size(); ++i ) {
    if ( i < 2 && diag->partons()[i]->hasAntiColour() )
      connectSource.push_back(make_pair(i,true));
    if ( i < 2 && diag->partons()[i]->hasColour() )
      connectSink.push_back(make_pair(i,false));
    if ( i > 1 && diag->partons()[i]->hasColour() )
      connectSource.push_back(make_pair(i,false));
    if ( i > 1 && diag->partons()[i]->hasAntiColour() )
      connectSink.push_back(make_pair(i,true));
  }

  assert(connectSource.size() == connectSink.size());

  list<list<list<pair<int,bool> > > > ret;

  do {

    vector<pair<int,bool> >::iterator source =
      connectSource.begin();
    vector<pair<int,bool> >::iterator sink =
      connectSink.begin();
    list<list<pair<int,bool> > > res;
    for ( ; source != connectSource.end(); ++source, ++sink ) {
      if ( source->first == sink->first ) {
	res.clear();
	break;
      }
      list<pair<int,bool> > line =
	colouredPath(*source,*sink,diag);
      if ( line.empty() ) {
	res.clear();
	break;
      }
      res.push_back(line);
    }

    if ( !res.empty() ) {

      // check, if all dressed properly
      vector<pair<int,int> > dressed((*diag).allPartons().size(),make_pair(0,0));
      for ( size_t p = 0; p < diag->allPartons().size(); ++p ) {
	if ( diag->allPartons()[p]->hasColour() &&
	     !diag->allPartons()[p]->hasAntiColour() )
	  dressed[p].first = 1;
	if ( diag->allPartons()[p]->hasAntiColour() &&
	     !diag->allPartons()[p]->hasColour() )
	  dressed[p].second = 1;
	if ( diag->allPartons()[p]->hasAntiColour() &&
	     diag->allPartons()[p]->hasColour() ) {
	  dressed[p].first = 1; dressed[p].second = 1;
	}
      }
      for ( list<list<pair<int,bool> > >::const_iterator l = res.begin();
	    l != res.end(); ++l ) {
	for ( list<pair<int,bool> >::const_iterator n = l->begin();
	      n != l->end(); ++n ) {
	  if ( !(n->second) )
	    dressed[n->first].first -= 1;
	  else
	    dressed[n->first].second -= 1;
	}
      }
      for ( vector<pair<int,int> >::const_iterator d = dressed.begin();
	    d != dressed.end(); ++d ) {
	if ( d->first != 0 || d->second != 0 ) {
	  res.clear();
	  break;
	}
      }

      if ( !res.empty() )
	ret.push_back(res);

    }

  } while ( std::next_permutation(connectSink.begin(),connectSink.end()) );

  return ret;

}

void ColourBasis::updateColourLines(Ptr<Tree2toNDiagram>::tcptr dd) {
  map<Ptr<Tree2toNDiagram>::tcptr,vector<string> >::const_iterator cl =
    theFlowMap.find(dd);
  assert(cl != theFlowMap.end());
  vector<ColourLines*> clines(cl->second.size());
  for ( size_t k = 0; k < cl->second.size(); ++k ) {
    if ( cl->second[k] == "" ) {
      clines[k] = 0;
      continue;
    }
    clines[k] = new ColourLines(cl->second[k]);
  }
  theColourLineMap[cl->first] = clines;
}

map<Ptr<Tree2toNDiagram>::tcptr,vector<ColourLines*> >&
ColourBasis::colourLineMap() {
  
  if ( !theColourLineMap.empty() )
    return theColourLineMap;

  for ( map<Ptr<Tree2toNDiagram>::tcptr,vector<string> >::const_iterator cl =
	  theFlowMap.begin(); cl != theFlowMap.end(); ++cl ) {
    vector<ColourLines*> clines(cl->second.size());
    for ( size_t k = 0; k < cl->second.size(); ++k ) {
      if ( cl->second[k] == "" ) {
	clines[k] = 0;
	continue;
      }
      clines[k] = new ColourLines(cl->second[k]);
    }
    theColourLineMap[cl->first] = clines;
  }

  return theColourLineMap;

}

Selector<const ColourLines *> ColourBasis::colourGeometries(tcDiagPtr diag,
							    const map<vector<int>,CVector>& amps) {
  Ptr<Tree2toNDiagram>::tcptr dd = 
    dynamic_ptr_cast<Ptr<Tree2toNDiagram>::tcptr>(diag);
  assert(dd && theFlowMap.find(dd) != theFlowMap.end());
  map<Ptr<Tree2toNDiagram>::tcptr,vector<ColourLines*> >::const_iterator colit =
    colourLineMap().find(dd);
  if ( colit == colourLineMap().end() ) {
    updateColourLines(dd);
    colit = colourLineMap().find(dd);
  }
  const vector<ColourLines*>& cl = colit->second;

  Selector<const ColourLines *> sel;
  size_t dim = amps.begin()->second.size();
  assert(dim == cl.size());
  double w = 0.;
  for ( size_t i = 0; i < dim; ++i ) {
    if ( !cl[i] )
      continue;
    w = 0.;
    for ( map<vector<int>,CVector>::const_iterator a = amps.begin();
	  a != amps.end(); ++a )
      w += real(conj((a->second)(i))*((a->second)(i)));
    if ( w > 0. )
      sel.insert(w,cl[i]);
  }
  assert(!sel.empty());
  return sel;
}

size_t ColourBasis::tensorIdFromFlow(tcDiagPtr diag, const ColourLines * flow) {

 Ptr<Tree2toNDiagram>::tcptr dd = 
    dynamic_ptr_cast<Ptr<Tree2toNDiagram>::tcptr>(diag);
  assert(dd && theFlowMap.find(dd) != theFlowMap.end());
  map<Ptr<Tree2toNDiagram>::tcptr,vector<ColourLines*> >::const_iterator colit =
    colourLineMap().find(dd);
  if ( colit == colourLineMap().end() ) {
    updateColourLines(dd);
    colit = colourLineMap().find(dd);
  }

  const vector<ColourLines*>& cl = colit->second;

  size_t res = 0;
  for ( ; res < cl.size(); ++res ) {
    if ( flow == cl[res] )
      break;
  }

  assert(res < cl.size());

  return res;

}

const symmetric_matrix<double,upper>& ColourBasis::scalarProducts(const cPDVector& sub) const {

  map<cPDVector,vector<PDT::Colour> >::const_iterator lit =
    theNormalOrderedLegs.find(sub);
  assert(lit != theNormalOrderedLegs.end());

  ScalarProductMap::const_iterator spit =
    theScalarProducts.find(lit->second);
  assert(spit != theScalarProducts.end());

  return spit->second;

}

const compressed_matrix<double>& ColourBasis::charge(const cPDVector& sub, size_t iIn) const {

  map<cPDVector,vector<PDT::Colour> >::const_iterator lit =
    theNormalOrderedLegs.find(sub);
  assert(lit != theNormalOrderedLegs.end());

  ChargeMap::const_iterator ct =
    theCharges.find(lit->second);
  assert(ct != theCharges.end());

  map<cPDVector,map<size_t,size_t> >::const_iterator trans
    = theIndexMap.find(sub);
  assert(trans != theIndexMap.end());
  size_t i = trans->second.find(iIn)->second;

  map<size_t,compressed_matrix<double> >::const_iterator cit
    = ct->second.find(i);
  assert(cit != ct->second.end());

  return cit->second;

}

const vector<pair<size_t,size_t> >& ColourBasis::chargeNonZero(const cPDVector& sub, size_t iIn) const {

  map<cPDVector,vector<PDT::Colour> >::const_iterator lit =
    theNormalOrderedLegs.find(sub);
  assert(lit != theNormalOrderedLegs.end());

  ChargeNonZeroMap::const_iterator ct =
    theChargeNonZeros.find(lit->second);
  assert(ct != theChargeNonZeros.end());

  map<cPDVector,map<size_t,size_t> >::const_iterator trans
    = theIndexMap.find(sub);
  assert(trans != theIndexMap.end());
  size_t i = trans->second.find(iIn)->second;

  map<size_t,vector<pair<size_t,size_t> > >::const_iterator cit
    = ct->second.find(i);
  assert(cit != ct->second.end());

  return cit->second;

}

const symmetric_matrix<double,upper>& ColourBasis::correlator(const cPDVector& sub,
							      const pair<size_t,size_t>& ijIn) const {

  map<cPDVector,vector<PDT::Colour> >::const_iterator lit =
    theNormalOrderedLegs.find(sub);
  assert(lit != theNormalOrderedLegs.end());

  CorrelatorMap::const_iterator cit =
    theCorrelators.find(lit->second);
  assert(cit != theCorrelators.end());

  map<cPDVector,map<size_t,size_t> >::const_iterator trans
    = theIndexMap.find(sub);
  assert(trans != theIndexMap.end());
  pair<size_t,size_t> ij(trans->second.find(ijIn.first)->second,
			 trans->second.find(ijIn.second)->second);
  if ( ij.first > ij.second )
    swap(ij.first,ij.second);

  map<pair<size_t,size_t>,symmetric_matrix<double,upper> >::const_iterator cijit
    = cit->second.find(ij);
  assert(cijit != cit->second.end());

  return cijit->second;

}

double ColourBasis::me2(const cPDVector& sub, 
			const map<vector<int>,CVector>& amps) const {

  const symmetric_matrix<double,upper>& sp = scalarProducts(sub);

  double res = 0.;

  for ( map<vector<int>,CVector>::const_iterator a = amps.begin();
	a != amps.end(); ++a ) {
    res += real(inner_prod(boost::numeric::ublas::conj(a->second),prod(sp,a->second)));
  }

  return res;

}

double ColourBasis::interference(const cPDVector& sub, 
				 const map<vector<int>,CVector>& amps1,
				 const map<vector<int>,CVector>& amps2) const {

  const symmetric_matrix<double,upper>& sp = scalarProducts(sub);

  double res = 0.;

  map<vector<int>,CVector>::const_iterator a = amps1.begin();
  map<vector<int>,CVector>::const_iterator b = amps2.begin();
  for ( ; a != amps1.end(); ++a, ++b ) {
    assert(a->first == b->first);
    res += 2.*real(inner_prod(boost::numeric::ublas::conj(a->second),prod(sp,b->second)));
  }

  assert(!isnan(res));

  return res;

}

double ColourBasis::colourCorrelatedME2(const pair<size_t,size_t>& ij,
					const cPDVector& sub, 
					const map<vector<int>,CVector>& amps) const {

  const symmetric_matrix<double,upper>& cij = correlator(sub,ij);

  double res = 0.;

  for ( map<vector<int>,CVector>::const_iterator a = amps.begin();
	a != amps.end(); ++a ) {
    res += real(inner_prod(boost::numeric::ublas::conj(a->second),prod(cij,a->second)));
  }

  return res;

}

Complex ColourBasis::interference(const cPDVector& sub, 
				  const CVector& left,
				  const CVector& right) const {

  const symmetric_matrix<double,upper>& sp = scalarProducts(sub);
  return inner_prod(boost::numeric::ublas::conj(left),prod(sp,right));

}

Complex ColourBasis::colourCorrelatedInterference(const pair<size_t,size_t>& ij,
						  const cPDVector& sub, 
						  const CVector& left,
						  const CVector& right) const {

  const symmetric_matrix<double,upper>& cij = correlator(sub,ij);
  return inner_prod(boost::numeric::ublas::conj(left),prod(cij,right));

}

double ColourBasis::me2(const cPDVector& sub, 
			const matrix<Complex>& amp) const {

  const symmetric_matrix<double,upper>& sp = scalarProducts(sub);

  double tr = 0;

  size_t n = amp.size1();

  for ( size_t i = 0; i < n; ++i ) {
    tr += real(inner_prod(row(sp,i),column(amp,i)));
  }

  return tr;

}

double ColourBasis::colourCorrelatedME2(const pair<size_t,size_t>& ij,
					const cPDVector& sub, 
					const matrix<Complex>& amp) const {

  const symmetric_matrix<double,upper>& cij = correlator(sub,ij);

  double tr = 0;

  size_t n = amp.size1();

  for ( size_t i = 0; i < n; ++i ) {
    tr += real(inner_prod(row(cij,i),column(amp,i)));
  }

  return tr;

}

struct pickColour {
  PDT::Colour operator()(tcPDPtr p) const {
    return p->iColour();
  }
};

vector<PDT::Colour> ColourBasis::projectColour(const cPDVector& sub) const {
  vector<PDT::Colour> res(sub.size());
  transform(sub.begin(),sub.end(),res.begin(),pickColour());
  return res;
}

vector<PDT::Colour> ColourBasis::normalOrder(const vector<PDT::Colour>& legs) const {
  vector<PDT::Colour> crosslegs = legs;
  if ( crosslegs[0] == PDT::Colour3 )
    crosslegs[0] = PDT::Colour3bar;
  else if ( crosslegs[0] == PDT::Colour3bar )
    crosslegs[0] = PDT::Colour3;
  if ( crosslegs[1] == PDT::Colour3 )
    crosslegs[1] = PDT::Colour3bar;
  else if ( crosslegs[1] == PDT::Colour3bar )
    crosslegs[1] = PDT::Colour3;
  int n3 = count_if(crosslegs.begin(),crosslegs.end(),matchRep(PDT::Colour3));
  int n8 = count_if(crosslegs.begin(),crosslegs.end(),matchRep(PDT::Colour8));
  vector<PDT::Colour> ordered(2*n3+n8,PDT::Colour8);
  int i = 0;
  while ( i < 2*n3 ) {
    ordered[i] = PDT::Colour3;
    ordered[i+1] = PDT::Colour3bar;
    i+=2;
  }
  return ordered;
}

string ColourBasis::file(const vector<PDT::Colour>& sub) const {

  string res = name() + "-";

  for ( vector<PDT::Colour>::const_iterator lit = sub.begin();
	lit != sub.end(); ++lit ) {
    if ( *lit == PDT::Colour3 )
      res += "3";
    if ( *lit == PDT::Colour3bar )
      res += "3bar";
    if ( *lit == PDT::Colour8 )
      res += "8";
  }

  if ( largeN() )
    res += "largeN";

  return res;

}

void ColourBasis::writeBasis(const string& prefix) const {

  if ( didWrite )
    return;

  set<vector<PDT::Colour> > legs;
  for ( map<cPDVector,vector<PDT::Colour> >::const_iterator lit
	  = theNormalOrderedLegs.begin(); lit != theNormalOrderedLegs.end(); ++lit ) {
    legs.insert(lit->second);
  }

  string searchPath = theSearchPath;

  if ( searchPath != "" )
    if ( *(--searchPath.end()) != '/' )
      searchPath += "/";

  for ( set<vector<PDT::Colour> >::const_iterator known = legs.begin();
	known != legs.end(); ++known ) {
    string fname = searchPath + prefix + file(*known) + ".cdat";
    if ( !( SamplerBase::runLevel() == SamplerBase::ReadMode ||
            SamplerBase::runLevel() == SamplerBase::BuildMode ) ) {
      ifstream check(fname.c_str());
      if ( check ) continue;
    }
    ofstream out(fname.c_str());
    if ( !out )
      throw Exception() << "ColourBasis: Failed to open "
			<< fname << " for storing colour basis information."
			<< Exception::runerror;
    out << setprecision(18);
    const symmetric_matrix<double,upper>& sp = 
      theScalarProducts.find(*known)->second;
    write(sp,out);
    if ( theCharges.find(*known) != theCharges.end() ) {
      out << "#charges\n";
      const map<size_t,compressed_matrix<double> >& tm =
	theCharges.find(*known)->second;
      const map<size_t,vector<pair<size_t,size_t> > >& tc =
	theChargeNonZeros.find(*known)->second;
      map<size_t,vector<pair<size_t,size_t> > >::const_iterator kc =
	tc.begin();
      for ( map<size_t,compressed_matrix<double> >::const_iterator k = tm.begin();
	    k != tm.end(); ++k, ++kc ) {
	out << k->first << "\n";
	write(k->second,out,kc->second);
      }
      const map<pair<size_t,size_t>,symmetric_matrix<double,upper> >& cm =
	theCorrelators.find(*known)->second;
      for ( map<pair<size_t,size_t>,symmetric_matrix<double,upper> >::const_iterator k =
	      cm.begin(); k != cm.end(); ++k ) {
	out << k->first.first << "\n" << k->first.second << "\n";
	write(k->second,out);
      }
    } else {
      out << "#nocharges\n";
    }
    out << flush;
  }

  didWrite = true;

}

bool ColourBasis::readBasis(const vector<PDT::Colour>& legs) {

  string searchPath = theSearchPath;

  if ( searchPath != "" )
    if ( *(--searchPath.end()) != '/' )
      searchPath += "/";

  string fname = searchPath + file(legs) + ".cdat";
  ifstream in(fname.c_str());
  if ( !in )
    return false;
  read(theScalarProducts[legs],in);
  string tag; in >> tag;
  if ( tag != "#nocharges" ) {
    for ( size_t k = 0; k < legs.size(); ++k ) {
      size_t i; in >> i;
      read(theCharges[legs][i],in,theChargeNonZeros[legs][i]);
    }
    for ( size_t k = 0; k < legs.size()*(legs.size()-1)/2; ++k ) {
      size_t i,j; in >> i >> j;
      read(theCorrelators[legs][make_pair(i,j)],in);
    }
  }

  readBasisDetails(legs);

  return true;

}

void ColourBasis::readBasis() {

  if ( didRead )
    return;

  string searchPath = theSearchPath;

  if ( searchPath != "" )
    if ( *(--searchPath.end()) != '/' )
      searchPath += "/";

  set<vector<PDT::Colour> > legs;
  for ( map<cPDVector,vector<PDT::Colour> >::const_iterator lit
	  = theNormalOrderedLegs.begin(); lit != theNormalOrderedLegs.end(); ++lit )
    legs.insert(lit->second);

  for ( set<vector<PDT::Colour> >::const_iterator known = legs.begin();
	known != legs.end(); ++known ) {
    if ( theScalarProducts.find(*known) != theScalarProducts.end() )
      continue;
    string fname = searchPath + file(*known) + ".cdat";
    if ( !readBasis(*known) )
      throw Exception() << "ColourBasis: Failed to open "
			<< fname << " for reading colour basis information."
			<< Exception::runerror;
  }

  didRead = true;

}

void ColourBasis::write(const symmetric_matrix<double,upper>& m, ostream& os) const {
  os << m.size1() << "\n";
  for ( size_t i = 0; i < m.size1(); ++i )
    for ( size_t j = i; j < m.size1(); ++j )
      os << m(i,j) << "\n";
  os << flush;
}

void ColourBasis::read(symmetric_matrix<double,upper>& m, istream& is) {
  size_t s; is >> s;
  m.resize(s);
  for ( size_t i = 0; i < m.size1(); ++i )
    for ( size_t j = i; j < m.size1(); ++j )
      is >> m(i,j);
}

void ColourBasis::write(const compressed_matrix<double>& m, ostream& os,
			const vector<pair<size_t,size_t> >& nonZeros) const {
  os << nonZeros.size() << "\n"
     << m.size1() << "\n"
     << m.size2() << "\n";
  for ( vector<pair<size_t,size_t> >::const_iterator nz = nonZeros.begin();
	nz != nonZeros.end(); ++nz )
    os << nz->first << "\n" << nz->second << "\n"
       << m(nz->first,nz->second) << "\n";
  os << flush;
}

void ColourBasis::read(compressed_matrix<double>& m, istream& is,
		       vector<pair<size_t,size_t> >& nonZeros) {
  size_t nonZero, size1, size2; 
  is >> nonZero >> size1 >> size2;
  nonZeros.resize(nonZero);
  m = compressed_matrix<double>(size1,size2,nonZero);
  for ( size_t k = 0; k < nonZero; ++k ) {
    size_t i,j; double val;
    is >> i >> j >> val;
    nonZeros[k] = make_pair(i,j);
    m(i,j) = val;
  }
}

void ColourBasis::doinit() {
  HandlerBase::doinit();
  if ( theSearchPath.empty() && factory() )
    theSearchPath = factory()->buildStorage();
  readBasis();
}

void ColourBasis::dofinish() {
  HandlerBase::dofinish();
  writeBasis();
}

void ColourBasis::doinitrun() {
  HandlerBase::doinitrun();
  if ( theSearchPath.empty() && factory() )
    theSearchPath = factory()->buildStorage();
  readBasis();
}

void ColourBasis::persistentOutput(PersistentOStream & os) const {
  os << theLargeN << theNormalOrderedLegs
     << theIndexMap << theFlowMap << theOrderingStringIdentifiers 
     << theOrderingIdentifiers << theFactory << theSearchPath;
  writeBasis();
}

void ColourBasis::persistentInput(PersistentIStream & is, int) {
  is >> theLargeN >> theNormalOrderedLegs
     >> theIndexMap >> theFlowMap >> theOrderingStringIdentifiers
     >> theOrderingIdentifiers >> theFactory >> theSearchPath;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeAbstractClass<ColourBasis,HandlerBase>
describeColourBasis("Herwig::ColourBasis", "Herwig.so");

void ColourBasis::Init() {

  static ClassDocumentation<ColourBasis> documentation
    ("ColourBasis is an interface to a colour basis "
     "implementation.");

  static Switch<ColourBasis,bool> interfaceLargeN
    ("LargeN",
     "Switch on or off large-N evaluation.",
     &ColourBasis::theLargeN, false, false, false);
  static SwitchOption interfaceLargeNOn
    (interfaceLargeN,
     "On",
     "Work in N=infinity",
     true);
  static SwitchOption interfaceLargeNOff
    (interfaceLargeN,
     "Off",
     "Work in N=3",
     false);

}

