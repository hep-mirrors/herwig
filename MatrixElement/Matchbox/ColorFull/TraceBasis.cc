// -*- C++ -*-
//
// TraceBasis.cc is a part of ColorFull
// Copyright (C) 2010-2011 Simon Platzer & Malin Sjodahl
//
// ColorFull is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TraceBasis class.
//

#include "TraceBasis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ColorFull;

TraceBasis::TraceBasis() {}

TraceBasis::~TraceBasis() {}

IBPtr TraceBasis::clone() const {
  return new_ptr(*this);
}

IBPtr TraceBasis::fullclone() const {
  return new_ptr(*this);
}

void TraceBasis::clear() {
  ColourBasis::clear();
  theBasisMap.clear();
  theScalarProducts.clear();
}

map<size_t,vector<vector<size_t> > > 
TraceBasis::basisList(const vector<PDT::Colour>& basisId) const {

  map<vector<PDT::Colour>,Trace_basis>::const_iterator bit =
    theBasisMap.find(basisId);

  map<size_t,vector<vector<size_t> > > res;

  const col_basis& cb = bit->second.cb;

  for ( size_t i = 0; i < cb.size(); ++i ) {

    vector<vector<size_t> > cstr;

    const Col_str& cs = cb.at(i).at(0);

    for ( size_t j = 0; j < cs.size(); ++j ) {

      const Quark_line& ql = cs.at(j);
      vector<size_t> qline;

      for ( size_t k = 0; k < ql.size(); ++k )
	qline.push_back(ql.at(k)-1);

      cstr.push_back(qline);

    }

    res[i] = cstr;

  }

  return res;

}

size_t TraceBasis::prepareBasis(const vector<PDT::Colour>& sub) {

  useMe();

  vector<PDT::Colour> mySub = normalOrder(sub);

  if ( theBasisMap.find(mySub) == theBasisMap.end() ) {

    int ng = count_if(mySub.begin(),mySub.end(),ColourBasis::matchRep(PDT::Colour8));
    int nq = (mySub.size() - ng)/2;

    Trace_basis basis;
    basis.create_basis(nq,ng);
    theBasisMap[mySub] = basis;

  }

  return theBasisMap[mySub].size();

}

void TraceBasis::readBasisDetails(const vector<PDT::Colour>& sub) {
  prepareBasis(sub);
}

double TraceBasis::scalarProduct(size_t i, size_t j,
				 const vector<PDT::Colour>& abBasis) const {

  if ( largeN() && i != j )
    return 0.;

  map<vector<PDT::Colour>,Trace_basis>::const_iterator bit =
    theBasisMap.find(abBasis);

  assert(bit != theBasisMap.end());

  const Trace_basis& Basis = bit->second;
  Col_str csi = Basis.cb.at(i).at(0);
  Col_str csj = Basis.cb.at(j).at(0);

  // Rename indices, and make string of new Col_strs, to use in map
  //pair<Col_str, Col_str> Css_new = Basis.rename_indices(csi,csj);
  Basis.rename_indices(csi,csj);
  Col_str Cs1 = csi;
  Col_str Cs2 = csj;
  ostringstream Cs_string;
  Cs_string << Cs1 << Cs2;

  map<string,Polynomial>::iterator pit = theScalarProducts.find(Cs_string.str());

  // Add result to map
  if ( pit == theScalarProducts.end() ) {
    Polynomial p = colorFunctions.scalar_product(Cs1, Cs2);
    theScalarProducts.insert(make_pair(Cs_string.str(), p));
    pit = theScalarProducts.find(Cs_string.str());
  }

  return
    largeN() ? 
    colorFunctions.double_num(colorFunctions.leading(pit->second)) :
    colorFunctions.double_num(pit->second);

}

double TraceBasis::tMatrixElement(size_t m, size_t i, size_t j,
				  const vector<PDT::Colour>& aBasis,
				  const vector<PDT::Colour>& bBasis) const {

  ++m;

  map<vector<PDT::Colour>,Trace_basis>::iterator ait =
    theBasisMap.find(aBasis);


  map<vector<PDT::Colour>,Trace_basis>::const_iterator bit =
    theBasisMap.find(bBasis);

  assert(bit != theBasisMap.end());
  assert(ait != theBasisMap.end());


  Trace_basis& ABasis = ait->second;
  const Trace_basis& BBasis = bit->second;

  pair<int,int> newNumbers = ABasis.new_vector_numbers(BBasis.cb.at(j).at(0),m);

  if ( (size_t) newNumbers.first == i )
    return 1.;

  if ( (size_t) newNumbers.second == i )
    return -1.;

  return 0.;

}

bool TraceBasis::colourConnected(const cPDVector& sub,
				 const vector<PDT::Colour>& basisId,
				 const pair<int,bool>& first,
				 const pair<int,bool>& second, 
				 size_t tensor) const {

  // get the basis
  map<vector<PDT::Colour>,Trace_basis>::const_iterator bit =
    theBasisMap.find(basisId);
  assert(bit != theBasisMap.end());
  const Trace_basis& basis = bit->second;

  // translate process to basis ids
  map<cPDVector,map<size_t,size_t> >::const_iterator trans
    = indexMap().find(sub);
  assert(trans != indexMap().end());

  int idColoured = first.second ? second.first : first.first;
  idColoured = trans->second.find(idColoured)->second;
  ++idColoured;
  int idAntiColoured = first.second ? first.first : second.first;
  idAntiColoured = trans->second.find(idAntiColoured)->second;
  ++idAntiColoured;

  const Col_str& cs = basis.cb.at(tensor).at(0);

  return cs.left_neighbor(idAntiColoured,idColoured);

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void TraceBasis::persistentOutput(PersistentOStream &) const {}

void TraceBasis::persistentInput(PersistentIStream & , int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<TraceBasis,Herwig::ColourBasis>
  describeTraceBasis("ColorFull::TraceBasis", 
		     "HwColorFull.so");

void TraceBasis::Init() {

  static ClassDocumentation<TraceBasis> documentation
    ("TraceBasis implements the trace colour basis.",
     "The colour algebra has been performed using ColorFull \\cite{Sjodahl:2014opa}",
     "%\\cite{Sjodahl:2014opa}\n"
     "\\bibitem{Sjodahl:2014opa}\n"
     "M.~Sjodahl,\n"
     "``ColorFull -- a C++ library for calculations in SU(Nc)color space,''\n"
     "arXiv:1412.3967 [hep-ph].\n"
     "%%CITATION = ARXIV:1412.3967;%%");

}

