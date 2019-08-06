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

TraceBasis::TraceBasis() { }

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


// m is emitting parton
// i is new vector number in new large aBasis
// j is old vector number in old small bBasis
// k is the new number of the emitter
// l is the number of the new gluon
// dict contains the map from old to new numbers of the partons not participating
double TraceBasis::tMatrixElement(size_t m, size_t i, size_t j,
		const vector<PDT::Colour>& aBasis,
		const vector<PDT::Colour>& bBasis,
		size_t k, size_t l,
		const map<size_t,size_t>& dict) const {

  // Call sMatrixElement if it's a gluon splitting
  if ( bBasis[m] == PDT::Colour8 && aBasis[k] != PDT::Colour8 ) 
    return sMatrixElement(m,i,j,aBasis,bBasis,k,l,dict);

	assert( dict.size()+1 == bBasis.size() );

	map<vector<PDT::Colour>,Trace_basis>::const_iterator ait =
			theBasisMap.find(aBasis);
	map<vector<PDT::Colour>,Trace_basis>::const_iterator bit =
			theBasisMap.find(bBasis);

	assert(bit != theBasisMap.end());
	assert(ait != theBasisMap.end());


	const Trace_basis& ABasis = ait->second;
	const Trace_basis& BBasis = bit->second;


	// Initial Col_amp, before emission
	Col_amp Ca1 = BBasis.at(j);

	int g_new, p_old;
	p_old = m+1;// ColorFull starts with parton number 1 (normally q or g)
	int g_old = aBasis.size();//Give the new g an index that isn't used 

	// If the emitted parton is not a gluon, this is because
	// we have backward evolution with a qqbar-> g
	if( aBasis[l] != PDT::Colour::Colour8 ){
		// The number of the new gluon is k+1
		g_new = k+1;
	} else{ // Standard case, l (+1) is the number of the new gluon
		g_new = l+1;
	}

	// Color structure after gluon split
	Col_amp Ca2 = colorFunctions.emit_gluon( Ca1, p_old, g_old);

	// The map should be the dict (containing the map of old non-involved partons)
	// + g_new mapped to itself and p_old mapped
	std::map<int, int> map;
	for ( size_t ii = 0; ii < bBasis.size(); ii++ ){
		if ( ii  != m ) {// exclude splitting gluon
			int parton = static_cast<int>( dict.at(ii) );
			map[ii+1] = parton+1;// Parton numbers one unit higher in ColorFull
		}
	}
	// New parton number mapping
	map[ g_old ] = g_new;
	// Old emitter should also be mapped
	if( aBasis[l] != PDT::Colour::Colour8 ){// Exceptional case
		map[ p_old ] = l+1;
	}
	else{ // Standard case
		map[ p_old ] = k+1;
	}
	// Check the size of the map
	assert( map.size() == aBasis.size() );

	// Color structure when partons have names in map
	Col_amp Ca3= colorFunctions.rename_partons( Ca2, map );

	// Check if the new color structure has a component for basis vector i in ABasis
	// Loop over Col_strs in Col_amp after split
	for ( uint Csi=0; Csi < Ca3.size(); Csi++ ){
		if( Ca3.at(Csi).cs == ABasis.at(i).at(0).cs ){// If col_strs are the same
			return colorFunctions.double_num( Ca3.at(Csi).Poly ); // Return Polynomial coefficient in front of Col_str
		};
	}

	// If no component corresponding to vector i is found
	return 0.;

}

// m is splitting gluon
// i is new vector number new large aBasis
// j is old vector number in old small bBasis
// k is the number of the new emitter after (q or qbar)
// l is the number of the "emission" (q or qbar)
// dict contains the map from old to new numbers of the partons not participating
double TraceBasis::sMatrixElement(size_t m, size_t i, size_t j,
		const vector<PDT::Colour>& aBasis,
		const vector<PDT::Colour>& bBasis,
		size_t k, size_t l,
		const map<size_t,size_t>& dict) const {

	// Check that dict has the right size (splitting gluon missing, hence +1)
	assert( dict.size()+1 == bBasis.size() );

	map<vector<PDT::Colour>,Trace_basis>::iterator ait =
			theBasisMap.find(aBasis);
	map<vector<PDT::Colour>,Trace_basis>::const_iterator bit =
			theBasisMap.find(bBasis);

	assert(bit != theBasisMap.end());
	assert(ait != theBasisMap.end());


	Trace_basis& ABasis = ait->second; // New basis
	const Trace_basis& BBasis = bit->second; // Old basis


	// Initial Col_amp, before split
	Col_amp Ca1 = BBasis.at(j);


	int g_old = m+1;// ColorFull starts with parton number 1 (normally q or g)
	int q_old, q_new, qbar_old, qbar_new;
	q_old = aBasis.size();// Give the q an index that isn't used
	qbar_old = aBasis.size()+1;// Give the qbar an index that isn't used
	if ( aBasis[k] == PDT::Colour::Colour3 && aBasis[l] == PDT::Colour::Colour3bar ) {
	  q_new = k+1;
	  qbar_new = l+1;
	} else if ( aBasis[l] == PDT::Colour::Colour3 && aBasis[k] == PDT::Colour::Colour3bar ) {
	  q_new = l+1;
	  qbar_new = k+1;
	} else {
	  assert(false);
	}

	// Color structure after gluon split
	// split_gluon also simplifies, so no polynomial factor in Qls
	Col_amp Ca2= colorFunctions.split_gluon( Ca1, g_old, q_old, qbar_old );

	// The map should be the dict (containing the map of old non-involved partons)
	std::map<int, int> map;
	for ( size_t ii = 0; ii < bBasis.size(); ii++ ){
		if ( ii  != m ) {// exclude splitting gluon
			int parton = static_cast<int>( dict.at(ii) );
			map[ii+1] = parton+1;// Parton numbers one unit higher in ColorFull
		}
	}

	map[q_old] = q_new;
	map[qbar_old] = qbar_new;

	assert( map.size() == aBasis.size() );

	// Color structure when partons have ColorFull default names
	Col_amp Ca3= colorFunctions.rename_partons( Ca2, map ); // does normal ordering as well

	// Check if the new color structure has a component for basis vector i in ABasis
	// Loop over Col_strs in Col_amp after split
	for ( uint Csi=0; Csi < Ca3.size(); Csi++ ){
		if( Ca3.at(Csi).cs == ABasis.at(i).at(0).cs ){// If col_strs are the same, note that they must be normal ordered
			return colorFunctions.double_num( Ca3.at(Csi).Poly ); // Return Polynomial coefficient in front of Col_str
		};
	}

	// If no component corresponding to vector i is found
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

map<size_t,size_t> TraceBasis::indexChange(const vector<PDT::Colour>& basis,
					   const size_t dim,
					   const map<size_t,size_t>& indPerm) const {
  // Change the map to indices starting at 1 (colorfull numbering of legs) from
  // starting at 0 (Herwig numbering of legs).
  map<int,int> iPerm;
  for ( map<size_t,size_t>::const_iterator it = indPerm.begin();
	it != indPerm.end(); it++ ) {
    iPerm[(it->first) + 1] = (it->second) + 1;
  }

  // Get the basis
  map<vector<PDT::Colour>,Trace_basis>::const_iterator bit =
    theBasisMap.find(basis);

  assert(bit != theBasisMap.end());

  const Trace_basis& Basis = bit->second;

  
  // Naive way: loop over every basis vector and rename the partons,
  //            then loop over the basis vectors to see which basis vector it became.
  //            This is sufficient for the current application, as it will
  //            only be used on the hard subprocess, which will never have
  //            a very large colour basis.
  map<size_t,size_t> indexMap;
  Col_amp Cai, Cai_re;
  Col_amp Caj;
  for ( size_t i = 0; i < dim; i++ ) {
    // Get the basis vector and change the leg numbering
    Cai = Basis.at(i);
    Cai_re = colorFunctions.rename_partons( Cai, iPerm );
    for ( size_t j = 0; j < dim; j++ ) {
      Caj = Basis.at(j);
      if ( Cai_re.at(0).cs == Caj.at(0).cs )
        indexMap[i] = j;
    }
  }
  // Check that the map was filled (every vector should have been
  // changed into a new vector)
  assert( indexMap.size() == dim );

  return indexMap;
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

