// -*- C++ -*-
/*
 *  Trace_type_basis.h
 *	Contains the declarations of the class Trace_type_basis, related types and operators.
 *  Created on: Aug 9, 2012
 *  Author: Malin Sjödahl
 */

#ifndef COLORFULL_Trace_type_basis_h
#define COLORFULL_Trace_type_basis_h

#include "Col_basis.h"


namespace ColorFull {

/// Trace_type_basis is used for the common features of
/// Trace_basis and Tree_level_gluon_basis which both inherit
/// from Trace_type_basis.
class Trace_type_basis:public Col_basis{

public:

	/// Default constructor.
	Trace_type_basis():Col_basis(){
		// The maximal number of quark_lines is
		// nq+n_g/2 for a Trace_basis and 1 for a Tree_level_gluon basis,
		// but at this time nq and n_g are not known.
		max_ql=0;
	}

	/// A function for decomposing the color amplitude ca in the basis,
	/// returning the result as a Polynomial.
	Poly_vec decompose( const Col_amp & Ca );


	/// Function for finding the new vector numbers in the new basis
	/// (this trace basis)
	/// after radiating a new gluon from the parton emitter.
	/// The old color structure is Cs, and after emission a linear combination
	/// of new basis vectors is obtained.
	/// For emission from a quark or an anti-quark there is only one resulting color
	/// structure, and -1 is returned in the place of the absent color structure.
	/// The second vector, where the new gluon is inserted before the emitter,
	/// comes with a minus sign in the new total amplitude.
	std::pair<int, int>  new_vector_numbers( const Col_str & Cs, int emitter ) ;


	/// This function is intended for tree-level processes with at most 2 qqbar-pairs.
	/// It finds the new vector numbers in the basis for n_p+1 partons
	/// after radiating a new gluon from the parton emitter.
	/// This function does not actually use the cb, but only calculates
	/// the basis vector number, which makes it much quicker than the general version.
	/// The old vector has number old_num, and there were, before emission,
	/// nq quarks (+ nq anti-quarks) and n_g-1 gluons, i.e. n_p=2 nq+ n_g-1.
	/// For emission from q or qbar there is only one resulting color structure,
	/// and -1 is returned in the place of the absent color structure.
	/// The second vector, where the new gluon is inserted before the emitter
	/// comes with a minus sign in the new total amplitude.
	/// The function has been explicitly tested against its sister function
	/// for initial states with 2 qqbar-pairs and up to 5 gluons,
	/// 1qqbar-pair and up to 7 gluons and 0 qqbar-pairs and up to 8 gluons.
	std::pair<int, int>  new_vector_numbers( int old_num, int emitter, int n_loop ) const;


protected:

	/// The maximal number of quark-lines allowed in the basis.
	/// This is used for constructing bases that only are valid
	/// up to a certain order in QCD, such that unused information
	/// need not be carried around.
	int max_ql;

	/// Function for finding the new vector number in the basis with n_p+1 partons (this basis)
	/// after inserting a new gluon with larger parton number at the place place.
	/// This function is only intended for the special cases of
	/// 0 qqbar-pairs and Tree_level_gluon_basis
	/// or Trace_basis and 1-2 qqbar-pairs at tree level, i.e. n_loop must be 0.
	/// This function doesn't actually use the cb, but only calculates the
	/// the basis vector number using find_parton.
	/// The old vector has number old_num, and there were, before emission
	/// nq quarks (+ nq anti-quarks) and n_g-1 gluons, i.e. n_p=2 nq + n_g-1.
	/// The function has been explicitly tested for initial states with 2qqbar-pairs and up to 5 gluons,
	/// 1qqbar-pair and up to 7 gluons, 0 qqbar-pairs and up to 8 gluons.
	int new_vector_number( int old_num, std::pair<int,int> place, int n_loop ) const;

	/// This function is only intended for special case of:
	/// 0 qqbar-pairs and Tree_level_gluon_basis
	/// or Trace_basis and 1-2 qqbar-pairs at tree level, i.e. n_loop must be 0.
	/// It locates the parton parton in the normal ordered basis,
	/// given the number of the vector vec_num, and the number of quarks and gluons in the basis.
	/// The function has been explicitly tested for Trace_basis with 1 qqbar-pair and up to 8 gluons
	/// and 2 qqbar-pairs and up to 7 gluons, and for Tree_level_gluon_bases with up to 9 gluons.
	/// The arguments n_quark and n_gluon has to be provided as it may be desirable to
	/// use the function with a different number of quarks and gluons than in the basis itself.
	std::pair<int,int> find_parton( int parton, int vec_num, int n_quark, int n_gluon, int n_loop ) const;

};

}

#endif /* COLBASIS_H_ */

