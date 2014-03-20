// -*- C++ -*-

/*
 * Poly_vec.h
 * Contains the declaration of the class Poly_vec related types and operators.
 * Created on: Aug 5, 2013
 * Author: Malin Sjodahl
 */

#ifndef COLORFULL_Poly_vec_h
#define COLORFULL_Poly_vec_h


#include "Polynomial.h"


namespace ColorFull {

/// To contain a vector of Polynomials.
typedef std::vector< Polynomial> poly_vec;


/// Class for containing vector of Polynomials, and functions
/// for Polynomial vectors.
class Poly_vec {


public:

	/// Default constructor, sets nothing.
	Poly_vec() {}

	/// Make a Poly_vec of a poly_vec.
	Poly_vec( poly_vec poly_v ){pv=poly_v;}

	/// To actually contain the polynomial information.
	poly_vec pv;

	/// Returning Polynomial at place i.
	const Polynomial& at( int i ) const { return pv.at(i); }

	/// Returning Polynomial at place i.
	Polynomial& at( int i ) { return pv.at(i); }

	/// Return the number of Polynomials in the vector,
	/// i.e., the size of the member pv.
	uint size() const { return pv.size(); }

	/// Erase information in vector.
	void clear() { pv.clear(); }

	/// Appends a Polynomial to data member pv.
	void push_back( Polynomial Poly ) { pv.push_back( Poly ); }

	/// Is the vector empty?
	bool empty() const { return pv.empty(); }

	/// Remove CF in the poly_vec member pv, i.e., replace CF by
	/// TR (Nc^2-1)/(2 Nc).
	void remove_CF();

	/// Normal order all Polynomials in the poly_vec member pv.
	/// (Uses the Polynomial.normal_order function.)
	void normal_order();

	/// Simplifies all polynomials in the poly_vec member pv.
	/// (Uses the simplify member function in Polynomial).
	void simplify();

	/// Conjugates the Poly_vec.
	void conjugate();

	/// Reads in a Polynomial vector of form {Poly1, Poly2,...}
	/// to the member pv from the file filename.
	/// Comments starting with # are allowed
	/// at the top of the file.
	void read_in_Poly_vec( std::string filename );

	/// Writes out the vector to the file filename.
	void write_out_Poly_vec( std::string filename ) const;

};

/// Operator << for poly_vec.
std::ostream& operator<<( std::ostream& out, const poly_vec & poly_v );

/// Operator << for Poly_vec.
std::ostream& operator<<( std::ostream& out, const Poly_vec & Pv );


}


#endif /* COLORFULL_Poly_vec_h */
