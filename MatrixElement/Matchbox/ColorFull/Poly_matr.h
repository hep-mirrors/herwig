// -*- C++ -*-

/*
 * Poly_matr.h
 * Contains the declaration of the class Poly_matr and related types and operators.
 * Created on: Aug 5, 2013
 * Author: Malin Sjodahl
 */

#ifndef COLORFULL_Poly_matr_h
#define COLORFULL_Poly_matr_h



#include "Poly_vec.h"


namespace ColorFull {

/// To contain a matrix of Polynomials, a vector of Poly_vec.
typedef std::vector< Poly_vec > poly_matr;


/// Class for containing a Polynomial matrix, and functions
/// for Polynomial matrices.
class Poly_matr {

public:

	/// Default constructor, leaves pm empty.
	Poly_matr() {}

	/// To actually contain the matrix of Polynomials.
	poly_matr pm;

	/// Returns the Poly_vec at place i.
	const Poly_vec& at( int i ) const {return pm.at(i);}

	/// Returns the Poly_vec at place i.
	Poly_vec& at( int i ) {return pm.at(i);}

	/// Returns the matrix element at i, j.
	Polynomial& at( int i, int j ) { return pm.at(i).pv.at(j);}

	/// Returns the matrix element at i, j.
	const Polynomial& at( int i, int j ) const { return pm.at(i).pv.at(j);}

	/// Is the matrix, stored in pm, empty?
	bool empty( ) const {return pm.empty();}

	/// Returns the size of the matrix, the number of Poly_vec's
	/// in the member pm.
	uint size( ) const {return pm.size();}

	/// Erases the matrix information.
	void clear()  {pm.clear();}

	/// Appends a Poly_vec to data member pm.
	void append( Poly_vec Pv ) {pm.push_back( Pv );}

	/// Remove CF in the poly_matr member pm, i.e., replace CF by
	/// TR (Nc^2-1)/Nc.
	void remove_CF();

	/// Normal orders all polynomials in the poly_matr member pm,
	/// (uses Polynomial.normal_order.)
	void normal_order();

	/// Simplifies all polynomials in the poly_matr member pm,
	/// (uses Polynomial.simplify.)
	void simplify();

	/// Conjugates the matrix.
	void conjugate();

	/// Reads in the matrix from the file filename.
	/// The file should be of the format
	/// {{Poly11,...,Poly1n},
	/// ...,
	/// {Polyn1,...,Polynn}},
	/// and may contain comment lines starting with # at the top.
	void read_in_Poly_matr( std::string filename );

	/// Writes out the matrix to the file filename.
	void write_out_Poly_matr( std::string filename ) const;

};

/// Operator << for poly_matr.
std::ostream& operator<<( std::ostream& out, const poly_matr & pm );

/// Operator << for poly_vec.
std::ostream& operator<<( std::ostream& out, const Poly_matr & Pm );

/// Define the operator == for Poly_matr,
/// each Poly_vec has to be identical.
bool operator==( const Poly_matr & Pm1, const Poly_matr & Pm2 );

/// Define the operator == for Poly_matr,
/// each Poly_vec has to be identical.
bool operator!=( const Poly_matr & Pm1, const Poly_matr & Pm2 );

}


#endif /* COLORFULL_Poly_matr_h */
