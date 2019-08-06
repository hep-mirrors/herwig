// -*- C++ -*-
/*
 * Quark_line.h
 *	Contains declaration of the class Quark_line and associated types and operators.
 *  Created on: Jul 7, 2010
 *  Author: Malin Sjodahl
 */

#ifndef COLORFULL_Quark_line_h
#define COLORFULL_Quark_line_h

#include "Polynomial.h"


namespace ColorFull {


/// Define a type to contain a quark-line with gluons attached,
/// the actual color information about how quarks are ordered in a
/// Quark_line.
typedef std::vector <int> quark_line;

/// A class to contain one quark-line with gluons attached
/// multiplying a Polynomial, Poly.
/// The partons are ordered as {q,g1,g2...gn,qbar} for an open quark-line
/// containing a q and a qbar, or (g1,g2...gn) for a closed
/// quark-line with only gluons.
/// The parton numbers are stored as components in a vector, a quark_line,
/// whereas the member open contains information about it the quark-line is
/// closed as a trace (true) or open, false.
class Quark_line {

public:

	/// Constructor used to set the color structure using a string.
	/// The string should be of form Polynomial*quark_line,
	/// used as "Quark_line Ql("5*TR*Nc^2 {1,6,7,2}");"
	/// for an open Quark_line with a quark with
	/// number 1, two gluons with number 6 and 7, and a qbar with number 2.
	/// For a closed Quark_line with 3 gluons the syntax is
	/// "Quark_line Ql("(1,2,3)");".
	/// The integers should be positive.
	/// The Polynomial should be in such a shape that it is readable by the
	/// Polynomial( std::string ) constructor.
	Quark_line( const std::string str );

	/// Default constructor.
	Quark_line();

	/// To actually contain the color information,
	/// in order {q, g1, ... gn, qbar} or (g1, ..., gn).
	quark_line ql;

	/// Polynomial factor, multiplying the quark_line.
	Polynomial Poly;

	/// Is the string open, with a q in the beginning and a qbar in the end, or not?
	bool open;

	/// Function for reading in the Quark_line from the file filename,
	/// uses Quark_line_of_str.
	void read_in_Quark_line( std::string filename );

	/// Function for writing out the Quark_line to a file
	/// with name filename.
	void write_out_Quark_line( std::string filename ) const;

	/// Returns the parton at place j.
	/// For closed quark_lines j may be between -size and 2*size.
	int at( int j ) const;

	/// The size of the quark_line.
	uint size() const{ return ql.size();}

	/// Erase information in quark_line ql.
	void clear() { ql.clear(); }

	/// Is the quark_line empty?
	bool empty() const { return ql.empty(); }

	/// If the quark_line is open, there is nothing to do,
	/// else order with smallest gluon index first
	/// (use that the trace is cyclic).
	void normal_order();

	/// To erase the parton at place i.
	void erase( int i );

	/// Conjugates the Quark_line by reversing the quark_line ql
	/// and conjugating the Polynomial Poly.
	void conjugate();

	/// Appends parton p to the Quark_line.
	void append( int p ) { ql.push_back( p ); }

	/// Appends a whole quark_line to the Quark_line.
	void append( const std::vector<int> & in_ql );

	/// Prepends parton p to the Quark_line.
	void prepend( int p );

	/// Prepends a whole quark_line to the Quark_line.
	void prepend( std::vector<int> in_ql );

	/// Inserting parton p at place j.
	void insert( int j, int p );

	/// Returns a Quark_line where the ql member is changed to contain
	/// only partons before place j.
	Quark_line before( int j ) const;

	/// Returns a Quark_line where the ql member is changed to contain
	/// only partons after place j.
	Quark_line after( int j ) const;

	/// Function for splitting a closed Quark_line into two Quark_lines.
	/// The gluons at j1 and j2 are removed in the split.
	/// May create 1-rings and 0-rings.
	std::pair<Quark_line, Quark_line> split_Quark_line( int j1, int j2 ) const;

	/// Function for finding the "smallest" Quark_line of Ql1 and Ql2,
	/// used for deciding which Quark_line should stand first while normal ordering.
	/// Does NOT first normal order the Quark_lines.
	/// If only one is open, that Quark_line should stand first.
	/// If both are open or both are closed, the longest Quark_line should stand first.
	/// If the size is the same, the Quark_line with smallest starting number should stand first.
	/// If the first number is the same, check the 2nd number, then the 3rd...
	/// 1 is returned if Ql1 should stand first, and 2 if Ql2 should stand first.
	/// If Ql1==Ql2, 0 is returned.
	int smallest( const Quark_line & Ql1 , const Quark_line & Ql2 ) const;

	/// Contracts neighboring gluons in the Quark_line starting at j,
	/// only intended for closed Quark_lines.
	void contract_neighboring_gluons( int j );

	/// Contracts neighboring gluons in a Quark_line starting at place 0,
	/// and checking all neighbors, only intended for closed Quark_lines.
	void contract_neighboring_gluons( );

	/// Contracts neighboring and next to neighboring gluons in the Quark_line,
	/// starting at place j (i.e. checking gluon j and j+2).
	/// Also looks for new neighbors, only intended for closed Quark_lines.
	void contract_next_neighboring_gluons( int j );

	/// Contracts neighboring and next to neighboring gluons in the Quark_line,
	/// starting with contracting neighbors, only intended for closed Quark_lines.
	void contract_next_neighboring_gluons( );


private:

	/// Function used to set the color structure using a string.
	/// The string should be of form Polynomial*quark_line,
	/// used as "Quark_line Ql("Polynomial {5,6,7}");"
	/// for an open Quark_line with a quark with
	/// number 5, a gluon with number 6 and a qbar with number 7, and as
	/// "Quark_line Ql("(5,6,7)");" for 3 gluons attached to a quark line.
	/// The Polynomial should be in such a shape that it's readable by the
	/// Polynomial( std::string ) constructor.
	void Quark_line_of_str( const std::string str );

	/// To make it easy to define a Quark_line using a string,
	/// used by Quark_line_of_str(std::string str)
	/// to set the quark_line member ql.
	void quark_line_of_str( std::string str );

};

/// Define the operator * for Quark_line and int.
/// The Polynomial of the Quark_line is multiplied with i.
Quark_line operator*( const Quark_line & Ql, const int i );

/// Define the operator * for Quark_line and int,
/// returns Ql*i.
Quark_line operator*( const int i, const Quark_line & Ql );

/// Define the operator * for Quark_line and cnum.
/// The Polynomial of the Quark_line is multiplied with c.
Quark_line operator*( const Quark_line & Ql, const cnum c );

/// Define the operator * for Quark_line and cnum,
/// returns Ql*c.
Quark_line operator*( const cnum c, const Quark_line & Ql );

/// Define the operator * for Quark_line and double.
/// The Polynomial of the Quark_line is multiplied with d.

Quark_line operator*( const Quark_line & Ql, const double d );
/// Define the operator * for Quark_line and double,

/// returns Ql*d.
Quark_line operator*( const double d, const Quark_line & Ql );

/// Define the operator * for Quark_line and Monomial.
/// The polynomial of the Quark_line is multipled with Mon.
Quark_line operator*( const Quark_line & Ql, const Monomial & Mon );

/// Define the operator * for Quark_line and Monomial.
/// returns Ql*Mon.
Quark_line operator*( const Monomial & Mon, const Quark_line & Ql );

/// Define the operator * for Quark_line and Polynomial.
/// The Polynomial of the Quark_line is multiplied with Poly.
Quark_line operator*( const Quark_line & Ql, const Polynomial & Poly );

/// Define the operator * for Quark_line and Polynomial
/// returns Ql*Poly.
Quark_line operator*( const Polynomial & Poly, const Quark_line & Ql );

/// Define the operator == for two Quark_lines
/// the quark_lines must be equal, the member variable open and the
/// Polynomials must be equal. Note that for Polynomials cf+Nc!=Nc+cf.
bool operator==( const Quark_line & Ql1, const Quark_line & Ql2 );

/// The negation of the Quark_line operator ==.
bool operator!=( const Quark_line & Ql1, const Quark_line & Ql2 );


/// Define the operator << for Quark_line.
std::ostream& operator<<( std::ostream& out, const Quark_line & Ql );


}// end namespace ColorFull

#endif /* COLORFULL_Quark_line_h */
