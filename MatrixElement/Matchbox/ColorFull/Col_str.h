// -*- C++ -*-
/*
 *  Col_str.h
 *	Contains declaration of the class Col_str and associated types and operators.
 *  Created on: Jul 7, 2010
 *  Author: Malin Sjodahl
 */

#ifndef COLORFULL_Col_str_h
#define COLORFULL_Col_str_h

#include "Quark_line.h"

namespace ColorFull {

/// For containing a vector (or list) of Quark_lines
/// the color information part of a Col_str.
/// The col_str is a product of Quark_lines,
/// contained in a vector of quark-lines.
typedef std::vector < Quark_line > col_str;


/// A class to contain ONE color structure, a direct product of Quark_lines,
/// multiplying a Polynomial, Poly.
/// The Quark_lines are stored as components in a vector, a col_str.
class Col_str {
public:

	/// Default constructor, leaves cs empty.
	Col_str(){};

	/// Constructor for setting the color structure using a string.
	/// Should be used as:
	/// "Col_str Cs("2*Nc*TR^(3) [{1,2,3,4}(5,6)(7,8)(9,10,11,12)]");",
	/// i.e. the argument should be a Polynomial * col_str.
	/// (The Polynomial should multiply the whole col_str,
	/// rather than a quark_line inside the [] brackets.)
	Col_str( const std::string str );

	/// Make a Col_str of a Quark_line.
	Col_str( Quark_line Ql ) {cs.push_back(Ql);}

	/// For containing the information about the color structure,
	/// a direct product of Quark_lines,
	/// contained in a vector of quark-lines.
	col_str cs;

	/// Polynomial factor multiplying the whole product of quark-lines.
	Polynomial Poly;

	/// Returns the Quark_line at place i.
	const Quark_line & at( int i ) const {return cs.at(i);}

	/// Returns the Quark_line at place i.
	Quark_line & at( int i ) {return cs.at(i);}

	/// Returns the parton at place j in Quark_line i.
	int at( int i, int j ) const;

	/// The size of the col_str
	uint size() const{ return cs.size(); }

	/// Is the col_str empty?
	bool empty() const { return cs.empty(); }

	/// Erase information in col_str.
	void clear() { cs.clear(); }

	/// Erases the Quark_line at place i.
	void erase( int i );

	/// Erases the parton at place i, j.
	void erase( int i, int j );

	/// Erases a parton at location place.
	void erase(std::pair<int, int> place);

	/// Appends a Quark_line to data member cs.
	void append( Quark_line Ql ) { cs.push_back( Ql ); }

	/// To insert the parton part_num in quark_line i
	/// at place j.
	void insert( int i, int j, int part_num );

	/// Append the content of a col_str to the cs of the Col_str.
	void append( col_str cs_in );

	/// Function for reading in the Col_str from the file filename.
	void read_in_Col_str( std::string filename );

	/// Function for writing out the Col_str to a file
	/// with name filename.
	void write_out_Col_str( std::string filename ) const;

	/// Locates the parton with number part_num in a Col_str.
	std::pair<int, int> find_parton( int part_num ) const;

	/// Function for telling if the partons p1 and p2 are neighbors.
	bool neighbor( int p1, int p2 ) const;

	/// Function for telling if parton p2 stands to the right of parton p1.
	bool right_neighbor( int p1, int p2 ) const;

	/// Function for telling if parton p2 stands to the left of parton p1.
	bool left_neighbor( int p1, int p2 ) const;

	/// Replaces the parton index old_ind with new_ind.
	void replace( int old_ind, int new_ind );

	/// Finds out if a parton is a quark, anti-quark or a gluon,
	/// returns "q", "qbar" or "g" respectively.
	/// This function does NOT loop over all partons, but assumes
	/// that the parton is a gluon if the Quark_line is closed,
	/// or if the Quark_line is open, but p cannot be found in the ends.
	std::string find_kind( int p ) const;

	/// Checks if the amplitude only has gluons, i.e. if all Quark_lines are closed.
	bool gluons_only() const;

	/// Counts the number of gluons in a Col_str.
	/// Counts all gluon indices, both free and contractable.
	int n_gluon() const;

	/// Counts the number of quarks (=number of anti-quarks) in a Col_str.
	/// Counts all quark indices, both free and contracted.
	int n_quark() const;

	/// Normal orders the Col_str by first
	/// normal order individual Quark_lines
	/// and then normal order different Quark_lines in the cs.
	/// For the ordering see the member function smallest in this
	/// class and in the Quark_line class.
	void normal_order();

	/// Finds out the "smallest" Col_str of two Col_strs, i.e.
	/// which Col_str should stand first in a normal ordered Col_amp or basis.
	/// Returns 1, if Cs1 should stand before Cs2
	/// and 2 if Cs2 should stand before Cs1.
	/// Both Col_strs have to be normal ordered for the result to be unique.
	/// The Col_strs are ordered by
	/// (1) number of Quark_lines
	/// (2) if the Quark_line at place 0,1,2... is open or not
	/// (3) the size of the Quark_line at place 1,2,3...
	/// (4) the parton numbers in the Quark_lines at place 1,2,3...,
	/// i.e. first the first parton in the first Quark_line is checked
	/// and last the last parton in the last Quark_line.
	/// The function returns 0 if Cs1=Cs2.
	int smallest( const Col_str & Cs1, const Col_str & Cs2 ) const;

	/// Returns the length of the longest Quark_line in the Col_str.
	int longest_quark_line() const;

	/// Removes Quark_lines with only one gluon as Tr(t^a)=0.
	void remove_1_rings();

	/// Removes Quark_lines without partons, equal to Nc (closed) or 1 (open).
	void remove_0_rings();

	/// Removes 0 and 1-rings,
	/// moves factors multiplying the individual Quark_lines to
	/// multiply the col_str instead (i.e., being stored in Poly)
	/// simplifies the Polynomial and normal orders the quark_lines.
	void simplify();

	/// Function for conjugating the Col_str by conjugating each Quark_line in cs,
	/// as well as the Polynomial Poly.
	void conjugate();

	/// Contracts neighboring and next to neighboring gluons in each
	/// Quark_line in the Col_str, starting with contracting neighbors.
	/// This function should only be used on Col_strs with only closed Quark_lines.
	void contract_next_neighboring_gluons( );

	/// Function for contracting gluon indices in closed Quark_lines with only 2 gluons.
	/// This removes the 2-ring, replaces one of the gluon indices,,
	/// and multiplies with a factor tr[t^a t^a]=TR (no sum),
	/// only intended for fully contractable Col_strs.
	void contract_2_rings( );

	/// Function for contracting quarks between two color structures Cs1 and Cs2.
	/// The result is stored in the Col_str itself.
	void contract_quarks( const Col_str Cs1, const Col_str Cs2 );


private:

	/// Function to tell which quark_line should stand first in normal order.
	/// The quark_lines at place i1 and i2 are compared.
	/// Returns i1 if i1 should stand first, i2 if i2 should stand first
	/// and i1 if the quark_lines are equal.
	int compare_quark_lines( int i1, int i2 ) const;

	/// Function to allow setting the color structure by using a string.
	/// Used by string constructor and by
	void Col_str_of_str( const std::string str  );

	/// Setting the col_str member cs
	/// (i.e. the non-Polynomial information) using a string,
	/// used by Col_str_of_str and read_in_Col_str.
	void col_str_of_str( std::string );

}; //end class Col_str


/// Define the operator == for two col_str's.
/// The col_str's must have equal length and
/// all Quark_lines must be the same
/// (i.e. have same Polynomial and same parton ordering).
/// The quark_lines are NOT normal ordered before comparison.
bool operator==(const col_str & cs1, const col_str & cs2);

/// Define the operator != for two col_str's.
/// Returns false if cs1==cs2 and false otherwise.
bool operator!=(const col_str & cs1, const col_str & cs2);

/// Define the operator << for col_str
std::ostream& operator<<( std::ostream& out, const col_str & cs );

/// Define the operator * for Col_str and int.
Col_str operator*( const Col_str & Cs, const int i );

/// Define the operator * for int and Col_str.
Col_str operator*( const int i, const Col_str & Cs );

/// Define the operator * for Col_str and double.
Col_str operator*( const Col_str & Cs, const double d );

/// Define the operator * for double and Col_str.
Col_str operator*( const double d, const Col_str & Cs );

/// Define the operator * for Col_str and cnum.
Col_str operator*( const Col_str & Cs, const cnum c);

/// Define the operator * for cnum and Col_str.
Col_str operator*( const cnum c, const Col_str & Cs );

/// Define the operator * for Col_str and Monomial.
Col_str operator*( const Col_str & Cs, const Monomial & Mon );

/// Define the operator * for Monomial and Col_str.
Col_str operator*( const Monomial & Mon, const Col_str & Cs);

/// Define the operator * for Col_str and Polynomial.
Col_str operator*( const Col_str & Cs, const Polynomial & Poly );

/// Define the operator * for Polynomial and Col_str.
Col_str operator*( const Polynomial & Poly, const Col_str & Cs );

/// Define the operator * for a Col_str and a Quark_line, adding Ql and multiplying Polynomial info.
Col_str operator*( const Col_str & Cs, const Quark_line & Ql);

/// Define the operator * for a Quark_line and a Col_str, adding Ql and multiplying Polynomial info.
Col_str operator*( const Quark_line & Ql, const Col_str & Cs );

/// Define the operator * for two Col_str's, adding Ql and multiplying Polynomial info.
Col_str operator*( const Col_str & Cs1, const Col_str & Cs2 );

/// Define the operator * for two Quark_lines. Clearly the result cannot be contained in
/// a Quark_line, but needs (at least) a Col_str.
Col_str operator*( const Quark_line & Ql1, const Quark_line & Ql2 );

/// Define the operator << for Col_str.
std::ostream& operator<<(std::ostream& out, const Col_str & Cs);

/// Define the operator == for two Col_str's
/// the col_strs must be equal, and the
/// Polynomials must be equal. Note that for Polynomials cf+Nc!=Nc+cf.
/// Function does NOT first normal order Col_strs.
bool operator==( const Col_str & Cs1, const Col_str & Cs2 );

/// Define the operator != for two Col_str's, the negation of ==.
/// Function does NOT first normal order Col_strs.
bool operator!=( const Col_str & Cs1, const Col_str & Cs2 );

}

#endif /* COLORFULL_Col_str_h */
