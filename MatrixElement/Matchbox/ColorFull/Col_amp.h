// -*- C++ -*-

/*
 * Col_amp.h
 *	Contains the declarations of the class Col_amp, related types and operators
 *  Created on: Jul 7, 2010
 *  Author: Malin Sjodahl
 */

#ifndef COLORFULL_Col_amp_h
#define COLORFULL_Col_amp_h

#include "Col_str.h"

namespace ColorFull {

/// Define a type to contain a linear combination of color structures,
/// a col_amp contains the actual color amplitude in a Col_amp.
typedef std::vector<Col_str> col_amp;


/// The full color amplitude is Scalar + Cs1+Cs2+Cs3...
/// Col_amp is a class to contain info on several Col_strs, a color amplitude.
class Col_amp {

public:

	/// Default constructor
	Col_amp() {Scalar = Scalar * 0;}

	/// Constructor taking a string as argument.
	/// The string should be of form
	/// Polynomial1* col_str1 + Polynomial2*col_str2,
	/// for example:
	/// Col_amp Ca("13*Nc*[(1, 3, 4, 2)] +2 TR 5 Nc^(-3) [(1, 4) (3,  2)]").
	/// (The Polynomials should multiply the whole col_strs in square brackets,
	/// rather than a quark_line inside the [] brackets.)
	Col_amp( const std::string str ){Col_amp_of_str( str );}

	/// Constructor converting a Col_str to a Col_amp.
	Col_amp( Col_str Cs ) {
		Scalar = Scalar * 0;
		ca.push_back(Cs);
	}

	/// To actually contain the info of the Col_strs, ca=Cs1+Cs2+Cs3+... .
	/// Technically the ca is a vector of Col_strs, a col_amp.
	col_amp ca;

	/// Scalar is Polynomial for collecting color factors when the color structure
	/// has been fully contracted.
	/// The full color amplitude is Scalar + Cs1+Cs2+Cs3..., the Polynomial should thus be non-zero
	/// only if all indices can be contracted.
	Polynomial Scalar;

	/// Reads in the Col_amp to the member ca from the file filename.
	/// (This is for reading in an actual color amplitude,
	/// nothing is read in the the Polynomial member scalar.)
	void read_in_Col_amp( std::string filename );

	/// Function for writing out the Col_amp to a file
	/// with name filename.
	void write_out_Col_amp( std::string filename ) const;

	/// Returns the Col_str at place i.
	const Col_str & at( int i ) const{ return ca.at(i); }

	/// Returns the Col_str at place i.
	Col_str & at( int i ) {return ca.at(i);}

	/// The size of the col_amp ca.
	uint size() const {return ca.size();}

	/// Is the col_amp empty?
	bool empty() const { return ca.empty(); }

	/// Erase information in col_amp.
	void clear() { ca.clear(); }

	/// Erases the Col_str at place i.
	void erase( int i );

	/// Appends the Col_strs in ca_in to the col_amp member ca.
	void append( col_amp ca_in );

	// Functions for probing the Col_amp

	/// Checks if the Col_amp only contains gluons, i.e., if all Quark_lines are closed.
	bool gluons_only() const;

	/// Returns the number of gluons in the Col_amp as the number of gluons in the first Col_str.
	/// Note that the other Col_strs could have a different number of (contracted) gluons.
	/// (Intended for tree-level Col_ams with only one Col_str.)
	int n_gluon() const;

	/// Returns the number of quarks  in the Col_amp as the number of quarks in the first Col_str.
	/// Note that the other Col_strs could have a different number of (contracted) quarks.
	/// (Intended for tree-level Col_ams with only one Col_str.)
	int n_quark() const;

	/// Returns the number of quarks in the Col_amp after checking that each Col_str
	/// has the same number of quarks.
	int n_quark_check() const;

	/// Returns the number of gluons in the Col_amp after checking that each Col_str
	/// has the same number of gluons.
	int n_gluon_check() const;

	/// Returns the length of the longest Quark_line in any Col_str.
	int longest_quark_line() const;

	// Functions for manipulating the Col_amp

	/// Remove Col_strs with quark_lines with just 1 gluon, they are 0 as Tr[t^a]=0.
	void remove_1_rings();

	/// Remove quark_lines with no gluons, they are N if closed, and defined to be 1 if open.
	void remove_0_rings();

	/// Removes empty Col_strs, an empty Col_str means that all indices have been contracted,
	/// so the Col_str is equal to its Polynomial, which is moved to the scalar part
	/// of the Col_amp.
	void remove_empty_Col_strs();

	/// Compares col_strs in a Col_amp, to collect similar col_strs,
	/// and only store once in ca.
	void collect_col_strs();

	/// Normal orders all col_strs in ca.
	void normal_order_col_strs();

	/// Normal orders the individual col_strs and then
	/// orders the Col_strs using the order defined in
	/// the Col_str function smallest.
	void normal_order();

	/// Function for simplifying an amplitude,
	/// removes 0 and 1-rings,
	/// compares col_strs,
	/// removes Col_strs multiplying 0 and
	/// simplifies Polynomials of the individual Col_strs.
	void simplify();

	/// Function for taking the conjugate of the Col_amp
	/// by conjugating each Col_str in ca and the
	/// Polynomial member Scalar.
	void conjugate();

	/// Contracts up to next to neighboring gluons in each Quark_line
	/// in each Col_str in each Col_amp, only intended for closed Quark_lines.
	void contract_next_neighboring_gluons( );

	/// Contract closed Quark_lines with only 2 gluons in
	/// each Quark_line in each Col_str in the Col_amp.
	/// This removes the 2-ring, replaces one of the gluon indices and
	/// multiplies with a factor tr[t^a t^a]=(1/2) (no sum).
	void contract_2_rings( );

	/// Function for contracting gluon indices within the Quark_lines.
	/// Checks only for ONE pair in each Quark_line, i.e.,
	/// if several gluon indices appear, only one pair is contracted
	/// in each Quark_line, only intended for closed Quark_lines.
	void contract_Quark_line_gluons( );

	/// Contracts one gluon, the first gluon in first Quark_line (in each Col_str),
	/// only intended for closed Quark_lines.
	void contract_a_gluon( );

	/// Function for contracting all gluon indices in a Col_amp,
	/// only intended for closed Quark_lines.
	void contract_all_gluons( );

	/// Function for contracting the (anti-)quarks in Ca1 with those
	/// in Ca2. The results is saved in this Col_amp.
	void contract_quarks( const Col_amp & Ca1, const Col_amp & Ca2 );

private:


	/// Contracts one gluon, the first gluon in first Quark_line, only intended
	/// for closed Quark_lines. This function is a member of Col_amp as
	/// the result is contained in this Col_amp.
	void contract_a_gluon( Col_str & Cs );

	/// Function for contracting gluon indices within the same Quark_line.
	/// Checks only for ONE pair, i.e. if several gluon indices appear
	/// within the Quark_line, only one pair is contracted.
	/// Only intended for closed Quark_lines.
	/// This function is a member of Col_amp as the result is a Col_amp.
	void contract_Quark_line_gluons( Quark_line & Ql );

	/// Function for contracting gluon indices within the same Quark_line.
	/// Checks only for ONE pair in each Quark_line in each Col_str,
	/// i.e. if several gluon indices appear within the Quark_line,
	/// only one pair is contracted.
	/// The function is only intended for closed Quark_lines.
	/// This function is a member of Col_amp as it saves the result in
	/// this Col_amp.
	void contract_Quark_line_gluons( Col_str & Cs );

	/// Function for contracting all gluon indices in a Col_str,
	/// assumes quarks already contracted and is only intended for
	/// closed Quark_lines. This function is a member of Col_amp as
	/// the result is contained in this Col_amp.
	void contract_all_gluons( Col_str & Cs );

	/// Converts a text string to a Col_amp,
	/// used by string constructor, and by read_in_Col_amp.
	/// The string should be of form
	/// Polynomial1* col_str1 + Polynomial2 col_str2,
	/// for example:
	/// Col_amp Ca("13*Nc*[(1,3,4,2)] +2 TR 5 Nc^(-3) [(1,4) (3,2)]").
	void Col_amp_of_str( const std::string str );

};

// Define operators involving Col_amp

/// Define the operator << for col_amp
std::ostream& operator<<( std::ostream& out, const col_amp & ca );

/// Define the operator << for Col_amp.
std::ostream& operator<<( std::ostream& out, const Col_amp & Ca );

/// Define the operator == for two Col_amps.
bool operator==( const Col_amp & Ca1, const Col_amp & Ca2 );

/// Define the operator != for two Col_amps.
bool operator!=( const Col_amp & Ca1, const Col_amp & Ca2 );

/// Define the operator + for Col_amp and Col_str, adds Col_str Cs to ca.
Col_amp operator+( const Col_amp & Ca, const Col_str & Cs );

/// Define the operator + for Col_str and Col_amp, adds Col_str Cs to ca.
Col_amp operator+( const Col_str & Cs, const Col_amp & Ca );

/// Define the operator + for two Col_amps.
/// Adds Scalars, and adds Col_strs.
Col_amp operator+( const Col_amp & Ca1, const Col_amp & Ca2 );

/// Define the operator += for two Col_amps.
Col_amp operator+=( Col_amp & Ca1, const Col_amp & Ca2 );

/// Define the operator - for two Col_amps.
/// Subtract Scalar of Ca2 from Ca1, and subtracts (=appends with minus sign) Col_strs.
Col_amp operator-( const Col_amp & Ca1, const Col_amp & Ca2 );

/// Define the operator * for Col_amps and integers.
Col_amp operator*( const Col_amp & Ca, const int i );

/// Define the operator * for integers number and Col_amp.
Col_amp operator*( const int i, const Col_amp & Ca );

/// Define the operator * for Col_amps and complex number.
Col_amp operator*( const Col_amp & Ca, const cnum c );

/// Define the operator * for complex number and Col_amp.
Col_amp operator*( const cnum c, const Col_amp & Ca );

/// Define the operator * for Col_amp and double.
Col_amp operator*( const Col_amp & Ca, const double d );

/// Define the operator * for double and Col_amp.
Col_amp operator*( const double d, const Col_amp & Ca );

/// Define the operator * for Col_amp and Monomial.
Col_amp operator*( const Col_amp & Ca, const Monomial & Mon );

/// Define the operator * for Monomial and Col_amp.
Col_amp operator*( const Monomial & Mon, const Col_amp & Ca );

/// Define the operator * for Col_amp and Polynomial.
Col_amp operator*( const Col_amp & Ca, const Polynomial & Poly );

/// Define the operator * for Monomial and Col_amp.
Col_amp operator*( const Polynomial & Poly, const  Col_amp & Ca );

/// Define the operator *= for Col_amp and Col_str.
Col_amp operator*( const Col_amp & Ca, const Col_str & Cs );

/// Define the operator * for Col_amp and Col_str.
Col_amp operator*( const Col_str & Cs, const Col_amp & Ca );

/// Define the operator * for Col_amps.
Col_amp operator*( const Col_amp & Ca1, const Col_amp & Ca2 );

/// Define the operator *= for two Col_amps.
Col_amp operator*=( Col_amp & Ca1, const Col_amp & Ca2 );

}// end namespace ColorFull

#endif /* COLORFULL_Col_amp_h */
