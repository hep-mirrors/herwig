// -*- C++ -*-
/*
 * Col_functions.h
 * Contains declarations of the class Col_str and associated types and operators
 * Author: Malin Sjodahl
 */

#ifndef COLORFULL_Col_functions_h
#define COLORFULL_Col_functions_h


#include "Col_amp.h"
#include "Poly_matr.h"
#include <list>
#include <map>
#include <memory>

namespace ColorFull {

using std::shared_ptr;

/// Library class containing functions for index contraction and
/// numerical evaluation.
/// This is where the parameters Nc, TR and CF are contained.
class Col_functions {

private:

/// The number of colors, used in numerical results,
/// changed by using set_Nc.
double Nc;

/// The trace convention Tr( t^a t^a )=TR (no sum).
/// The normalization of the SU(Nc) generators, to be used in numerical evaluation,
/// changed by using set_TR.
/// The value 1/2 corresponds to the Gell-Mann normalization.
double TR;

/// The value of CF=TR*Nc-TR/Nc, changed by using set_CF.
/// Note that CF can be changed independently of Nc.
double CF;

/// While evaluating leading terms one may want to keep the full value of CF for
/// TR(Nc^2-1)/Nc, or only keep the leading Nc term =TR*Nc (default).
/// full_CF is used by the Polynomial version of leading
/// (and hence also Poly_vec and Poly_matr versions etc).
/// The leading functions replaces CF by TR*Nc if full_CF is false (default)
/// while evaluating the leading terms.
/// If full_CF is true, CF is replaced by TR(Nc^2-1)/Nc.
/// Clearly this affects the result of subsequent numerical evaluation.
/// In the Col_basis class (and derived) the matrix version of leading
/// is used to evaluate scalar product matrices.
bool full_CF;


public:

/// Default constructor.
Col_functions()
  : Nc(3.0), TR(0.5), CF(4.0/3.0), full_CF( false )  {}

/// Set the number of colors.
/// The value of CF is adjusted accordingly.
void set_Nc( double n) {
  Nc = n;
  CF = TR*(Nc*Nc-1.)/Nc;
}

/// Set the normalization of the generators.
/// The value of CF is adjusted accordingly.
void set_TR( double tr) {
  CF *= tr/TR;
  TR = tr;
}

/// Set the value of CF.
/// The value of Nc is NOT adjusted accordingly.
void set_CF( double cf) {
  CF = cf;
}

/// Switch on/off full_CF.
void set_full_CF( bool is_full ) { full_CF = is_full; }

/// Returns the number of colors.
double get_Nc() const { return Nc; }

/// Returns the normalization of the generators,
/// tr(t^a t^b)=TR*delta^{a,b}.
double get_TR() const { return TR; }

/// Returns the value of CF.
double get_CF() const { return CF; }

/// Returns true, if full CF is used.
bool get_full_CF() const { return full_CF; }


/****************** Functions for leading terms *************************/
// The functions called leading(...) depend on the variable full_CF.
// As it would be messy to let each Polynomial carry around
// its own full_CF, these functions are kept here.

/// Function for finding the leading power of Nc in a Poly_vec,
/// i.e., the power of Nc plus the power of CF.
int  leading_Nc_pow( const Polynomial & Poly ) const;

/// Function for finding the leading power of Nc in a Poly_vec.
int  leading_Nc_pow( const Poly_vec & Pv ) const;

/*
/// Function for finding the leading power of Nc in a
/// vector of pointers to Polynomials.
int  leading_Nc_pow( const std::vector< shared_ptr<Polynomial> > & Pvp) const;
*/

/// Takes the leading Nc terms of a Polynonmial, i.e. Monomials with highest
/// power of Nc+CF. If full_CF is false (default), CF is replaced by TR Nc.
/// If full_CF is true CF is replaced by TR(Nc^2-1)/Nc.
Polynomial leading( const Polynomial & Poly ) const;

/// Take the leading part of a Poly_vec.
/// Keeps only Monomials with maximal power of CF plus Nc,
/// uses leading( const Polynomial & Poly).
/// If full_CF is false (default), CF is replaced by TR Nc.
/// If full_CF is true CF is replaced by TR(Nc^2-1)/Nc.
/// Note that taking the leading terms of a Poly_vec is not
/// the same as taking the leading terms in each Polynomial.
// Used only by Poly_matr version of leading
Poly_vec leading( const Poly_vec & Pv ) const;

/// Takes the leading part of a matrix of Polynomials,
/// keeping only those with maximal power of CF plus Nc.
/// If full_CF is false (default), CF is replaced by TR Nc.
/// If full_CF is true CF is replaced by TR(Nc^2-1)/Nc.
/// Note that taking the leading terms of a Poly_matr is not
/// the same as taking the leading terms in each Poly_vec.
// Used only once in Col_basis
Poly_matr leading( const Poly_matr & Pm ) const;

/*
/// Take the leading part of a Poly_vec, given a vector of pointers to the Polynomials.
/// Keeps only Monomials with maximal power of CF plus Nc.
// Currently never used
Poly_vec leading( const std::vector<shared_ptr<Polynomial> > & Pvp) const;

/// Take the leading part of a Poly_matr, given a vector of vector of pointers to the Polynomials.
/// Loops over Monomials in all Polynomials
/// and keeps only those with maximal power of CF plus Nc.
// used only by scalar_product_matrix_mem in Col_functions
dmatr leading( const std::vector< std::vector< shared_ptr<Polynomial> > > & Pm ) const;
*/

/********************* Functions for numerical evaluation *************************/
// These functions has to be kept in Col_functions class as they need numerical
// values for evaluation. Letting each Polynomial carry around its own Nc etc.
// would be messy.


/// Numerically evaluates a Monomial using the Nc, TR and CF data members.
cnum cnum_num( const Monomial & Mon ) const;

/// Numerically evaluates a Polynomial, using the Nc, TR and CF data members.
cnum cnum_num( const Polynomial & Poly ) const;

/// Numerically evaluates a Poly_vec (vector of Polynomial),
/// using cnum_num (Polynomial).
cvec cnum_num( const Poly_vec & Pv ) const;

/// Numerically evaluates a Poly_matr (vector of Poly_vec),
/// using cnum_num( Poly_vec ) for each Poly_vec.
cmatr cnum_num( const Poly_matr & Pm ) const;

/// Numerically evaluates a Monomial to a double using the Nc, TR and CF data members.
double double_num( const Monomial & Mon ) const;

/// Numerically evaluates a Polynomial to a double using the Nc, TR and CF data members.
double double_num( const Polynomial & Poly ) const;

/// Numerically evaluates a Poly_vec (vector of Polynomial)
/// using the Nc, TR and CF data members.
dvec double_num( const Poly_vec & Pv ) const;

/// Numerically evaluates a Poly_matr (vector of Poly_vec),
/// using the Nc, TR and CF data members.
dmatr double_num( const Poly_matr & Pm ) const;

/*
/// Returns a double vector. The argument is a vector of pointers to Polynomials.
dvec double_num( const std::vector<std::shared_ptr<Polynomial> > & Pv ) const;

/// Returns a double matrix. The argument is a vector of vector of pointers to Polynomials.
// (not used 14 06 08)
dmatr double_num( const std::vector<std::vector<std::shared_ptr<Polynomial> > > & Pm ) const;

/// To take the numerical value of a map.
// (not used 14 06 08)
std::map< std::string, double > double_num( std::map< std::string, std::shared_ptr<Polynomial> > mem_map ) const;

/// To take the numerical value of a map.
// (not used 14 06 08)
std::map< std::string, double > double_num( std::map< std::string, Polynomial > mem_map ) const;
*/

/// Numerically evaluates a Polynomial using the value of the data member Nc,
/// and stores in the format of a Polynomial with only one term with only a numerical part.
Polynomial Polynomial_cnum_num( const Polynomial & Poly ) const;

/// Numerically evaluates a Poly_vec (vector of Polynomial)
/// and stores in the form of a Poly_vec, uses polynomial_cnum_num( Pv.at( p ) ).
/// for each Polynomial.
Poly_vec Poly_vec_cnum_num( const Poly_vec & Pv ) const;

/// Numerically evaluates a Poly_matr (vector of Poly_vec)
/// and stores in the form of a Poly_matr.
Poly_matr Poly_matr_cnum_num( const Poly_matr & Pm ) const;


/****************** Functions for scalar products *************************/

/// Function for calculating the scalar products between Col_amps.
/// Does not add implicit state in the gluons only case.
Polynomial scalar_product( const Col_amp & Ca1, const Col_amp & Ca2 ) const;

/// Function for calculating the scalar product between two Col_strs.
/// Does not add implicit state in the gluons only case.
Polynomial scalar_product( const Col_str & Cs1, const Col_str & Cs2 ) const;


/****************** Functions for gluon emission exchange, and splitting *************/

/// Function for emitting a gluon from a Col_str.
/// When the gluon is inserted before the emitter in a Quark_line,
/// the amplitude comes with a minus sign.
Col_amp emit_gluon( const Col_str & in_Col_str, int emitter, int g_new ) const;

/// Function for emitting a gluon from a Col_amp.
/// When the gluon is inserted before the emitter in a Quark_line,
/// the amplitude comes with a minus sign.
Col_amp emit_gluon( const Col_amp & Ca_in, int emitter, int g_new ) const;

/// Function for splitting the gluon g_old in a Col_str to a qqbar pair.
Col_amp split_gluon( const Col_str & in_Col_str, int g_old, int q_new, int qbar_new ) const;

/// Function for splitting the gluon g_old in a Col_amp to a qqbar pair.
Col_amp split_gluon( const Col_amp & in_Col_amp, int g_old, int q_new, int qbar_new ) const;


/// Function for exchanging a gluon between the partons p1 and p2 in the Col_str Cs.
/// When the gluon is inserted before the emitter in a Quark_line,
/// the amplitude comes with a minus sign.
Col_amp exchange_gluon( const Col_str & Cs, int p1, int p2 ) const;

/// Function for exchanging a gluon between two partons p1 and p2 in the Col_amp Ca.
/// When the gluon is inserted before the emitter in a Quark_line,
/// the amplitude comes with a minus sign.
/// (There is no special treatment of the glons only cases.)
Col_amp exchange_gluon( const Col_amp & Ca, int p1, int p2 ) const;

/// Calculates < M | T_i T_j | M >, the "color correlator"
/// relevant for coherent gluon emission of gluon g_new from
/// parton i and parton j, or gluon exchange between i and j.
/// The Ca should thus be | M >,
/// and i and j are the partons involved in the emission (exchange).
/// (The incoming amplitude is what it is, there is no special
/// treatment of gluons only cases.)
Polynomial color_correlator( const Col_amp Ca, int i, int j ) const;


/********************* Other functions *************************/

// As dvec and dmatr are not classes some read and write functions
// are contained here.

/// Reads in a numerical vector and save it as a double vector, dvec.
/// The file should be of the format
/// {d11,...,d1n},
/// and may contain comment lines starting with # at the top.
dvec read_in_dvec( std::string filename ) const;

/// Reads in a numerical matrix and save it as a double matrix, dmatr.
/// The file should be of the format
/// {{d11,...,d1n},
/// ...,
/// {dn1,...,dnn}},
/// and may contain comment lines starting with # at the top.
dmatr read_in_dmatr( std::string filename ) const;

/// Function for writing out a numerical vector,
/// to the file filename.
void write_out_dvec( const dvec & dv, std::string filename ) const;

/// Writes out the double version of a (scalar product) matrix
/// to the file filename.
void write_out_dmatr( const dmatr & matr, std::string filename ) const;

/// The factorial of an int, 0! is defined as 1.
int factorial( int i ) const;

/// Function that finds the default parton numbers for a Col_str.
/// The default numbers are 1,...,N_parton, where quarks have the first
/// odd numbers, anti-quarks have the first even numbers, and gluons have
/// subsequent numbers. The intended usage is after gluon splitting, where
/// the gluon g_old has split into q_new and qbar_new, and the parton numbers
/// before splitting are assumed to be default.
std::map<int, int> default_parton_numbers( const Col_str &, int g_old, int q_new, int qbar_new ) const;


/// Function that renames the partons in a Col_str using a map where, in each pair,
/// the first number is to be replaced by the second. (The Col_functions member
/// function default_parton_numbers returns a map where the partons are given
/// default numbers.)
Col_str rename_partons( const Col_str &, const std::map <int, int> replacements ) const;


/// Function that renames the partons in a Col_amp using a map where, in each pair,
/// the first number is to be replaced by the second. (The Col_functions member
/// function default_parton_numbers returns a map where the partons are given
/// default numbers.)
Col_amp rename_partons( const Col_amp &, const std::map <int, int> replacements ) const;


}; // end class Col_functions


/////////////////////// DECLEARING OPERATORS /////////////////////

/// Defining + operator to be able to erase elements at place
std::list<int>::iterator operator+( std::list<int>::iterator x, int n );

std::list<int>::iterator operator-( std::list<int>::iterator x, int n );

/// Defining + operator to be able to erase elements at place
std::list< Quark_line >::iterator operator+( std::list < Quark_line >::iterator x, int n );

col_str::iterator operator-( col_str::iterator x, int n );

/// Define the operator << for vector of int.
std::ostream& operator<<( std::ostream& out, const std::vector<int> & vec );

/// Define the operator << for cvec.
std::ostream& operator<<( std::ostream& out, const cvec & cv );

/// Define the operator << for dvec.
std::ostream& operator<<( std::ostream& out, const dvec & dv );

/// Define the operator << for cmatr.
std::ostream& operator<<( std::ostream& out, const cmatr & cm );

/// Define the operator << for dmatr.
std::ostream& operator<<( std::ostream& out, const dmatr & matr );

/// Define the operator << for std::pair<int, int>.
std::ostream& operator<<( std::ostream& out, std::pair<int, int> pair );
} // end namespace ColorFull

#endif /* COLORFULL_Col_functions_h */
