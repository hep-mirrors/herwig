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
#include <boost/shared_ptr.hpp>

namespace ColorFull {

using boost::shared_ptr;


/// Library class containing functions for index contraction and
/// numerical evaluation.
/// This is where the parameters Nc, TR and CF are contained.
class Col_functions {

private:

/// The number of colors, used in numerical results. Change here for different value.
double Nc;

/// The trace convention Tr( t^a t^a )=TR (no sum).
/// The normalization of the SU(Nc) generators, to be used in numerical evaluation.
/// The value 1/2 corresponds to the Gell-Mann normalization.
double TR;

/// The value of CF=TR (Nc^2-1)/(Nc).
/// Note that CF can be changed independently of Nc.
double CF;

/// While evaluating leading terms one may want to keep the full value of CF for
/// Nc=3, (TR(Nc^2-1)/Nc), or only keep the leading Nc term =TR*Nc (default).
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

/// Get the number of colors.
double get_Nc() const { return Nc; }

// Returns the normalization of the generators,
/// tr(t^a t^b)=TR*delta^{a,b}.
double get_TR() const { return TR; }

/// Get the number of colors.
double get_CF() const { return CF; }

/// Return true, if full CF is used
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

/// Function for finding the leading power of Nc in a
/// vector of pointer to Polynomials.
int  leading_Nc_pow( std::vector< boost::shared_ptr<Polynomial> > Pvp) const;

/// Takes the leading Nc terms of a Polynonmial, i.e. Monomials with highest
/// power of Nc+CF. If full_CF is false (default), CF is replaced by TR Nc.
/// If full_CF is true CF is replaced by TR(Nc^2-1)/Nc.
Polynomial leading( const Polynomial & Poly ) const;

/// Take the leading part of a Poly_vec.
/// Keeps only Monomials with maximal power of CF + Nc.
/// Uses leading( const Polynomial & Poly).
/// If full_CF is false (default), CF is replaced by TR Nc.
/// If full_CF is true CF is replaced by TR(Nc^2-1)/Nc.
/// Note that taking the leading terms of a Poly_vec is not
/// the same as taking the leading terms in each Polynomial.
// Used only by Poly_matr version of leading
Poly_vec leading( const Poly_vec & Pv ) const;

/// Takes the leading part of a matrix of Polynomials,
/// keeping only those with maximal power of CF + Nc.
/// If full_CF is false (default), CF is replaced by TR Nc.
/// If full_CF is true CF is replaced by TR(Nc^2-1)/Nc.
/// Note that taking the leading terms of a Poly_matr is not
/// the same as taking the leading terms in each Poly_vec.
// Used only once in Col_basis
Poly_matr leading( const Poly_matr & Pm ) const;

/// Take the leading part of a Poly_vec, given a vector of pointers to the Polynomials.
/// Keeps only Monomials with maximal power of CF + Nc.
// Currently never used
Poly_vec leading( std::vector<boost::shared_ptr<Polynomial> > Pvp) const;

/// Take the leading part of a Poly_matr, given a vector of vector of pointers to the Polynomials.
/// Loops over Monomials in all Polynomials
/// and keeps only those with maximal power of CF + Nc.
// used only by scalar_product_matrix_mem in Col_functions
dmatr leading( std::vector< std::vector< boost::shared_ptr<Polynomial> > > Pm ) const;

/// To keep only leading terms in a map.
// used only on Col_functions by scalar product_matrix_mem_2 and radiation_amplitude_matrix
std::map< std::string, Polynomial > leading( std::map< std::string, Polynomial > mem_map ) const;


/********************* Functions for numerical evaluation *************************/
// These functions has to be kept in Col_functions class as they need numerical
// values for evaluation. Letting each Polynomial carry around its own Nc etc.
// would be messy.

/// To take the numerical value of a map.
std::map< std::string, double > double_num( std::map< std::string, boost::shared_ptr<Polynomial> > mem_map ) const;

/// Numerically evaluates a Monomial using the Nc and CF variables;
cnum cnum_num( const Monomial & Mon ) const;

// Numerically evaluates a Monomial to a double.
double double_num( const Monomial & Mon ) const;

/// Numerically evaluates a Polynomial, using the CF and Nc variables.
cnum cnum_num( const Polynomial & Poly ) const;

/// Numerically evaluates a Polynomial using the data members Nc, CF and TR.
double double_num( const Polynomial & Poly ) const;

/// To take the numerical value of a map.
std::map< std::string, double > double_num(std::map< std::string, Polynomial > mem_map) const;

/// Numerically evaluates a Poly_vec (vector of Polynomial) for Nc=3.
dvec double_num( const Poly_vec & Pv ) const;

/// Returns a double value. The argument is a vector of pointers to Polynomials.
dvec double_num( const std::vector<boost::shared_ptr<Polynomial> > & Pv ) const;

/// Returns a double value. The argument is a vector of vector of pointers to Polynomials.
dmatr double_num( const std::vector<std::vector<boost::shared_ptr<Polynomial> > > & Pm ) const;

/// Numerically evaluates a Polynomial for Nc=3,
/// and stores in the format of a Polynomial with only one term with only a numerical part.
Polynomial Polynomial_cnum_num( const Polynomial & Poly ) const;

/// Numerically evaluates a Poly_vec (vector of Polynomial),
/// using cnum_num (Polynomial).
cvec cnum_num( const Poly_vec & Pv ) const;

/// Numerically evaluates a Poly_vec (vector of Polynomial)
/// and stores in the form of a Poly_vec, uses polynomial_cnum_num( Pv.at( p ) ).
/// for each Polynomial.
Poly_vec Poly_vec_cnum_num( const Poly_vec & Pv ) const;

/// Numerically evaluates a Poly_matr (vector of Poly_vec)
/// and stores in the form of a Poly_matr.
Poly_matr Poly_matr_cnum_num( const Poly_matr & Pm ) const;

/// Numerically evaluates a Poly_matr (vector of Poly_vec),
/// using cnum_num( Poly_vec ) for each Poly_vec.
cmatr cnum_num( const Poly_matr & Pm ) const;

/// Numerically evaluates a Poly_matr (vector of Poly_vec).
dmatr double_num( const Poly_matr & Pm ) const;

/****************** Functions for scalar products *************************/

/// Function for calculating scalar products between Col_amps.
/// Does not add implicit state in the gluons only case.
Polynomial scalar_product( const Col_amp & Ca1, const Col_amp & Ca2 ) const;

/// Function for calculating scalar product between two Col_strs.
/// Does not add implicit state in the gluons only case.
Polynomial scalar_product( const Col_str & Cs1, const Col_str & Cs2 ) const;


/****************** Functions for gluon emission and exchange *************/

/// Function for emitting a gluon from a Col_str.
/// When the gluon is inserted before the emitter in a Quark_line,
/// the amplitude comes with a minus sign.
Col_amp emit_gluon( const Col_str & in_Col_str, int emitter, int g_new ) const;

/// Function for emitting a gluon from a Col_amp.
/// When the gluon is inserted before the emitter in a Quark_line,
/// the amplitude comes with a minus sign.
Col_amp emit_gluon( const Col_amp & Ca_in, int emitter, int g_new ) const;

/// Function for exchanging a gluon between the partons p1 and p2 in the Col_str Cs.
/// When the gluon is inserted before the emitter in a Quark_line,
/// the amplitude comes with a minus sign.
Col_amp exchange_gluon( const Col_str & Cs, int p1, int p2 ) const;

/// Function for exchanging a gluon between two partons p1 and p2 in the Col_amp Ca.
/// When the gluon is inserted before the emitter in a Quark_line,
/// the amplitude comes with a minus sign.
/// (The incoming amplitude is what it is, there is no special
/// treatment of glons only cases.)
Col_amp exchange_gluon( const Col_amp & Ca, int p1, int p2 ) const;

/// Calculates < M | T^(i) T^(j) | M >, the "color correlator"
/// relevant for coherent gluon emission of gluon g_new from
/// parton i and parton j, or gluon exchange between i and j.
/// The Ca should thus be | M >, g_new should be a unique dummy index,
/// and i and j are the partons involved in the emission (exchange).
/// (The incoming amplitude is what it is, there is no special
/// treatment of gluons only cases.)
Polynomial color_correlator( const Col_amp Ca, int i, int j, int g_new ) const;


/********************* Other functions *************************/

// As dvec and dmatr are not classes some read and write functions
// are contained here.

/// Read in a numerical matrix from filename and save it as a double matrix, dmatr.
/// The file should be in the format
/// {d11,...,d1n}.
dvec read_in_dvec( std::string filename ) const;

/// Read in a numerical matrix from filename and save it as a double matrix, dmatr.
/// The file should be in the format
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

}; // end class Col_functions


/////////////////////// DECLEARING OPERATORS /////////////////////

/// Defining + operator to be able to erase elements at place
std::list<int>::iterator operator+( std::list<int>::iterator x, int n );

std::list<int>::iterator operator-( std::list<int>::iterator x, int n );

/// Defining + operator to be able to erase elements at place
std::list< Quark_line >::iterator operator+( std::list < Quark_line >::iterator x, int n );

col_str::iterator operator-( col_str::iterator x, int n );

/// Define the operator << for vector of int.
std::ostream& operator<<( std::ostream& out, std::vector<int> vec );

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
