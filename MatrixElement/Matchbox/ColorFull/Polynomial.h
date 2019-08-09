// -*- C++ -*-
/*
 * Polynomial.h
 *	Contains declaration of the class Polynomial and associated types and operators.
 *  Created on: Jul 7, 2010
 *  Author: Malin Sjodahl
 */

#ifndef COLORFULL_Polynomial_h
#define COLORFULL_Polynomial_h

#include "Monomial.h"

namespace ColorFull {

/// For containing the info (as opposed to the functions) of a Polynomial.
/// The Polynomial is essentially a sum of Monomials, contained in a
/// vector of Monomials.
/// Note that a Monomial "is" TR^a Nc^b CF^c* int_part*cnum_part where int_part is an
/// integer factor and cnum_part a complex numerical factor.
/// An empty polynomial is defined as 1.
typedef std::vector < Monomial > polynomial;


/// For containing a Polynomial (in Nc, CF and TR), as a sum of Monomials.
/// Note that a Monomial "is" TR^a*Nc^b*CF^c*int_part*cnum_part where int_part is an
/// integer factor, cnum_part a complex numerical factor and a,b, and c integers
/// (not necessarily positive).
/// An empty Polynomial is defined as 1.
class Polynomial {
public:

	/// Default constructor, leaves polynomial empty=1.
	Polynomial(){};

	/// Constructor allow setting the polynomial by using a string.
	/// Should be used as for example "Polynomial Poly("(-20*TR^(5))/Nc + 28*Nc*TR^(5) - 10*Nc^3*TR^(5)")".
	/// The Momomials should be separated by + or -, see also the
	/// Monomial string constructor.
	Polynomial( const std::string str );

	/// Constructor allowing setting the Polynomial using a double.
	/// The Polynomial gets one Monomial where the real part of
	/// cnum_part gets the value of dnum.
	Polynomial( double dnum );

	/// Constructor allowing setting the Polynomial using an int.
	/// The Polynomial gets one Monomial where int_part
	/// has the value num.
	Polynomial( int num );

	/// Contains the polynomial, a sum of Monomials, as an std::vector
	/// of Monomials.
	/// An empty Polynomial is defined as 1, to get 0, multiply with 0.
	polynomial poly;

	/// Function for reading in the Polynomial from the file filename,
	/// uses Polynomial_of_str.
	void read_in_Polynomial( std::string filename );

	/// Function for writing out the Polynomial to a file
	/// with name filename.
	void write_out_Polynomial( std::string filename ) const;

	/// Returns Monomial at place i.
	const Monomial& at( int i ) const {return poly.at(i);}

	/// Returns Monomial at place i.
	Monomial& at( int i ) {return poly.at(i);}

	/// Returns the number of terms in the Polynomial.
	int size() const {return poly.size();}

	/// Adding a Monomial term.
	void append( const Monomial Mon ) {poly.push_back(Mon);}

	/// Erases the Monomial at place i.
	void erase( int i ) {poly.erase(poly.begin() + i);}

	/// Is the polynomial empty?
	bool empty() const {return poly.empty();}

	/// Erases info in polynomial.
	void clear() {poly.clear();}

	/// Take complex conjugate of the polynomial.
	/// Note that this changes the Polynomial itself.
	void conjugate();

	/// Collects terms with same power of TR and Nc and Cf.
	void simplify();

	/// Replaces CF with TR*Nc-TR/Nc.
	void remove_CF();

	/// Orders terms in Polynomial in a unique form,
	/// first according to pow_Nc+pow_CF, then according to pow_Nc (for same pow_Nc+pow_CF)
	/// then according to int_part*num, then according to int_part, and finally according to pow_TR.
	void normal_order();

private:
	/// The factorial of an int.
	int factorial (int i) const;

	/// Function for setting the Polynomial using a string,
	/// used by string constructor.
	/// The Momomials should be separated by + or -, see also the
	/// Monomial string constructor.
	void Polynomial_of_str( const std::string str );

};

/// Define the operator << for Polynomial.
std::ostream& operator<<(std::ostream& out, const Polynomial & Poly);

/// Define the operator << for polynomial.
std::ostream& operator<<(std::ostream& out, const polynomial & poly);

/// Define the operator == for Polynomial
/// By definition, each Monomial has to be identical,
/// as order matters 1+2 is not 2+1.
bool operator==( const Polynomial & Poly1, const Polynomial & Poly2);

/// Operator != for Polynomial. Returns false if Poly1==Poly2, and
/// true otherwise.
bool operator!=( const Polynomial & Poly1, const Polynomial & Poly2 );

/// Operator + for Polynomial and a single Monomial.
Polynomial operator+(const Polynomial & Poly, const Monomial & Mon );

/// Operator + for a single Monomial and a Polynomial,
/// returns Poly+Mon.
Polynomial operator+(const Monomial & Mon, const Polynomial & Poly );

/// Define the operator += for Polynomials.
/// The Monomial is appended unless Mon.int_part=0.
Polynomial operator+=( Polynomial & Poly, const Monomial & Mon );

/// Define the operator += for Polynomials.
/// This operator appends the Monomials of Poly2 to Poly1.
Polynomial operator+=( Polynomial & Poly1, const Polynomial & Poly2 );

/// Operator - for Polynomial and a single Monomial.
/// Changes sign of Monomial by changing int_part and returns Poly+(-Mon).
Polynomial operator-(const Polynomial & Poly, const Monomial & Mon);

/// Operator - for Polynomial and a single Monomial
/// changes the sign of Poly by changing int_part,
/// and then adds Mon to (-Poly).
Polynomial operator-(const Monomial & Mon, const Polynomial & Poly);

/// Operator + for Polynomials, appends Momomails from Poly2 to Poly1.
/// If one of the Polynomials is empty=1, this is compensated for by first
/// appending a Monomial=1 to the empty Polynomial.
Polynomial operator+( const Polynomial & Poly1, const Polynomial & Poly2 );

/// Operator - for Polynomials,
/// changes sign of Poly2, and then returns Poly1 + (-Poly2).
Polynomial operator-(const Polynomial & Poly1, const Polynomial & Poly2);

/// Operator * for Polynomial and int,
/// loops over Monomials and multiplies int_part with i.
/// If Poly is empty=1, a default Monomial=1 is first
/// appended to Poly.
Polynomial operator*( const Polynomial & Poly, int i );

/// Operator * for int and Polynomial.
/// Returns Poly*i.
Polynomial operator*( int i, const Polynomial & Poly );

/// Operator * for Polynomial and cnum,
/// loops over Monomials and multiplies cnum_partwith c.
/// If Poly is empty=1, a default Monomial=1 is first
/// appended to poly.
Polynomial operator*( const Polynomial & Poly, const cnum c );

/// Operator * for cnum and Polynomial, returns Poly*c.
Polynomial operator*( const cnum c, const Polynomial & Poly );

/// Operator * for Polynomial and double,
/// loops over Monomials and multiplies cnum_partwith d.
/// If Poly is empty=1, a default Monomial=1 is first
/// appended to poly.
Polynomial operator*( const Polynomial & Poly, const double d );

/// Operator * for Polynomial and double, returns Poly*d.
Polynomial operator*( const double d, const Polynomial & Poly );

/// Operator * for Polynomial and Monomial, loops over Monomials in
/// Poly, and multiplies each with Mon. If Poly is empty=1,
/// a Polynomial with one Monomial=Mon is returned.
Polynomial operator*( const Polynomial & Poly, const Monomial & Mon );

/// Operator * for Monomial and Polynomial,
/// returns Poly*Mon.
Polynomial operator*( const Monomial & Mon, const Polynomial & Poly );

/// Operator * for Polynomials. (For details, see code.)
Polynomial operator*( const Polynomial & Poly1, const Polynomial & Poly2 );

/// Operator *= for Polynomial and int.
/// Multiplies Poly with i.
Polynomial operator*=( Polynomial & Poly, const int i );

/// Operator *= for Polynomial and double.
/// Multiplies Poly with d.
Polynomial operator*=( Polynomial & Poly, const double d );

/// Operator *= for Polynomial and cnum.
/// Multiplies Poly with c.
Polynomial operator*=( Polynomial & Poly, const cnum c );

/// Operator *= for Polynomial and Monomial.
/// Multiplies Poly with Mon.
Polynomial operator*=( Polynomial & Poly, const Monomial & Mon );

/// Operator *= for Polynomial and Polynomial.
/// Multiplies Poly1 with Poly2.
Polynomial operator*=( Polynomial & Poly1, const Polynomial & Poly2 );

}


#endif /* COLORFULL_Polynomial_h */
