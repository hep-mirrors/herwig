// -*- C++ -*-
/*
 * Monomial.h
 *	Contains declaration of the class Monomial and associated types and operators
 *  Created on: Jul 7, 2010
 *  Author: Malin Sjodahl
 */

#ifndef COLORFULL_Monomial_h
#define COLORFULL_Monomial_h

#include "types.h"


namespace ColorFull {


/// A class to contain the factor of form TR^a*Nc^b*CF^c*int_part*cnum_part,
/// where the powers a, b and c may be negative.
/// A default Monomial is defined to be 1, and has int_part and cnum_part=1.
/// A 0-Monomial has int_part=0.
/// A polynomial is a sum of Monomials.
class Monomial {
public:

	/// Default constructor sets int_part and cnum_part=1, and Pow_Nc=pow_TR=pow_CF=0.
	Monomial(){
		pow_TR=pow_Nc=pow_CF=0;
		int_part=1;
		cnum_part=1.0;
	}

	/// Constructor using a double.
	/// The cnum_part member is set to contain the value.
	Monomial( double dnum ){
	    pow_TR=pow_Nc=pow_CF=0;
	    int_part=1;
		cnum_part.real(dnum);
	}

	/// Constructor using an int.
	/// The int_part member is set to contain the value.
	Monomial( int num ){
	    pow_TR=pow_Nc=pow_CF=0;
	    int_part=num;
	    cnum_part=1.0;
	}

	/// Constructor taking a string as argument.
	/// The argument should be of the form given in form (for example)
	/// -(20*TR^5)/Nc or -20 TR^(5)/Nc or 20 / TR^(-5)Nc^(1) CF^(3).
	/// NOTE: All spaces and * are ignored, except in  "*(-1)" and *-1, which
	/// is understood as (*-1).
	/// EVERYTHING standing after / is
	/// divided with, whereas everything standing before is multiplied with.
	/// Parentheses are ignored unless they appear in powers,
	/// i.e, directly after ^.
	/// No spaces are allowed inside the powers.
	/// If the string contains no info or is empty the Monomial is put to 1,
	/// pow_TR = pow_Nc = pow_CF = 0, int_part = 1, cnum_part = 1.0.
	/// (Expanded Mathematica 8 expressions are in this form.)
	Monomial( std::string str );

	/// Power of TR in Monomial.
	int pow_TR;

	/// Power of the number of colors.
	int pow_Nc;

	/// Power of CF=TR (Nc^2-1)/(Nc)
	int pow_CF;

	/// Integer multiplying monomial, can be 0.
	int int_part;

	/// Complex number multiplying monomial.
	cnum cnum_part;

	/// Take the complex conjugate.
	/// Note that this changes the Monomial itself.
	void conjugate() {
		cnum_part=conj( cnum_part );
	}

	/// Function for reading in the Monomial from the file filename,
	/// uses Monomial_of_str.
	void read_in_Monomial( std::string filename );


private:

	/// Function for makinga a Monomial from a string.
	/// The argument should be of the form given in form (for example)
	/// -(20*TR^5)/Nc or -20 TR^(5)/Nc or 20 / TR^(-5)Nc^(1) CF^3.
	/// NOTE: All spaces and * are ignored, except in  "*(-1)" and *-1, which
	/// is understood as (*-1).
	/// EVERYTHING standing after / is
	/// divided with, whereas everything standing before is multiplied with.
	/// Parentheses are ignored unless they appear in powers,
	/// i.e, directly after ^.
	/// No spaces are allowed inside the powers.
	/// If the string contains no info or is empty the Monomial is put to 1,
	/// pow_TR = pow_Nc = pow_CF = 0, int_part = 1, cnum_part = 1.0.
	void Monomial_of_str( std::string str );


};


/// Define the operator << for Monomial
std::ostream& operator<<( std::ostream& out, const Monomial & Mon );

/// Operator * for Monomial and int.
/// The int_part member is multiplied by i, whereas other
/// members are kept constant.
Monomial operator*( const Monomial & Mon,  const int i );

/// Operator * for int and Monomial.
/// Returns Mon*i.
Monomial operator*( const int i, const Monomial & Mon );

/// Operator * for Monomial and cnum. The member Mon.cnum_part
/// is multiplied by c, whereas other members are kept
/// the same.
Monomial operator*( const Monomial & Mon, const cnum c );

/// Operator * for cnum and Monomial, returns Mon*c.
Monomial operator*( const cnum c, const Monomial & Mon );

/// Operator * for Monomial and double. The member Mon.cnum_part
/// is multiplied by c, whereas other members are kept
/// the same.
Monomial operator*( const Monomial & Mon, const double d );

/// Operator * for double and Monomial, returns Mon*d.
Monomial operator*( const double d, const Monomial & Mon );

/// Operator * for Monomials. The powers, pow_TR, pow_Nc
/// and pow_CF are added, and the numbers int_part and mon
/// are multiplied.
Monomial operator*( const Monomial & Mon1, const Monomial & Mon2 );

/// Operator == for Monomials, all parts must be equal.
/// For the numerical part, an accuracy is used for the ratio.
bool operator==( const Monomial & Mon1, const Monomial & Mon2 );

/// Define the operator != Monomial. Returns false if Mon1==Mon2,
/// and true otherwise.
bool operator!=( const Monomial & Mon1, const Monomial & Mon2 );

/// Operator to find the "smallest" of two Monomials.
/// The Monomials are ordered first according to pow_Nc+pow_CF,
/// then according to pow_Nc (for same pow_Nc+pow_CF)
/// then according to int_part*abs(cnum_part), then according to int_part, and finally according to pow_TR.
/// NOTE: this ordering does not agree with ordering according to numerical value.
/// NOTE: If the Monomials are equal, Mon1 is not smaller so false will be returned.
bool operator<( const Monomial & Mon1, const Monomial & Mon2 );

}

#endif /* COLORFULL_Monomial_h */
