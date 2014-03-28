// -*- C++ -*-
/* Polynomial.cc
 *	Contains definition of the class Polynomial and associated types and operators
 *  Created on: Jul 7, 2010
 *  Author: Malin Sjodahl
 */
#include "Polynomial.h"
#include "parameters.h"
#include <limits>
#include <cassert>
#include <fstream>
#include <iostream>


namespace ColorFull {

Polynomial::Polynomial( double dnum ) {

	Monomial Mon;
	Mon.cnum_part.real(dnum);
	poly.push_back(Mon);

}


Polynomial::Polynomial( int num ) {
	Monomial Mon;
	Mon.int_part= num ;
	poly.push_back(Mon);
}


Polynomial::Polynomial( const std::string str ) {
	Polynomial_of_str( str );
}


void Polynomial::Polynomial_of_str( const std::string str ) {

	// Special case of empty string
	if( str.empty() ) {
		polynomial empty;
		poly=empty;
		return;
	}

	uint i = 0, j=0;
	std::string Mon_str;

	int left_brackets=0, right_brackets=0;
	// Check that left and right brackets match up
	while (j < str.size()) {
		if(str.at(j)=='(') left_brackets++;
		if(str.at(j)==')') right_brackets++;
		j++;
	}
	if(left_brackets != right_brackets){
		std::cerr << "Polynomial::Polynomial_of_str: The brackets in the polynomial\"" << str <<"\" do not seem to match up. There were "
				<< left_brackets <<" left bracket(s) and "<< right_brackets << " right bracket(s)." << std::endl;
		assert( 0 );
	}

	// Skip potential spaces
	while (i < str.size() && (str.at(i) == ' ')) i++;

	// Skip initial left bracket(s) surrounding polynomial
	// there are left brackets-right brackets such brackets,
	// we many have ((TR+Nc))
	int found_left_brackets=0;

	// if we have a surplus of left_brackets before we have a + or -
	// then (with some exceptions) it's an overall bracket
	int right_before_plus=0, left_before_plus=0;
	j=0;
	bool new_term=false;
	while ( j<str.size() and !new_term ){

		// we should allow brackets like (2TR)+4
		// but remove overall brackets, like (2 TR+4)
		// a + or - can means that we may have reached a new term
		// and should stop looking for overall brackets
		if(str.at(j)=='-'){
			// we can allow - inside powers, i.e. after ^(, ^(- or ^-
			if( j>0 ){
				// if (- continue, (this is not a new term)
				if(str.at(j-1)=='('){}
				// if ^- continue, (this is not a new term)
				else if(str.at(j-1)=='^'){}
				// if *- continue, (this is not a new term)
				else if(str.at(j-1)=='*'){}
				else new_term=true;
			}
			else new_term=true;
		}

		// a plus sign could mean that it's a new term
		if(str.at(j)=='+'){
			// we can allow - inside powers, i.e. after ^(, ^(- or ^-
			if( j>0 ){
				// if (+ continue, (this is not a new term)
				if(str.at(j-1)=='('){}
				// if ^+ continue, (this is not a new term)
				else if(str.at(j-1)=='^'){}
				// if *+ continue, (this is not a new term)
				else if(str.at(j-1)=='*'){}
				else new_term=true;
			}
			else new_term=true;
		}

		if(str.at(j)== '(') left_before_plus++;
		if(str.at(j)== ')') right_before_plus++;
		j++;
	}

	int extra_left_brackets =  left_before_plus-right_before_plus;

	while (i < str.size() and ( str.at(i)=='(' or str.at(i)==' ')  and found_left_brackets < extra_left_brackets ) {
		if(str.at(i)== '(') found_left_brackets++;
		i++;
	}

	int Moncount=0;

	// Look for Monomials until the end of the string
	while (i < str.size()) {
		Mon_str.clear();
		Moncount++;

		// Check for one + or - for each Monomial
		while (i < str.size() && (str.at(i) == '+' or str.at(i) == '-' )){
			Mon_str.push_back(str.at(i));
			i++;
		}

		// Look for more factors in same Monomial
		while (i < str.size()) {

			// There should be no surplus of ) in the Monomial
			left_brackets=right_brackets=0;

			// Check for Monomials
			// Everything but + and - should go to Monomail
			// Read until another + or - while making sure we don't get final )'s
			if ( i < str.size() ){
				if( str.at(i) != '+' and str.at(i) != '-' ){
					Mon_str.push_back(str.at(i));
				i++;
				}
			}

			// we can allow - after (, "(-"for example inside powers
			if ( i < str.size() and  str.at(i) == '-'){
				// removed i-1>= 0 as always true, no check that read inside
				if(( str.at(i-1)=='(')){
					Mon_str.push_back(str.at(i));
					i++;
				}
			}

			// We can also allow for - after*, for example 0.5*-1
			if (i < str.size() and  str.at(i) == '-' and ( (i>= 2) and str.at(i-1)=='*'))  {
				// read in sign
				Mon_str.push_back(str.at(i));
				i++; // skip sign
				// Keep reading in while numbers
				while (i< str.size() and  (str.at(i) == '0' or str.at(i) == '1' or str.at(i) == '2'
						or str.at(i) == '3' or str.at(i) == '4'
								or str.at(i) == '5' or str.at(i) == '6'
										or str.at(i) == '7' or str.at(i) == '8'
												or str.at(i) == '9')){
					Mon_str.push_back(str.at(i));
					i++;
				}
			}

			// We can also allow for - after ( as in (-1)
			if (i < str.size() and str.at(i) == '-' and (i>= 3) and str.at(i-1)=='(') {
				// read in sign
				Mon_str.push_back(str.at(i));
				i++; // skip sign

				// Keep reading in while numbers
				while (i< str.size() and  (str.at(i) == '0' or str.at(i) == '1' or str.at(i) == '2'
						or str.at(i) == '3' or str.at(i) == '4'
								or str.at(i) == '5' or str.at(i) == '6'
										or str.at(i) == '7' or str.at(i) == '8'
												or str.at(i) == '9')){
					Mon_str.push_back(str.at(i));
					i++;
				}

				// Skip spaces
				while ( i< str.size() and str.at(i) == ' ') i++;

				// Look for closing )
				if(str.at(i)==')') {
					// read in )
					Mon_str.push_back(str.at(i));
					i++;}
				else {
					std::cerr << "Polynomial::Polynomial_of_str: Expects final ) in *(number)"
							<< " got " << str.at(i) << std::endl;
					assert( 0 );
				}
			}// Stop looking for (-

			// Break if we have + or - and not a special case
			if ( i == str.size() ) break;
			if ( str.at(i)== '+' ) break;
			if (str.at(i) == '-'
					and !(
					(i < str.size() and  str.at(i) == '-'
							and ( ( i-3 > 0 and str.at(i-1)=='(' and str.at(i-2)=='^')
							or ( i-2> 0 and  str.at(i-1)=='^' )	))
					// We can also allow for - after*, for example 0.5*-1
					or (i < str.size() and  str.at(i) == '-' and ( (i-2> 0) and str.at(i-1)=='*'))
					// We can also allow for - after 0.5*(-1)
					or (i < str.size() and str.at(i) == '-' and (i-3> 0) and str.at(i-1)=='(' and str.at(i-2)=='*')) ) break;
		} // end of while (i < str.size()), i.e. now keep looking for factors in same Monomial
		// Then make a Monomial out of the str


		// Check that left and right brackets match up
		left_brackets=right_brackets=0;
		for ( uint j=0; j< Mon_str.size(); j++ ){
			if(Mon_str.at(j)=='(') left_brackets++;
			if(Mon_str.at(j)==')') right_brackets++;
		}

		int brack_diff=right_brackets-left_brackets;
		// Remove (-brackets in the end until they match ), i.e. brack_diff==0.
		uint j=1;
		while( ( Mon_str.at( Mon_str.size()-j )==')' or Mon_str.at( Mon_str.size()-j )==' ' or Mon_str.at( Mon_str.size()-j )=='*') and brack_diff > 0 ){
			if( Mon_str.at( Mon_str.size()-j )==')' ) {
				brack_diff--;
				Mon_str.at( Mon_str.size()-j )= ' ';
			}
			j++;
		}

		poly.push_back(Monomial(Mon_str));

	} // end of while (i < str.size())

}


void Polynomial::read_in_Polynomial( std::string filename) {

	// Read in file
	std::ifstream fin(filename.c_str());

	// Check that file exists
	if( !fin ){
		std::cerr << "Polynomial::read_in_Polynomial: The file "
				<< filename << " could not be opened." << std::endl;
		assert( 0 );
	}

	// Erase current information
	poly.clear();

	// Copy info from file to string
	std::string str((std::istreambuf_iterator<char>(fin)), std::istreambuf_iterator<char>());
	//str.erase(str.size()-1);
	Polynomial_of_str( str );
}


void Polynomial::write_out_Polynomial( std::string filename ) const {

	std::ofstream outfile(filename.c_str());

	if ( !outfile )
	std::cerr << "Polynomial::write_out_Polynomial: Cannot write out diagonal scalar products as the file \""
		<< filename.c_str() << "\" could not be opened. (Does the directory exist? Consider creating the directory.)" << std::endl;

	outfile << *this;
}


void Polynomial::conjugate(){
	// Loop over Monomials and complex conjugate them
	for(uint i=0; i< poly.size(); i++) poly.at(i).conjugate();

}


void Polynomial::simplify() {

	// If empty Polynomial or only one Monomial, do nothing
	if (poly.size() <= 1)
		return;

	// Collect different Monomials, keep 0th Monomial
	polynomial terms;
	terms.push_back(poly.at(0));
	poly.erase(poly.begin());

	// Move terms from Polynomial to unique set of terms
	// as long as terms remain in Polynomial
	while (poly.size() > 0) {
		bool was_found = false;

		// Check if term already there
		for (uint i = 0; (i < terms.size() && !was_found); i++) {
			// If same powers of TR, Nc and CF
			if (poly.at(0).pow_TR == terms.at(i).pow_TR && poly.at(0).pow_Nc
					== terms.at(i).pow_Nc && poly.at(0).pow_CF
					== terms.at(i).pow_CF) {
				was_found = true;

				// Change the int_part (if no numerical factor)
				// or the num factor of term to be added

				// If already stored term, or part to be added, has a numeric part,
				// store info in numeric part
				if ( abs(poly.at(0).cnum_part-1.0) > accuracy or abs( terms.at(i).cnum_part -1.0 ) > accuracy) {

					// Float to add to num, put numerical factor, and int_part,
					// and sign here
					cnum c_to_add;
					c_to_add = poly.at(0).cnum_part;
					c_to_add *= poly.at(0).int_part;

					// Move int-part of existing to numeric part
					terms.at(i).cnum_part *= terms.at(i).int_part;
					terms.at(i).int_part = 1;

					// Add cnum from the term
					terms.at(i).cnum_part += c_to_add;
				} else { //if no numeric part change the int-part instead
					terms.at(i).int_part += poly.at(0).int_part;
				}
			}
			// If, after this the i:th term is 0, remove it,
			// to avoid carrying around 0's
			if( ( terms.at(i).int_part==0 or terms.at(i).cnum_part==0.0 )  and terms.size()>1 ){
				terms.erase(terms.begin()+i);
			}

		}
		// If term not already there, add term
		if ( !was_found and poly.at(0).int_part!=0 )
			terms.push_back(poly.at(0));

		// Erase the term in poly
		poly.erase(poly.begin());
	}

	for (uint i = 0; i < terms.size() ; i++) {
		// If term not already there, add term
		if ((terms.at(i).int_part == 0 or terms.at(i).cnum_part == 0.0)
				and terms.size() > 1) {
			terms.erase(terms.begin() + i);
		}
	}
	poly = terms;
}


void Polynomial::remove_CF() {

	// If empty Polynomial, do nothing
	if (poly.size() == 0)
		return;
	else{//non-empty polynomial
		polynomial poly_res; // to contain the result
		// Loop over terms in Polynomial and replace CF
		for (uint i = 0; i < poly.size(); i++) {

			// If the term contains no CF, simply keep in terms
			if( poly.at(i).pow_CF == 0 )
				poly_res.push_back( poly.at(i) );

			// replace CF if raced to positive power
			// (unless the int part is 0, in which case nothing should be added)
			else if( poly.at(i).pow_CF > 0 and poly.at(i).int_part != 0){

				// Factors given by binomial theorem
				// one term for each power (Nc)^power_Nc (-1/Nc)^(pow_CF-power_Nc)
				// i.e. TR ( (Nc)^power_Nc (-1/Nc)^(pow_CF-power_Nc) )
				// =TR Nc^(power_Nc-(pow_CF-power_Nc))==TR Nc^(2*power_Nc-pow_CF)
				Monomial const_fact=poly.at(i); // To contain the factor that was multiplying CF
				const_fact.pow_CF=0; // replace CF
				const_fact.pow_TR+=poly.at(i).pow_CF; // we always get a factor TR
				for ( int power_Nc = 0; power_Nc <= poly.at(i).pow_CF; power_Nc++) {

					Monomial term; // one term in (a+b)^n

					// The corresponding binomial coefficient
					int bin_coef=factorial( poly.at(i).pow_CF) /factorial( poly.at(i).pow_CF - power_Nc )/factorial( power_Nc );

					term.pow_Nc=2*power_Nc-poly.at(i).pow_CF;
					term=bin_coef*term*int(pow( -1.0, int( poly.at(i).pow_CF-power_Nc) ));// fix sign

					// if the term is not 0, add it
					if( term.int_part!=0 && const_fact.int_part!=0 )
						poly_res.push_back( term*const_fact );
				}// end looping over binomial factors
			}// end if pow_CF>0


			else if ( poly.at(i).pow_CF < 0 ) {
				std::cerr << "Polynomial::remove_CF(): Warning: cannot replace negative powers of CF. Leaving Monomial term << "
						<< poly.at(i) << " as it is." << std::endl;
				poly_res.push_back( poly.at(i) );
			}

		}// end looping over Monomials in poly

		// If, after this, poly_res is still empty, this is because nothing was added
		// and the polynomial should be 0 (can for example happen for 0*CF)
		if ( poly_res.empty() ) poly_res.push_back(Monomial("0"));
		poly=poly_res;
	}
}


void Polynomial::normal_order() {

	// To contain the Polynomial in order
	polynomial poly_ordered;

	// Order the different Monomials in the Polynomial.
	// Do this by moving the Monomials one by one to poly_ordered
	while (poly.size() > 0) {
		// Next monomial to put in place (the last Monomial in poly)
		Monomial Mon_next = poly.at( poly.size() - 1 );

		// Then insert the Mon_next among the ordered Monomials
		// Count how many steps left Mon_next should be moved in poly_ordered
		uint steps_left = 0;
		while ( (steps_left < (poly_ordered.size())) && (!(Mon_next < poly_ordered.at(poly_ordered.size()-1-steps_left ))) ==1) {
			steps_left++;
		}

		// Insert the Monomial in the right place among the ordered Col_strs
		polynomial::iterator it=poly_ordered.end()-steps_left;
		poly_ordered.insert( it, Mon_next);

		// Erase the Monomial from poly
		poly.erase(  poly.end()-1 ) ;
	}
	poly=poly_ordered;
}


int Polynomial::factorial( int i ) const{
	if(i<0) {
		std::cerr << "Polynomial::factorial: intended for int >=0, was " << i << std::endl;
		std::cerr.flush();
		assert( 0 );
	}
	if (i==0) return 1;
	return factorial(i-1)*i; // Recursive call
}


std::ostream& operator<<( std::ostream& out, const Polynomial & Poly ) {

	if (Poly.size() == 0)
		out << "1";
	else if (Poly.size() == 1) out << Poly.at(0);
	else {
		if (Poly.size() > 1)
			out << "(";

		for (int i = 0; i < Poly.size(); i++) {
			out << Poly.at(i);
			if (i != Poly.size() - 1)
				out << " + ";
		}

		if (Poly.size() > 1)
			out << ")";
	}
	return out;
}


std::ostream& operator<<(std::ostream& out, const polynomial & poly) {

	if (poly.size() == 0)
		out << "1";
	else if (poly.size() == 1) out << poly.at(0);
	else {
		if (poly.size() > 1)
			out << "(";
		for ( uint i = 0; i < poly.size(); i++ ) {
			out << poly.at(i);
			if (i != poly.size() - 1)
				out << " + ";
		}
		if (poly.size() > 1)
			out << ")";
	}
	return out;
}


bool operator==( const Polynomial & Poly1, const Polynomial & Poly2){

	if( Poly1.size() != Poly2.size() ) return false;

	// All terms should be the same in order for Polynomials to be the same
	for (int i=0; i < Poly1.size(); i++ ){
		if( Poly1.at(i) != Poly2.at(i ) ) return false;
	}

	return true;
}


bool operator!=( const Polynomial & Poly1, const Polynomial & Poly2){
	if( Poly1==Poly2) return false;
	else return true;
}


Polynomial operator+(const Polynomial & Poly, const Monomial & Mon){

	Polynomial out_Poly(Poly);
	// If original Polynomial was empty=1, push_back dummy Monomial=1
	// (such that after the addition the Polynomial has indeed two terms)
	if(out_Poly.empty()){
		Monomial dummy_Mon;
		out_Poly.push_back(dummy_Mon);
	}

	// Add term unless it's 0
	if(Mon.int_part != 0) out_Poly.push_back( Mon );

	return out_Poly;
}


Polynomial operator-( const Polynomial & Poly, const  Monomial & Mon){

	Monomial out_Mon(Mon);
	// Change sing of Monomial
	out_Mon.int_part*=(-1);
	// Add term unless it's 0
	return Poly+out_Mon;
}


Polynomial operator-(const Monomial & Mon, const Polynomial & Poly){

	Polynomial out_Poly(Poly);

	// If Poly is empty=1, add a default Monomial=1, which can change sign
	if (out_Poly.empty()) {
		Monomial dummy_Mon;
		out_Poly.push_back(dummy_Mon);
	}

	// Change sing of Polynomial, by looping over terms
	for (uint m=0; m < out_Poly.poly.size(); m++)
		out_Poly.poly.at(m).int_part*=(-1);

	// Add Polynomial and Monomial
	return out_Poly+Mon;
}


Polynomial operator+(const Monomial & Mon, const Polynomial & Poly ){

	return Poly + Mon;
}


Polynomial operator+=( Polynomial & Poly, const Monomial & Mon ) {

	// If original Poly was empty=1, push_back dummy Monomial=1
	// to make the Poly 1.
	if ( Poly.empty() ) {
		Monomial dummy_Mon;
		Poly.push_back( dummy_Mon );
	}

	if(Mon.int_part != 0) Poly.push_back( Mon );

	return Poly;
}


Polynomial operator+=( Polynomial & Poly1, const Polynomial & Poly2 ) {

	if (&Poly1 == &Poly2){
		std::cerr << "operator+=: Polynomials need to have different address, both arguments " << Poly1
				<< "."<<  std::endl;
		assert(0);
	}

	// If original Poly2 was empty=1, push_back dummy Monomial=1
	// add a dummy Monomial to Poly1
	if ( Poly2.empty() ) {
		Monomial dummy_Mon;
		Poly1.push_back( dummy_Mon );
	}
	// otherwise add Monomials from Poly2 to Poly 1
	else for (int i = 0; i < Poly2.size(); i++) {
		Poly1.push_back( Poly2.at(i) );
	}

	return Poly1;
}


Polynomial operator+( const Polynomial & Poly1, const Polynomial & Poly2) {

	Polynomial out_Poly1(Poly1);
	Polynomial out_Poly2(Poly2);

	// If original Poly1 was empty=1, push_back dummy Monomial=1
	// after adding empty Monomial the Polynomial is still 1
	if (out_Poly1.empty()) {
		Monomial dummy_Mon;
		out_Poly1.push_back(dummy_Mon);
	}
	// If original Poly2 was empty=1, push_back dummy Monomial=1
	// after adding empty Monomial the Polynomial is still 1
	if (out_Poly2.empty()) {
		Monomial dummy_Mon;
		out_Poly2.push_back(dummy_Mon);
	}

	// Add terms from Poly2 to Poly1
	for (int i = 0; i < out_Poly2.size(); i++) {
		out_Poly1 = out_Poly1 + out_Poly2.at(i);
	}

	return out_Poly1;
}


Polynomial operator-( const Polynomial & Poly1, const Polynomial & Poly2 ) {

	Polynomial out_Poly2(Poly2);

	// If Poly2 is empty=1 make it contain a Monomial,
	// as a Monomial by construction is 1
	if(out_Poly2.empty()) {
		Monomial dummy_Mon;
		out_Poly2.push_back(dummy_Mon);
	}

	// Change sing of Poly2, by looping over terms
	for (uint m=0; m < out_Poly2.poly.size(); m++)
		out_Poly2.poly.at(m).int_part*=(-1);

	// Now the subtraction is equal to an addition
	// if Poly1 is empty this is taken care of in + operator
	return Poly1+out_Poly2;
}


Polynomial operator*( const Polynomial & Poly, const int in ){

	Polynomial out_Poly(Poly);
	// If initially empty Polynomial=1, append empty Monomial=1, to have something to multiply with
	Monomial dummy_Mon;
	if (out_Poly.empty())
		out_Poly.push_back(dummy_Mon);

	// Loop over terms in Polynomial
	for (int i=0; i< out_Poly.size(); i++){
		out_Poly.at(i).int_part=out_Poly.at(i).int_part*in;
	}
	return  out_Poly;
}


Polynomial operator*( const int i, const Polynomial & Poly ) {
	return Poly * i;
}


Polynomial operator*( const Polynomial & Poly, const cnum c ) {

	Polynomial out_Poly(Poly);

	// If initially empty Polynomial=1, append empty Monomial=1, to have something to multiply with
	Monomial dummy_Mon;
	if (out_Poly.empty())
		out_Poly.push_back(dummy_Mon);

	// Loop over terms in Polynomial
	for (int i = 0; i < out_Poly.size(); i++) {
		out_Poly.poly.at(i).cnum_part = out_Poly.at(i).cnum_part * c;
	}
	return out_Poly;
}
Polynomial operator*( const cnum c, const Polynomial & Poly ) {
	return Poly * c;
}


Polynomial operator*( const Polynomial & Poly, const double d ) {

	Polynomial out_Poly(Poly);

	// If initially empty Polynomial=1, append empty Monomial=1, to have something to multiply with
	Monomial dummy_Mon;
	if (out_Poly.empty())
		out_Poly.push_back(dummy_Mon);

	// Loop over terms in Polynomial
	for (int i = 0; i < out_Poly.size(); i++) {
		out_Poly.poly.at(i).cnum_part = out_Poly.at(i).cnum_part * d;
	}
	return out_Poly;
}
Polynomial operator*(const double d, const Polynomial & Poly) {
	return Poly * d;
}


Polynomial operator*( const Polynomial & Poly, const Monomial & Mon ){

	// To contain the result
	Polynomial Poly_res;

	// Special case that the Polynomial is empty=1
	if( Poly.empty() ) {
		Poly_res.push_back(Mon);
		return Poly_res;
	}
	else{
		for (int i1=0; i1< Poly.size(); i1++){
			Poly_res.push_back( Poly.at(i1)*Mon );
		}
	}
	return Poly_res;
}
Polynomial operator*( const Monomial & Mon, const Polynomial & Poly ){

	return  Poly*Mon;
}


Polynomial operator*( const Polynomial & Poly1, const Polynomial & Poly2 ){

	// To contain the result
	Polynomial Poly_res;

	// Special case that the vectors are empty=1, needed for special treatment of sum of 0-monomials
	// Poly_res is also empty, and thus =1, 1*1=1
	if(Poly1.empty() && Poly2.empty()){
		return Poly_res;
	}

	// If one Poly is empty, =1, return other
	if(!Poly1.empty() && Poly2.empty() ) return Poly1;
	else if( Poly1.empty() && !Poly2.empty()) return Poly2;
	// If both vectors contain info
	else {
		for (int i1=0; i1< Poly1.size(); i1++){
			for (int i2=0; i2< Poly2.size(); i2++){
				// Add terms if non of them is 0, else just don't add
				if( Poly1.at(i1).int_part!=0 && Poly2.at(i2).int_part!=0  )
					Poly_res.push_back(Poly1.at(i1)*Poly2.at(i2));
			}
		}
		// If, by now, the Poly_res is empty that's because all terms were 0
		if( Poly_res.empty() ){
			Monomial Mon_tmp;
			Mon_tmp.int_part=0;
			Poly_res.push_back( Mon_tmp );
		}
	}

	return Poly_res;
}

}// end namespace ColorFull



