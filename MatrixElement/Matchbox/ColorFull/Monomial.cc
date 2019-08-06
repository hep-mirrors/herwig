// -*- C++ -*-
/*
 * Monomial.cc
 *	Contains definition of the class Monomial and associated types and operators
 *  Created on: Jul 7, 2010
 *  Author: Malin Sjodahl
 */
#include "Monomial.h"
#include "parameters.h"
#include <cassert>
#include <fstream>
#include <iostream>

namespace ColorFull {


Monomial::Monomial( std::string str ){
	Monomial_of_str( str );
}


void Monomial::Monomial_of_str( std::string str ) {
	//std::cout << "Monomial::Monomial: got string \"" << str << "\"" << std::endl;

	// Strings to contain various parts
	std::string num_str, Nc_pow_str, TR_pow_str, CF_pow_str, int_part_str;
	// Is it in denominator?
	bool denominator=0;

	// Start with setting to default
	pow_TR = pow_Nc = pow_CF = 0;
	int_part = 1;
	cnum_part = 1.0;

	// Check that left and right brackets match up
	int left_brackets=0, right_brackets=0;
	uint j=0;
	while (j < str.size()) {
		if(str.at(j)=='(') left_brackets++;
		if(str.at(j)==')') right_brackets++;
		j++;
	}
	if(left_brackets != right_brackets){
		std::cerr << "Monomial::Monomial_of_str: The brackets in the monomial\"" << str <<"\" do not seem to match up. There were "
				<< left_brackets <<" left bracket(s) and "<< right_brackets << " right bracket(s)." << std::endl;
		assert( 0 );
	}

	// Check that only allowed signs
	j = 0;
	while (j < str.size()) {

		if (!(str.at(j) == '+' or str.at(j) == '-' or str.at(j) == '.' or str.at(j) == '\n'
				or str.at(j) == '*' or str.at(j) == '/' or str.at(j) == '^'
						or str.at(j) == 'T' or str.at(j) == 'R' or str.at(j) == 'N'
								or str.at(j) == 'c' or str.at(j) == 'C' or str.at(j) == 'F'
										or str.at(j) == '('
												or str.at(j) == ')' or str.at(j) == ' ' or str.at(j) == '0'
														or str.at(j) == '1' or str.at(j) == '2' or str.at(j) == '3'
																or str.at(j) == '4' or str.at(j) == '5' or str.at(j) == '6'
																		or str.at(j) == '7' or str.at(j) == '8' or str.at(j) == '9')) {
			std::cerr
			<< "Monomial::Monomial_of_str: A disallowed sign encountered in string for Monomial constructor: "
			<< str.at(j) << std::endl;
			assert( 0 );

		}

		num_str.push_back( str.at(j) );
		j++;
	}

	// Read the string, starting from 0th element
	uint i = 0;

	// We may have to skip some chars containing spaces
	while (i < str.size() && (str.at(i) == ' ' or str.at(i) == '\n'  or str.at(i) == '*' or str.at(i) == '(' or str.at(i) == ')')) i++;
	// Check plus or minus
	if (i < str.size() and str.at(i) == '+') {
		i++;
	}
	// Pick up the sign
	if (i < str.size()  and str.at(i) == '-' ) {
		int_part = int_part * (-1);
		i++;
	}

	while (i < str.size()) {

		// Clear strings as we may loop around many times
		num_str.clear();
		Nc_pow_str.clear();
		TR_pow_str.clear();
		CF_pow_str.clear();

		// We may have to skip some chars containing spaces and *
		while (i < str.size() && (str.at(i) == ' ' or str.at(i) == '\n' or str.at(i) == '*' or str.at(i) == '(' or str.at(i) == ')')) i++;

		// look for denominator "/"
		if(i< str.size() and  str.at(i) == '/') {

			if (denominator) {
				std::cerr << "Monomial::Monomial_of_str: Can only handle one / "
						<< "but got " << str;
				assert( 0 );
			};
			denominator=true;
			i++;
		}

		// If a number is encountered which is not *-1 or *(-1)
		while (i< str.size() and ( str.at(i) == '.' or str.at(i) == '0' or str.at(i) == '1'
				or str.at(i) == '2' or str.at(i) == '3' or str.at(i) == '4'
						or str.at(i) == '5' or str.at(i) == '6' or str.at(i) == '7'
								or str.at(i) == '8' or str.at(i) == '9')  ) {
			num_str.push_back(str.at(i));
			i++;
		}

		// If the num_str contains something, save
		if (!num_str.empty()) {
			double num_fac=1.0;

			// Make a number of the string
			std::istringstream num_str_st(num_str);
			// Set the num member variable
			num_str_st >> num_fac;

			if(denominator) cnum_part=cnum_part/num_fac;
			else { // factor is in numerator, try to put in int-part if it looks like an int
				int num_fac_int=static_cast<int>(floor(num_fac + 0.50));
				if( std::abs(num_fac_int-num_fac) < accuracy )
					int_part= int_part*num_fac_int;
				else // put factor in numeric part
					cnum_part=cnum_part*num_fac;
			}
		}

		// If TR is encountered (T should always be followed by R)
		if (i< str.size() and (str.at(i) == 'T' and (i==str.size()-1  or str.at(i + 1) != 'R')) ){
			std::cerr << "Monomial::Monomial_of_str: Got a string containing T, but T was not followed by R (as in TR). " << std::endl;
			assert( 0 );
		}
		if (i< str.size() and (str.at(i) == 'T' and str.at(i + 1) == 'R') ) {
			i++;
			i++;
			// If there is no ^, the power of TR is 1
			int pow_cont=1;

			// If ^, start reading in power of TR
			if (i < str.size() and str.at(i) == '^') {
				i++;

				// If we have something in a bracket
				if(i < str.size() and  str.at(i)=='('){
					i++;
					// keep reading until end-bracket
					while(i < str.size() and str.at(i)!=')'){
						TR_pow_str.push_back(str.at(i));
						i++;
					}
					// Skip the )
					if(i < str.size() and  str.at(i)==')') i++;
				}
				// otherwise, as long as we have numbers
				else while (i< str.size() and (str.at(i) == '0' or str.at(i) == '1' or str.at(i) == '2'
						or str.at(i) == '3' or str.at(i) == '4'
								or str.at(i) == '5' or str.at(i) == '6'
										or str.at(i) == '7' or str.at(i) == '8'
												or str.at(i) == '9' or str.at(i) == '-')) {
					TR_pow_str.push_back(str.at(i));
					i++;
				}
				// Make a number of the string
				std::istringstream TR_pow_str_st(TR_pow_str);
				// Set the TR_pow member variable
				TR_pow_str_st >> pow_cont;
			}

			if(denominator) pow_cont=-pow_cont;
			pow_TR=pow_TR + pow_cont;

		}

		// If Nc is encountered
		if (i < str.size() and (str.at(i) == 'N')) {
			if( i+1< str.size() and str.at(i + 1) == 'c') i++; // get to c
			else{// allow only Nc
				std::cerr << "Monomial::Monomial_of_str: got a string containing N, " << str <<", but N was not followed by c (as in Nc). " << std::endl;
				assert( 0 );
			}
			i++; // get to next sign

			int Nc_pow = 1;
			//if(i < str.size())  std::cout << "Monomial::Monomial: " <<  str.at(i) << endl;
			// If ^, start reading in power of Nc
			if (i< str.size() and str.at(i) == '^') {
				i++; // compensate for ^

				// If we have something in a bracket
				if(i < str.size() and  str.at(i)=='('){
					i++;
					// keep reading until end-bracket
					while(i < str.size() and  str.at(i)!=')'){
						Nc_pow_str.push_back(str.at(i));
						i++;
					}
					// Skip the )
					if(i < str.size() and str.at(i)==')') i++;
				}
				// as long as we have numbers or -
				else while (i< str.size() and  (str.at(i) == '0' or str.at(i) == '1' or str.at(i) == '2'
						or str.at(i) == '3' or str.at(i) == '4'
								or str.at(i) == '5' or str.at(i) == '6'
										or str.at(i) == '7' or str.at(i) == '8'
												or str.at(i) == '9' or str.at(i) == '-')){
					Nc_pow_str.push_back(str.at(i));
					i++;
				}
				// Make a number of the string
				std::istringstream Nc_pow_str_st(Nc_pow_str);
				Nc_pow_str_st >>  Nc_pow ;
			}
			if (denominator) pow_Nc=pow_Nc-Nc_pow;
			else pow_Nc=pow_Nc + Nc_pow;
		}

		// If CF is encountered
		if (i < str.size() and (str.at(i) == 'C')) {
			std::cout.flush();
			if( i+1< str.size() and str.at(i + 1) == 'F') i++; // get to f
			i++; // get to next sign

			int CF_pow = 1;
			std::cout.flush();
			// If ^, start reading in power of CF
			if (i< str.size() and str.at(i) == '^') {
				i++; // compensate for ^

				// If we have something in a bracket
				if(i < str.size() and  str.at(i)=='('){
					i++; // Skip the (
					// keep reading until end-bracket
					while(i < str.size() and  str.at(i)!=')'){
						CF_pow_str.push_back(str.at(i));
						i++;
					}
					// Here we may need to skip some spaces as well
					if(i < str.size() and str.at(i)==' ') i++;
					// Skip the )
					if(i < str.size() and str.at(i)==')') i++;
				}
				// as long as we have numbers
				else while (i< str.size() and  (str.at(i) == '0' or str.at(i) == '1' or str.at(i) == '2'
						or str.at(i) == '3' or str.at(i) == '4'
								or str.at(i) == '5' or str.at(i) == '6'
										or str.at(i) == '7' or str.at(i) == '8'
												or str.at(i) == '9' or str.at(i) == '-')){
					CF_pow_str.push_back(str.at(i));
					i++;
				}
				// Make a number of the string
				std::istringstream CF_pow_str_st(CF_pow_str);
				CF_pow_str_st >>  CF_pow ;
			}
			if (denominator) pow_CF=pow_CF-CF_pow;
			else pow_CF=pow_CF+ CF_pow;
		}

		// If *-1
		if ((i < str.size() and str.at(i) == '-') and  ( (i-1>0) and str.at(i-1)=='*')  ) {
			i++; //skip the -
			std::cout.flush();

			// Keep reading in while numbers
			while (i< str.size() and  (str.at(i) == '0' or str.at(i) == '1' or str.at(i) == '2'
					or str.at(i) == '3' or str.at(i) == '4'
							or str.at(i) == '5' or str.at(i) == '6'
									or str.at(i) == '7' or str.at(i) == '8'
											or str.at(i) == '9')){
				int_part_str.push_back(str.at(i));
				i++;
			}
			int n_fac=1;
			std::istringstream n_str(int_part_str);
			n_str >> n_fac;
			int_part=int_part * n_fac*(-1);
		}
		// If *(-1)
		if ((i < str.size() and str.at(i) == '-') and (  ( ( (i-2> 0) and str.at(i-1)=='(') and  str.at(i-2)=='*') )) {
			i++; //skip the -
			i++; //skip the (
			std::cout.flush();

			// Keep reading in while numbers
			while (i< str.size() and  (str.at(i) == '0' or str.at(i) == '1' or str.at(i) == '2'
					or str.at(i) == '3' or str.at(i) == '4'
							or str.at(i) == '5' or str.at(i) == '6'
									or str.at(i) == '7' or str.at(i) == '8'
											or str.at(i) == '9')){
				int_part_str.push_back(str.at(i));
				i++;
			}
			int n_fac=1;
			std::istringstream n_str(int_part_str);
			n_str >> n_fac;
			int_part=int_part * n_fac*(-1);
			// Skip potential white spaces
			while (i< str.size() and  str.at(i) == ' ') i++;
			// Make sure we have a closing )
			if(str.at(i)==')') i++;
			else {
				std::cerr << "Monomial::Monomial_of_str: Expects final ) in *(number)"
						<< " got " << str.at(i) << std::endl;
				assert( 0 );
			}
		}

		// We may have to skip some spaces
		while (i < str.size() && str.at(i) == ' ') {i++;}
		std::cout.flush();

	}
}// end Monomial_of_str


void Monomial::read_in_Monomial( std::string filename ) {

	// Read in file
	std::ifstream fin( filename.c_str() );

	// Check that file exists
	if( !fin ){
		std::cerr << "Monomial::read_in_Monomial: The file "
				<< filename << " could not be opened." << std::endl;
		assert( 0 );
	}

	// Copy info from file to string
	std::string str((std::istreambuf_iterator<char>(fin)), std::istreambuf_iterator<char>());

	// Skip lines starting with #
	while( str.at(0)== '#' ){
		while (str.at(0) != '\n'){
			str.erase(str.begin());
		}
		// erase endl sign(s)
		while(str.at(0)== '\n'){
			str.erase(str.begin());
		}
	}

	// Remove endl chars at the end of the file
	while( str.at(str.size()-1) == '\n' ) str.erase(str.size()-1);
	Monomial_of_str( str );
}


void Monomial::write_out_Monomial( std::string filename ) const {

	std::ofstream outfile(filename.c_str());

	if ( !outfile )
	std::cerr << "Monomial::write_out_Monomial: Cannot write out Monomial as the file \""
		<< filename.c_str() << "\" could not be opened. (Does the directory exist? Consider creating the directory.)" << std::endl;

	outfile << *this;
}

std::ostream& operator<<( std::ostream& out, const Monomial & Mon ){

	// If multiplied by 0, the whole Monomial is 0
	if (Mon.int_part == 0 ) out << "0";
	else{
		// If real, print only real part
		if( imag( Mon.cnum_part )==0 ) {
			// If both int and numeric parts are 1, print 1
			if( Mon.int_part==1 && Mon.cnum_part.real() ==1 ) out << "1";
			// If int part is 1, print numeric part if >0 (minus signs gives problems)
			else if( Mon.int_part==1 && Mon.cnum_part.real() > 1 ) out << Mon.cnum_part.real();
			// If numeric part is 1, print int part if >0 (minus signs gives problems)
			else if( Mon.int_part > 1 && Mon.cnum_part.real() ==1 ) out << Mon.int_part;
			else out << real( Mon.cnum_part ) << "*" << Mon.int_part;
		}
		// else print full complex number
		else out << Mon.cnum_part  << "*" << Mon.int_part;

		if( Mon.pow_TR !=0 ) {
			if( Mon.pow_TR ==1 ) out << " TR";
			else out << " TR^" <<"(" << Mon.pow_TR <<")";
		}
		if( Mon.pow_Nc !=0 ){
			if( Mon.pow_Nc ==1 ) out << " Nc";
			else out <<" Nc^" << "(" << Mon.pow_Nc <<")";
		}
		if( Mon.pow_CF !=0 ) {
			if( Mon.pow_CF ==1 ) out << " CF";
			else out << " CF^" << "(" << Mon.pow_CF <<")";
		}
	}
	return out;
}


Monomial operator*( const Monomial & Mon, const int i ){
	Monomial out_Mon(Mon);
	out_Mon.int_part=out_Mon.int_part*i;
	return out_Mon;
}


Monomial operator*( const int i, const Monomial & Mon){
	return Mon*i;
}

Monomial operator*=( Monomial & Mon, const int i ){
	return Mon.int_part*=i;
}


Monomial operator*( const Monomial & Mon, const cnum c ){
	Monomial out_Mon(Mon);
	out_Mon.cnum_part=out_Mon.cnum_part*c;
	return out_Mon;
}


Monomial operator*( const cnum c, const Monomial & Mon ){
  return Mon*c;
}

Monomial operator*=( Monomial & Mon, const cnum c ){
	Mon.cnum_part*=c;
	return Mon;
}

Monomial operator*(const Monomial & Mon, const double d){
		Monomial out_Mon(Mon);
		out_Mon.cnum_part=out_Mon.cnum_part*d;
		return out_Mon;
}

Monomial operator*(const double d, const Monomial & Mon){
	return Mon*d;
}

Monomial operator*=( Monomial & Mon, const double d ){
	Mon.cnum_part*=d;
	return Mon;
}

Monomial operator*(const Monomial & Mon1, const Monomial & Mon2){
  Monomial Mon_out;

  // Adding powers to get total power of TR, Nc and CF
  Mon_out.pow_TR = Mon1.pow_TR + Mon2.pow_TR;
  Mon_out.pow_Nc = Mon1.pow_Nc + Mon2.pow_Nc;
  Mon_out.pow_CF = Mon1.pow_CF + Mon2.pow_CF;

  // Multiplying factors
  Mon_out.int_part=Mon1.int_part * Mon2.int_part;
  Mon_out.cnum_part=Mon1.cnum_part * Mon2.cnum_part;
  return Mon_out;
}

Monomial operator*=( Monomial & Mon1, const Monomial & Mon2){

  // Adding powers to get total power of TR, Nc and CF
  Mon1.pow_TR += Mon2.pow_TR;
  Mon1.pow_Nc += Mon2.pow_Nc;
  Mon1.pow_CF += Mon2.pow_CF;

  // Multiplying factors
  Mon1.int_part *= Mon2.int_part;
  Mon1.cnum_part *= Mon2.cnum_part;
  return Mon1;
}


bool operator==(const Monomial & Mon1 , const Monomial & Mon2){

  // If both are 0
  if(Mon1.int_part == 0  && Mon2.int_part == 0) return true;

  // all factors should be same in order for Monomial's to be same
  if(Mon1.pow_TR != Mon2.pow_TR ) return false;
  if(Mon1.pow_Nc != Mon2.pow_Nc ) return false;
  if(Mon1.pow_CF != Mon2.pow_CF ) return false;
  if(Mon1.int_part != Mon2.int_part ) return false;

  // For comparing numerical part use an accuracy, as numbers may be small, compare ratio
  const double r1 = real(Mon1.cnum_part);
  const double r2 = real(Mon2.cnum_part);
  const double i1 = imag(Mon1.cnum_part);
  const double i2 = imag(Mon2.cnum_part);

  if( r2 != 0 and i2 != 0 ) {
    if( r1/r2  < (1-accuracy) or r1/r2  > (1+accuracy) ) return false;
    if( i1/i2  < (1-accuracy) or i1/i2  > (1+accuracy) ) return false;
  }
  else if( r2 == 0 and i2 != 0 )  { 
    if( r1 != 0 ) return false;
    if( i1/i2  < (1-accuracy) or i1/i2  > (1+accuracy) ) return false;
  }
  else if ( r2 != 0 and i2 == 0 ) { 
    if( r1/r2  < (1-accuracy) or r1/r2  > (1+accuracy) ) return false;
    if( i1 != 0 ) return false;
  }
  else { // both r2,i2 are zero
    if( r1 != 0 ) return false;
    if( i1 != 0 ) return false;
  }

  return true;
}


bool operator!=(const Monomial & Mon1 , const Monomial & Mon2){
  if( Mon1==Mon2) return false;
  else return true;
}


bool operator<(const Monomial & Mon1, const Monomial & Mon2) {

	// If different Nc+CF power
	if (Mon1.pow_Nc + Mon1.pow_CF < Mon2.pow_Nc + Mon2.pow_CF)
		return true;
	else if (Mon1.pow_Nc + Mon1.pow_CF > Mon2.pow_Nc + Mon2.pow_CF)
		return false;
	else {	  // same total Nc power
		// if different powers of N
		if (Mon1.pow_Nc < Mon2.pow_Nc)
			return true;
		else if (Mon1.pow_Nc > Mon2.pow_Nc)
			return false;
		else { // same pow_Nc+pow_CF and same pow_Nc
			// order according to cnum_part*int_part
			if ( (Mon1.int_part)*(abs(Mon1.cnum_part)) < (Mon2.int_part)*(abs(Mon2.cnum_part)))
				return true;
			else if ( (Mon1.int_part)*(abs(Mon1.cnum_part)) > (Mon2.int_part)*(abs(Mon2.cnum_part)))
				return false;
			else { // same pow_Nc+pow_CF and same pow_Nc, and same numerical cnum_part*int_part
				// order according to int_part
				if (Mon1.int_part < Mon2.int_part)
					return true;
				else if (Mon1.int_part < Mon2.int_part)
					return false;
				else {// same pow_Nc+pow_CF and same pow_Nc, and same numerical cnum_part*int_part, and same int_part
					// order according to pow_TR
					if (Mon1.pow_TR < Mon2.pow_TR)
						return true;
					else if (Mon1.pow_TR < Mon2.pow_TR)
						return false;
					else { // The Monomials are equal
						return false;
					}
				}
			}
		}

	}
	std::cerr << "Monomial::operator<: Could not decide on Monomial ordering of "
			<< Mon1 << " and " << Mon2 << ". This should not happen, report bug."
			<< std::endl;
	assert( 0 );
	return false;

}

} // end namespace ColorFull
