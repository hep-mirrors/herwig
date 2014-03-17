// -*- C++ -*-

/*
 * Poly_vec.cc
 * Contains the definition of the class Poly_vec related operators.
 * Created on: Aug 5, 2013
 * Author: Malin Sjodahl
 */

#include "Poly_vec.h"
#include <cassert>
#include <fstream>
#include <iostream>

namespace ColorFull{


std::ostream& operator<<(std::ostream& out, const poly_vec & poly_v ){
	out <<"{";
	// Loop over entries
	for( uint i=0; i< poly_v.size(); i++ ){
		out << poly_v.at(i);
		// If not last element print ","
		if (i<poly_v.size()-1 ) out << ", ";
	}
	out <<"}";
	return out;
}


std::ostream& operator<<(std::ostream& out, const Poly_vec & Pv ){
	out << Pv.pv;
	return out;
}


void Poly_vec::remove_CF()  {

	// Loop over entries in vector (Polynomials)
	for ( uint i = 0; i < pv.size(); i++ ) {
		// Remove CF in each term
		pv.at(i).remove_CF();
		//std::cout << "Col_functions::remove_CF: Pv_res.at(i)" << Pv_res.at(i) << std::endl;
	}
}


void Poly_vec::normal_order()  {

	// Loop over entries in vector (Polynomials)
	for ( uint i = 0; i < pv.size(); i++ ) {
		// Normal order each term
		pv.at(i).normal_order();
	}
}


void Poly_vec::simplify( ) {

	// Loop over entries in vector (poly_vec)
	for ( uint i = 0; i < pv.size(); i++ ) {
		// Remove CF in each term
		pv.at(i).simplify();
	}
}


void Poly_vec::conjugate( ) {

	for ( uint i=0; i < pv.size(); i++ ){
		pv.at(i).conjugate();
	}
}


void Poly_vec::read_in_Poly_vec( std::string filename ) {

	// Read in file
	std::ifstream fin(filename.c_str() );

	// Check that file exists
	if( !fin ){
		std::cerr << "Poly_vec::read_in_Poly_vec: The file "
				<< filename << " could not be opened." << std::endl;
		assert(0);
	}

	// Copy info from file to string
	std::string str((std::istreambuf_iterator<char>(fin)),
			std::istreambuf_iterator<char>());

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

	// First char in file should be '{'
	if ( str.at(0) != '{' ) {
		std::cerr
		<< "Poly_vec::read_in_Poly_vec: First char in vector data after comments should be '{', it was: "
		<< str.at(0) << std::endl;
		assert(0);
	}

	// Row to contain numbers
	Poly_vec row;

	// Read the string, starting from 0th element
	unsigned int i = 0;
	while ( i < str.size() - 2 ) {

		// To contain the Polynomial string
		std::string Poly_str;
		Poly_str.clear();

		// We may have to skip some spaces, end-lines and {
		while (i < str.size() - 2 && (str.at(i) == ' ' or str.at(i) == '\n' or str.at(i) == '{'))
			i++;

		// Keep reading the number while not ',' or '}'
		while (i < str.size() - 2 && (str.at(i) != ',' && str.at(i) != '}')) {
			Poly_str.push_back(str.at(i));
			i++;
		}

		Polynomial Poly( Poly_str );
		pv.push_back( Poly );

		// If we have a new row
		if (i < str.size() - 2 && str.at(i) == '}') { i++;}
		i++;
	}
}


void Poly_vec::write_out_Poly_vec( std::string filename ) const {
	std::ofstream outfile( filename.c_str() );

	if ( !outfile )
	std::cerr << "Poly_vec::write_out_Poly_vec: Cannot write out vector of Polynomials as the file \""
		<< filename.c_str() << "\" could not be opened. (Does the directory exist? Consider creating the directory.)" << std::endl;


	outfile << pv;
}

} // end namespace ColorFull




