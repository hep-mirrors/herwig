// -*- C++ -*-

/*
 * Poly_matr.cc
 * Contains the definition of the class Poly_matr and related operators.
 * Created on: Aug 5, 2013
 * Author: Malin Sjodahl
 */

#include "Poly_matr.h"
#include <cassert>
#include <fstream>
#include <iostream>


namespace ColorFull{


std::ostream& operator<<(std::ostream& out, const poly_matr & pm){

	out <<"{" <<std::endl;
	// Loop over rows
	for(uint i=0; i< pm.size(); i++ ){
		out <<"{";
		// Loop over columns
		for(uint j=0; j< pm.at(i).size(); j++ ){
			// Print element
			std::cout.width( 20 );
			std::ostringstream outstr;
			outstr << pm.at(i).at(j);
			// If not last element print ","
			if (j<pm.at(i).size()-1 ) outstr << ",";
			out << outstr.str();
			//out << pm.at(i).at(j) << "";
		}
		out <<"}";
		// If not last row, print ","
		if (i<pm.at(i).size()-1 ) out << ",";
		out << std::endl;
	}
	out <<"}" <<std::endl;
	return out;
}


void Poly_matr::remove_CF()  {

	// Loop over vectors matrix
	for ( uint i = 0; i < pm.size(); i++ ) {
		// Remove CF in each term
		pm.at(i).remove_CF();
	}
}


std::ostream& operator<<(std::ostream& out, const Poly_matr & Pm){
	out << Pm.pm;
	return out;
}


bool operator==( const Poly_matr & Pm1, const Poly_matr & Pm2 ){

	if( Pm1.size() != Pm2.size() ) return false;

	// All terms should be the same
	for ( uint i=0; i < Pm1.size(); i++ ){
		if( Pm1.at(i) != Pm2.at(i ) ) return false;
	}

	return true;
}

bool operator!=(  const Poly_matr & Pm1, const Poly_matr & Pm2 ){
	if( Pm1==Pm2 ) return false;
	else return true;
}



void Poly_matr::normal_order( ) {

	// Loop over vectors in matrix
	for ( uint i = 0; i < pm.size(); i++ ) {
		// Remove CF in each Poly_vec
		pm.at(i).normal_order() ;
	}
}


void Poly_matr::simplify( ) {

	// Loop over vectors in matrix
	for ( uint i = 0; i < pm.size(); i++ ) {
		// Remove CF in each term
		pm.at(i).simplify();
	}

}


void Poly_matr::conjugate( ) {

  for ( uint i=0; i < pm.size(); i++ ){
	  pm.at(i).conjugate();
  }
}


void Poly_matr::read_in_Poly_matr( std::string filename ) {

	// Read in file
	std::ifstream fin( filename.c_str() );

	// Check that file exists
	if( !fin ){
		std::cerr << "Poly_matr::read_in_Poly_matr: The file "
				<< filename << " could not be opened." << std::endl;
		assert(0);
	}

	// Copy info from file to string
	std::string str((std::istreambuf_iterator<char>(fin)),
			std::istreambuf_iterator<char>());

	// Skip lines starting with #
	while(str.at(0)== '#'){
		while (str.at(0) != '\n'){
			str.erase(str.begin());
		}
		// erase endl sign(s)
		while(str.at(0)== '\n'){
			str.erase(str.begin());
		}
	}

	// First char in file should be '{'
	if (str.at(0) != '{') {
		std::cerr
		<< "Poly_matr::read_in_Poly_matr: First char in matrix data after comments file should be '{', it was: "
		<< str.at(0) << std::endl;
		assert( 0 );
	}

	// Row to contain numbers
	Poly_vec row;

	// Read the string, starting from 0th element
	uint i = 0;
	while (i < str.size() - 2) {

		// To contain the Polynomial string
		std::string Poly_str;
		Poly_str.clear();

		// We may have to skip some spaces, end-lines and {
		while (i < str.size() - 2 && (str.at(i) == ' ' or str.at(i) == '\n' or str.at(i) == '{'))
			i++;

		// Read next Polynomial, until , or {
		// Keep reading the number while not ',' or '}'
		while (i < str.size() - 2 && (str.at(i) != ',' && str.at(i) != '}')) {
			Poly_str.push_back(str.at(i));
			i++;
		}

		Polynomial Poly( Poly_str );
		row.append( Poly );

		// If we have a new row
		if (i < str.size() - 2 && str.at(i) == '}') {
			// Save row in matrix, and empty row
			pm.push_back(row);
			row.clear();

			// We may have to skip some chars
			while (i< str.size()-2 &&(str.at(i) == ',' or str.at(i) == '}' or str.at(i) == ' ' or str.at(i) == '\n' ) ) i++;
		}
		i++;
	}
}

void Poly_matr::write_out_Poly_matr( std::string filename ) const {
	std::ofstream outfile( filename.c_str() );

	if ( !outfile )
	std::cerr << "Poly_matr::write_out_Poly_matr: Cannot write out Polynomial matrix as the file \""
		<< filename.c_str() << "\" could not be opened. (Does the directory exist? Consider creating the directory.)" << std::endl;
	outfile << pm;
}

}// end namespace ColorFull
