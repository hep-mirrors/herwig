/*
 * Orthogonal_basis.cc
 * Contains the definitions of the class Orthogonal_basis, related types and operators.
 *  Created on: May 25, 2013
 *      Author: Malin Sjodahl
 */

#include "Orthogonal_basis.h"
#include "parameters.h"
#include <cassert>
#include <fstream>
#include <iostream>


namespace ColorFull {


void Orthogonal_basis::scalar_product_matrix(){
	if( ng+nq>5 ){
		std::cout << "Orthogonal_basis::scalar_product_matrix: nq+n_g0=" << nq+ng << " is large, consider using numerical and/or memory version.  "  << std::endl;
		std::cout.flush();
	}
	return scalar_product_matrix( true, true, false );
}


void Orthogonal_basis::diagonal_scalar_product_matrix( bool save_P_diagonal_spm, bool save_d_diagonal_spm, bool use_mem ){

	// If the diagonal_P_spm and diagonal_d_spm have already been calculated, erase them
	diagonal_P_spm.clear();
	diagonal_d_spm.clear();

	if(  (cb.size()==0 ) ) {
		std::cout << "Orthogonal_basis::diagonal_scalar_product_matrix: There are no basis vectors in this basis, consider using read_in_basis." << std::endl;
		std::cout.flush();
		return ;
	}


	// Check that all Polynomials are real
	// (check doesn't cover the case of complex Polynomial in a Quark_line)
	// but if an entry is not real this will be discovered as the d_spm is calculated
	for ( uint cbi=0; cbi < cb.size(); cbi++){
		for ( uint j=0; j< cb.at(cbi).size(); j++){
			if( imag (Col_fun.cnum_num(  cb.at(cbi).at(j).Poly  )) > accuracy ){
				std::cerr << "Orthogonal_basis::diagonal_scalar_product_matrix: ColorFull expects real Polynomial multiplying the color structure, but the Polynomial\n"
						<< cb.at(cbi).at(j).Poly  << std::endl
						<< " appearing in front of the Col_str " << cb.at(cbi).at(j).cs
						<< " in basis vector number " << cbi << " is not real."<< std::endl << std::endl;
				std::cerr.flush();
				assert( 0 );
			}
		}
	}

	// For remembering already calculated topologies
	std::map<std::string, std::shared_ptr<Polynomial> > mem_map;

	// Loop over basis vectors in Basis
	for( uint i=0; i < cb.size(); i++){

		// To contain the result of vector i square
		Polynomial iiRes;
		iiRes=iiRes*0;

		if(use_mem){

			// Loop over Col_strs in first Ca
			for( uint Ca1i=0; Ca1i< cb.at(i).size(); Ca1i++){

				// Loop over Col_strs in second Ca
				for( uint Ca2i=0; Ca2i< cb.at(i).size(); Ca2i++){

					// To contain the contribution to the ii-th entry if memoization is used
					std::shared_ptr<Polynomial> iiEntry_contr;

					// Rename indices, and make string of new col_strs, to use in map
					// strip off Polynomial information to make the map minimal
					Col_str Cs1, Cs2;
					Cs1.append(cb.at(i).at(Ca1i).cs);
					Cs2.append(cb.at(i).at(Ca2i).cs);

					rename_indices( Cs1, Cs2 );

					std::ostringstream Cs_string;
					Cs_string << Cs1 << Cs2;

					// Calculate element contribution to the ijth element in scalar product matrix
					// If this has scalar product has occurred before, reuse old value
					if (mem_map.count( Cs_string.str() ) > 0) {
						iiEntry_contr = mem_map[Cs_string.str()];
					}
					// Otherwise, calculate value and save topology
					else {
						Polynomial p = Col_fun.scalar_product(Cs1, Cs2);
						iiEntry_contr = shared_ptr<Polynomial> (new Polynomial(p));
						mem_map[Cs_string.str()] = iiEntry_contr;

					}

					// Sum up all the contributions to one Polynomial, recall to multiply with Polynomials
					Polynomial iiEntry_contr_poly=(*iiEntry_contr)* cb.at(i).at(Ca1i).Poly*cb.at(i).at(Ca2i).Poly;


					if( !save_P_diagonal_spm ) {
						double iiEntry_contr_d= Col_fun. double_num( iiEntry_contr_poly );
						Monomial Mon( iiEntry_contr_d );
						iiRes+= Mon;
					}
					else iiRes += iiEntry_contr_poly;


				} // end looping over Cs in first Ca
			} // end looping over Cs in 2nd Ca
			// If Polynomial result is wanted, simplify
			if ( save_P_diagonal_spm ) iiRes.simplify();

			// Otherwise convert to numerical
			else {
				double num_iiRes= Col_fun. double_num( iiRes );
				iiRes.clear();
				Monomial Mon( num_iiRes );
				iiRes.push_back(Mon);
			}

			diagonal_P_spm.push_back( iiRes );

		} // end if ( use_mem )
		else{
			// Calculate element ij in scalar product matrix
			Polynomial iiRes=Col_fun.scalar_product(cb.at(i), cb.at(i));

			iiRes.simplify();
			diagonal_P_spm.push_back( iiRes );

		} //end if (not mem)

	}// end looping over i

	if ( save_d_diagonal_spm ){
		diagonal_d_spm=Col_fun.double_num(diagonal_P_spm);
	}

	if (! save_P_diagonal_spm){
		diagonal_P_spm.clear();
	}

	return;
}


void Orthogonal_basis::scalar_product_matrix( bool save_P_spm, bool save_d_spm, bool use_mem ) {

	// If the P_spm and d_spm have already been calculated, erase them
	P_spm.clear();
	d_spm.clear();

	// First calculate diagonal version
	diagonal_scalar_product_matrix( save_P_spm, save_d_spm, use_mem );

	// Then copy content to matrix
	if( save_P_spm  ){
		// First empty P_spm
		P_spm.clear();

		Polynomial Zero;
		Zero=0*Zero;
		for (uint i=0; i< diagonal_P_spm.size(); i++ ){
			Poly_vec rowi;
			for (uint j=0; j< diagonal_P_spm.size(); j++ ){
				if(i==j){// diagonal entries
					rowi.push_back( diagonal_P_spm.at(i));
				}
				else{// non-diagonal parts
					rowi.push_back(Zero);
				}
			}
			P_spm.push_back(rowi);

		}
	}
	if( save_d_spm  ){
		// First empty P_spm
		d_spm.clear();
		for (uint i=0; i< diagonal_d_spm.size(); i++ ){
			dvec rowi;
			for (uint j=0; j< diagonal_d_spm.size(); j++ ){
				if(i==j){// diagonal entries
					rowi.push_back( diagonal_d_spm.at(i));
				}
				else{// non-diagonal parts
					rowi.push_back( 0 );
				}
			}
			d_spm.push_back(rowi);
		}
	}

}


Poly_vec Orthogonal_basis::decompose( const Col_amp & Ca ) {

	// Check that we have a basis
	if(cb.size()==0){
		std::cerr << "Orthogonal_basis::decompose: The basis vector cb is empty consider reading in basis." << std::endl;
		assert( 0 );
	}

	// Check that quark and gluon content agree with that in Col_amp
	else if(Ca.size()>0 ){
		if(Ca.at(0).n_quark() != nq)
		{
			std::cerr << "Orthogonal_basis::decompose: The number of quarks in the argument Col_amp, " <<  Ca.at(0).n_quark()
					<< ", does not fit the number of quarks in the basis "
					<< nq << std::endl;}
		if(Ca.at(0).n_gluon() != ng )
		{
			std::cerr << "Orthogonal_basis::decompose: The number of gluons in the argument Col_amp " << Ca.at(0).n_gluon()
						<< " does not fit the number of gluons in the basis "
						<< nq << std::endl;
		}

	}

	// This version of decompose (as opposed to trace type versions)
	// need scalar product matrix
	// If P_spm is not empty, but has wrong size
	if(  P_spm.size() != cb.size() and  !P_spm.empty()  ){
		std::cerr << "Orthogonal_basis::decompose: The size of the scalar product matrix and the basis do not agree." << std::endl;
	assert( 0 );
	}
	// If diagonal_P_spm is not empty, but has wrong size
	if( diagonal_P_spm.size() != cb.size() and !diagonal_P_spm.empty() ){
		std::cerr << "Orthogonal_basis::decompose: The size of the diagonal scalar product matrix and the basis do not agree." << std::endl;
	assert( 0 );
	}
	// If both are empty, calculate diagonal version
	if( P_spm.empty() and diagonal_P_spm.empty() ){
		diagonal_scalar_product_matrix(true, true, true);
	}

	// Use matrix information if diagonal spm is empty
	if( diagonal_P_spm.empty() ){
		for ( uint i=0; i< P_spm.size(); i++ ){
			diagonal_P_spm.push_back( P_spm.at(i).at(i) );
		}
		diagonal_d_spm=Col_fun.double_num(diagonal_P_spm);
	}

	// To contain the decomposed vector
	Poly_vec Decv;

	// Loop over all vectors in cb and calculate projection
	for (uint i = 0; i < cb.size(); i++) {
		double inv_norm=1.0/diagonal_d_spm.at(i);
		//Decv.at(i)= Col_fun.scalar_product( cb.at(i) , Ca )*inv_norm;
		Decv.push_back( Col_fun.scalar_product( cb.at(i) , Ca )*inv_norm );
	}

	return Decv;
}


std::string Orthogonal_basis::diagonal_spm_file_name(const bool leading, const bool poly ) const{

	// First construct filename
	std::ostringstream ss;
	std::string filename;
	ss << "ColorResults";
	ss << '/';
	// CF as in ColorFull
	ss << "CF_";
	// Prefix according to basis type
	ss << 	"OB_";
	// diagonal version
	ss << 	"diagonal_";
	// Polynomial or numerical matrix?
	if( poly )ss << "P_";
	else ss << "d_";
	ss << "spm_q";
	ss << nq;
	ss << "_g";
	ss << ng;
	// is the result leading, and if so, how?
	if ( leading ) ss << "_l";
	if ( Col_fun.get_full_CF() ) ss << "_cff";
	else ss << "_cfl";
	if( Col_fun.get_Nc() != 3 ){
		ss << "_Nc_";
		ss << Col_fun.get_Nc();
	}
	if(Col_fun.get_TR() != 0.5 ){
		ss << "_TR_";
		ss << Col_fun.get_TR();
	}
	filename=ss.str();
	return filename;
}


void Orthogonal_basis::write_out_diagonal_d_spm( std::string filename ) const{

	Col_fun.write_out_dvec( diagonal_d_spm, filename );

}


void Orthogonal_basis::write_out_diagonal_d_spm( ) const{

	std::string filename = diagonal_spm_file_name( false, false );
	write_out_diagonal_d_spm( filename );
}


void Orthogonal_basis::write_out_diagonal_P_spm( std::string filename ) const{

	diagonal_P_spm.write_out_Poly_vec( filename );
}


void Orthogonal_basis::write_out_diagonal_P_spm( ) const{

	std::string filename = diagonal_spm_file_name( false, true );
	write_out_diagonal_P_spm( filename );
}


Polynomial Orthogonal_basis::scalar_product( const Col_amp & Ca1, const Col_amp & Ca2 )  {

	// Check that we have a basis
	if(cb.size()==0){
		std::cerr << "Orthogonal_basis::scalar_product: The basis vector cb is empty consider using create_basis or read_in_basis." << std::endl;
		assert( 0 );
	}

	// Check that size of P_spm is consistent (if calculated)
	if( (P_spm.size() != cb.size()) and  P_spm.size() !=0 ) {
		std::cerr << "Orthogonal_basis::scalar_product: Size of scalar product matrix P_spm and color basis cb do not agree." << std::endl;
		assert( 0 );
	}

	// Check that size of diagonal_P_spm is consistent (if calculated)
	if( (diagonal_P_spm.size() != cb.size()) and  diagonal_P_spm.size() !=0 ) {
		std::cerr << "Orthogonal_basis::scalar_product: Size of diagonal_P_spm and color basis cb do not agree." << std::endl;
		assert( 0 );
	}

	// Check if at least one of P_spm and diagonal_P_spm exist,
	// If non exist, calculate diagonal
	if( P_spm.empty() and diagonal_P_spm.empty() ){
		diagonal_scalar_product_matrix(true, true, true);
	}
	// If matrix form, but not diagonal exist, copy diagonal entries
	else if(diagonal_P_spm.empty() ){
		for( uint i=0; i< P_spm.size(); i++ ){
			diagonal_P_spm.push_back(P_spm.at(i).at(i));
		}
	}

	// To contain the resulting Polynomial
	Polynomial Poly_res;
	Poly_res=Poly_res*0;

	// Decompose the Col_amps
	Poly_vec Polyv1=decompose(Ca1);
	Polyv1.conjugate();
	Poly_vec Polyv2=decompose( Ca2 );

	// Then add contributions
	for (uint m1=0; m1< cb.size(); m1++){
		// Diagonal terms
		Poly_res=Poly_res+Polyv1.at(m1) *Polyv2.at(m1) *diagonal_P_spm.at(m1);
	}
	return Poly_res;
}


cnum Orthogonal_basis::scalar_product_num( const Col_amp & Ca1, const Col_amp & Ca2 )  {

	// Check that we have a basis
	if(cb.size()==0){
		std::cerr << "Orthogonal_basis::scalar_product_num: The basis vector cb is empty consider using create_basis or read_in_basis." << std::endl;
		assert( 0 );
	}

	// Check that size of d_spm is consistent (if calculated)
	if( (d_spm.size() != cb.size()) and  d_spm.size() !=0 ) {
		std::cerr << "Orthogonal_basis::scalar_product_num: Size of scalar product matrix d_spm and color basis cb do not agree." << std::endl;
		assert( 0 );
	}

	// Check that size of diagonal_d_spm is consistent (if calculated)
	if( (diagonal_d_spm.size() != cb.size()) and  diagonal_d_spm.size() !=0 ) {
		std::cerr << "Orthogonal_basis::scalar_product_num: Size of diagonal_d_spm and color basis cb do not agree." << std::endl;
		assert( 0 );
	}

	// Check if at least one of d_spm and diagonal_d_spm exist,
	// If non exist, calculate diagonal version
	if( d_spm.empty() and diagonal_d_spm.empty() ){
		diagonal_scalar_product_matrix(false, true, true);
	}
	// If matrix form, but not diagonal exist, copy diagonal entries
	else if(diagonal_d_spm.empty() ){
		for( uint i=0; i< d_spm.size(); i++ ){
			diagonal_d_spm.push_back(d_spm.at(i).at(i));
		}
	}

	// To contain the result
	double res=0;

	// Decompose the Col_amps
	Poly_vec Polyv1=decompose(Ca1);
	Polyv1.conjugate();
	Poly_vec Polyv2=decompose( Ca2 );
	dvec v1=Col_fun.double_num( Polyv1 );
	dvec v2=Col_fun.double_num( Polyv2 );


	// Then add contributions
	for (uint m1=0; m1< cb.size(); m1++){
		// Diagonal terms
		res=res+v1.at(m1) *v2.at(m1) *diagonal_d_spm.at(m1);
	}
	return res;
}


cnum Orthogonal_basis::scalar_product_num( const cvec & v1, const cvec & v2) {

	if(v1.size()!= v2.size()){
		std::cerr << "Orthogonal_basis::scalar_product_num: Size of first vector "
		<< v1.size() << " does not agree with size of second vector "
		<< v2.size() << std::endl;
		assert( 0 );
	}

	// Check that size of d_spm is consistent (if calculated)
	if( ( d_spm.size() != cb.size()) and  d_spm.size() !=0 ) {
		std::cerr << "Orthogonal_basis::scalar_product_num: Size of scalar product matrix d_spm and color basis cb do not agree." << std::endl;
		assert( 0 );
	}

	// Check that size of diagonal_d_spm is consistent (if calculated)
	if( ( diagonal_d_spm.size() != cb.size()) and  diagonal_d_spm.size() !=0 ) {
		std::cerr << "Orthogonal_basis::scalar_product_num: Size of diagonal_d_spm and color basis cb do not agree." << std::endl;
		assert( 0 );
	}

	// Check if at least one of d_spm and diagonal_d_spm exist,
	// If non exist, calculate diagonal version
	if( d_spm.empty() and diagonal_d_spm.empty() ){
		diagonal_scalar_product_matrix(false, true, true);
	}
	// If matrix form, but not diagonal exist, copy diagonal entries
	else if(diagonal_d_spm.empty() ){
		for( uint i=0; i< d_spm.size(); i++ ){
			diagonal_d_spm.push_back(d_spm.at(i).at(i));
		}
	}

	uint basis_size=v1.size();

	// To contain the result
	cnum res=0;

	// Add contributions, diagonal parts only
	for (uint m1=0; m1< basis_size; m1++){
		res += conj( v1.at(m1) ) *v2.at(m1) *diagonal_d_spm.at(m1);
	}

	return res;
}


void Orthogonal_basis::write_out_diagonal_spm( const dvec & dv, const bool leading ) const{

	std::string filename = diagonal_spm_file_name( leading, false );
	std::ofstream outfile(filename.c_str());

	if ( !outfile )
	std::cerr << "Orthogonal_basis::write_out_diagonal_spm: Cannot write out diagonal scalar products as the file \""
		<< filename.c_str() << "\" could not be opened. (Does the directory exist? Consider creating the directory.)" << std::endl;


	outfile << dv;
}


void Orthogonal_basis::write_out_diagonal_spm( const Poly_vec & Pv, const bool leading ) const{

	std::string filename = diagonal_spm_file_name( leading, true );
	std::ofstream outfile(filename.c_str());
	outfile << Pv;
}


} //end namespace ColorFull
