// -*- C++ -*-
/*
 * Col_basis.cc
 * Contains definition of the base class Col_basis and associated types and operators.
 * Created on: Aug 9, 2012
 * Author: Malin Sjodahl
 */

#include "Col_basis.h"
#include "parameters.h"
#include <cassert>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <iomanip>


namespace ColorFull {


void Col_basis::append( col_basis cb_in ) {
	for (uint m = 0; m < cb_in.size(); m++) {
		cb.push_back( cb_in.at(m) );
	}
}

void Col_basis::scalar_product_matrix_no_mem(){
	if(ng+nq>6){
		std::cout << "Col_basis::scalar_product_matrix: nq+ng=" << nq+ng << " is large, consider using numerical and/or memory version.  "  << std::endl;
		std::cout.flush();
	}
	return scalar_product_matrix( true, true, false );
}


void Col_basis::scalar_product_matrix(){
	return scalar_product_matrix( true, true, true );
}


void Col_basis::scalar_product_matrix_num_no_mem(){
	return scalar_product_matrix( false, true, false );
}


void Col_basis::scalar_product_matrix_num(){
	return scalar_product_matrix( false, true, true );
}

void Col_basis::leading_scalar_product_matrix(){

	if( cb.size()==0 ) {
		std::cerr << "Col_basis::leading_scalar_product_matrix: There are no basis vectors in this basis, consider using create_basis or read_in_Col_basis." << std::endl;
		std::cerr.flush();
		return ;
	}

	if( P_spm.empty() ) {
		scalar_product_matrix( true, true, true );
	}

	leading_P_spm = Col_fun.leading( P_spm );
	leading_d_spm = Col_fun.double_num( leading_P_spm );

}


std::string Col_basis::basis_file_name() const{

	// First construct filename
	std::ostringstream ss;
	std::string filename;
	ss << "ColorResults";
	ss << '/';
	// ColorFull
	ss << "CF_";
	// Basis Type
	if( trace_basis )  ss << "TB_";
	else if ( tree_level_gluon_basis ) ss << "TGB_";
	else if ( orthogonal_basis ) ss << "OB_";
	else  ss << "CB_";
	ss << "q_";
	ss << nq;
	ss << "_g_";
	ss << ng;
	// If Nc is not 3, append Nc info
	if( Col_fun.get_Nc() != 3 ){
		ss << "_Nc_";
		ss << Col_fun.get_Nc();
	}
	// If TR is not 1/2, append TR info
	if( Col_fun.get_TR() != 0.5 ){
		ss << "_TR_";
		ss << Col_fun.get_TR();
	}

	filename=ss.str();
	return filename.c_str();
}


std::string Col_basis::spm_file_name(const bool leading, const bool poly ) const{

	// First construct filename
	std::ostringstream ss;
	std::string filename;
	ss << "ColorResults";
	ss << '/';
	// CF as in ColorFull
	ss << "CF_";
	// Prefix according to basis type
	if( trace_basis ) ss << "TB_";
	else if ( tree_level_gluon_basis ) ss << 	"TGB_";
	else if ( orthogonal_basis ) ss << 	"OB_";
	else ss << "CB_";
	// Polynomial or numerical matrix?
	if( poly )ss << "P_";
	else ss << "d_";
	ss << "spm_q";
	ss << nq;
	ss << "_g";
	ss << ng;
	if ( leading ) ss << "_l";
	if ( Col_fun.get_full_CF() ) ss << "_cff";
	else ss << "_cfl";

	if(Col_fun.get_Nc() != 3 ){
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


void Col_basis::write_out_Col_basis() const{
	write_out_Col_basis( basis_file_name() );
}


void Col_basis::write_out_Col_basis( std::string filename ) const {

	if ((cb.size() == 0)) {
		std::cerr
		<< "Col_basis::write_out_Col_basis(filename): There are no basis vectors in this basis, consider using create_basis or read_in_Col_basis."
		<< std::endl;
		return;
	}

	std::ofstream outfile(filename.c_str());

	if ( !outfile )
	std::cerr << "Col_basis::write_out_Col_basis: Cannot write out basis as the file \""
		<< filename.c_str() << "\" could not be opened. (Does the directory exist? Consider creating the directory.)" << std::endl;


	outfile << std::setprecision(16);

	for ( uint m = 0; m < cb.size(); m++ ) {
		outfile << m << "      " << cb.at(m) << std::endl;
	}
	outfile.flush();
}

std::ostream& Col_basis::write_out_Col_basis_to_stream( std::ostream& out ) const{

	if(  (cb.size()==0 ) ) {
		std::cerr << "Col_basis::write_out_Col_basis(): There are no basis vectors in this basis, consider using create_basis or read_in_Col_basis." << std::endl;
		std::cerr.flush();
	}

	for (uint m = 0; m < cb.size(); m++) {
		out << m << "      "<< cb.at(m) << std::endl;
	}
	return out;
}

void Col_basis::write_out_d_spm( std::string filename ) const{

	Col_fun.write_out_dmatr( d_spm, filename );

}


void Col_basis::write_out_d_spm( ) const{

	std::string filename = spm_file_name( false, false );
	write_out_d_spm( filename );
}


void Col_basis::write_out_P_spm( std::string filename ) const{

	P_spm.write_out_Poly_matr( filename );
}


void Col_basis::write_out_P_spm( ) const{

	std::string filename = spm_file_name( false, true );
	write_out_P_spm( filename );
}


void Col_basis::write_out_leading_d_spm( std::string filename ) const{

	Col_fun.write_out_dmatr( leading_d_spm, filename );
}


void Col_basis::write_out_leading_d_spm( ) const{

	std::string filename = spm_file_name( true, false );
	write_out_leading_d_spm( filename );
}


void Col_basis::write_out_leading_P_spm( std::string filename ) const{

	leading_P_spm.write_out_Poly_matr( filename );
}


void Col_basis::write_out_leading_P_spm( ) const{

	std::string filename = spm_file_name( true, true );

	write_out_leading_P_spm( filename );
}


void Col_basis::read_in_Col_basis( std::string filename ) {

	// Read in file
	std::ifstream fin(filename.c_str());

	// Check that file exists
	if( !fin ){
		std::cerr << "Col_basis::read_in_Col_basis: The file "
				<< filename << " could not be opened." << std::endl;
		assert( 0 );
	}

	// Erase current information
	cb.clear();

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

	Col_basis_of_str( str );

}


void Col_basis::read_in_d_spm( std::string filename){

	d_spm= Col_fun.read_in_dmatr( filename );

}


void Col_basis::read_in_d_spm( ){

	d_spm=Col_fun.read_in_dmatr( spm_file_name( false, false ).c_str() );

}


void Col_basis::read_in_leading_d_spm( std::string filename){

	leading_d_spm=Col_fun.read_in_dmatr( filename );

}


void Col_basis::read_in_leading_d_spm( ){

	leading_d_spm=Col_fun.read_in_dmatr( spm_file_name( true, false ).c_str() );
}


void Col_basis::read_in_P_spm( std::string filename ){

	P_spm.read_in_Poly_matr( filename );
}


void Col_basis::read_in_P_spm( ){

	P_spm.read_in_Poly_matr( spm_file_name( false, true) );

}


void Col_basis::read_in_leading_P_spm( std::string filename){

	leading_P_spm.read_in_Poly_matr( filename );

}


void Col_basis::read_in_leading_P_spm( ){

	leading_P_spm.read_in_Poly_matr( spm_file_name( true, true) );
}


int Col_basis::n_quark_check() const {

	if (empty()) return 0;
	int nq=cb.at(0).n_quark_check();
	for( uint n=0; n< cb.size(); n++ )
		if( cb.at(n).n_quark_check() != nq) {
			std::cerr << "Col_basis::n_quark_check: The Col_amps in " << cb << " have differently many quarks." << std::endl;
		}

	return nq;
}


int Col_basis::n_gluon_check() const {

	if ( empty()) return 0;
	int nq=cb.at(0).n_gluon_check();
	for( uint n=0; n< cb.size(); n++ )
		if( cb.at(n).n_gluon_check() != nq) {
			std::cerr << "Col_basis::n_gluon_check: The Col_amps in " << cb << " have differently many gluons." << std::endl;
		}
	return nq;
}


void Col_basis::simplify(){

	for( uint i=0; i< cb.size(); i++ ){
		cb.at(i).simplify();
	}
}


Poly_vec Col_basis::decompose( const Col_amp & Ca ){

	// This line is to avoid compiler warnings
	// Ca should be an argument despite that it's not used
	(void) Ca;

	std::cerr << "Col_basis::decompose: This function is not implemented for the Col_basis class. Try using a derived class (such as Trace_basis). " << std::endl;

	Poly_vec Pv;
	assert( 0 );

	return Pv;
}


Col_amp Col_basis::exchange_gluon( uint vec, int p1, int p2 ) {

	if(  (cb.size()==0 ) ) {
		std::cerr << "Col_basis::exchange_gluon: There are no basis vectors in this basis, consider using create_basis or read_in_Col_basis." << std::endl;
		std::cerr.flush();
		assert( 0 );
	}

	// Check that the basis vector exists
	if ( vec >= cb.size() ) {
		std::cerr << "Col_basis::exchange_gluon: Basis vector number "<< vec
				<< " does not exist, as the basis only have " << cb.size()
				<< " basis vectors." << std::endl;
		assert( 0 );
	}

	return Col_fun.exchange_gluon( cb.at(vec), p1, p2);

}

Poly_matr Col_basis::color_gamma( int p1, int p2 ) {

	if(  (cb.size()==0 ) ) {
		std::cerr << "Col_basis::color_gamma: There are no basis vectors in this basis, consider using create_basis or read_in_Col_basis." << std::endl;
		std::cerr.flush();
		assert( 0 );
	}

	// Function not available for Tree_level_gluon basis
	if( tree_level_gluon_basis ) {
		std::cerr << "Col_basis::color_gamma: This function is not available for Tree_level_gluon_basis "
				<< "as the vector resulting after gluon exchange contains vectors which are not in the basis. "
				<< "Consider using Trace_basis."<< std::endl;
		std::cerr.flush();
		assert( 0 );
	}

	// Function should only be used for Orthogonal_basis and Trace_basis
	if( ( !trace_basis ) and (! orthogonal_basis) ) {
		std::cerr << "Col_basis::color_gamma: This function is only implemented for Trace_basis and Orthogonal_basis. "
				<< "If your basis is not orthogonal, the result will not be correct. "
				<< "If your basis is orthogonal, consider using Orthogonal_basis."
				<< std::endl;
		std::cerr.flush();
	}

	// To contain the resulting matrix
	Poly_matr gamma_res;
	// Fill gamma_res with 0s
	Polynomial Zero;
	Zero=Zero*0;
	for( uint i=0; i < cb.size(); i++ ){
		Poly_vec rowi;
		for( uint j=0; j < cb.size(); j++ ){
			rowi.append(Zero);
		}
		gamma_res.append(rowi);
	}

	// To contain the Col_amp after gluon exchange
	Col_amp Ca_ae;

	// First exchange a gluon between partons p1 and p2 in the vector v1
	for ( uint vi=0; vi < cb.size(); vi++ ){

		Ca_ae=exchange_gluon( vi, p1, p2 );

		// Then decompose the result into the basis,
		Poly_vec col_vi = decompose( Ca_ae );

		// This gives the vi:th column in color_gamma
		// Loop over rows in that column
		for ( uint vj=0; vj < cb.size(); vj++ ){
			gamma_res.at(vj).pv.at(vi) = col_vi.at(vj);
		}
	}


	return gamma_res;
}


void Col_basis::Col_basis_of_str( std::string str ){

	if (str.size() == 0) {
		std::cerr << "Col_basis::Col_basis_of_str: The basis string is empty. This should not happen." << std::endl;
		assert( 0 );
	}

	// Skip lines starting with #
	while(str.at(0)== '#'){
		while (str.at(0) != '\n'){
			str.erase( str.begin() );
		}
		// erase endl sign(s)
		while(str.at(0)== '\n'){
			str.erase(str.begin());
		}
	}

	// First char in string should be '0'
	if (str.at(0) != '0') {
		std::cerr
		<< "Col_basis::Col_basis_of_str: First char in basis data file after comments should be '0', as in basis vector 0 as in vector number 0, but it was: "
		<< str.at(0) << std::endl;
		assert( 0 );
	}

	// Read the string, starting from 0th element
	uint i = 0;
	while (i < str.size() - 2) {
		//i += 1; testing to move down 130917

		// We may have some white spaces
		while (i< str.size()-2 &&(str.at(i) == ' ')) i++;

		// After white spaces there should be a number, the vector number
		while (i< str.size()-2 &&(str.at(i) == '0' or str.at(i) == '1' or str.at(i) == '2' or str.at(i) == '3' or str.at(i) == '4' or str.at(i) == '5' or str.at(i) == '6' or str.at(i) == '7' or str.at(i) == '8' or str.at(i) == '9') ) i++;

		// We may have some more white spaces
		while (i< str.size()-2 &&(str.at(i) == ' ')) i++;


		// String to make a Col_amp of
		std::string Ca_str;
		Ca_str.clear();

		// Keep reading the Ca_str while not '\n which marks new vector
		//while ( i< str.size()-1 && (str.at(i) != '\n') )
		while ( i< str.size() && (str.at(i) != '\n') ){
			Ca_str.push_back(str.at(i));
			i++;
		}

		// Make the Col_amp (basis vector)
		Col_amp Ca( Ca_str );

		// Simplify, for example move Polynomials to multiply Col_strs, rather than Quark-lines
		Ca.simplify();

		append(Ca);
		i++;
	}

	// Check that the basis is not empty
	if( empty()){
		std::cerr << "Col_basis::Col_basis_of_str: The basis is empty. " << std::endl;
		assert( 0 );
	}

	// Check that first vector is not empty
	if(cb.at(0).ca.empty()){
		std::cerr << "Col_basis::Col_basis_of_str: The first Col_amp (vector) in the basis is empty. " << std::endl;
		assert( 0 );
	}

	// Check that first col_str in first vector is not empty
	if(cb.at(0).ca.at(0).cs.empty()){
		std::cerr << "Col_basis::Col_basis_of_str: The first Col_str in the first basis vector is empty. " << std::endl;
		assert( 0 );
	}

	// Check that first quark_line in the first col_str in first vector is not empty
	if(cb.at(0).ca.at(0).cs.at(0).ql.empty()){
		std::cerr << "Col_basis::Col_basis_of_str: The first Quark_line in the first Col_str in the first basis vector is empty. " << std::endl;
		assert( 0 );
	}

	// Find the number of quarks and gluons
	nq=n_quark_check();
	ng=n_gluon_check();

	// Test if it's a trace basis
	if( trace_basis ){
		bool is_trace_basis=true;
		for ( uint j=0; j< cb.size(); j++ ){
			if( cb.at(j).size() !=1 ) is_trace_basis=false;

			if(! is_trace_basis ){
				std::cerr << "Col_basis::Col_basis_of_str: You are trying to read in a non-trace basis to a trace basis object. "
						<<"The length of the Col_amp " << cb.at(j) << " is not 1."
						<< std::endl;
				assert( 0 );
			}
		}
	}
}


void Col_basis::read_in_Col_basis( ) {

	read_in_Col_basis( basis_file_name() );
}


void Col_basis::scalar_product_matrix( bool save_P_spm, bool save_d_spm, bool use_mem ) {

	if(  (cb.size()==0 ) ) {
		std::cerr << "Col_basis::scalar_product_matrix: There are no basis vectors in this basis, consider using create_basis or read_in_Col_basis." << std::endl;
		std::cerr.flush();
		return ;
	}

	// Simplify the basis vectors
	simplify();

	// If the P_spm and d_spm have already been calculated, erase them
	if( !P_spm.empty() ) P_spm.clear();
	if( !d_spm.empty() ) d_spm.clear();

	// Check that all Polynomials are real
	// (check doesn't cover the case of complex Poly in a Quark_line)
	// but if an entry is not real this will be discovered as the d_spm is calculated
	for ( uint cbi=0; cbi < cb.size(); cbi++){
		for ( uint j=0; j< cb.at(cbi).size(); j++){
			if( imag (Col_fun.cnum_num(  cb.at(cbi).at(j).Poly  )) > accuracy ){
				std::cerr << "Col_basis::scalar_product_matrix: ColorFull expects real Polynomial multiplying the color structure, but the Polynomial\n"
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

		// To contain a row of the resulting scalar product matrix
		Poly_vec rowi;

		// Loop over basis vectors in Basis, normally over all vectors,
		// but for orthogonal_basis only for i==j
		uint minj=0, maxj=cb.size();

		for( uint j=minj; j< maxj; j++){
			Polynomial ijRes;
			ijRes=ijRes*0;
			if(use_mem){

				// Loop over Col_strs in first Ca
				for( uint Ca1i=0; Ca1i< cb.at(i).size(); Ca1i++){

					// Loop over Col_strs in second Ca
					for( uint Ca2i=0; Ca2i< cb.at(j).size(); Ca2i++){

						// To contain the contribution to the ij-th entry if memoization is used
						std::shared_ptr<Polynomial> ijEntry_contr;

						// Rename indices, and make string of new col_strs, to use in map
						// strip off Polynomial information to make the map minimal
						Col_str Cs1, Cs2;
						Cs1.append(cb.at(i).at(Ca1i).cs);
						Cs2.append(cb.at(j).at(Ca2i).cs);

						rename_indices( Cs1, Cs2 );

						std::ostringstream Cs_string;
						Cs_string << Cs1 << Cs2;

						// Calculate element contribution to the ij:th element in scalar product matrix
						// If this has scalar product has occurred before, reuse old value
						if (mem_map.count( Cs_string.str() ) > 0) {
							ijEntry_contr = mem_map[Cs_string.str()];
						}
						// Otherwise, calculate value and save topology
						else {
							if( !tree_level_gluon_basis ){
								Polynomial p = Col_fun.scalar_product(Cs1, Cs2);
								ijEntry_contr = std::shared_ptr<Polynomial> (new Polynomial(p));
								mem_map[Cs_string.str()] = ijEntry_contr;

							}
							else if( tree_level_gluon_basis ){
								int sign = (ng % 2 ? -1 : 1);
								Col_str Cs2_conj=Cs2;
								Cs2_conj.conjugate();
								Polynomial P = 2* Col_fun. scalar_product( Cs1, Cs2 )
										+ sign*2*Col_fun. scalar_product( Cs1, Cs2_conj );
								ijEntry_contr = std::shared_ptr<Polynomial> (new Polynomial(P));
								mem_map[Cs_string.str()] = ijEntry_contr;
							}
						}
						// Sum up all the contributions to one Polynomial, recall to multiply with Polynomials
						Polynomial ijEntry_contr_poly=(*ijEntry_contr)* cb.at(i).at(Ca1i).Poly*cb.at(j).at(Ca2i).Poly;

						// If Polynomial version is not needed convert all algebraic information
						// to numerical in the vectors
						if( !save_P_spm ) {
							double ijEntry_contr_d= Col_fun.double_num( ijEntry_contr_poly );
							Monomial Mon( ijEntry_contr_d );
							ijRes+= Mon;
						}
						else ijRes += ijEntry_contr_poly;

					} // end looping over Cs in first Ca
				} // end looping over Cs in 2nd Ca

				// If Polynomial result is wanted, simplify
				if ( save_P_spm ) ijRes.simplify();

				// Otherwise convert to numerical
				else {
					double num_ijRes= Col_fun.double_num(ijRes);
					ijRes.clear();
					Monomial Mon( num_ijRes );
					ijRes.append(Mon);
				}

				rowi.append( ijRes );

			} // end if ( use_mem )
			else{
				// Calculate element ij in scalar product matrix
				ijRes=ij_entry( i, j );
				rowi.append( ijRes );

				// Convert to numerical if not saving Polynomial version
				if ( ! save_P_spm ) {
					double num_ijRes= Col_fun.double_num(ijRes);
					ijRes.clear();
					Monomial Mon( num_ijRes );
					ijRes.append(Mon);
				}

			} //end if (not mem)
		}
		P_spm.append(rowi);
	}

	// For saving and checking symmetry of numerical version
	// d_spm has to be calculated even if it's not saved to facilitate
	// symmetry check
	d_spm= Col_fun.double_num( P_spm );

	// Making consistency checks (needs double version)
	check_spm();

	// Erasing the double spm
	if( ! save_d_spm ) {
		d_spm.clear();
	}
	// Erasing the Polynomial spm
	if( !save_P_spm ) {
		P_spm.clear();
	}

	return;
}


Polynomial Col_basis::ij_entry( const int i, const int j ) const{

	Polynomial ijEntry;
	ijEntry=Col_fun.scalar_product(cb.at(i), cb.at(j));
	ijEntry.simplify();
	return ijEntry;
}


bool Col_basis::check_symmetry( const dmatr & matr ) const {

	if( matr.size()==0 ){
		std::cout << "Col_basis::check_symmetry( dmatr ): The numerical matrix is empty..." << std::endl;
	}
	else{
		//std::cout << "Col_basis::check_symmetry( dmatr ): Numerically verifying symmetry of matrix...";
	}
	// Verifying that the matrix is symmetric to accuracy "accuracy",
	// if not sym is put to false
	bool sym=true;
	// Loop over basis vectors in Basis
	for (uint i = 0; i < d_spm.size(); i++) {
		Poly_vec rowi;
		// Loop over basis vectors in Basis
		for (uint j = 0; j <=i; j++) {
			if ((fabs(d_spm.at(i).at(j) / d_spm.at(j).at(i) - 1.0) > accuracy) && (d_spm.at(i).at(j) > accuracy) && (d_spm.at(j).at(i) > accuracy)) {
				sym=false;
				std::cerr
				<< "Col_basis::check_symmetry( dmatr ): Error, the resulting scalar product matrix is not symmetric. \n "
				<< "Element " << i << "," << j << ": "
				<< d_spm.at(i).at(j) << ", Element " << j << "," << i
				<< ": " << d_spm.at(j).at(i) << std::endl
				<< "This indicates an error in calculation of scalar products. " << std::endl;
			}
		}
	}
	//if( matr.size() !=0  and sym ) std::cout << " done." << std::endl;
	return sym;
}


bool Col_basis::check_diagonal( const dmatr & matr ) const{

	if( matr.size()==0 ){
		std::cerr << "Col_basis::check_diagonal( dmatr ): The numerical matrix is empty...";
		std::cout.flush();
	}
	else{
		//std::cout << "Col_basis::check_diagonal( dmatr ): Numerically checking if the matrix is diagonal...";
		//std::cout.flush();
	}

	// Verifying that the matrix is diagonal
	// if not sym is put to false
	bool diag=true;
	// Loop over basis vectors in Basis
	for (uint i = 0; i < matr.size(); i++) {

		Poly_vec rowi;
		// Loop over basis vectors in Basis
		for (uint j = 0; j <=i; j++) {
			// If non-diagonal elements are sufficiently large write warning
			if ( std::abs( matr.at(i).at(j) ) > accuracy and i!=j ) {
				diag=false;
				if( !(trace_basis or tree_level_gluon_basis) ){
					std::cout
					<< "Col_basis::check_diagonal( matr ): Warning, the matrix is not diagonal. \n "
					<< "Element " << i << "," << j << ": "
					<< matr.at(i).at(j) << std::endl;
				}
			}
		}
	}
	if (! diag ){
		std::cout << "Col_basis::check_diagonal: the matrix is not diagonal." << std::endl;
	}
	//else if (diag) std::cout << " done." << std::endl;

	return diag;
}


void Col_basis::check_spm() const{

	// Verifying that the matrix is symmetric to accuracy
	bool sym=check_symmetry(d_spm);

	if( !sym ) {
		std::cerr <<"Col_basis::check_spm(): scalar product matrix not symmetric. Please report bug." << std::endl;
		std::cerr.flush();
		assert( 0 );
	}

	// Verifying that the leading terms only sit on the diagonal
	bool diagonal=true;
	if( !leading_d_spm.empty() ){
		diagonal=check_diagonal(leading_d_spm);
	}
	if(!diagonal and ( trace_basis or tree_level_gluon_basis ) ) {
		std::cout <<"Col_basis::check_spm(): Leading terms appear of the diagonal. This should not happen in a trace type basis."
				<< " For numerical bases it can appear to happen as powers of Nc may hide in numerical constants." << std::endl;
		std::cout.flush();
	}
}


void Col_basis::rename_indices(Col_str & Cs1, Col_str & Cs2) const{

	/*
     This is a two-step process:
     1) Compute new indices
     2) Replace indices

     The new indices are stored in a vector named new_indices, such that
     old_index should be replaced by new_indices[old_index]
	 */

	std::vector<int> new_indices;
	uint n_indices_total = 2*Cs1.n_quark() + Cs1.n_gluon();
	new_indices.resize(n_indices_total+1);
	assert(2*Cs1.n_quark() + Cs1.n_gluon() == 2*Cs1.n_quark() + Cs1.n_gluon());
	int new_ind=0; // The new index
	uint Nql= Cs1.cs.size();
	// Loop over Quark_lines in Cs1
	for(uint i=0; i< Nql; i++ ){
		Quark_line & Ql_i = Cs1.cs.at(i);
		uint Nind = Ql_i.ql.size();
		// Loop over indices in the Quark_line
		for(uint j=0; j< Nind; j++ ){
			new_ind++;
			// The old index in Cs1, to be replaced also in Cs2
			uint old_ind=Ql_i.ql.at(j);
			// Replace index in Cs1
			Ql_i.ql.at(j) = new_ind;
			// Insert into replacement map
			// If the assert fails, we have not allocated enough space, and will fail.
			assert(old_ind <= n_indices_total);
			new_indices[old_ind] = new_ind;
		}
	}
	// Loop over Quark_lines in Cs2
	Nql= Cs2.cs.size();
	for(uint i=0; i< Nql; i++ ){
		Quark_line & Ql_i = Cs2.cs.at(i);
		uint Nind = Ql_i.ql.size();
		// Loop over indices in the Quark_line
		for(uint j=0; j< Nind; j++ ){
			// Replace index in Cs2
			Ql_i.ql.at(j) = new_indices[Ql_i.ql.at(j)];
		}
	}
}


Polynomial Col_basis::scalar_product( const Col_amp & Ca1, const Col_amp & Ca2 )  {


	// Check that we have a basis
	if(cb.size()==0){
		std::cerr << "Col_basis::scalar_product: The basis vector cb is empty, consider using create_basis or read in basis." << std::endl;
		assert( 0 );
	}

	// Check that the Polynomial scalar product matrix is calculated, if not calculate it
	if( P_spm.size() != cb.size() ) {
		std::cerr << "Col_basis::scalar_product: This function uses the scalar product matrix which has not yet been calculated." << std::endl;
		assert( 0 );
	}

	// To contain the resulting Polynomial
	Polynomial Poly_res;
	Poly_res=Poly_res*0;

	// Decompose the Col_amps
	Poly_vec Polyv1=decompose(Ca1);
	Polyv1.conjugate();
	Poly_vec Polyv2=decompose( Ca2 );

	// Then add contributions
	for ( uint m1=0; m1< cb.size(); m1++ ){
		// Diagonal terms
		Poly_res+= Polyv1.at(m1) *Polyv2.at(m1) *P_spm.at(m1).at(m1);
		for (uint m2=0; m2< m1; m2++){
			// Other terms, use symmetry of scalar product matrix
			Poly_res+= ( Polyv1.at(m1) *Polyv2.at(m2)+Polyv1.at(m2) *Polyv2.at(m1) ) *P_spm.at(m1).at(m2);
		}
	}
	return Poly_res;
}


cnum Col_basis::scalar_product_num( const Col_amp & Ca1, const Col_amp & Ca2 )  {


	// Check that we have a basis
	if(cb.size()==0){
		std::cerr << "Col_basis::scalar_product_num: The basis vector cb is empty consider using create_basis or read in basis." << std::endl;
		assert( 0 );
	}

	// Check that the Polynomial scalar product matrix is calculated, if not calculate it
	if( d_spm.size() != cb.size() ) {
		std::cerr << "Col_basis::scalar_product_num: This function uses the numerical scalar product matrix which has not yet been calculated." << std::endl;
		assert( 0 );
	}
	// To contain the resulting Polynomial
	cnum res=0;
	uint basis_size=cb.size();

	// Decompose the Col_amps
	Poly_vec Polyv1, Polyv2;
	Polyv1=decompose(Ca1);
	if( &Ca1==&Ca2 ) Polyv2=Polyv1;
	else Polyv2=decompose(Ca2);

	// Conjugate first arg
	Polyv1.conjugate();

	// Make ordinary vectors to speed up multiplication
	cvec v1, v2;
	v1.reserve( basis_size );
	for (uint i=0; i< basis_size; i++){
		v1.push_back( Col_fun.cnum_num( Polyv1.at(i) ) );
	}

	if (false) v2=v1;
	else{
		v2.reserve(basis_size);
		for (uint i=0; i< basis_size; i++) v2.push_back( Col_fun.cnum_num( Polyv2.at(i) ) );
	}

	// Then add contributions
	cnum v1m1, v2m1;
	dvec row;
	for (uint m1=0; m1< basis_size; m1++){
		v1m1=v1.at(m1);
		v2m1=v2.at(m1);
		row=d_spm.at( m1 );
		// Diagonal terms
		res=res+ v1m1 *v2.at(m1) *row.at(m1);
		for (uint m2=0; m2< m1; m2++){
			res += (v1m1*v2.at(m2) +v1.at(m2)*v2m1 ) *row.at(m2);

		}
	}

	return res;

}


cnum Col_basis::scalar_product_num( const cvec & v1, const cvec & v2) {

	if( v1.size()!= v2.size() ){
		std::cerr << "Col_basis::scalar_product_num: Size of first vector "
				<< v1.size() << " does not agree with size of second vector "
				<< v2.size() << std::endl;
		assert( 0 );
	}
	if( v1.size()!= d_spm.size() ){
		std::cerr << "Col_basis::scalar_product_num: Size of vectors "
				<< v1.size() << " does not agree with size of d_spm matrix "
				<< d_spm.size() << std::endl;
		assert( 0 );
	}

	uint basis_size=v1.size();

	// To contain the result
	cnum res=0;
	cnum tmp=0, tmp2=0;

	// Add contributions
	cnum v1m1, v2m1;
	for (uint m1=0; m1< basis_size; m1++){
		v1m1= conj( v1.at(m1) );
		v2m1=v2.at(m1);
		const dvec &row = d_spm.at(m1);
		// Diagonal parts
		res += v1m1 *v2.at(m1) *row.at(m1);
		tmp = 0;
		tmp2=0;
		for (uint m2=0; m2< m1; m2++){
			tmp += (v1m1*v2[m2]+v2m1*conj(v1[m2])) *row[m2];
		}
		res+=tmp;
	}

	return res;
}


cnum Col_basis::scalar_product_num_diagonal( const cvec & v1, const cvec & v2 ) {

	if(v1.size()!= v2.size()){
		std::cerr << "Col_basis::scalar_product_num_diagonal: Size of first vector "
				<< v1.size() << " does not agree with size of second vector "
				<< v2.size() << std::endl;
		assert( 0 );
	}
	if(v1.size()!= d_spm.size()){
		std::cerr << "Col_basis::scalar_product_num_diagonal: Size of vectors "
				<< v1.size() << " do not agree with size of d_spm matrix "
				<< d_spm.size() << std::endl;
		assert( 0 );
	}

	// To contain the result
	cnum res=0;

	// Check dimension
	uint basis_size=v1.size();

	// Then add contributions
	for (uint m1=0; m1< basis_size; m1++){
		res=res+ conj( v1.at(m1) ) *v2.at(m1) *d_spm.at(m1).at(m1);
	}

	return res;
}


std::ostream& operator<<( std::ostream& out, const Col_basis & Cb ) {

	return Cb.write_out_Col_basis_to_stream( out );
}


std::ostream& operator<<( std::ostream& out, const col_basis & cb ) {
	int max = cb.size();
	for (int i = 0; i < max; i++) {
		out <<",  " << cb.at(i);
	}
	return out;
}

}
