// -*- C++ -*-
/*
 * Tree_level_gluon_basis.h
 * Contains definition of the class Tree_level_gluon_basis and associated types and operators.
 * Created on: Aug 9, 2012
 * Author: Malin Sjodahl
 */

#include "Tree_level_gluon_basis.h"
#include "parameters.h"
#include <cassert>
#include <fstream>
#include <iostream>


namespace ColorFull {


void Tree_level_gluon_basis::create_basis( int n_gluon ) {

	// Setting basis variable
	ng=n_gluon;

	// Remove a potentially already calculated basis
	cb.clear();
	// The Col_amp containing the basis
	Col_amp Ca_basis;

	// Create the basis using the old function for a maximal number of loops
	Ca_basis = create_trace_basis(ng);

	// Sort the resulting Col_amp into the Col_basis cb
	for ( uint i = 0; i < Ca_basis.ca.size(); i++ ) {
		// A Col_amp to contain a basis vector
		Col_amp Ca_vec;
		Ca_vec.ca.push_back(Ca_basis.ca.at(i));
		cb.push_back(Ca_vec);
	}
}


void Tree_level_gluon_basis::read_in_Col_basis( std::string filename ){


	// First read in basis as normally

	// Read in file
	std::ifstream fin(filename.c_str());

	// Check that file exists
	if( !fin ){
		std::cerr << "Tree_level_gluon_basis::read_in_Col_basis: The file "
				<< filename << " could not be opened." << std::endl;
		assert( 0 );
	}

	// Erase current information
	cb.clear();

	// Copy info from file to string
	std::string str((std::istreambuf_iterator<char>(fin)), std::istreambuf_iterator<char>());

	Col_basis_of_str( str );

	// Then check that it's really a tree level gluon basis
	// Check that length of each Col_amp is 2
	for( uint bv=0; bv < cb.size(); bv++ ){
		if( cb.at(bv).size()!=2 ){
			std::cerr << "Tree_level_gluon_basis::read_in_Col_basis: The basis read in from file " << filename <<
					" has basis vectors with length > 2, and is thus not a tree level gluon basis."
					<< std::endl;
			assert( 0 );
		}

		// Check that the first term is the conjugate of the second
		// Check col_str exactly
		Col_str Cs_conj=cb.at(bv).at(1) ;
		Cs_conj.conjugate();
		Cs_conj.simplify();
		if( cb.at(bv).at(0).cs !=  Cs_conj.cs ){
			std::cerr << "Tree_level_gluon_basis::read_in_Col_basis: The basis read in from file " << filename
					<< " has basis vector " << bv << " with non self-conjugate color structure "
					<< cb.at(bv)
					<<". The col_str " <<  cb.at(bv).at(0).cs << " is not the same as " << Cs_conj.cs << std::endl;
			assert( 0 );
		}
		// Check Polynomial numerically
		if( std::abs(  Col_fun. cnum_num(cb.at(bv).at(0).Poly - (pow(-1.0, ng))* Cs_conj.Poly) ) > accuracy ){
			std::cerr << "Tree_level_gluon_basis::read_in_Col_basis: The basis read in from file " << filename
					<< " has basis vector " << bv << " which is not real due to multiplying Polynomial."
					<< std::endl;
			assert( 0 );
		}

		// Then, remove implicit part
		cb.at(bv).erase(1);
	}
}


void Tree_level_gluon_basis::write_out_Col_basis( ) const{

	write_out_Col_basis( basis_file_name() );

}


void Tree_level_gluon_basis::read_in_Col_basis( ) {

	read_in_Col_basis( basis_file_name() );

}


void Tree_level_gluon_basis::write_out_Col_basis( std::string filename ) const{

	if(  (cb.size()==0 ) ) {
		std::cout << "Tree_level_gluon_basis::write_out_Col_basis(string): There are no basis vectors in this basis, consider using create_basis or read_in_basis." << std::endl;
		std::cout.flush();
		return ;
	}

	std::ofstream outfile(filename.c_str());

	if ( !outfile )
	std::cerr << "Tree_level_gluon_basis::write_out_Col_basis: Cannot write out basis as the file \""
		<< filename.c_str() << "\" could not be opened. (Does the directory exist? Consider creating the directory.)" << std::endl;


	int sign=pow(-1,cb.at(0).n_gluon());
	for (uint m = 0; m < cb.size(); m++) {
		outfile << m << "      "<< cb.at(m);
		if(sign==1 ) outfile << " + ";
		else outfile << " - ";
		Col_amp Ca_conj=cb.at(m);
		Ca_conj.conjugate();
		Ca_conj.normal_order();
		outfile << Ca_conj << std::endl;
	}
	outfile.flush();
}


std::ostream&  Tree_level_gluon_basis::write_out_Col_basis_to_stream( std::ostream&  out ) const{

	if(  (cb.size()==0 ) ) {
		std::cout << "Tree_level_gluon_basis::write_out_Col_basis(): There are no basis vectors in this basis, consider using create_basis." << std::endl;
		std::cout.flush();
	}

	int sign=pow(-1,cb.at(0).n_gluon());
	for (uint m = 0; m < cb.size(); m++) {
		std::cout << m << "      "<< cb.at(m);
		if(sign==1 ) std::cout << " + ";
		else std::cout << " - ";
		Col_amp Ca_conj=cb.at(m);
		Ca_conj.conjugate();
		Ca_conj.normal_order();
		out << Ca_conj << std::endl;
	}
	return out;
}


Polynomial Tree_level_gluon_basis::ij_entry( const int i, const int j ) const{

	// Loop over basis vectors in Basis
	Polynomial ijEntry;
	uint Ng=cb.at(i).n_gluon();

	// The sign of the interference, (-1)^Ng
	int sign=(Ng % 2 ? -1:1);

	Col_amp Cbi_conj=cb.at(i);
	Cbi_conj.conjugate();
	ijEntry=2*Col_fun.scalar_product( cb.at(i), cb.at(j) )
			+sign*2*Col_fun.scalar_product( cb.at(j), Cbi_conj );
	ijEntry.simplify();
	return ijEntry;
}


Col_amp Tree_level_gluon_basis::create_trace_basis( int n_g ) const {

	// To contain the resulting basis
	Col_amp Basis;

	// There has to be at least two gluons
	if (n_g <= 1) {
		std::cerr
				<< "Tree_level_gluon_basis::create_trace_basis: For 0 quarks there is no basis with only "
				<< n_g << " gluons" << std::endl;
		assert( 0 );
	}

	// If 2 or more gluons, build from the 2-gluon basis
	else if (n_g >= 2) {
		Col_str Cs_OnlyState("[(1,2)]");
		Col_amp Ca_tmp;
		Ca_tmp.ca.push_back(Cs_OnlyState);

		Basis = Ca_tmp;
		// For 2 gluons, the work is done
		if (n_g == 2)
			return Basis;
	}

	// Then, add the gluons one at the time
	// As there are only gluons the generation should start from gluon 3,
	// otherwise from 2*n_q+1;
	for (int g_new = 3; g_new <= n_g; g_new++) {
		// If only gluons start from the 2-gluon state, so add gluon 3
		Basis = add_one_gluon(Basis, g_new);
	}

	// Normal order the Col_str's
	Basis.normal_order();
	return Basis;

}


Col_amp Tree_level_gluon_basis::add_one_gluon( const Col_str & Cs, int g_new ) const {

	// For storing the new basis
	Col_amp New_tensors;

	// Add the new gluon in all possible ways to the old Color structure
	// Loop over the Quark_lines
	for (uint ql = 0; ql < Cs.cs.size(); ql++) {
		// The old Quark_line, before insertion of the new gluon index
		Quark_line Old_Ql = Cs.cs.at(ql);
		Col_str New_tensor = Cs;

		// Special case of insertion of a gluon in a 2-ring in a gluons only basis
		// in this case only "half" the basis states for rings with >= 3 gluons are
		// generated. The rest are obtained by taking the indices in anti-cyclic order
		// i.e. by complex conjugating
		if ((Old_Ql.ql.size() == 2 && Cs.gluons_only())) {
			// Insert the new gluon after the existing gluons
			Quark_line New_Ql = Old_Ql;
			New_Ql.ql.push_back(g_new);
			// Replace the old Col_str with the new and add to the new basis states
			New_tensor.cs.at(ql) = New_Ql;
			New_tensors = New_tensors + New_tensor;
		} else { // ordinary case

			// Loop over (potential) insertion places in the Quark_lines, starting from the end
			for (int j = Old_Ql.ql.size(); j > 0; j--) {

				// Special treatment of last place, insert here only if the ring is open
				// (the gluon index cannot take the place of the a quark index)
				if (Old_Ql.open && j == static_cast<int>(Old_Ql.ql.size()))
					j--;

				Quark_line New_Ql = Old_Ql;
				quark_line::iterator it = New_Ql.ql.begin() + j;
				New_Ql.ql.insert(it, g_new);

				// Replace the old Col_str with the new and add to the new basis states
				New_tensor.cs.at(ql) = New_Ql;
				New_tensors = New_tensors + New_tensor;
			}
		}
	}

	return New_tensors;
}


Col_amp Tree_level_gluon_basis::add_one_gluon(const Col_amp & Old_basis, int g_new) const {

	// For storing the new basis
	Col_amp New_bas;

	// Add the new gluon to each of the previous color tensors
	for (uint t = 0; t < Old_basis.ca.size(); t++) {
		New_bas = New_bas + add_one_gluon(Old_basis.ca.at(t), g_new);
	}

	return New_bas;
}

} // end namespace ColorFull

