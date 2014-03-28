// -*- C++ -*-
/*
 *  Col_str.cc
 *	Contains definition of the class Col_str and associated types and operators.
 *  Created on: Jul 7, 2010
 *  Author: Malin Sjodahl
 */

#include "Col_str.h"
#include <cassert>
#include <fstream>
#include <iostream>


namespace ColorFull {

Col_str::Col_str( const std::string str ) {
	Col_str_of_str( str );
}


void Col_str::Col_str_of_str( const std::string str ) {

	// First split the string into Col_str and Polynomial part
	uint j=0;

	// Check that left and right normal brackets match up
	j=0;
	int left_brackets=0,right_brackets=0;
	while (j < str.size()) {
		if(str.at(j)=='(') left_brackets++;
		if(str.at(j)==')') right_brackets++;
		j++;
	}
	if(left_brackets != right_brackets){
		std::cerr << "Col_str::Col_str_of_str: The normal brackets, (), in the Col_str\"" << str <<"\" do not seem to match up. There were "
				<< left_brackets <<" left bracket(s) and "<< right_brackets << " right bracket(s)." << std::endl;
		assert( 0 );
	}

	// Check that left and right curly brackets match up
	left_brackets=0,right_brackets=0;
	j=0;
	while (j < str.size()) {
		if(str.at(j)=='{') left_brackets++;
		if(str.at(j)=='}') right_brackets++;
		j++;
	}
	if(left_brackets != right_brackets){
		std::cerr << "Col_str::Col_str_of_str: The curly brackets in the Col_str\"" << str <<"\" do not seem to match up. There were "
				<< left_brackets <<" left bracket(s) and "<< right_brackets << " right bracket(s)." << std::endl;
		assert( 0 );
	}

	// Check that left and right [] brackets match up
	j=0;
	left_brackets=0, right_brackets=0;
	while (j < str.size()) {
		if(str.at(j)=='[') left_brackets++;
		if(str.at(j)==']') right_brackets++;
		j++;
	}
	if(left_brackets != right_brackets){
		std::cerr << "Col_str::Col_str_of_str: The square brackets, [], in the string \"" << str <<"\" do not seem to match up. There were "
				<< left_brackets <<" left bracket(s) and "<< right_brackets << " right bracket(s)." << std::endl;
		assert( 0 );
	}


	if( left_brackets != 1){
		std::cerr << "Col_str::Col_str_of_str: Found " << left_brackets << " squared, [], left brackets in the string " << str
				<<" but there should be 1;" << std::endl;
		assert( 0 );
	}

	// Read in the Polynomial until a [ is found
	j=0;
	std::string  Poly_string;
	Poly_string.empty();
	while( j<str.size() and str.at(j)!='['){
		Poly_string.push_back(str.at(j));
		j++;
	}

	// Then read in the col_str
	std::string  col_string;
	col_string.empty();
	while( j< str.size() and str.at(j)!=']'){
		col_string.push_back( str.at(j) );
		j++;
	}

	// Get the final ], In this way the final sign will be ]
	col_string.push_back( str.at(j) );

	Poly=Polynomial( Poly_string );

	col_str_of_str( col_string );
}


void Col_str::col_str_of_str(std::string str) {
	// First sign should be '['
	if (str.at(0) != '[') {
		std::cerr
		<< "Col_str::col_str_of_str: First char in col_str should be '[', it was: "
		<< str.at(0) << std::endl;
		std::cerr.flush();
		assert( 0 );
	}

	uint i = 0; // set i to 0 to start at position 1 after adding 1
	while (i < str.size() - 2) {
		// To contain the argument to the Quark_line constructor
		std::string Ql_arg;
		Ql_arg.clear();

		// Increase i with 1, to compensate for ')'
		i += 1;

		// Make argument to the Quark_line constructor
		// Keep reading while not ')' or '}'
		while (str.at(i) != ')' && str.at(i) != '}') {
			Ql_arg.push_back(str.at(i));
			i++;
			//cout << "Col_str::Col_str(str): For " << i << " The str is: " << Ql_arg << std::endl;
		}
		// Append last char ')' or '}'
		Ql_arg.push_back(str.at(i));
		Quark_line Ql(Ql_arg);
		cs.push_back(Ql_arg);
	}

	// Last sign should be ']'
	if (str.at(str.size() - 1) != ']') {
		std::cerr
		<< "Col_str::col_str_of_str: Last char in col_str should be ']', it was: "
		<< str.at(str.size() - 1) << std::endl;
		std::cerr.flush();
		assert( 0 );
	}
}


void Col_str::read_in_Col_str( std::string filename ) {

	// Read in file
	std::ifstream fin(filename.c_str() );

	// Check that file exists
	if( !fin ){
		std::cerr << "Col_str::read_in_Col_str: The file "
				<< filename << " could not be opened." << std::endl;
		assert( 0 );
	}

	// Copy info from file to string
	std::string str((std::istreambuf_iterator<char>(fin)), std::istreambuf_iterator<char>());
	Col_str_of_str( str );
}


void Col_str::write_out_Col_str( std::string filename ) const {

	if ((cs.size() == 0)) {
		std::cout
		<< "Col_str::write_out_Col_str: The Col_str is empty."
		<< std::endl;
		std::cout.flush();
		return;
	}

	std::ofstream outfile(filename.c_str());


	if ( !outfile )
	std::cerr << "Col_str::write_out_Col_str: Cannot write out Col_str as the file \""
		<< filename.c_str() << "\" could not be opened. (Does the directory exist? Consider creating the directory.)" << std::endl;


	outfile << *this;
}


int Col_str::at( int i, int j ) const {
	if (i < 0) {
		std::cout << "Col_str::at: First argument <0\n";
		std::cerr.flush();
		assert( 0 );
	} else if (i >= static_cast<int>(cs.size()) ) {
		std::cerr << "Col_str::at: First argument > size -1\n";
		std::cerr.flush();
		assert( 0 );
	}
	if (j < 0) {
		std::cerr << "Col_str::at: Second argument <0 \n";
		std::cerr.flush();
		assert( 0 );
	} else if (j >= static_cast<int>(cs.at(i).ql.size())) {
		std::cerr << "Col_str::at: Second argument > size -1\n";
		std::cerr.flush();
		assert( 0 );
	}

	return cs.at(i).ql.at(j);
}


void Col_str::erase( int i ) {
	cs.erase(cs.begin() + i);
}


void Col_str::erase( int i, int j ) {
	cs.at(i).ql.erase(cs.at(i).ql.begin() + j);
}

  
  void Col_str::erase( std::pair<int, int> place ) {
    
    erase(place.first, place.second);
}


void Col_str::insert( int i, int j, int part_num ) {
	// Checking that location is within range
	if (i < 0) {
		std::cerr << "Col_str::insert: First argument <0\n";
		std::cerr.flush();
		assert( 0 );
	} else if (i >= static_cast<int> (cs.size()) ) {
		std::cerr << "Col_str::insert: First argument > size -1\n";
		std::cerr.flush();
		assert( 0 );
	}
	if (j < 0) {
		std::cerr << "Col_str::insert: Second argument <0\n";
		std::cerr.flush();
		assert( 0 );
	} else if (j > static_cast<int>( cs.at(i).ql.size()) ) {
		std::cerr << "Col_str::insert: Second argument > size, was " << j << std::endl;
		std::cerr.flush();
		assert( 0 );
	}
	// Inserting element at place given by iterator itj
	quark_line::iterator itj = cs.at(i).ql.begin() + j;
	cs.at(i).ql.insert(itj, part_num);
}


void Col_str::append( col_str cs_in ) {
	for (uint j = 0; j < cs_in.size(); j++) {
		cs.push_back(cs_in.at(j));
	}
}


std::pair<int, int> Col_str::find_parton( int part_num ) const {

	// Loop over all Quark_lines
	for (uint i = 0; i < cs.size(); i++) {
		// Loop over all places in the Quark_lines
		for (uint j = 0; j < cs.at(i).ql.size(); j++) {
			if (cs.at(i).ql.at(j) == part_num) {
				// cout << j<< "\n";
			  return std::make_pair(i, j);
			}
		}
	}
	// Assert parton found
	std::cerr
	<< "Col_str::find_parton: The function find_parton did not find the parton "
	<< part_num << "in \n" << cs;
	std::cerr.flush();
	assert( 0 );
	return std::make_pair(-1, -1);
}


bool Col_str::neighbor(int p1,int p2) const{

	// The places
  std::pair<int, int> place1 = find_parton(p1);
  std::pair<int, int> place2 = find_parton(p2);

	// First make sure the partons are in the same Quark_line
	if( place1.first!=place2.first ) return false;
	if( place2.second==place1.second + 1 or place2.second==place1.second - 1) return true;

	return false;
}


bool Col_str::right_neighbor(int p1,int p2) const{

  return left_neighbor(p2,p1);

}


bool Col_str::left_neighbor(int p1,int p2) const{

	// The places
  std::pair<int, int> place1 = find_parton(p1);
  std::pair<int, int> place2 = find_parton(p2);

	// First make sure the partons are in the same Quark_line
	if( place1.first!=place2.first ) return false;
	bool closed = !at(place1.first).open;
	int length = at(place1.first).size();
	if ( !(place1.second-place2.second == 1 ||
	       (place1.second-place2.second == 1 - length && closed)) ) return false;

	return true;
}


void Col_str::replace(int old_ind, int new_ind) {

	// First, locate the index to replace
  std::pair<int, int> place = find_parton(old_ind);

	// Then erase the index at that place
	erase(place.first, place.second);

	// Then insert the new index
	insert(place.first, place.second, new_ind);
}


std::string Col_str::find_kind( int part_num ) const {

	// Locate the parton in the Col_str
  std::pair<int, int> place = find_parton(part_num);

	// Check if the quark-line is closed
	// If the parton is in a closed quark line it's a g
	if (!at(place.first).open) {
		return "g";
	}
	// If the parton is first in an open quark-line it's a q
	// (the first relevant place is place 1)
	else if (place.second == 0) {
		return "q";
	}
	// If the parton is last in an open quark-line it's a qbar
	// To find location, subtract the number of characters it takes to write part_num
	else if (place.second == static_cast<int> (cs.at(place.first).ql.size()) - 1) {
		return "qbar";
	}
	// add warning if not found?
	// If the parton wasn't found, it must be a gluon
	else {
		return "g";
	}
}


bool Col_str::gluons_only() const {
	// Loop over Quark_lines
	for (uint i = 0; i < cs.size(); i++) {
		if (cs.at(i).open)
			return false;
	}
	// If no Ql was open, the Cs has gluons only
	return true;
}


int Col_str::n_gluon() const {
	int ng = 0;

	// Loop over Quark_lines
	for (uint i = 0; i < cs.size(); i++) {
		// If the ql is closed, all partons are gluons
		// if it is open, all -2 are gluons
		ng = ng+at(i).ql.size();
		if (at(i).open) ng = ng - 2;
	}
	return ng;
}


int Col_str::n_quark() const {
	int nq = 0;
	// Loop over Quark_lines
	for (uint i = 0; i < cs.size(); i++) {
		// If the ql is closed, there is no gluon
		// If it is open, there is one quark (and one anti-quark)
		if (at(i).open) nq++;
	}
	return nq;
}


void Col_str::normal_order() {

	// All individual quark_lines should be normal_ordered
	for (uint i = 0; i < cs.size(); i++) {
		cs.at(i).normal_order();
	}

	// Loop over the ql's which are candidates to move
	for (uint i = 1; i < cs.size(); i++) {

		// How many steps to the left should the ql be moved?
		int steps_left = 0;
		while (steps_left <= static_cast<int>(i) - 1 && compare_quark_lines(i, i - steps_left - 1) == static_cast<int>(i))
			steps_left++;

		if (steps_left != 0) {
			// Insert ql in new place
			cs.insert(cs.begin() + i - steps_left, cs.at(i));
			// Eras in old
			cs.erase(cs.begin() + i + 1);
		}

	}
}


int Col_str::smallest( const Col_str & Cs1, const Col_str & Cs2 ) const{

	// First order Col_strs according to how many Quark_lines they have
	// The Col_str with fewest Quark_lines is "smallest"
	if ( Cs1.size() < Cs2.size() ) return 1;
	else if( Cs2.size() < Cs1.size() ) return 2;

	// First judge depending on if the Qls are open or not
	// open ql's are "smaller"
	for( uint i=0; i< std::min( Cs1.cs.size(), Cs2.cs.size() ); i++ ){
		// The "smallest" Ql at place i
		if (Cs1.cs.at(i).open &&  !Cs2.cs.at(i).open) return 1;
		else if (Cs2.cs.at(i).open &&  !Cs1.cs.at(i).open) return 2;
	}

	// Then, judge depending on size of the Ql's
	// longer qls are "smaller", should stand first
	for( uint i=0; i< std::min( Cs1.cs.size(), Cs2.cs.size() ); i++ ){
		// The "longest" Ql at place i
		if (     Cs1.cs.at(i).ql.size() > Cs2.cs.at(i).ql.size() ) return 1;
		else if (Cs2.cs.at(i).ql.size() > Cs1.cs.at(i).ql.size() ) return 2;
	}

	// Then, loop over the Quark_lines to see which ql is "smaller", taking index ordering into account
	// First different index decides, ql with smallest index is smaller
	for( uint i=0; i< std::min( Cs1.cs.size(), Cs2.cs.size() ); i++ ){
		// The "smallest" Ql at place i
		int OneOrTwo=Cs1.cs.at(i).smallest( Cs1.cs.at(i), Cs2.cs.at(i) );

		if     ( OneOrTwo==1 ) return 1;
		else if( OneOrTwo==2 ) return 2;
	}

	//If the col_str's are identical, return 0
	return 0;
}


int Col_str::longest_quark_line() const{

	int length=0;
	for ( uint j = 0; j < cs.size(); j++ ) {
		if( static_cast<int> ( at(j).size() )> length ) length=static_cast<int> ( at(j).size() );
	}

	return length;
}


void Col_str::remove_1_rings() {
	// Loop over quarl_lines
	for (uint j = 0; j < cs.size(); j++) {
		// If a quark_line contains only one gluon, replace whole Col_str with a 0-monomial
		if (cs.at(j).ql.size() == 1) {
			// If the quark_line is closed the Col_str is 0
			if (!cs.at(j).open) {
				// Remove irrelevant Polynomial and col_str
				// The col_str is now multiplying 0
				cs.clear();
				Poly.clear();
				Monomial Mon0;
				Mon0.int_part = 0;
				Poly.push_back(Mon0);
			}

			// If the Ql is open and has only one element, something is wrong
			else if (cs.at(j).open) {
				std::cerr
				<< "Col_str::remove_1_rings: An open quark_line cannot have only one parton, but it had in \n"
				<< cs << std::endl;
				std::cerr.flush();
				assert( 0 );
			}
		}
	}
}


void Col_str::remove_0_rings() {
	// Loop over Quark_lines
	for (int j = 0; j < static_cast<int>(cs.size()); j++) {
		// If a quark_line contains no gluons and is closed
		// it is equal to Nc, move factor Nc to the Polynomial of the Col_str
		if (cs.at(j).ql.size() == 0) {
			// Move color factor of the Quark_line to the Polynomial of the Col_str
			Poly = Poly * cs.at(j).Poly;
			// Multiply with Nc if the quark_line is closed
			if (!cs.at(j).open) {
				Monomial Mon_tmp;
				Mon_tmp.pow_Nc = 1;
				Poly = Poly * Mon_tmp;
			}
			// If the ql is open and has 0 elements,
			// it is defined as 1 and can be removed
			// Erase the Ql
			erase(j);

			// In order to check all elements, decrease j when Ql removed
			// j may get to -1, so can not be seen as uint
			j--;
		}
	}
}


void Col_str::simplify() {
	remove_1_rings();
	remove_0_rings();

	// Move factors multiplying the individual Ql's to multiply
	// the Col_str instead
	for (uint i = 0; i < cs.size(); i++) {
		Poly = Poly * cs.at(i).Poly;
		cs.at(i).Poly.clear();
	}

	// Simplify Polynomial of Col_str
	Poly.simplify();

	// Normal order
	normal_order();
}


void Col_str::conjugate(){

	// Conjugating Polynomial
	Poly.conjugate();

	// Take conjugate of cs by conjugating each Quark_line
	for (uint i=0; i < cs.size(); i++ ){
		cs.at(i).conjugate();
	}
}


int Col_str::compare_quark_lines( int i1, int i2 ) const {

	int OneOrTwo= cs.at(i1).smallest( cs.at(i1), cs.at(i2) );
	if ( OneOrTwo==1 ) return i1;
	else if ( OneOrTwo==2 ) return i2;
	// If ql's are equal, return i1
	else if ( OneOrTwo==0 ) return i1;
	else{
		std::cerr << "Col_str::compare_quark_lines: cannot decide on ordering of quark_lines "
				<< cs.at(i1) << " and " << cs.at(i2);
		return 0;
	}

}

void Col_str::contract_2_rings( ) {

	// Contract gluons from 2-ring, this gives only one term
	// For storing place of two-ring
	std::vector<int> place1;
	// For storing place of gluon with same index as second gluon in two-ring
	// i.e. the gluon to be replaced with the first index in the 2-ring
	std::vector<int> place2;
	// The second index in two-ring, to be removed
	int the_g;

	// Loop over Quark_lines to find a two-ring
	for (uint i = 0; i < cs.size(); i++) {
		place1.clear();
		place2.clear();
		// Search for 2-rings
		// A gluon 2-ring has length 2 and is closed
		if ( cs.at(i).ql.size() == 2 && !cs.at(i).open ) {

			// Save place of two-ring
			place1.push_back(i);
			// Pick second gluon (for no good reason)
			// this index should be replaced by first index
			// in the other place where it occurs, and the 2-ring should be removed
			place1.push_back(1);

			the_g = at(place1.at(0), place1.at(1));

			// Now, in all Quark_lines, look for same the_g until found
			for ( uint i2 = 0; i2 < size(); i2++ ) {
				for (uint j2 = 0; j2 < at(i2).ql.size(); j2++) {
					// If the_g was found at a DIFFERENT place, (not same g again)
					if ((i2 != i or j2 != 1) && at(i2, j2) == the_g) {
						place2.push_back(i2);
						place2.push_back(j2);
					}
				}
			}

			// Check that the gluon index was found again
			if (place2.empty()) {
				std::cerr << "Col_str:contract_2_rings: Only found the index " << the_g
						<< " once in " << *this << std::endl;
				assert( 0 );
			}

			// If the_g was found twice in the same ql
			if (place1.at(0) == place2.at(0)) {
				// Keep the Monomial
				// Multiply with Nc Mon from the contraction
				Monomial Mon_tmp;
				Mon_tmp.pow_Nc = 1;
				Mon_tmp.pow_CF = 1;
				// Multiply the Poly of the Col_str with the Poly of the ql
				// and the color factor from the contraction
				Poly = Poly * cs.at(place1.at(0)).Poly * Mon_tmp;
				// Erase the ql
				erase(place1.at(0));
			}
			else{
				// That index should be changed to the index of the first gluon
				// in the 2-ring
				cs.at(place2.at(0)).ql.at(place2.at(1))
							= at(place1.at(0), 0);
				// Multiply the Poly of the Col_str with the Poly of the ql
				// and the color factor from the contraction
				Monomial Mon_tmp;
				Mon_tmp.pow_TR = 1;
				Poly = Poly * Mon_tmp*cs.at(place1.at(0)).Poly;
				// The two ring should be removed
				cs.erase( cs.begin() + i);
			}
			i--; // if we found a 2-ring, we also erased it
		} // end if we found a 2-ring
	}// end looping over Quark_lines
	// One-rings may have been created
	remove_1_rings();

	return;
}

void Col_str::contract_quarks( const Col_str Cs1, const Col_str Cs2 ) {

	if( !cs.empty() or Poly.size()!=0 ){
		std::cerr
		<< "Col_str::contract_quarks(Cs1,Cs2): This member function "
		<< "stores the result from contracting quarks in the Col_str itself. "
		<< "It therefore expects an empty initially Col_str, but it was:" << *this << std::endl;
	}

	std::vector<int> q_place;
	std::vector<int> q_place2;

	// The conjugate of Cs1
	Col_str conj_Cs1 = Cs1;
	conj_Cs1.conjugate();

	// The total color structure
	*this = conj_Cs1*Cs2;

	// Count how many quarks should be contracted
	int n_q = n_quark();

	// As long as there are quark_lines left to contract
	while (n_q > 0) {
		// Find first quark in Cs1 by looping over Quark_lines
		for (int i = 0; (n_q>0 && i <  static_cast<int>( size()) ); i++) {
			// Check if the quark-line is open, in which case it has a q
			if ( cs.at(i).open ) {
				// The first quark is located and has position
				q_place.clear();
				q_place.push_back(i);
				q_place.push_back(0);
				// and number
				int q = at(q_place.at(0), q_place.at(1));

				// Locate same quark a second time
				// Loop over Quark_lines
				q_place2.clear();
				int i2 = i + 1; // Quark_line of second occurrence
				while (q_place2.empty()) { // As long as quark not found a second time
					if ( cs.at(i2).at( cs.at(i2).ql.size() - 1) == q) {// If quark found, store place
						q_place2.push_back(i2);
						q_place2.push_back( cs.at(i2).ql.size() - 1);
					}
					i2++;
				}
				if (q_place2.empty()) {
					std::cerr << "Col_functions::contract_quarks(Cs1, Cs2): Found q " << q
							<< " only once in " << *this << std::endl;
				}

				// Prepare new Quark_line
				// to be inserted at the place of found open Quark_line
				Quark_line new_Quark_line;
				Quark_line part2_new_Quark_line;
				// The first part of the new Quark_line should be the Quark_line
				// containing q in the conjugate
				new_Quark_line = cs.at(q_place2.at(0));

				// Erasing q in the end
				new_Quark_line.ql.erase(--new_Quark_line.ql.end());
				part2_new_Quark_line = cs.at(q_place.at(0));

				// Erasing the q in the beginning of second part
				part2_new_Quark_line.ql.erase(part2_new_Quark_line.ql.begin());

				new_Quark_line.append(part2_new_Quark_line.ql);

				// So far we have not included the Polynomial of part2_new_Quark_line
				new_Quark_line.Poly=new_Quark_line.Poly*part2_new_Quark_line.Poly;

				// If the first q index and the last qbar index in the new
				// Quark_line is the same (and the Quark_line is "open"), the indices
				// should be removed and the Quark_line should be closed
				if (new_Quark_line.ql.at(0) == new_Quark_line.ql.at(
						new_Quark_line.ql.size() - 1) && new_Quark_line.open) {
					// The string is closed
					new_Quark_line.open = false;
					// Remove last and first index
					new_Quark_line.ql.erase(--new_Quark_line.ql.end(),
							new_Quark_line.ql.end());
					new_Quark_line.ql.erase(new_Quark_line.ql.begin());
				}

				// Inserting new Quark_line in the place of the old
				cs.at(i) = new_Quark_line;

				// Remove quark_line with q in Cs
				cs.erase(( cs.begin() + q_place2.at(0)));
				i=-1; // reset to keep looking from the beginning in the new Cs (i will be increased to 0)
			}// end of if (open)

			n_q = n_quark();

		} // end of for, loop over quark_lines
	} // end while (n_q > 0)
	return;
}



void  Col_str::contract_next_neighboring_gluons( ) {

  // Loop over Quark_lines and remove neighboring and next to neighboring
  // gluon indices
  for( uint i=0; i < cs.size(); i++ ){

    // Create a Quark_line to use as argument
    Quark_line Ql= at(i);

    Ql.contract_next_neighboring_gluons( );

    // Collect Polynomial in Cs Polynomial
    Poly=Poly*Ql.Poly;
    // and put powers in Ql to 0
    Ql.Poly.clear();

    cs.insert( cs.begin() +i, Ql);

    // Erase old version of Quark_line (now at place i+1)
    cs.erase( cs.begin() +i + 1 );
  }
  // Remove 1 and 0-rings, simplify Poly and normal order
  simplify();

  return;
}

bool operator==(const col_str & cs1, const col_str & cs2){

	// col_str's must have equal length
	if( cs1.size() != cs2.size() ) return false;

	// Individual ql's be be equal
	for ( uint i=0; i< cs1.size(); i++){
		if(cs1.at(i).ql != cs2.at(i).ql ) return false; // Equal color structure
		if(cs1.at(i).Poly != cs2.at(i).Poly ) return false; // Equal polynomial structure
	}
	//If all ql's equal the cs are considered equal
	return true;
}


bool operator!=(const col_str & cs1, const col_str & cs2){
	if(cs1==cs2) return false;
	else return true;
}


std::ostream& operator<<(std::ostream& out, const col_str & cs) {
	int max = cs.size();
	if (max == 0)
		out << "[]";
	else {
		out << "[";
		for (int i = 0; i < max - 1; i++) {
			out << cs.at(i);
		}
		out << cs.at(max - 1) << "]"; //<< std::endl;
	}
	return out;
}


Col_str operator*( const Col_str & Cs, const int i){
	Col_str Cs_res=Cs;
	Cs_res.Poly=Cs.Poly*i;
	return Cs_res;
}
// Define the operator * for int and Col_str
Col_str operator*( const int i, const Col_str & Cs){
	Col_str Cs_res=Cs;
	Cs_res.Poly=Cs.Poly*i;
	return Cs_res;
}


Col_str operator*( const Col_str & Cs, const double d){
	// Multiply Polynomial of Col_str with the double
	Col_str Cs_res=Cs;
	Cs_res.Poly=Cs.Poly*d;
	return Cs_res;
}
Col_str operator*( const double d, const Col_str & Cs){
	Col_str Cs_res=Cs;
	Cs_res.Poly=Cs.Poly*d;
	return Cs_res;
}



Col_str operator*(const Col_str & Cs, const cnum c){
	Col_str Cs_res=Cs;
	Cs_res.Poly=Cs.Poly*c;
	return Cs_res;
	return Cs;
}
Col_str operator*(const cnum c, const Col_str & Cs){
	Col_str Cs_res=Cs;
	Cs_res.Poly=Cs.Poly*c;
	return Cs_res;
}


Col_str operator*(const Col_str & Cs, const Monomial & Mon){
	Col_str Cs_res=Cs;
	Cs_res.Poly=Cs.Poly*Mon;
	return Cs_res;
}
Col_str operator*( Monomial & Mon, const Col_str & Cs){
	Col_str Cs_res=Cs;
	Cs_res.Poly=Cs.Poly*Mon;
	return Cs_res;
}


Col_str operator*( const Col_str & Cs, const Polynomial & Poly ){
	Col_str Cs_res=Cs;
	Cs_res.Poly=Cs.Poly*Poly;
	return Cs_res;
}
Col_str operator*( const Polynomial & Poly, const Col_str & Cs ){
	Col_str Cs_res=Cs;
	Cs_res.Poly=Cs.Poly*Poly;
	return Cs_res;
}


Col_str operator*( const Col_str & Cs, const Quark_line & Ql){

	// Col_str to return
	Col_str out_Cs=(Cs);

	// Quark_line copy
	Quark_line out_Ql(Ql);

	//  Multiplying Polynomials
	out_Cs.Poly = out_Ql.Poly*Cs.Poly;

	// Remove Polynomial info from Quark_line (put to 1),
	// as this info is stored in Polynomial of Col_str instead
	out_Ql.Poly.clear();

	// Append Quark_line (without multiplicative polynomial) to out_Col_str
	out_Cs.cs.push_back( out_Ql );

	return out_Cs;
}
Col_str operator*( const Quark_line & Ql , const Col_str & Cs){
	return Cs*Ql;
}


Col_str operator*( const Col_str & Cs1, const Col_str & Cs2 ){

	// Col_str to return
	Col_str out_Cs=Cs1;

	//  Multiplying Polynomials
	out_Cs.Poly = Cs1.Poly*Cs2.Poly;

	// Looping over Quark_lines in Cs2 to add to Cs1
	for (uint i=0; i < Cs2.cs.size(); i++ ){
		// Append i:th Quark_line at i:th place in out_Col_str
		out_Cs.cs.push_back( Cs2.cs.at(i) );
	}

	return out_Cs;
}


std::ostream& operator<<(std::ostream& out, const Col_str & Cs){

	// Write out polynomial if it is not 1, or if the cs has no structure
	Polynomial Poly1; // For comparison, the default Polynomial=1

	if( Cs.Poly!=Poly1 or Cs.cs.size()==0 ) out << Cs.Poly;
	out << Cs.cs;
	return out;
}


bool operator==(const Col_str & Cs1, const Col_str & Cs2){

	// col_str's be equal
	if( Cs1.cs != Cs2.cs ) return false;

	// Poly's must be equal
	if( Cs1.Poly != Cs2.Poly ) return false;

	//If all col_strs and all Polynomials are equal the Col_strs are considered equal
	return true;
}


bool operator!=( const Col_str & Cs1, const Col_str & Cs2){
	if( Cs1==Cs2 ) return false;
	else return true;
}

} //end namespace ColorFull
