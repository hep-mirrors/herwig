// -*- C++ -*-
/*
 * Col_amp.cc
 *	Contains the definition of the class Col_amp and related operators
 *  Author: Malin Sjodahl
 */

#include "Col_amp.h"
#include <cassert>
#include <fstream>
#include <iostream>

namespace ColorFull {

std::ostream& operator<<(std::ostream& out, const col_amp & ca) {
	int max = ca.size();
	for (int i = 0; i < max; i++) {
		if(i != 0) out <<" + ";
		out << ca.at(i) ;
	}
	return out;
}


void Col_amp::normal_order_col_strs() {
	for (uint m = 0; m < ca.size(); m++) {
		ca.at(m).normal_order();
	}
}


void Col_amp::normal_order() {

	// First normal order the Col_strs
	normal_order_col_strs();

	// To contain the Col_strs in order
	col_amp ca_ordered;

	// Order the different Col_strs in the Col_amp
	// Do this by moving the Cs one by one to ca_ordered
	while ( ca.size() > 0 ) {
		// Next Cs to put in place (the last Cs in ca)
		Col_str Cs_next = ca.at( ca.size() - 1 );

		// Then insert the Cs among the ordered Col_strs
		// Count how many steps left Cs_next should be moved in ca_ordered
		uint steps_left = 0;
		while ( (steps_left < (ca_ordered.size())) && Cs_next.smallest( Cs_next, ca_ordered.at(ca_ordered.size()-1-steps_left )) ==1) {
			steps_left++;
		}

		// Insert the Cs in the right place among the ordered Col_strs
		col_amp::iterator it=ca_ordered.end()-steps_left;
		ca_ordered.insert( it, Cs_next);

		// Erase the Cs from ca
		ca.erase(  ca.end()-1 ) ;

	}
	ca=ca_ordered;
}


void Col_amp::erase( int i ) {
	ca.erase(ca.begin() + i);
}


void Col_amp::append(col_amp ca_in) {
	for (uint m = 0; m < ca_in.size(); m++) {
		ca.push_back(ca_in.at(m));
	}
}


void Col_amp::read_in_Col_amp( std::string filename) {

	// Read in file
	std::ifstream fin( filename.c_str() );

	// Check that file exists
	if( !fin ){
		std::cerr << "Col_amp::read_in_Col_amp: The file "
				<< filename << " could not be opened." << std::endl;
		assert( 0 );
	}

	// Erase current information
	ca.clear();
	Scalar.clear();

	// Copy info from file to string
	std::string str((std::istreambuf_iterator<char>(fin)), std::istreambuf_iterator<char>());
	str.erase(str.size()-1);
	Col_amp_of_str( str );
}

void Col_amp::write_out_Col_amp( std::string filename ) const {

	if ((ca.size() == 0)) {
		std::cout
		<< "Col_amp::write_out_Col_amp: The Col_amp is empty."
		<< std::endl;
		std::cout.flush();
		return;
	}

	std::ofstream outfile(filename.c_str());


	if ( !outfile )
	std::cerr << "Col_amp::write_out_Col_amp: Cannot write out Col_amp as the file \""
		<< filename.c_str() << "\" could not be opened. (Does the directory exist? Consider creating the directory.)" << std::endl;


	outfile << *this;
}


void Col_amp::Col_amp_of_str( const std::string str ){

	// First split the string into Col_str and Polynomial part
	uint j=0;

	// Check that left and right normal brackets match up
	int left_brackets=0,right_brackets=0;
	while (j < str.size()) {
		if(str.at(j)=='(') left_brackets++;
		if(str.at(j)==')') right_brackets++;
		j++;
	}
	if(left_brackets != right_brackets){
		std::cerr << "Col_amp::Col_amp_of_str: The normal brackets, (), in the Col_amp \"" << str <<"\" do not seem to match up. There were "
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
		std::cerr << "Col_amp::Col_amp_of_str: The curly brackets in the Col_amp \"" << str <<"\" do not seem to match up. There were "
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
		std::cerr << "Col_amp::Col_amp_of_str: The square brackets, [], in the string \"" << str <<"\" do not seem to match up. There were "
				<< left_brackets <<" left bracket(s) and "<< right_brackets << " right bracket(s)." << std::endl;
		assert( 0 );
	}

	if(left_brackets < 1){
		std::cerr << "Col_amp::Col_amp_of_str: The Col_amp string constructor requires at least one Col_str, i.e."
				<< " at least one set of square brackets, []. The string \"" << str << "\" contains "
				<< left_brackets  << std::endl;
		assert( 0 );
	}

	// Find various Col_strs, read until final ]
	j=0;
	std::string Cs_string;
	while (j< str.size() ){

		// Skip some white spaces
		while (j< str.size() and str.at(j)==' ') j++;

		Cs_string.clear();
		while (j< str.size() and str.at(j)!=']'){
			Cs_string.push_back( str.at(j) );
			j++;
		}
		// get fina l ]
		if( j< str.size() and str.at(j)==']' ) Cs_string.push_back( str.at(j) );

		std::cout.flush();
		j++;
		if( !Cs_string.empty() ) ca.push_back(  Col_str(Cs_string) );
	}

	Scalar=Scalar*0;
}


void Col_amp::remove_1_rings() {
	// Loop Col_strs, if empty, remove
	for (uint m = 0; m < ca.size(); m++) {
		// Try removing 1-rings in Col_str
		// If there was a 1-ring, the col_str will multiply 0 after Quark_line.remove_1_rings()
		ca.at(m).remove_1_rings();
		// If the result is multiplying only 0 (in the form of a monomial=0), remove the Col_str
		if (ca.at(m).Poly.size() == 1 && ca.at(m).Poly.at(0).int_part == 0)
			erase(m);
	}
}


void Col_amp::remove_0_rings() {
	// Loop Col_strs, and remove empty Quark_lines
	for (uint m = 0; m < ca.size(); m++) {
		// Try removing 0-rings in Col_str
		ca.at(m).remove_0_rings();
		// If the result is empty (no Quark_lines), keep the Polynomial (in scalar),
		// but remove that Col_str
		if (ca.at(m).cs.empty()) {
			Scalar = Scalar + ca.at(m).Poly;
			erase(m);
		}
	}
}


void Col_amp::remove_empty_Col_strs() {
	// Loop Col_strs
	for (uint m = 0; m < ca.size(); m++) {
		if ( ca.at(m).cs.empty() ){
			Scalar=Scalar+ca.at(m).Poly;
			erase(m);
		}
	}
}


bool Col_amp::gluons_only() const {
	// Loop over Col_strs
	for (uint m = 0; m < ca.size(); m++) {
		if (!ca.at(m).gluons_only())
			return false;
	}
	// If all Cs had gluons only
	return true;
}


int Col_amp::n_gluon() const {
	if	(ca.size()>0) return ca.at(0).n_gluon();
	else{
		std::cerr << "Col_amp::n_gluon(): ca has no Col_str " << std::endl;
		return 0;
	}
}


int Col_amp::n_quark() const {
	if	(ca.size()>0) return ca.at(0).n_quark();
	else{
		std::cerr << "Col_amp::n_quark(): ca has no Col_str " << std::endl;
		return 0;
	}
}


int Col_amp::n_gluon_check() const {
	int ng=n_gluon();
	for(uint n=0; n< ca.size(); n++)
		if( ca.at(n).n_gluon() != ng) {
			std::cerr << "Col_amp::n_gluon_check: The Col_strs in " << ca << " have differently many gluons." << std::endl;

		}
	return ng;
}


int Col_amp::n_quark_check() const {
	int nq=n_quark();
	for(uint n=0; n< ca.size(); n++)
		if( ca.at(n).n_quark() != nq) {
			std::cerr << "Col_amp::n_quark_check: The Col_strs in " << ca << " have differently many quarks." << std::endl;
		}
	return nq;
}


int Col_amp::longest_quark_line() const{

	int length=0;
	for ( uint j = 0; j < ca.size(); j++ ) {
		if( static_cast<int> ( at(j).longest_quark_line() )> length ) length=static_cast<int> ( at(j).longest_quark_line() );
	}

	return length;
}

void Col_amp::collect_col_strs() {

	normal_order_col_strs();

	// Special case of empty Col_amp, do nothing
	if (ca.empty())
		return;

	// Resulting col_amp
	col_amp ca_out;

	// Move first Cs to ca_out
	ca_out.push_back(ca.at(0));
	ca.erase(ca.begin());

	// Move elements (starting from the beginning) in ca to ca_out
	// and sort as long as elements remain to sort
	while (!ca.empty()) {

		// Compare next =0th col_str to all terms in ca_out
		bool was_found = false;
		for (uint i = 0; (i < ca_out.size() && !was_found); i++) {

			if (ca.at(0).cs == ca_out.at(i).cs) {
				was_found = true;
				// Add Polynomial multiplying the Monomial to Polynomial multiplying
				// existing term
				ca_out.at(i).Poly = ca_out.at(i).Poly + ca.at(0).Poly;
				// See if Polynomial can be simplified
				//ca_out.at(i).Poly=simplify( ca_out.at(i).Poly );
			}
		}
		// If the Col_str wasn't represented add info
		if (!was_found)
			ca_out.push_back(ca.at(0));

		// Erase Cs after saving info
		ca.erase(ca.begin());
	}
	ca = ca_out;
}


void Col_amp::simplify() {
	remove_1_rings();
	remove_0_rings();

	// Simplify scalar factor
	Scalar.simplify();
	// Collect similar col_str's
	collect_col_strs();

	// Loop over Col_strs
	for (uint m = 0; (m < ca.size() ); m++) {

		// Simplify Col_strs
		ca.at(m).simplify();

		// If there's only one term in the Polynomial, and the term is 0
		// remove Col_str
		if (ca.at(m).Poly.size() == 1 && ca.at(m).Poly.at(0).int_part == 0){
			erase(m);
		}
	}
}


void Col_amp::conjugate( ) {

	Scalar.conjugate();

	// Take conjugate by reversing order of all partons in each Quark_line
	// Loop over Col_strs
	for (uint i=0; i < ca.size(); i++ ){
		ca.at(i).conjugate();
	}
}


void Col_amp::contract_next_neighboring_gluons(  )  {
	// Loop over Col_str's
	for (uint m = 0; m < ca.size(); m++) {
		ca.at(m).contract_next_neighboring_gluons( );
	}
	return;
}

void Col_amp::contract_2_rings( ){
	// Loop over Col_str's
	for (uint m = 0; m < size(); m++) {
		at(m).contract_2_rings( );
	}
	return;
}


void Col_amp::contract_Quark_line_gluons( Quark_line & Ql )  {

	// Check that all quark indices have been removed
	if (Ql.open){
		std::cerr
		<< "Col_amp::contract_Quark_line_gluons(Ql): all quark indices were not contracted in " << Ql << std::endl;
	}

	if( !ca.empty() or (Scalar.size()==!1 or Scalar.at(0).int_part!=0 ) ){
		std::cerr
		<< "Col_amp::contract_Quark_line_gluons(Ql): This member function "
		<< "stores the result from contracting the Quark_line in the Col_amp itself. "
		<< "It therefore expects an empty initially Col_amp, but it was:" << *this << std::endl;
	}

	// If the Quark_line was empty, return a Col_amp with info in Col_str
	if (Ql.empty()) {
		Col_str Cs;
		// Move Polynomial to Col_str
		Cs.Poly = Ql.Poly;
		Ql.Poly.clear();
		Cs.cs.push_back(Ql);
		ca.push_back(Cs);
		return;
	}

	// To keep track if partner was found or not
	bool found_pair = false;

	// Take first term and look for partner
	for (uint j1 = 0; j1 < Ql.ql.size() - 1; j1++) {
		int the_g = Ql.ql.at(j1);
		// Look for same g
		for (uint j2 = j1 + 1; j2 < Ql.ql.size(); j2++) {
			// if same g found
			if (Ql.ql.at(j2) == the_g) {

				// Special case that the gluons are neighbors
				if (j2 == j1 + 1 or (!Ql.open && j1 == 0 && j2 == Ql.ql.size()
						- 1)) {
					Ql.contract_neighboring_gluons( j1 );
					// Make a Col_str out of the Quark_line
					Col_str Cs;
					Cs.cs.push_back(Ql);
					ca.push_back(Cs);
					simplify();
					return;
				}
				// Special case that the gluons are next to neighbors
				else if (j2 == j1 + 2 //normal case
						or (!Ql.open && j1 == 0 && j2 == Ql.ql.size()- 2) 	//first and second last
				or (!Ql.open && j1 == 1 && j2 == Ql.ql.size()- 1) //second and last
				) {
					if (j2 == j1 + 2)
						Ql.contract_next_neighboring_gluons( j1 );
					if ((!Ql.open && j1 == 0 && j2 == Ql.ql.size() - 2))
						Ql.contract_next_neighboring_gluons( Ql.ql.size() - 2);
					if( (!Ql.open && j1 == 1 && j2 == Ql.ql.size()- 1) )
						Ql.contract_next_neighboring_gluons( Ql.ql.size() - 1);

					// Make a Col_str out of the Quark_line
					Col_str Cs;
					Cs.cs.push_back(Ql);
					ca.push_back(Cs);
					simplify();
					return;
				}
				// normal case, not neighbors or next to neighbors
				else{
					// Make a Col_str out of the Quark_line
					Col_str Cs;

					// Move color factor to multiply Cs instead of Ql
					Cs.Poly = Ql.Poly;
					Ql.Poly.clear();
					Cs.cs.push_back(Ql);

					// Make two copies to store both terms
					ca.push_back(Cs);
					ca.push_back(Cs);

					std::pair <Quark_line, Quark_line> Ql_parts=ca.at(0).cs.at(0).split_Quark_line( j1, j2 );
					Col_str Cs_Ql_parts;
					Cs_Ql_parts.cs.push_back(Ql_parts.first);
					Cs_Ql_parts.cs.push_back(Ql_parts.second);

					ca.at(0)=Cs_Ql_parts;

					// This is multiplying a factor TR
					Monomial Mon_tmp;
					Mon_tmp.pow_TR = 1;

					ca.at(0).Poly = ca.at(0).Poly * Mon_tmp;
					// The split can generate new (next to) neighbors
					ca.at(0).contract_next_neighboring_gluons( );

					// The second Nc suppressed term, obtained by just removing gluon
					// remove g
					quark_line::iterator it1 = ca.at(1).cs.at(0).ql.begin() + j1;
					quark_line::iterator it2 = ca.at(1).cs.at(0).ql.begin() + j2;
					ca.at(1).cs.at(0).ql.erase(it2);
					ca.at(1).cs.at(0).ql.erase(it1);
					// Multiply with -TR/(Nc)
					Mon_tmp.pow_TR = 1;
					Mon_tmp.pow_Nc = -1;
					Mon_tmp.int_part = -1;

					ca.at(1).Poly = ca.at(1).Poly * Mon_tmp;

					// A pair was found, use for stopping
					found_pair = true;
				}
			}
			if (found_pair)
				break; // Stop j2 loop
		}
		if (found_pair)
			break; // Stop j1 loop
	}

	// If a pair was never found return the info of the original Ql
	// but with the Polynomial moved to the Col_str
	if(!found_pair){
		Col_str Cs;
		Cs.Poly=Ql.Poly;
		Ql.Poly.clear();
		Cs.cs.push_back(Ql);
		ca.push_back(Cs);
	}

	// Remove possible 0 and 1-rings
	remove_1_rings();
	remove_0_rings();

	return;
}


void Col_amp::contract_Quark_line_gluons( Col_str & Cs ) {

	if( !ca.empty() or (Scalar.size()==!1 or Scalar.at(0).int_part!=0 ) ){
		std::cerr
		<< "Col_amp::contract_Quark_line_gluons(Cs): This member function "
		<< "stores the result from contracting the Quark_line in the Col_amp itself. "
		<< "It therefore expects an empty initially Col_amp, but it was:" << *this << std::endl;
	}

	Cs.remove_1_rings();
	Cs.remove_0_rings();

	// If the Cs is empty there is no color structure
	if( Cs.empty() ) {
		ca.push_back(Cs);
		return;
	}

	// Make sure to get Poly, both from Cs and Ql.at(0)
	// First the result of the contracted Quakr_line is stored in this Col_amp
	contract_Quark_line_gluons( Cs.at(0) );
	// if the Col_str had a Polynomial we multiply with it here
	*this=*this*Cs.Poly;

	// Loop over Quark_lines
	for( uint i=1; i < Cs.size(); i++ ){

		// If the Quark_line is short do nothing with it as no gluons
		// can not be neighbors or next to neighbors
		if ( Cs.at(i).size()< 6 ){
			*this=*this*Col_str( Cs.at(i) );
		}
		else{
			Col_amp contracted_Ql;
			contracted_Ql.contract_Quark_line_gluons( Cs.at(i) );
			*this*=contracted_Ql;
		}
	}

	simplify();

	return;
}

void Col_amp::contract_Quark_line_gluons( ) {

	// If all the Col_strs are short nothing can be gained
	// compared to contract_next_neighboring gluons, do nothing
	// If there is no remaining color structure, do nothing
	if ( longest_quark_line() < 6 ) return;

	// Make a copy of the old ca part (from this Ca)
	Col_amp Ca_copy;
	Ca_copy.ca=ca;

	// Keep the scalar part, but replace the ca part
	ca.clear();

	for ( uint i = 0; i < Ca_copy.size(); i++ ) {
		Col_amp Ca_from_Cs;
		Ca_from_Cs.contract_Quark_line_gluons( Ca_copy.at(i) );
		*this+=Ca_from_Cs;
	}
}

void Col_amp::contract_quarks( const Col_amp Ca1, const Col_amp Ca2 ) {

	if( !ca.empty() or (Scalar.size()==!1 or Scalar.at(0).int_part!=0 ) ){
		std::cerr
		<< "Col_amp::contract_quarks(Ca1, Ca2): This member function "
		<< "stores the result from contracting quarks in the Col_amp itself. "
		<< "It therefore expects an empty initially Col_amp, but it was:" << *this << std::endl;
	}


	if(Ca1.empty()){
		std::cerr << "Col_amp::contract_quarks: Expects non-empty Col_amps, got first argument "
				<< Ca1 << std::endl;
		assert(0);
	}
	if(Ca2.empty()){
		std::cerr << "Col_amp::contract_quarks: Expects non-empty Col_amps, got second argument "
				<< Ca2 << std::endl;
		assert(0);
	}

	//Col_amp Ca_res;
	Col_amp Ca1_copy=Ca1;
	Col_amp Ca2_copy=Ca2;

	// Make sure the Col_strs are not empty "[]"=1, as all indices contracted
	Ca1_copy.remove_empty_Col_strs();
	Ca1_copy.remove_empty_Col_strs();

	// Loop over Col_strs, and contract quarks between all possible combinations
	// Loop over Col_strs in Ca1
	for(uint m1=0; m1 < Ca1_copy.ca.size(); m1++ ){
		// Loop over Col_strs in Ca2
		for(uint m2=0; m2 < Ca2_copy.ca.size(); m2++ ){
			Col_str Cs_tmp;
			Cs_tmp.contract_quarks( Ca1_copy.ca.at(m1), Ca2_copy.ca.at(m2));
			ca.push_back( Cs_tmp );
		}
	}

	return;
}

void Col_amp::contract_a_gluon( Col_str & Cs )  {

	if( Cs.n_quark()!=0 ){
		std::cerr << "Col_amp::contract_a_gluon(Cs): Expects Col_str with gluons only, got Cs" <<
				std::endl;
	}


	if( !ca.empty() or (Scalar.size()==!1 or Scalar.at(0).int_part!=0 ) ){
		std::cerr
		<< "Col_amp::contract_Quark_line_gluons(Cs): This member function "
		<< "stores the result from contracting the Quark_line in the Col_amp itself. "
		<< "It therefore expects an empty initially Col_amp, but it was:" << *this << std::endl;
	}

	//If the Col_str is empty, or the first ql is empty return it as Col_amp
	if ( Cs.empty() or Cs.cs.at(0).empty() ) {
		ca.push_back(Cs);
		return;
	}

	// Pick first gluon
	int the_g = Cs.at(0, 0);
	std::vector<int> place;
	place.push_back(0);
	place.push_back(0);

	// For storing place of second gluon
	std::vector<int> place2;

	// Locate the same gluon
	for (uint i2 = 0; i2 < Cs.cs.size(); i2++) {
		for (uint j2 = 0; j2 < Cs.cs.at(i2).ql.size(); j2++) {
			// Make sure it is a different gluon
			if ((i2 != 0 or j2 != 0) && Cs.at(i2, j2) == the_g) {
				place2.push_back(i2);
				place2.push_back(j2);
			}
		}
	}

	// If the gluon was not found, something went wrong
	if (place2.empty()) {
		std::cerr << "Col_functions::contract_a_gluon: The gluon " << the_g
				<< " was only found once " << Cs;
	    std::cerr.flush();
		assert( 0 );
	}

	// If the gluon was found in same Ql, use contract_Ql_gluons
	if (place.at(0) == place2.at(0)) {
		contract_Quark_line_gluons( Cs );
		return;
	}
	// If the removed gluon was part of a two-ring
	else if( Cs.cs.at(place.at(0)).ql.size()==2 or Cs.cs.at(place2.at(0)).ql.size()==2){
		Cs.contract_2_rings( );
		ca.push_back(Cs);
		return;
	}
	else
	{
	// The result is a sum of two terms, to contain these
	Col_str Cs1 = Cs;
	Col_str Cs2;

	// The Nc suppressed term, obtained by just removing the g
	// location of first and second g to remove
	quark_line::iterator it1 = Cs1.cs.at(place.at(0)).ql.begin();
	quark_line::iterator it2 = Cs1.cs.at(place2.at(0)).ql.begin() + place2.at(1);
	// Remove the gluon in both places
	Cs1.cs.at(place2.at(0)).ql.erase(it2);
	Cs1.cs.at(place.at(0)).ql.erase(it1);


	// The non-suppressed term, obtained by joining the ql's (where the g is removed)
	Cs2 = Cs1;
	// Construct the contracted ql by
	// 1, taking the first part of the 2nd ql, by erasing everything after the erased gluon
	Quark_line ql_first_part= Cs.cs.at(place2.at(0)).before( place2.at(1) );
	// 2, the in between part=first ql, the contracted gluon is already erased
	Quark_line middle_part=Cs2.cs.at(place.at(0));
	// 3, Last part equal to second part of 2nd ql
	Quark_line last_part= Cs.cs.at(place2.at(0)).after( place2.at(1) );
	// The new combined Ql, to replace the 2nd involved Qualk_line
	Quark_line new_Ql=ql_first_part;
	new_Ql.append( middle_part.ql );
	new_Ql.append( last_part.ql );
	// Copy Polynomial info, both from first and 2nd involved Ql
	new_Ql.Poly=Cs.cs.at(place2.at(0)).Poly*Cs.cs.at(place.at(0)).Poly;
	// Replace the first involved Ql with the constructed Ql
	Cs2.cs.at(0)=new_Ql;

	// erase moved ql
	col_str::iterator it = Cs2.cs.begin() + place2.at(0);
	Cs2.cs.erase(it);

	// Multiply with TR for non-suppressed term
	Monomial Mon_tmp;
	Mon_tmp.pow_TR = 1;
	Cs2.Poly = Cs2.Poly * Mon_tmp;

	// Multiply with -TR/Nc for suppressed term
	Mon_tmp.pow_Nc = -1;
	Mon_tmp.int_part = -1;
	Cs1.Poly = Cs1.Poly * Mon_tmp;

	// Look for simple simplifications in Cs1
	// This has to be done after Cs1 is used to construct Cs2
	// After the removal there may be new next to neighbors in the affected ql's
	// First affected Ql
	// look forward at nn
	Cs1.cs.at(place.at(0)).contract_neighboring_gluons(place.at(0));
	// look backward at nn
	Cs1.cs.at(place.at(0)).contract_neighboring_gluons( place.at(0) - 1);
	// look forward at next nei
	Cs1.cs.at(place.at(0)).contract_next_neighboring_gluons( place.at(0) );
	// look backward at next nei
	Cs1.cs.at(place.at(0)).contract_next_neighboring_gluons( place.at(0) - 2);
	// Second affected ql
	// look forward at nn
	Cs1.cs.at(place2.at(0)).contract_neighboring_gluons( place2.at(0));
	// look backward at nn
	Cs1.cs.at(place2.at(0)).contract_neighboring_gluons( place2.at(0) - 1);
	// look forward at next nei
	Cs1.cs.at(place2.at(0)).contract_next_neighboring_gluons( place2.at(0) );
	// look backward at next nei
	Cs1.cs.at(place2.at(0)).contract_next_neighboring_gluons( place2.at(0) - 2);

	ca.push_back(Cs1);
	ca.push_back(Cs2);

	// remove rings with just one gluon (=0)
	remove_1_rings();
	// remove rings with no partons (=Nc, if closed)
	remove_0_rings();
	}
	return;

}


void Col_amp::contract_a_gluon( ) {

  // If the Ca has only an empty Col_str or a Col_str with an empty Quark_line, do nothing
  if( size()==1 && ( at(0).empty() or at(0).cs.at(0).empty())) return ;


  // Copy ca to copy, not that the ca is in the copy and the Scalar is still in this Ca
  Col_amp Ca_copy;
  Ca_copy.ca=ca;
  ca.clear();

  // Make one contraction in each Col_str
	for ( uint i=0; i< Ca_copy.ca.size(); i++ ) {

			// Contract one gluon in Col_str
			Col_amp Ca_part;
			// The copy has no Scalar, so adding does not double count Scalar
			Ca_part.contract_a_gluon( Ca_copy.ca.at(i) );
			// ... and add resulting Col_amp to Ca_res
			*this+=Ca_part;
	}
	return;
}


void Col_amp::contract_all_gluons( Col_str & Cs )  {

	if( !ca.empty() or (Scalar.size()==!1 or Scalar.at(0).int_part!=0 ) ){
		std::cerr
		<< "Col_amp::contract_all_gluons(Cs): This member function "
		<< "stores the result from contracting the Quark_line in the Col_amp itself"
		<< "It therefore expects an empty initially Col_amp, but it was:" << *this << std::endl;
	}

	ca.push_back( Cs );
	contract_all_gluons();

	return ;
}


void Col_amp::contract_all_gluons( ) {
	//std::cout << "Col_functions::contract_all_gluons: incoming " << Ca << std::endl;


	// Make sure the Col_strs are not empty "[]"=1, as all indices contracted
	remove_empty_Col_strs();
	remove_empty_Col_strs();

	// Check that all quark indices have been removed
	if ( !gluons_only() ) {
		std::cerr << "Col_amp::contract_all_gluons(Ca): Error, all quark indices were not contracted in " << *this << std::endl;
		std::cerr.flush();
		assert( 0 );
	}

	// First normal order and look for simple simplification
	simplify();

	Col_amp Ca_old; // For comparison
	// Contract remaining gluons as long as result keep changing
		// Start with contractions only giving one term
		while (Ca_old != *this) {

			Ca_old = *this;

			// Contract gluons from 2-ring as this gives only one term
			// and may result in new neighbors or next to neighbors
			contract_2_rings( );

			//std::cout << "Col_functions::contract_all_gluons: after 2 rings " << Ca << std::endl;


			// Remove neighboring and next to neighboring gluon indices
			// in each quark_line as this gives only (at most) one term
			contract_next_neighboring_gluons( );
			//std::cout << "Col_functions::contract_all_gluons: after nn " << Ca << std::endl;


			// Contract gluons within same quark_line
			// This will give >1 term, but at least rings will be shorter and shorter
			contract_Quark_line_gluons(  );
			//std::cout << "Col_functions::contract_all_gluons: after Qlg " << Ca << std::endl;


			// To maximize number of CF's
			contract_next_neighboring_gluons( );
			//std::cout << "Col_functions::contract_all_gluons: after nn " << Ca << std::endl;

			// Make "arbitrary" gluon contraction when no simple or Ql contraction remains
			// In each Col_str, contract the first gluon in the first Ql
			contract_a_gluon( );

			// Look if same Col_str is represented more than once in the Col_amp and simplify
			simplify();

		} //end of while (Ca keeps changing)

		return;
}


std::ostream& operator<<(std::ostream& out, const Col_amp & Ca) {
	int max = Ca.ca.size();

	// Write out Scalar part if it is not obviously 0
	Polynomial Poly0; // For comparison, the default Polynomial=1
	Poly0=Poly0*0;
	if( ! (Ca.Scalar==Poly0) ) {
		out << Ca.Scalar << " + ";
	}

	// Print color structure
	if (max == 0)
		out << "{[]}";
	else {
		for (int i = 0; i < max - 1; i++) {
			out << Ca.ca.at(i) << " + ";
		}
		out << Ca.ca.at(max - 1);
	}
	return out;
}


bool operator==( const Col_amp & Ca1, const Col_amp & Ca2 ) {

	// The Ca's should have equal length
	if (Ca1.ca.size() != Ca2.ca.size())
		return false;

	// All Col_strs should be the same
	for (uint i = 0; i < Ca1.ca.size(); i++) {
		if (Ca1.ca.at(i) != Ca2.ca.at(i))
			return false;
	}
	// The scalar part should be the same
	if (Ca1.Scalar!=Ca2.Scalar) return false;
	return true;
}


bool operator!=( const Col_amp & Ca1,const Col_amp & Ca2 ) {
	if (Ca1 == Ca2)
		return false;
	else
		return true;
}


Col_amp operator+( const Col_amp & Ca, const Col_str & Cs ){

	Col_amp Ca_out(Ca);

	// Add Col_str Cs to Ca
	Ca_out.ca.push_back(Cs);

	return Ca_out;
}


Col_amp operator+( const Col_str & Cs, const Col_amp & Ca ){

	return Ca+Cs;
}


Col_amp operator+( const Col_amp & Ca1, const Col_amp & Ca2 ){

	// To contin the result
	Col_amp Ca_out;

	// Add Scalars
	Ca_out.Scalar = Ca1.Scalar +Ca2.Scalar;

	// Add col_strs of Ca1 and Ca2 to Ca_out
	Ca_out.append(Ca1.ca);
	Ca_out.append(Ca2.ca);

	return Ca_out;
}


Col_amp operator+=( Col_amp & Ca1, const Col_amp & Ca2 ){

	// Add Scalars
	Ca1.Scalar+=Ca2.Scalar;

	// Add col_strs of Ca1 and Ca2 to Ca_out
	Ca1.append(Ca2.ca);

	return Ca1;
}

Col_amp operator-( const Col_amp & Ca1, const Col_amp & Ca2 ) {

	Col_amp Ca2_out(Ca2);

	// 	Change sign of Ca2 in Polynomial
	Ca2_out.Scalar=Ca2_out.Scalar*(-1);

	// Loop over Col_strs to change sign for each term
	for( uint m=0; m < Ca2_out.ca.size(); m++ ) {
		Ca2_out.ca.at(m).Poly=Ca2_out.ca.at(m).Poly*(-1);
	}

	// Now the subtraction is equal to an addition
	return Ca1+Ca2_out;
}


Col_amp operator*( const Col_amp & Ca, const int i ){

	Col_amp Ca_res(Ca);
	// Multiply scalar factor (Polynomial)
	Ca_res.Scalar=Ca_res.Scalar*i;

	// Loop over Col_strs and multiply Polynomials
	for(uint m=0; m< Ca_res.ca.size(); m++){
		// Multiply Polynomial of Col_str with the int
		Ca_res.ca.at(m).Poly=Ca_res.ca.at(m).Poly*i;
	}
	return Ca_res;
}
Col_amp operator*( const int i, const Col_amp & Ca ) {

	return Ca*i;
}


Col_amp operator*(const Col_amp & Ca, const cnum c){
	Col_amp Ca_res = Ca;
	// Multiply scalar factor (Polynomial)
	Ca_res.Scalar = Ca_res.Scalar * c;

	// Loop over Col_strs and multiply Polynomials
	for (uint m = 0; m < Ca_res.ca.size(); m++) {
		// Multiply Polynomial of Col_str with the complex number
		Ca_res.ca.at(m).Poly = Ca_res.ca.at(m).Poly * c;
	}
	return Ca_res;
}
Col_amp operator*( const cnum c, const Col_amp & Ca ){

	return Ca*c;
}


Col_amp operator*( const Col_amp & Ca, const double d ){
	Col_amp Ca_res = Ca;
	// Multiply scalar factor (Polynomial)
	Ca_res.Scalar = Ca_res.Scalar * d;
	// Loop over Col_strs and multiply Polynomials
	for (uint m = 0; m < Ca_res.ca.size(); m++) {
		// Multiply Polynomial of Col_str with the complex number
		Ca_res.ca.at(m).Poly = Ca_res.ca.at(m).Poly * d;
	}
	return Ca_res;
}
Col_amp operator*( const double d, const Col_amp & Ca ) {

	return Ca*d;
}


Col_amp operator*( const Col_amp & Ca, const Monomial & Mon ){

	// To contain result
	Col_amp Ca_out(Ca);

	// Multiply scalar factor (Polynomial)
	Ca_out.Scalar=Ca_out.Scalar*Mon;

	// Loop over Col_strs and multiply Polynomials
	for(uint m=0; m< Ca_out.ca.size(); m++){
		// Multiply Polynomial of Col_str with the complex number
		Ca_out.ca.at(m).Poly=Ca_out.ca.at(m).Poly*Mon;
	}
	return Ca_out;
}
Col_amp operator*(const Monomial & Mon, const Col_amp & Ca){
	return Ca*Mon;
}


Col_amp operator*(const Col_amp & Ca, const Polynomial & Poly){

	Col_amp Ca_out(Ca);

	// Multiply scalar factor (Polynomial)
	Ca_out.Scalar=Ca_out.Scalar*Poly;
	// Loop over Col_strs and multiply Polynomials
	for(uint m=0; m< Ca_out.ca.size(); m++){

		// Multiply Polynomial of Col_str with the Polynomial
		Ca_out.ca.at(m).Poly=Ca_out.ca.at(m).Poly*Poly;
	}
	return Ca_out;
}
Col_amp operator*( const Polynomial & Poly, const Col_amp & Ca ){
	return Ca*Poly;
}


Col_amp operator*( const Col_amp & Ca, const Col_str & Cs ){

	Col_amp Ca_res;

	// The Scalar part of Ca gives rise to a Col_str=Cs*Scalar when multiplied with Cs
	// but this should only be kept if the scalar Polynomial is non-zero
	if(! (Ca.Scalar.size()==1 and Ca.Scalar.at(0).int_part==0) ) {
		Col_str Cs_tmp=	Ca.Scalar*Cs;
		Ca_res=Ca_res+Cs_tmp;
	}

	// Multiply each Col_str in Col_amp with Cs
	for(uint m=0; m< Ca.ca.size(); m++){
		// Multiply Col_str in Col_amp with the Col_str
		Col_str Cs_tmp=	Ca.ca.at(m)*Cs;
		Ca_res=Ca_res+Cs_tmp;
	}

	return Ca_res;
}
Col_amp operator*( const Col_str & Cs, const Col_amp & Ca ){

	return Ca*Cs;
}


Col_amp operator*(const Col_amp & Ca1, const Col_amp & Ca2){

	Col_amp Ca_res;

	Ca_res=Ca2*Ca1.Scalar;

	// Multiply each Col_str in Col_amp Ca1 with Ca2 and add to result
	for( uint m=0; m< Ca1.ca.size(); m++){
		// Multiply Col_str in Col_amp with the Col_str
		Ca_res+=Ca1.ca.at(m)*Ca2;
	}

	return Ca_res;
}

Col_amp operator*=( Col_amp & Ca1, const Col_amp & Ca2){

	Col_amp Ca_res;
	Col_amp Ca1_old=Ca1;

	Ca1=Ca2*Ca1_old.Scalar;

	// Multiply each Col_str in Col_amp Ca1 with Ca2 and add to result
	for( uint m=0; m< Ca1_old.ca.size(); m++){
		// Multiply Col_str in Col_amp with the Col_str
		Ca1+=Ca1_old.ca.at(m)*Ca2;
	}

	return Ca1;
}

}// end namespace ColorFull
