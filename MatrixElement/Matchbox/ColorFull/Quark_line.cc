// -*- C++ -*-
/*
 * Quark_line.cc
 *	Contains definitions of the class Quark_line and associated types and operators.
 *  Created on: Jul 7, 2010
 *  Author: Malin Sjodahl
 */

#include "Quark_line.h"
#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>


namespace ColorFull {

Quark_line::Quark_line( const std::string str ) {
	Quark_line_of_str( str );
}


Quark_line::Quark_line() {
	// member set just to avoid complaints
	open=true;
}


int Quark_line::at( int j ) const{
	int size = ql.size();
	if (0 <= j and j < size)
		return ql.at(j);
	else if ( !open and (size <= j && j < 2 * size) )
		return ql.at(j - ql.size());
	else if ( !open and (-size < j and j < 0) )
		return ql.at(j + size);

	std::cerr << "Quark_line::at(j): j=" << j << " is out of range" << std::endl;
	std::cerr.flush();
	assert(0);
	return 0;

}


void Quark_line::Quark_line_of_str( const std::string str ){

	int j=0;

	// Check that the string contains a comma
	bool contains_comma=false;
	while (j < (int)str.size()) {
		if(str.at(j)==',') contains_comma=true;
		j++;
	}
	if( !contains_comma ){
		std::cerr << "Quark_line::Quark_line_of_str: The string " << str
				<< " contains no comma. Each quark_line must contain a comma since (g1)=0"
				<< " and for each quark, there is an anti-quark.";
		assert(0);

	}

	// Check that left and right brackets match up
	j=0;
	int left_brackets=0,right_brackets=0;
	while (j < (int)str.size()) {
		if(str.at(j)=='(') left_brackets++;
		if(str.at(j)==')') right_brackets++;
		j++;
	}
	if(left_brackets != right_brackets){
		std::cerr << "Quark_line::Quark_line_of_str: The brackets in the Quark_line\"" << str <<"\" do not seem to match up. There were "
				<< left_brackets <<" left bracket(s) and "<< right_brackets << " right bracket(s)." << std::endl;
		assert( 0 );
	}

	// Check that left and right curly brackets match up
	left_brackets=0,right_brackets=0;
	j=0;
	while (j < (int)str.size()) {
		if(str.at(j)=='{') left_brackets++;
		if(str.at(j)=='}') right_brackets++;
		j++;
	}
	if(left_brackets != right_brackets){
		std::cerr << "Quark_line::Quark_line_of_str: The curly brackets in the Quark_line\"" << str <<"\" do not seem to match up. There were "
				<< left_brackets <<" left bracket(s) and "<< right_brackets << " right bracket(s)." << std::endl;
		assert( 0 );
	}

	if( left_brackets > 1){
		std::cerr << "Quark_line::Quark_line_of_str: Found " << left_brackets << " curly left brackets in the string " << str
				<<" but there should be at most 1;" << std::endl;
		assert( 0 );
	}

	// The string should be split into a part containing the quark_line
	int i=0;
	std::string Poly_str, ql_str;
	Poly_str.clear();
	ql_str.clear();

	// If there were curly brackets just read until one is found
	if( left_brackets>0 ) {
		while( i< (int)str.size() and str.at(i)!='{'){
			Poly_str.push_back(str.at(i));

			i++;
		}
	}
	// If no curly brackets, look for ,
	else {
		// See if the Col_str contains ",",
		// then we know that we should look backwards
		// until the ( is found
		int first_comma=0;
		j=0;
		while ( j < (int)str.size() and first_comma==0 ) {
			if(str.at(j)==',') first_comma=j;
			j++;
		}

		if (first_comma>0){
			j=first_comma;
			while (j>=0 and str.at(j)!= '(') {
				j--;
			}
		}
		// We should have found a (
		if (j==0 and str.at(j)!= '('){
			std::cerr << "Quark_line::Quark_line_of_str: Found no starting bracket in string" << str << std::endl;
			assert( 0 );
		}

		// Now, the quark_line starting ( is at place j
		i=0;
		while(i<j){
			Poly_str.push_back(str.at(i));
			i++;
		}
	}

	// The rest of the string contains the quark_line
	while( i < (int)str.size() ){
		ql_str.push_back(str.at(i));
		i++;
	}

	// Assigning info
	Poly=Polynomial(Poly_str);
	quark_line_of_str(ql_str);

}


void Quark_line::read_in_Quark_line( std::string filename ) {

	// Read in file
	std::ifstream fin( filename.c_str() );

	// Check that file exists
	if( !fin ){
		std::cerr << "Quark_line::read_in_Quark_line: The file "
				<< filename << " could not be opened." << std::endl;
		assert( 0 );

	}

	// Erase current information
	ql.clear();
	Poly.clear();

	// Copy info from file to string
	std::string str((std::istreambuf_iterator<char>(fin)), std::istreambuf_iterator<char>());
	Quark_line_of_str( str );
}


//uint Quark_line::size() const{ return ql.size();}

void Quark_line::quark_line_of_str(std::string str) {
	//std::cout << "quark_line_of_str, got \"" << str.c_str() << "\""<< endl;

	// First char tells if the string is open or not
	if (str.at(0) == '{')
		open = true;
	else if (str.at(0) == '(')
		open = false;
	else {

		std::cerr << "Quark_line::quark_line_of_str: First and last char in Quark_line should be '{}' for open or '()' as in closed. It was: "
				<< str.at(0) << std::endl;
		std::cerr.flush();
		assert(0);
	}

	// Next parton to be added (as string)
	std::string parton_str;
	// Next parton to be added (as int)
	int parton_num;

	// To keep track of the place
	uint i = 1;
	// Loop over characters in string, last char is '}' or ')'
	while (i < str.size() - 1) {
		// When a ',' is found, write content in parton_str to the
		// quark_line and
		if (str.at(i) == ',') {
			std::istringstream parton_str_st(parton_str);
			parton_str_st >> parton_num;
			ql.push_back(parton_num);
			i++;
			// Empty parton string to prepare for new content
			parton_str.clear();
		}
		while (i < str.size() - 1 && str.at(i) != ',') {
			parton_str.push_back(str.at(i));
			i++;
		}
	}
	// Adding last term
	std::istringstream parton_str_st(parton_str);
	parton_str_st >> parton_num;
	ql.push_back(parton_num);

	if (str.at(str.size() - 1) != ')' && str.at(str.size() - 1) != '}') {
		std::cerr
		<< "Quark_line::quark_line_of_str:  Last char in Quark_line should match '{' for open or '(' as for closed. It was: "
		<< str.at(str.size() - 1) << std::endl;
		std::cerr.flush();
		assert(0);
	}
}


void Quark_line::write_out_Quark_line( std::string filename ) const {

	if ((ql.size() == 0)) {
		std::cout
		<< "Quark_line::write_out_Quark_line: The Quark_line is empty."
		<< std::endl;
		std::cout.flush();
		return;
	}

	std::ofstream outfile(filename.c_str());

	if ( !outfile )
	std::cerr << "Quark_line::write_out_Quark_line: Cannot write out Quark_line as the file \""
		<< filename.c_str() << "\" could not be opened. (Does the directory exist? Consider creating the directory.)" << std::endl;

	outfile << *this;
}


void Quark_line::normal_order() {
	// Do nothing if the quark_line is empty
	if (ql.empty())
		return;

	if (!open) {
		// For storing smallest gluon number (should be first in normal ordering)
		int smallest = at(0);
		int place = 0;

		for (uint j = 1; j < ql.size(); j++) {
			if (at(j) < smallest) {
				smallest = at(j);
				place = j;
			}

			// If indices are equal, check next
			else if (at(j) == smallest) {
				// Keep checking next until a difference is found
				// or all elements are compared
				int k = 1;
				while (at(j + k) == at(place + k) && k < static_cast<int>(ql.size()))
					k++;
				// Now, at place j+k, the indices are different
				// Change place if smaler
				if (at(j + k) < at(place + k)) {
					place = j;
				}
			}
		}
		// Not the place of the lowest index is located and the quark_lins should
		// start with that index
		// Move elements before place of lowest number to the end
		for (int j = 0; j < place; j++) {
			ql.push_back(at(0));
			ql.erase(ql.begin());
		}
	}
}


void Quark_line::conjugate( ) {

	Poly.conjugate();

	// To contain the conjugated quark_line
	quark_line ql_new;

	// Take conjugate by reversing order of all partons in the Quark_line
	for (int i=ql.size()-1; i>=0;){
		ql_new.push_back( ql.at(i) );
		i--;
	}

	ql=ql_new;
}

// Returns a Quark_line with the the elements before j in a Quark_line
Quark_line Quark_line::before(  int j ) const{
	//cout << "before( Quark_line Ql, int j ): the size was " << static_cast<int>(Ql.ql.size());
	//cout << "for the Quark_line " << Ql << endl;

	Quark_line Ql_copy=*this;

	if( j<0 or j> static_cast<int>( ql.size()-1) ){
		std::cerr << "Quark_line::before( int j ): the size was " << static_cast<int>( ql.size() );
		std::cerr << "for the quark_line " << *this << std::endl;
		std::cerr.flush();
		std::cerr <<"before: the argument must be >=0 and inside vector, was " << j << std::endl;
		std::cerr.flush();
		assert(0);
	}

	// In first part erase everything after and including j
	quark_line::iterator it1=Ql_copy.ql.begin()+j;
	quark_line::iterator it2=Ql_copy.ql.end();


	Ql_copy.ql.erase(it1, it2);

	return Ql_copy;
}


Quark_line Quark_line::after( int j ) const{

	Quark_line Ql_copy=*this;

	if(j<0 or j> static_cast<int>( ql.size()-1)){
		std::cerr <<"Quark_line::after: the argument must be >=0 and inside vector, was " << j << std::endl;
		std::cerr.flush();
		assert(0);
	}

	// In first part erase everything after and including j
	quark_line::iterator it1=Ql_copy.ql.begin();
	quark_line::iterator it2=Ql_copy.ql.begin()+j+1;
	Ql_copy.ql.erase(it1, it2);

	return Ql_copy;
}


void Quark_line::erase( int i ){
	ql.erase( ql.begin()+i );
}


void Quark_line::append( int p ){
	ql.insert( ql.end(), p );
}


void Quark_line::append(  const std::vector<int> & in_ql ){
	for (uint j=0; j<in_ql.size(); j++ ){
		ql.push_back( in_ql.at(j) );
	}
}


void Quark_line::prepend( int p){
	ql.insert( ql.begin(), p );
}


void Quark_line::prepend( std::vector<int> in_ql ){
	for (uint j=0; j<ql.size(); j++ ){
		in_ql.push_back(ql.at(j));
	}
	ql=in_ql;
}


void Quark_line::insert( int j, int p ){
	if( j< static_cast<int> ( ql.size() ) and j>=0 )
	{
		quark_line::iterator itj = ql.begin() + j;
		ql.insert(itj, p);
	}
	else if( j>= static_cast<int> ( ql.size() )) {
		std::cerr << "Quark_line::insert: The size of ql is " << ql.size()
									<< ", so no parton can be inserted at place " << j
									<< " in " << *this << std::endl;
		assert(0);
	}
	else if( j<0 ) {
		std::cerr << "Quark_line::insert: Can not insert at place " << j
				<< " < 0 in " << *this << std::endl;
		assert(0);
	}
}


std::pair< Quark_line, Quark_line > Quark_line::split_Quark_line( int j1, int j2 ) const{

	if(open){
		std::cerr <<"Quark_line::split_Quark_line: expects a closed quark_line" << std::endl;
		std::cerr.flush();
		assert(0);
	}
	Quark_line Ql1=*this;
	Quark_line Ql2=*this;

	// In first part erase everything between j1 and j2
	quark_line::iterator it1=Ql1.ql.begin()+j1;
	quark_line::iterator it2=Ql1.ql.begin()+j2;
	Ql1.ql.erase(it1, it2+1);

	// In second part, erase from j2...
	it1=Ql2.ql.begin()+j1;
	it2=Ql2.ql.begin()+j2;
	Ql2.ql.erase( it2, Ql2.ql.end() );
	// ... and before j1
	Ql2.ql.erase( Ql2.ql.begin(), it1+1 );

	return std::make_pair(Ql1,Ql2);

}

int Quark_line::smallest(const Quark_line & Ql1, const Quark_line & Ql2) const {

	// If only one is open, that Ql should stand first
	if (Ql1.open && !Ql2.open)
		return 1;
	else if (Ql2.open && !Ql1.open)
		return 2;

	// If both are open or both are closed, the longest should stand first
	else if (Ql1.ql.size() > Ql2.ql.size())
		return 1;
	else if (Ql2.ql.size() > Ql1.ql.size())
		return 2;

	// If the size is the same, the Ql with smallest starting number should stand first
	// If the first number is the same, check the 2nd number, then the 3rd...
	else {
		// Loop over places in the Ql's
		for (uint j = 0; j < Ql1.ql.size(); j++) {
			if (Ql1.ql.at(j) < Ql2.ql.at(j))
				return 1;
			if (Ql2.ql.at(j) < Ql1.ql.at(j))
				return 2;
		}
		// If the Ql's are equal, return 0
		return 0;
	}
}


void Quark_line::contract_neighboring_gluons(int j) {

	if (ql.empty())
		return;

	if (open) {
		std::cerr
				<< "Quark_line::contract_neighboring_gluons(j): Expects a closed Quark_line, got "
				<< *this << std::endl;
		assert(0);
	}

	// If asked for -1st element, check last
	if (j == -1)
		j = size() - 1;
	// If asked for last+1 element, check first
	if (j == static_cast<int>(ql.size() - 1))
		j = 0;

	// Keep contracting gluons as long as:
	// there are at least two gluons to contract
	// neighbors are equal
	// and j is not the last parton in the ql
	while (j < static_cast<int>(size() - 1) && (size()) >= 2
			&& (at(j)) == at(j + 1)) {

		quark_line::iterator it1 = ql.begin() + j;
		quark_line::iterator it2 = ql.begin() + j + 2;

		// Removing neighboring gluons
		ql.erase(it1, it2);

		// This should increase the power of CF
		Monomial Mon_tmp;
		Mon_tmp.pow_CF = 1;
		Poly = Poly * Mon_tmp;
		if (j > 2)
			j = j - 2;
	}

	// If j is the last parton and the Quark_line is closed and first=last
	// and there are at least 2 gluons
	while (static_cast<int>(size()) >= 2 && j == static_cast<int>(size() - 1)
			&& !open && at(0) == at(size() - 1)) {

		quark_line::iterator it1 = ql.begin();
		quark_line::iterator it2 = ql.end() - 1;
		ql.erase(it2);
		ql.erase(it1);

		// This should increase the power of CF
		Monomial Mon_tmp;
		Mon_tmp.pow_CF = 1;
		Poly = Poly * Mon_tmp;

		// If now, all gluons are removed and there are no quarks,
		// increase Nc and open the ql (if not already open)
		// An open quark-line with no partons is defined as 1
		if (empty() && !open) {
			Monomial Mon_tmp;
			Mon_tmp.pow_Nc = 1;
			Poly = Poly * Mon_tmp;
			open = true;
		}
		j = j - 2;
	}
}

void  Quark_line::contract_neighboring_gluons( ){

	if( empty() ) return;

	if( open ){
		std::cerr << "uark_line::contract_neighboring_gluons( ): Expects a closed Quark_line, got "
				<< *this << std::endl;
		assert(0);
	}

	// Counting powers of CF,
	// Should be increases once for every gluon that is contracted
	// initially put to dummy value, minimal int
	int CF_pow_old=std::numeric_limits<int>::min();
	int CF_pow_new=0;

	// Keep looping over elements as long as CF_pow changes
	while( CF_pow_old!= CF_pow_new ){
		CF_pow_old=CF_pow_new;

		// Loop over elements in the Quark_line, starting at first parton
		// and ending with last
		for (uint j=0; j < size(); j++){
			contract_neighboring_gluons( j );
		}
		// To see if CF_pow changed, arbitrarily compare 0th element in Polynomial (if existing)
		// (all pow_CF changes)
		if( !Poly.empty() )
			CF_pow_new=Poly.at(0).pow_CF; // put CF_pow_new to current pow
	}

	return ;
}


void Quark_line::contract_next_neighboring_gluons( int j ) {

	if( ql.empty()) return;

	if( open ){
		std::cerr << "Quark_line::contract_next_neighboring_gluons( j): Expects a closed Quark_line, got "
				<< *this << std::endl;
		assert(0);
	}

	// If asked for -1st element, check last
	if(j==-1) j=ql.size()-1;
	// If asked for -2nd element, check second last
	if(j==-2) j=ql.size()-2;

	// Keep track of changes in number of gluons
	// add dummy number '1' to run at least once
	int old_size=ql.size()+1;

	// Loop over index in quark_line as long as the size keep changing,
	// or at least once
	while( old_size!=  static_cast<int>(ql.size()) &&  (ql.size())>=4){
		old_size=  (ql.size());
		// There are three cases:
		// the normal case, and, if the Quark_line is closed,
		// 2nd last =first (0th) and  last =2nd (1st)
		bool normal_case=(ql.size()>=4 && static_cast<uint>(j)< ql.size()-2 && ql.at(j)==ql.at(j+2));
		bool second_last_case=( ql.size()>=4 && !open && static_cast<uint>(j)==ql.size()-2 &&  ql.at(j)==ql.at(0) );
		bool last_case=(ql.size()>=4 && !open && static_cast<uint>(j)==ql.size()-1 &&  ql.at(j)==ql.at(1));
		// See if next to neighbors are equal
		// The point of having a while is to check again if one pair was found
		while(normal_case or second_last_case or last_case) {

			// Erase the gluons (at places depending on size)
			quark_line::iterator it1;
			quark_line::iterator it2;
			if(normal_case){
				it1=ql.begin()+j;
				it2=ql.begin()+j+2;
			}
			else if(second_last_case){
				it1=ql.begin();
				it2=ql.end()-2;
			}
			else if(last_case){
				it1=ql.begin()+1;
				it2=ql.end()-1;
			}
			// Erase right most gluon first
			ql.erase(it2);
			ql.erase(it1);

			// The surviving term comes with -TR/(Nc)
			Monomial Mon_tmp;
			Mon_tmp.pow_TR = 1;
			Mon_tmp.int_part=-1;
			Mon_tmp.pow_Nc=-1;
			Poly=Poly*Mon_tmp;

			// Check if new neighbors are equal
			uint old_size2=ql.size();
			// Check if inbetween gluon is now neighbor with itself to the right
			contract_neighboring_gluons( j );
			// Check if inbetween gluon is now neighbor with itself to the left
			contract_neighboring_gluons( j-1 );
			// if a neighboring pair was removed
			if(  (ql.size()) != old_size2 ) j=j-2;

			// Two gluons have been erased so to keep looking forward in the
			// Quark_line j must not be increased
			j--;

			normal_case=( ql.size()>=4 && static_cast<uint>(j)< ql.size()-2 && ql.at(j)==ql.at(j+2));
			second_last_case=( ql.size()>=4 && !open && j== static_cast<int>(ql.size())-2 && ql.at(j)==ql.at(0));
			last_case=( ql.size()>=4 && !open && j== static_cast<int>(ql.size())-1 &&  ql.at(j)==ql.at(1));
		}
	}
	return ;
}

void Quark_line::contract_next_neighboring_gluons(  ) {

	if( empty() ) return ;

	if( open ){
		std::cerr << "Quark_line::contract_next_neighboring_gluons: Expects a closed Quark_line, got "
				<< *this << std::endl;
		assert(0);
	}

	// Start with looking for neighbors
	contract_neighboring_gluons();

	// Keep track of changes in number of gluons
	// add dummy number '1' to run at least once
	int old_size= size()+1;

	// Loop over index in quark-line as long as the size keep changing,
	// or at least once
	while( old_size!=  static_cast<int>( size()) &&  ( size())>=4){
		// For seeing if size has changed
		old_size= size();
		for (int j=0; j<  static_cast<int>( size()); j++){
			// Contract next to neighbors starting at place j
			contract_next_neighboring_gluons( j );
		}
	}
	return ;
}




Quark_line operator*( const Quark_line & Ql, const int i){
	Quark_line Ql_out(Ql);
	Ql_out.Poly=Ql_out.Poly*i;
	return Ql_out;
}
Quark_line operator*(const int i, const Quark_line & Ql){
	return Ql*i;
}


Quark_line operator*(const Quark_line & Ql, const cnum c){
	Quark_line Ql_out(Ql);
	Ql_out.Poly=Ql_out.Poly*c;
	return Ql_out;
}
Quark_line operator*(const cnum c, const Quark_line & Ql){

	return Ql*c;
}


Quark_line operator*(const Quark_line & Ql, const double d){
	Quark_line Ql_out(Ql);
	Ql_out.Poly=Ql_out.Poly*d;
	return Ql_out;
}
Quark_line operator*(const double d, const Quark_line & Ql){

	return Ql*d;
}


Quark_line operator*(const Quark_line & Ql, const Monomial & Mon){

	Quark_line Ql_out(Ql);
	Ql_out.Poly=Ql_out.Poly*Mon;
	return Ql_out;
}
Quark_line operator*(const Monomial & Mon, const Quark_line & Ql){

	return Ql*Mon;
}


Quark_line operator*(const Quark_line & Ql, const Polynomial & Poly){

	Quark_line Ql_out(Ql);
	Ql_out.Poly=Ql_out.Poly*Poly;
	return Ql_out;
}
Quark_line operator*(const Polynomial & Poly, const Quark_line & Ql){

	return Ql*Poly;
}


std::ostream& operator<<(std::ostream& out, const Quark_line & Ql) {
	int max = Ql.ql.size();

	// Write out Polynomial if it is not 1, of if the ql has no structure
	Polynomial Poly1; // For comparison, the default Polynomial=1
	if (Ql.Poly != Poly1 or Ql.ql.size() == 0)
		out << Ql.Poly;

	if (max == 0 && Ql.open)
		out << "{}";
	else if (max == 0 && !Ql.open)
		out << "()";
	else if (max > 0 && Ql.open)
		out << "{";
	else if (max > 0 && !Ql.open)
		out << "(";
	for (int i = 0; i < max - 1; i++) {
		out << Ql.ql.at(i) << ",";
	}
	if (max > 0)
		out << Ql.ql.at(max - 1);
	if (max > 0 && Ql.open)
		out << "}";
	else if (max > 0 && !Ql.open)
		out << ")";

	return out;
}


} // end namespace ColorFull
