// -*- C++ -*-
/*
 * Trace_type_basis.cc
 * Contains definition of the class Trace_basis and associated types and operators.
 * Created on: Aug 9, 2012
 * Author: Malin Sjödahl
 */

#include "Trace_type_basis.h"
#include <cassert>
#include <iostream>


namespace ColorFull {


Poly_vec Trace_type_basis::decompose( const Col_amp & Ca ) {

	Col_amp Ca_copy=Ca;
	Ca_copy.simplify();

	// To contain the decomposed vector
	Poly_vec Decv;
	Polynomial Zero;
	Zero*=0;
	if(cb.size()==0){
		std::cerr << "Trace_type_basis::decompose: The basis vector cb is empty consider using create_basis or read in basis." << std::endl;
		assert( 0 );
	}
	else if(Ca_copy.size()>0 ){
		if(Ca_copy.at(0).n_quark() != nq)
		{
			std::cerr << "Trace_type_basis::decompose: The number of quarks in the argument Col_amp, " <<  Ca
					<< ", does not fit the number of quarks in the basis "
					<< nq << "."<< std::endl;}
		if(Ca_copy.at(0).n_gluon() != ng )
		{
			std::cerr << "Trace_type_basis::decompose: The number of gluons in the argument Col_amp, " << Ca
					<< ", does not fit the number of gluons in the basis "
					<< ng << "."<< std::endl;
		}
	}

	// Initially set all components of the vectors to 0
	for (uint m2 = 0; m2 < cb.size(); m2++){
		Decv.append(Zero);
	}

	// Loop over Cs in Ca and check which basis vector they equal
	for (uint m1 = 0; m1 < Ca_copy.ca.size(); m1++) {

		for ( uint m2 = 0; m2 < cb.size(); m2++ ) {
			bool found=false;

			// For a trace type basis it is enough to compare to see if the basis vector
			// is equal to the Col_str in the Ca.
			if (Ca_copy.ca.at(m1).cs == cb.at(m2).at(0).cs) {
				found=true;
				Decv.at(m2) += Ca_copy.ca.at(m1).Poly;
				Decv.at(m2).simplify();
			}

			if(found) break; // The tensor is ONE of the tensors in the Basis
		}
	}

	return Decv;
}

std::pair<int, int>  Trace_type_basis::new_vector_numbers( const Col_str & Cs, int emitter ) {

	// Number of new gluon
	int new_g=Cs.n_gluon()+2*Cs.n_quark()+1;

	// First check that basis is not empty
	if( cb.empty() ){
		std::cerr << "Trace_type_basis::new_vector_numbers: The basis has no vectors, "
					<< "consider using create_basis or read_in_basis." << std::endl;
		assert( 0 );
	}
	// Then check that number of quarks is the same as for the Cs
	if( nq != Cs.n_quark() ){
		std::cerr << "Trace_type_basis::new_vector_numbers: The number of quarks in the (new) basis, " << nq << " is not the same "
					<< "as the number of quarks in Cs, " << Cs.n_quark() << ", in Cs."<< std::endl;
		assert( 0 );
	}
	// Then check that number of gluons is one more than for the Cs
	if( ng != Cs.n_gluon()+1 ){
		std::cerr << "Trace_type_basis::new_vector_numbers: The number of gluons in the (new) basis, " << ng
					<< ", is not one plus the number of gluons in Cs << "
					<< Cs.n_gluon()<< std::endl;
		assert( 0 );
	}

	// To contain the numbers of the new non-zero vectors
	int plus_comp=-1;
	int minus_comp=-1;

	// The vector after emission
	Col_amp Ca=Col_fun.emit_gluon(Cs, emitter, new_g);
	Poly_vec vec=decompose( Ca );

	// Loop over vector entries to identify the non-zero ones
	for ( uint veci= 0; veci < cb.size(); veci++) {
		if( Col_fun. double_num(vec.at(veci))!=0 ){
			if(Col_fun.double_num(vec.at(veci)) > 0){plus_comp=veci;}
			if(Col_fun.double_num(vec.at(veci)) < 0){minus_comp=veci;}
		}
	}

	std::pair<int, int> res= std::make_pair(plus_comp,minus_comp);

	return res;
}


std::pair<int, int> Trace_type_basis::new_vector_numbers( int old_num, int emitter, int n_loop ) const{

	int n_g_old=ng-1;

	if( !(nq==0 or nq==1 or nq==2) or n_loop!=0 ){
		std::cerr << "Trace_type_basis:new_vector_numbers(int, int, int): Function only intended for special case of 0-2 q qbar pair at tree level. For the general case use the general version." << std::endl;
		assert( 0 );
	}

	// Find the place of the parton
	std::pair<int,int> parton_place=find_parton( emitter, old_num, nq, n_g_old, n_loop);

	// If the q or qbar is radiating, there is only one term
	if(parton_place.second==0 && nq>=1)// If the q is radiating
		return std::make_pair(new_vector_number(old_num, std::make_pair(parton_place.first,parton_place.second+1),  n_loop),-1);

	else if(parton_place.second==n_g_old+1 && nq==1)
		return std::make_pair(-1,new_vector_number(old_num, parton_place, n_loop));

	// If the emitter is a g there are two different terms
	std::pair<int,int> first_insertion_place, second_insertion_place;
	if( nq==0 && parton_place.second!=0) {
		first_insertion_place=std::make_pair(parton_place.first, parton_place.second + 1);
		second_insertion_place=std::make_pair(parton_place.first, parton_place.second );
	}
	// Special case that the emitter is standing first
	else if(nq==0 && parton_place.second==0){
		first_insertion_place=std::make_pair(0,  1);
		second_insertion_place=std::make_pair(0, n_g_old );
	}
	else if(nq==1){ // special cases already excluded
		first_insertion_place=std::make_pair(parton_place.first, parton_place.second+1);
		second_insertion_place=std::make_pair(parton_place.first, parton_place.second);
	}
	else if(nq==2){ // special case of quark already excluded

		// If the emitter is an anti-quark
		if(emitter==2 or emitter==4){
			return std::make_pair(-1,new_vector_number(old_num, parton_place, n_loop));

		}

		// Emitter is a gluon, and is thus not standing in the ends
		first_insertion_place=std::make_pair(parton_place.first, parton_place.second+1);
		second_insertion_place=std::make_pair(parton_place.first, parton_place.second);
	}

	int first_tens = new_vector_number(old_num, first_insertion_place, n_loop);
	int second_tens = new_vector_number(old_num, second_insertion_place, n_loop);

	return std::make_pair( first_tens, second_tens );
}



int Trace_type_basis::new_vector_number( int old_num, std::pair<int,int> place, int n_loop) const{

	int n_g_old=ng-1;

	std::cout.flush();
	if( !( nq==0 or nq==1 or nq==2 ) ){
		std::cerr << "Trace_type_basis::new_vector_number: Function only intended for special case of 1 or 2 open quark-lines." << std::endl;
		std::cerr.flush();
		assert( 0 );
	}
	if( nq==0 and !tree_level_gluon_basis ){
		std::cerr << "Trace_type_basis:new_vector_number: For 0 qqbar-pairs this function is only available for Tree_level_gluon_basis." << std::endl;
		assert( 0 );
	}
	if( nq>0 and !trace_basis ){
		std::cerr << "Trace_type_basis:new_vector_number: The basis type should be Trace_basis for processes with quarks, not Tree_level_gluon_basis." << std::endl;
		assert( 0 );
	}

	if(place.second == 0 && nq>0 ){
		std::cerr << "Trace_type_basis::new_vector_number: Cannot insert gluon at place 0, reserved for quark." << std::endl;
		std::cerr.flush();
		assert( 0 );
	} else if(place.second==n_g_old+2 && nq==1){
		std::cerr << "Trace_type_basis::new_vector_number: Cannot insert gluon at place n_g_old+2, reserved for qbar." << std::endl;
		std::cerr.flush();
		assert( 0 );}
	else if(nq==2 && (0 > place.second or place.second > n_g_old+1+1) ){
		std::cerr << "Trace_type_basis::new_vector_number: Cannot insert gluon outside range 1...n_g_old." << std::endl;
		std::cerr.flush();
		assert( 0 );
	}
	else if( 0 > place.second or place.second > n_g_old+1+nq){
		std::cerr << "Trace_type_basis::new_vector_number: Cannot insert gluon outside range 1...n_g_old ." << std::endl;
		std::cerr.flush();
		assert( 0 );
	}
	// If asked to insert a gluon at place 0, it should be inserted beyond the last place instead
	if( nq==0 && place.second==0 ){
		place.second=n_g_old;
	}

	// In the gluons only case we need to keep track of the next position relative to gluon 2
	int place2=0;

	int fnext=old_num; // To contain the remaining part of tensor number
	// If there are 2 quarks, there is an effective scalar place in {1,g1...gm,2}{3,gm+1...4}
	// obtained by crossing 2}{3 out. This place is given by
	// (place in the first ql) if it is in the first ql
	// (number of gluons in first ql) + place in second if it's in the second ql

	// If nq=2 the color structure is of type {q1, g1....gm,q2}{q3,gm+1...gn,q4}
	// This is almost as in the case of only one Ql {q1, g1....gm,gm+1...gn,q4}
	// with a breakpoint "q2}{q3" inserted or q4}{q3.
	// To find the tensor number, first locate the breakpoint.
	// This can be done by noting that there are 2*2*n_g_old! options for each breakpoints,
	// where the factor 2*2 comes from the options for the quarks
	// we thus have
	int break_after=0;
	int fsplit=0; // Contribution from where the ql is split
	int fq1=0; // Contribution from first quark
	if (nq == 2) {

		// To contain the number of gluon in the first ql
		// First the number of gluons in 2nd ql
		// For each gluon in 2nd ql there was
		// 2 (from nq)*2 (from n_qbar)*factorial(n_g_old) (from the gluons)
		// tensors with lower numbers from other breaks
		// Note that for all other splitting the factor from the first q is relevant
		break_after = old_num /( 2*2*Col_fun.factorial(n_g_old) ) ;


		fnext= old_num % (2*2*Col_fun.factorial(n_g_old));
		// The factor from the split
		fsplit=old_num-fnext;
		break_after=n_g_old-break_after;

		// Check if the first quark is 1 or 3
		// The factor from the first quark
		fq1=(fnext/( 2*Col_fun.factorial(n_g_old) ))*2*Col_fun.factorial(n_g_old);
		// Special case of 2 quark lines of equal length, then there is no option for q1 (it is 1)
		if( break_after==n_g_old-break_after ) fq1=0; // 1 is first, no factor from first quark
		fnext= fnext % (2*Col_fun.factorial(n_g_old));

		// First treat the special case when the ql's switch place
		// This can happen when
		// 1) The 2nd ql gets longer than the first such as
		// 		emission from 6 in {1,5,2}{3,6,4} -> {3,6,7,4}{1,5,2}
		// 2) The 2nd is one shorter than the first but has the q=1 standing first
		// 		{3,5,6,2}{1,7,4} ->{1,8,7,4}{3,5,6,2}
		if( ((break_after==n_g_old-break_after) && place.first==1) or ((break_after==n_g_old+1-break_after) && fq1>0 && place.first==1)) {

			// Locate all partons in old ql
			std::vector< std::pair<int, std::pair<int,int> > > parton_and_place;
			for( int p=1; p<=n_g_old+2*nq; p++ ){
				std::pair<int,int> p_place=find_parton(p, old_num, nq, n_g_old, n_loop);
				make_pair(p,p_place);
				// note that parton p stands at place p-1
				parton_and_place.push_back( make_pair(p,p_place) ) ;
			}

			// To contain new tensor number
			int new_num=0;

			// First calculate the number from the split
			// After the split the number of gluons in the qls are
			// break_after in the 2nd ql
			// This is true for both special cases of form
			// {1,5,2}{3,6,4} -> {3,6,7,4}{1,5,2}
			// 	{3,5,6,2}{1,7,4} ->{1,8,7,4}{3,5,6,2}
			new_num=(break_after)*2*2*Col_fun.factorial(n_g_old+1);

			// Then calculate the number from the first quark
			// There is a factor 2 from the anti-quark and a factor (n_g_old+1)! from the gluons
			// if the first factor is 3, which it is in swapings of form
			// {1,5,2}{3,6,4} -> {3,6,7,4}{1,5,2}
			if( (break_after==n_g_old-break_after) && place.first==1 )	new_num=new_num+2*Col_fun.factorial(n_g_old+1);
			// In cases of form {3,5,6,2}{1,7,4} ->{1,8,7,4}{3,5,6,2}
			// the factor from the first quark is 0, and


			// For the anti-quark qbar=2, check if it stands in first or 2nd ql
			// If it stands in first ql, there is no contribution to the tensor numbers
			// otherwise there is a contribution from all orders of gluons in 2nd ql
			if( parton_and_place.at(2-1).second.first==0) new_num+=Col_fun.factorial(break_after);

			// The 2nd anti-quark and 2nd quark gives no contribution

			// Then calculate the contribution to the new tensor number from the gluons
			// Loop over gluons
			for( int p=5; p<=n_g_old+2*nq; p++ ){

				// The NEW scalar place of the gluon
				int p_place=parton_and_place.at(p-1).second.second; // contribution from pos in ql
				// If a gluon was standing in the first ql it is afterwards standing in the 2nd
				// All gluons in the new first (old 2nd ql) are therefore standing in front
				// they are n_g_old-break_after +1 (the +1 is added below)
				if(parton_and_place.at(p-1).second.first==0) p_place+=(n_g_old-break_after);
				// If the new gluon was inserted before the parton, then the p_place should be increased with one
				// This happens when the parton stands in the new 2nd ql= old first ql  (always)
				// or when the gluon is in the new first ql, but after the emitter
				if( parton_and_place.at(p-1).second.first==0 or place.second <=parton_and_place.at(p-1).second.second ) p_place++;// ?? is this right

				// See how many smaller stand to the right
				int smaller_right_of_p=0;
				// Check all gluons with smaller number than p to see how many stand to the right
				for( int ps=5; ps<p; ps++ ){
					int ps_place=parton_and_place.at(ps-1).second.second; // Contribution from pos in ql
					if(parton_and_place.at(ps-1).second.first==0) ps_place+=(n_g_old-break_after); // From first ql
					if( parton_and_place.at(ps-1).second.first==0 or place.second <=parton_and_place.at(ps-1).second.second ) ps_place++; // From new g
					// If the smaller parton ps stand to the right of the parton p
					if( ps_place>p_place ) smaller_right_of_p++;
				}

				// Increase the new tensor number for partons to the right
				if( smaller_right_of_p>0 ){
					// If the gluon stands in the first ql (former 2nd), there is also a factor of 2 from the
					// choice of anti-quarks
					// All smaller_right_of_p ways of replacing p with a smaller number contributes a factor
					// Col_fun.factorial(gluons to the right), and there are smaller_right_of_p options
					if(parton_and_place.at(p-1).second.first==1) new_num+=2*smaller_right_of_p*Col_fun.factorial(n_g_old+1-p_place);
					else new_num+=smaller_right_of_p*Col_fun.factorial(n_g_old+1-p_place);
				}

			} // End of loop over gluons

			// Add contribution from the new gluon. The gluon number is larger than every other gluon number,
			// so all gluons standing to the right are smaller
			// In the new first ql the new gluon is at place place.second
			// and there are in total n_g_old-break_after+1 gluons, so the total number of gluons to the right is
			// break_after+(n_g_old-break_after+1) -place.second
			int n_right_new=(n_g_old+1)-place.second;
			// The overall contribution also has a factor 2 from the qbar choice
			// The factor 2 is from the qbar, Col_fun.factorial(n_right_new) from ordering gluons
			new_num+=2*n_right_new*Col_fun.factorial(n_right_new);

			return new_num;
		} // End of hard case
	}


	if ( nq==0 ) place2=find_parton(2, old_num, nq, n_g_old, n_loop).second;
	// Initially we had q=1, g1,....gn, qbar=2
	// After the insertion at place place we have q=1, g1, ...g(place-1), g_new, g(place), ...gn, qbar=2
	// The number of the new vector can be seen as a sum of 3 factors:
	// f1 = the contribution from q=1, g1, ...g(place-1)
	// fnew= the contribution from the newly inserted gluon
	// f2= the contribution from g(place), ...gn, qbar=2 (this contribution is the same as before the insertion)
	// Of these contributions fnew is easiest to calculate it is simply:
	int fnew=0;
	// If one ql
	if( nq==1 ) fnew=(n_g_old - place.second + 1)*Col_fun.factorial(n_g_old - place.second + 1);
	// Special case of place after last gluon
	else if( nq==0 && place.second == n_g_old+1 ) fnew=0;
	// If one closed ql and the parton is to the right of 2
	else if( nq==0 && place.second > place2 && place.second < n_g_old+1) fnew=(n_g_old - place.second )*Col_fun.factorial(n_g_old - place.second );
	// If the parton is before 2
	else if( nq==0 && place2 >= place.second ) fnew=(n_g_old - place.second )*(Col_fun.factorial(n_g_old - place.second )/2); //Rel order 2 3 fixed

	// The factor f1 is trickier to obtain. Its initial value (before emission) is
	// sum_{i=1}^(place-1) n_sri (n_p-i)!
	// where n_sri denotes the number of smaller partons standing to the rigt at place i
	// these factors are not known, but can be calculated.
	// Once these factors are known the final value of f1 can also be obtained from
	// sum_{i=1}^(place-1) n_sri (n_p-i+1)!
	// where the extra 1 in the factorial is present as after the emission there is
	// one more parton standing to the right.

	int f1new=0; // To contain the new f1
	int f1old=0; // To contain the old f1

	int scalar_place=0; // The scalar place
	if( place.first!=0 )  scalar_place=break_after;
	scalar_place+=place.second;

	if( nq==2 ) {
		fnew=(n_g_old + 1- scalar_place )*Col_fun.factorial(n_g_old + 1- scalar_place );

		// If g inserted in first ql, there is an extra factor of 2 from interchanging quarks
		if( place.first==0 ) fnew*=2;

	} // end of if( nq==2 )

	// The last relevant parton for calculating the contribution to f1
	int last_rel=scalar_place-1;

	for (int i=1; i<=last_rel; i++){
		//cout << "i " << i << endl;
		// number of partons to the right
		int n_r=n_g_old-i;
		if( nq==0 ) n_r--;

		// nsri is the number of smaller partons to the right
		int nsri=0;
		if (nq==1 ) {
			// Note integer division, how many times do we find the factorial
			nsri=fnext/(Col_fun.factorial(n_r));
			// The contribution to the old f1
			f1old+=nsri*Col_fun.factorial(n_r);
			// The next effective f, to use for finding next nsri
			fnext=fnext%(Col_fun.factorial(n_r));
			// The contribution to the new f1
			f1new+=nsri*(Col_fun.factorial(n_r+1));
		}
		else if (nq==2) {
			int q_part=1;
			// Keeps track of a possible multiplicative factor from first qbar option
			if ( i <= break_after) q_part=nq;

			// Note integer division, how many times do we find the factorial
			nsri=fnext/(q_part*Col_fun.factorial(n_r));
			// The contribution to the old f1
			f1old+=nsri*q_part*Col_fun.factorial(n_r);
			// The next effective f, to use for finding next nsri
			fnext=fnext%(q_part*Col_fun.factorial(n_r));
			// The contribution to the new f1
			f1new+=nsri*(q_part*Col_fun.factorial(n_r+1));

			// Take care of factor from first qbar "when we reach it"
			// if( i == break_after && place.second!=break_after+1 && break_after=n_g_old-break_after) {
			if( i == break_after ) {
				// Note integer division, how many times do we find the factorial
				nsri=fnext/(Col_fun.factorial(n_r));
				int f1old_before_qbar=f1old;
				// The contribution to the old f1
				f1old+=nsri*Col_fun.factorial(n_r);
				// The next effective f, to use for finding next nsri
				fnext=fnext%(Col_fun.factorial(n_r));
				// if the gluon is inserted just before the qbar the factor from
				// the qbar should be taken care of in the f2 instead
				// (f1old still has to be recalculated as it is used for finding f2 which is unchanged)
				// If the g is NOT inserted just before the qbar
				// i.e. if it is either inserted in 2nd ql or at other place than break_after+1
				if ( !(scalar_place == break_after+1) or place.first ==1){
					// The contribution to the new f1 from the anti-quark
					//f1new+=nsri*(factorial(n_r+1));
					// There are (n_r+1)! options to have 2 before 4
					// and we should only have a contribution if 4 is before 2
					f1new+=nsri*(Col_fun.factorial(n_r+1));
				}
				else { // If the gluon is inserted just before the qbar
					// The factor should be accounted for in fnew instead
					// If the second qbar was 4, this should be compensated for in fnew
					if( f1old_before_qbar!=f1old && place.first ==0){

						// The factor is the number of ways of ordering the gluons in 2nd ql
						// The q and qbar are fixed and do not contribute
						// There is no factor if no gluons to the right ??
						// There are n_r gluons in the second ql (also after emission)
						fnew=fnew+Col_fun.factorial(n_r);
					}
				}
			}

		} // end if (nq==2)

		// If i is after the place of 2 we have the standard case
		else if(nq==0 && i >= place2) {
			nsri=fnext/(Col_fun.factorial(n_r));
			f1old+=nsri*Col_fun.factorial(n_r);
			fnext=fnext%(Col_fun.factorial(n_r));
			f1new+=nsri*(Col_fun.factorial(n_r+1));
		}
		// If i is before the place of 2 ordering of 2 and 3 irrelevant
		else if(nq==0 &&  i < place2) {

			nsri=fnext/(Col_fun.factorial(n_r)/2);
			f1old+=nsri*(Col_fun.factorial(n_r)/2);
			fnext=fnext%(Col_fun.factorial(n_r)/2);
			f1new+=nsri*(Col_fun.factorial(n_r+1)/2);
		}
	}

	// Once f1 is known then f2 is simply,
	int f2=0;
	if( nq==2 ) {
		// For 2 qqbar pairs we have to take care of factor from split
		f2=old_num-f1old-fsplit-fq1;
	}
	else f2=old_num-f1old;

	// and the contribution from f2 is the same as before (irrespectively of position relative 2)

	// In the case of 2 quarks, if there was a contribution to the tensor number from the split
	// this has to be added back

	// The contribution was (n_g_old-break_after)*(2*factorial(n_g_old))
	// The special case of ql swapping should already have been taken care of
	if( nq==2 ) {

		// Calculating new factor from the split
		int fsplitnew=0;

		// New factor from split after insertion in first ql
		if(place.first==0) {
			// The length of 2nd ql is then still (n_g_old-break_after)
			// and for each break there are 2*factorial(n_g_old+1) options
			// The lengths of first and 2nd ql are never equal as 1st was at least as long as 2nd
			// -> factor 2 from option of first q
			fsplitnew=(n_g_old-break_after)*2*2*Col_fun.factorial(n_g_old+1);
		}

		// If the insertion was in the second ql (of new length n_g_old-break_after+1)
		if(place.first==1) {
			// The length of the 2nd ql is now (n_g_old+1-break_after)
			// all gluon orders matter, and the qbar matters
			// Also the factor fomr the split, coming from all "lower splits"
			// has a factor 2 as the q matters for the lower splits
			fsplitnew=(n_g_old+1-break_after)*2*2*Col_fun.factorial(n_g_old+1); // 1*2*2!

		}

		int fq1new=fq1*(n_g_old+1);

		// If after emission the 2nd ql has the same length as the first
		if(place.first==1 && break_after==n_g_old+1-break_after){
			fq1new=fq1new/2; // How can this be right?? There should be no factor from first quark then
		}

		f1new=fsplitnew+fq1new+f1new;

	} // End of nq==2

	return f1new+fnew+f2;
}


std::pair<int,int> Trace_type_basis::find_parton( int parton, int vec_num, int n_quark, int n_gluon, int n_loop ) const{
	if( !(n_quark==1 or n_quark==0 or n_quark==2) or n_loop >0 ){
		std::cerr << "Trace_type_basis:find_parton: Function only intended for special case of 0-2 qqbar pairs, and a tree level bases." << std::endl;
		assert( 0 );
	}
	if( n_quark==0 and !tree_level_gluon_basis ){
		std::cerr << "Trace_type_basis:find_parton: For 0 qqbar-pairs this function is only available for Tree_level_gluon_basis." << std::endl;
		assert( 0 );
	}
	if( n_quark>0 and !trace_basis ){
		std::cerr << "Trace_type_basis:find_parton: The basis type should be Trace_basis for processes with quarks, not Tree_level_gluon_basis." << std::endl;
		assert( 0 );
	}

	if( parton==1 && n_quark == 1) return std::make_pair( 0, 0 ); // q=1 or g=1 has position 0 if one ql
	if( parton==1 && n_quark == 0) return std::make_pair( 0, 0 ); // q=1 or g=1 has position 0 if one ql
	if( parton==2 && n_quark==1 ) return std::make_pair( 0, n_gluon+1 ); // qbar has last pos

	// If n_quark=2 the color structure is of type {q1, g1....gm,q2}{q3,gm+1...gn,q4}
	// This is almost as in the case of only one Ql {q1, g1....gm,gm+1...gn,q4}
	// with a breakpoint "q2}{q3" inserted or q4}{q3.
	// To find the tensor number, first locate the breakpoint.
	// This can be done by noting that there are 2*n_gluon! options for each breakpoints,
	// where the factor 2 comes from the 2 options for the first anti-quark
	// we thus have
	int break_after=n_gluon;
	if (n_quark == 2) {
		// There are n_quark!n_quark!n_gluon! assignments for each break
		break_after = break_after-(vec_num) / (2*2*Col_fun.factorial(n_gluon));
		vec_num= vec_num -(n_gluon-break_after)*(2*2*Col_fun.factorial(n_gluon));

		// Check if first parton is 1 or 3, there are 2 options for each from the qbar being 2 or 4
		bool first3=false;
		if (  ( vec_num / ( 2*Col_fun.factorial(n_gluon)) ) ) first3=true;

		vec_num= vec_num % (2*Col_fun.factorial(n_gluon));

		// Special cases of asked for quark
		if(parton == 1){
			if(!first3) return std::make_pair(0,0);
			else return std::make_pair(1,0);
		}
		if(parton == 3){
			if(first3) return std::make_pair(0,0);
			else return std::make_pair(1,0);
		}
	}

	// The number of partons gluons with lower parton number than the parton
	// (and possibly standing to the right, hence -1)
	int lower=parton-2-n_quark;
	if ( n_quark==2 ) lower=parton -2*n_quark-1; // quarks don't count

	// To contain the number of gluons standing to the right of the parton under consideration
	// and being smaller than 2, noting to the right is smaller than 2, hence 0
	// needed for special case of gluons only
	int lower2=0;

	// Number of partons to the right, when checking the various positions
	int n_r=0;
	// For quark lines(s) all gluons except the one under consideration matter
	if ( n_quark>0 ) n_r=n_gluon-1;
	// For a gluon line we have to subtract 1 extra as we start read at place 1 anyway
	else if ( n_quark==0 ) n_r=n_gluon-2;

	// Was the gluon 2 found or not, needed for special case of gluons only
	bool found_2=false;

	// The resulting quark_line number (can change only in case 2 q qbar pairs)
	int ql_num=0;

	// Loop over possible places of the gluon, start with checking place 1
	// n_l_r contains the number of partons lower than the parton at the position under consideration
	int n_l_r=0; // To make compiler not warn
	int looped_over=0;
	while(true){
		if( n_quark==1 ) {
			n_l_r=vec_num/Col_fun.factorial(n_r); // The result of integer division
			vec_num=vec_num%Col_fun.factorial(n_r); // The rest after integer division
		}

		if( n_quark==2 ) {

			// If both qbars to right extra factor 2!
			if( break_after > looped_over) {

				n_l_r=vec_num/(Col_fun.factorial(n_r)*2); // The result of integer division
				vec_num=vec_num%(Col_fun.factorial(n_r)*2); // The rest after integer division
			}
			else{
				// Check if it is time to compensate for first qbar and second q
				if ( break_after==looped_over ){

					// If we have already looped over all gluons
					if(n_r==-1) n_r=0;
					// Special cases if asked for first qbar = 2 or 4
					if( parton==2){
						// All gluons which are not looped over are to right
						if (vec_num/(Col_fun.factorial(n_r+1))==0) return std::make_pair(0,break_after+1);
						else return std::make_pair(1,n_gluon-break_after+1);
					}
					if( parton==4){
						if (vec_num/(Col_fun.factorial(n_r+1))==1) return std::make_pair(0,break_after+1);
						else return std::make_pair(1,n_gluon-break_after+1);
					}

					vec_num=vec_num%(Col_fun.factorial(n_r+1)); // compensated for first qbar


					// After looking at first qbar the ql_num should be 1
					ql_num=1;
				}

				// Otherwise the situation is as for gluons only
				n_l_r=vec_num/(Col_fun.factorial(n_r)); // The result of integer division
				vec_num=vec_num%Col_fun.factorial(n_r); // The rest after integer division
			}
		}
		else if (n_quark==0){
			if( found_2 ){
				n_l_r=vec_num/Col_fun.factorial(n_r);
				vec_num=vec_num%Col_fun.factorial(n_r);
			}
			else{
				// First check if 2 was at the position
				n_l_r=vec_num/(Col_fun.factorial(n_r)); // If 2 was at the place, then all orders matter
				// If 2 was found
				if( n_l_r==lower2 ) {
					found_2=true;
					vec_num=vec_num%(Col_fun.factorial(n_r));				}
				else{ // Recalculate knowing that we did not find 2
					n_l_r=vec_num/(Col_fun.factorial(n_r)/2); // There is only one relative order of gluon 2 and 3
					vec_num=vec_num%(Col_fun.factorial(n_r)/2); // The rest after integer division
				}
			}
		}
		if( n_l_r==lower ) break; // Number smaller to right = number smaller than parton to right
		if( n_l_r < lower ) lower--;
		n_r--;
		looped_over++; // needed to compare with break in 2q case
	}

	// and then the number in the quark_line
	int pos_num=n_gluon-n_r; // n_quark==1 case
	if( n_quark==0) pos_num=n_gluon-n_r-1;
	else if( n_quark==2 && ql_num==0) pos_num=n_gluon-n_r;
	else if(n_quark==2 && ql_num==1) pos_num=n_gluon-break_after-n_r;

	return std::make_pair( ql_num, pos_num );
}



} // end namespace ColorFull
