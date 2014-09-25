// -*- C++ -*-
/*
 * Trace_basis.cc
 * Contains the definitions of the class Trace_basis, related types and operators
 * Created on: Aug 9, 2012
 * Author: Malin Sjodahl
 */

#include "Trace_basis.h"
#include <cassert>
#include <iostream>


namespace ColorFull {


void Trace_basis::create_basis(int n_quark, int n_gluon){

	// Create the basis with an upper bound for n_loop to
	// keep all basis vectors
	// Maximal number of "loops"
	int n_l=n_gluon/2; // Integer division
	create_basis( n_quark, n_gluon, n_l);
}


void Trace_basis::create_basis(int n_quark, int n_gluon, int n_loop){

	// Setting basis variables
	nq=n_quark;
	ng=n_gluon;

	// If only to finite "loop" the maximal number of qls change
	max_ql=std::max(nq,1)+n_loop;

	// Remove a potentially already calculated basis
	cb.clear();
	// The Col_amp containing the basis
	Col_amp Ca_basis;

	// Create the basis using the old function for a maximal number of loops
	Ca_basis=create_trace_basis(nq, ng, n_loop);

	// Sort the resulting Col_amp into the Col_basis cb
	for ( uint i=0; i < Ca_basis.ca.size(); i++ ){
		// A Col_amp to contain a basis vector
		Col_amp Ca_vec;
		Ca_vec.ca.push_back(Ca_basis.ca.at(i));
		cb.push_back(Ca_vec);
	}
}


Col_amp Trace_basis::create_trace_basis( int n_quark, int n_g, int n_loop ) const {

	// To contain the resulting basis
	Col_amp Basis;

	// First connect all q and qbars in all n_quark! ways
	// This gives a basis for the case of n_quark quarks and n_quarkbar anti-quarks
	if( n_quark != 0)  Basis= connect_quarks( n_quark );

	// If there were no quarks
	if( n_quark==0 ){
		// There has to be at least two gluons
		if( n_g <= 1 ) {
			std::cerr << "Trace_basis::create_trace_basis: For 0 quarks there is no basis with " << n_g << " gluons" << std::endl;
			assert( 0 );
		}
		// If 2 or more gluons, build from the 2-gluon basis
		else if( n_g>=2 ){
			Col_str Cs_OnlyState("[(1,2)]");
			Col_amp Ca_tmp;
			Ca_tmp.ca.push_back(Cs_OnlyState);

			Basis= Ca_tmp;
			// For 2 gluons, the work is done
			if( n_g==2 ) return Basis;
			//cout << Basis;
		}
	}

	// Then, add the gluons one at the time
	// If there are only gluons the generation should start from gluon 3,
	// otherwise from 2*n_quark+1;
	for( int g_new=std::max(3, 2*n_quark+1); g_new <= 2*n_quark+ n_g; g_new++){
		// If only gluons start from the 2-gluon state, so add gluon 3
		Basis=add_one_gluon(Basis, n_quark ,n_g, g_new, n_loop);
		}

	// Remove basis vectors with too many Quark_lines
	uint max_Ql=std::max(1, n_quark)+n_loop;
	for(int i=Basis.size()-1; i>-1; i--){
		if (Basis.ca.at(i).size() > max_Ql ) Basis.erase(i);
	}

	// Normal order the Col_str's
	Basis.normal_order();

	return Basis;

}


Col_amp Trace_basis::connect_quarks( int n_quark ) const{

	Col_amp Basis;

	// If n_quark=1, there is only one state
	if( n_quark==1 )	{
		Col_str Cs_tmp("[{1,2}]");
		Basis=Basis+Cs_tmp;
	}
	else{ // For more than one q qbar pair, build up states iteratively
		// Start from the state with one less q qbar state
		Col_amp Old_bas=connect_quarks( n_quark-1 );

		// Loop over old basis vectors
		for (uint old_v = 0; old_v < Old_bas.ca.size(); old_v++) {

			// Some, 1*(n_quark-1)!, new states are obtained by just adding q_new qbar_new to the old states
			// Construct the new Ql with i.e. {q_new, qbar_new}
			Quark_line Ql_new;
			Ql_new.open = true;
			Ql_new.ql.push_back(2 * n_quark - 1);
			Ql_new.ql.push_back(2 * n_quark);
			Col_str New_state;
			New_state.cs.push_back(Ql_new);
			New_state.append( Old_bas.ca.at(old_v).cs );
			Basis=Basis+New_state;
		}

		// The other (n_quark-1)*(n_quark-1)! states are obtained by combining the new qbar with any old quark
		// Loop over old basis states
		for (uint old_v = 0; old_v < Old_bas.ca.size(); old_v++) {
			// Loop over old qbar indices
			for (int qbar_old = 2; qbar_old < 2 * n_quark; qbar_old += 2) {
				// To contain the new basis vector
				Col_str New_state = Old_bas.ca.at(old_v);

				// Replace the index qbar_old with 2*nq
				New_state.replace(qbar_old, 2 * n_quark);

				// Create the Quark_line {q_new, qbar_old}
				Quark_line Ql_new;
				Ql_new.open = true;
				Ql_new.ql.push_back(2 * n_quark - 1);
				Ql_new.ql.push_back(qbar_old);

				// Add the Col_str from  {q_new, qbar_old} to the new basis vector
				New_state.cs.push_back(Ql_new);

				// Add the new basis tensor to the basis
				Basis = Basis + New_state;
			}
		}
	}

	return Basis;
}


Col_amp Trace_basis::add_one_gluon( const Col_str & Cs, int g_new, int ) const {

	// For storing the new basis
	Col_amp New_tensors;

	// Add the new gluon in all possible ways to the old Color structure
	// Loop over the Quark_lines
	for (uint ql = 0; ql < Cs.cs.size(); ql++) {

		// The old Quark_line, before insertion of the new gluon index
		Quark_line Old_Ql = Cs.cs.at(ql);
		Col_str New_tensor = Cs;

		// Loop over (potential) insertion places in the Quark_lines, starting from the end
		for (int j = Old_Ql.ql.size(); j > 0; j--) {

			// Special treatment of last place, insert here only if the ring is open
			// (the gluon index cannot take the place of the a quark index)
			if (Old_Ql.open && j == static_cast<int>(Old_Ql.ql.size()) ) j--;

			Quark_line New_Ql = Old_Ql;
			quark_line::iterator it = New_Ql.ql.begin() + j;
			New_Ql.ql.insert(it, g_new);

			// Replace the old Col_str with the new and add to the new basis states
			New_tensor.cs.at(ql) = New_Ql;
			New_tensors = New_tensors + New_tensor;
		}

	}
	return New_tensors;
}


Col_amp Trace_basis::add_one_gluon(const Col_amp & Old_basis, int n_q, int, int g_new, int n_loop) const{

	// For storing the new basis
	Col_amp New_bas;

	// Add the new gluon to each of the previous color tensors
	for (uint t = 0; t < Old_basis.ca.size(); t++) {
		New_bas = New_bas + add_one_gluon(Old_basis.ca.at(t), g_new, n_loop);
	}

	// Create new states with 2-rings formed from taking the new gluon and any old gluon
	// these states can be created by taking the ng-2 basis but use the indices
	// 1...i_ex-1, i_ex+1 ng-1
	// for the gluons, that is exclude one index (to be combined with the index ng in a 2-ring)
	// First, create the ng-2 basis here, must have at least 4 gluons
	// Do this only is not tree-level mode, i.e. n_loop>0
	if ( n_loop>0 && (( n_q>0 && g_new-2*n_q>=2 ) or (n_q==0 && g_new >= 4) )) {
		Col_amp Bas_n_g_minus_2 = create_trace_basis(n_q, g_new-n_q*2-2, n_loop);

		// Loop over indices to group with the new gluon index
		for (int g_ex = 2 * n_q + 1; (g_ex <= g_new - 1 ); g_ex++) {

			// Create the 2-ring with g_ex and g_new
			Quark_line Ql2_ring;
			Ql2_ring.open = false;
			Ql2_ring.append(g_ex);
			Ql2_ring.append(g_new);

			// Loop over Col_str's=basis vectors
			for (uint bv = 0; bv < Bas_n_g_minus_2.ca.size(); bv++) {

				Col_str Cs_new; // To contain the new basis vector
				Cs_new = Bas_n_g_minus_2.ca.at(bv);

				// Loop over the Quark_lines in the basis vectors
				for (uint ql = 0; ql < Bas_n_g_minus_2.ca.at(bv).cs.size(); ql++) {
					// Loop over positions in the quark_line
					for (uint pos = 0; pos < Bas_n_g_minus_2.ca.at(bv).cs.at(ql).ql.size(); pos++)

						// Jump over index g_ex by increasing the index g_ex, and all indices above with 1
						if (Bas_n_g_minus_2.ca.at(bv).cs.at(ql).ql.at(pos)
								>= g_ex) {
							// Change the index if >= g_ex
							Cs_new.cs.at(ql).ql.at(pos)++;
						}
				}
				Cs_new.cs.push_back(Ql2_ring);
				New_bas = New_bas + Cs_new;
			}
		}
	}
	return New_bas;
}



}// end namespace ColorFull {


