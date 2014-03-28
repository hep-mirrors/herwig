// -*- C++ -*-
/*
 * Trace_basis.h
 * Contains the declarations of the class Trace_basis, related types and operators
 * Created on: Aug 9, 2012
 * Author: Malin Sjodahl
 */

#ifndef COLORFULL_Trace_basis_h
#define COLORFULL_Trace_basis_h

#include "Trace_type_basis.h"

namespace ColorFull {


/// In a trace basis the basis vectors are closed or open quark-lines
/// or products of close and open quark-lines.
class Trace_basis:public Trace_type_basis {
public:

	/// Default constructor.
	Trace_basis():Trace_type_basis(){
		initialize();
	}

	/// Constructor for creating a trace basis for n_quark qqbar-pairs
	/// and n_gluon gluons.
	/// (Note: For electroweak interactions more color structures may
	/// be needed.)
	Trace_basis( int n_quark, int n_gluon ){
		initialize();
		create_basis( n_quark, n_gluon );
	}

	/// Constructor for creating a trace basis for n_quark qqbar-pairs
	/// and ng gluons, keeping only those color structures that
	/// can appear to order n_loop in QCD.
	/// (Note: For electroweak interactions more color structures may
	/// be needed.)
	Trace_basis( int n_quark, int n_gluon, int n_loop ){
		initialize();
		create_basis( n_quark, n_gluon, n_loop );
	}

	/// Little helper function, called by all constructors.
	void initialize(){
		nq=0;
		ng=0;
		tree_level_gluon_basis = false;
		orthogonal_basis = false;
		trace_basis = true;
	}

	/******************** Functions for basis creation **********************/

	/// Creates a trace basis with basis vectors saved in the cb member.
	/// Keeps all possible basis vectors, i.e., the basis is valid to any order
	/// in perturbation theory.
	void create_basis( int n_q, int n_g );

	/// Creates a trace basis with basis vectors saved in the cb member.
	/// Keeps only basis vectors needed up to n_loop in pure QCD.
	void create_basis( int n_q, int n_g, int n_loop );


private:

	/******************** Internal function for basis creation **********************/

	/// Function for creating a basis with n_q=n_qbar
	/// quarks and ng gluons to order n_loop,
	/// i.e. each Col_str is a product of at most
	/// max(1+nq)+n_loop Quark_lines.
	Col_amp create_trace_basis( int n_q, int n_g, int n_loop ) const;

	/// Connect the n_q quarks in all n_q! ways.
	/// This function is used when creating a basis.
	Col_amp connect_quarks( int n_quark ) const;

	/// Compute the new basis states coming from
	/// one old vector when one gluon, g_new, is added.
	Col_amp add_one_gluon( const Col_str & Cs, int g_new, int n_loop ) const;

	/// Compute the basis if one gluon is added to the old basis Old_basis.
	/// If n_loop==0, only tree level states are constructed.
	Col_amp add_one_gluon( const Col_amp & Old_basis, int n_q, int n_g, int g_new, int n_loop ) const;

};

} /* namespace ColorFull */
#endif /* COLORFULL_Trace_basis_h */
