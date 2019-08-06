// -*- C++ -*-
/*
 *  Trace_basis.h
 *	Contains the declarations of the class Tree_level_gluon_basis, related types and operators
 *  Created on: Aug 9, 2012
 *  Author: Malin Sjodahl
 */

#ifndef COLORFULL_Tree_level_gluon_basis_h
#define COLORFULL_Tree_level_gluon_basis_h

#include "Trace_type_basis.h"

namespace ColorFull {


/// This class is for containing tree level gluon bases,
/// i.e. bases of form Tr(t^a t^b....t^z) +/- Tr(t^z .... t^b t^a),
/// where the sign is given by (-1)^n_g.
/// Technically only one of the quark-lines are carried
/// around in the cb member, the other is implicit.
/// This speeds up calculations.
class Tree_level_gluon_basis:public Trace_type_basis {
public:

	/// Default constructor.
	Tree_level_gluon_basis():Trace_type_basis(){
		initialize();
	}

	/// Constructor for creating a tree level gluon basis with n_g gluons.
	Tree_level_gluon_basis( int n_g ):Trace_type_basis(){
		initialize();
		create_basis( n_g );
	}

	/******************** Basis creation **********************/
	// Special functions are needed as half the color structure is implicit

	/// Creates a basis with basis vectors saved in the cb member.
	/// Each basis vector is a sum of two traces,
	/// of form Tr(t^a t^b....t^z) +/- Tr(t^z .... t^b t^a).
	/// The charge conjugated trace
	/// is implicit, and only one trace is actually
	/// carried around.
	void create_basis( int n_g );

	/******************** Basis reading and writing **********************/

	/// Function for reading in the basis from a file.
	/// The file should contain the whole basis,
	/// including the charge conjugated part.
	void read_in_Col_basis( std::string filename );

	/// Function reading in the basis from default name
	/// (see basis_file_name in the Col_basis class).
	/// The full basis, including the charge conjugated
	/// part should be contained in the file. (This is to simplify
	/// comparison with other programs, such as ColorMath.)
	void read_in_Col_basis( );

	/// Function for writing out the basis to default name,
	/// (see basis_file_name in the Col_basis class).
	/// The full basis, including the charge conjugated
	/// part is written out. (This is to simplify
	/// comparison with other programs, such as ColorMath.)
	void write_out_Col_basis( ) const;

	/// Function for writing out the basis to filename,
	/// (see basis_file_name in the Col_basis class).
	/// The full basis, including the charge conjugated
	/// part is written out. (This is to simplify
	/// comparison with other programs, such as ColorMath.)
	void write_out_Col_basis( std::string filename ) const;


protected:

	/// Function for writing out the basis in a human readable
	/// format to an ostream.
	virtual std::ostream&  write_out_Col_basis_to_stream( std::ostream&  out ) const;


private:

	/// Little helper function, called by constructors.
	void initialize(){
		tree_level_gluon_basis = true;
		max_ql=1;
	}

	/// Calculate element ij in scalar product matrix
	/// using the implicit presence of a charge conjugated Col_str.
	Polynomial ij_entry( const int i, const int j ) const;

	/******************** Basis creation **********************/

	/// Function for creating a basis with n_g gluons.
	Col_amp create_trace_basis( int n_g ) const;

	/// Computes the new basis states coming from one old Col_str
	/// when one gluon, g_new, is added.
	/// g_new is the number of the new gluon.
	Col_amp add_one_gluon( const Col_str & Cs, int g_new ) const;

	/// Compute the basis if one gluon is added to the old basis Old_basis.
	/// (For a Trace_type_basis, the basis is can be contained in a
	/// Col_amp, and this is the case for Old_basis).
	/// g_new is the number of the new gluon.
	/// Used by create_basis.
	Col_amp add_one_gluon( const Col_amp & Old_basis, int g_new ) const;

};

} /* namespace ColorFull */
#endif /* COLORFULL_Tree_level_gluon_basis_h */
