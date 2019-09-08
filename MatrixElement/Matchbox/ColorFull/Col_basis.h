// -*- C++ -*-
/*
 * Col_basis.h
 * Contains the declarations of the class Col_basis, related types and operators.
 * Created on: Aug 9, 2012
 * Author: Malin Sjodahl
 */

#ifndef COLORFULL_Col_basis_h
#define COLORFULL_Col_basis_h


#include "Col_functions.h"


namespace ColorFull {

/// Define a type to store all basis vectors.
typedef std::vector<Col_amp> col_basis;

/// To contain a color basis, where each basis vector is a Col_amp.
/// Technically the color information is contained
/// in the member cb, which is a vector of Col_amps, a col_basis.
class Col_basis {

public:

	/// Default constructor
	Col_basis(){
		nq=0;
		ng=0;
		trace_basis = false;
		tree_level_gluon_basis = false;
		orthogonal_basis = false;
	}

	/// Destructor.
	virtual ~Col_basis(){}

	/// The number of qqbar-pairs, initially set to 0 and (possibly)
	/// changed with create_basis, or while reading in the basis.
	int nq;

	/// The number of gluons, initially set to 0 and (possibly)
	/// changed with create_basis, or while reading in the basis.
	int ng;

	/// To actually contain the info about the basis vectors cb= vector1, vector2... .
	/// Technically cb is a col_basis, a vector of Col_amps.
	col_basis cb;

	/// To contain the Polynomial version of the scalar product matrix.
	Poly_matr P_spm;

	/// To contain the Polynomial version of the leading part of the scalar product matrix.
	Poly_matr leading_P_spm;

	/// To contain the double version of the scalar product matrix.
	dmatr d_spm;

	/// To contain the double version of the leading part of the scalar product matrix.
	dmatr leading_d_spm;

	/// To contain the set of Col_functions used.
	Col_functions Col_fun;

	/// Is it a Trace_basis?
	bool  is_Trace_basis() const {return trace_basis;}

	/// Is it an Orthogonal_basis?
	bool  is_Orthogonal_basis() const {return orthogonal_basis;}

	/// Is it a Tree_level_gluon_basis?
	bool  is_Tree_level_gluon_basis() const {return tree_level_gluon_basis;}

	/// Returns the number of basis vectors.
	uint size() const {return cb.size();}

	/// Returns the Col_amp (basis vector) at place i.
	const Col_amp & at( const int i ) const{return cb.at(i);}

	/// Returns the Col_amp (basis vector) at place i.
	Col_amp & at( const int i ) {return cb.at(i);}

	/// Is the col_basis empty?
	bool empty() const { return cb.empty(); }

	/// Erase the basis, stored in cb.
	void clear() { cb.clear(); }

	/// Appends a Col_amp to the basis, stored in cb.
	void append( Col_amp Ca ) { cb.push_back( Ca ); }

	/// Appends the Col_amps in cb_in to the col_basis member cb.
	void append( col_basis cb_in );


	/******************** Functions for scalar products **********************/


	/// Function for calculating the scalar products matrix.
	/// This function loops over all basis vectors and stores the
	/// value of the scalar product between basis vector
	/// i and basis vector j in the i,j -entry in P_spm and d_spm.
	/// The calculation is done using memoization.
	/// The symmetry of the scalar product matrix is not used
	/// for the calculation, instead it is checked that the
	/// resulting matrix is indeed symmetric.
	void scalar_product_matrix();

	/// This function works as scalar_product_matrix, but does the
	/// calculation numerically. It hence only calculates d_spm.
	void scalar_product_matrix_num();

	/// This function works like scalar_product_matrix, but does
	/// not use memoization.
	virtual void scalar_product_matrix_no_mem();

	/// This function works like scalar_product_matrix_num, but does
	/// not use memoization.
	void scalar_product_matrix_num_no_mem();

	/// Finds the leading Nc scalar product matrices,
	/// leading_P_spm and leading_d_spm.
	/// If the polynomial scalar product matrix, P_spm has
	/// been calculated, P_spm is used, otherwise P_spm is first calculated
	/// and the leading Nc limit is then taken of P_spm.
	void leading_scalar_product_matrix();

	/// Function for calculating scalar products algebraically
	/// using the basis and the scalar product matrix (Poly_matr) in the basis.
	/// (Does add implicit conjugated part for Tree_level_gluon_basis,
	/// as these terms are contained in the matrix of scalar products.)
	virtual Polynomial scalar_product( const Col_amp & Ca1, const Col_amp & Ca2 );

	/// Function for calculating scalar products numerically, knowing the basis
	/// and the scalar product matrix in numerical form.
	/// (Does add implicit conjugated part for Tree_level_gluon_basis,
	/// as these terms are contained in the matrix of scalar products.)
	virtual cnum scalar_product_num( const Col_amp & Ca1, const Col_amp & Ca2 );

	/// Calculates the scalar product between numerical (complex) amplitudes v1, V2
	/// using the numerical scalar product matrix, d_spm.
	/// The vectors thus needs to be expressed in the basis contained in cb.
	/// (Does add implicit conjugated part for Tree_level_gluon_basis,
	/// as these terms are contained in the matrix of scalar products.)
	virtual cnum scalar_product_num( const cvec & v1, const cvec & v2 );

	/// Calculates the scalar product between numerical (complex) amplitudes v1, v2
	/// using the numerical scalar product matrix, d_spm.
	/// Assumes that there are only diagonal contributions.
	/// This is useful for calculations in leading Nc limit.
	/// (Does add implicit conjugated part for Tree_level_gluon_basis,
	/// as these terms are contained in the matrix of scalar products.)
	cnum scalar_product_num_diagonal( const cvec & v1, const cvec & v2 );


	/******************** Functions for reading and writing **********************/

	// Functions for naming files

	/// Returns a standard filename, used for writing out
	/// the basis to a file.
	std::string basis_file_name() const;

	/// Returns a standard filename, used for writing out
	/// scalar product matrices.
	/// If leading is true, "_l" is appended to the filename.
	/// If "poly" is true "P_" is added to the filename, and if it is
	/// false "d_", as in double, is added to the filename.
	std::string spm_file_name( const bool leading,  const bool poly) const;

	/// Function for reading in the basis from the file filename.
	/// The basis should be in human readable format, of form:<br>
	/// 0  [{1,3,4,2}]<br>
	/// 1  [{1,4,3,2}]<br>
	/// 2  [{1,2}(3,4)]<br>
	/// i.e. first basis vector number 0,1,2..., then
	/// the Col_amp corresponding to the basis vector in question.
	/// The Col_amps may consist of several Col_strs, for example<br>
	/// 0  [(1,2,3,4)]+[(1,4,3,2)]<br>
	/// 1  [(1,2,4,3)]+[(1,3,4,2)]<br>
	///	2  [(1,3,4,2)]+[(1,2,4,3)]<br>
	/// and each Col_str may also contain a Polynomial.
	/// (The Polynomial should multiply the whole col_str,
	/// rather than a quark_line inside the []-brackets.)
	virtual void read_in_Col_basis( std::string filename );

	/// Function for reading in the basis from default filename
	/// (see basis_file_name).
	virtual void read_in_Col_basis();

	/// Read in a numerical matrix from a file filename
	/// (see spm_file_name) and save it as a double matrix, dmatr,
	/// in the member variable d_spm.
	/// The file should be in the format<br>
	/// {{d11,...,d1n},<br>
	/// ...,<br>
	/// {dn1,....dnn}},<br>
	/// and may contain comment lines starting with # at the top.
	void read_in_d_spm( std::string filename );

	/// Read in a numerical matrix from a file with default filename
	/// (see spm_file_name) and save it as a double matrix, dmatr,
	/// in the member variable d_spm.
	/// The file should be in the format
	/// {{d11,...,d1n},<br>
	/// ...,<br>
	/// {dn1,....dnn}},<br>
	/// and may contain comment lines starting with # at the top.
	void read_in_d_spm( );

	/// Read in a numerical matrix from the file filename and save it as a double matrix, dmatr,
	/// in the member variable leading_d_spm.
	/// The file should be in the format<br>
	/// {{d11,0,...,0},<br>
	/// ...,<br>
	/// {0,0,....dnn}},<br>
	/// and may contain comment lines starting with # at the top.
	void read_in_leading_d_spm( std::string filename );

	/// Read in a numerical matrix from a file with default filename (see spm_file_name)
	/// and save it as a double matrix, dmatr, in the member variable leading_d_spm.
	/// The file should be in the format<br>
	/// {{d11,0,...,0},<br>
	/// ...,<br>
	/// {0,0,....dnn}},<br>
	/// and may contain comment lines starting with # at the top.
	void read_in_leading_d_spm( );

	/// Read in a Polynomial matrix from the file filename and save it as a Poly_matr
	/// in the member variable P_spm.
	/// The file should be in the format<br>
	/// {{Poly11,...,Poly1n}, <br>
	/// ...,<br>
	/// {Polyn1,...,Polynn}},<br>
	/// and may contain comment lines starting with # at the top.
	void read_in_P_spm( std::string filename );

	/// Read in a Polynomial matrix from a file with default filename (see spm_file_name)
	/// and save it as a Poly_matr in the member variable P_spm.
	/// The file should be in the format<br>
	/// {{Poly11,...,Poly1n},<br>
	/// ...,<br>
	/// {Polyn1,...,Polynn}},<br>
	/// and may contain comment lines starting with # at the top.
	void read_in_P_spm( );

	/// Read in a Polynomial matrix from a file with default filename
	/// (see spm_file_name)  and save it as a Poly_matr
	/// in the member variable leading_P_spm.
	/// The file should be in the format<br>
	/// {{Poly11,0...,0},<br>
	/// ...,<br>
	/// {0,...,Polynn}},<br>
	/// and may contain comment lines starting with # at the top.
	void read_in_leading_P_spm( std::string filename );

	/// Reads in a Polynomial matrix from default filename (see spm_file_name)
	/// and save it as a Poly_matr in the member variable leading_P_spm.
	/// The file should be in the format<br>
	/// {{Poly11,0...,0},<br>
	/// ...,<br>
	/// {0,...,Polynn}},<br>
	/// and may contain comment lines starting with # at the top.
	void read_in_leading_P_spm(  );


	// Functions for writing

	/// Function for writing out the basis to a file with name filename.
	virtual void write_out_Col_basis( std::string filename) const;

	/// Function for writing out the basis to a file
	/// with default name (see basis_file_name).
	virtual void write_out_Col_basis() const;

	/// Writes out d_spm to the file filename.
	void write_out_d_spm( std::string filename ) const;

	/// Writes out d_spm to the standard filename, see spm_file_name.
	void write_out_d_spm( ) const;

	/// Writes out P_spm to the file filename.
	void write_out_P_spm( std::string filename ) const;

	/// Writes out P_spm to the standard filename, see spm_file_name.
	void write_out_P_spm( ) const;

	/// Writes out leading_d_spm to the file filename.
	void write_out_leading_d_spm( std::string filename ) const;

	/// Writes out leading_d_spm to the standard filename, see spm_file_name.
	void write_out_leading_d_spm( ) const;

	/// Writes out leading_P_spm to the file filename.
	void write_out_leading_P_spm( std::string filename ) const;

	/// Writes out leading_P_spm to the standard filename, see spm_file_name.
	void write_out_leading_P_spm( ) const;

	/******************** Other functions **********************/

	/// Simplifies all the basis vectors by using simplify on the
	/// individual Col_amps in the basis.
	void simplify();

	/// Each type of color basis has to implement a function for decomposing
	/// an amplitude in the color basis.
	virtual Poly_vec decompose( const Col_amp & Ca );

	/// Function for finding the resulting Col_amp after exchanging
	/// a gluon between parton p1 and parton p2 in the
	/// basis vector vec.
	Col_amp exchange_gluon( uint vec, int p1, int p2 );

	/// Function for calculating the color structure part of the soft anomalous
	/// dimension matrix. First calculates the effect of gluon exchange on a
	/// basis vector and then decomposes the result into the basis.
	/// For this to work the basis must clearly contain all resulting
	/// basis vectors, meaning for example that it can not be
	/// used for Tree_level_gluon_basis.
	/// The function is only available for the Trace_basis
	/// and the Orthogonal_basis classes.
	/// The ij-component of the resulting matrix gives the amplitude
	/// for ending up in component i if starting in component j.
	Poly_matr color_gamma( int p1, int p2 );

	/// Returns the number of quarks in the Col_basis after
	/// checking that each Col_str in each Col_amp
	/// has the same number of quarks.
	int n_quark_check() const;

	/// Returns the number of gluons in the Col_basis after
	/// checking that each Col_str in each Col_amp
	/// has the same number of gluons.
	int n_gluon_check() const;

	/// A function to rename the indices in two Col_strs, such that in the first
	/// they are called 1,2,3..., and in the second the relative order is kept.
	void rename_indices( Col_str & Cs1, Col_str & Cs2 ) const;

	/******************** Internal functions **********************/

protected:

	bool trace_basis;
	bool tree_level_gluon_basis;
	bool orthogonal_basis;

	/// The underlying function for calculation of scalar product
	/// matrices. Calculation depends on the arguments:
	/// save_P_spm, save_d_spm and use_mem.
	/// If save_P_spm is true the Polynomial scalar product
	/// matrix is saved to P_spm.
	/// If save_d_spm is true, the numerical scalar product matrix
	/// is saved to d_spm.
	/// The leading_d_spm is only calculated and written out if save_P_spm is true.
	/// If use_mem is true memoization is used
	/// in order to calculate a color topology only once.
	virtual void scalar_product_matrix( bool save_P_spm, bool save_d_spm, bool use_mem );

	/// Calculates element i,j in scalar product matrix using the scalar product.
	virtual Polynomial ij_entry( const int i, const int j ) const;

	/// Checking that a numerical (double) matrix is diagonal.
	/// Used for the leading version of the scalar product matrix.
	/// Returns true if the matrix is diagonal and false otherwise.
	bool check_diagonal( const dmatr & matr ) const;

	/// Checking that a numerical (double) matrix (the scalar product matrix)
	/// is symmetric.
	/// Returns true if the matrix is symmetric and false otherwise.
	bool check_symmetry( const dmatr & matr ) const;

	/// Makes consistency checks on the scalar product matrix.
	void check_spm() const;

	/// Converts a text string (in a file) to a basis,
	/// used by the read_in_basis functions.
	void Col_basis_of_str( std::string str );

	/// Function for writing out the basis in a human readable
	/// format to an ostream.
	virtual std::ostream&  write_out_Col_basis_to_stream( std::ostream&  out ) const;

	/// The operator << for Col_basis operator must be able to access the
	/// write_out_Col_basis_to_stream operator.
	friend std::ostream& operator<<( std::ostream& out, const Col_basis & Cb );

};// end class Col_basis


/// Define the operator << for col_basis.
std::ostream& operator<<( std::ostream& out, const col_basis & cb );


/// Define the operator << for Col_basis.
std::ostream& operator<<( std::ostream& out, const Col_basis & Cb );

}// end namespace ColorFull


#endif /* COLBASIS_H_ */

