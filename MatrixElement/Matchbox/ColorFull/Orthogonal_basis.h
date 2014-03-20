/*
 * Orthogonal_basis.h
 * Contains the declarations of the class Orthogonal_basis, related types and operators.
 * Created on: May 25, 2013
 * Author: Malin Sjodahl
 */

#ifndef COLORFULL_Orthogonal_basis_h
#define COLORFULL_Orthogonal_basis_h


#include "Col_basis.h"


namespace ColorFull {

/// This class is for containing orthogonal bases,
/// i.e. bases where the scalar product matrix is diagonal.
/// ColorFull has (currently) no functionality for creating
/// orthogonal bases, but orthogonal bases
/// can be read in using read_in _Col_basis.
class Orthogonal_basis:public Col_basis {
public:

	/// Default constructor.
	Orthogonal_basis():Col_basis(){
		orthogonal_basis = true;
	}

	/// To contain information about scalar products as a dvec,
	/// entry i is the square of vector i.
	dvec diagonal_d_spm;

	/// To contain information about scalar products as a Poly_vec,
	/// i.e., entry i is the square of vector i.
	Poly_vec diagonal_P_spm;


	/******************** Functions for scalar products **********************/

	/// Calculates the scalar product matrix assuming the basis to be orthogonal.
	/// Calculates both the double (d_spm) and the Polynomial (P_spm) matrices and
	/// saves to default file names.
	void scalar_product_matrix();

	/// Calculates the diagonal entries in the scalar product matrix,
	/// and (depending on arguments), saves them to the member variables
	/// diagonal_P_spm and diagonal_d_spm.
	/// This function is used by the Orthogonal_basis version of
	/// scalar_product_matrix.
	void diagonal_scalar_product_matrix( bool save_P_diagonal_spm, bool save_d_diagonal_spm, bool use_mem );

	/// The decomposition of a Col_amp in an orthogonal basis is done
	/// by calculating scalar products and dividing out the norm.
	/// The norm is evaluated numerically.
	Poly_vec decompose( const Col_amp & Ca );

	/// Function for calculating scalar products
	/// given the information about the basis and the scalar product matrix in the basis.
	/// The Col_amps are first decomposed using decompose,
	/// and then squared using the scalar product matrix P_spm.
	/// An orthogonal scalar product matrix is assumed.
	Polynomial scalar_product( const Col_amp & Ca1, const Col_amp & Ca2 );

	/// Function for calculating scalar products
	/// given the information about the basis
	/// and the scalar product matrix in numerical form.
	/// The Col_amps are first decomposed using decompose,
	/// and then squared using diagonal_d_spm.
	/// For Orthogonal_basis, an orthogonal scalar product matrix is assumed.
	cnum scalar_product_num( const Col_amp & Ca1, const Col_amp & Ca2 );

	/// Calculates the scalar product between decomposed amplitudes v1, V2
	/// using the diagonal_d_spm diagonal numerical scalar product matrix.
	/// The vectors needs to be expressed in the basis contained in cb,
	/// i.e., the decomposition has to be known.
	cnum scalar_product_num( const cvec & v1, const cvec & v2 );


	/******************** Functions for reading and writing **********************/

	/// Creates a default filename for writing out diagonal scalar products.
	/// The boolean variable leading should be true if the name is for a leading
	/// Nc variable. The filename is then modified accordingly.
	std::string diagonal_spm_file_name( const bool leading, const bool poly ) const;

	/// Writes out diagonal_d_spm to file filename.
	void write_out_diagonal_d_spm( std::string filename ) const;

	/// Writes out diagonal_d_spm to the standard filename, see diagonal_spm_file_name.
	void write_out_diagonal_d_spm( ) const;

	/// Writes out diagonal_P_spm to file filename.
	void write_out_diagonal_P_spm( std::string filename ) const;

	/// Writes out diagonal_P_spm to the standard filename, see diagonal_spm_file_name.
	void write_out_diagonal_P_spm( ) const;

	/// Function for writing out a dvec (the diagonal scalar products, diagonal_d_spm)
	/// to a file with standard filename given by diagonal_spm_file_name.
	/// The boolean variable leading should be true if dv only has leading
	/// Nc contributions. The filename is then modified accordingly.
	void write_out_diagonal_spm( const dvec & dv, const bool leading ) const;

	/// Function for writing out a Poly_vec (the diagonal scalar products, diagonal_d_spm)
	/// to a file with standard filename given by diagonal_spm_file_name.
	/// The boolean variable leading should be true if Poly_vec only has leading
	/// Nc contributions. The filename is then modified accordingly.
	void write_out_diagonal_spm( const Poly_vec & pv, const bool leading ) const;

private:

	/// Function for calculating the scalar products matrix.
	/// This function uses diagonal_scalar_product_matrix
	/// and saved the value of the diagonal scalar products between basis vector
	/// i and basis vector j in the i,j -entry in
	/// P_spm (if save P_spm is true) and d_spm (if save_d_spm is true).
	/// If use_mem is true, memoization is used.
	void scalar_product_matrix( bool save_P_spm, bool save_d_spm, bool use_mem );


};// end class Orthogonal_basis

}// end namespace ColorFull


#endif /* COLORFULL_Orthogonal_basis_h */
