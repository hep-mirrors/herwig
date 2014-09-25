/*
 * Col_functions.cc
 *	Contains definitions of functions and operators used for treating the color structure
 *  Created on: Jul 8, 2010
 *      Author: malin
 */

#include "Col_functions.h"
#include "parameters.h"
#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>

namespace ColorFull {

using boost::shared_ptr;

int  Col_functions::leading_Nc_pow( const Polynomial & Poly ) const{
	if( Poly.poly.empty()) {
		if( double_num(Poly)!=0 ) return 0;
		else return std::numeric_limits<int>::min();
	}
	int leading_pow=Poly.at(0).pow_Nc + Poly.at(0).pow_CF;

	// Loop over non-zero Monomials
	for (uint i=0; i<Poly.poly.size(); i++){
		if( Poly.at(i).pow_Nc + Poly.at(i).pow_CF > leading_pow && Poly.at(i).int_part!=0 )
			leading_pow=Poly.at(i).pow_Nc + Poly.at(i).pow_CF;
	}
	return leading_pow;
}

int  Col_functions::leading_Nc_pow( const Poly_vec & Pv ) const{

	int leading_pow=leading_Nc_pow( Pv.at(0) );

	// Loop over Polynomials
	for (uint p=0; p<Pv.size(); p++){
		if ( leading_Nc_pow( Pv.at(p) ) > leading_pow )
			leading_pow=leading_Nc_pow( Pv.at(p) );
	}

	return leading_pow;
}

int  Col_functions::leading_Nc_pow( std::vector< boost::shared_ptr<Polynomial> > Pvp) const{
	int leading_pow=leading_Nc_pow( *(Pvp.at(0)) );

	// Loop over Polynomials
	for (uint p=1; p<Pvp.size(); p++){
		if ( leading_Nc_pow( *(Pvp.at(p)) ) > leading_pow )
			leading_pow=leading_Nc_pow( *(Pvp.at(p)) );
	}

	return leading_pow;
}

Polynomial Col_functions::leading( const Polynomial & Poly ) const{


	// If the Polynomial is empty=1, or it it has only one term,
	// just return the Poly itself
	if( Poly.empty() ) return Poly;

	// Candidate for highest power of Nc+CF
	int cand_pow = Poly.at(0).pow_Nc + Poly.at(0).pow_CF;

	// For containing the leading Nc Polynomial
	Polynomial Leading_pol;
	Leading_pol.push_back(Poly.at(0));

	// If only one term, keep that term
	if( Poly.size()==1 ) Leading_pol=Poly;

	// If more than one term, look for higher powers
	// Looks for Monomials with highest Nc power
	if( Poly.size()>1 ){// found bug for Polynomials with one term 13 07 12
		for (int k = 1; k < Poly.size(); k++) {

			// If power higher
			if ((Poly.at(k).pow_Nc + Poly.at(k).pow_CF) > cand_pow) {
				cand_pow = Poly.at(k).pow_Nc + Poly.at(k).pow_CF;
				Leading_pol.clear();
				Leading_pol = Leading_pol * 0;
				Leading_pol = Leading_pol + Poly.at(k);
			}
			// If power equal
			else if (Poly.at(k).pow_Nc + Poly.at(k).pow_CF == cand_pow) {
				Leading_pol = Leading_pol + Poly.at(k);
			}
		}
	}
	// now Leading_pol contains terms with maximal power of CF plus Nc

	// If the variable full_CF is false (default), CF should be replaced with TR*Nc.
	// If numerical evaluation is done after this, CF will thus be evaluated to TR Nc.
	// If full_cf is true, CF will be replaced by  CF(Nc).
	if ( ! full_CF ) {
		// Loop over terms and replace CF with the Nc -> infinity limit of CF =TR*Nc
		for (int i = 0; i < Leading_pol.size(); i++) {
			int pow_CF = Leading_pol.at(i).pow_CF;
			Leading_pol.at(i).pow_Nc += pow_CF;
			Leading_pol.at(i).pow_TR += pow_CF;
			Leading_pol.at(i).pow_CF = 0;
		}
	}

	Leading_pol.simplify();
	return Leading_pol;
}

Poly_vec Col_functions::leading( const Poly_vec & Pv ) const{

	// To contain the result
	Poly_vec Pv_res;

	// Take the leading part of each Polynomial
	// after this all Monomials in the SAME Polynomial have the same Nc power
	for (uint i = 0; i < Pv.size(); i++) {
		// Take the leading terms in each component
		Pv_res.push_back( leading( Pv.at(i) ) );
	}

	// Find the leading power
	int leading_pow = leading_Nc_pow(Pv_res);

	// Loop over entries and keep only those with maximal power
	for (uint p = 0; p < Pv_res.size(); p++) {
		if (leading_Nc_pow(Pv_res.at(p)) != leading_pow) {
			Pv_res.at(p) = 0 * Pv_res.at(p);
			Pv_res.at(p).simplify();
		}
	}

	return Pv_res;
}

Poly_matr Col_functions::leading( const Poly_matr & Pm )  const{

	// To contain the result
	Poly_matr Pm_res;
	Pm_res.pm.reserve(Pm.pm.size());
	// Loop over Poly_vecs, and take the leading part of each Poly_vec
	for (uint i = 0; i < Pm.size(); i++) {
		// Take the leading terms in each Poly_vec
		Poly_vec Pvl= leading( Pm.at(i) );
		Pm_res.push_back( Pvl );
	}
	int cand_pow = leading_Nc_pow( Pm.at(0).pv );

	// Loop over Poly_vecs to locate highest power
	for (uint k = 0; k < Pm_res.size(); k++) {
	  int pow_at_k = leading_Nc_pow(Pm_res.at(k).pv);
		if (pow_at_k > cand_pow)
			cand_pow = pow_at_k;
	}

	// In all Polynomials, in all Poly_vec's,
	// put all elements which doesn't have highest power to 0
	for (uint k = 0; k < Pm_res.size(); k++) {
		// In each Poly_vec, loop over all Polynomials
		for (uint p = 0; p < Pm_res.at(k).size(); p++) {
			if (leading_Nc_pow(Pm_res.at(k).at(p)) != cand_pow) {
				//Pm_res.at(k).pv.at(p) = Pm_res.at(k).pv.at(p) * 0;
				Pm_res.at( k,p )=Pm_res.at(k).at(p) * 0;
				Pm_res.at(k,p).simplify();
			}
		}
	}

	return Pm_res;
}


Poly_vec Col_functions::leading( std::vector<boost::shared_ptr<Polynomial> >  Pvp )  const{
	Poly_vec Pv_res;

	// Loop over entries in vector (Poly_vecs), and take the leading part of each Polynomial
	// after this each Monomial has the same Nc+CF-power
	for (uint i = 0; i < Pvp.size(); i++) {
		// Take the leading terms in each component
		boost::shared_ptr<Polynomial> the_pointer=Pvp.at(i);
		Pv_res.push_back( leading( (*the_pointer) ) );
	}

	// Find the leading power
	int leading_pow = leading_Nc_pow(Pv_res);

	// Loop over entries and keep only those with maximal power
	for (uint p = 0; p < Pv_res.size(); p++) {
		if (leading_Nc_pow(Pv_res.at(p)) != leading_pow) {
			Pv_res.at(p) = 0 * Pv_res.at(p);
			Pv_res.at(p).simplify();
		}
	}

	return Pv_res;
}


dmatr Col_functions::leading( std::vector< std::vector< boost::shared_ptr<Polynomial> > >  Ppm ) const {

	// To contain the result as a matrix of double
	dmatr res;
	res.reserve( Ppm.size() );

	int cand_pow = leading_Nc_pow( Ppm.at(0) );

	// Loop over Poly_vecs to locate highest power
	for (uint k = 0; k < Ppm.size(); k++) {
		int the_pow=leading_Nc_pow(Ppm.at(k));
		if ( the_pow > cand_pow )
			cand_pow = the_pow;
	}

	// In all Polynomials, in all Poly_vec's,
	// put all elements which doesn't have highest power to 0
	for (uint k = 0; k < Ppm.size(); k++) {
		dvec dummy;
		res.push_back(dummy); // Not to run out of range
		// In each Poly_vec, loop over all Polynomials
		for (uint p = 0; p < Ppm.at(k).size(); p++) {
			Polynomial the_poly=*Ppm.at(k).at(p);
			leading_Nc_pow( the_poly );
			if ( leading_Nc_pow( the_poly ) != cand_pow ) {res.at(k).push_back(0);}
			else res.at(k).push_back( double_num(leading(the_poly)) );
		}
	}

	return res;
}

std::map< std::string, Polynomial > Col_functions::leading( std::map< std::string, Polynomial > mem_map  )  const{

	// To contain (string, leading(Polynomial))
	std::map< std::string, Polynomial > res;

	// Iterator for looping
	std::map< std::string, Polynomial >::iterator iter=mem_map.begin();

    // Find the highest power of Nc+CF
    int pow_cand=leading_Nc_pow( (iter->second) );
    for( iter =  mem_map.begin(); iter !=  mem_map.end(); ++iter ) {
    	Polynomial Poly= (iter->second);
    	if( leading_Nc_pow(Poly) > pow_cand ) pow_cand=leading_Nc_pow(Poly);
    }

    for( iter =  mem_map.begin(); iter !=  mem_map.end(); ++iter ) {

    	// Insert pair of string and the leading versions of the Polynomial
    	Polynomial lead_Poly= leading( (iter->second) );
    	if ( leading_Nc_pow(lead_Poly) != pow_cand ) lead_Poly=lead_Poly*0;

    	res.insert( make_pair( iter->first, lead_Poly ) );
    } // After this only the leading part of each Polynomial contributes

	return res;
}

cnum Col_functions::cnum_num( const Monomial & Mon ) const {

	cnum res=Mon.cnum_part*static_cast<double>(Mon.int_part);
	for ( int ix=0; ix<Mon.pow_Nc; ix++ )
		res *= Nc;
	for ( int ix=0; ix<Mon.pow_CF; ix++ )
		res *= CF;
	for ( int ix=0; ix<Mon.pow_TR; ix++ )
		res *= TR;
	for ( int ix=0; ix>Mon.pow_Nc; ix-- )
		res /= static_cast<double>(Nc);
	for ( int ix=0; ix>Mon.pow_CF; ix-- )
		res /= static_cast<double>(CF);
	for ( int ix=0; ix>Mon.pow_TR; ix-- )
		res /= static_cast<double>(TR);
	//cnum res=pow(Nc, Mon.pow_Nc)*pow(CF, Mon.pow_CF)*pow(TR, Mon.pow_TR)*Mon.int_part* Mon.cnum_part;

	return res;
}

double  Col_functions::double_num( const Monomial & Mon ) const{

	double im=imag( cnum_num(Mon) );
	double re=real( cnum_num(Mon) );

	// Warn if the complex number has significant imaginary parts
	if( im!=0.0 and re!=0.0 ) {
		double ratio =im/re;
		if( ratio > accuracy )
			std::cerr << "Col_functions::double_num(Mon): Warning keeping only real part of complex number, the ratio im/re was " << ratio << std::endl;
	}
	else if ( im!=0.0 and re==0.0 ){
					std::cerr << "Col_functions::double_num(Mon): Warning keeping only real part of complex number, imaginary part was " << im << std::endl;
	}

	return re;
}

cnum Col_functions::cnum_num( const Polynomial & Poly ) const {

	// An empty Polynomial has numerical value 1
	if( Poly.empty() ) return 1.0;

	cnum res = 0;

	// Add contributions from Monomials
	for ( int k = 0; k < Poly.size(); k++ ) {
		cnum part_k;
		part_k=cnum_num( Poly.at(k) );
		res = res + part_k;
	}

	return res;
}

double Col_functions::double_num( const Polynomial & Poly ) const{
	double im=imag( cnum_num(Poly) );
	double re=real( cnum_num(Poly) );

	// Warn if the complex number has significant imaginary parts
	if( im!=0.0 and re!=0.0 ) {
		double ratio =im/re;
		if( ratio > accuracy )
			std::cerr << "Col_functions::double_num: Warning keeping only real part of complex number, the ratio im/re was " << ratio << std::endl;
	}
	else if ( im!=0.0 and re==0.0 ){
		std::cerr << "Col_functions::double_num: Warning keeping only real part of complex number, imaginary part was " << im << std::endl;
	}

	return re;
}


Col_amp Col_functions::emit_gluon( const Col_str & in_Col_str, int emitter, int g_new ) const{

  // Locate the emitter in the Col_str
  std::pair<int, int> place=in_Col_str.find_parton(emitter);

  // Find what kind, q qbar or g, the emitter is
  std::string kind=in_Col_str.find_kind(emitter);

  // Defining Col_str to return (two needed in case of g)
  Col_str out_Col_str1=in_Col_str;
  Col_str out_Col_str2=in_Col_str;
  // Defining Col_amp to return
  Col_amp out_Col_amp;

  // If the emitter is a quark
  if( kind =="q" ){
    // Add new parton index at second first place in relevant quark_line
    out_Col_str1.insert( place.first, place.second+1, g_new );
    out_Col_amp.ca.push_back( out_Col_str1 );
  }
  // If the emitter is an anti-quark
  else if( kind =="qbar" ){
    // Add new gluon before the qbar
    out_Col_str1.insert( place.first, place.second, g_new );
    // Change sign
    out_Col_str1.Poly=out_Col_str1.Poly*(-1);
    // Making a Col_amp to return
    out_Col_amp.ca.push_back( out_Col_str1 );
  }
  // If the emitter is a g
  else if( kind =="g" ){
    // If the emitter is a gluon a new gluon line should be inserted in two
    // different places, before and after the emitter

    // -sign when inserting the gluon before
    // Change sign Poly for first term
    Monomial Mon_tmp;
    Mon_tmp.int_part=-1;
    out_Col_str1.Poly=out_Col_str1.Poly*Mon_tmp;

    // Inserting the gluon before
    out_Col_str1.insert( place.first, place.second, g_new );
    // Inserting the gluon after
    out_Col_str2.insert( place.first, place.second+1, g_new );
    // Appending result to out_Col_amp
    out_Col_amp.ca.push_back( out_Col_str1 );
    out_Col_amp.ca.push_back( out_Col_str2 );
    // Normal order
    out_Col_amp.normal_order();
  }
  out_Col_amp.simplify();

  return out_Col_amp;
}


Col_amp Col_functions::emit_gluon( const Col_amp & Ca_in, int emitter, int g_new ) const{

	Col_amp Ca_out;

	// Emit from each Col_str, and append to new Col_amp
	for( uint m=0; m< Ca_in.ca.size(); m++ ){
		Col_amp part_m=emit_gluon(Ca_in.ca.at(m), emitter, g_new);
		Ca_out.append(part_m.ca);
	}

	return Ca_out;
}


// See my general color structure paper
Col_amp Col_functions::exchange_gluon( const Col_str & Cs, int p1, int p2 )  const{

	Col_str Cs_copy=Cs;

	// Find out kind and location of partons
	std::string kind1 = Cs_copy.find_kind(p1);
	std::string kind2 = Cs_copy.find_kind(p2);
	std::pair<int, int> place1 = Cs_copy.find_parton(p1);
	std::pair<int, int> place2 = Cs_copy.find_parton(p2);



	// To contain result
	Col_amp Ca;

	// Make sure the col_f in multiplying participating Ql's are 0
	Polynomial Poly1; // For comparing Poly is 1, all powers are 0
	if (!(Poly1 == Cs_copy.cs.at(place1.first).Poly)) {
		// Move factor of Quark_line to Col_str
		Cs_copy.Poly = Cs_copy.Poly * Cs_copy.cs.at(place1.first).Poly;
		Cs_copy.cs.at(place1.first).Poly = Poly1;
	}
	if (!(Cs_copy.cs.at(place2.first).Poly == Poly1)) {
		// Move factor of Quark_line to Col_str
		Cs_copy.Poly = Cs_copy.Poly * Cs_copy.cs.at(place2.first).Poly;
		Cs_copy.cs.at(place2.first).Poly = Poly1;
	}

	//	sign
	if ( place1.first == place2.first && place1.second == place2.second ) {

		Monomial mon;
		// For qq or qbar qbar
		// If the qs are the same, the result is just CF*old for q and qbar
		// there are 0 minus signs for quarks and 2 for anti-quarks
		if (kind1 == "q" or kind1 == "qbar")
			mon.pow_CF = 1;
		// for gluons
		else{
			mon.pow_Nc = 1;
			// There is a factor TR from the contracted gluon,
			// The sign can be checked by writing out four terms
			// corresponding to inserting the new gluon
			// first/first, first/last, last/first, last/last
			mon.pow_TR+=1;
			mon.int_part*=2;
		}
		Cs_copy.Poly = Cs_copy.Poly * mon;
		Ca.ca.push_back( Cs_copy );
		return Ca;
	}

	// Exchange between qq or qbar qbar
	if ((kind1 == "q" && kind2 == "q") or (kind1 == "qbar" && kind2 == "qbar")) {

		Col_str Cs2 = Cs_copy;

		// The q's normally sit on different quark_lines

		// The leading Nc-tems has the q's exchanged
		Cs_copy.cs.at(place1.first).ql.at(place1.second) = p2;
		Cs_copy.cs.at(place2.first).ql.at(place2.second) = p1;
		// Multiply with TR
		Monomial Mon;
		Mon.pow_TR = 1;

		Cs_copy.Poly = Cs_copy.Poly * Mon;

		// The suppressed term is a copy of the old term,
		// but multiplied with -TR /Nc
		Mon.pow_TR = 1;
		Mon.int_part = -1;
		Mon.pow_Nc = -1;
		Cs2.Poly = Cs2.Poly * Mon;

		Ca.ca.push_back(Cs_copy);
		Ca.ca.push_back(Cs2);
	}

	// Exchange between q qbar
	if ((kind1 == "q" && kind2 == "qbar") or (kind1 == "qbar" && kind2 == "q")) {

	  std::pair<int, int> q_place;
	  std::pair<int, int> qbar_place;
		int the_q=-1;
		int the_qbar=-1;

		if ((kind1 == "q" && kind2 == "qbar")) {
			q_place = place1;
			qbar_place = place2;
			the_q = p1;
			the_qbar = p2;
		} else if ((kind1 == "qbar" && kind2 == "q")) {
			q_place = place2;
			qbar_place = place1;
			the_q = p2;
			the_qbar = p1;
		}

		Col_str Cs2 = Cs_copy;

		// If q and qbar are not part of same ql
		if (q_place.first != qbar_place.first) {
			// The non-suppressed term, obtained by erasing q qbar and joining rest
			Cs_copy.erase(q_place);
			Cs_copy.erase(qbar_place);

			// In the ql of the qbar, append the ql of the q, then erase the ql of the q
			Cs_copy.cs.at(qbar_place.first).append(Cs_copy.cs.at(q_place.first).ql);
			// and erase the ql of the q
			Cs_copy.erase(q_place.first);
		}
		// If q and qbar sit on same ql
		else {
			// The non-suppressed term, obtained by erasing q qbar and joining rest to a ring
			Cs_copy.erase(qbar_place);
			Cs_copy.erase(q_place);
			Cs_copy.cs.at(q_place.first).open = false;
			// If only one quark remains the term should be 0
		}

		// q and qbar must be found
		assert( the_q!=-1);
		assert( the_qbar!=-1);

		// Construct Quark_line={q,qbar}, and append it
		Quark_line Ql;
		Ql.open = true;
		Ql.push_back( the_q );
		Ql.push_back( the_qbar );
		Cs_copy.push_back(Ql);

		// Multiply with TR
		Monomial Mon;
		Mon.pow_TR = 1;
		Cs_copy.Poly = Cs_copy.Poly * Mon;

		// The suppressed term is a copy of the old term,
		// but multiplied with -TR / Nc
		Mon.pow_TR=1;
		Mon.int_part = -1;
		Mon.pow_Nc = -1;
		Cs2.Poly = Cs2.Poly * Mon;

		// Construct  the resulting color_amplitude and return
		// Make sure there is not only one or 0 gluons in a closed ql
		// Append Cs to ca if >=2 gluons
		// if open, append
		if (Cs_copy.cs.at(q_place.first).open)
			Ca.ca.push_back(Cs_copy);
		// if closed but >= 2 gluons, append
		else if (Cs_copy.cs.at(q_place.first).ql.size() >= 2 && !Cs_copy.cs.at(
				q_place.first).open)
			Ca.ca.push_back(Cs_copy);
		// If 1 gluon, the term is 0 and is not appended

		// If 0 gluons the empty colsed ql is equivalent to multiplying with Nc, obtained after simplification
		else if (Cs_copy.cs.at(q_place.first).empty()
				&& !Cs_copy.cs.at(q_place.first).open) {
			// return 2nd term * (1-Nc^2)
			Ca.ca.push_back(Cs_copy);
			Ca.ca.push_back(Cs2);
			Ca.simplify();
			return Ca;
		}
		Ca.ca.push_back(Cs2);
	}

	// Exchange between q g
	if ((kind1 == "q" && kind2 == "g") or (kind1 == "g" && kind2 == "q")) {

	  std::pair<int, int> q_place;
	  std::pair<int, int> g_place;
		int the_q=0;
		int the_g=0;

		if (kind1 == "q" && kind2 == "g") {
			q_place = place1;
			g_place = place2;
			the_q = p1;
			the_g = p2;
		} else if (kind1 == "g" && kind2 == "q") {
			q_place = place2;
			g_place = place1;
			the_q = p2;
			the_g = p1;
		}

		// For containing Quark_lines
		Quark_line Ql11;
		Quark_line Ql12;
		Quark_line Ql21;
		Quark_line Ql22;

		Col_str Cs1 = Cs_copy;
		Col_str Cs2 = Cs_copy;

		// If the q and the g are not part of same ql
		if (g_place.first != q_place.first) {

			// Take split the Ql with the gluon after the gluon
			Quark_line Qlg_start = Cs_copy.cs.at(g_place.first).before( g_place.second );

			Quark_line Qlg_end = Cs_copy.cs.at(g_place.first).after( g_place.second );
			Quark_line Qlq_end = Cs_copy.cs.at(q_place.first).after( q_place.second );

			// First part, + sign
			// First part, first Ql
			Ql11 = Qlg_end;
			Ql11.prepend(the_g);
			Ql11.prepend(the_q);

			// First part, second Ql
			Ql12 = Qlg_start;
			Ql12.append(Qlq_end.ql);

			// If the ql of the gluon was closed account for this by gluing qls together
			if (!Qlg_start.open) {
				Ql11.append(Ql12.ql);
				Ql11.open = true;
			};

			// Second part, first Ql
			Ql21 = Qlg_end;
			Ql21.prepend(the_q);

			// Second part, second Ql
			Ql22 = Qlg_start;
			Ql22.append(the_g);
			Ql22.append(Qlq_end.ql);

			// If the ql of the gluon was closed account for this by gluing qls together
			if (!Qlg_start.open) {
				Ql21.append(Ql22.ql);
				Ql21.open = true;
			};

			// Replace Ql's in the Col_str
			Cs1.cs.at(q_place.first) = Ql11;
			Cs2.cs.at(q_place.first) = Ql21;

			// There is only one ql if the ql of the gluon was closed
			if (Qlg_start.open) {
				Cs1.cs.at(g_place.first) = Ql12;
				Cs2.cs.at(g_place.first) = Ql22;
			} else { // Erase 2nd ql
				Cs1.erase(g_place.first);
				Cs2.erase(g_place.first);
			}

		}
		// If q and g are part of same ql
		else {

			Quark_line Ql1 = Cs_copy.cs.at(g_place.first).before( g_place.second );
			Ql1.ql.erase(Ql1.ql.begin());
			Quark_line Ql2 = Cs_copy.cs.at(g_place.first).after( g_place.second);

			// First term, First ql
			Ql11 = Ql1;
			Ql11.open = false;

			// First term, 2nd ql
			Ql12 = Ql2;
			Ql12.prepend(the_g);
			Ql12.prepend(the_q);

			// 2nd term, 1st ql
			Ql21 = Ql1;
			Ql21.append(the_g);
			Ql21.open = false;

			// 2nd term 2nd ql
			Ql22 = Ql2;
			Ql22.prepend(the_q);

			// Replace Ql's in the Col_str
			Cs1.cs.at(q_place.first) = Ql11;
			Cs1.cs.insert(Cs1.cs.begin() + q_place.first + 1, Ql12);

			// Replace Ql's in the Col_str
			Cs2.cs.at(q_place.first) = Ql21;
			Cs2.cs.insert(Cs2.cs.begin() + q_place.first + 1, Ql22);

		}

		// Multiply first part with TR
		Monomial Mon_tmp;
		Mon_tmp.pow_TR = 1;
		Cs1.Poly = Cs1.Poly * Mon_tmp;

		// Multiply second part with -TR
		Monomial Mon_tmp2;
		Mon_tmp2.pow_TR = 1;
		Mon_tmp2.int_part = -1;
		Cs2.Poly = Cs2.Poly * Mon_tmp2;


		// Add to result
		Ca.ca.push_back(Cs1);
		Ca.ca.push_back(Cs2);
	}// end exchange between q g


	// Exchange between qbar g (q and g in my paper)
	if ((kind1 == "qbar" && kind2 == "g") or (kind1 == "g" && kind2 == "qbar")) {

	  std::pair<int, int> qbar_place;
	  std::pair<int, int> g_place;
		int the_qbar=0;
		int the_g=0;

		if (kind1 == "qbar" && kind2 == "g") {
			qbar_place = place1;
			g_place = place2;
			the_qbar = p1;
			the_g = p2;
		} else if (kind1 == "g" && kind2 == "qbar") {
			qbar_place = place2;
			g_place = place1;
			the_qbar = p2;
			the_g = p1;
		}

		// For containing Quark_lines
		Quark_line Ql11, Ql12, Ql21, Ql22;

		// For containing resulting Col_str's
		Col_str Cs1 = Cs_copy;
		Col_str Cs2 = Cs_copy;

		// If ql and g on different quark_lines
		if (g_place.first != qbar_place.first) {
			// Split ql's into parts
			Quark_line Qlg_start = Cs_copy.cs.at(g_place.first).before( g_place.second) ;
			Quark_line Qlg_end = Cs_copy.cs.at(g_place.first).after( g_place.second );
			Quark_line Qlqbar_start =  Cs_copy.cs.at( qbar_place.first ).before(qbar_place.second);

			// First part, + sign, first Ql
			Ql11 = Qlqbar_start;
			Ql11.append(the_g);
			Ql11.append(Qlg_end.ql);

			// First part,  + sign, second Ql
			Ql12 = Qlg_start;
			Ql12.append(the_qbar);

			// If the ql of the gluon was closed account for this by gluing qls together
			if (!Qlg_start.open) {
				Ql11.append(Ql12.ql);
				Ql11.open = true;
			};

			// Second part, - sign, first Ql
			Ql21 = Qlqbar_start;
			Ql21.append(Qlg_end.ql);

			// Second part, - sign,  second Ql
			Ql22 = Qlg_start;
			Ql22.append(the_g);
			Ql22.append(the_qbar);

			// If the ql of the gluon was closed account for this by gluing qls together
			if (!Qlg_start.open) {
				Ql21.append(Ql22.ql);
				Ql21.open = true;
			};

			// Replace Ql's in the Col_str
			Cs1 = Cs_copy;
			Cs2 = Cs_copy;
			if (Qlg_start.open) {
				Cs1.cs.at(qbar_place.first) = Ql11;
				Cs1.cs.at(g_place.first) = Ql12;
				Cs2.cs.at(qbar_place.first) = Ql21;
				Cs2.cs.at(g_place.first) = Ql22;
			} else { // If the ql of the g was closed, there is only one ql, erase 2nd ql
				Cs1.cs.at(qbar_place.first) = Ql11;
				Cs1.erase(g_place.first);

				Cs2.cs.at(qbar_place.first) = Ql21;
				Cs2.erase(g_place.first);
			}
		}
		// If qbar and g are part of same quark_line
		else {

			Quark_line Ql1 = Cs_copy.cs.at(g_place.first).before(g_place.second );
			Quark_line Ql2 = Cs_copy.cs.at(g_place.first).after( g_place.second );
			Ql2.ql.erase(Ql2.ql.end() - 1);

			// First term, First ql
			Ql11 = Ql1;
			Ql11.append(the_qbar);

			// First term, 2nd ql
			Ql12 = Ql2;
			Ql12.prepend(the_g);
			Ql12.open = false;

			// 2nd term, 1st ql
			Ql21 = Ql1;
			Ql21.append(the_g);
			Ql21.append(the_qbar);

			// 2nd term 2nd ql
			Ql22 = Ql2;
			Ql22.open = false;

			// Replace Ql's in the Col_str
			Cs1 = Cs_copy;
			Cs1.cs.at(qbar_place.first) = Ql11;
			Cs1.cs.insert(Cs1.cs.begin() + qbar_place.first + 1, Ql12);

			// Replace Ql's in the Col_str
			Cs2 = Cs_copy;
			Cs2.cs.at(qbar_place.first) = Ql21;
			Cs2.cs.insert(Cs2.cs.begin() + qbar_place.first + 1, Ql22);
		}
		// Multiply with TR
		Monomial Mon_tmp;
		Mon_tmp.pow_TR = 1;
		Cs1.Poly = Cs1.Poly * Mon_tmp;

		// Multiply with -TR
		Monomial Mon_tmp2;
		Mon_tmp2.pow_TR = 1;
		Mon_tmp2.int_part = -1;
		Cs2.Poly = Cs2.Poly * Mon_tmp2;

		Ca.ca.push_back(Cs1);
		Ca.ca.push_back(Cs2);

	}

	// Exchange between g g
	if (kind1 == "g" && kind2 == "g") {

		// For containing resulting Col_str's
		Col_str Cs1 = Cs_copy;
		Col_str Cs2 = Cs_copy;
		Col_str Cs3 = Cs_copy;
		Col_str Cs4 = Cs_copy;

		// For containing Quark_lines
		Quark_line Ql11, Ql12, Ql21, Ql22, Ql31, Ql32, Ql41, Ql42;

		// If g's on different Ql's
		if (place1.first != place2.first) {

			// Split both ql's into parts
			Quark_line Ql1_start = Cs_copy.cs.at(place1.first).before( place1.second );
			Quark_line Ql1_end = Cs_copy.cs.at(place1.first).after( place1.second );
			Quark_line Ql2_start = Cs_copy.cs.at(place2.first).before( place2.second );
			Quark_line Ql2_end = Cs_copy.cs.at(place2.first).after( place2.second );

			// First term, both quarks on 2nd ql, First ql
			Ql11 = Ql1_start;
			Ql11.append(Ql2_end.ql);
			// First term, 2nd ql
			Ql12 = Ql2_start;
			Ql12.append(p2);
			Ql12.append(p1);
			Ql12.append(Ql1_end.ql);

			// 2nd term, g1 on ql1, g2 on ql2, 1st ql
			Ql21 = Ql1_start;
			Ql21.append(p1);
			Ql21.append(Ql2_end.ql);
			// 2nd term, 2nd ql
			Ql22 = Ql2_start;
			Ql22.append(p2);
			Ql22.append(Ql1_end.ql);

			// 3rd term, g2 on ql1, g1 on ql2,1st ql
			Ql31 = Ql1_start;
			Ql31.append(p2);
			Ql31.append(Ql2_end.ql);
			// 3rd term, 2nd ql
			Ql32 = Ql2_start;
			Ql32.append(p1);
			Ql32.append(Ql1_end.ql);

			// last term, both quarks on 1st ql,1st ql
			Ql41 = Ql1_start;
			Ql41.append(p1);
			Ql41.append(p2);
			Ql41.append(Ql2_end.ql);

			// last term, 2nd ql
			Ql42 = Ql2_start;
			Ql42.append(Ql1_end.ql);

			// Replace Ql's in the Col_str
			// If both Ql's are open
			if (Cs_copy.cs.at(place1.first).open && Cs_copy.cs.at(place2.first).open) {
				Cs1.cs.at(place1.first) = Ql11;
				Cs1.cs.at(place2.first) = Ql12;

				Cs2.cs.at(place1.first) = Ql21;
				Cs2.cs.at(place2.first) = Ql22;

				Cs3.cs.at(place1.first) = Ql31;
				Cs3.cs.at(place2.first) = Ql32;

				Cs4.cs.at(place1.first) = Ql41;
				Cs4.cs.at(place2.first) = Ql42;
			}
			// If string where closed glue together ends
			// If first ql closed
			else if (!Cs_copy.cs.at(place1.first).open) {
				Ql12.append(Ql11.ql);
				Ql22.append(Ql21.ql);
				Ql32.append(Ql31.ql);
				Ql42.append(Ql41.ql);

				//The result should be open/closed if the ql2 was open/closed
				Ql12.open = Cs_copy.cs.at(place2.first).open;
				Ql22.open = Cs_copy.cs.at(place2.first).open;
				Ql32.open = Cs_copy.cs.at(place2.first).open;
				Ql42.open = Cs_copy.cs.at(place2.first).open;

				//The reslut contaied int Qlx2 sould be used for the Color structure
				// Replace Ql's in the Col_str
				Cs1.cs.at(place2.first) = Ql12;
				Cs2.cs.at(place2.first) = Ql22;
				Cs3.cs.at(place2.first) = Ql32;
				Cs4.cs.at(place2.first) = Ql42;

				Cs1.erase(place1.first);
				Cs2.erase(place1.first);
				Cs3.erase(place1.first);
				Cs4.erase(place1.first);
			}
			// if only 2nd ql closed
			if (!Cs_copy.cs.at(place2.first).open && Cs_copy.cs.at(place1.first).open) {
				Ql11.append(Ql12.ql);
				Ql21.append(Ql22.ql);
				Ql31.append(Ql32.ql);
				Ql41.append(Ql42.ql);

				// The result should be open
				Ql11.open = true;
				Ql21.open = true;
				Ql21.open = true;
				Ql21.open = true;

				//The result contained int Qlx1 should be used for the Color structure
				// Replace Ql's in the Col_str
				Cs1.cs.at(place2.first) = Ql11;
				Cs2.cs.at(place2.first) = Ql21;
				Cs3.cs.at(place2.first) = Ql31;
				Cs4.cs.at(place2.first) = Ql41;

				Cs1.erase(place1.first);
				Cs2.erase(place1.first);
				Cs3.erase(place1.first);
				Cs4.erase(place1.first);

			}
		}
		// If g's on same ql's
		if (place1.first == place2.first) {

		  Quark_line Ql_start = Cs_copy.cs.at(place1.first).before( std::min(place2.second, place1.second));
			// the part between the g's
			Quark_line Ql_between = Cs_copy.cs.at(place1.first).before(std::max(
					place2.second, place1.second));
			Ql_between = Ql_between.after( std::min(place2.second, place1.second) );

			// the end after 2nd g
			Quark_line Ql_end = Cs_copy.cs.at(place1.first).after( std::max(place2.second,
					place1.second));

			// we have to know what g is first to put them back in right order
			int first_g, last_g;
			if( place1.second < place2.second ) {
				first_g=p1;
				last_g=p2;
			}
			else {
				first_g=p2;
				last_g=p1;
			}


			// First term, both gluons inbetween
			Ql11 = Ql_start;
			Ql11.append(Ql_end.ql);
			// First term, part between quarks
			Ql12 = Ql_between;
			Ql12.append( last_g );
			Ql12.append( first_g );
			Ql12.open = false;

			// 2nd term, g1 after first part, g2 after 2nd
			Ql21 = Ql_start;
			Ql21.append( first_g );
			Ql21.append(Ql_end.ql);
			// First term, part between quarks
			Ql22 = Ql_between;
			Ql22.append( last_g );
			Ql22.open = false;

			// 3rd term, g2 after first part, g1 after 2nd
			Ql31 = Ql_start;
			Ql31.append( last_g);
			Ql31.append(Ql_end.ql);
			// First term, part between quarks
			Ql32 = Ql_between;
			Ql32.append( first_g );
			Ql32.open = false;

			// 4th term, g2 after first part, g1 after 2nd
			Ql41 = Ql_start;
			Ql41.append( first_g );
			Ql41.append( last_g );
			Ql41.append(Ql_end.ql);
			// First term, part between quarks
			Ql42 = Ql_between;
			Ql42.open = false;

			Cs1.cs.at(place1.first) = Ql11;
			Cs1.cs.insert(Cs1.cs.begin() + place1.first + 1, Ql12);
			Cs2.cs.at(place1.first) = Ql21;
			Cs2.cs.insert(Cs2.cs.begin() + place1.first + 1, Ql22);
			Cs3.cs.at(place1.first) = Ql31;
			Cs3.cs.insert(Cs3.cs.begin() + place1.first + 1, Ql32);
			Cs4.cs.at(place1.first) = Ql41;
			Cs4.cs.insert(Cs4.cs.begin() + place1.first + 1, Ql42);

		}

		// Multiplying with TR where appropriate
		Monomial Mon_tmp;
		Mon_tmp.pow_TR = 1;

		Cs2.Poly = Cs2.Poly * Mon_tmp;
		Cs3.Poly = Cs3.Poly * Mon_tmp;

		// Multiplying with -TR where appropriate
		Mon_tmp.pow_TR = 1;

		Mon_tmp.int_part = -1;
		Cs1.Poly = Cs1.Poly * Mon_tmp;
		Cs4.Poly = Cs4.Poly * Mon_tmp;

		Ca.ca.push_back(Cs1);
		Ca.ca.push_back(Cs2);
		Ca.ca.push_back(Cs3);
		Ca.ca.push_back(Cs4);

	}

	// Simplify Col_amp
	Ca.simplify();

	return Ca;

}


Col_amp Col_functions::exchange_gluon( const  Col_amp & Ca, int p1, int p2 ) const{

  // Final col_amp to return
  Col_amp Ca_out;

  // Exchange in each Col_str, and append to new Col_amp
  for(uint m=0; m< Ca.ca.size(); m++ ){
    // Exchange in Col_str m
    Col_amp part_m=exchange_gluon( Ca.ca.at(m), p1, p2 );
    // Append result
    Ca_out.append(part_m.ca);
  }

  return Ca_out;
}

/*
Col_str Col_functions::contract_quarks( Col_str Cs1, Col_str Cs2 )  const{

	std::vector<int> q_place;
	std::vector<int> q_place2;

	// The conjugate of Cs1
	Col_str conj_Cs1 = Cs1;
	conj_Cs1.conjugate();

	// The total color structure
	Col_str Cs = conj_Cs1*Cs2;

	// Count how many quarks should be contracted
	int n_q = Cs.n_quark();

	// As long as there are quark_lines left to contract
	while (n_q > 0) {
		// Find first quark in Cs1 by looping over Quark_lines
		for (int i = 0; (n_q>0 && i <  static_cast<int>(Cs.cs.size()) ); i++) {
			// Check if the quark-line is open, in which case it has a q
			if (Cs.cs.at(i).open) {
				// The first quark is located and has position
				q_place.clear();
				q_place.push_back(i);
				q_place.push_back(0);
				// and number
				int q = Cs.at(q_place.at(0), q_place.at(1));

				// Locate same quark a second time
				// Loop over Quark_lines
				q_place2.clear();
				int i2 = i + 1; // Quark_line of second occurrence
				while (q_place2.empty()) { // As long as quark not found a second time
					if (Cs.cs.at(i2).at(Cs.cs.at(i2).ql.size() - 1) == q) {// If quark found, store place
						q_place2.push_back(i2);
						q_place2.push_back(Cs.cs.at(i2).ql.size() - 1);
					}
					i2++;
				}
				if (q_place2.empty()) {
					std::cerr << "Col_functions::contract_quarks(Cs1, Cs2): Found q " << q
							<< " only once in " << Cs << std::endl;
				}

				// Prepare new Quark_line
				// to be inserted at the place of found open Quark_line
				Quark_line new_Quark_line;
				Quark_line part2_new_Quark_line;
				// The first part of the new Quark_line should be the Quark_line
				// containing q in the conjugate
				new_Quark_line = Cs.cs.at(q_place2.at(0));

				// Erasing q in the end
				new_Quark_line.ql.erase(--new_Quark_line.ql.end());
				part2_new_Quark_line = Cs.cs.at(q_place.at(0));

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
				Cs.cs.at(i) = new_Quark_line;

				// Remove quark_line with q in Cs
				Cs.cs.erase((Cs.cs.begin() + q_place2.at(0)));
				i=-1; // reset to keep looking from the beginning in the new Cs (i will be increased to 0)
			}// end of if (open)

			n_q = Cs.n_quark();

		} // end of for, loop over quark_lines

	}

	return Cs;

}


Col_amp Col_functions::contract_quarks( Col_amp Ca1, Col_amp Ca2 ) const{

	if(Ca1.empty()){
		std::cerr << "Col_functions::contract_quarks: Expects non-empty Col_amps, got first argument "
				<< Ca1 << std::endl;
		assert(0);
	}
	if(Ca2.empty()){
		std::cerr << "Col_functions::contract_quarks: Expects non-empty Col_amps, got second argument "
				<< Ca2 << std::endl;
		assert(0);
	}

	Col_amp Ca_res;

	// Make sure the Col_strs are not empty "[]"=1, as all indices contracted
	Ca1.remove_empty_Col_strs();
	Ca2.remove_empty_Col_strs();

	// Loop over Col_strs, and contract quarks between all possible combinations
	  // Loop over Col_strs in Ca1
	  for(uint m1=0; m1 < Ca1.ca.size(); m1++ ){
	    // Loop over Col_strs in Ca2
	    for(uint m2=0; m2 < Ca2.ca.size(); m2++ ){
	    	Col_str Cs_tmp;
	    	Cs_tmp.contract_quarks( Ca1.ca.at(m1), Ca2.ca.at(m2));
	      Ca_res.ca.push_back( Cs_tmp );
	    }
	  }

	  return Ca_res;
}
*/

std::map< std::string, double > Col_functions::double_num( std::map< std::string, Polynomial > mem_map )  const{

	std::map< std::string, double > res;

	// Loop over entries
	std::map< std::string, Polynomial >::iterator iter=mem_map.begin();

	for( iter =  mem_map.begin(); iter !=  mem_map.end(); ++iter) {
		// Insert pair of string and the leading versions of the Polynomial
		double_num( (iter->second) );
		res.insert(std::make_pair( iter->first, double_num((iter->second)) ));
	}

	return res;
}


Polynomial Col_functions::Polynomial_cnum_num( const Polynomial & Poly ) const{

	// Store content in numerical part of the Monomial
	Monomial Mon;
	Mon.cnum_part = cnum_num(Poly);

	Polynomial res;
	res = res * Mon;

	return res;
}


cvec Col_functions::cnum_num( const Poly_vec & Pv )  const{

	// To contain the numerical result
	cvec res;

	// Loop over Polynomials in the vector, and add the numerical value
	// to the vector to return
	for (uint p = 0; p < Pv.size(); p++) {
		res.push_back(cnum_num(Pv.at(p)));
	}
	return res;
}


dvec Col_functions::double_num( const Poly_vec & Pv )  const{

	// To contain the numerical result
	dvec res;

	// Loop over Polynomials in the vector, and add the numerical value
	// to the vector to return
	for (uint p = 0; p < Pv.size(); p++) {
		res.push_back(double_num( Pv.at(p)) );
	}
	return res;
}


dmatr Col_functions::double_num( const Poly_matr & Pm )  const{
	// To contain the numerical result
	dmatr res;

	// Loop over Poly_vecs in the vector, and add the numerical value
	// to the vector to return
	for (uint pv = 0; pv < Pm.size(); pv++) {
		res.push_back(double_num( Pm.at(pv).pv ) );
	}
	return res;
}


dvec Col_functions::double_num( const std::vector<boost::shared_ptr<Polynomial> > & Pv )  const{

	// To contain the numerical result
	dvec res;

	// Loop over Polynomials in the vector, and add the numerical value
	// to the vector to return
	for (uint p = 0; p < Pv.size(); p++) {
		boost::shared_ptr<Polynomial>  the_pointer=Pv.at(p);
		// Want double of the polynomial which the pointer points at
		res.push_back( double_num(*the_pointer) );
	}
	return res;
}

Poly_vec Col_functions::Poly_vec_cnum_num( const Poly_vec & Pv)  const{

	// To contain the result
	Poly_vec res;

	// Loop over Polynomials in the vector, put each Polynomial to its numerical value
	for (uint p = 0; p < Pv.size(); p++) {
		res.push_back( Polynomial_cnum_num( Pv.at( p ) ) );
	}
	return res;
}


Poly_matr Col_functions::Poly_matr_cnum_num( const Poly_matr & Pm ) const {

	// To contain the result
	Poly_matr res_matr;

	// Loop over Polynomials in the matrix
	// and change each Polynomial to its numerical version
	for (uint v = 0; v < Pm.size(); v++) {
		res_matr.push_back( Poly_vec_cnum_num( Pm.at( v ).pv ));
	}

	return res_matr;
}


cmatr Col_functions::cnum_num( const Poly_matr & Pm )  const{

	// To contain the numerical result
	cmatr  res;

	// Loop over Polynomials in the vector, and add the numerical value
	// to the vector to return
	for (uint v = 0; v < Pm.size(); v++) {
		res.push_back(cnum_num( Pm.at(v).pv ));
	}
	return res;
}


dmatr Col_functions::double_num( const std::vector<std::vector<boost::shared_ptr<Polynomial> > > & Pm )  const{

	// To contain the numerical result
	dmatr  res;

	// Loop over Polynomials in the vector, and add the numerical value
	// to the vector to return
	for (uint v = 0; v < Pm.size(); v++) {
		res.push_back(double_num(Pm.at(v)));
	}
	return res;
}


Polynomial Col_functions::scalar_product( const Col_amp & Ca1 , const Col_amp & Ca2 ) const{
	//std::cout << "Col_functions::scalar_product, incoming Ca1 " << Ca1 <<" and Ca2 " << Ca2 << std::endl;

	if( !Ca1.Scalar.empty() and ( cnum_num(Ca1.Scalar).real()!=0 or cnum_num(Ca1.Scalar).imag()!=0 ) ){
		std::cerr << "Col_functions::scalar_product(Ca1,Ca2): "
				<< "Expects Col_amps with empty Scalar parts, but the Scalar of the first Col_amp was " <<
				Ca1.Scalar << std::endl;
		assert(0);
	}
	if( !Ca2.Scalar.empty() and ( cnum_num(Ca2.Scalar).real()!=0 or cnum_num(Ca2.Scalar).imag()!=0 ) ){
		std::cerr << "Col_functions::scalar_product(Ca1,Ca2): "
				<< "Expects Col_amps with empty Scalar parts, but the Scalar of the second Col_amp was " <<
				Ca2.Scalar << std::endl;
		assert(0);
	}
	// To contain the result
	Col_amp Ca_res;

	// Contract the quarks
	//Ca_res = contract_quarks(Ca1, Ca2);
	Ca_res.contract_quarks( Ca1, Ca2 );
	//std::cout << "Col_functions::scalar_product, contracted quarks " << Ca_res << std::endl;

	// Look for simple simplifications
	Ca_res.simplify();
	//Ca_res.simplify(); // why twice

	//std::cout << "Col_functions::scalar_product, simplified " << Ca_res << std::endl;

	// Contract the gluons
	Ca_res.contract_all_gluons();

	if (!Ca_res.empty()) {
		std::cerr << "Col_functions::scalar_product: terminating due to non-contracted indices."
				<< std::endl;
		std::cerr << "The Col_amp is " << Ca_res << std::endl;
	    std::cerr.flush();
		assert( 0 );
		return Ca_res.Scalar;
	}
	else
		return Ca_res.Scalar;
}


Polynomial  Col_functions::scalar_product( const Col_str & Cs1, const Col_str & Cs2 ) const{

	Col_str Cs_tmp;

	// Contract the quarks
	Cs_tmp.contract_quarks( Cs1, Cs2 );

	Col_amp Ca_tmp(Cs_tmp);

	// Contract the gluons
	Ca_tmp.contract_all_gluons();

	if ( !Ca_tmp.empty() ){
		std::cerr << "Col_functions::scalar_product: terminating due to non-contracted quark indices." <<std::endl;
		std::cerr << "The col_amp is " << Ca_tmp << std::endl;
		assert( 0 );
	}

	return Ca_tmp.Scalar;
}


Polynomial Col_functions::color_correlator( const Col_amp Ca, int p1, int p2, int g_new ) const{

	// The amplitudes after emission
	Col_amp Cai = emit_gluon(Ca, p1, g_new);
	Col_amp Caj = emit_gluon(Ca, p2, g_new);

	Polynomial res=0;

	res = scalar_product( Cai, Caj );

	return res;
}


int Col_functions::factorial( int i ) const{
	if(i<0) {
		std::cerr << "Col_functions::factorial: intended for int >=0, argument was " << i << std::endl;
		std::cerr.flush();
		assert( 0 );
	}
	if (i==0) return 1;
	return factorial(i-1)*i; // Recursive call
}


void Col_functions::write_out_dvec( const dvec & dv, std::string filename ) const {

	std::ofstream outfile(filename.c_str() );
	outfile << dv;
}


dmatr Col_functions::read_in_dmatr( std::string filename ) const {


	// Read in file
	std::ifstream fin(filename.c_str());

	// Check that file exists
	if( !fin ){
		std::cerr << "Col_functions::read_in_dmatr: The file "
				<< filename << " could not be opened." << std::endl;
		assert( 0 );

	}

	// Copy info from file to string
	std::string str((std::istreambuf_iterator<char>(fin)), std::istreambuf_iterator<char>());

	// Skip lines starting with #
	while(str.at(0)== '#'){
		while (str.at(0) != '\n'){
			str.erase(str.begin());
		}
		// erase endl sign(s)
		while(str.at(0)== '\n'){
			str.erase(str.begin());
		}
	}

	// First char in file should be '{'
	if (str.at(0) != '{') {
		std::cerr
		<< "Col_functions::read_in_dmatr: First char in matrix data file after comments should be '{', it was: "
		<< str.at(0) << std::endl;
		assert( 0 );
	}

	// Check that only allowed characters
	uint j = 0;
	while (j < str.size()) {

		if (!(str.at(j) == '+' or str.at(j) == '-' or str.at(j) == '.'
				or str.at(j) == '{' or str.at(j) == '}' or str.at(j) == '\n'
						or str.at(j) == ',' or str.at(j) == ' ' or str.at(j) == '0'
								or str.at(j) == '1' or str.at(j) == '2' or str.at(j) == '3'
										or str.at(j) == '4' or str.at(j) == '5' or str.at(j) == '6'
												or str.at(j) == '7' or str.at(j) == '8' or str.at(j) == '9')) {
			std::cerr
			<< "Col_functions::read_in_dmatr: A disallowed characters encountered in string for dmatr: "
			<< str.at(j) << ", in file " << filename <<  std::endl;
			std::cerr << "Col_functions::read_in_dmatr expects a numerical matrix." << std::endl;
			assert( 0 );
		}
		j++;
	}


	// Row to contain numbers
	dvec row;

	// To contain matrix of scalar products
	dmatr matr;

	// Read the string, starting from 0th element
	uint i = 0;
	while (i < str.size() - 2) {
		i += 1;

		// We may have to skip some chars
		while (i< str.size()-2 &&(str.at(i) == ',' or str.at(i) == '}' or str.at(i) == ' ' or str.at(i) == '\n' or str.at(i) == ' ' or str.at(i) == '{') ) i++;

		// String to make a number of, and double to contain number
		std::string num_str;
		num_str.clear();
		double num;

		// Keep reading the number while not ',' or '}'
		while ( i< str.size()-2 && (str.at(i) != ',' && str.at(i) != '}') )
		{
			num_str.push_back(str.at(i));
			i++;
		}

		// num_str contains the string to make a number of
		std::istringstream parton_str_st( num_str );
		parton_str_st >> num;

		// Add number to vector
		row.push_back(num);

		// If we have a new row
		if( i< str.size()-2 && str.at(i)=='}'){
			// Save row in matrix, and empty row
			matr.push_back(row);
			row.clear();

			// We may have to skip some chars
			while (i< str.size()-2 &&(str.at(i) == ',' or str.at(i) == '}' or str.at(i) == ' ' or str.at(i) == '\n' ) ) {
				i++;
			}
		}
		// Otherwise just keep on reading the next number in row
	}

	return matr;
}


dvec Col_functions::read_in_dvec( std::string filename ) const {

	// Read in file
	std::ifstream fin(filename.c_str());

	// Check that file exists
	if( !fin ){
		std::cerr << "Col_functions::read_in_dvec: The file "
				<< filename << " could not be opened." << std::endl;
		assert( 0 );
	}

	// Copy info from file to string
	std::string str((std::istreambuf_iterator<char>(fin)), std::istreambuf_iterator<char>());

	// Skip lines starting with #
	while(str.at(0)== '#'){
		while (str.at(0) != '\n'){
			str.erase(str.begin());
		}
		// erase endl sign(s)
		while(str.at(0)== '\n'){
			str.erase(str.begin());
		}
	}

	// First char in file should be '{'
	if (str.at(0) != '{') {
		std::cerr
		<< "Col_functions::read_in_dvec: First char in matrix data file after comments should be '{', it was: "
		<< str.at(0) << std::endl;
		assert( 0 );
	}


	// Check that only allowed characters
	uint j = 0;
	while (j < str.size()) {

		if (!(str.at(j) == '+' or str.at(j) == '-' or str.at(j) == '.'
				or str.at(j) == '{' or str.at(j) == '}' or str.at(j) == '\n'
						or str.at(j) == ',' or str.at(j) == ' ' or str.at(j) == '0'
								or str.at(j) == '1' or str.at(j) == '2' or str.at(j) == '3'
										or str.at(j) == '4' or str.at(j) == '5' or str.at(j) == '6'
												or str.at(j) == '7' or str.at(j) == '8' or str.at(j) == '9')) {
			std::cerr
			<< "Col_functions::read_in_dvec: A disallowed character encountered in string for dmatr: "
			<< str.at(j) << ", in file " << filename <<  std::endl;
			std::cerr << "Col_functions::read_in_dvec expects a numerical matrix." << std::endl;
			assert( 0 );
		}
		j++;
	}

	// Row to contain numbers
	dvec row;

	// To contain matrix of scalar products
	dmatr matr;

	// Read the string, starting from 0th element
	unsigned int i = 0;
	while (i < str.size() - 1) {
		i += 1;

		// We may have to skip some chars
		while (i< str.size()-2 &&(str.at(i) == ',' or str.at(i) == '}' or str.at(i) == ' ' or str.at(i) == '\n' or str.at(i) == ' ' or str.at(i) == '{') ) i++;

		// String to make a number of, and double to contain number
		std::string num_str;
		num_str.clear();
		double num;

		// Keep reading the number while not ',' or '}'
		while ( i< str.size()-1 && (str.at(i) != ',' && str.at(i) != '}') )
		{
			num_str.push_back(str.at(i));
			i++;
		}
		// now, at(i), there is either , or }

		// num_str contains the string to make a number of
		std::istringstream num_str_st( num_str );
		num_str_st >> num;

		// Add number to vector
		row.push_back(num);

		// Skip signs in end and make sure not to enter loop one extra time
		while( i<str.size() and ( str.at(i)==' ' or str.at(i)=='\n' or str.at(i)=='}' ) ) i++;
	}

	return row;

}


void Col_functions::write_out_dmatr( const dmatr & matr, std::string filename ) const {
	std::ofstream outfile(filename.c_str());
	outfile << matr;
}


std::list<int>::iterator operator+( std::list<int>::iterator x, int n ) {

  while(n>0){
    x++;
    n--;
  }
  while(n<0){
    x--;
    n++;
  }
  return x;
}

std::list<int>::iterator operator-( std::list<int>::iterator x, int n ) {

  while(n>0){
    x--;
    n--;
  }
  while(n<0){
    x++;
    n++;
  }
  return x;
}


std::list< Quark_line  >::iterator operator+( std::list < Quark_line >::iterator x, int n ){

  while(n>0){
    x++;
    n--;
  }
  while(n<0){
    x--;
    n++;
    }
  return x;
}

col_str::iterator operator-( col_str::iterator x, int n ){

  while(n>0){
    x--;
    n--;
  }
  while(n<0){
    x++;
    n++;
  }
  return x;
}


std::ostream& operator<<( std::ostream& out, std::vector<int> vec ){
  int max=vec.size();
  if(max==0)  out <<"{}";
  else{
    out <<"{";
    for (int i=0; i<max-1; i++){
      out << vec.at(i) << ",";
    }
    out << vec.at(max-1) <<"}";
  }
  return out;
}


std::ostream& operator<<( std::ostream& out, const cvec & cv ) {

	out << "{";
	// Loop over entries
	for (uint i = 0; i < cv.size(); i++) {
		// Print element
		std::cout.width(6);
		std::ostringstream outstr;
		outstr << cv.at(i);
		out << outstr.str();
		// If not last element print ","
		if (i < (cv.size() -1 )) out << ", ";
	}
	out << "}";
	return out;
}


std::ostream& operator<<(std::ostream& out, const dvec & dv) {

	out << "{";
	// Loop over entries
	for (uint i = 0; i < dv.size(); i++) {
		// Print element
		std::cout.width(6);
		std::ostringstream outstr;
		outstr << dv.at(i);
		out << outstr.str();
		// If not last element print ","
		if (i < (dv.size() -1 )) out << ", ";
	}
	out << "}";
	return out;
}


std::ostream& operator<<( std::ostream& out, const cmatr & cm ){

	out <<"{" << std::endl;
	// Loop over rows
	for(uint i=0; i< cm.size(); i++ ){
		out <<"{";
		// Loop over columns
		for(uint j=0; j< cm.at(i).size(); j++ ){
			// Print element
			std::cout.width( 6 );
			std::ostringstream outstr;
			outstr << cm.at(i).at(j);
			// If not last element print ","
			if (j<cm.at(i).size()-1 ) outstr << ",";
			out << outstr.str();
			//out << Poly_m.at(i).at(j) << "";
		}
		out <<"}";
		// If not last row, print ","
		if (i<cm.at(i).size()-1 ) out << ",";
		out << std::endl;
	}
	out <<"}" << std::endl;
  return out;
}


std::ostream& operator<<( std::ostream& out, const dmatr & matr ){

  out <<"{" << std::endl;
  // Loop over rows
  for(uint i=0; i< matr.size(); i++ ){
    out <<"{";
    // Loop over columns
    for(uint j=0; j< matr.at(i).size(); j++ ){
      // Print element
      std::ostringstream outstr;
      outstr.width( 20 );
      outstr.precision(16);
      // If the result is larger than accuracy print it out
      if( std::abs( matr.at(i).at(j) )> accuracy ){
    	  outstr << std::fixed << matr.at(i).at(j);
      }
      // otherwise is should probably be 0
      else
		  outstr << std::fixed << 0;
      // If not last element print ","
      if (j<matr.at(i).size()-1 ) outstr << ",";
      out << outstr.str();
    }
    out <<"}";
    // If not last row, print ","
    if (i<matr.at(i).size()-1 ) out << ",";
    out << std::endl;
  }
  out <<"}" <<std::endl;
  return out;
}


std::ostream& operator<<( std::ostream& out, std::pair<int, int> pair ) {

	out << "(";
	out << pair.first;
	out << ", ";
	out << pair.second;
	out << ")";

	return out;
}


} //end namespace ColorFull
