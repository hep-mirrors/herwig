// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CheckId class.
//

#include "CheckId.h"
#include "Herwig++/Utilities/HwDebug.h"

using namespace Herwig;
// using namespace Pythia7;


long CheckId::diquarkId(const long id1, const long id2) {
  if ( id1 * id2 < 0  ||  ! CheckId::isQuark(id1)  ||  ! CheckId::isQuark(id2) ) return 0;
  long ida = abs(id1);
  long idb = abs(id2);
  if (ida < idb) {
    ida = abs(id2);
    idb = abs(id1);
  }
  long idnew = ida*1000 + idb*100 + 1;
  // Diquarks made of quarks of the same type: uu, dd, ss, cc, bb, 
  // have spin 1, and therefore the less significant digit (which
  // corresponds to 2*J+1) is 3 rather than 1 as all other Diquarks.
  if (id1 == id2) idnew += 2;
  return id1 > 0 ? idnew : -idnew ;
}


bool CheckId::canBeBaryon(const long id1, const long id2, const long id3) {
  bool result = false;
  if (id3 == 0) {
    // In this case, to be a baryon, one component must be a (anti-)quark and
    // the other a (anti-)diquark; furthermore, only:
    //   quark - diquark   or  anti-quark - anti-diquark  
    // (that means same sign for both components) are allowed.
    result = id1*id2 > 0  &&
      ( ( isQuark(id1) && isDiquark(id2) ) || ( isDiquark(id1) && isQuark(id2) ) );
  } else {
    // In this case, to be a baryon, all three components must be (anti-)quarks
    // and with the same sign.
    result = id1*id2 > 0  &&  id2*id3 > 0  && 
      isQuark(id1) && isQuark(id2) && isQuark(id3);
  }
  return result;    
}


bool CheckId::hasBeauty(const long id1, const long id2, const long id3) {
  bool result = false;
  if ( id2 == 0  &&  id3 == 0 ) {
    result = ( abs(id1) == ParticleID::b                            ||   // quark
	       isDiquarkWithB(id1)                                  ||   // diquark
	       isMeson(id1) && (abs(id1)/100)%10  == ParticleID::b  ||   // meson
	       isBaryon(id1) && (abs(id1)/1000)%10 == ParticleID::b );   // baryon
  } else {
    result = ( abs(id1) == ParticleID::b  ||  isDiquarkWithB(id1)  || 
	       abs(id2) == ParticleID::b  ||  isDiquarkWithB(id2)  || 
	       abs(id3) == ParticleID::b  ||  isDiquarkWithB(id3) ); 
  }
  return result;
}  


bool CheckId::hasCharm(const long id1, const long id2, const long id3) {
  bool result = false;
  if ( id2 == 0  &&  id3 == 0 ) {
    result = ( abs(id1) == ParticleID::c                                || // quark
	       isDiquarkWithC(id1)                                      || // diquark
	       isMeson(id1) && ( (abs(id1)/100)%10   == ParticleID::c  ||  // meson no b
				 (abs(id1)/10)%10    == ParticleID::c ) || // meson with b
	       isBaryon(id1) && ( (abs(id1)/1000)%10 == ParticleID::c  ||  // baryon no b
				  (abs(id1)/100)%10  == ParticleID::c  ||  // baryon with b
				  (abs(id1)/10)%10   == ParticleID::c ) ); // baryon with b
  } else {
    result = ( abs(id1) == ParticleID::c  ||  isDiquarkWithC(id1)  || 
	       abs(id2) == ParticleID::c  ||  isDiquarkWithC(id2)  || 
	       abs(id3) == ParticleID::c  ||  isDiquarkWithC(id3) ); 
  }
  return result;
}  


bool CheckId::hasStrangeness(const long id1, const long id2, const long id3) {
  // If used with a single meson it returns always false: you have to use
  // the other overloaded method in the case of a meson.
  bool result = false;
  if ( id2 == 0  &&  id3 == 0 ) {
    result = ( abs(id1) == ParticleID::s                                || // quark
	       isDiquarkWithS(id1)                                      || // diquark
	       isBaryon(id1) && ( (abs(id1)/1000)%10 == ParticleID::s  || // baryon no c/b
				  (abs(id1)/100)%10  == ParticleID::s  || // baryon with c/b
				  (abs(id1)/10)%10   == ParticleID::s ) );// baryon with c/b
  } else {
    result = ( abs(id1) == ParticleID::s  ||  isDiquarkWithS(id1)  || 
	       abs(id2) == ParticleID::s  ||  isDiquarkWithS(id2)  || 
	       abs(id3) == ParticleID::s  ||  isDiquarkWithS(id3) ); 
  }  
  return result;         
}


bool CheckId::hasStrangeness(const long id, const double rnd) {
  // If used with anything different than  a meson it returns always false: 
  // you have to use the other overloaded method in these situations. 
  bool result = false;
  if ( isMeson(id) ) {  
    long id1 = (abs(id)/100)%10;
    long id2 = (abs(id)/10)%10;    
    if ( id1 != id2  &&  ( id1 == ParticleID::s  ||  id2 == ParticleID::s ) ) {
      result = true;  // Simple case, no mixing.
    } else if ( id1 == id2  &&  ( id1 == ParticleID::s  ||  id1 == 2 ) ) {

      // Consider now the most difficult case: mesons which are a
      // u ubar , d dbar , s sbar  admixtures.
      // First initialize the Octet-Singlet isoscalar mixing angles in degrees
      // (taken from Fortran Herwig subroutine HWIGIN and then used in subroutine 
      //  HWURES (where the tables for hadrons are also prepared.)
      double idealAngleMix = atan( 1.0 / sqrt(2.0) ) * 180.0 / acos(-1.0);
      double etamix = -23.0;          //  eta - eta'
      double phimix = +36.0;          //  phi - omega
      double h1mix  = idealAngleMix;  //  h_1(1380) - h_1(1170)
      double f0mix  = idealAngleMix;  //  missing - f_0(1370) 
      double f1mix  = idealAngleMix;  //  f_1(1420) - f_1(1285)
      double f2mix  = +26.0;          //  f'_2 - f_2
      // double omhmix = idealAngleMix;  //  missing - omega(1650)     NOT FOUND in Pythia7
      // double et2mix = idealAngleMix;  //  eta_2(1645) - eta_2(1870) NOT FOUND in Pythia7
      // double ph3mix = +28.0;          //  phi_3 - omega_3           NOT FOUND in Pythia7
      
      int order = 1;   // 1 for the first of the pair; 2 for the second one
      double angleMix  = 999.9;
      switch ( abs(id) ) {
	
      case ParticleID::eta :      angleMix = etamix; break;
      case ParticleID::etaprime : angleMix = etamix; order = 2; break;
	
      case ParticleID::phi :      angleMix = phimix; break;
      case ParticleID::omega :    angleMix = phimix; order = 2; break;
	
      case ParticleID::hprime_1 : angleMix = h1mix; break;
      case ParticleID::h_1 :      angleMix = h1mix; order = 2; break;
	
      case ParticleID::f_0 :      angleMix = f0mix; order = 2; break;
	
      case ParticleID::fprime_1 : angleMix = f1mix; break;
      case ParticleID::f_1 :      angleMix = f1mix; order = 2; break;
	
      case ParticleID::fprime_2 : angleMix = f2mix; break;
      case ParticleID::f_2 :      angleMix = f2mix; order = 2; break;
	
      }
      
      double prob = 1.0/3.0;
      if ( fabs( angleMix ) <= 360.0 ) {
	prob = probabilityMixing(angleMix,order);
      }	  
      // Because prob is the probability of the (|uubar>+|ddbar>)/sqrt(2)
      // component in  eta  then the probability for the |ssbar> component
      // is 1 - prob. Therefore, in order to consider  eta  as with strangeness
      // the random number (flat between 0 and 1) must be greater than prob.
      if (rnd > prob) result = true;  

      // Debugging
      if ( HERWIG_DEBUG_LEVEL >= 100 ) {    
	cerr << "CheckId::hasStrangeness : id=" << id << "  angleMix=" << angleMix 
	     << "  Prob(|ssbar>)=" << 1.0-prob << endl;
      }
      
    } // end if ( is a mixture of u ubar, d dbar, s sbar)
    
  } // end if ( isMesons(id) )
  
  return result;
         
}


double CheckId::probabilityMixing(const double angleMix, const int order) {
  // Calculate mixing weight for (|uubar>+|ddbar>)/sqrt(2) component.
  // Note: the convention used is (we are calling "eta" the first of the
  //       two pair (i.e. order = 1), and "eta'" the second one (order = 2)):
  //          eta  = cos(angleMix)*|etaOctet> - sin(angleMix)*|etaSinglet>
  //          eta' = sin(angleMix)*|etaOctet> + cos(angleMix)*|etaSinglet>
  //       with:
  //          etaSinglet = ( |uubar> + |ddbar> + |ssbar> ) / sqrt(3)
  //          etaOctect  = ( |uubar> + |ddbar> - 2*|ssbar> ) / sqrt(6)
  //       In deriving the formula below (taken from Fortran Herwig):
  // 
  //          probability for (|uubar>+|ddbar>)/sqrt(2) in eta =
  //            cos( mixAngle + atan( sqrt(2) )**2
  // 
  //       the following trigonometric identity has been used:
  // 
  //          ( cos(alpha) - sqrt(2)*sin(alpha) ) / sqrt(3)  =
  //          cos( alpha + atan( sqrt(2) )
  // 
  //       which can be easily proved by using:
  //          cos(alpha + beta) = cos(alpha)*cos(beta) - sin(alpha)*sin(beta)
  //          cos(beta) = 1 / sqrt( 1 + tan(beta)**2 )
  //          sin(beta) = tan(beta) / sqrt( 1 + tan(beta)**2 )
  //       where in our case tan(beta) = sqrt(2).
  //       Special cases for the mixing angle are the following:
  //
  //         1) mixAngle = atan(1/sqrt(2)) = 35.26 degrees (ideal mixing angle)
  //              |ssbar> component: 100% in eta , 0% in eta'
  //
  //         2) mixAngle = 0  (eta pure Octect, eta' pure Singlet)
  //              |ssbar> component: 2/3 in eta , 1/3 in eta' 
  //
  //         3) mixAngle = -atan(sqrt(2)) = -54.7 degrees 
  //              |ssbar> component: 0% in eta , 100% in eta'  

  static double pi = 3.14159265358979323846;
  double prob = 1.0;
  if ( order > 0 ) {
    if (order == 1)      prob = sqr( cos( angleMix*pi/180.0 + atan( sqrt(2.0) ) ) );
    else if (order == 2) prob = sqr( sin( angleMix*pi/180.0 + atan( sqrt(2.0) ) ) );
  }
  return prob;
}


