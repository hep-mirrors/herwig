// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HadronsSelector class.
//

#include "HadronsSelector.h"
#include "Pythia7/Interface/ClassDocumentation.h"
#include "Pythia7/Interface/Parameter.h"
#include "Pythia7/Persistency/PersistentOStream.h"
#include "Pythia7/Persistency/PersistentIStream.h"
#include "Pythia7/PDT/ParticleData.h"
#include "Pythia7/PDT/EnumParticles.h"
#include "Pythia7/Repository/EventGenerator.h"
#include "Herwig++/Utilities/HwDebug.h"
#include "Herwig++/Utilities/CheckId.h"
#include "Herwig++/Utilities/Kinematics.h"

using namespace Herwig;
// using namespace Pythia7;


HadronsSelector::~HadronsSelector() {}


void HadronsSelector::persistentOutput(PersistentOStream & os) const {
  os << _PwtDquark << _PwtUquark << _PwtSquark 
     << _PwtCquark << _PwtBquark << _PwtDIquark
     << _SngWt << _DecWt << _vecHad << _Pwt << _locHad; 
}


void HadronsSelector::persistentInput(PersistentIStream & is, int) {
  is >> _PwtDquark >> _PwtUquark >> _PwtSquark 
     >> _PwtCquark >> _PwtBquark >> _PwtDIquark
     >> _SngWt >> _DecWt >> _vecHad >> _Pwt >> _locHad;
}


ClassDescription<HadronsSelector> HadronsSelector::initHadronsSelector;
// Definition of the static class description member.


void HadronsSelector::Init() {

  static ClassDocumentation<HadronsSelector> documentation
    ("This class selects a proper pair (or single) hadrons...");

  static Parameter<HadronsSelector,double>
    interfacePwtDquark("PwtDquark","Weight for choosing a quark D",
		       &HadronsSelector::_PwtDquark, 0, 1.0, 0.0, 10.0);
  static Parameter<HadronsSelector,double>
    interfacePwtUquark("PwtUquark","Weight for choosing a quark U",
		       &HadronsSelector::_PwtUquark, 0, 1.0, 0.0, 10.0);
  static Parameter<HadronsSelector,double>
    interfacePwtSquark("PwtSquark","Weight for choosing a quark S",
		       &HadronsSelector::_PwtSquark, 0, 1.0, 0.0, 10.0);
  static Parameter<HadronsSelector,double>
    interfacePwtCquark("PwtCquark","Weight for choosing a quark C",
		       &HadronsSelector::_PwtCquark, 0, 1.0, 0.0, 10.0);
  static Parameter<HadronsSelector,double>
    interfacePwtBquark("PwtBquark","Weight for choosing a quark B",
		       &HadronsSelector::_PwtBquark, 0, 1.0, 0.0, 10.0);
  static Parameter<HadronsSelector,double>
    interfacePwtDIquark("PwtDIquark","Weight for choosing a DIquark",
			&HadronsSelector::_PwtDIquark, 0, 1.0, 0.0, 10.0);
  static Parameter<HadronsSelector,double>
    interfaceSngWt("SngWt","Weight for singlet baryons",
                  &HadronsSelector::_SngWt, 0, 1.0, 0.0, 10.0);
  static Parameter<HadronsSelector,double>
    interfaceDecWt("DecWt","Weight for decuplet baryons",
                  &HadronsSelector::_DecWt, 0, 1.0, 0.0, 10.0);

}


// ------------- PUBLIC MAIN METHODS ----------------------------

long HadronsSelector::
lightestHadron(const long id1, const long id2, const long id3) const {
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Hadronization ) {
    safetyCheck(id1,id2,id3);
  }
  long lightest = 0;
  if (id3 == 0) {    // The method assumes id3 == 0   
    int flav1 = convertIdToFlavour(id1);
    int flav2 = convertIdToFlavour(id2);
    if ( flav1 != 0  &&  flav2 != 0 ) {
      lightest = _vecHad[ _locHad[flav1][flav2].lightest ].id;
      int sign = signHadron(id1,id2,lightest);      
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Hadronization ) {  // Safety check, usually skipped.  
	if ( lightest == 0  ||  sign == 0 ) {  
	  generator()->logWarning( Exception("HadronsSelector::lightestHadron "
					     "***Zero lightest or sign*** ",
					     Exception::warning) );
	  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
	    generator()->log() << "         ===>" 
			       << "  id1=" << id1 << "  id2=" << id2 
			       << "  lightest=" << lightest << "  sign=" << sign 
			       << endl << endl;
	  }
	}
      } 
      lightest *= sign;
    } else {
      generator()->logWarning( Exception("HadronsSelector::lightestHadron "
					 "***Zero flav1 or flav2*** ",
					 Exception::warning) );
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
	generator()->log() << "         ===>" 
			   << "  id1=" << id1 << "  id2=" << id2 
			   << "  flav1=" << flav1 << "  flav2=" << flav2 
			   << endl << endl;
      }
    }
  }
  return lightest;
}


Energy HadronsSelector::
massLightestHadron(const long id1, const long id2, const long id3) const {
  long id = lightestHadron(id1,id2,id3);
  if ( ! getParticleData(id) ) {
    generator()->logWarning( Exception("HadronsSelector::massLightestHadron "
				       "***Lightest Hadron not found*** ", 
				       Exception::warning) );
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
      generator()->log() << "         ===>" 
			 << "  id1=" << id1 << "  id2=" << id2 << "  id3=" << id3 
			 << "  idHad=" << id << endl << endl;
    }
    return Energy();
  }
  return getParticleData(id)->mass();
}


pair<long,long> HadronsSelector::
lightestHadronsPair(const long id1, const long id2, const long id3) const {

  // This method has been implemented simply by calling twice the method  
  // lightestHadron, once with (id1.-idPartner) , and once with (id2,idPartner)
  // where  idPartner (with sign defined properly) is either the id of
  //  d  or  u . In fact, the idea is that whatever the flavour of id1 
  // and id2, no matter if (anti-)quark or (anti-)diquark, the lightest
  // pair of hadrons containing flavour id1 and id2 will have either 
  // flavour d or u, being the lightest quarks.
  // The method returns the pair (0,0) if anything goes wrong. 
  
  long idHad1 = 0, idHad2 = 0;
  
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Hadronization ) {
    safetyCheck(id1,id2,id3);
  }

  Charge ch = Charge();
  ch = getParticleData(id1)->charge() + getParticleData(id2)->charge();
  if ( id3 ) ch += getParticleData(id3)->charge();
  // Because charges are implemented as double, to protect against limited 
  // numerical precision, when comparing charges of particles w.r.t. the sum
  // of the charges of their components, an unphysical discrepancy of less
  // than a tenth of the smallest physical charge ( d quark: e/3) is allowed.
  Charge delta_ch = fabs( getParticleData(ParticleID::d)->charge() ) / 10.0; 
  
  // The method assumes id3==0 (otherwise we don't know how to proceed: a 
  // possible, trivial way would be to randomly select two of the three 
  // (anti-)quarks and treat them as a (anti-)diquark, reducing the problem
  // to two components as treated below.
  // In the normal (two components) situation, the strategy is the following:
  // treat in the same way the two possibilities:  (d dbar)  (i=0) and  
  // (u ubar)  (i=1)  as the pair quark-antiquark necessary to form a
  // pair of hadrons containing the input flavour  id1  and  id2; finally,
  // select the one that produces the lightest pair of hadrons, compatible
  // with the charge conservation contraint.
  if ( id3 == 0 ) {
    
    long vIdHad1[2] = {0,0}, vIdHad2[2] = {0,0};
    bool vOk[2] = {false, false};
    Energy vMassPair[2] = { Energy(), Energy() };
    
    for (int i = 0; i < 2; i++) {
	long idPartner = 0;
	if (i == 0) {
	  idPartner = ParticleID::d;
	} else {
	  idPartner = ParticleID::u;
	}
	
	// Change sign to idPartner (transform it into a anti-quark) if it is not
	// possible to form a meson or a baryon.
	if ( ! CheckId::canBeMeson(id1,-idPartner)  &&  
	     ! CheckId::canBeBaryon(id1,-idPartner) ) {  
	  idPartner *= -1;
	} 
	
	vIdHad1[i] = lightestHadron(id1, -idPartner);
	vIdHad2[i] = lightestHadron(id2, idPartner);

	// Test the charge conservation.
	if ( getParticleData(vIdHad1[i])  &&  getParticleData(vIdHad2[i])  &&
	     fabs( getParticleData(vIdHad1[i])->charge() + 
		   getParticleData(vIdHad2[i])->charge() - ch ) < delta_ch ) {
	  vOk[i] = true;
	  vMassPair[i] = getParticleData(vIdHad1[i])->mass() + 
	    getParticleData(vIdHad2[i])->mass();
	}  
    } // end of the if loop
    
    // Take the lightest pair compatible with charge conservation.
    if ( vOk[0]  &&  ( ! vOk[1]  ||  vMassPair[0] <= vMassPair[1] ) ) {      
	idHad1 = vIdHad1[0];
	idHad2 = vIdHad2[0];
    } else if ( vOk[1]  &&  ( ! vOk[0]  ||  vMassPair[1] < vMassPair[0] ) ) {      
	idHad1 = vIdHad1[1];
	idHad2 = vIdHad2[1];
    } else {
	// Sanity check (normally skipped) for debugging. 
	if ( HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Hadronization ) {
	  generator()->logWarning( Exception("HadronsSelector::lightestHadronsPair "
					     "***Not found lightest pair*** ",
					     Exception::warning) );
	  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
	    generator()->log() << "         ===>"
			       << " id1=" << id1 << " id2=" << id2 << " ch=" << ch << endl 
			       << " \t candidate " << vIdHad1[0] << " " << vIdHad2[0]
			       << " m=" << vMassPair[0] << " ch="
			       << ( getParticleData(vIdHad1[0]) && getParticleData(vIdHad2[0]) ?
				    getParticleData(vIdHad1[0])->charge() + 
				    getParticleData(vIdHad2[0])->charge() : -888 )
			       << " ok=" << vOk[0] << endl
			       << " \t candidate " << vIdHad1[1] << " " << vIdHad2[1]
			       << " m=" << vMassPair[1] << " ch="
			       << ( getParticleData(vIdHad1[1]) && getParticleData(vIdHad2[1]) ?
				    getParticleData(vIdHad1[1])->charge() + 
				    getParticleData(vIdHad2[1])->charge() : -888 )
			       << " ok=" << vOk[1] << endl
			       << endl << endl;
	  }      
	}
    }
    
  } // end of if (id3 == 0)
  
  return pair<long,long>(idHad1,idHad2);   
}


Energy HadronsSelector::
massLightestHadronsPair(const long id1, const long id2, const long id3) const {
  pair<long,long> pairId = lightestHadronsPair(id1,id2,id3);
  if ( ! getParticleData( pairId.first )  ||  ! getParticleData( pairId.second ) ) {
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Hadronization ) {
      generator()->logWarning( Exception("HadronsSelector::masslightestHadronsPair "
					 "***Lightest Hadron Pair not found*** ",
					 Exception::warning) );
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
	generator()->log() << "         ===>" << " pairId.first=" << pairId.first 
			   << " pairId.second=" << pairId.second << endl << endl;
      }      
    }
    return Energy();
  } 
  return ( getParticleData( pairId.first )->mass() + 
	   getParticleData( pairId.second )->mass() ); 
}


pair<long,long> HadronsSelector::
chooseHadronsPair(const Energy cluMass, const long id1, const long id2, const long id3) 
  throw(Veto, Stop, Exception) {

  // Kupco's method is used, rather than Brian's original one 
  // still in used in Fortran Herwig 6.3.
  // The idea is to build on the fly a table of all possible pairs
  // of hadrons (Had1,Had2) (that we can call "cluster decay channels")
  // which are kinematically above threshold  and have flavour 
  // Had1=(id1,-idQ), Had2=(idQ,id2), where idQ is the id of: 
  //    ---  d, u, s, c, b   
  //                        if either id1 or id2 is a diquark;      
  //    ---  d, u, s, c, b, dd, ud, uu, sd, su, ss, 
  //                        cd, cu, cs, cc, bd, bu, bs, bc, bb
  //                        if both id1 and id2 are quarks.
  // The weight associated with each channel is given by the product
  // of: the phase space available including the spin factor 2*J+1, 
  //     the constant weight factor for chosen idQ, 
  //     the octet-singlet isoscalar mixing factor, and finally 
  //     the singlet-decuplet weight factor.
     
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Hadronization ) {
    safetyCheck(id1,id2,id3);
  }

  pair<long,long> hadPair = pair<long,long>(0,0);
 
  if (id3 == 0) {    // The method assums id3 == 0
    
    // Define Kupco's table of possible cluster decay channels. 
    // The first element is [0] (not, [1], as vecHad).
    const int NmaxClusterDecayChannel = 3000;  // Should be enough!
    struct ClusterDecayChannel {
      long   idQ;     // id of the created flavour of the cluster decay channel.         
      long   idHad1;  // id of the first hadron of the cluster decay channel.              
      long   idHad2;  // id of the second hadron of the cluster decay channel.
      Energy weight;  // weight of the cluster decay channel.
    } kupcoTable[NmaxClusterDecayChannel]; 
    
    // If both id1 and id2 are not diquarks, then the possible Q values 
    // (of the Q Qbar pair from the vacuum) are assumed to be: 
    //    u, d, s, c, b, ud, dd, uu, sd, su, ss, 
    //                   cd, cu, cs, cc, bd, bu, bs, bc, bb
    // otherwise only u, d, s, c, b. 
    int numConsideredQuarksDiquarks = NumQuarksDiquarks;
    if ( CheckId::isDiquark(id1)  ||  CheckId::isDiquark(id2) ) {
      numConsideredQuarksDiquarks = B;
    }   
    
    // Fill Kupco's cluster decay channel table (starting from [0])
    // with all channels above threshold.
    int numChan = 0; // it counts the number of cluster decay channel.
    Energy maxWeight = Energy();  // store the maximum weight on the table.
    int flav1 = convertIdToFlavour(id1);
    int flav2 = convertIdToFlavour(id2);
    if ( flav1 != 0  &&  flav2 !=0 ) {
      for (int i=D; i <= numConsideredQuarksDiquarks; ++i) {
	if ( cluMass > ( _vecHad[ _locHad[flav1][i].lightest ].mass +  
			 _vecHad[ _locHad[flav2][i].lightest ].mass ) ) {         
	  // Loop over all hadron pairs with given flavour.
	  for (int iH1 = _locHad[flav1][i].first; iH1 <= _locHad[flav1][i].last; ++iH1) {
	    for (int iH2 = _locHad[flav2][i].first; iH2 <= _locHad[flav2][i].last; ++iH2) {
	      Energy pCmStar = Kinematics::pstarTwoBodyDecay( cluMass, _vecHad[iH1].mass,
							      _vecHad[iH2].mass );
	      if ( pCmStar > Energy() ) {  // the two hadrons are above threshold.
		if ( numChan >= NmaxClusterDecayChannel-1 ) { // send warning and skip channel.
		  generator()->logWarning( Exception("HadronsSelector::chooseHadronsPair "
						     "***Exceeded NmaxClusterDecayChannel*** ",
						     Exception::warning) );
		} else {
		  Energy weight = _Pwt[i] * pCmStar * 
		    _vecHad[iH1].overallWeight * _vecHad[iH2].overallWeight; 
		  if ( weight > maxWeight ) maxWeight = weight;
		  int signQ = 0;
		  long idQ = convertFlavourToId(i);
		  if ( ( CheckId::canBeMeson(id1,-idQ) || 
			 CheckId::canBeBaryon(id1,-idQ) ) &&   
		       ( CheckId::canBeMeson(idQ,id2)  || 
			 CheckId::canBeBaryon(idQ,id2)  ) ) {  
		    signQ = +1;
		  } else if ( ( CheckId::canBeMeson(id1,idQ)  || 
				CheckId::canBeBaryon(id1,idQ)  )  &&   
			      ( CheckId::canBeMeson(-idQ,id2) || 
				CheckId::canBeBaryon(-idQ,id2) ) ) {  
		    signQ = -1;
		  }
		  int signHad1 = 0, signHad2 = 0;
		  if ( signQ != 0 ) {
		    signHad1 = signHadron( id1, -signQ*idQ, _vecHad[iH1].id );
		    signHad2 = signHadron( id2,  signQ*idQ, _vecHad[iH2].id );
		  }
		  if ( signHad1 != 0  &&  signHad2 != 0 ) {
		    kupcoTable[numChan].idQ    = signQ * idQ;
		    kupcoTable[numChan].idHad1 = signHad1 * _vecHad[iH1].id;
		    kupcoTable[numChan].idHad2 = signHad2 * _vecHad[iH2].id;
		    kupcoTable[numChan].weight = weight;
		    numChan++;
		  } else {
		    generator()->logWarning( Exception("HadronsSelector::chooseHadronsPair "
						       "***Inconsistent Hadron*** ",
						       Exception::warning) );
		    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
		      generator()->log() << "         ===>"
					 << "  i=" << i << "  idQ=" << idQ << "  signQ=" << signQ 
					 << endl << "             "
					 << "  id1=" << id1 << "   " << -signQ*i << "   idHad1=" 
					 << _vecHad[iH1].id  << "  signHad1=" << signHad1 
					 << endl  << "             "
					 << "  id2=" << id2 << "   " << signQ*i  << "   idHad2=" 
					 << _vecHad[iH2].id  << "  signHad2=" << signHad2 << endl  
					 << endl;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    
    // Choose one decay channel.
    if ( numChan > 0 ) {
      int numTries = 0;
      int iChan;
      do {
	iChan = irnd(0,numChan);  // draw a flat random integer number [0,numChan[
      } while ( maxWeight*rnd() > kupcoTable[iChan].weight  &&  ++numTries < MaxNumTries );
      // Take the lightest pair if too many attempts failed.
      if ( numTries >= MaxNumTries ) {
	hadPair = lightestHadronsPair(id1,id2);
	generator()->logWarning( Exception("HadronsSelector::chooseHadronsPair"
					   "***Too many failed attempts, take lightest pair*** ", 
					   Exception::warning) );
      } else {
	hadPair = pair<long,long>( kupcoTable[iChan].idHad1 , kupcoTable[iChan].idHad2 );
      }
    }
    
    // Safety checks and debugging infos, usually skipped.
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Hadronization ) {
      if ( numChan == 0 ) {
	generator()->logWarning( Exception("HadronsSelector::chooseHadronsPair "
					   "***numChan = 0  in kupcoTable*** ",
					   Exception::warning) );
	if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
	  generator()->log() << "         ===>" 
			     << " cluMass=" << cluMass 
			     << "   id1=" << id1 << "   id2=" << id2 << "   id3=" << id3 
			     << endl << endl;
	}
      } else if ( maxWeight < 1.0e-9*GeV ) {
	generator()->logWarning( Exception("HadronsSelector::chooseHadronsPair "
					   "***very low maxWeight  in kupcoTable*** ",
					   Exception::warning) );
	if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
	  generator()->log() << "         ===>" 
			     << " cluMass=" << cluMass 
			     << "   id1=" << id1 << "   id2=" << id2 << "   id3=" << id3 
			     << endl << endl;
	}
      }
      if ( HERWIG_DEBUG_LEVEL == HwDebug::extreme_Hadronization ) {    
	generator()->log() << "\t --- Kupco's Table --- : cluMass=" << cluMass 
			   << "  id1=" << id1 << "  id2=" << id2 << endl
			   << "\t numChan=" << numChan << "   maxWeight=" << maxWeight 
			   << endl;
	for (int ich=0; ich < numChan; ++ich) {
	  generator()->log() << "\t \t" << ich << "   idQ=" << kupcoTable[ich].idQ
			     << "   idHad1=" << kupcoTable[ich].idHad1 
			     << "   idHad2=" << kupcoTable[ich].idHad2
			     << "   weight=" << kupcoTable[ich].weight 
			     << endl;
	}
      }      
    }
    
  } // end if (id3 == 0)
  
  return hadPair;
  
}


// ------------- PRIVATE ACCESSORY METHODS ----------------------------

void HadronsSelector::safetyCheck(const long id1, const long id2, const long id3) const
  throw(Veto, Stop, Exception) {
  // Check if the components can form a meson or a baryon.
  if ( ! CheckId::canBeMeson(id1,id2)  &&  ! CheckId::canBeBaryon(id1,id2,id3) ) {  
    generator()->logWarning( Exception("HadronsSelector::safetyCheck "
				       "***The cluster is inconsistent*** ",
				       Exception::warning) );
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
      generator()->log() << "         ===>" 
			 << " id1=" << id1 << " id2=" << id2 << " id3=" << id3 
			 << endl << endl;
    }
  }
}


int HadronsSelector::convertIdToFlavour(const long id) const {
  int flavour = 0;
  if ( abs(id) >= D  &&  abs(id) <= B ) { 
    flavour = abs(id);
  } else {
    switch ( abs(id) ) {
    case 1103: flavour = DD;  break;
    case 2101:
    case 2103: flavour = UD;  break;  
    case 2203: flavour = UU;  break;
    case 3101:
    case 3103: flavour = SD;  break;
    case 3201:
    case 3203: flavour = SU;  break;
    case 3303: flavour = SS;  break;
    case 4101:
    case 4103: flavour = CD;  break;
    case 4201:
    case 4203: flavour = CU;  break;
    case 4301:
    case 4303: flavour = CS;  break;
    case 4403: flavour = CC;  break;
    case 5101:
    case 5103: flavour = BD;  break;
    case 5201:
    case 5203: flavour = BU;  break;
    case 5301:
    case 5303: flavour = BS;  break;
    case 5401:
    case 5403: flavour = BC;  break;
    case 5503: flavour = BB;  break;
    }
  }
  return flavour;
}


long HadronsSelector::convertFlavourToId(const int flavour) const {
  // Notice that we are always converting diquarks to the 
  // spin 0 states, when available, never to spin 1.
  // This is just a simplifications that hopefully should
  // be harmless (if not, to change is quite straighforward,
  // for example randomly assigned 50% to spin 0 and 50% to spin 1).
  long id = 0;
  switch ( flavour ) {
  case D:   id  = 1;     break;
  case U:   id  = 2;     break;
  case S:   id  = 3;     break;
  case C:   id  = 4;     break;
  case B:   id  = 5;     break;
  case DD:  id  = 1103;  break;
  case UD:  id  = 2101;  break;
  case UU:  id  = 2203;  break;
  case SD:  id  = 3101;  break;
  case SU:  id  = 3201;  break;
  case SS:  id  = 3303;  break;
  case CD:  id  = 4101;  break;
  case CU:  id  = 4201;  break;
  case CS:  id  = 4301;  break;
  case CC:  id  = 4403;  break;
  case BD:  id  = 5101;  break;
  case BU:  id  = 5201;  break;
  case BS:  id  = 5301;  break;
  case BC:  id  = 5401;  break;
  case BB:  id  = 5503;  break;
  }
  return id;
}


int HadronsSelector::signHadron(const int idQ1, const int idQ2, const int idHad) const {

  // This method receives in input three PDG ids, whose the
  // first two have proper signs (corresponding to particles, id > 0, 
  // or antiparticles, id < 0 ), whereas the third one must
  // be always positive (particle not antiparticle),
  // corresponding to:
  //  --- quark-antiquark, or antiquark-quark, or
  //      quark-diquark, or diquark-quark, or
  //      antiquark-antidiquark, or antidiquark-antiquark
  //      for the first two input (idQ1, idQ2);
  //  --- meson or baryon for the third input (idHad): 
  // The method returns:
  //  --- + 1  if the two partons (idQ1, idQ2) are exactly
  //           the constituents for the hadron idHad;
  //  --- - 1  if the two partons (idQ1, idQ2) are exactly
  //           the constituents for the anti-hadron -idHad;
  //  --- + 0  otherwise.
  // The method it is therefore useful to decide the
  // sign of the id of the produced hadron as appeared 
  // in the vector _vecHad (where only hadron idHad > 0 are present)  
  // given the two constituent partons.

  // Safety check, usually skipped.
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Hadronization ) {
    if ( ! ( getParticleData(idQ1)  &&  getParticleData(idQ2)  &&  getParticleData(idHad)  &&
	     ( CheckId::canBeMeson(idQ1,idQ2)  && CheckId::isMeson(idHad)  ) ||
	     ( CheckId::canBeBaryon(idQ1,idQ2) && CheckId::isBaryon(idHad) ) ) ) {
      generator()->logWarning( Exception("HadronsSelector::signHadron "
					 "***The input is inconsistent*** ",
					 Exception::warning) );
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
	generator()->log() << "         ===>" 
			   << " idQ1=" << idQ1 << "  idQ2=" << idQ2 << "  idHad=" << idHad 
			   << endl << endl;
      }
    }
  }

  int sign = 0;
  if ( getParticleData(idQ1)  &&  getParticleData(idQ2)  &&  getParticleData(idHad) ) {
    Charge chargeIn  = getParticleData(idQ1)->charge() + getParticleData(idQ2)->charge();
    Charge chargeOut = getParticleData(idHad)->charge();
    Charge delta_ch  = fabs( getParticleData(D)->charge() ) / 10.0; 
    if ( fabs(chargeIn-chargeOut) < delta_ch  &&  fabs(chargeIn+chargeOut) > delta_ch ) {
      sign = +1;
    } else if ( fabs(chargeIn-chargeOut) > delta_ch  &&  fabs(chargeIn+chargeOut) < delta_ch ) {
      sign = -1;
    } else if ( fabs(chargeIn) < delta_ch  &&  fabs(chargeOut) < delta_ch ) {  
      // In the case of same null charge, there are four cases:
      //  i) K0-like mesons, B0-like mesons, Bs-like mesons
      //     the PDG convention is to consider them "antiparticle" (idHad < 0) 
      //     if the "dominant" (heavier) flavour (respectively, s, b)
      //     is a quark (idQ > 0): for instance, B0s = (b, sbar) has id < 0
      //     Remember that there is an important exception for K0L (id=130) and
      //     K0S (id=310): they don't have antiparticles, therefore idHad > 0
      //     always. We use below the fact that K0L and K0S are the unique
      //     hadrons having 0 the first (less significant) digit of their id.
      //  2) D0-like mesons: the PDG convention is to consider them "particle"
      //     (idHad > 0) if the charm flavour is carried by a c: (c,ubar) has id>0
      //  3) the remaining mesons should not have antiparticle, therefore their
      //     sign is always positive.
      //  4) for baryons, that is when one of idQ1 and idQ2 is a (anti-) quark and 
      //     the other one is a (anti-) diquark the sign is negative when both
      //     constituents are "anti", that is both with id < 0; positive otherwise.
      if ( CheckId::isQuark(idQ1) && CheckId::isQuark(idQ2) ) {  // meson
	int idQa = abs(idQ1), idQb = abs(idQ2), dominant = idQ2;
	if (idQa > idQb) {
	  idQa = abs(idQ2); idQb = abs(idQ1); dominant = idQ1;
	}
	if ( ( idQa==D && idQb==S ) || ( idQa==D && idQb==B ) || ( idQa==S && idQb==B ) ) {
	  if ( dominant < 0  ||  idHad % 10 == 0 ) {  // idHad%10 is zero for K0L,K0S
	    sign = +1;
	  } else if ( dominant > 0 ) {
	    sign = -1;
	  }
	} else if ( idQa==U && idQb==C ) {
	  if ( dominant > 0 ) {
	    sign = +1;
	  } else if ( dominant < 0 ) { 
	    sign = -1;
	  } 
	} else if ( idQa==idQb ) {
	  sign = +1;
	}
      } else if ( CheckId::isDiquark(idQ1) || CheckId::isDiquark(idQ2) ) { // baryon
	if ( idQ1 > 0  &&  idQ2 > 0 ) {
	  sign = +1;
	} else if ( idQ1 < 0  &&  idQ2 < 0 ) {
	  sign = -1;
	}
      }
    }
  }
  return sign;
}


// ------------- INITIALIZATION  ----------------------------

void HadronsSelector::initialize() {

  // This is the main method to initialize the hadron data
  // (mainly the weights associated to each hadron, taking
  //  into account its spin, eventual isoscalar-octect
  //  mixing, singlet-decuplet factor).
  // Plenty of debugging information is also store in the
  // log file.  

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
    generator()->log() << "  --- BEGIN INITIALIZE --- " << endl;
  }      

  // Initialize  _Pwt  and  _RepWt as is done in Fortran Herwig.
  // Note on the factor 0.5 below for diquarks of different quarks:
  //      the rationale for that factor can be explained by taking
  //      as a prototype example that of exact SU(2) limit, in 
  //      which:  m(up) = m(down)   and 
  //              m(proton) = m(neutron) = m(Delta)     
  //      and we will have equal numbers of u and d quarks produced.
  //      Suppose that we weight 1 the diquarks made of the same 
  //      quark and 1/2 those made of different quarks, the fractions
  //      of u and d baryons (p, n, Delta) we get are the following:
  //        --- Delta++  : 1 possibility only  u uu  with weight 1
  //        --- Delta-   : 1 possibility only  d dd  with weight 1
  //        --- p,Delta+ : 2 possibilities:    u ud  with weight 1/2
  //                                           d uu  with weight 1
  //        --- n,Delta0 : 2 possibilities:    d ud  with weight 1/2
  //                                           u dd  with weight 1
  //      In the latter two cases, we have to take into account the 
  //      fact that  p  and  n  have spin 1/2 whereas  Delta+  and  Delta0
  //      have spin 3/2 therefore from phase space we get a double weight 
  //      for  Delta+  and  Delta0  relative to  p  and  n  respectively.
  //      Therefore the relative amount of these baryons that is
  //      produced is the following:
  //       # p = # n = ( 1/2 + 1 ) * 1/3 = 1/2
  //       # Delta++ = # Delta- = 1 = ( 1/2 + 1) * 2/3 # Delta+ = # Delta0
  //      which is correct, and therefore the weight 1/2 for the
  //      diquarks of different types of quarks is justified (at least
  //      in this limit of exact SU(2) ).

  _Pwt[0]  = 0.0; // not used.
  _Pwt[D]  = _PwtDquark;
  _Pwt[U]  = _PwtUquark;
  _Pwt[S]  = _PwtSquark;
  _Pwt[C]  = _PwtCquark;
  _Pwt[B]  = _PwtBquark;
  _Pwt[DD] =       _PwtDIquark * _PwtDquark * _PwtDquark;
  _Pwt[UD] = 0.5 * _PwtDIquark * _PwtUquark * _PwtDquark;
  _Pwt[UU] =       _PwtDIquark * _PwtUquark * _PwtUquark;
  _Pwt[SD] = 0.5 * _PwtDIquark * _PwtSquark * _PwtDquark;
  _Pwt[SU] = 0.5 * _PwtDIquark * _PwtSquark * _PwtUquark;
  _Pwt[SS] =       _PwtDIquark * _PwtSquark * _PwtSquark;
  _Pwt[CD] = 0.5 * _PwtDIquark * _PwtCquark * _PwtDquark;
  _Pwt[CU] = 0.5 * _PwtDIquark * _PwtCquark * _PwtUquark;
  _Pwt[CS] = 0.5 * _PwtDIquark * _PwtCquark * _PwtSquark;
  _Pwt[CC] =       _PwtDIquark * _PwtCquark * _PwtCquark;
  _Pwt[BD] = 0.5 * _PwtDIquark * _PwtBquark * _PwtDquark;
  _Pwt[BU] = 0.5 * _PwtDIquark * _PwtBquark * _PwtUquark;
  _Pwt[BS] = 0.5 * _PwtDIquark * _PwtBquark * _PwtSquark;
  _Pwt[BC] = 0.5 * _PwtDIquark * _PwtBquark * _PwtCquark;
  _Pwt[BB] =       _PwtDIquark * _PwtBquark * _PwtBquark;

  double _Repwt[Lmax][Jmax][Nmax];  // Weights for excited mesons;  
  for (int l = 0; l < Lmax; ++l ) {
    for (int j = 0; j < Jmax; ++j) {
      for (int n = 0; n < Nmax; ++n) {
	_Repwt[l][j][n] = 1.0;
      }
    }
  }

  // Initialize to the default value.  
  for (int i=0; i < NumHadrons; ++i) {
    _vecHad[i].id            = 0;
    _vecHad[i].ptrData       = tPDPtr();
    _vecHad[i].swtef         = 1.0;
    _vecHad[i].wt            = 1.0;
    _vecHad[i].overallWeight = 0.0;           // This is important to skip uknown or wrong input.
    _vecHad[i].mass          = 0.999999e+9*GeV;   // Just a crazy large value!
  }
  for (int i=0; i <= NumQuarksDiquarks; ++i) {
    for (int j=0; j <= NumQuarksDiquarks; ++j) {
      _locHad[i][j].first    = -997;          // Just crazy large negative values,
      _locHad[i][j].last     = -998;          // furthermore wrongly ordered!
      _locHad[i][j].lightest = -999;
    } 
  }
 
  int countHadrons = fillDataHadrons();

  // Calculate the overallWeight, and the other remaining stuff.
  // Notice that if the id is < 0, it will be forced to be positive, because
  // we are following the convention to have only hadrons ( id > 0 ), 
  // not anti-hadrons ( id < 0 ) in the vector _vecHad.
  // We use the PDG id convention to identify:
  //  --- Singlet (Lambda-like) baryons:  J=1/2 and the two less
  //      significant digits associated to quark flavours are ordered
  //      differently than usual: the less significant of the two is
  //      bigger than the more significant one;
  //  --- Decuplet baryons: J=3/2
  //  --- Angular or Radial excited mesons: L, or J, or N greater than 0 

  multiset<long> unknownIds; // Set of hadrons not found in Pythia7
  for (int i=1; i <= countHadrons; ++i) {
    if ( _vecHad[i].id < 0 ) _vecHad[i].id *= -1; // force positive id.
    int nj = _vecHad[i].id % 10;
    if ( nj == 0 ) nj = 1;  // K0L and K0S only have nj == 0 
    if ( nj == 2  ||  nj == 4 ) {     // Baryon : J = 1/2 or 3/2
      if ( nj == 2 ) {
	int nq3 = ( _vecHad[i].id / 10 )  % 10;
        int nq2 = ( _vecHad[i].id / 100 ) % 10;
        if ( nq2 < nq3 ) {
	  _vecHad[i].swtef = pow( _SngWt, 2 );   // Singlet (Lambda-like) baryon
	}
      } else {                    
	 _vecHad[i].swtef = pow( _DecWt, 2 );    // Decuplet baryon
      }
    } else if ( 2*(nj/2) != nj ) {   // Meson
      int j  = (nj - 1) / 2;                     // Total angular momentum
      int nl = ( _vecHad[i].id / 10000 )  % 10;  // related to Orbital angular momentum l
      int l  = -999;  
      int n  = ( _vecHad[i].id / 100000 ) % 10;  // Radial excitation
      if ( j == 0 ) {
	l = nl;
      } else if ( nl == 0 ) {
	l = j - 1;
      } else if ( ( nl == 1) || ( nl == 2) ) {
	l = j;
      } else if ( nl == 3 ) {
	l = j + 1;
      }
      if (  ( l > 0  ||  j > 0  ||  n > 0 )  &&
            l >= 0  &&  l < Lmax  &&  j < Jmax  &&  n < Nmax ) {
	_vecHad[i].swtef = pow( _Repwt[l][j][n], 2 );  // Angular or Radial excited meson
      }    
    }
    _vecHad[i].ptrData = getParticleData( _vecHad[i].id ); 
    if ( _vecHad[i].ptrData != tPDPtr() ) {
      _vecHad[i].overallWeight = nj * _vecHad[i].wt * _vecHad[i].swtef;
      _vecHad[i].mass          = _vecHad[i].ptrData->mass();
    } else {
      unknownIds.insert( _vecHad[i].id );
    }
  }
 
  // Find the lightest hadron position for each class.
  for (int i=D; i <= B; ++i) {
    for (int j=i; j <= NumQuarksDiquarks; ++j) {
      int lightestPos = _locHad[i][j].first;
      for (int k = lightestPos+1; k <= _locHad[i][j].last; ++k) {
	if ( _vecHad[k].mass < _vecHad[lightestPos].mass ) lightestPos = k;
      } 
      _locHad[i][j].lightest = lightestPos;
    } 
  }

  // Symmetrize the table _locHad, just for the easy of use.
  for (int i=D; i <= B; ++i) {
    for (int j=i+1; j <= NumQuarksDiquarks; ++j) {
      if ( i != j) _locHad[j][i] = _locHad[i][j];
    }
  }

  // Sanity checks and debugging infos (normally skipped)
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Hadronization ) {
    if ( unknownIds.size() != 0 ) {
      generator()->logWarning( Exception("HadronsSelector::initialize "
					 "***Hadrons UNKNOWN is Pythia7: they will be skipped! *** ",
					 Exception::warning) );
      for ( multiset<long>::const_iterator iter = unknownIds.begin();
	    iter != unknownIds.end(); ++iter) {
	generator()->log() << "\t UNKNOWN  id = " << *iter << endl; 
      }
    }
    for (int i=1; i <= countHadrons; ++i) {
      if ( ( unknownIds.find( _vecHad[i].id ) == unknownIds.end() ) &&
           ( _vecHad[i].id == 0   ||  _vecHad[i].ptrData == tPDPtr()  || 
	     _vecHad[i].wt < 0.0  ||  _vecHad[i].swtef < 0.0  ||  _vecHad[i].overallWeight < 0  ||
	     ( ! CheckId::isMeson( _vecHad[i].id )  &&  ! CheckId::isBaryon( _vecHad[i].id ) ) ) ) { 
	generator()->logWarning( Exception("HadronsSelector::initialize "
					   "***Wrong input data: it will be ignored!*** ",
					   Exception::warning) );
	if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
	  generator()->log() << "         ===>" << "  i=" << i << "  id=" << _vecHad[i].id 
			     << "  ptrData=" << _vecHad[i].ptrData
                             << "  wt=" << _vecHad[i].wt << "  swtef=" << _vecHad[i].swtef
	                     << "  overallWeight=" << _vecHad[i].overallWeight
                             << "  mass=" << _vecHad[i].mass << endl << endl;
	}
      } 
    }
    map<long,int> allIds; // To check the number of times each hadron appear in _vecHad
    for (int i=D; i <= B; ++i) {
      for (int j=i; j <= NumQuarksDiquarks; ++j) {
	if ( (_locHad[i][j].first < 0  &&  _locHad[j][i].first < 0) || 
	     (_locHad[i][j].last  < 0  &&  _locHad[j][i].last  < 0) || 
	     (_locHad[i][j].last  < _locHad[i][j].first  &&  
	      _locHad[j][i].last  < _locHad[j][i].first) ) {
	  generator()->log() << "\t NO ENTRY : "  << " i=" << i << " j=" << j 
			     << " first=" << _locHad[i][j].first 
			     << " last =" << _locHad[i][j].last << endl;
	} else {
	  // Check that the same hadron is not repeated more than ones in each
          // hadron class; furthermore, in the overall hadron list (_vecHad)
          // each hadron should not appear more than 3 times (for mesons, this
          // can happen only for Octet-Singlet isoscalar mixing, like eta-eta',
          // which contribute to classes DD, UU, SS; for baryons, this can 
          // happen for those made by three different quark flavours: for
          // example for Lambda0: D SU, U SD, S UD.
	  map<long,int> classIds; 
	  for (int k=_locHad[i][j].first; k <= _locHad[i][j].last; ++k) {
	    if ( classIds.find( _vecHad[k].id ) != classIds.end() ) {
              ( classIds.find( _vecHad[k].id ) )->second++; 
	    } else {
	      classIds.insert( pair<long,int>(_vecHad[k].id,1) );
	      if ( allIds.find( _vecHad[k].id ) != allIds.end() ) {
		( allIds.find( _vecHad[k].id ) )->second++; 
	      } else {
		allIds.insert( pair<long,int>(_vecHad[k].id,1) );
	      }
	    }
	  }
	  for ( map<long,int>::const_iterator iter = classIds.begin(); iter != classIds.end(); ++iter) {
	    if ( iter->second > 1 ) {
	      generator()->logWarning( Exception("HadronsSelector::initialize "
						 "***Hadron repeated in the same class*** ", 
						 Exception::warning) );
	      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
		generator()->log() << "         ===>" << "  i=" << i << "  j=" << j 
                                   << "  id=" << iter->first << "  times=" << iter->second 
				   << endl << endl;
	      }
	    }
	  }
	}  
      } 
    }
    // Check now for repetitions in the all _vecHad 
    for ( map<long,int>::const_iterator iter = allIds.begin(); iter != allIds.end(); ++iter) {
      if ( iter->second > 3 ) {
	generator()->logWarning( Exception("HadronsSelector::initialize "
					   "***Hadron repeated overall too many times*** ", 
					   Exception::warning) );
	if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
	  generator()->log() << "         ===>" << "  id=" << iter->first 
			     << "  times=" << iter->second << endl << endl;
	}
      }
    }
    // General debugging infos
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
      generator()->log() << "   Total number of hadrons in the table: " << countHadrons << endl;
      for (int i=D; i <= B; ++i) {
	for (int j=i; j <= NumQuarksDiquarks; ++j) {
	  generator()->log() << "   -----------  i=" << i << "   j=" << j << " -------- from " 
			     << _locHad[i][j].first << "  to " << _locHad[i][j].last 
			     << "  lightest " << _locHad[i][j].lightest << endl;
	  for (int k=_locHad[i][j].first; k <= _locHad[i][j].last; ++k) {
	    generator()->log() << "\t   " << k << "   " 
			       << (_vecHad[k].ptrData ? _vecHad[k].ptrData->PDGName() : "UNKNOWN")
                               << "  id=" << _vecHad[k].id 
			       << "  wt=" << _vecHad[k].wt << "  swtef=" << _vecHad[k].swtef
			       << "  overallWeight=" << _vecHad[k].overallWeight 
                               << "  mass=" << _vecHad[k].mass << endl;
	  }
	}
      }
      generator()->log() << "  --- END INITIALIZE --- " << endl << endl;
    }      
  }
}


int HadronsSelector::fillDataHadrons() {

  // Here the hadrons data must be filled.

  // Octet-Singlet isoscalar mixing angles in degrees 
  // (taken from Fortran Herwig subroutine HWIGIN)
  const double idealAngleMix = atan( 1.0 / sqrt(2.0) ) * 180.0 / pi;
  const double etamix = -23.0;          //  eta - eta'
  const double phimix = +36.0;          //  phi - omega
  const double h1mix  = idealAngleMix;  //  h_1(1380) - h_1(1170)
  const double f0mix  = idealAngleMix;  //  missing - f_0(1370) 
  const double f1mix  = idealAngleMix;  //  f_1(1420) - f_1(1285)
  const double f2mix  = +26.0;          //  f'_2 - f_2

  //***LOOKHERE*** The following three mixing cases are at the 
  //               moment comment out because the corresponding
  //               particles are not (yet) present in Pythia7/PDT/EnumParticles.h.
  //               The PDG Ids and masses of these particles are
  //               (according to Herwig 6.3) the following: 
  //                  omega(1650) : id = 30223 ,  m = 1649 MeV 
  //                  eta_2(1645) : id = 10225 ,  m = 1632 MeV 
  //                  eta_2(1870) : id = 10335 ,  m = 1854 MeV 
  //                  phi_3       : id = 337   ,  m = 1854 MeV 
  //                  omega_3     : id = 227   ,  m = 1667 MeV 
  // const double omhmix = idealAngleMix;  //  missing - omega(1650)     
  // const double et2mix = idealAngleMix;  //  eta_2(1645) - eta_2(1870) 
  // const double ph3mix = +28.0;          //  phi_3 - omega_3          

  double angleMix = 0.0;  
  int order = 0;         

  int i = 1;  // _vecHad[0] is left empty; start from 1

  //==================== Begin filling table of hadrons ===================
  //
  //           90 Classes: 15 Meson classes + 75 Baryon classes
  // Note:
  //      --- List of hadrons taken from PDG 2000.
  //          Many of these hadrons are not yet present in Pythia7: 
  //          you can equivalently comment out the relevant lines 
  //          of code, or use (uncomment) them: in the latter case
  //          you will receive a warning message and zero overall weights
  //          will be given to them (so they will never be produced).
  //      --- Only mesons of type   (X, Y)   with  X <= Y  and
  //          only baryons of type  (X, YZ)  with  Z <= Y  (whatever X)
  //          (where X, Y, Z can be  U, D, S, C, B; notice the order 
  //                  D < U < S < C < B )
  //          are filled (directly or indirectly, see below) in the table 
  //          in this method (the symmetric ones will be filled automatically
  //          by the method initialize() ).
  //      --- For mesons, the  U U  class of mesons is automatically copied 
  //          from the  D D  ones.
  //      --- For baryons, the  "X Y Z"  class of baryons, where two of X, Y, Z 
  //          are "combined" into a diquark, is automatically copied from the 
  //          class  "L MN"  where:
  //                   L = MIN ( X, Y, Z)
  //                   M = MAX ( X, Y, Z)
  //                   N = the remaining one  
  //          Example: (U SD) and (S UD) are copied from (D, SU)
  //          Example: (S,UU) is copied from (U, SU).
  //      --- For Octet-Singlet isoscalar mixing, they should be placed
  //          (as shown below for eta-eta') at the beginning of the  D D
  //          and  S S  classes of mesons (not also  U U  because this is
  //          copied from  D D : see above).   
  //          We use the method   probabilityMixing(angleMix,order)
  //          where:  angleMix  is the mixing angle in degree;
  //                  order  is  1 for the first  (eta-like) and 
  //                             2 for the second (eta'-like)
  //          the return value of the method is the probability of the                    
  //          (|uubar>+|ddbar>)/sqrt(2) component. Therefore, in 
  //          D D  class we use 1/2 such probability (the other half
  //          is for  U U), and for  S S  class we use 1 - such probability. 
  //      --- The order in which the various hadrons of each class are
  //          filled in is not important: in particular, the lightest one
  //          can be placed anywhere, because its location is explicitly
  //          determined at the end by the method initialize().
  //

  // ---------------------- Mesons : 15 classes ---------------------------

  //---------------------- D D Mesons -------------------------------------
  _locHad[D][D].first = i;  // Just once at the beginning.
  
  //--- Octet-Singlet isoscalar mixing ---

  // *** eta - eta' ***                     --- First eta ---- 
  angleMix = etamix;                  // Mixing angle in degree for eta-eta'  
  order    = 1;                       // The first (eta)
  _vecHad[i].wt = 0.5 * CheckId::probabilityMixing(angleMix,order);
  _vecHad[i++].id = 221;  // eta  
  //                                        --- Second eta' ---
  order    = 2;                       // The second (eta')
  _vecHad[i].wt = 0.5 * CheckId::probabilityMixing(angleMix,order);
  _vecHad[i++].id = 331;  // eta'  

  // *** phi - omega ***                    --- First phi ---- 
  angleMix = phimix;                  // Mixing angle in degree for phi-omega  
  order    = 1;                       // The first (phi)
  _vecHad[i].wt = 0.5 * CheckId::probabilityMixing(angleMix,order);
  _vecHad[i++].id = 333;  // phi 
  //                                        --- Second omega ---
  order    = 2;                       // The second (omega)
  _vecHad[i].wt = 0.5 * CheckId::probabilityMixing(angleMix,order);
  _vecHad[i++].id = 223;  // omega  

  // *** h_1(1380) - h_1(1170) ***          --- First h_1(1380) ---- 
  angleMix = h1mix ;                  // Mixing angle in degree for h_1(1380)-h_1(1170)  
  order    = 1;                       // The first (h_1(1380))
  _vecHad[i].wt = 0.5 * CheckId::probabilityMixing(angleMix,order);
  _vecHad[i++].id = 10333;  // h_1(1380) 
  //                                        --- Second h_1(1170) ---
  order    = 2;                       // The second (h_1(1170))
  _vecHad[i].wt = 0.5 * CheckId::probabilityMixing(angleMix,order);
  _vecHad[i++].id = 10223;  // h_1(1170)  

  // *** missing - f_0(1370) *** 
  angleMix = f0mix;                   // Mixing angle in degree for missing-f_0(1370)  
  order    = 2;                       // The second (f_0(1370))
  _vecHad[i].wt = 0.5 * CheckId::probabilityMixing(angleMix,order);
  _vecHad[i++].id = 10221;  // f_0(1370)  

  // *** f_1(1420) - f_1(1285) ***          --- First f_1(1420) ---- 
  angleMix = f1mix;                   // Mixing angle in degree for f_1(1420)-f_1(1285)  
  order    = 1;                       // The first (f_1(1420))
  _vecHad[i].wt = 0.5 * CheckId::probabilityMixing(angleMix,order);
  _vecHad[i++].id = 20333;  // f_1(1420) 
  //                                        --- Second f_1(1285) ---
  order    = 2;                       // The second (f_1(1285))
  _vecHad[i].wt = 0.5 * CheckId::probabilityMixing(angleMix,order);
  _vecHad[i++].id = 20223;  // f_1(1285)  

  // *** f'_2 - f_2 ***                     --- First f'_2 ---- 
  angleMix = f2mix;                   // Mixing angle in degree for f'_2-f_2  
  order    = 1;                       // The first (f'_2)
  _vecHad[i].wt = 0.5 * CheckId::probabilityMixing(angleMix,order);
  _vecHad[i++].id = 335;  // f'_2 
  //                                        --- Second f_2 ---
  order    = 2;                       // The second (f_2)
  _vecHad[i].wt = 0.5 * CheckId::probabilityMixing(angleMix,order);
  _vecHad[i++].id = 225;  // f_2  

  // *** missing - omega(1650)  ***             <---  ***LOOKHERE*** NOT FOUND in Pythia7 
  // angleMix = omhmix;                  // Mixing angle in degree for missing-omega(1650)  
  // order    = 2;                       // The second (omega(1650))
  // _vecHad[i].wt = 0.5 * CheckId::probabilityMixing(angleMix,order);
  // _vecHad[i++].id = 30223;  // omega(1650)  

  // *** eta_2(1645) - eta_2(1870) ***          <---  ***LOOKHERE*** NOT FOUND in Pythia7 
  // angleMix = et2mix;                  // Mixing angle in degree for eta_2(1645)-eta_2(1870)  
  // order    = 1;                       // The first (eta_2(1645))
  // _vecHad[i].wt = 0.5 * CheckId::probabilityMixing(angleMix,order);
  // _vecHad[i++].id = 10225;  // eta_2(1645)  
  // //                                        --- Second eta_2(1870) ---
  // order    = 2;                       // The second (eta_2(1870))
  // _vecHad[i].wt = 0.5 * CheckId::probabilityMixing(angleMix,order);
  // _vecHad[i++].id = 10335;  // eta_2(1870)  

  // *** phi_3 - omega_3     ***                <---  ***LOOKHERE*** NOT FOUND in Pythia7 
  // angleMix = ph3mix;                  // Mixing angle in degree for phi_3-omega_3  
  // order    = 1;                       // The first (phi_3)
  // _vecHad[i].wt = 0.5 * CheckId::probabilityMixing(angleMix,order);
  // _vecHad[i++].id = 337;  // phi_3 
  // //                                        --- Second omega_3 ---
  // order    = 2;                       // The second (omega_3)
  // _vecHad[i].wt = 0.5 * CheckId::probabilityMixing(angleMix,order);
  // _vecHad[i++].id = 227;  // omega_3  

  //--- d dbar ,  u ubar ,  s sbar  mesons not included above ----
  //    Maybe there is something better to do than simply assumed 1/3
  //    for each of:  d dbar ,  u ubar ,  s sbar  as done below.
  //    (for example  idealAngleMix  could be assumed, but the problem is
  //     then to define the proper pair of partner scalars...)
  //    Similar thing for the  S S  class (see below).
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id = 9000221;  // f0(400-1200)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id = 9010221;  // f0(980)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id =  100221;  // eta(1295)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id =  100331;  // eta(1440)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id = 9020221;  // f0(1500)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id =   10331;  // f0(1710)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id =  200221;  // eta(1760) 
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id = 9030221;  // f0(2020)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id = 9040221;  // f0(2060)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id = 9050221;  // f0(2200)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id = 9060221;  // eta(2225)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id =  100223;  // omega(1420)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id = 9000223;  // f1(1510)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id =  100333;  // phi(1680)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id = 9000225;  // f2(1430) 
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id = 9010225;  // f2(1565) 
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id = 9020225;  // f2(1640)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id =  100225;  // f2(1810) 
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id = 9030225;  // f2(1950)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id =  100335;  // f2(2010)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id = 9040225;  // f2(2150)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id = 9050225;  // f2(2300)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id = 9060225;  // f2(2340)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id =     229;  // f4(2050)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id = 9000339;  // fj(2220) 
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id = 9000229;  // f4(2300)  

  //---  d dbar ,  u ubar  mesons ---

  _vecHad[i].wt = 0.5;   _vecHad[i++].id =     111;  // pi0 
  _vecHad[i].wt = 0.5;   _vecHad[i++].id = 9000111;  // a0(980)0 
  _vecHad[i].wt = 0.5;   _vecHad[i++].id =  100111;  // pi(1300)0 
  _vecHad[i].wt = 0.5;   _vecHad[i++].id =   10111;  // a0(1450)0 
  _vecHad[i].wt = 0.5;   _vecHad[i++].id =  200111;  // pi(1800)0 
  _vecHad[i].wt = 0.5;   _vecHad[i++].id =     113;  // rho(770)0 
  _vecHad[i].wt = 0.5;   _vecHad[i++].id =   10113;  // b1(1235)0
  _vecHad[i].wt = 0.5;   _vecHad[i++].id =   20113;  // a1(1260)0
  _vecHad[i].wt = 0.5;   _vecHad[i++].id = 9000113;  // pi1(1400)0 
  _vecHad[i].wt = 0.5;   _vecHad[i++].id =  100113;  // rho(1450)0 
  _vecHad[i].wt = 0.5;   _vecHad[i++].id = 9010113;  // pi1(1600)0 
  _vecHad[i].wt = 0.5;   _vecHad[i++].id = 9020113;  // a1(1640)0
  _vecHad[i].wt = 0.5;   _vecHad[i++].id =   30113;  // rho(1700)0 
  _vecHad[i].wt = 0.5;   _vecHad[i++].id = 9030113;  // rho(2150)0 
  _vecHad[i].wt = 0.5;   _vecHad[i++].id =     115;  // a2(1320)0 
  _vecHad[i].wt = 0.5;   _vecHad[i++].id = 9000115;  // a2(1660)0 
  _vecHad[i].wt = 0.5;   _vecHad[i++].id =   10115;  // pi2(1670)0
  _vecHad[i].wt = 0.5;   _vecHad[i++].id = 9010115;  // a2(1750)0
  _vecHad[i].wt = 0.5;   _vecHad[i++].id = 9020115;  // pi2(2100)0
  _vecHad[i].wt = 0.5;   _vecHad[i++].id =     117;  // rho3(1690)0
  _vecHad[i].wt = 0.5;   _vecHad[i++].id = 9000117;  // rho3(2250)0
  _vecHad[i].wt = 0.5;   _vecHad[i++].id =     119;  // a4(2040)0

  _locHad[D][D].last = i - 1;     // Just once at the end.

  //---------------------- D U Mesons -------------------------------------
  _locHad[D][U].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id =     211;  // pi+ 
  _vecHad[i++].id = 9000211;  // a0(980)+
  _vecHad[i++].id =  100211;  // pi(1300)+
  _vecHad[i++].id =   10211;  // a0(1450)+
  _vecHad[i++].id =  200211;  // pi(1800)+
  _vecHad[i++].id =     213;  // rho(770)+
  _vecHad[i++].id =   10213;  // b1(1235)+
  _vecHad[i++].id =   20213;  // a1(1260)+
  _vecHad[i++].id = 9000213;  // pi1(1400)+
  _vecHad[i++].id =  100213;  // rho(1450)+
  _vecHad[i++].id = 9010213;  // pi1(1600)+
  _vecHad[i++].id = 9020213;  // a1(1640)+
  _vecHad[i++].id =   30213;  // rho(1700)+
  _vecHad[i++].id = 9030213;  // rho(2150)+
  _vecHad[i++].id =     215;  // a2(1320)+
  _vecHad[i++].id = 9000215;  // a2(1660)+
  _vecHad[i++].id =   10215;  // pi2(1670)+
  _vecHad[i++].id = 9010215;  // a2(1750)+
  _vecHad[i++].id = 9020215;  // pi2(2100)+
  _vecHad[i++].id =     217;  // rho3(1690)+
  _vecHad[i++].id = 9000217;  // rho3(2250)+
  _vecHad[i++].id =     219;  // a4(2040)+

  _locHad[D][U].last = i - 1;     // Just once at the end.

  //---------------------- D S Mesons -------------------------------------
  _locHad[D][S].first = i;  // Just once at the beginning.
  
  // NOTE: neither K0L (id=130) nor K0S (id=310): only K0, K0bar matters
  //       as far as the cluster decays is concerned.

  _vecHad[i++].id =     311;  // K0
  _vecHad[i++].id =   10311;  // K0*(1430)0
  _vecHad[i++].id =  100311;  // K(1460)0 
  _vecHad[i++].id =  200311;  // K(1830)0 
  _vecHad[i++].id = 9000311;  // K0*(1950)0 
  _vecHad[i++].id =     313;  // K*(892)0 
  _vecHad[i++].id =   10313;  // K1(1270)0
  _vecHad[i++].id =   20313;  // K1(1400)0
  _vecHad[i++].id =  100313;  // K*(1410)0 
  _vecHad[i++].id = 9000313;  // K1(1650)0
  _vecHad[i++].id =   30313;  // K*(1680)0 
  _vecHad[i++].id =     315;  // K2*(1430)0 
  _vecHad[i++].id = 9000315;  // K2(1580)0 
  _vecHad[i++].id =   10315;  // K2(1770)0 
  _vecHad[i++].id =   20315;  // K2(1820)0 
  _vecHad[i++].id =  100315;  // K2*(1980)0 
  _vecHad[i++].id = 9010315;  // K2(2250)0 
  _vecHad[i++].id =     317;  // K3(1780)0 
  _vecHad[i++].id = 9010317;  // K1(2320)0 
  _vecHad[i++].id =     319;  // K4*(2045)0 
  _vecHad[i++].id = 9000319;  // K4(2500)0 

  _locHad[D][S].last = i - 1;     // Just once at the end.

  //---------------------- D C Mesons -------------------------------------
  _locHad[D][C].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id =     411;  // D+ 
  _vecHad[i++].id =   10411;  // D0*+
  _vecHad[i++].id =     413;  // D*(2010)+ 
  _vecHad[i++].id =   10413;  // D1(2420)+
  _vecHad[i++].id =   20413;  // D1(H)+   
  _vecHad[i++].id =     415;  // D2*(2460)+ 

  _locHad[D][C].last = i - 1;     // Just once at the end.

  //---------------------- D B Mesons -------------------------------------
  _locHad[D][B].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id =     511;  // B0
  _vecHad[i++].id =   10511;  // B0*0
  _vecHad[i++].id =     513;  // B*0
  _vecHad[i++].id =   10513;  // B1(L)0
  _vecHad[i++].id =   20513;  // B1(H)0
  _vecHad[i++].id =     515;  // B2*0

  _locHad[D][B].last = i - 1;     // Just once at the end.

  //---------------------- U U Mesons : Skip ------------------------------
  _locHad[U][U].first = _locHad[D][D].first;
  _locHad[U][U].last  = _locHad[D][D].last;

  //---------------------- U S Mesons -------------------------------------
  _locHad[U][S].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id =     321;  // K+
  _vecHad[i++].id =   10321;  // K0*(1430)+ 
  _vecHad[i++].id =  100321;  // K(1460)+
  _vecHad[i++].id =  200321;  // K(1830)+
  _vecHad[i++].id = 9000321;  // K0*(1950)+
  _vecHad[i++].id =     323;  // K*(892)+
  _vecHad[i++].id =   10323;  // K1(1270)+
  _vecHad[i++].id =   20323;  // K1(1400)+
  _vecHad[i++].id =  100323;  // K*(1410)+
  _vecHad[i++].id = 9000323;  // K1(1650)+
  _vecHad[i++].id =   30323;  // K*(1680)+
  _vecHad[i++].id =     325;  // K2*(1430)+
  _vecHad[i++].id = 9000325;  // K2(1580)+
  _vecHad[i++].id =   10325;  // K2(1770)+
  _vecHad[i++].id =   20325;  // K2(1820)+
  _vecHad[i++].id =  100325;  // K2*(1980)+
  _vecHad[i++].id = 9010325;  // K2(2250)+
  _vecHad[i++].id =     327;  // K3(1780)+
  _vecHad[i++].id = 9010327;  // K1(2320)+
  _vecHad[i++].id =     329;  // K4*(2045)+
  _vecHad[i++].id = 9000329;  // K4(2500)+

  _locHad[U][S].last = i - 1;     // Just once at the end.

  //---------------------- U C Mesons -------------------------------------
  _locHad[U][C].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id =     421;  // D0 
  _vecHad[i++].id =   10421;  // D0*0
  _vecHad[i++].id =     423;  // D*(2007)0 
  _vecHad[i++].id =   10423;  // D1(2420)0
  _vecHad[i++].id =   20423;  // D1(H)0   
  _vecHad[i++].id =     425;  // D2*(2460)0 

  _locHad[U][C].last = i - 1;     // Just once at the end.

  //---------------------- U B Mesons -------------------------------------
  _locHad[U][B].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id =     521;  // B+
  _vecHad[i++].id =   10521;  // B0*+ 
  _vecHad[i++].id =     523;  // B*+
  _vecHad[i++].id =   10523;  // B1(L)+
  _vecHad[i++].id =   20523;  // B1(H)+
  _vecHad[i++].id =     525;  // B2*+

  _locHad[U][B].last = i - 1;     // Just once at the end.

  //---------------------- S S Mesons -------------------------------------
  _locHad[S][S].first = i;  // Just once at the beginning.
  
  //--- Octet-Singlet isoscalar mixing ---

  // *** eta - eta' ***                     --- First eta ---- 
  angleMix = etamix;                  // Mixing angle in degree for eta-eta'  
  order    = 1;                       // The first (eta)
  _vecHad[i].wt = 1.0 - CheckId::probabilityMixing(angleMix,order);
  _vecHad[i++].id = 221;  // eta  
  //                                        --- Second eta' ---
  order    = 2;                       // The second (eta')
  _vecHad[i].wt = 1.0 - CheckId::probabilityMixing(angleMix,order);
  _vecHad[i++].id = 331;  // eta'      

  // *** phi - omega ***                    --- First phi ---- 
  angleMix = phimix;                  // Mixing angle in degree for phi-omega  
  order    = 1;                       // The first (phi)
  _vecHad[i].wt = 1.0 - CheckId::probabilityMixing(angleMix,order);
  _vecHad[i++].id = 333;  // phi  
  //                                        --- Second omega ---
  order    = 2;                       // The second (omega)
  _vecHad[i].wt = 1.0 - CheckId::probabilityMixing(angleMix,order);
  _vecHad[i++].id = 223;  // omega      

  // *** h_1(1380) - h_1(1170) ***          --- First h_1(1380) ---- 
  angleMix = h1mix;                   // Mixing angle in degree for h_1(1380)-h_1(1170)  
  order    = 1;                       // The first (h_1(1380))
  _vecHad[i].wt = 1.0 - CheckId::probabilityMixing(angleMix,order);
  _vecHad[i++].id = 10333;  // h_1(1380)  
  //                                        --- Second h_1(1170) ---
  order    = 2;                       // The second (h_1(1170))
  _vecHad[i].wt = 1.0 - CheckId::probabilityMixing(angleMix,order);
  _vecHad[i++].id = 10223;  // h_1(1170)      

  // *** missing - f_0(1370) *** 
  angleMix = f0mix;                   // Mixing angle in degree for phi-f_0(1370)  
  order    = 2;                       // The second (f_0(1370))
  _vecHad[i].wt = 1.0 - CheckId::probabilityMixing(angleMix,order);
  _vecHad[i++].id = 10221;  // f_0(1370)      

  // *** f_1(1420) - f_1(1285) ***          --- First f_1(1420) ---- 
  angleMix = f1mix;                   // Mixing angle in degree for f_1(1420)-f_1(1285)  
  order    = 1;                       // The first (f_1(1420))
  _vecHad[i].wt = 1.0 - CheckId::probabilityMixing(angleMix,order);
  _vecHad[i++].id = 20333;  // f_1(1420)  
  //                                        --- Second f_1(1285) ---
  order    = 2;                       // The second (f_1(1285))
  _vecHad[i].wt = 1.0 - CheckId::probabilityMixing(angleMix,order);
  _vecHad[i++].id = 20223;  // f_1(1285)      

  // *** f'_2 - f_2 ***                     --- First f'_2 ---- 
  angleMix = f2mix;                   // Mixing angle in degree for f'_2-f_2  
  order    = 1;                       // The first (f'_2)
  _vecHad[i].wt = 1.0 - CheckId::probabilityMixing(angleMix,order);
  _vecHad[i++].id = 335;  // f'_2  
  //                                        --- Second f_2 ---
  order    = 2;                       // The second (f_2)
  _vecHad[i].wt = 1.0 - CheckId::probabilityMixing(angleMix,order);
  _vecHad[i++].id = 225;  // f_2      

  // *** missing - omega(1650)  ***             <---  ***LOOKHERE*** NOT FOUND in Pythia7 
  // angleMix = omhmix;                  // Mixing angle in degree for missing-omega(1650)  
  // order    = 2;                       // The second (omega(1650))
  // _vecHad[i].wt = 1.0 - CheckId::probabilityMixing(angleMix,order);
  // _vecHad[i++].id = 30223;  // omega(1650)  

  // *** eta_2(1645) - eta_2(1870) ***          <---  ***LOOKHERE*** NOT FOUND in Pythia7 
  // angleMix = et2mix;                  // Mixing angle in degree for eta_2(1645)-eta_2(1870)  
  // order    = 1;                       // The first (eta_2(1645))
  // _vecHad[i].wt = 1.0 - CheckId::probabilityMixing(angleMix,order);
  // _vecHad[i++].id = 10225;  // eta_2(1645)  
  // //                                        --- Second eta_2(1870) ---
  // order    = 2;                       // The second (eta_2(1870))
  // _vecHad[i].wt = 1.0 - CheckId::probabilityMixing(angleMix,order);
  // _vecHad[i++].id = 10335;  // eta_2(1870)  

  // *** phi_3 - omega_3     ***                <---  ***LOOKHERE*** NOT FOUND in Pythia7 
  // angleMix = ph3mix;                  // Mixing angle in degree for phi_3-omega_3  
  // order    = 1;                       // The first (phi_3)
  // _vecHad[i].wt = 1.0 - CheckId::probabilityMixing(angleMix,order);
  // _vecHad[i++].id = 337;  // phi_3 
  // //                                        --- Second omega_3 ---
  // order    = 2;                       // The second (omega_3)
  // _vecHad[i].wt = 1.0 - CheckId::probabilityMixing(angleMix,order);
  // _vecHad[i++].id = 227;  // omega_3  

  //--- d dbar ,  u ubar ,  s sbar  mesons not included above ----
  //    Maybe there is something better to do than simply assumed 1/3
  //    for each of:  d dbar ,  u ubar ,  s sbar  as done below.
  //    (for example  idealAngleMix  could be assumed, but the problem is
  //     then to define the proper pair of partner scalars...)
  //    Similar thing for the  U U  class (see above).
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id = 9000221;  // f0(400-1200)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id = 9010221;  // f0(980)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id =  100221;  // eta(1295)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id =  100331;  // eta(1440)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id = 9020221;  // f0(1500)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id =   10331;  // f0(1710)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id =  200221;  // eta(1760) 
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id = 9030221;  // f0(2020)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id = 9040221;  // f0(2060)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id = 9050221;  // f0(2200)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id = 9060221;  // eta(2225)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id =  100223;  // omega(1420)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id = 9000223;  // f1(1510)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id =  100333;  // phi(1680)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id = 9000225;  // f2(1430) 
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id = 9010225;  // f2(1565) 
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id = 9020225;  // f2(1640)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id =  100225;  // f2(1810) 
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id = 9030225;  // f2(1950)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id =  100335;  // f2(2010)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id = 9040225;  // f2(2150)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id = 9050225;  // f2(2300)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id = 9060225;  // f2(2340)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id =     229;  // f4(2050)  
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id = 9000339;  // fj(2220) 
  // _vecHad[i].wt = 0.333333;   _vecHad[i++].id = 9000229;  // f4(2300)  

  //--- pure s sbar not mixed mesons  ---

  // _vecHad[i++].id = xyz;  // ??? 

  _locHad[S][S].last = i - 1;     // Just once at the end.

  //---------------------- S C Mesons -------------------------------------
  _locHad[S][C].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id =     431;  // Ds+ 
  _vecHad[i++].id =   10431;  // Ds0*+
  _vecHad[i++].id =     433;  // Ds*+ 
  _vecHad[i++].id =   10433;  // Ds1(2536)+
  _vecHad[i++].id =   20433;  // Ds1(H)+
  _vecHad[i++].id =     435;  // Ds2*+ 

  _locHad[S][C].last = i - 1;     // Just once at the end.

  //---------------------- S B Mesons -------------------------------------
  _locHad[S][B].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id =     531;  // Bs0
  _vecHad[i++].id =   10531;  // Bs0*0
  _vecHad[i++].id =     533;  // Bs*0
  _vecHad[i++].id =   10533;  // Bs1(L)0
  _vecHad[i++].id =   20533;  // Bs1(H)0 
  _vecHad[i++].id =     535;  // Bs2*0

  _locHad[S][B].last = i - 1;     // Just once at the end.

  //---------------------- C C Mesons -------------------------------------
  _locHad[C][C].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id =     441;  // eta_c(1S) 
  _vecHad[i++].id =   10441;  // chi_c0(1P)
  _vecHad[i++].id =  100441;  // eta_c(2S) 
  _vecHad[i++].id =     443;  // J/psi(1S) 
  _vecHad[i++].id =   10443;  // h_c(1P)
  _vecHad[i++].id =   20443;  // chi_c1(1P)
  _vecHad[i++].id =  100443;  // psi(2S) 
  _vecHad[i++].id =   30443;  // psi(3770)
  _vecHad[i++].id = 9000443;  // psi(4040) 
  _vecHad[i++].id = 9010443;  // psi(4160) 
  _vecHad[i++].id = 9020443;  // psi(4415) 
  _vecHad[i++].id =     445;  // chi_c2(1P) 
  _vecHad[i++].id = 9000445;  // psi(3836) 

  _locHad[C][C].last = i - 1;     // Just once at the end.

  //---------------------- C B Mesons -------------------------------------
  _locHad[C][B].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id =     541;  // Bc+
  _vecHad[i++].id =   10541;  // Bc0*+
  _vecHad[i++].id =     543;  // Bc*+
  _vecHad[i++].id =   10543;  // Bc1(L)+
  _vecHad[i++].id =   20543;  // Bc1(H)+
  _vecHad[i++].id =     545;  // Bc2*+

  _locHad[C][B].last = i - 1;     // Just once at the end.

  //---------------------- B B Mesons -------------------------------------
  _locHad[B][B].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id =     551;  // eta_b(1S) 
  _vecHad[i++].id =   10551;  // chi_b0(1P) 
  _vecHad[i++].id =  100551;  // eta_b(2S)
  _vecHad[i++].id =  110551;  // chi_b0(2P) 
  _vecHad[i++].id =  200551;  // eta_b(3S)
  _vecHad[i++].id =  210551;  // chi_b0(3P) 
  _vecHad[i++].id =     553;  // upsilon(1S)
  _vecHad[i++].id =   10553;  // h_b(1P) 
  _vecHad[i++].id =   20553;  // chi_b1(1P) 
  _vecHad[i++].id =   30553;  // upsilon1(1D) 
  _vecHad[i++].id =  100553;  // upsilon(2S) 
  _vecHad[i++].id =  110553;  // h_b(2P)
  _vecHad[i++].id =  120553;  // chi_b1(2P) 
  _vecHad[i++].id =  130553;  // upsilon1(2D) 
  _vecHad[i++].id =  200553;  // upsilon(3S) 
  _vecHad[i++].id =  210553;  // h_b(3P) 
  _vecHad[i++].id =  220553;  // chi_b1(3P) 
  _vecHad[i++].id =  300553;  // upsilon(4S) 
  _vecHad[i++].id = 9000553;  // upsilon(10860) 
  _vecHad[i++].id = 9010553;  // upsilon(11020)
  _vecHad[i++].id =     555;  // chi_b2(1P)
  _vecHad[i++].id =   10555;  // eta_b2(1D) 
  _vecHad[i++].id =   20555;  // upsilon2(1D)
  _vecHad[i++].id =  100555;  // chi_b2(2P) 
  _vecHad[i++].id =  110555;  // eta_b2(2D) 
  _vecHad[i++].id =  120555;  // upsilon2(2D) 
  _vecHad[i++].id =  200555;  // chi_b2(3P)
  _vecHad[i++].id =     557;  // upsilon3(1D)
  _vecHad[i++].id =  100557;  // upsilon3(2D) 

  _locHad[B][B].last = i - 1;     // Just once at the end.


  // ---------------------- Baryons : 75 classes ----------------------------
  // No excited baryons!

  //---------------------- D DD Baryons -------------------------------------
  _locHad[D][DD].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id = 1114;  // Delta- 

  _locHad[D][DD].last = i - 1;     // Just once at the end.

  //---------------------- D UD Baryons -------------------------------------
  _locHad[D][UD].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id = 2112;  // n0 
  _vecHad[i++].id = 2114;  // Delta0 

  _locHad[D][UD].last = i - 1;     // Just once at the end.

  //---------------------- D UU Baryons -------------------------------------
  _locHad[D][UU].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id = 2212;  // p+ 
  _vecHad[i++].id = 2214;  // Delta+ 

  _locHad[D][UU].last = i - 1;     // Just once at the end.

  //---------------------- D SD Baryons -------------------------------------
  _locHad[D][SD].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id = 3112;  // Sigma-
  _vecHad[i++].id = 3114;  // Sigma*-

  _locHad[D][SD].last = i - 1;     // Just once at the end.

  //---------------------- D SU Baryons -------------------------------------
  _locHad[D][SU].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id = 3122;  // Lambda0 
  _vecHad[i++].id = 3212;  // Sigma0 
  _vecHad[i++].id = 3214;  // Sigma*0 

  _locHad[D][SU].last = i - 1;     // Just once at the end.

  //---------------------- D SS Baryons -------------------------------------
  _locHad[D][SS].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id = 3312;  // Xi- 
  _vecHad[i++].id = 3314;  // Xi*- 

  _locHad[D][SS].last = i - 1;     // Just once at the end.

  //---------------------- D CD Baryons -------------------------------------
  _locHad[D][CD].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id = 4112;  // Sigma_c0
  _vecHad[i++].id = 4114;  // Sigma_c*0

  _locHad[D][CD].last = i - 1;     // Just once at the end.

  //---------------------- D CU Baryons  -------------------------------
  _locHad[D][CU].first = i;  // Just once at the beginning. 
  
  _vecHad[i++].id = 4122;  // Lambda_c+ 
  _vecHad[i++].id = 4212;  // Sigma_c+ 
  _vecHad[i++].id = 4214;  // Sigma_c*+ 

  _locHad[D][CU].last = i - 1;     // Just once at the end.

  //---------------------- D CS Baryons -------------------------------------
  _locHad[D][CS].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id = 4132;  // Xi_c0 
  _vecHad[i++].id = 4312;  // Xi_c'0 
  _vecHad[i++].id = 4314;  // Xi_c*0 

  _locHad[D][CS].last = i - 1;     // Just once at the end.

  //---------------------- D CC Baryons -------------------------------------
  _locHad[D][CC].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id = 4412;  // Xi_cc+
  _vecHad[i++].id = 4414;  // Xi_cc*+

  _locHad[D][CC].last = i - 1;     // Just once at the end.

  //---------------------- D BD Baryons -------------------------------------
  _locHad[D][BD].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id = 5112;  // Sigma_b- 
  _vecHad[i++].id = 5114;  // Sigma_b*- 

  _locHad[D][BD].last = i - 1;     // Just once at the end.

  //---------------------- D BU Baryons -------------------------------------
  _locHad[D][BU].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id = 5122;  // Lambda_b0 
  _vecHad[i++].id = 5212;  // Sigma_b0 
  _vecHad[i++].id = 5214;  // Sigma_b*0 

  _locHad[D][BU].last = i - 1;     // Just once at the end.

  //---------------------- D BS Baryons -------------------------------------
  _locHad[D][BS].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id = 5132;  // Xi_b- 
  _vecHad[i++].id = 5312;  // Xi_b'- 
  _vecHad[i++].id = 5314;  // Xi_b*- 

  _locHad[D][BS].last = i - 1;     // Just once at the end.

  //---------------------- D BC Baryons -------------------------------------
  _locHad[D][BC].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id = 5142;  // Xi_bc0
  _vecHad[i++].id = 5412;  // Xi_bc'0
  _vecHad[i++].id = 5414;  // Xi_bc*0

  _locHad[D][BC].last = i - 1;     // Just once at the end.

  //---------------------- D BB Baryons -------------------------------------
  _locHad[D][BB].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id = 5512;  // Xi_bb-
  _vecHad[i++].id = 5514;  // Xi_bb*-

  _locHad[D][BB].last = i - 1;     // Just once at the end.

  //---------------------- U DD Baryons : Skip ------------------------------
  _locHad[U][DD].first = _locHad[D][UD].first;
  _locHad[U][DD].last  = _locHad[D][UD].last;

  //---------------------- U UD Baryons : Skip ------------------------------
  _locHad[U][UD].first = _locHad[D][UU].first;
  _locHad[U][UD].last  = _locHad[D][UU].last;

  //---------------------- U UU Baryons -------------------------------------
  _locHad[U][UU].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id = 2224;  // Delta++ 

  _locHad[U][UU].last = i - 1;     // Just once at the end.

  //---------------------- U SD Baryons : Skip ------------------------------
  _locHad[U][SD].first = _locHad[D][SU].first; 
  _locHad[U][SD].last  = _locHad[D][SU].last;

  //---------------------- U SU Baryons -------------------------------------
  _locHad[U][SU].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id = 3222;  // Sigma+ 
  _vecHad[i++].id = 3224;  // Sigma*+ 

  _locHad[U][SU].last = i - 1;     // Just once at the end.

  //---------------------- U SS Baryons -------------------------------------
  _locHad[U][SS].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id = 3322;  // Xi0 
  _vecHad[i++].id = 3324;  // Xi*0 

  _locHad[U][SS].last = i - 1;     // Just once at the end.

  //---------------------- U CD Baryons : Skip ------------------------------
  _locHad[U][CD].first = _locHad[D][CU].first; 
  _locHad[U][CD].last  = _locHad[D][CU].last;

  //---------------------- U CU Baryons -------------------------------------
  _locHad[U][CU].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id = 4222;  // Sigma_c++ 
  _vecHad[i++].id = 4224;  // Sigma_c*++ 

  _locHad[U][CU].last = i - 1;     // Just once at the end.

  //---------------------- U CS Baryons -------------------------------------
  _locHad[U][CS].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id = 4232;  // Xi_c+ 
  _vecHad[i++].id = 4322;  // Xi_c'+ 
  _vecHad[i++].id = 4324;  // Xi_c*+ 

  _locHad[U][CS].last = i - 1;     // Just once at the end.

  //---------------------- U CC Baryons -------------------------------------
  _locHad[U][CC].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id = 4422;  // Xi_cc++
  _vecHad[i++].id = 4424;  // Xi_cc*++

  _locHad[U][CC].last = i - 1;     // Just once at the end.

  //---------------------- U BD Baryons : Skip ------------------------------
  _locHad[U][BD].first = _locHad[D][BU].first; 
  _locHad[U][BD].last  = _locHad[D][BU].last;

  //---------------------- U BU Baryons -------------------------------------
  _locHad[U][BU].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id = 5222;  // Sigma_b+ 
  _vecHad[i++].id = 5224;  // Sigma_b*+ 

  _locHad[U][BU].last = i - 1;     // Just once at the end.

  //---------------------- U BS Baryons -------------------------------------
  _locHad[U][BS].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id = 5232;  // Xi_b0 
  _vecHad[i++].id = 5322;  // Xi_b'0 
  _vecHad[i++].id = 5324;  // Xi_b*0 

  _locHad[U][BS].last = i - 1;     // Just once at the end.

  //---------------------- U BC Baryons -------------------------------------
  _locHad[U][BC].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id = 5242;  // Xi_bc+
  _vecHad[i++].id = 5422;  // Xi_bc'+
  _vecHad[i++].id = 5424;  // Xi_bc*+

  _locHad[U][BC].last = i - 1;     // Just once at the end.

  //---------------------- U BB Baryons -------------------------------------
  _locHad[U][BB].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id = 5522;  // Xi_bb0
  _vecHad[i++].id = 5524;  // Xi_bb*0

  _locHad[U][BB].last = i - 1;     // Just once at the end.

  //---------------------- S DD Baryons : Skip ------------------------------
  _locHad[S][DD].first = _locHad[D][SD].first;
  _locHad[S][DD].last  = _locHad[D][SD].last;

  //---------------------- S UD Baryons : Skip ------------------------------
  _locHad[S][UD].first = _locHad[D][SU].first;
  _locHad[S][UD].last  = _locHad[D][SU].last;

  //---------------------- S UU Baryons : Skip ------------------------------
  _locHad[S][UU].first = _locHad[U][SU].first; 
  _locHad[S][UU].last  = _locHad[U][SU].last;

  //---------------------- S SD Baryons : Skip ------------------------------
  _locHad[S][SD].first = _locHad[D][SS].first; 
  _locHad[S][SD].last  = _locHad[D][SS].last;

  //---------------------- S SU Baryons : Skip ------------------------------
  _locHad[S][SU].first = _locHad[U][SS].first; 
  _locHad[S][SU].last  = _locHad[U][SS].last;

  //---------------------- S SS Baryons -------------------------------------
  _locHad[S][SS].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id = 3334;  // Omega- 

  _locHad[S][SS].last = i - 1;     // Just once at the end.

  //---------------------- S CD Baryons : Skip ------------------------------
  _locHad[S][CD].first = _locHad[D][CS].first; 
  _locHad[S][CD].last  = _locHad[D][CS].last;

  //---------------------- S CU Baryons : Skip ------------------------------
  _locHad[S][CU].first = _locHad[U][CS].first; 
  _locHad[S][CU].last  = _locHad[U][CS].last;

  //---------------------- S CS Baryons -------------------------------------
  _locHad[S][CS].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id = 4332;  // Omega_c0 
  _vecHad[i++].id = 4334;  // Omega_c*0 

  _locHad[S][CS].last = i - 1;     // Just once at the end.

  //---------------------- S CC Baryons -------------------------------------
  _locHad[S][CC].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id = 4432;  // Omega_cc+
  _vecHad[i++].id = 4434;  // Omega_cc*+

  _locHad[S][CC].last = i - 1;     // Just once at the end.

  //---------------------- S BD Baryons : Skip ------------------------------
  _locHad[S][BD].first = _locHad[D][BS].first; 
  _locHad[S][BD].last  = _locHad[D][BS].last;

  //---------------------- S BU Baryons : Skip ------------------------------
  _locHad[S][BU].first = _locHad[U][BS].first; 
  _locHad[S][BU].last  = _locHad[U][BS].last;

  //---------------------- S BS Baryons -------------------------------------
  _locHad[S][BS].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id = 5332;  // Omega_b- 
  _vecHad[i++].id = 5334;  // Omega_b*- 

  _locHad[S][BS].last = i - 1;     // Just once at the end.

  //---------------------- S BC Baryons -------------------------------------
  _locHad[S][BC].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id = 5342;  // Omega_bc0
  _vecHad[i++].id = 5432;  // Omega_bc'0
  _vecHad[i++].id = 5434;  // Omega_bc*0

  _locHad[S][BC].last = i - 1;     // Just once at the end.

  //---------------------- S BB Baryons -------------------------------------
  _locHad[S][BB].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id = 5532;  // Omega_bb-
  _vecHad[i++].id = 5534;  // Omega_bb*-

  _locHad[S][BB].last = i - 1;     // Just once at the end.

  //---------------------- C DD Baryons : Skip ------------------------------
  _locHad[C][DD].first = _locHad[D][CD].first; 
  _locHad[C][DD].last  = _locHad[D][CD].last;

  //---------------------- C UD Baryons : Skip ------------------------------
  _locHad[C][UD].first = _locHad[D][CU].first; 
  _locHad[C][UD].last  = _locHad[D][CU].last;

  //---------------------- C UU Baryons : Skip ------------------------------
  _locHad[C][UU].first = _locHad[U][CU].first; 
  _locHad[C][UU].last  = _locHad[U][CU].last;

  //---------------------- C SD Baryons : Skip ------------------------------
  _locHad[C][SD].first = _locHad[D][CS].first; 
  _locHad[C][SD].last  = _locHad[D][CS].last;

  //---------------------- C SU Baryons : Skip ------------------------------
  _locHad[C][SU].first = _locHad[U][CS].first; 
  _locHad[C][SU].last  = _locHad[U][CS].last;

  //---------------------- C SS Baryons : Skip ------------------------------
  _locHad[C][SS].first = _locHad[S][CS].first; 
  _locHad[C][SS].last  = _locHad[S][CS].last;

  //---------------------- C CD Baryons : Skip ------------------------------
  _locHad[C][CD].first = _locHad[D][CC].first; 
  _locHad[C][CD].last  = _locHad[D][CC].last;

  //---------------------- C CU Baryons : Skip ------------------------------
  _locHad[C][CU].first = _locHad[U][CC].first; 
  _locHad[C][CU].last  = _locHad[U][CC].last;

  //---------------------- C CS Baryons : Skip ------------------------------
  _locHad[C][CS].first = _locHad[S][CC].first; 
  _locHad[C][CS].last  = _locHad[S][CC].last;

  //---------------------- C CC Baryons -------------------------------------
  _locHad[C][CC].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id = 4444;  // Omega_ccc++

  _locHad[C][CC].last = i - 1;     // Just once at the end.

  //---------------------- C BD Baryons : Skip ------------------------------
  _locHad[C][BD].first = _locHad[D][BC].first; 
  _locHad[C][BD].last  = _locHad[D][BC].last;

  //---------------------- C BU Baryons : Skip ------------------------------
  _locHad[C][BU].first = _locHad[U][BC].first; 
  _locHad[C][BU].last  = _locHad[U][BC].last;

  //---------------------- C BS Baryons : Skip ------------------------------
  _locHad[C][BS].first = _locHad[S][BC].first; 
  _locHad[C][BS].last  = _locHad[S][BC].last;

  //---------------------- C BC Baryons -------------------------------------
  _locHad[C][BC].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id = 5442;  // Omega_bcc+
  _vecHad[i++].id = 5444;  // Omega_bcc*+

  _locHad[C][BC].last = i - 1;     // Just once at the end.

  //---------------------- C BB Baryons -------------------------------------
  _locHad[C][BB].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id = 5542;  // Omega_bbc0 
  _vecHad[i++].id = 5544;  // Omega_bbc*0 

  _locHad[C][BB].last = i - 1;     // Just once at the end.

  //---------------------- B DD Baryons : Skip ------------------------------
  _locHad[B][DD].first = _locHad[D][BD].first; 
  _locHad[B][DD].last  = _locHad[D][BD].last;

  //---------------------- B UD Baryons : Skip ------------------------------
  _locHad[B][UD].first = _locHad[D][BU].first; 
  _locHad[B][UD].last  = _locHad[D][BU].last;

  //---------------------- B UU Baryons : Skip ------------------------------
  _locHad[B][UU].first = _locHad[U][BU].first; 
  _locHad[B][UU].last  = _locHad[U][BU].last;

  //---------------------- B SD Baryons : Skip ------------------------------
  _locHad[B][SD].first = _locHad[D][BS].first; 
  _locHad[B][SD].last  = _locHad[D][BS].last;

  //---------------------- B SU Baryons : Skip ------------------------------
  _locHad[B][SU].first = _locHad[U][BS].first; 
  _locHad[B][SU].last  = _locHad[U][BS].last;

  //---------------------- B SS Baryons : Skip ------------------------------
  _locHad[B][SS].first = _locHad[S][BS].first; 
  _locHad[B][SS].last  = _locHad[S][BS].last;

  //---------------------- B CD Baryons : Skip ------------------------------
  _locHad[B][CD].first = _locHad[D][BC].first; 
  _locHad[B][CD].last  = _locHad[D][BC].last;

  //---------------------- B CU Baryons : Skip ------------------------------
  _locHad[B][CU].first = _locHad[U][BC].first; 
  _locHad[B][CU].last  = _locHad[U][BC].last;

  //---------------------- B CS Baryons : Skip ------------------------------
  _locHad[B][CS].first = _locHad[S][BC].first; 
  _locHad[B][CS].last  = _locHad[S][BC].last;

  //---------------------- B CC Baryons : Skip ------------------------------
  _locHad[B][CC].first = _locHad[C][BC].first; 
  _locHad[B][CC].last  = _locHad[C][BC].last;

  //---------------------- B BD Baryons : Skip ------------------------------
  _locHad[B][BD].first = _locHad[D][BB].first; 
  _locHad[B][BD].last  = _locHad[D][BB].last;

  //---------------------- B BU Baryons : Skip ------------------------------
  _locHad[B][BU].first = _locHad[U][BB].first; 
  _locHad[B][BU].last  = _locHad[U][BB].last;

  //---------------------- B BS Baryons : Skip ------------------------------
  _locHad[B][BS].first = _locHad[S][BB].first; 
  _locHad[B][BS].last  = _locHad[S][BB].last;

  //---------------------- B BC Baryons : Skip ------------------------------
  _locHad[B][BC].first = _locHad[C][BB].first; 
  _locHad[B][BC].last  = _locHad[C][BB].last;

  //---------------------- B BB Baryons -------------------------------------
  _locHad[B][BB].first = i;  // Just once at the beginning.
  
  _vecHad[i++].id = 5554;  // Omega_bbb-

  _locHad[B][BB].last = i - 1;     // Just once at the end.


  //===================== End filling table of hadrons ======================

  i--;  // Now i is really the number of hadrons.

  return i;

}


// -------------------------------------------------------
// ------------- END OF FILE  ----------------------------
// -------------------------------------------------------
