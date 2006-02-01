

// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HadronsSelector class.
//

#include "HadronSelector.h"
#include <ThePEG/Interface/ClassDocumentation.h>
#include <ThePEG/Interface/Parameter.h>
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/PDT/ParticleData.h>
#include <ThePEG/PDT/EnumParticles.h>
#include <ThePEG/Repository/EventGenerator.h>
#include "Herwig++/Utilities/HwDebug.h"
#include "Herwig++/Utilities/CheckId.h"
#include "Herwig++/Utilities/Kinematics.h"

using namespace Herwig;

HadronSelector::~HadronSelector() {}


void HadronSelector::persistentOutput(PersistentOStream & os) const {
  os << _PwtDquark << _PwtUquark << _PwtSquark 
     << _PwtCquark << _PwtBquark << _PwtDIquark
     << _SngWt << _DecWt << _Pwt << _ClusterDKMode; 
  for(int i = D; i<=B; i++)
    for(int j = i; j<BB; j++)
      os << _table[i][j];
}


void HadronSelector::persistentInput(PersistentIStream & is, int) {
  is >> _PwtDquark >> _PwtUquark >> _PwtSquark 
     >> _PwtCquark >> _PwtBquark >> _PwtDIquark
     >> _SngWt >> _DecWt >> _Pwt >> _ClusterDKMode;
  for(int i = D; i<=B; i++) {
    for(int j = i; j<BB; j++) {
      is >> _table[i][j];
      if(i != j) _table[j][i] = _table[i][j];
    }
  }
}


ClassDescription<HadronSelector> HadronSelector::initHadronSelector;
// Definition of the static class description member.


void HadronSelector::Init() {

  static ClassDocumentation<HadronSelector> documentation
    ("This class selects a proper pair (or single) hadrons...");

  static Parameter<HadronSelector,double>
    interfacePwtDquark("PwtDquark","Weight for choosing a quark D",
		       &HadronSelector::_PwtDquark, 0, 1.0, 0.0, 10.0,false,false,false);
  static Parameter<HadronSelector,double>
    interfacePwtUquark("PwtUquark","Weight for choosing a quark U",
		       &HadronSelector::_PwtUquark, 0, 1.0, 0.0, 10.0,false,false,false);
  static Parameter<HadronSelector,double>
    interfacePwtSquark("PwtSquark","Weight for choosing a quark S",
		       &HadronSelector::_PwtSquark, 0, 1.0, 0.0, 10.0,false,false,false);
  static Parameter<HadronSelector,double>
    interfacePwtCquark("PwtCquark","Weight for choosing a quark C",
		       &HadronSelector::_PwtCquark, 0, 1.0, 0.0, 10.0,false,false,false);
  static Parameter<HadronSelector,double>
    interfacePwtBquark("PwtBquark","Weight for choosing a quark B",
		       &HadronSelector::_PwtBquark, 0, 1.0, 0.0, 10.0,false,false,false);
  static Parameter<HadronSelector,double>
    interfacePwtDIquark("PwtDIquark","Weight for choosing a DIquark",
			&HadronSelector::_PwtDIquark, 0, 1.0, 0.0, 100.0,false,false,false);
  static Parameter<HadronSelector,double>
    interfaceSngWt("SngWt","Weight for singlet baryons",
                  &HadronSelector::_SngWt, 0, 1.0, 0.0, 10.0,false,false,false);
  static Parameter<HadronSelector,double>
    interfaceDecWt("DecWt","Weight for decuplet baryons",
                  &HadronSelector::_DecWt, 0, 1.0, 0.0, 10.0,false,false,false);

  static Parameter<HadronSelector,int>
    interfaceCDKM("DKMode","Decay mode for Clusters...hw64=0,kupco,default",
                  &HadronSelector::_ClusterDKMode, 0, 0, 0, 2,false,false,false);

  static Parameter<HadronSelector,int>
    interfaceTrials("Trial","Debug param for trials...pions=1,spin<3,<4",
		  &HadronSelector::_trial, 0, 0, 0, 5,false,false,false);
}


// ------------- PUBLIC MAIN METHODS ----------------------------

long HadronSelector::
lightestHadron(const long id1, const long id2, const long id3) const {
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Hadronization ) {
    safetyCheck(id1,id2,id3);
  }
  long lightest = 0;
  if (id3 == 0) {    // The method assumes id3 == 0   
    int flav1 = convertIdToFlavour(id1);
    int flav2 = convertIdToFlavour(id2);
    if ( flav1 != -1  &&  flav2 != -1 ) {
      lightest = _table[flav1][flav2].begin()->id;
      int sign = signHadron(id1,id2,lightest);
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Hadronization ) {  // Safety check, usually skipped.  
	if ( lightest == 0  ||  sign == 0 ) {  
	  generator()->logWarning( Exception("HadronSelector::lightestHadron "
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
      generator()->logWarning( Exception("HadronSelector::lightestHadron "
					 "***-1 flav1 or flav2*** ",
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


Energy HadronSelector::
massLightestHadron(const long id1, const long id2, const long id3) const {
  Energy mass=-123456789.;
  if (id3 == 0) {    // The method assumes id3 == 0   
    int flav1 = convertIdToFlavour(id1);
    int flav2 = convertIdToFlavour(id2);
    if ( flav1 != -1  &&  flav2 != -1 ) {
      mass = _table[flav1][flav2].begin()->mass;
    } else {
      generator()->logWarning( Exception("HadronSelector::massLightestHadron "
					 "***Zero flav1 or flav2*** ",
					 Exception::warning) );
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
	generator()->log() << "         ===>" 
			   << "  id1=" << id1 << "  id2=" << id2 
			   << "  flav1=" << flav1 << "  flav2=" << flav2 
			   << endl << endl;
      }
    }
  } else {
    throw Exception() << "HadronSelector::massLightestHadron() could not set a mass."
		      << Exception::abortnow;
  }
  return mass;
}

Energy HadronSelector::massLightestBaryonPair(const long id1, 
					      const long id2) const {
  // Make sure that we don't have any diquarks as input, return arbitrarily
  // large value if we do
  static const Energy numericMax = 1e10*GeV;
  if(abs(id1) > ParticleID::t || abs(id2) > ParticleID::t) return numericMax;
  int f1 = convertIdToFlavour(id1);
  int f2 = convertIdToFlavour(id2);
  int currentLow = DD;
  Energy currentSum = numericMax; 
  for(int i = DD; i<=BB; i++) {
    if(_table[f1][i].size() == 0 || _table[i][f2].size() == 0) continue; 
    Energy s = _table[f1][i].begin()->mass + _table[i][f2].begin()->mass;
    if(currentSum > s) {
      currentSum = s;
      currentLow = i;
    }
  }
  return currentSum;
}

pair<long,long> HadronSelector::
lightestHadronPair(const long id1, const long id2, const long id3) const {

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
	  generator()->logWarning( Exception("HadronSelector::lightestHadronsPair "
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


Energy HadronSelector::
massLightestHadronPair(const long id1, const long id2, const long id3) const {
  pair<long,long> pairId = lightestHadronPair(id1,id2,id3);
  if ( ! getParticleData( pairId.first )  ||  ! getParticleData( pairId.second ) ) {
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Hadronization ) {
      generator()->logWarning( Exception("HadronSelector::masslightestHadronsPair "
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


pair<long,long> HadronSelector::
chooseHadronPair(const Energy cluMass, const long id1, const long id2, 
		 const long id3) throw(Veto, Stop, Exception) {

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
     
  if(HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Hadronization) {
    safetyCheck(id1,id2,id3);
  }

  int maxFlav = BB;
  if(id1 >= DD || id2 >= DD) maxFlav = B;
  pair<long,long> lighthad = lightestHadronPair(id1,id2);
  if(lighthad.first == 0 || lighthad.second == 0) 
     cout << "We have 0's! First id = " << id1 << " second = " << id2 << endl;
  Energy PCMax = Kinematics::pstarTwoBodyDecay(cluMass, 
				  getParticleData(lighthad.first)->mass(),
				  getParticleData(lighthad.second)->mass());
  if(_ClusterDKMode == 0) return hw64(cluMass, id1, id2, PCMax, maxFlav);
  else if(_ClusterDKMode==1) return kupco(cluMass, id1, id2, PCMax, maxFlav);
  else return hwpp(cluMass, id1, id2, PCMax, maxFlav); 
}

pair<long,long> HadronSelector::hw64(const Energy cluMass, const long id1, 
				     const long id2, Energy PCMax, int maxFlav)
{
  pair<long,long> hadPair = pair<long,long>(0,0);
  long had1,had2;
  had1 = had2 = 0;
  int flav1 = convertIdToFlavour(id1);
  int flav2 = convertIdToFlavour(id2);
  //maxFlav = S;
  int mFlav;
  int ntry = 0;
  int idQ = -123456789;
  const int nmax = 5000;
  Energy p;
  /*for(mFlav = S; mFlav<=maxFlav; mFlav++) {
     if(cluMass < _table[flav1][mFlav].begin()->mass + 
	_table[mFlav][flav2].begin()->mass) break;
  } */
  mFlav = maxFlav;  
  do {
    int flav = int(rnd()*double(mFlav));
    KupcoData::iterator it1,it2;
    /*if(cluMass < _table[flav1][flav].begin()->mass + 
                 _table[flav][flav2].begin()->mass) {
      // do nothing
      } else*/
    if(_Pwt[flav] > rnd()) {
      int num1,num2;
      idQ = convertFlavourToId(flav);
      do { 
	num1 = int(double(_table[flav1][flav].size())*rnd()); 
	for(it1 = _table[flav1][flav].begin(); num1; ) { num1--; it1++; }
      } while(it1 != _table[flav1][flav].end() && it1->overallWeight < rnd());
      if(it1!=_table[flav1][flav].end()) had1 = it1->id;
      else had1 = 0;
      do {
	num2 = int(double(_table[flav][flav2].size())*rnd());
	for(it2 = _table[flav][flav2].begin(); num2; ) { num2--; it2++; }
      } while(it2 != _table[flav][flav2].end() && it2->overallWeight < rnd());
      if(it2!=_table[flav][flav2].end()) had2 = it2->id;
      else had2 = 0;  
    } else {
      throw Exception() << "HadronSelector::hw64() did not initialize idQ"
			<< Exception::abortnow;
    }
    
    if(had1 && had2) {
      p = Kinematics::pstarTwoBodyDecay(cluMass, it1->mass, it2->mass);
      if(p/PCMax < rnd()) { had1 = 0; had2 = 0; ntry++;}
    } 
  } while((had1 == 0 || had2 == 0) && ntry < nmax);
  if(ntry >= nmax) return lightestHadronPair(id1,id2);
  else hadPair = pair<long,long>(had1,had2);
  int signQ = 0, signHad1 = -12345678, signHad2 = -12345678;
  if((CheckId::canBeMeson(id1,-idQ) || CheckId::canBeBaryon(id1,-idQ)) &&
     (CheckId::canBeMeson(idQ,id2)  || CheckId::canBeBaryon(idQ,id2))) 
	  signQ = +1;
  if((CheckId::canBeMeson(id1,idQ ) || CheckId::canBeBaryon(id1,idQ)) &&
     (CheckId::canBeMeson(-idQ,id2) || CheckId::canBeBaryon(-idQ,id2)))
	  signQ = -1;
  if(signQ) {
    signHad1 = signHadron(id1, -signQ*idQ, had1);
    signHad2 = signHadron(id2, signQ*idQ, had2);
  } else {
    throw Exception() << "HadronSelector::hw64() did not initialize signHad1/2"
		      << Exception::abortnow;
  }
  hadPair = pair<long,long>(signHad1*had1, signHad2*had2);
  return hadPair;
}

pair<long,long> HadronSelector::kupco(const Energy cluMass, const long id1,
				      const long id2, Energy PCMax, 
				      int maxFlav) { 
  multiset<Kupco> weights;
  pair<long,long> hadPair = pair<long,long>(0,0);

  // If both id1 and id2 are not diquarks, then the possible Q values 
  // (of the Q Qbar pair from the vacuum) are assumed to be: 
  //    u, d, s, c, b, ud, dd, uu, sd, su, ss, 
  //                   cd, cu, cs, cc, bd, bu, bs, bc, bb
  // otherwise only u, d, s, c, b. 
  //int numFlavs = BB;
  //if(CheckId::isDiquark(id1) || CheckId::isDiquark(id2)) numFlavs = B;
  
  // Fill Kupco's cluster decay channel table with all channels 
  // above threshold.
  Energy sumWeight = 0.0;  // store the sum weight on the table.
  //Energy maxWeight = Energy();
  Energy p, weight;
  int flav1 = convertIdToFlavour(id1);
  int flav2 = convertIdToFlavour(id2);
  weights.clear();
  if(flav1 != -1  &&  flav2 !=-1 ) {
    for(int i=D; i <= maxFlav; ++i) {
      if(cluMass > (_table[flav1][i].begin()->mass + 
		    _table[i][flav2].begin()->mass)) { 
	// Loop over all hadron pairs with given flavour.
	for(KupcoData::iterator H1 = _table[flav1][i].begin(); 
	    H1 != _table[flav1][i].end(); H1++) {
	  for(KupcoData::iterator H2 = _table[i][flav2].begin();
	      H2 != _table[i][flav2].end(); H2++) {
	    if(cluMass < H1->mass + H2->mass) break;
	    p = Kinematics::pstarTwoBodyDecay(cluMass, H1->mass, H2->mass );
	    if(p > Energy()) {  // the hadrons are above threshold.
	      weight = _Pwt[i]*p*H1->overallWeight*H2->overallWeight;
	      //if(weight > maxWeight || maxWeight==0.) maxWeight = weight;
	      //sumWeight += weight;
	      int signQ = 0;
	      long idQ = convertFlavourToId(i);
#define m(a,c) CheckId::canBeMeson(a,c)
#define b(a,c) CheckId::canBeBaryon(a,c)
	      if((m(id1,-idQ) || b(id1,-idQ)) && 
		 (m(idQ,id2) || b(idQ,id2)))
		signQ = +1;
	      else if((m(id1,idQ) || b(id1,idQ)) && 
		      (m(-idQ,id2) || b(-idQ,id2)))
		signQ = -1;
#undef m
#undef b
	      int signHad1 = 0, signHad2 = 0;
	      if(signQ != 0) {
		signHad1 = signHadron(id1, -signQ*idQ, H1->id);
		signHad2 = signHadron(id2,  signQ*idQ, H2->id);
	      }
	      if(signHad1 != 0  &&  signHad2 != 0) {
		Kupco a;
		a.idQ = signQ * idQ;
		a.idHad1 = signHad1 * H1->id;
		a.idHad2 = signHad2 * H2->id;
		a.weight = weight;
		sumWeight += weight;
		weights.insert(a);
	      } else {
		ostringstream s;
		s << "HadronSelector::chooseHadronsPair " 
		  << "***Inconsistent Hadron*** signQ " << signQ 
		  << ", " << signHad1 << ", " << signHad2;
		generator()->logWarning(Exception(s.str(),
						  Exception::warning));
	      }
	    } // if pCmstar > 
	  } // for Had2
	} // for Had1
      } // if clumass
    } // for i
  } // if flav1 != 0
  
  // Choose one decay channel.
  if ( weights.size() > 0 ) {
    //unsigned int iChan;      
    // print the full Kupco table for 
    /* if (cluMass > 499*GeV && generator()->currentEventNumber() == 1 
       && HERWIG_DEBUG_LEVEL == 66) {
       Energy M[6], P[6], pcm; 
       double sumw[6];
       int k, fl1, fl2;
       for (k=0; k<6; k++) sumw[k] = 0.0;
       M[0] = .5*GeV; M[1] = 1*GeV; M[2] = 2*GeV; M[3] = 4*GeV; 
       M[4] = 16*GeV; M[5] = 32*GeV; 
       for(iChan = 0; iChan<weights.size(); iChan++) {
       fl1 = weights[iChan].idHad1;
       fl2 = weights[iChan].idHad2;
       pcm = Kinematics::
       pstarTwoBodyDecay(cluMass, getParticleData(fl1)->mass(), 
       getParticleData(fl2)->mass());
       for (k=0; k<6; k++) {
       P[k] = Kinematics::
       pstarTwoBodyDecay(M[k], getParticleData(fl1)->mass(), 
       getParticleData(fl2)->mass());
       if (pcm > 0) sumw[k] += weights[iChan].weight/pcm*P[k];	
       }	  
       }
       for(iChan = 0; iChan<weights.size(); iChan++) {
       fl1 = weights[iChan].idHad1;
       fl2 = weights[iChan].idHad2;
       cout << setw(3) << iChan 
       << setw(6) << weights[iChan].idQ
       << setw(9) << weights[iChan].idHad1
       << setw(9) << weights[iChan].idHad2
       << setw(15) << getParticleData(weights[iChan].idHad1)->PDGName()
       << setw(15) << getParticleData(weights[iChan].idHad2)->PDGName();
       for (k=0; k<6; k++) {
       pcm = Kinematics::
       pstarTwoBodyDecay(cluMass, getParticleData(fl1)->mass(), 
       getParticleData(fl2)->mass());
       P[k] = Kinematics::
       pstarTwoBodyDecay(M[k], getParticleData(fl1)->mass(), 
       getParticleData(fl2)->mass());
       if (pcm > 0 && sumw[k] > 0) {
       cout << setw(10) << setprecision(4) 
       << weights[iChan].weight/pcm*P[k]/sumw[k];	
       } else {
       if (pcm > 0) cout << setw(8) << setprecision(5) << sumw[k];
       else cout << setw(8) << setprecision(5) << "-1";
       }
       }
       cout << endl;
       }
       }*/
    multiset<Kupco>::iterator it;
    double sumWt = 0.0;
    double r;
    //sumWeight = 0.0;
    //for(it = weights.begin(); it != weights.end(); it++) sumWeight+=it->weight;
    r = rnd(sumWeight);
    for(it = weights.begin(); it != weights.end(); it++) {
      sumWt += it->weight;
      if(r<=sumWt) break;
    }
    if(it == weights.end()) {
      ostringstream os;
      os << "HadronSelector::chooseHadronsPair ***No channel found! "
	 << "Setting to first*** size of list is " << weights.size() 
	 << " Sum of weights is " << sumWeight << " r is " << r << " sumwt is "
	 << sumWt;
      generator()->logWarning(Exception(os.str(), Exception::warning));
      it = weights.begin();
    }
    hadPair = pair<long,long>(it->idHad1, it->idHad2);
  }
  return hadPair;
}
pair<long,long> HadronSelector::hwpp(const Energy cluMass, const long id1,
			  	     const long id2, Energy PCMax, 
				     int maxFlav) { 
  multiset<Kupco> weights;
  pair<long,long> hadPair = pair<long,long>(0,0);
  
  // Fill Kupco's cluster decay channel table with all channels 
  // above threshold.
  Energy sumWeight = 0.0;  // store the sum weight on the table.
  Energy p, weight;
  int flav1 = convertIdToFlavour(id1);
  int flav2 = convertIdToFlavour(id2);
  weights.clear();
  int startFlav = D;
  // Choose the meson sector if between 0 and 1/(1+B), otherwise choose 
  // baryon sector
  //Energy a = massLightestBaryonPair(id1,id2);
  if(cluMass > massLightestBaryonPair(id1,id2)) {
    double r = rnd();
    if(r < 1./(1.+pwtDIquark())) { startFlav = D; maxFlav = B; }
    else { startFlav = DD; maxFlav = BB; }
  } else { startFlav = D; maxFlav = B; }

  if(flav1 == -1  ||  flav2 ==-1 ) return pair<long,long>(0,0);

  //cout << "startflav = " << startFlav << " and maxFlav = " << maxFlav << " cluMass = " << cluMass << " a is " << a << endl;
  for(int i=startFlav; i <= maxFlav; ++i) {
    if(cluMass > (_table[flav1][i].begin()->mass + 
		  _table[i][flav2].begin()->mass)) { 
      // Loop over all hadron pairs with given flavour.
      for(KupcoData::iterator H1 = _table[flav1][i].begin(); 
	  H1 != _table[flav1][i].end(); H1++) {
	for(KupcoData::iterator H2 = _table[i][flav2].begin();
	    H2 != _table[i][flav2].end(); H2++) {
	  if(cluMass < H1->mass + H2->mass) break;
	  p = Kinematics::pstarTwoBodyDecay(cluMass, H1->mass, H2->mass );
	  if(p > Energy()) {  // the hadrons are above threshold.
	    weight = _Pwt[i]*p*H1->overallWeight*H2->overallWeight;
	    int signQ = 0;
	    long idQ = convertFlavourToId(i);
#define m(a,c) CheckId::canBeMeson(a,c)
#define b(a,c) CheckId::canBeBaryon(a,c)
	    if((m(id1,-idQ) || b(id1,-idQ)) && 
	       (m(idQ,id2) || b(idQ,id2)))
	      signQ = +1;
	    else if((m(id1,idQ) || b(id1,idQ)) && 
		    (m(-idQ,id2) || b(-idQ,id2)))
	      signQ = -1;
#undef m
#undef b
	    int signHad1 = 0, signHad2 = 0;
	    if(signQ != 0) {
	      signHad1 = signHadron(id1, -signQ*idQ, H1->id);
	      signHad2 = signHadron(id2,  signQ*idQ, H2->id);
	    }
	    if(signHad1 != 0  &&  signHad2 != 0) {
	      Kupco a;
	      a.idQ = signQ * idQ;
	      a.idHad1 = signHad1 * H1->id;
	      a.idHad2 = signHad2 * H2->id;
	      a.weight = weight;
	      sumWeight += weight;
	      weights.insert(a);
	      //cout << "adding weight for " << a.idHad1 << "x" << a.idHad2 
	      //   << " = " << weight << endl;
	    } else {
	      ostringstream s;
	      s << "HadronSelector::chooseHadronsPair " 
		<< "***Inconsistent Hadron*** signQ " << signQ 
		<< ", " << signHad1 << ", " << signHad2;
	      generator()->logWarning(Exception(s.str(),
						Exception::warning));
	    }
	  } // if pCmstar > 
	} // for Had2
      } // for Had1
    } // if clumass
  } // for i

  // Choose one decay channel.
  if ( weights.size() > 0 ) {
    multiset<Kupco>::iterator it;
    double sumWt = 0.0;
    double r;
    r = rnd(sumWeight);
    for(it = weights.begin(); it != weights.end(); it++) {
      //cout << "Looking at " << it->idHad1 << "x" << it->idHad2 << " = " 
      //   << it->weight << endl;
      sumWt += it->weight;
      if(r<=sumWt) break;
    }
    if(it == weights.end()) {
      ostringstream os;
      os << "HadronSelector::chooseHadronsPair ***No channel found! "
	 << "Setting to first*** size of list is " << weights.size() 
	 << " Sum of weights is " << sumWeight << " r is " << r << " sumwt is "
	 << sumWt;
      generator()->logWarning(Exception(os.str(), Exception::warning));
      it = weights.begin();
    }
    hadPair = pair<long,long>(it->idHad1, it->idHad2);
  } 
  return hadPair;
}

// ------------- PRIVATE ACCESSORY METHODS ----------------------------

void HadronSelector::safetyCheck(const long id1, const long id2, const long id3) const
  throw(Veto, Stop, Exception) {
  // Check if the components can form a meson or a baryon.
  if ( ! CheckId::canBeMeson(id1,id2)  &&  ! CheckId::canBeBaryon(id1,id2,id3) ) {  
    generator()->logWarning( Exception("HadronSelector::safetyCheck "
				       "***The cluster is inconsistent*** ",
				       Exception::warning) );
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
      generator()->log() << "         ===>" 
			 << " id1=" << id1 << " id2=" << id2 << " id3=" << id3 
			 << endl << endl;
    }
  }
}


int HadronSelector::convertIdToFlavour(const long id) const {
  int flavour = -1;
  // Flavours are (starting at 0) D U S C B UU UD DD US DS SS
  if ( abs(id) >= 1  &&  abs(id) <= 5 ) { 
    flavour = abs(id)-1;
  } else {
    switch ( abs(id) ) {
    case 1103: flavour = DD;  break;
    case 2101:
    case 2103: flavour = DU;  break;  
    case 2203: flavour = UU;  break;
    case 3101:
    case 3103: flavour = DS;  break;
    case 3201:
    case 3203: flavour = US;  break;
    case 3303: flavour = SS;  break;
    case 4101:
    case 4103: flavour = DC;  break;
    case 4201:
    case 4203: flavour = UC;  break;
    case 4301:
    case 4303: flavour = SC;  break;
    case 4403: flavour = CC;  break;
    case 5101:
    case 5103: flavour = DB;  break;
    case 5201:
    case 5203: flavour = UB;  break;
    case 5301:
    case 5303: flavour = SB;  break;
    case 5401:
    case 5403: flavour = CB;  break;
    case 5503: flavour = BB;  break;
    }
  }
  return flavour;
}


long HadronSelector::convertFlavourToId(const int flavour) const {
  // Notice that we are always converting diquarks to the 
  // spin 0 states, when available, never to spin 1.
  // This is just a simplifications that hopefully should
  // be harmless (if not, to change is quite straighforward,
  // for example randomly assigned 50% to spin 0 and 50% to spin 1).
  long id = 0;
  switch ( flavour ) {
  case D:   id  = ParticleID::d;     break;
  case U:   id  = ParticleID::u;     break;
  case S:   id  = ParticleID::s;     break;
  case C:   id  = ParticleID::c;     break;
  case B:   id  = ParticleID::b;     break;
  case DD:  id  = ParticleID::dd_1;  break;
  case DU:  id  = ParticleID::ud_0;  break;
  case UU:  id  = ParticleID::uu_1;  break;
  case DS:  id  = ParticleID::sd_0;  break;
  case US:  id  = ParticleID::su_0;  break;
  case SS:  id  = ParticleID::ss_1;  break;
  case DC:  id  = ParticleID::cd_1;  break;
  case UC:  id  = ParticleID::cu_0;  break;
  case SC:  id  = ParticleID::cs_0;  break;
  case CC:  id  = ParticleID::cc_1;  break;
  case DB:  id  = ParticleID::bd_1;  break;
  case UB:  id  = ParticleID::bu_0;  break;
  case SB:  id  = ParticleID::bs_0;  break;
  case CB:  id  = ParticleID::bc_0;  break;
  case BB:  id  = ParticleID::bb_1;  break;
  }
  return id;
}


int HadronSelector::signHadron(const int idQ1, const int idQ2, 
				const int idHad) const {

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
  if(HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Hadronization) {
    if(!(getParticleData(idQ1) && getParticleData(idQ2) && 
	 getParticleData(idHad)  &&
	 (CheckId::canBeMeson(idQ1,idQ2) && CheckId::isMeson(idHad)) ||
	 (CheckId::canBeBaryon(idQ1,idQ2) && CheckId::isBaryon(idHad)))) {
      generator()->logWarning(Exception("HadronSelector::signHadron "
				    "***The input is inconsistent*** ",
				    Exception::warning));
      if(HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization) {    
	generator()->log() << "         ===>" << " idQ1=" << idQ1 << "  idQ2=" 
			   << idQ2 << "  idHad=" << idHad << endl << endl;
      }
    }
  }
  int sign = 0;
  if(getParticleData(idQ1) && getParticleData(idQ2) && getParticleData(idHad)){
    Charge chargeIn  = getParticleData(idQ1)->charge() + 
                       getParticleData(idQ2)->charge();
    Charge chargeOut = getParticleData(idHad)->charge();
    // Delta charge is the charge of d/10
    Charge delta_ch  = fabs(getParticleData(ParticleID::d)->charge())/10.0; 
    if(fabs(chargeIn-chargeOut) < delta_ch && 
       fabs(chargeIn+chargeOut) > delta_ch) {
      sign = +1;
    } else if(fabs(chargeIn-chargeOut) > delta_ch &&  
              fabs(chargeIn+chargeOut) < delta_ch ) {
      sign = -1;
    } else if (fabs(chargeIn) < delta_ch && fabs(chargeOut) < delta_ch ) {  
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
      if(CheckId::isQuark(idQ1) && CheckId::isQuark(idQ2)) {  // meson
	int idQa = abs(idQ1), idQb = abs(idQ2), dominant = idQ2;
	if(idQa > idQb) {
	  idQa = abs(idQ2); idQb = abs(idQ1); dominant = idQ1;
	}
	if((idQa==ParticleID::d && idQb==ParticleID::s) ||
	   (idQa==ParticleID::d && idQb==ParticleID::b) ||
	   (idQa==ParticleID::s && idQb==ParticleID::b)) {
	  if (dominant < 0 || idHad%10 == 0) { // idHad%10 is zero for K0L,K0S
	    sign = +1;
	  } else if(dominant > 0) {
	    sign = -1;
	  }
	} else if(idQa==ParticleID::u && idQb==ParticleID::c) {
	  if(dominant > 0) {
	    sign = +1;
	  } else if(dominant < 0) { 
	    sign = -1;
	  } 
	} else if(idQa==idQb) {
	  sign = +1;
	}
      } else if(CheckId::isDiquark(idQ1) || CheckId::isDiquark(idQ2)) { 
	// baryon
	if(idQ1 > 0 && idQ2 > 0) {
	  sign = +1;
	} else if(idQ1 < 0 && idQ2 < 0) {
	  sign = -1;
	}
      }
    } 
  }
  return sign;
}


// ------------- INITIALIZATION  ----------------------------

void HadronSelector::initialize() {

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

  //_Pwt[0]  = 0.0; // not used.
  _Pwt[D]  = _PwtDquark;
  _Pwt[U]  = _PwtUquark;
  _Pwt[S]  = _PwtSquark;
  _Pwt[C]  = _PwtCquark;
  _Pwt[B]  = _PwtBquark;
  _Pwt[DD] =       _PwtDIquark * _PwtDquark * _PwtDquark;
  _Pwt[DU] = 0.5 * _PwtDIquark * _PwtUquark * _PwtDquark;
  //_Pwt[DU] = _PwtDIquark * _PwtUquark * _PwtDquark;
  _Pwt[UU] =       _PwtDIquark * _PwtUquark * _PwtUquark;
  _Pwt[DS] = 0.5 * _PwtDIquark * _PwtSquark * _PwtDquark;
  //_Pwt[DS] = _PwtDIquark * _PwtSquark * _PwtDquark;
  _Pwt[US] = 0.5 * _PwtDIquark * _PwtSquark * _PwtUquark;
  //_Pwt[US] = _PwtDIquark * _PwtSquark * _PwtUquark;
  _Pwt[SS] =       _PwtDIquark * _PwtSquark * _PwtSquark;
  // Commenting out heavy di-quark weights
  _Pwt[DC] = 0.0;
  _Pwt[UC] = 0.0;
  _Pwt[SC] = 0.0;
  _Pwt[CC] = 0.0;
  _Pwt[DB] = 0.0;
  _Pwt[UB] = 0.0;
  _Pwt[SB] = 0.0;
  _Pwt[CB] = 0.0;
  _Pwt[BB] = 0.0;
  //_Pwt[DC] = 0.5 * _PwtDIquark * _PwtCquark * _PwtDquark;
  //_Pwt[UC] = 0.5 * _PwtDIquark * _PwtCquark * _PwtUquark;
  //_Pwt[SC] = 0.5 * _PwtDIquark * _PwtCquark * _PwtSquark;
  //_Pwt[CC] =       _PwtDIquark * _PwtCquark * _PwtCquark;
  //_Pwt[DB] = 0.5 * _PwtDIquark * _PwtBquark * _PwtDquark;
  //_Pwt[UB] = 0.5 * _PwtDIquark * _PwtBquark * _PwtUquark;
  //_Pwt[SB] = 0.5 * _PwtDIquark * _PwtBquark * _PwtSquark;
  //_Pwt[CB] = 0.5 * _PwtDIquark * _PwtBquark * _PwtCquark;
  //_Pwt[BB] =       _PwtDIquark * _PwtBquark * _PwtBquark; 
  double pmax = 0.0;
  for(int pw = D; pw<=BB; pw++) if(pmax<_Pwt[pw]) pmax = _Pwt[pw];
  for(int pw = D; pw<=BB; pw++) _Pwt[pw]/=pmax;
  for (int l = 0; l < Lmax; ++l ) {
    for (int j = 0; j < Jmax; ++j) {
      for (int n = 0; n < Nmax; ++n) {
	_Repwt[l][j][n] = 1.0;
      }
    }
  }
  
  fillHadronData();

  // General debugging infos
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
    for (int i=D; i <= B; ++i) {
      for (int j=i; j <= NumFlavs; ++j) {
	generator()->log() << "Hadrons for " << i << " and " << j << endl;
	int num = 1;
	for (KupcoData::iterator it = _table[i][j].begin(); it!=_table[i][j].end(); it++) {
	  generator()->log() << "\t   " << num++ << "   " 
			     << (it->ptrData ? it->ptrData->PDGName() : "UNKNOWN")
			     << "  id=" << it->id 
			     << "  wt=" << it->wt << "  swtef=" << it->swtef
			     << "  overallWeight=" << it->overallWeight 
			     << "  mass=" << it->mass << endl;
	}
      }
    }
    generator()->log() << "  --- END INITIALIZE --- " << endl << endl;
  }      
}


void HadronSelector::fillHadronData() {
  // Here the hadrons data must be filled.
  map<int,int> flavToIndex;
  map<int,string> itoS;
  // Set this simple map up to save some hassle later
  itoS[flavToIndex[1] = D] = "d";
  itoS[flavToIndex[2] = U] = "u";
  itoS[flavToIndex[3] = S] = "s";   
  itoS[flavToIndex[4] = C] = "c";
  itoS[flavToIndex[5] = B] = "b";
  itoS[flavToIndex[22] = UU] = "uu";
  itoS[flavToIndex[12] = DU] = "du";
  itoS[flavToIndex[11] = DD] = "dd";
  itoS[flavToIndex[23] = US] = "us";
  itoS[flavToIndex[13] = DS] = "ds";
  itoS[flavToIndex[33] = SS] = "ss";
  itoS[flavToIndex[14] = DC] = "dc";
  itoS[flavToIndex[15] = DB] = "db";
  itoS[flavToIndex[24] = UC] = "uc";
  itoS[flavToIndex[25] = UB] = "ub";
  itoS[flavToIndex[34] = SC] = "sc";
  itoS[flavToIndex[35] = SB] = "sb";
  itoS[flavToIndex[44] = CC] = "cc";
  itoS[flavToIndex[45] = CB] = "cb";
  itoS[flavToIndex[55] = BB] = "bb";
  ParticleMap particles = generator()->particles();
  double maxdd,maxss,maxrest;
  maxdd = maxss = maxrest = 0.0;
  for(ParticleMap::iterator it=particles.begin(); it!=particles.end(); it++) {
    // Don't include non-hadrons
    if(it->second->id() < 100) continue;
    // Only include those with 2J+1 less than...5
    if(_trial==2 && it->second->id()%10 >= 5) continue;
    // Only include those with 2J+1 less than...7
    if(_trial==3 && it->second->id()%10 >= 7) continue;
    // Only include pions
    if(_trial==1 && it->second->id()!=111 && it->second->id()!=211) continue;

    int flav1, flav2;
    // Get the flavours
    int x4 = (it->second->id()/1000)%10; 
    int x3 = (it->second->id()/100 )%10;
    int x2 = (it->second->id()/10  )%10;
    if(x3 == 0 || x2 == 0) {
      continue; // Skip non-hadrons (susy particles, etc...)
    } else if(x4 == 0) { 
      flav1 = x3; 
      flav2 = x2; 
    } else {
      if(x3 > x2) flav1 = 10*x2 + x3;
      else flav1 = 10*x3 + x2;
      flav2 = x4;
    }
    // Get spin
    int nj = it->second->id()%10;
    // if j is 0 and it isn't the K0S, then it is an abnormal state, ignore it
    if(nj==0) continue;// && it->second->id() != 310) continue;
    //if(nj==0) nj = 1;

    HadronInfo a;
    a.id = it->second->id();
    a.ptrData = it->second;
    a.swtef = specialWeight(it->second->id());
    a.mass = it->second->mass();

    if(flav1 == flav2 && (flav1 == 3 || flav1 == 2)) {
      // load up ssbar> uubar> ddbar> admixture states
      if(mixingState(it->second->id()))
	a.wt = mixingStateWeight(it->second->id());
      else a.wt = 0.33;
      if(_ClusterDKMode != 0) 
	a.overallWeight = 1./sqrt(2.) * a.wt * a.swtef * nj;
      else 
	a.overallWeight = a.wt * nj;
      _table[D][D].insert(a);
      _table[U][U].insert(a);
      if(_ClusterDKMode == 0 && a.overallWeight > maxdd) 
	maxdd = a.overallWeight;
      if(_ClusterDKMode != 0) a.wt = 1.-sqrt(2.)*a.wt;
      else a.wt = 1-a.wt;
      if(a.wt > 0) {
	a.overallWeight = a.wt * a.swtef * nj;
	_table[S][S].insert(a);
	if(_ClusterDKMode == 0 && a.overallWeight > maxss) 
	  maxss = a.overallWeight;
      }
    } else if(flav1 == 1 && flav2 == 1) {
      // load up ddbar> uubar> admixtures
      if(_ClusterDKMode != 0) a.overallWeight = 1./sqrt(2.) * a.wt * nj;
      else a.overallWeight = a.wt * nj;
      _table[D][D].insert(a);
      _table[U][U].insert(a);
      if(_ClusterDKMode == 0 && a.overallWeight > maxdd) 
	maxdd = a.overallWeight;
    } else if((flav1 == 1 && flav2 == 11) || (flav1 == 11 && flav2 == 1) ||
	      (flav1 == 2 && flav2 == 22) || (flav1 == 22 && flav2 == 2) ||
	      (flav1 == 3 && flav2 == 33) || (flav1 == 33 && flav2 == 3)) {
      if(_ClusterDKMode != 0) a.overallWeight = 1.5*a.wt*a.swtef*nj;
      else a.overallWeight = a.wt * nj;
      _table[flavToIndex[flav1]][flavToIndex[flav2]].insert(a);
      _table[flavToIndex[flav2]][flavToIndex[flav1]].insert(a);
      if(_ClusterDKMode == 0 && a.overallWeight > maxrest)
	maxrest = a.overallWeight;
    } else {
      if(_ClusterDKMode != 0) a.overallWeight = a.wt * a.swtef * nj;
      else a.overallWeight = a.wt * nj;
      _table[flavToIndex[flav1]][flavToIndex[flav2]].insert(a);
      if(flav1 != flav2)
	_table[flavToIndex[flav2]][flavToIndex[flav1]].insert(a);
      if(_ClusterDKMode == 0 && a.overallWeight > maxrest) 
	maxrest = a.overallWeight;
    }
  }
  // Account for identical combos of diquark/quarks and symmetrical elements
  // e.g. U UD = D UU
  for(int i = D; i<=C; i++) {
    for(int j = DD; j<=BB; j++) {
      int k, l, sub;
      if(j >= DB) {
	k = B;
	sub = DB;
      } else if(j >= DC) {
	k = C;
	sub = DC;
      } else if(j >= DS) {
	k = S;
	sub = DS;
      } else if(j >= DU) {
	k = U;
	sub = DU;
      } else if(j == D) {
	k = U;
	sub = DU;
      } else continue;
      if(j-sub > i) l = flavToIndex[10*(i+1)+(j-sub+1)];
      else l = flavToIndex[10*(j-sub+1) + (i+1)];
      // Only do ones that haven't been filled
      if(k < DD && l < DD) continue;
      if(_table[i][j].size() == 0) _table[i][j] = _table[j][i] = _table[k][l];
    }
  }
  if(_ClusterDKMode == 0) {
    for(int i = D; i<BB; i++) {
      for(int j = D; j<BB; j++) {
	for(KupcoData::iterator it = _table[i][j].begin(); 
	    it!=_table[i][j].end(); it++) {
	  if(i==j && (i==D || i==U)) it->rescale(1./maxdd);
	  else if(i==j && i==S) it->rescale(1./maxss);
	  else it->rescale(1./maxrest);
	}
      }
    }
  }
  if(HERWIG_DEBUG_LEVEL >= 66) {
    ofstream kupco1;
    kupco1.open("kupco_table/table_size.dat");
    kupco1 << setw(4) << "";
    for(int i = 0; i<20; i++)
      kupco1 << setw(4) << itoS[i];
    kupco1 << endl;
    //  kupco1 << "\td\tu\ts\tc\tb\tuu\tud\tdd\tus\tds\tss\n";
    for(int i = 0; i<20; i++) {
      kupco1 << setw(4) << itoS[i] << setw(4);
      for(int j = 0; j<20; j++) {
	kupco1 << _table[i][j].size();
	if(j != 19) kupco1 << setw(4);
	else kupco1 << endl;
      }
    }
    kupco1.close();
    for(int i = D; i<=B; i++) {
      for(int j = i; j<=BB; j++) {
	ofstream flavOut;
	string s = "kupco_table/" + itoS[i] + itoS[j] + ".dat";
	flavOut.open(s.c_str());
	for(KupcoData::iterator it = _table[i][j].begin(); 
	    it != _table[i][j].end(); it++)
	  flavOut << setw(10) << it->ptrData->PDGName() << "\t" 
		<< setw(10) << it->overallWeight << "\t"
		<< setw(10) << it->mass << "\t" << setw(8) << it->id << endl;
	flavOut.close();
      }
    }
  }
}

double HadronSelector::specialWeight(long id) {
  int nj = id % 10;
  if(nj == 0) nj = 1;  // K0L and K0S only have nj == 0 
  if(nj == 2 || nj == 4) {     // Baryon : J = 1/2 or 3/2
    if(nj == 2) {
      int nq3 = (id/10 )%10;
      int nq2 = (id/100)%10;
      if(nq2 < nq3) return  sqr(_SngWt);   // Singlet (Lambda-like) baryon
    } else return sqr(_DecWt);    // Decuplet baryon
  } else if(nj % 2 == 1) {   // Meson
    int j  = (nj - 1) / 2;                     // Total angular momentum
    int nl = (id/10000 )%10;  // related to Orbital angular momentum l
    int l  = -999;  
    int n  = (id/100000)%10;  // Radial excitation
    if(j == 0) l = nl;
    else if(nl == 0) l = j - 1;
    else if((nl == 1) || (nl == 2)) l = j;
    else if(nl == 3) l = j + 1;
    if((l||j||n) && l>=0  &&  l<Lmax  &&  j<Jmax  &&  n<Nmax) {
      return sqr(_Repwt[l][j][n]);  // Angular or Radial excited meson
    }    
  }
  return 1.0;
}

bool HadronSelector::mixingState(long id) {
  switch(id) {
  case ParticleID::eta:
  case ParticleID::etaprime: 
  case ParticleID::omega:
  case ParticleID::phi:
  case ParticleID::h_1:
  case ParticleID::hprime_1:
  case ParticleID::f_0:
  case ParticleID::fprime_1:
  case ParticleID::f_1:
  case ParticleID::fprime_2:
  case ParticleID::f_2: 
    return true;
  default:
    return false;
  }
  return false;
}

double HadronSelector::mixingStateWeight(long id) {

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
  //               particles are not (yet) present in ThePEG/PDT/EnumParticles.h.
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
  switch(id) {
  case ParticleID::eta:      return 0.5*CheckId::probabilityMixing(etamix,1);
  case ParticleID::etaprime: return 0.5*CheckId::probabilityMixing(etamix,2);
  case ParticleID::omega:    return 0.5*CheckId::probabilityMixing(phimix,2);
  case ParticleID::phi:      return 0.5*CheckId::probabilityMixing(phimix,1);
  case ParticleID::h_1:      return 0.5*CheckId::probabilityMixing(h1mix,2);
  case ParticleID::hprime_1: return 0.5*CheckId::probabilityMixing(h1mix,1); 
  case ParticleID::f_0:      return 0.5*CheckId::probabilityMixing(f0mix,2);
  case ParticleID::fprime_1: return 0.5*CheckId::probabilityMixing(f1mix,1);
  case ParticleID::f_1:      return 0.5*CheckId::probabilityMixing(f1mix,2);
  case ParticleID::fprime_2: return 0.5*CheckId::probabilityMixing(f2mix,1);
  case ParticleID::f_2:      return 0.5*CheckId::probabilityMixing(f2mix,2);
  default:                   return 1./3.;
  }
  return 1.0;
}

