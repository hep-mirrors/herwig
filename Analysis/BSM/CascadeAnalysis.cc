// -*- C++ -*-
//
// CascadeAnalysis.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CascadeAnalysis class.
//

#include "CascadeAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "Herwig++/Models/Susy/MSSM.h"

using namespace Herwig;

CascadeAnalysis::CascadeAnalysis() : 
  theModel(0), theMassNormalize(true), theNBins(70), theResonances(7,0), 
  theMassMaxima(6), theX(0.), theY(0.), theZ(0.), theYc(0.), theZc(0.), 
  theNChain(6,0), theFractions(6,0), theAlpha(0.) {}

void CascadeAnalysis::analyze(tEventPtr event, long, int, int) {
  ParticleVector::const_iterator pit = 
    event->primaryCollision()->primarySubProcess()->outgoing().begin();
  ParticleVector::const_iterator pend = 
    event->primaryCollision()->primarySubProcess()->outgoing().end();

  for( ; pit != pend; ++pit ) {
    if( abs((*pit)->id()) == theResonances[0] ) {
      tPPtr parent = showeredProduct(*pit);
      ParticleVector products = parent->children();
      if( products.size() != 2 ) continue;
      if( !products[0]->data().charged() || 
	  !products[1]->data().charged() )
	analyseNeutralChain(parent, products);
      else
	analyseChargedChain(parent, products);
    }
    else if(abs((*pit)->id()) == theResonances[0] - 1 ) {
      tPPtr parent = showeredProduct(*pit);
      ParticleVector products = parent->children();
      if( products.size() != 2 ) continue;
      analyseChargedChain(parent, products);
    }
    else {}
  }

}

void CascadeAnalysis::persistentOutput(PersistentOStream & os) const {
  os << theModel  << theMassNormalize << theNBins;
}

void CascadeAnalysis::persistentInput(PersistentIStream & is, int) {
  is >> theModel  >> theMassNormalize >> theNBins;
}

ClassDescription<CascadeAnalysis> CascadeAnalysis::initCascadeAnalysis;
// Definition of the static class description member.

void CascadeAnalysis::Init() {

  static ClassDocumentation<CascadeAnalysis> documentation
    ("There is no documentation for the CascadeAnalysis class");
  
  static Switch<CascadeAnalysis,unsigned int> interfaceModelName
    ("ModelName",
     "The model under study",
     &CascadeAnalysis::theModel, 0, false, false);
  static SwitchOption interfaceModelNameMSSM
    (interfaceModelName,
     "MSSM",
     "The minimal-supersymmetric-standard-model",
     0);
  static SwitchOption interfaceModelNameMUED
    (interfaceModelName,
     "MUED",
     "The minimal universal extra dimensions",
     1);
  
  static Switch<CascadeAnalysis,bool> interfaceNormalizeMass
    ("NormalizeMass",
     "Normalize the mass to the maximum",
     &CascadeAnalysis::theMassNormalize, true, false, false);
  static SwitchOption interfaceNormalizeMassNormalize
    (interfaceNormalizeMass,
     "Yes",
     "Normalize the masses",
     true);
  static SwitchOption interfaceNormalizeMassNoNormalize
    (interfaceNormalizeMass,
     "No",
     "Do not normalize the masses",
     false);

  static Parameter<CascadeAnalysis,unsigned int> interfaceNumberOfBins
    ("NumberOfBins",
     "The number of bins for each histogram",
     &CascadeAnalysis::theNBins, 70, 0, 1000,
     false, false, Interface::limited);
}

void CascadeAnalysis::
analyseNeutralChain(tPPtr inpart, const ParticleVector & products) {
  if( products.empty() ) return;
//   //D->C,c
  tPPtr partC, quark;
  long tid0(abs(products[0]->id())), tid1(abs(products[1]->id()));
  if( tid0 == theResonances[1] ) { 
    partC = products[0];
    quark = products[1];
  } else if( tid1 == theResonances[1] ) { 
    partC = products[1];
    quark = products[0];
  } else return;
  //C->B,b
  ParticleVector prd = showeredProduct(partC)->children();
  if( prd.size() != 2 ) return;
  tPPtr partB, ln;
  tid0 = abs(prd[0]->id());
  tid1 = abs(prd[1]->id());
  if( ( tid0 == ParticleID::muminus &&
        tid1 == theResonances[2] + 2) ||
      ( tid0 == ParticleID::eminus &&
        tid1 == theResonances[2]) ) { 
    ln = prd[0];
    partB = prd[1];    
  } else if( (tid1 == ParticleID::muminus &&
	      tid0 == theResonances[2] + 2) ||
	     (tid1 == ParticleID::eminus &&
	      tid0 == theResonances[2]) ) {
    ln = prd[1];
    partB = prd[0];    
  } else return;
    //B->A,a
  prd = showeredProduct(partB)->children();
  if( prd.size() != 2 ) return;
  tid0 = abs(prd[0]->id());
  tid1 = abs(prd[1]->id());
  tPPtr lf, partA;
  if( (tid0 == ParticleID::muminus &&
       tid1 == theResonances[3]) ||
      (tid0 == ParticleID::eminus &&
       tid1 == theResonances[3]) ) { 
    lf = prd[0];
    partA = prd[1];    
  } else if( (tid1 == ParticleID::muminus &&
	      tid0 == theResonances[3]) ||
	     (tid1 == ParticleID::eminus &&
	      tid0 == theResonances[3]) ) { 
    lf = prd[1];
    partA = prd[0];    
  } else return;
  
  // Invariant masses normalized appropriately according to user options
  double qln = sqrt( (quark->momentum() + ln->momentum()).m2() )/theMassMaxima[0];
  double qlf = sqrt( (quark->momentum() + lf->momentum()).m2() )/theMassMaxima[1];
  double lnlf = sqrt( (ln->momentum() + lf->momentum()).m2() )/theMassMaxima[2];
  
  if( inpart->id() > 0 ) {
    ++theNChain[0];
    if( ln->id() > 0 ) {
      *theHistograms["qlnm"] += qln;
      *theHistograms["qlm"] += qln*theMassMaxima[0]/theMassMaxima[3];
    }
    else { 
      *theHistograms["qlnp"] += qln;
      *theHistograms["qlp"] += qln*theMassMaxima[0]/theMassMaxima[3];
    }

    if( lf->id() > 0 ) {
      *theHistograms["qlfm"] += qlf;
      *theHistograms["qlm"] += qlf*theMassMaxima[1]/theMassMaxima[3];
    }
    else { 
      *theHistograms["qlfp"] += qlf;
      *theHistograms["qlp"] += qlf*theMassMaxima[1]/theMassMaxima[3];
    }
  }
  else {
    ++theNChain[1];
    if( ln->id() > 0 ) {
      *theHistograms["qblnm"] += qln;
      *theHistograms["qlm"] += qln*theMassMaxima[0]/theMassMaxima[3];
    }
    else { 
      *theHistograms["qblnp"] += qln;
      *theHistograms["qlp"] += qln*theMassMaxima[0]/theMassMaxima[3];
    }

    if( lf->id() > 0 ) {
      *theHistograms["qblfm"] += qlf;
      *theHistograms["qlm"] += qlf*theMassMaxima[1]/theMassMaxima[3];
    }
    else { 
      *theHistograms["qblfp"] += qlf;
      *theHistograms["qlp"] += qlf*theMassMaxima[1]/theMassMaxima[3];
    }
  }

  *theHistograms["ll"] += lnlf;
}

void CascadeAnalysis::
analyseChargedChain(tPPtr, const ParticleVector & products) {
  if( products.empty() ) return;
  tPPtr partC, quark;
  long tid0(abs(products[0]->id())), tid1(abs(products[1]->id()));
  if( tid0 == theResonances[4] ) { 
    partC = products[0];
    quark = products[1];
  } else if( tid1 == theResonances[4] ) { 
    partC = products[1];
    quark = products[0];
  } else return;
  if( abs(quark->id()) > 4 ) return;
  //C->W,A
  ParticleVector prd = showeredProduct(partC)->children();
  if( prd.size() != 2 ) return;
  tPPtr wboson, partA;
  tid0 = abs(prd[0]->id());
  tid1 = abs(prd[1]->id());
  if( tid0 == theResonances[5] && tid1 == theResonances[6] ) {
    wboson = prd[0];
    partA = prd[1];    
  } else if( tid0 == theResonances[6] && tid1 == theResonances[5] ) {
    wboson = prd[1];
    partA = prd[0];    
  } else return;
  prd = showeredProduct(wboson)->children();
  if( prd.size() != 2 ) return;
  tPPtr lepton;
  tid0 = abs(prd[0]->id());
  tid1 = abs(prd[1]->id());
  if( tid0 == 11 || tid0 == 13 ) lepton = prd[0];
  else if( tid1 == 11 || tid1 == 13 ) lepton = prd[1];
  else return;

  Energy mlq = sqrt( (quark->momentum() + lepton->momentum() ).m2());

  if( abs(quark->id()) % 2 == 0 ) {
    *theHistograms["W-p1"] += mlq/theMassMaxima[4];
    if( lepton->id() > 0 ) {
      *theHistograms["W-lm"] += mlq/theMassMaxima[4];
      ++theNChain[2];
    }
    else {
      *theHistograms["W-lp"] += mlq/theMassMaxima[4];
      ++theNChain[5];
    }
  }
  else {
    *theHistograms["W-p2"] += mlq/theMassMaxima[5];
    if( lepton->id() > 0 ) {
      *theHistograms["W-lm"] += mlq/theMassMaxima[5];
      ++theNChain[3];
    }
    else {
      *theHistograms["W-lp"] += mlq/theMassMaxima[5];
      ++theNChain[4];
    }
  }

}

tPPtr CascadeAnalysis::showeredProduct(tPPtr p) const {
  if( isLastInShower(*p) ) return p;
  long pid = p->id();
  ParticleVector prds = p->children();
  ParticleVector::size_type n = p->children().size();
  for( ParticleVector::size_type i = 0; i < n; ++i ) {
    tPPtr child = prds[i];
    if( child->id() == pid ) {
      p = showeredProduct(child);
      break;
    }
  }
  return p;
}

void CascadeAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();

  if( theModel == 0 ) {
    //neutralino
    theResonances[0] = 1000002;
    theResonances[1] = 1000023;
    theResonances[2] = 2000011;
    theResonances[3] = 1000022;
    //chargino
    theResonances[4] = 1000037;
    theResonances[5] = 24;
    theResonances[6] = 1000022;
  }
  else {
    theResonances[0] = 5100002;
    theResonances[1] = 5100023;
    theResonances[2] = 6100011;
    theResonances[3] = 5100022;
    //W_1
    theResonances[4] = 5100024;
    theResonances[5] = 24;
    theResonances[6] = 5100022;
  }

  Energy2 mD2 = sqr(getParticleData(theResonances[0])->mass());
  Energy2 mC2 = sqr(getParticleData(theResonances[1])->mass());
  Energy2 mB2 = sqr(getParticleData(theResonances[2])->mass());
  Energy2 mA2 = sqr(getParticleData(theResonances[3])->mass());
  theX = mC2/mD2;
  theY = mB2/mC2;
  theZ = mA2/mB2;
  // Maximum invariant mass
  double mln_max = sqrt( (mD2 - mC2)*(mC2 - mB2)/mC2 )/GeV;
  double mlf_max = sqrt( (mD2 - mC2)*( mB2 - mA2)/mB2 )/GeV;
  double mll_max = sqrt( (mC2 - mB2)*(mB2 - mA2)/mB2 )/GeV;
  
  mC2 = sqr(getParticleData(theResonances[4])->mass());
  mB2 = sqr(getParticleData(theResonances[5])->mass());
  mA2 = sqr(getParticleData(theResonances[6])->mass());
  theYc = mB2/mC2;
  theZc = mA2/mC2;

  // neutralino histograms
  string name[15] = { "qlnm", "qlnp", "qblnm", "qblnp", "qlfm", "qlfp",
		      "qblfm", "qblfp", "qlm", "qlp", "ll",
		      "W-p1", "W-p2", "W-lm", "W-lp" };
  for( unsigned int i = 0; i < 11; ++i ) {
    double limit(0.);
    if( i <= 3 ) limit = theMassNormalize ? 1.01 : mln_max;
    else if( i <= 7 && i >= 4 ) limit = theMassNormalize ? 1.01 : mlf_max;
    else if( i == 8 || i == 9 ) 
      limit = theMassNormalize ? 1.01 : max( mln_max, mlf_max);
    else limit = theMassNormalize ? 1.01 : mll_max;
    
    theHistograms[ name[i] ] = new_ptr( Histogram(0., limit, theNBins) );
  }
  //chargino histograms
  double limit(0.);
  if( theMassNormalize ) 
    limit = 2.;
  else {
    double x = mB2/sqr( getParticleData(theResonances[0])->mass() );
    Energy2 mqlmax = mB2*(1. - x)*
      ( (1. + theYc - theZc) 
	+ sqrt( sqr(1. + theYc - theZc) 
		- 4.*theZc ) )/2./x;
    limit = sqrt(mqlmax)/GeV;
  }
  for( unsigned int i = 11; i < 15; ++i ) 
    theHistograms[ name[i] ] = new_ptr( Histogram(0., limit, theNBins) );
    
  if( theMassNormalize ) {
    theMassMaxima[0] = mln_max*GeV;
    theMassMaxima[1] = mlf_max*GeV;
    theMassMaxima[2] = mll_max*GeV;
    theMassMaxima[3] = max( mln_max, mlf_max)*GeV;
    double x = mC2/sqr( getParticleData(theResonances[0] - 1)->mass() );
    theMassMaxima[4] = 0.5*sqrt( mC2*(1. - x)/x );
    x = mC2/sqr( getParticleData(theResonances[0])->mass() );
    theMassMaxima[5] = 0.5*sqrt( mC2*(1. - x)/x );
  }
  else 
    theMassMaxima.resize(6, 1.*GeV);
  
  //for charged distributions
  if( theModel == 0 ) {
    //need coupling factors 
    tMSSMPtr mssm = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
    MixingMatrix neut = *mssm->neutralinoMix();
    MixingMatrix uch = *mssm->charginoUMix();
    MixingMatrix vch  = *mssm->charginoVMix();
    unsigned int eic = ( theResonances[4] == 1000024 ) ? 0 : 1;
    Complex oL = neut(0,1)*conj(vch(eic,0)) - (neut(0,3)*conj(vch(eic,1))/sqrt(2));
    Complex oR = conj(neut(0,1))*uch(eic,0) + (conj(neut(0,2))*uch(eic,1)/sqrt(2));
    //common variables
    theAlpha = (oR - oL).real()/(oL + oR).real();
  }

}

void CascadeAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  
  //neutralino chain fractions
  long Nt = theNChain[0] + theNChain[1];
  theFractions[0] = (double)theNChain[0]/(double)Nt;
  theFractions[1] = (double)theNChain[1]/(double)Nt;

  //chargino/W chain fractions
  Nt = theNChain[2] + theNChain[3];
  theFractions[2] = (double)theNChain[2]/(double)Nt;
  theFractions[3] = (double)theNChain[3]/(double)Nt;
  Nt = theNChain[4] + theNChain[5];
  theFractions[4] = (double)theNChain[4]/(double)Nt;
  theFractions[5] = (double)theNChain[5]/(double)Nt;

  string filename = CurrentGenerator::current().filename();
  if( theModel == 0 ) filename += string("-mssm_spin.top");
  else filename += string("-mued_spin.top");
  
  ofstream os(filename.c_str());
  
  map<string, HistogramPtr>::const_iterator miter = theHistograms.begin();
  map<string, HistogramPtr>::const_iterator mend = theHistograms.end();
  using namespace HistogramOptions;
  for( ; miter != mend; ++miter ) {
    string title = (*miter).first;
    if( (*miter).second->globalStatistics().numberOfPoints() > 0 ) {
      (*miter).second->topdrawOutput(os, Frame, "BLUE", title);
      
      if( theMassNormalize ) {
	if( title.find("W") != string::npos ) {
	  for( double m = 0.0; m < 2.; m += 0.01 )
	    os << m << "\t\t\t" << chargedChain(title, m) << '\n';
	  os << "JOIN BLACK\n";
	}
	else {
	  for( double m = 0.0; m < 1.01; m += 0.01 )
	    os << m << "\t\t\t" << neutralChain(title, m) << '\n';
	  os << "JOIN BLACK\n";
	}
      }
    }
    
  }
  
}

double CascadeAnalysis::neutralChain(string s, double m) {
  if( m > 1. ) return 0.;
  if( theModel == 0 ) {
    if( s == "qlm" ) {
      double n1 = theMassMaxima[3]/theMassMaxima[0];
      double n2 = theMassMaxima[3]/theMassMaxima[1];
      double val = 
	theFractions[0]*( neutralChainSUSY("qlnm", m*n1)*n1 
			  + neutralChainSUSY("qlfm", m*n2)*n2 ) 
	+ theFractions[1]*( neutralChainSUSY("qblnm", m*n1)*n1 
			    + neutralChainSUSY("qblfm", m*n2)*n2 ) ;
      return 0.5*val;
    } 
    else if( s ==  "qlp" ) {
      double n1 = theMassMaxima[3]/theMassMaxima[0];
      double n2 = theMassMaxima[3]/theMassMaxima[1];
      double val = 
	theFractions[0]*( neutralChainSUSY("qlnp", m*n1)*n1 
			  + neutralChainSUSY("qlfp", m*n2)*n2 ) 
	+ theFractions[1]*( neutralChainSUSY("qblnp", m*n1)*n1 
			    + neutralChainSUSY("qblfp", m*n2)*n2 ) ;
      return 0.5*val;
    }
    else return neutralChainSUSY(s, m);
  }
  else {
    if( s == "qlm" ) {
      double n1 = theMassMaxima[3]/theMassMaxima[0];
      double n2 = theMassMaxima[3]/theMassMaxima[1];
      double val = 
	theFractions[0]*( neutralChainUED("qlnm", m*n1)*n1 
			  + neutralChainUED("qlfm", m*n2)*n2 ) 
	+ theFractions[1]*( neutralChainUED("qblnm", m*n1)*n1 
			    + neutralChainUED("qblfm", m*n2)*n2 ) ;
      return 0.5*val;
    } 
    else if( s ==  "qlp" ) {
      double n1 = theMassMaxima[3]/theMassMaxima[0];
      double n2 = theMassMaxima[3]/theMassMaxima[1];
      double val = 
	theFractions[0]*( neutralChainUED("qlnp", m*n1)*n1 
			  + neutralChainUED("qlfp", m*n2)*n2 ) 
	+ theFractions[1]*( neutralChainUED("qblnp", m*n1)*n1 
			    + neutralChainUED("qblfp", m*n2)*n2 ) ;
      return 0.5*val;
    }
    else return neutralChainUED(s, m);  
  }
  
}

double CascadeAnalysis::neutralChainSUSY(string s, double m) {
  if( m > 1. ) return 0.;

  double y(0.);
  //near process 1
  if( s == "qlnp" || s == "qblnm" ) y = 4.*pow(m, 3);
  //near process 2
  else if( s == "qlnm" || s == "qblnp" ) y = 4.*m*(1. - m*m);
  //far process 1
  else if( s == "qlfm" || s == "qblfp" ) {
    if(m <= sqrt(theY) )
      y = 1. - theY + log(theY);
    else
      y = 1. - m*m + 2.*log(m);
    
    y *= -4.*m/sqr(1. - theY);  
  }
  //far process 2
  else if( s == "qlfp" || s == "qblfm" ) {
    if( m <= sqrt(theY) )
      y = 1. - theY + theY*log(theY);
    else
      y = 1. - m*m + 2.*theY*log(m);
    
    y *= 4.*m/sqr(1. - theY);
  }
  else if( s == "ll" ) y = 2.*m;
  else y = 0.;
  
  return y;
}

double CascadeAnalysis::neutralChainUED(string s, double m) {
  if( m > 1. ) return 0.;
 
  double y(0.);
  //near process 1
  if( s == "qlnp" || s == "qblnm" ) {
    y = theY + 4.*(1. - theY + theX*theY)*m*m 
      - 4.*(1. - theX)*(1.- theY)*pow(m,4);
    y *= 6.*m/(1. + 2.*theX)/(2. + theY);
  }
  //near process 2
  else if( s == "qlnm" || s == "qblnp" ) {
    y = 4.*theX + theY + 4.*(1. - 2.*theX - theY + theX*theY)*m*m 
      - 4.*(1. - theX)*(1.- theY)*pow(m,4);
    y *= 6.*m/(1. + 2.*theX)/(2. + theY);
  }
  //far process 1
  else if(  s == "qlfm" || s == "qblfp" ) {
    if( m < sqrt(theY) ) {
      y = (1. - theY)*(4.*theX - theY 
		       + 2.*(2. + 3.*theY - 2.*theX*(5. + theY))*theZ
		    - 4.*m*m*(2. - 3.*theX)*(1. - 2.*theZ))
	- ( theY*(1. - 2.*(4. + theY)*theZ) 
	    + 4.*theX*(2.*theZ - theY*(1. - 4.*theZ))
	    + 4.*m*m*(1. + theY - theX*(2. + theY))*(1. - 2.*theZ) )*log(theY);
    }
    else {
      y = (1. - m*m)*( 4.*theX*(1. + 2.*theY - 5.*theZ - 6.*theY*theZ) - 5.*theY 
		       + 2.*(2. + 9.*theY)*theZ - 4.*m*m*(1. - theX)*(1. - theZ) )
	- 2.*( theY*(1. - 2.*(4. + theY)*theZ) 
	       + 4.*theX*(2.*theZ - theY*(1. - 4.*theZ))
	       + 4.*m*m*(1. + theY - theX*(2. + theY))*(1. - 2.*theZ) )*log(m);
    }
    y *= 12.*m/(1. + 2.*theX)/(2. + theY)/(1. + 2.*theZ)/sqr(1. - theY);
  }
  //far process 2
  else if(  s == "qlfp" || s == "qblfm" ) {
    if( m < sqrt(theY) ) {
      y = (1. - theY)*(2.*(2. + 2.*theX*(1. - theY) + 3.*theY)*theZ - theY
		      - 4.*m*m*(2. - theX)*(1. - 2.*theZ) )
	- ( theY*(1. - 2.*(4. + theY)*theZ) 
	    + 4.*m*m*(1. + (1. - theX)*theY)*(1. - 2.*theZ) )*log(theY);
    }
    else {
      y = (1. - m*m)*( 4.*(1. + theX)*theZ - theY*(5. - 18.*theZ + 8.*theX*theZ)
			 - 4.*m*m*(1. - theX)*(1 - theZ) )
	- 2.*( theY*(1. - 2.*(4. + theY)*theZ) 
	       + 4.*m*m*(1. + (1. - theX)*theY)*(1. - 2*theZ) )*log(m);
    
    }
    y *= 12.*m/(1. + 2.*theX)/(2. + theY)/(1. + 2.*theZ)/sqr(1. - theY);
  }
  else if( s == "ll" ) {
    y = theY + 4.*theZ + (2. - theY)*(1. - 2.*theZ)*sqr(m);
    y *= 4.*m/(2. + theY)/(1. + 2.*theZ);
  } 
  else y = 0.;

  return y;
}

double CascadeAnalysis::chargedChain(string s, double m) {
  if( theModel != 0 ) return 0.;

  double k1 = 1 + theYc - theZc;
  double k2 = sqrt(k1*k1 - 4.*theYc);
  double k1mk2 = k1 - k2;
  double k1pk2 = k1 + k2;

  if( sqr(m) > 2.*k1pk2 ) return 0.;
  double al2 = sqr(theAlpha);
  double brac = sqr(1 + theAlpha);
  double den = 32.*k2*( (1. + al2)*(k1*k1 + 2*theYc - 3.*k1*theYc) 
 			- 6.*theYc*sqrt(theZc)*(1. - al2) );

  double y(0.);
  //process-1
  if( s == "W-p1" ) {
    if( sqr(m) <= 2.*k1mk2 )
      y = 16.*( k2*( brac*k1 + 4.*brac*theYc 
		     + m*m*(1. + 4.*theAlpha + al2 + (al2 - 1.)*sqrt(theZc)) )
		- ( 2.*theAlpha*k1*m*m 
		    + (2.*brac*(1. + k1) 
		       + (1. + al2)*m*m)*theYc)*log(k1pk2/k1mk2) );
    else 
      y = 8.*k1*k1*brac + (1. + (theAlpha - 6.)*theAlpha)*pow(m, 4)
	+8.*k1*( brac*(k2 + 4.*theYc) 
		 + m*m*( sqrt(theZc) - 2. 
			 - theAlpha*(2. + theAlpha*(2 + sqrt(theZc))) ) )
	+ 8.*m*m*( k2*(1. + 4.*theAlpha + al2 - (1. - al2)*sqrt(theZc)) 
		   - 2.*brac*theYc )
	+ 16.*theYc*( brac*(5. + 2.*k2) - 4.*(1. - al2)*sqrt(theZc))
	- 16.*( 2.*theAlpha*k1*m*m 
		+ (2.*brac*(1. + k1) + (1. + al2)*m*m)*theYc)*log(2*k1pk2/m/m);
    
    y *= 3.*m/den;
  }
  //process-2
  else if( s == "W-p2" ) {
    if( sqr(m) <= 2.*k1mk2 )
      y = -16.*( k2*(m*m + 4*theYc - k1*sqr(1. - theAlpha) 
		     + theAlpha*(m*m*(4. + theAlpha) + 4.*theAlpha*theYc)
		     - m*m*sqrt(theZc)*(1. - al2) )
		 - (theYc*(2. + m*m - 4.*sqrt(theZc)) 
		    + 2.*theAlpha*(k1*m*m + 2*theYc)
		    + theYc*(2. + m*m + 4.*sqrt(theZc))*al2)*log(k1pk2/k1mk2));
    else
      y = 8.*k1*k1*sqr(1.-theAlpha) + pow(m,4)*(1. + theAlpha*(6. + theAlpha))
	-16.*theYc*(3. + 2.*k2 - 4.*sqrt(theZc) 
		+ theAlpha*(10. + theAlpha*(3. + 2.*k2 + 4.*sqrt(theZc))))
	+ 8.*k1*( k2*sqr(1 - theAlpha) 
		  + 2.*m*m*theAlpha - m*m*sqrt(theZc)*(1. - al2) 
		  - 4.*theYc*(1. + al2)) 
	- 8.*m*m*(-2*theYc*(1.+al2) + k2*(1. + 4.*theAlpha + al2 
					  + sqrt(theZc)*(-1. + al2)))
	+ 16.*(theYc*(2 + m*m - 4*sqrt(theZc)) + 2.*(k1*m*m + 2.*theYc)*theAlpha 
	       +theYc*(2. + m*m + 4.*sqrt(theZc))*al2)*log(2.*k1pk2/m/m); 

    y *= 3.*m/den;
  }
  //observable lepton
  else if( s == "W-lm" )
    y = theFractions[2]*chargedChain("W-p1", m) 
      + theFractions[3]*chargedChain("W-p2", m);
  //observable anti-lepton
  else if( s == "W-lp" )
    y = theFractions[4]*chargedChain("W-p2", m) 
      + theFractions[5]*chargedChain("W-p1", m);
  else return 0.;

  return y;
}
