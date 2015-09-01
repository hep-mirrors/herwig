// -*- C++ -*-
//
// UEDBase.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDBase class.
//

#include "UEDBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Repository/Repository.h" 
#include "ThePEG/Repository/CurrentGenerator.h"

using namespace Herwig;

UEDBase::UEDBase() : theRadCorr(true), theInvRadius(500.*GeV), 
		     theLambdaR(20.), theMbarH(), theSinThetaOne(0.),
		     theVeV(246.*GeV), includeSMMass_(true), fixedCouplings_(false), includeGaugeMixing_(true)
{}

void UEDBase::doinit() {
  readDecays(false);
  BSMModel::doinit();
  //level-1 masses and mixing angle
  calculateKKMasses(1);
  writeSpectrum();
  //add the level-1 vertices.
  addVertex(theF1F1Z0Vertex);
  addVertex(theF1F1G0Vertex);
  addVertex(theF1F0G1Vertex);
  addVertex(theG1G1G0Vertex);
  addVertex(theG0G0G1G1Vertex);
  addVertex(theF1F1P0Vertex);
  addVertex(theF1F1W0Vertex);
  addVertex(theF1F0W1Vertex);
  addVertex(theF1F0H1Vertex);
  addVertex(theP0H1H1Vertex);
  addVertex(theZ0H1H1Vertex);
  addVertex(theW0A1H1Vertex);
  addVertex(theZ0A1h1Vertex);
  addVertex(theW0W1W1Vertex);
  readDecays(true);
  if(decayFile()=="") return;
  decayRead();
}

void UEDBase::persistentOutput(PersistentOStream & os) const {
  os << theRadCorr << ounit(theInvRadius, GeV) << theLambdaR 
     << theF1F1Z0Vertex << theF1F1G0Vertex << theF1F0G1Vertex
     << theG1G1G0Vertex << theG0G0G1G1Vertex << theF1F1P0Vertex
     << theF1F1W0Vertex << theF1F0W1Vertex << theF1F0H1Vertex 
     << theP0H1H1Vertex << theZ0H1H1Vertex << theW0A1H1Vertex 
     << theZ0A1h1Vertex << theW0W1W1Vertex << ounit(theVeV,GeV) 
     << ounit(theMbarH, GeV) << theSinThetaOne << includeSMMass_
     << fixedCouplings_ << includeGaugeMixing_;
}

void UEDBase::persistentInput(PersistentIStream & is, int) {
  is >> theRadCorr >> iunit(theInvRadius, GeV) >> theLambdaR
     >> theF1F1Z0Vertex >> theF1F1G0Vertex >> theF1F0G1Vertex
     >> theG1G1G0Vertex >> theG0G0G1G1Vertex >> theF1F1P0Vertex
     >> theF1F1W0Vertex >> theF1F0W1Vertex >> theF1F0H1Vertex 
     >> theP0H1H1Vertex >> theZ0H1H1Vertex >> theW0A1H1Vertex 
     >> theZ0A1h1Vertex >> theW0W1W1Vertex >> iunit(theVeV,GeV) 
     >> iunit(theMbarH, GeV) >> theSinThetaOne >> includeSMMass_
     >> fixedCouplings_ >> includeGaugeMixing_;
}

ClassDescription<UEDBase> UEDBase::initUEDBase;
// Definition of the static class description member.

void UEDBase::Init() {

  static ClassDocumentation<UEDBase> documentation
    ("This class implements/stores the necessary information for the simulation"
     " of a Universal Extra Dimensions model.",
     "Universal extra dimensions model based on \\cite{Cheng:2002iz,Appelquist:2000nn}.",
     "%\\cite{Cheng:2002iz}\n"
     "\\bibitem{Cheng:2002iz}\n"
     "  H.~C.~Cheng, K.~T.~Matchev and M.~Schmaltz,\n"
     "  ``Radiative corrections to Kaluza-Klein masses,''\n"
     "  Phys.\\ Rev.\\  D {\\bf 66}, 036005 (2002)\n"
     "  [arXiv:hep-ph/0204342].\n"
     "  %%CITATION = PHRVA,D66,036005;%%\n"
     "%\\cite{Appelquist:2000nn}\n"
     "\\bibitem{Appelquist:2000nn}\n"
     "  T.~Appelquist, H.~C.~Cheng and B.~A.~Dobrescu,\n"
     "  ``Bounds on universal extra dimensions,''\n"
     "  Phys.\\ Rev.\\  D {\\bf 64}, 035002 (2001)\n"
     "  [arXiv:hep-ph/0012100].\n"
     "  %%CITATION = PHRVA,D64,035002;%%\n"
   );

  static Switch<UEDBase,bool> interfaceRadiativeCorrections
    ("RadiativeCorrections",
     "Calculate the radiative corrections to the masses",
     &UEDBase::theRadCorr, true, false, false);
  static SwitchOption interfaceRadiativeCorrectionsYes
    (interfaceRadiativeCorrections,
     "Yes",
     "Calculate the radiative corrections to the masses",
     true);
  static SwitchOption interfaceRadiativeCorrectionsNo
    (interfaceRadiativeCorrections,
     "No",
     "Leave the masses of the KK particles as n/R",
     false);

  static Parameter<UEDBase,Energy> interfaceInverseRadius
    ("InverseRadius",
     "The inverse radius of the compactified dimension ",
     &UEDBase::theInvRadius, GeV, 500.*GeV, ZERO, ZERO,
     true, false, Interface::nolimits);

  static Parameter<UEDBase,double> interfaceLambdaR
    ("LambdaR",
     "The product of the cut-off scale  and the radius of compactification",
     &UEDBase::theLambdaR, 20.0, 0.0, 0,
     false, false, Interface::lowerlim);

    static Parameter<UEDBase,Energy> interfaceBoundaryMass
    ("HiggsBoundaryMass",
     "The boundary mass for the Higgs",
     &UEDBase::theMbarH, GeV, ZERO, ZERO, ZERO,
     false, false, Interface::lowerlim);

  static Parameter<UEDBase,Energy> interfaceVeV
    ("HiggsVEV",
     "The vacuum expectation value of the Higgs field",
     &UEDBase::theVeV, GeV, 246.*GeV, ZERO, ZERO,
     true, false, Interface::nolimits);
    
  static Reference<UEDBase,Helicity::AbstractFFVVertex> interfaceF1F1Z
    ("Vertex/F1F1Z",
     "The F1F1Z UED Vertex",
     &UEDBase::theF1F1Z0Vertex, false, false, true, false, false);

  static Reference<UEDBase,Helicity::AbstractFFVVertex> interfaceF1F1G0
    ("Vertex/F1F1G0",
     "The F1F1G UED Vertex",
     &UEDBase::theF1F1G0Vertex, false, false, true, false, false);

  static Reference<UEDBase,Helicity::AbstractFFVVertex> interfaceF1F0G1
    ("Vertex/F1F0G1",
     "The F1F0G0 UED Vertex",
     &UEDBase::theF1F0G1Vertex, false, false, true, false, false);

  static Reference<UEDBase,Helicity::AbstractVVVVertex> interfaceG1G1G0
    ("Vertex/G1G1G0",
     "The G1G1G0 UED Vertex",
     &UEDBase::theG1G1G0Vertex, false, false, true, false, false);

  static Reference<UEDBase,Helicity::AbstractVVVVVertex> interfaceG0G0G1G1
    ("Vertex/G0G0G1G1",
     "The G0G0G1G1 UED Vertex",
     &UEDBase::theG0G0G1G1Vertex, false, false, true, false, false);

  static Reference<UEDBase,Helicity::AbstractFFVVertex> interfaceF1F1P
    ("Vertex/F1F1P",
     "The F1F1P UED Vertex",
     &UEDBase::theF1F1P0Vertex, false, false, true, false, false);

  static Reference<UEDBase,Helicity::AbstractFFVVertex> interfaceF1F1W
    ("Vertex/F1F1W",
     "The F1F1W UED Vertex",
     &UEDBase::theF1F1W0Vertex, false, false, true, false, false);

  static Reference<UEDBase,Helicity::AbstractFFVVertex> interfaceF1F0W1
    ("Vertex/F1F0W1",
     "The F1F0W1 UED Vertex",
     &UEDBase::theF1F0W1Vertex, false, false, true, false, false);

  static Reference<UEDBase,Helicity::AbstractFFSVertex> interfaceF1F0H1
    ("Vertex/F1F0H1",
     "The F1F0H1 UED Vertex",
     &UEDBase::theF1F0H1Vertex, false, false, true, false, false);

  static Reference<UEDBase,Helicity::AbstractVSSVertex> interfaceP0H1H1
    ("Vertex/P0H1H1",
     "The P0H1H1 UED Vertex",
     &UEDBase::theP0H1H1Vertex, false, false, true, false, false);

  static Reference<UEDBase,Helicity::AbstractVSSVertex> interfaceZ0H1H1
    ("Vertex/Z0H1H1",
     "The Z0H1H1 UED Vertex",
     &UEDBase::theZ0H1H1Vertex, false, false, true, false, false);

  static Reference<UEDBase,Helicity::AbstractVSSVertex> interfaceW0A1H1
    ("Vertex/W0A1H1",
     "The W0A1H1 UED Vertex",
     &UEDBase::theW0A1H1Vertex, false, false, true, false, false);

  static Reference<UEDBase,Helicity::AbstractVSSVertex> interfaceZ0A1h1
    ("Vertex/Z0A1h1",
     "The W0A1H1 UED Vertex",
     &UEDBase::theZ0A1h1Vertex, false, false, true, false, false);

  static Reference<UEDBase,Helicity::AbstractVVVVertex> interfaceW0W1W1
    ("Vertex/W0W1W1",
     "The W0W1W1 UED Vertex",
     &UEDBase::theW0W1W1Vertex, false, false, true, false, false);

  static Switch<UEDBase,bool> interfaceIncludeSMMass
    ("IncludeSMMass",
     "Whether or not to include the SM mass in the calculation of the masses of the KK states.",
     &UEDBase::includeSMMass_, true, false, false);
  static SwitchOption interfaceIncludeSMMassYes
    (interfaceIncludeSMMass,
     "Yes",
     "Include them",
     true);
  static SwitchOption interfaceIncludeSMMassNo
    (interfaceIncludeSMMass,
     "No",
     "Don't include them",
     false);

  static Switch<UEDBase,bool> interfaceFixedCouplings
    ("FixedCouplings",
     "Use fixed or running couplings to calculate the masses.",
     &UEDBase::fixedCouplings_, false, false, false);
  static SwitchOption interfaceFixedCouplingsYes
    (interfaceFixedCouplings,
     "Yes",
     "Use fixed couplings",
     true);
  static SwitchOption interfaceFixedCouplingsNo
    (interfaceFixedCouplings,
     "No",
     "Use running couplings",
     false);

  static Switch<UEDBase,bool> interfaceIncludeGaugeMixing
    ("IncludeGaugeMixing",
     "Whether or not to include mixing between the KK photon"
     " and Z in the vertices, always included in the mass",
     &UEDBase::includeGaugeMixing_, true, false, false);
  static SwitchOption interfaceIncludeGaugeMixingYes
    (interfaceIncludeGaugeMixing,
     "Yes",
     "Include the mixing",
     true);
  static SwitchOption interfaceIncludeGaugeMixingNo
    (interfaceIncludeGaugeMixing,
     "No",
     "Don't include the mixing",
     false);

}

void UEDBase::calculateKKMasses(const unsigned int n) {
  useMe();
  if(n == 0)
    throw InitException() << "UEDBase::resetKKMasses - "
			  << "Trying to reset masses with KK number == 0!"
			  << Exception::warning;
    if(theRadCorr) {
      fermionMasses(n);
      bosonMasses(n);
    }
    else {
      cerr << 
	"Warning: Radiative corrections to particle masses have been "
	"turned off.\n  The masses will be set to (n/R + m_sm)^1/2 and "
	"the spectrum will be\n  highly degenerate so that no decays "
	"will occur.\n  This is only meant to be used for debugging "
	"purposes.\n";
      //set masses to tree level for each kk mode
      long level1 = 5000000 + n*100000;
      long level2 = 6000000 + n*100000;
      Energy2 ndmass2 = sqr(n*theInvRadius);
      for ( int i = 1; i < 38; ++i ) {
	if(i == 7 || i == 17) i += 4;
	if(i == 26) i += 10;
	Energy kkmass = sqrt( ndmass2 + sqr(getParticleData(i)->mass()) );
	resetMass(level1 + i, kkmass);
	if( i < 7 || i == 11 || i == 13 || i == 15 )
	  resetMass(level2 + i, kkmass);
      }
    }

}

void UEDBase::bosonMasses(const unsigned int n) {
  // Common constants
  const Energy2 invRad2 = theInvRadius*theInvRadius;
  const double g_em2 = fixedCouplings_ ? 
    4.*Constants::pi*alphaEMMZ() : 4.*Constants::pi*alphaEM(invRad2);
  const double g_s2 = fixedCouplings_ ? 
    4.*Constants::pi*alphaS()    : 4.*Constants::pi*alphaS(invRad2);
  const double g_W2 = g_em2/sin2ThetaW();
  const double g_P2 = g_em2/(1-sin2ThetaW());
  //Should probably use a function  to calculate zeta.
  const double zeta3 = 1.20206;
  const Energy2 nmass2 = sqr(n*theInvRadius); 
  const double pi2 = sqr(Constants::pi);
  const double norm = 1./16./pi2;
  const double nnlogLR = n*n*log(theLambdaR);
  long level = 5000000 + n*100000;
  //gluon
  Energy2 deltaGB = g_s2*invRad2*norm*(23.*nnlogLR - 3.*zeta3/2./pi2 );
  resetMass(level + 21, sqrt(nmass2 + deltaGB));

  //W+/-
  Energy2 deltaGW = g_W2*invRad2*norm*( 15.*nnlogLR - 5.*zeta3/2./pi2 );

  Energy2 mw2 = sqr(getParticleData(24)->mass());
  resetMass(level + 24, sqrt(mw2 + nmass2 + deltaGW));

  //Z and gamma are a mixture of Bn and W3n
  deltaGB = -g_P2*invRad2*norm*( 39.*zeta3/2./pi2 + nnlogLR/3. );
  Energy2 mz2 = sqr(getParticleData(23)->mass());
  Energy2 fp = 0.5*(mz2 + deltaGB + deltaGW + 2.*nmass2);
  Energy2 sp = 0.5*sqrt( sqr(deltaGB - deltaGW - 2.*mw2 + mz2)
			 - 4.*mw2*(mw2 - mz2) );
  resetMass(level + 22, sqrt(fp - sp));
  resetMass(level + 23, sqrt(fp + sp));
  //mixing angle will now depend on both Z* and gamma* mass
  //Derived expression:
  // 
  // cos^2_theta_N = ( (n/R)^2 + delta_GW + mw^2 - m_gam*^2)/(m_z*^2 - m_gam*^2)
  //
  if(includeGaugeMixing_) {
    double cn2 = (nmass2 + deltaGW + mw2 - fp + sp)/2./sp;
    double sn = sqrt(1. - cn2);
    theMixingAngles.insert(make_pair(n, sn));
    if( n == 1 ) theSinThetaOne = sn;
  }
  else {
    theMixingAngles.insert(make_pair(n,0.));
    if( n == 1 ) theSinThetaOne = 0.;
  }
  //scalars
  Energy2 mh2 = sqr(getParticleData(25)->mass());
  double lambda_H = mh2/sqr(theVeV);
  deltaGB = nnlogLR*norm*invRad2*(3.*g_W2 + (3.*g_P2/2.) - 2.*lambda_H) 
    + sqr(theMbarH);
  //H0
  Energy2 new_m2 = nmass2 + deltaGB;
  resetMass(level + 25, sqrt( mh2 + new_m2 ));
  //A0
  resetMass(level + 36, sqrt( mz2 + new_m2 ));
  //H+
  resetMass(level + 37, sqrt( mw2 + new_m2 ));
}

void UEDBase::fermionMasses(const unsigned int n) {
  const Energy2 invRad2 = theInvRadius*theInvRadius;
  const double g_em2 = fixedCouplings_ ? 
    4.*Constants::pi*alphaEMMZ() : 4.*Constants::pi*alphaEM(invRad2);
  const double g_s2 = fixedCouplings_ ? 
    4.*Constants::pi*alphaS() : 4.*Constants::pi*alphaS(invRad2); 
  const double g_W2 = g_em2/sin2ThetaW();
  const double g_P2 = g_em2/(1-sin2ThetaW());
  const Energy nmass = n*theInvRadius;
  const Energy norm = 
    nmass*log(theLambdaR)/16./Constants::pi/Constants::pi;
  const Energy topMass = getParticleData(6)->mass();
  const double ht = sqrt(2)*topMass/theVeV;
  //doublets
  Energy deltaL = norm*(6.*g_s2 + (27.*g_W2/8.) + (g_P2/8.));
  Energy deltaQ = deltaL;
  Energy2 shift = sqr(nmass + deltaL);
  long level = 5000000 + n*100000;
  for(long i = 1; i < 17; ++i) {
    if(i == 5)  {
      i += 6;
      deltaL = norm*( (27.*g_W2/8.) + (9.*g_P2/8.) );
      shift = sqr(nmass + deltaL); 
    }
    Energy2 new_m2 = includeSMMass_ ? sqr(getParticleData(i)->mass()) + shift : shift;
    resetMass(level + i, sqrt(new_m2));
  }
  //singlet shifts
  const Energy  deltaU = norm*(6.*g_s2 + 2.*g_P2);
  const Energy  deltaD = norm*(6.*g_s2 + 0.5*g_P2);
  const Energy2 shiftU = sqr(nmass + deltaU);
  const Energy2 shiftD = sqr(nmass + deltaD);
  
  //Top quarks seperately as they have different corrections
  const Energy2 mt2 = sqr(topMass);
  const Energy delta_Q3 = -3.*ht*ht*norm/2.;
  const Energy deltaTD = deltaQ + delta_Q3;
  const Energy deltaTS = deltaU + 2.*delta_Q3;
  Energy second_term = 0.5*sqrt( sqr(2.*nmass + deltaTS + deltaTD) + 4.*mt2 );
  //doublet
  resetMass(level           + 6, abs(0.5*(deltaTD - deltaTS) - second_term) );
  //singlet
  resetMass(level + 1000000 + 6, 0.5*(deltaTD - deltaTS) + second_term);
  //Bottom quarks
  const Energy2 mb2 = sqr(getParticleData(5)->mass());
  const Energy deltaBS = deltaD;
  second_term = 0.5*sqrt( sqr(2.*nmass + deltaBS + deltaTD) + 4.*mb2 );
  //doublet
  resetMass(level + 1000000 + 5, abs(0.5*(deltaTD - deltaBS) - second_term) );
  //singlet
  resetMass(level           + 5, 0.5*(deltaTD - deltaBS) + second_term);
  // others
  //lepton 
  Energy delta = 9.*norm*g_P2/2.;
  shift = sqr(nmass + delta);

  level += 1000000;
  for(long i = 1; i < 17; ) {
    if(i == 5) i += 6;
    Energy2 smMass2(sqr(getParticleData(i)->mass()));
    if(i < 6) {
      Energy2 new_m2 = includeSMMass_ ? smMass2 : ZERO;
      if( i % 2 == 0) new_m2 = shiftU;
      else            new_m2 = shiftD;
      resetMass(level + i, sqrt(new_m2));
      ++i;
    }
    else {
      if(includeSMMass_)
	resetMass(level + i, sqrt(smMass2 + shift));
      else
	resetMass(level + i, sqrt(shift));
      i += 2;
    }
  }
}

void UEDBase::resetMass(long id, Energy mass) {
  theMasses.push_back(make_pair(id, mass));
  StandardModel::resetMass(id,mass);
}

void UEDBase::writeSpectrum() {
  sort(theMasses.begin(), theMasses.end(), lowerMass);
  ostream & ofs = CurrentGenerator::current().misc();
  ofs << "# MUED Model Particle Spectrum\n"
      << "# R^-1: " << theInvRadius/GeV << " GeV\n"
      << "# Lambda * R: " << theLambdaR << "\n"
      << "# Higgs Mass: " << getParticleData(25)->mass()/GeV << " GeV\n";
  ofs << "#\n# ID\t\t\tMass(GeV)\n";
  while (!theMasses.empty()) {
    IDMassPair tmp = theMasses.back();
    tcPDPtr data = getParticleData(tmp.first);
    ofs << tmp.first << "\t\t\t" << tmp.second/GeV << "\t\t" << (data? data->PDGName() : "") 
	<< endl;
    theMasses.pop_back();
  }
  ofs << "#\n";
}

double UEDBase::sinThetaN(const unsigned int n) const {
  WAMap::const_iterator pos = theMixingAngles.find(n);
  if(pos != theMixingAngles.end())
    return pos->second;
  else {
    throw Exception() << "UEDBase::sinThetaN() - A mixing angle has "
		      << "been requested for a level that does not "
		      << "exist. Check that the radiative corrections "
		      << "for the " << n << "th level have been "
		      << "calculated." << Exception::warning;
    return 0.0;
  }
}
