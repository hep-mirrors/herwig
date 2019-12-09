// -*- C++ -*-
//
// MEHiggsPair.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2009-2019 The Herwig Collaboration
//
// Herwig is licenaaaced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEHiggsPair class.
//

#include "MEHiggsPair.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Cuts/Cuts.h"
#include "Herwig/Utilities/Maths.h"
#include <fstream>
#include <cmath>

namespace Herwig {


}

using namespace Herwig;


MEHiggsPair::MEHiggsPair()
  : _selfcoupling(1.0),_hhHcoupling(1.0), _process(0), _mh(), _wh(), _cH(0.0), _c6(0.0), _ct1(0.0),  _cb1(0.0), _ct2(0.0),  _cb2(0.0), _cg1(0.0), _cg2(0.0), _EFTScale(1000*GeV),  _fixedalphaS(0), _alphasfixedvalue(0.1), _alphascale(90.0*GeV), _basescale(125.0*GeV), _fixedscale(0), _scalemultiplier(1.0), _yH(1.0), _yh(1.0), _ybH(1.0), _ybh(1.0)   {
  massOption(vector<unsigned int>(2,0));
}

int MEHiggsPair::nDim() const {
  return 5;
}


void MEHiggsPair::doinit() {

  if(_process == 5 || _process == 6) { cerr << "HH and hH production not implemented yet, please choose an hh production subprocess" << endl; exit(1); } 
  if(_process == 4) { cerr << "WARNING: In the EFT, the option SelfCoupling is not used." << endl; }

  higgs(getParticleData(ParticleID::h0));
  PDPtr wboson = getParticleData(ParticleID::Wplus);
  _Wmass = wboson->mass();
  PDPtr top = getParticleData(ParticleID::t);
  _topmass = top->mass();
  PDPtr zboson = getParticleData(ParticleID::Z0);
  _zmass = zboson->mass();
  PDPtr bottom = getParticleData(ParticleID::b);
  _bottommass = bottom->mass();
  PDPtr heavyH = getParticleData(35);
  _heavyHmass = heavyH->mass();
  _heavyHwidth = heavyH->width();

  // higgs stuff
  _mh = _higgs->mass();
  _m1 = _higgs->mass();
  _m2 = _higgs->mass();
  _wh = _higgs->width();
  if(_higgs->massGenerator()) {
    _hmass=dynamic_ptr_cast<GenericMassGeneratorPtr>(_higgs->massGenerator());
  }
  HwMEBase::doinit();
}

void MEHiggsPair::rebind(const TranslationMap & trans) {
  HwMEBase::rebind(trans);
}

IVector MEHiggsPair::getReferences() {
  IVector ret = HwMEBase::getReferences();
  return ret;
}

IBPtr MEHiggsPair::clone() const {
    return new_ptr(*this);
  }
  
IBPtr MEHiggsPair::fullclone() const {
    return new_ptr(*this);
  }

void MEHiggsPair::persistentOutput(PersistentOStream & os) const {
  os << _selfcoupling << _hhHcoupling << _process << _higgs << ounit(_topmass,GeV) << ounit(_bottommass,GeV) << ounit(_zmass,GeV) << ounit(_m1,GeV) << ounit(_m2,GeV) << ounit(_mh,GeV) << ounit(_wh,GeV) << _hmass << _cH << _c6 << _ct1 << _cb1 << _ct2 << _cb2 << _cg1 << _cg2 << ounit(_EFTScale,GeV) << _fixedalphaS << _alphasfixedvalue << ounit(_heavyHmass,GeV) << ounit(_heavyHwidth,GeV) << ounit(_basescale,GeV) << ounit(_alphascale,GeV) << _fixedscale << _scalemultiplier << _yH << _yh << _ybH << _ybh;
}

void MEHiggsPair::persistentInput(PersistentIStream & is, int) {
  is >> _selfcoupling >> _hhHcoupling >> _process >> _higgs >> iunit(_topmass,GeV) >> iunit(_bottommass,GeV) >> iunit(_zmass,GeV) >> iunit(_m1,GeV) >> iunit(_m2,GeV) >> iunit(_mh,GeV) >> iunit(_wh,GeV) >> _hmass >> _cH >> _c6 >> _ct1 >> _cb1 >> _ct2 >> _cb2 >> _cg1 >> _cg2 >> iunit(_EFTScale,GeV) >>  _fixedalphaS >> _alphasfixedvalue >> iunit(_heavyHmass,GeV) >> iunit(_heavyHwidth,GeV) >> iunit(_basescale,GeV) >> iunit(_alphascale,GeV) >> _fixedscale >> _scalemultiplier >> _yH >> _yh >> _ybH >> _ybh;

}

Energy2 MEHiggsPair::scale() const {
  if(_fixedscale==0 || _fixedscale==5 ) { return sqr(_scalemultiplier)*sHat(); } 
  else if (_fixedscale==1 || _fixedscale==2 || _fixedscale==3) { return (sqr(_scalemultiplier)*sqr(_basescale)); }
}

// Definition of the static class description member.
ClassDescription<MEHiggsPair> MEHiggsPair::initMEHiggsPair;

void MEHiggsPair::Init() {

  static ClassDocumentation<MEHiggsPair> documentation
    ("The MEHiggsPair class implements the gg -> hh processes in hadron-hadron"
     " collisions");

 static Parameter<MEHiggsPair, double> interfaceSelfCoupling
    ("SelfCoupling",
     "Multiplier for the SM Higgs triple coupling",
     &MEHiggsPair::_selfcoupling, 1.0, -1000000.0, 1000000.0,
     false, false, Interface::limited);

 static Parameter<MEHiggsPair, double> interfacecH
    ("cH",
     "Effective theory coefficient cH",
     &MEHiggsPair::_cH, 0.0, -1000000.0, 1000000.0,
     false, false, Interface::limited);

 static Parameter<MEHiggsPair, double> interfacec6
    ("c6",
     "Effective theory coefficient c6",
     &MEHiggsPair::_c6, 0.0, -1000000.0, 1000000.0,
     false, false, Interface::limited);

 static Parameter<MEHiggsPair, double> interfacect1
    ("ct1",
     "Effective theory coefficient ct1",
     &MEHiggsPair::_ct1, 0.0, -1000000.0, 1000000.0,
     false, false, Interface::limited);

 static Parameter<MEHiggsPair, double> interfacecb1
    ("cb1",
     "Effective theory coefficient fb1",
     &MEHiggsPair::_cb1, 0.0, -1000000.0, 1000000.0,
     false, false, Interface::limited);

 static Parameter<MEHiggsPair, double> interfacect2
    ("ct2",
     "Effective theory coefficient ft2",
     &MEHiggsPair::_ct2, 0.0, -1000000.0, 1000000.0,
     false, false, Interface::limited);

 static Parameter<MEHiggsPair, double> interfacecb2
    ("cb2",
     "Effective theory coefficient cb2",
     &MEHiggsPair::_cb2, 0.0, -1000000.0, 1000000.0,
     false, false, Interface::limited);

 static Parameter<MEHiggsPair, double> interfacecg1
    ("cg1",
     "Effective theory coefficient cg1",
     &MEHiggsPair::_cg1, 0.0, -1000000.0, 1000000.0,
     false, false, Interface::limited);

static Parameter<MEHiggsPair, double> interfacecg2
    ("cg2",
     "Effective theory coefficient cg2",
     &MEHiggsPair::_cg2, 0.0, -1000000.0, 1000000.0,
     false, false, Interface::limited);

  static Parameter<MEHiggsPair, double> interfacecyytHeavy
    ("ytH",
     "Multiplier for Heavy Higgs-top Yukawa",
     &MEHiggsPair::_yH, 1.0, -1000000.0, 1000000.0,
     false, false, Interface::limited);

  static Parameter<MEHiggsPair, double> interfacecyytLight
    ("yth",
     "Multiplier for Light Higgs-top Yukawa",
     &MEHiggsPair::_yh, 1.0, -1000000.0, 1000000.0,
     false, false, Interface::limited);

   static Parameter<MEHiggsPair, double> interfacecyybHeavy
    ("ybH",
     "Multiplier for Heavy Higgs-top Yukawa",
     &MEHiggsPair::_ybH, 1.0, -1000000.0, 1000000.0,
     false, false, Interface::limited);

  static Parameter<MEHiggsPair, double> interfacecyybLight
    ("ybh",
     "Multiplier for Light Higgs-top Yukawa",
     &MEHiggsPair::_ybh, 1.0, -1000000.0, 1000000.0,
     false, false, Interface::limited);



 static Parameter<MEHiggsPair, Energy> interfaceEFTScale
    ("EFTScale",
     "Effective theory coefficient f1",
     &MEHiggsPair::_EFTScale, GeV, 1000.0*GeV, 1.0*GeV, 10000000*GeV,
     false, false, Interface::limited);


  static Switch<MEHiggsPair,unsigned int> interfaceFixedAlphaS
    ("FixedAlphaS",
     "Which implementation to use",
     &MEHiggsPair::_fixedalphaS, 0, false, false);
  static SwitchOption interfaceFixedAlphaSNo
    (interfaceFixedAlphaS,
     "No",
     "Don't fix the scale of alphaS.",
     0);
  static SwitchOption interfaceFixedAlphaSScale
    (interfaceFixedAlphaS,
     "Scale",
     "Fix the scale of alphaS.",
     1);
  static SwitchOption interfaceFixedAlphaSAlphaS
    (interfaceFixedAlphaS,
     "AlphaS",
     "Fix alphaS directly via the AlphaS interface.",
     2);


  static Parameter<MEHiggsPair, Energy> interfaceBaseScale
    ("BaseScale",
     "Base scale if fixed",
     &MEHiggsPair::_basescale, GeV, 125.0*GeV, 0.0*GeV, 10000000.0*GeV,
     false, false, Interface::limited);

 
  static Switch<MEHiggsPair,unsigned int> interfaceFixedScale
    ("FixedScale",
     "Whether to use fixed base scale or not",
     &MEHiggsPair::_fixedscale, 0, false, false);
  static SwitchOption interfaceFixedScaleNo
    (interfaceFixedScale,
     "sHat",
     "Don't fix the base scale, use sHat().",
     0);
  static SwitchOption interfaceFixedScaleFixed
    (interfaceFixedScale,
     "Fixed",
     "Fix the scale of the process to chosen (base scale)^2.",
     2);
  static SwitchOption interfaceFixedScaleBaseQuad
    (interfaceFixedScale,
     "FixedQuadPt",
     "Fix the scale of the process to (base scale)^2 + pt^2.",
     1);
  static SwitchOption interfaceFixedScaleBaseLin
    (interfaceFixedScale,
     "FixedLinPt",
     "Fix the scale of the process to (base scale + pt)^2.",
     3);
  static SwitchOption interfaceFixedScaleMHH
    (interfaceFixedScale,
     "MHH",
     "Don't fix the base scale, use MHH.",
     5);

  static Parameter<MEHiggsPair, double> interfaceScaleMultiplier
    ("ScaleMultiplier",
     "Multiplier scale",
     &MEHiggsPair::_scalemultiplier, 1.0, 0.0, 100000.0,
     false, false, Interface::limited);

  static Parameter<MEHiggsPair, Energy> interfaceAlphaSScale
    ("AlphaSScale",
     "Scale for alphaS if fixed",
     &MEHiggsPair::_alphascale, GeV, 100.0*GeV, 0.0*GeV, 10000000.0*GeV,
     false, false, Interface::limited);

  static Parameter<MEHiggsPair, double> interfacehHHCoupling
    ("hhHCoupling",
     "Multiplier for the hh-H triple coupling",
     &MEHiggsPair::_hhHcoupling, 1.0, -10.0, 10.0,
     false, false, Interface::limited);
 
 //choose whether to include the triangle, box or both
  static Switch<MEHiggsPair,unsigned int> interfaceProcess
    ("Process",
     "Which subprocesses to include",
     &MEHiggsPair::_process, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Include all SM gg->hh subprocesses",
     0);
  static SwitchOption interfaceProcessTriangle
    (interfaceProcess,
     "ggToTriangleTohh",
     "Include only SM gg->hh triangle subprocess",
     1);
  static SwitchOption interfaceProcessBox
    (interfaceProcess,
     "ggToBoxTohh",
     "Include only SM gg->hh box subprocess",
     2);
  static SwitchOption interfaceProcessHhh
    (interfaceProcess,
     "ggToHTohh",
     "Include all gg->hh subprocess, with heavy Higgs",
     3);
  static SwitchOption interfaceProcessEFT
    (interfaceProcess,
     "EFT",
     "Include all SM gg->hh subprocesses in the EFT",
     4);
  static SwitchOption interfaceProcessHH
    (interfaceProcess,
     "ggToHH",
     "Include all gg->HH subprocesses",
     5);
  static SwitchOption interfaceProcessHh
    (interfaceProcess,
     "ggToHh",
     "Include only gg->hH subprocesses",
     6);

 

 static Parameter<MEHiggsPair, double> interfaceAlphaS
    ("AlphaS",
     "The value of AlphaS if FixedAlphaS is AlphaS, option 1",
     &MEHiggsPair::_alphasfixedvalue, 0.1, 0., 100000.0,
     false, false, Interface::limited);

 
}

CrossSection MEHiggsPair::dSigHatDR() const {
  using Constants::pi;
  return me2()/(16.0*pi*sqr(sHat())/GeV/GeV)*sqr(hbarc);
}



Selector<MEBase::DiagramIndex>
MEHiggsPair::diagrams(const DiagramVector & diags) const {
  // select the diagram, this is easy for us as we have already done it
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    sel.insert(1.0, i);
  }
  return sel;
}

//diagrams with box or triangle
void MEHiggsPair::getDiagrams() const {
  // get the particle data objects
  PDPtr gluon = getParticleData(ParticleID::g);
  PDPtr boxon = getParticleData(99926);
  PDPtr triangon = getParticleData(99927);
  //cout<<  _higgs->mass() << endl;
  if(_process==0||_process==1||_process==3||_process==4) {  
    add(new_ptr((Tree2toNDiagram(2),gluon,gluon,
		 1,triangon,3,higgs(),higgs(),-1)));    
  }
  if(_process==0||_process==2||_process==3||_process==4) {
    add(new_ptr((Tree2toNDiagram(2),gluon,gluon,
		 1,boxon,3,higgs(),higgs(),-2)));    
  }
   
}

//both diagrams possess the same colour structure
Selector<const ColourLines *>
MEHiggsPair::colourGeometries(tcDiagPtr diag) const {
  // colour lines for gg to hh
  static const ColourLines cgghh("1 -2, -1 2");

  Selector<const ColourLines *> sel;
  switch(abs(diag->id())) {
  case 1: 
    sel.insert(1.0, &cgghh);
    break;
  case 2:
    sel.insert(1.0, &cgghh);
    break;
  }
  return sel;
}

//the matrix element squared with the appropriate factors
double MEHiggsPair::me2() const {
  using Constants::pi;
  double mesq(0.);

  //get the masses and Mandestam variables
  Energy mh(_mh);
  double Mc = _mh/GeV;
  double Md = _mh/GeV;
  double s = sHat()/GeV/GeV;
  double t = tHat()/GeV/GeV;
  double u = uHat()/GeV/GeV;

  //calculate the matrix element from the form factors
  mesq = MATRIX(s, t, u, Mc, Md);
  //QCD factors
  double ALPS=SM().alphaS(scale());
  if(_fixedalphaS==1) { ALPS = _alphasfixedvalue; }
  double QCD=sqr(ALPS);
  //XLAM = beta*S
  double XLAM=sqrt(sqr(s+sqr(Mc)-sqr(Md))-4.*s*sqr(Mc));
  //the jacobian is jacobian()*XLAM
  mesq *= jacobian()*XLAM; 
  mesq *= QCD / (2048.*sqr(pi));
  return mesq;    
}

bool MEHiggsPair::generateKinematics(const double * r) { 
  //C--FUNCTION FOR PARTONIC CROSS SECTION
  using Constants::pi;

  //the form factors
  Complex A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2;

  //phase space cut, masses, and Mandelstam variables
  double EPS=1E-8;
  double S = sHat()/GeV/GeV;
  double T = tHat()/GeV/GeV;
  double U = uHat()/GeV/GeV;
  double M1 = _m1/GeV;
  double M2 = _m2/GeV;
  double mh = M1; 
  double WHAT = sqrt(S);

  //generate the random variables
  vector<double> Y;
  int N=2;
  Y.resize(N);
  for(int I=0; I < N; I++) { 
    Y[I]=EPS+(1.-2.*EPS)*r[I];
  }
  jacobian(1.);
  double DJAC = 1-2*EPS;
  double V=Y[0]/2.;

  double XLAM=sqrt(sqr(S+sqr(M1)-sqr(M2))-4.*S*sqr(M1));
  double BET=XLAM/(S+sqr(M1)-sqr(M2));
  //sample T1 and U1
  double T1=-(S+sqr(M1)-sqr(M2))*((1.+BET)/2.-BET*V);
  double U1=-(S+sqr(M1)-sqr(M2))*((1.-BET)/2.+BET*V);
  T=T1+sqr(M1);
  U=U1+sqr(M1);

  //calculate the pT squared
  double PT2=(T1*U1-S*sqr(M1))/S;

  //if identical particles, divide by two
  DJAC=DJAC/2.;

  jacobian(DJAC*jacobian());
  meMomenta()[2].setMass(mh*GeV);
  meMomenta()[3].setMass(mh*GeV);
  
  if(WHAT<2*mh) {
    jacobian(0.);
    return false;
  }
  Energy q = ZERO;
  Energy pt = sqrt(PT2)*GeV;
  
  //calculate the energy/theta/phi given the sqrt(S) = WHAT and pT
  q = sqrt(sqr(WHAT)/4 - sqr(mh)) * GeV;
  double cth = sqrt( 1 - sqr(pt)/sqr(q) );
  double phi = UseRandom::rnd()*2.0*pi;

  //set up the momenta in the COM 
  meMomenta()[2].setX(pt*sin(phi));
  meMomenta()[2].setY(pt*cos(phi));
  meMomenta()[2].setZ(q*cth);
  
  meMomenta()[3].setX(-pt*sin(phi));
  meMomenta()[3].setY(-pt*cos(phi));
  meMomenta()[3].setZ(-q*cth);

  meMomenta()[2].rescaleEnergy();
  meMomenta()[3].rescaleEnergy();

  return true;
}

void MEHiggsPair::setKinematics() {
  HwMEBase::setKinematics();
}



vector<Complex> MEHiggsPair::iniscal(double AMQ, double S, double T,double U, double M1, double M2) const {
//INITIALIZATION OF SCALAR INTEGRALS
//taken from hpair

  Complex C0AB,C0AC,C0AD,C0BC,C0BD,C0CD,D0ABC,D0BAC,D0ACB;
  double EPM = 1.E-6;
  double c12[1] = {EPM};
  double c3[1] = {S};
  double c456[1] = {AMQ};
  double zerod[1] = {0.};
  double sqrm1[1] = {sqr(M1)};
  double sqrm2[1] = {sqr(M2)};
  double ss[1] = {S};
  double tt[1] = {T};
  double uu[1] = {U};

  Complex DQ2=Complex(sqr(AMQ),0.);

  //get the scalar integrals using the associated hpair fortran functions
  C0AB = c03_(c12, c12, c3, c456, c456, c456)*DQ2;
  C0AC = c03_(c12,sqrm1,tt, c456,c456,c456)*DQ2;
  C0AD = c03_(c12,sqrm2,uu,c456,c456,c456)*DQ2;
  C0BC = c03_(c12,sqrm1,uu,c456,c456,c456)*DQ2;
  C0BD = c03_(c12,sqrm2,tt,c456,c456,c456)*DQ2;
  C0CD = c03_(sqrm1,sqrm2,ss,c456,c456,c456)*DQ2;

  D0ABC = sqr(DQ2)*d04_(zerod,zerod,sqrm1,sqrm2,ss,uu,c456,c456,c456,c456);
  D0BAC = sqr(DQ2)*d04_(zerod,zerod,sqrm1,sqrm2,ss,tt,c456,c456,c456,c456);
  D0ACB = sqr(DQ2)*d04_(zerod,sqrm1,zerod,sqrm2,tt,uu,c456,c456,c456,c456);
  

  //push them back into a vector
  vector<Complex> FormRes;
  FormRes.push_back(C0AB);
  FormRes.push_back(C0AC);
  FormRes.push_back(C0AD);
  FormRes.push_back(C0BC);
  FormRes.push_back(C0BD);
  FormRes.push_back(C0CD);
  FormRes.push_back(D0ABC);
  FormRes.push_back(D0BAC);
  FormRes.push_back(D0ACB);
  
  return FormRes;
}


//calculate the LO matrix element squared as given by hpair, constants are taken care of by me2()
//both top and bottom loops are considered.
double MEHiggsPair::MATRIX(double S, double T,double U, double M1, double M2) const {
  //C--LO MATRIX ELEMENT FOR GG -> HH
  using Constants::pi;

  Complex A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2;
  Complex F1,F2,PROH, PROHEAVY;
  
  //bottom and top masses
  double AMT = _topmass/GeV;
  double AMB = _bottommass/GeV;
  
  //W mass
  double MW = _Wmass/GeV;

  //higgs width and mass
  double FACH=1; 
  double GAMH= _wh/GeV; 
  double AMH = _mh/GeV;

  //vacuum expectation value, self-coupling, top and bottom Yukawas
  double V =  1./sqrt(sqrt(2) * SM().fermiConstant()*GeV*GeV);
  double GHHH =  _selfcoupling * 3.*sqr(AMH) / V;
  double GHT=_yh*AMT/V;
  double GHB=_ybh*AMB/V;

  double ALPS=SM().alphaS(scale());
  if(_fixedalphaS==1) { ALPS = _alphasfixedvalue; }
  double eftscale = _EFTScale/GeV;

  //Higgs propagator
  PROH = Complex(S-sqr(AMH),AMH*GAMH*FACH);

  //EFT-modified self-coupling
  double GHHH_EFT = 0;
  if(_process == 4) { 
    GHHH_EFT = (3.*sqr(AMH) / V) * ( - (3./2.) * _cH + _c6 ) ;  
  }
  
  //EFT-modified top/bottom Yukawas
  double GHT_EFT = 0, GHB_EFT = 0;
  if(_process == 4) { 
    GHT_EFT = AMT/V * ( _ct1 - 0.5 * _cH );
    GHB_EFT = AMB/V * ( _cb1 - 0.5 * _cH );
  }
  //Heavy Higgs
  double AMHEAVY = _heavyHmass/GeV;
  double GAMHEAVY = _heavyHwidth/GeV;
  double GHTHEAVY=_yH*AMT/V;
  double GHBHEAVY=_ybH*AMB/V;


  double FACHEAVY = 1.;
  PROHEAVY = Complex(S-sqr(AMHEAVY),AMHEAVY*GAMHEAVY*FACHEAVY);
  double GHhh = _hhHcoupling * 3.*sqr(AMHEAVY) / V;  // Hhh coupling (for _process == 3)
  
  //calculation of scalar integrals C_ij and D_ijk
  vector<Complex> ScalFacs;
  ScalFacs = iniscal(AMT, S, T, U, M1, M2);
  double amq[1] = {AMT};
  double m1[1] = {M1};
  double m2[1] = {M2};
  double ss[1] = {S};
  double tt[1] = {T};
  double uu[1] = {U};

  Complex C0AB[1] = {ScalFacs[0]};
  Complex C0AC[1] = {ScalFacs[1]};
  Complex C0AD[1] = {ScalFacs[2]};
  Complex C0BC[1] = {ScalFacs[3]};
  Complex C0BD[1] = {ScalFacs[4]};
  Complex C0CD[1] = {ScalFacs[5]};
  Complex D0ABC[1] = {ScalFacs[6]};
  Complex D0BAC[1] = {ScalFacs[7]};
  Complex D0ACB[1] = {ScalFacs[8]};
    //form factors for top quark contributions
  formfac_(amq, ss, tt, uu, m1, m2, C0AB, C0AC, C0AD, C0BC, C0BD, C0CD, D0ABC, D0BAC, D0ACB);

  // this factor needs to be divided out to go from SM EFT to DIM-6 EFT
  double CH = ALPS / (8. * pi);
  
  F1 = Complex(0.,0.);
  F2 = Complex(0.,0.);
  //triangle
  if(_process==0||_process==1||_process==3||_process==4) { 
    F1 = F1 + (form_.H1*(GHT*GHHH/PROH)); //(1)
  }
  //triangle: EFT
  if(_process==4) { 
    F1 = F1 + (form_.H1*(GHT*GHHH_EFT+GHT_EFT*GHHH)/PROH);
  }


  //box
  if(_process==0||_process==2||_process==3||_process==4) { 
    F1 = F1 + (form_.HH1*GHT*GHT); //(2)
    F2 = F2 + (form_.HH2*GHT*GHT);
  }
  //box: EFT
  if(_process==4) { 
    F1 = F1 + 2. * (form_.HH1*GHT_EFT*GHT);
    F2 = F2 + 2. * (form_.HH2*GHT_EFT*GHT);
  }

  //additional EFT from gg -> h -> hh diagram:
  if(_process==4) { 
    F1 = F1 + _cg1 * (2.*S/V)*2*(GHHH/PROH); 
    /* divide (1) by form_.H1*AMT/(2*S) to remove F_triangle
     * multiply by 2
     * and multiply by cg1
     */
  }
  //additional EFT from gg -> hh diagram:
  if(_process==4) { 
    F1 = F1 + _cg2 * (2.*S/sqr(V))*2.;
     /* divide (2) by form_.HH1*sqr(AMT)/(2*S) to remove F_box
      * multiply by 2 
      * and multiply by cg2
      */
  }
  //EFT: new t-tbar-h-h diagram
  if(_process==4) { 
    double GttHH = (3.*_ct2 - _cH) *AMT/(2.*sqr(V)); // we have already replaced yf/sqrt(2) = mf/V
    F1 = F1 + (form_.H1*2.*GttHH); //factor of 2 from tthh feynman rule
  }
  //TESTING BELOW

  /* cout << "form_.H1 = " << form_.H1 << " form_.HH1 = " << form_.HH1 << endl; 
  cout << "Ftri = " << real(form_.H1)*AMT/(2*S) << endl; 
  cout << "Fbox = " << real(form_.HH1)*sqr(AMT)/(2*S) << endl; */

  //form factors for bottom quark contributions
  ScalFacs = iniscal(AMB, S, T, U, M1, M2);
  amq[0] = AMB;
  C0AB[0] = ScalFacs[0];
  C0AC[0] = ScalFacs[1];
  C0AD[0] = ScalFacs[2];
  C0BC[0] = ScalFacs[3];
  C0BD[0] = ScalFacs[4];
  C0CD[0] = ScalFacs[5];
  D0ABC[0] = ScalFacs[6];
  D0BAC[0] = ScalFacs[7];
  D0ACB[0] = ScalFacs[8];
  formfac_(amq, ss, tt, uu, m1, m2, C0AB, C0AC, C0AD, C0BC, C0BD, C0CD, D0ABC, D0BAC, D0ACB);
  //triangle 
  if(_process==0||_process==1||_process==3||_process==4) { 
    F1 = F1 + (form_.H1*(GHB*GHHH/PROH));
  }
  //triangle: EFT
  if(_process==4) { 
    F1 = F1 + (form_.H1*(GHB*GHHH_EFT+GHB_EFT*GHHH)/PROH);
  }

  //box
  if(_process==0||_process==2||_process==3||_process==4) { 
    F1 = F1 + (form_.HH1*GHB*GHB);
    F2 = F2 + (form_.HH2*GHB*GHB);
  }
  //box: EFT
  if(_process==4) { 
    F1 = F1 + 2. * (form_.HH1*GHB_EFT*GHB);
    F2 = F2 + 2. * (form_.HH2*GHB_EFT*GHB);
  }

  //EFT: new b-bbar-h-h diagram
  if(_process==4) { 
    double GbbHH = (3.*_cb2 - _cH) *AMB/(2.*sqr(V));// we have already replaced yf/sqrt(2) = mf/V
    F1 = F1 + (form_.H1*2.*GbbHH);
  }
  
  if(_process==3) { 
    m1[0] = AMH;
    m2[0] = AMHEAVY;
    ScalFacs = iniscal(AMT, S, T, U, M1, M2);
    amq[0] = AMT;
    C0AB[0] = ScalFacs[0];
    C0AC[0] = ScalFacs[1];
    C0AD[0] = ScalFacs[2];
    C0BC[0] = ScalFacs[3];
    C0BD[0] = ScalFacs[4];
    C0CD[0] = ScalFacs[5];
    D0ABC[0] = ScalFacs[6];
    D0BAC[0] = ScalFacs[7];
    D0ACB[0] = ScalFacs[8];
    formfac_(amq, ss, tt, uu, m1, m2, C0AB, C0AC, C0AD, C0BC, C0BD, C0CD, D0ABC, D0BAC, D0ACB);
    F1 = F1 + (form_.H1*(GHTHEAVY*GHhh/PROHEAVY));
    ScalFacs = iniscal(AMB, S, T, U, M1, M2);
    amq[0] = AMB;
    C0AB[0] = ScalFacs[0];
    C0AC[0] = ScalFacs[1];
    C0AD[0] = ScalFacs[2];
    C0BC[0] = ScalFacs[3];
    C0BD[0] = ScalFacs[4];
    C0CD[0] = ScalFacs[5];
    D0ABC[0] = ScalFacs[6];
    D0BAC[0] = ScalFacs[7];
    D0ACB[0] = ScalFacs[8];
    formfac_(amq, ss, tt, uu, m1, m2, C0AB, C0AC, C0AD, C0BC, C0BD, C0CD, D0ABC, D0BAC, D0ACB);
    F1 = F1 + (form_.H1*(GHBHEAVY*GHhh/PROHEAVY));
  }
  
  //square
  double DMAT = 2. * (norm(F1) + norm(F2));  
  return DMAT;
}


 
