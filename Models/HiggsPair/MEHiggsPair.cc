// -*- C++ -*-
//
// MEHiggsPair.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2009-2011 The Herwig Collaboration
//
// Herwig++ is licenaaaced under version 2 of the GPL, see COPYING for details.
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
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Cuts/Cuts.h"
#include "Herwig++/Utilities/Maths.h"
#include <fstream>
#include <cmath>

namespace Herwig {


}

using namespace Herwig;


MEHiggsPair::MEHiggsPair()
  : _selfcoupling(1.0),_hhHcoupling(1.0), _process(0), _mh(), _wh() {
  massOption(vector<unsigned int>(2,0));
}

int MEHiggsPair::nDim() const {
  return 5;
}


void MEHiggsPair::doinit() {

  if(_process == 4 || _process == 5) { cerr << "HH and hH production not implemented yet, please choose an hh production subprocess" << endl; exit(1); } 

  higgs(getParticleData(ParticleID::h0));
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
  os << _selfcoupling << _hhHcoupling << _process << _higgs << ounit(_topmass,GeV) << ounit(_bottommass,GeV) << ounit(_zmass,GeV) << ounit(_m1,GeV) << ounit(_m2,GeV) << ounit(_mh,GeV) << ounit(_wh,GeV) << _hmass;
}

void MEHiggsPair::persistentInput(PersistentIStream & is, int) {
  is >> _selfcoupling >> _hhHcoupling >> _process >> _higgs >> iunit(_topmass,GeV) >> iunit(_bottommass,GeV) >> iunit(_zmass,GeV) >> iunit(_m1,GeV) >> iunit(_m2,GeV) >> iunit(_mh,GeV) >> iunit(_wh,GeV) >> _hmass;

}

Energy2 MEHiggsPair::scale() const {
   return sHat();
}

// Definition of the static class description member.
ClassDescription<MEHiggsPair> MEHiggsPair::initMEHiggsPair;

void MEHiggsPair::Init() {

  static ClassDocumentation<MEHiggsPair> documentation
    ("The MEHiggsPair class implements the transplanckian 2->2 processes in hadron-hadron"
     " collisions");

 static Parameter<MEHiggsPair, double> interfaceSelfCoupling
    ("SelfCoupling",
     "Multiplier for the SM Higgs triple coupling",
     &MEHiggsPair::_selfcoupling, 1.0, -10.0, 10.0,
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
  static SwitchOption interfaceProcessHH
    (interfaceProcess,
     "ggToHH",
     "Include all gg->HH subprocesses",
     4);
  static SwitchOption interfaceProcessHh
    (interfaceProcess,
     "ggToHh",
     "Include only gg->hH subprocesses",
     5);
}

CrossSection MEHiggsPair::dSigHatDR() const {
  using Constants::pi;
  return me2()/(16.0*sqr(pi)*sHat())*sqr(hbarc);
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
  if(_process==0||_process==1||_process==3) {  
    add(new_ptr((Tree2toNDiagram(2),gluon,gluon,
		 1,triangon,3,higgs(),higgs(),-1)));    
  }
  if(_process==0||_process==2||_process==3) {
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
  double mesq = 0.;

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
  double QCD=sqr(ALPS);
  //XLAM = beta*S
  double XLAM=sqrt(sqr(s+sqr(Mc)-sqr(Md))-4.*s*sqr(Mc));
  //the jacobian is jacobian()*XLAM
  mesq *= jacobian()*XLAM; 
  mesq *= QCD * SM().fermiConstant()*GeV*GeV / (128 * sqrt(2) * sqr(pi)  );

  return mesq/s;    
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
 
  Complex A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2;
  Complex F1,F2,PROH, PROHEAVY;
  
  //bottom and top masses
  double AMT = _topmass/GeV;
  double AMB = _bottommass/GeV;

  //higs width and mass
  double FACH=1; 
  double GAMH= _wh/GeV; 
  double AMH = _mh/GeV;

  //vacuum expectation value, self-coupling, top and bottom Yukawas
  double V =  1./sqrt(sqrt(2) * SM().fermiConstant()*GeV*GeV);
  double GHHH =  _selfcoupling * 3.*sqr(AMH) / V;
  double GHT=AMT/V;
  double GHB=AMB/V;



  //Higgs propagator
  PROH = Complex(S-sqr(AMH),AMH*GAMH*FACH);

  //Heavy Higgs
  double AMHEAVY = _heavyHmass/GeV;
  double GAMHEAVY = _heavyHwidth/GeV;			      
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
  
  F1 = Complex(0.,0.);
  F2 = Complex(0.,0.);
  //triangle
  if(_process==0||_process==1||_process==3) { 
    F1 = F1 + AMT*(form_.H1*(GHT*GHHH/PROH));
  }
  //box
  if(_process==0||_process==2||_process==3) { 
    F1 = F1 + AMT*(form_.HH1*GHT*GHT);
    F2 = F2 + AMT*(form_.HH2*GHT*GHT);
  }

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
  if(_process==0||_process==1||_process==3) { 
    F1 = F1 + AMB*(form_.H1*(GHB*GHHH/PROH));
  }
  //box
  if(_process==0||_process==2||_process==3) { 
    F1 = F1 + AMB*(form_.HH1*GHB*GHB);
    F2 = F2 + AMB*(form_.HH2*GHB*GHB);
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
    F1 = F1 + AMT*(form_.H1*(GHT*GHhh/PROHEAVY));
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
    F1 = F1 + AMB*(form_.H1*(GHT*GHhh/PROHEAVY));
  }
  
  //square
  double DMAT = 2. * (norm(F1) + norm(F2));  
  return DMAT;
}


 
