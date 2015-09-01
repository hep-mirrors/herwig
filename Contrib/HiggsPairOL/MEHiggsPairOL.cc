// -*- C++ -*-
//
// MEHiggsPairOL.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2009-2011 The Herwig Collaboration
//
// Herwig is licenaaaced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEHiggsPairOL class.
//

#include "MEHiggsPairOL.h"
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


MEHiggsPairOL::MEHiggsPairOL()
  : _selfcoupling(1.0), _process(0), _mh(), _wh(), _fixedalphaS(0), _alphasfixedvalue(0.123601), _alphascale(100.*GeV), _fixedscale(0),
    _implementation(0),
    _includeWidths(0), _includeBquark(1), _basescale(125.*GeV), _scalemultiplier(1.0),
    Mass_E(0.0),
    Mass_M(0.0),
    Mass_L(0.0),
    Mass_T(100.),
    Mass_U(0.0),
    Mass_C(0.0),
    Mass_D(0.0),
    Mass_S(0.0),
    Mass_B(5.),
    Mass_Z(91.),
    Mass_W(80.4),
    Mass_H(125.),
    Width_C(0.0),
    Width_B(0.0),
    Width_T(0.0),
    Width_W(0.),
    Width_Z(0.),
    Width_H(0.0),
    Coupl_Alpha_QED(1./132.348905),
    Coupl_Alpha_QCD(0.12360931820752918),
    _id(0), 
    _idtri(0), 
    _idbox(0),
    _idint(0)
{
  massOption(vector<unsigned int>(2,0));
}

int MEHiggsPairOL::nDim() const {
  return 5;
}

void MEHiggsPairOL::doinitrun() {
  //  if(_implementation == 0 || _implementation == 2) { InitOpenLoops(); }
}

void MEHiggsPairOL::doinit() {
  _theModel = generator()->standardModel();
  tcHiggsPairPtr hwHiggsPair=dynamic_ptr_cast<tcHiggsPairPtr>(_theModel);
  if(hwHiggsPair){
    _selfcoupling=hwHiggsPair->selfcoupling();
    _process=hwHiggsPair->process();
    _includeBquark=hwHiggsPair->includeBquark();
    _includeWidths=hwHiggsPair->includeWidths();
    _implementation=hwHiggsPair->implementation();
    _fixedscale=hwHiggsPair->fixedscale();
    _alphascale=hwHiggsPair->alphascale();
    _fixedalphaS=hwHiggsPair->fixedalphaS();
    _alphasfixedvalue = hwHiggsPair->alphasfixedvalue();
    _basescale=hwHiggsPair->basescale();
    _scalemultiplier=hwHiggsPair->scalemultiplier();
    if(_implementation == 0 || _implementation == 2) { 
      _id =  hwHiggsPair->id();
      _idtri = hwHiggsPair->idtri();
      _idbox = hwHiggsPair->idbox();
      _idint = hwHiggsPair->idint();
    }
  }

  if(_implementation == 0 && _process != 0) { cerr << "The OpenLoops implementation currently only contains SM production" << endl; exit(1); }

  if(_process == 4 || _process == 5) { cerr << "HH and hH production not implemented yet, please choose an hh production subprocess" << endl; exit(1); } 

  higgs(getParticleData(ParticleID::h0));
  PDPtr top = getParticleData(ParticleID::t);
  _topmass = top->mass();
  //  cout << "_topmass = " << _topmass << endl;
  PDPtr zboson = getParticleData(ParticleID::Z0);
  PDPtr wboson = getParticleData(24);

  _zmass = zboson->mass();
  PDPtr bottom = getParticleData(ParticleID::b);
  _bottommass = bottom->mass();
  PDPtr heavyH = getParticleData(35);
  _heavyHmass = heavyH->mass();
  _heavyHwidth = heavyH->width();
  
  Mass_E = 0.0;
  Mass_M = 0.0;
  Mass_L = 0.0;

  Mass_T = top->mass()/GeV;
  Mass_U = 0.0;
  Mass_C = 0.0;
  Mass_D = 0.0;
  Mass_S=  0.0;
  Mass_B = bottom->mass()/GeV;
  if(_includeBquark == 0) { Mass_B = 0; cout << "Not including B quark mass." << endl; }

  Mass_Z = zboson->mass()/GeV;
  Mass_W = wboson->mass()/GeV;
  Mass_H = _higgs->mass()/GeV;

  Width_C = 0.0;
  Width_B = 0.0;
  Width_T = 0.0;
  Width_W = 0.0;
  Width_Z = 0.0;
  Width_H = 0.0;
  if(_includeWidths == 1 && _implementation == 0) {
    Width_T = top->width()/GeV;
    Width_W = wboson->width()/GeV;
    Width_Z = zboson->width()/GeV;
    Width_H = _higgs->width()/GeV;
  } 
  if(_includeWidths == 1 && _implementation != 0) {
    cout << "Warning: top quark width effects are only included in the OpenLoops implementation" << endl;
  }

  Coupl_Alpha_QED= SM().alphaEM(sqr(2*Mass_H*GeV)); //change to correct scale.
  Coupl_Alpha_QCD = 0.12360931820752918;//SM().alphaS(sqr(2*Mass_H*GeV)); //change to correct scale
  // cout << "init init init" << endl;
  //  if(_implementation == 0 || _implementation == 2) { InitOpenLoops(); }
  //  cout << "init init init 2" << endl;

  _mh = _higgs->mass();
  _m1 = _higgs->mass();
  _m2 = _higgs->mass();
  _wh = _higgs->width();

  if(_higgs->massGenerator()) {
    _hmass=dynamic_ptr_cast<GenericMassGeneratorPtr>(_higgs->massGenerator());
  }

  HwMEBase::doinit();
}

void MEHiggsPairOL::rebind(const TranslationMap & trans) {
  HwMEBase::rebind(trans);
}

IVector MEHiggsPairOL::getReferences() {
  IVector ret = HwMEBase::getReferences();
  return ret;
}

IBPtr MEHiggsPairOL::clone() const {
    return new_ptr(*this);
  }
  
IBPtr MEHiggsPairOL::fullclone() const {
    return new_ptr(*this);
  }

void MEHiggsPairOL::persistentOutput(PersistentOStream & os) const {
  os <<_theModel << _selfcoupling  << _process << _higgs << ounit(_topmass,GeV) << ounit(_bottommass,GeV) << ounit(_zmass,GeV) << ounit(_m1,GeV) << ounit(_m2,GeV) << ounit(_mh,GeV) << ounit(_wh,GeV) << _hmass << Mass_E << Mass_M << Mass_L << Mass_T << Mass_U << Mass_C << Mass_D << Mass_S << Mass_B << Mass_Z << Mass_W << Mass_H << Width_C << Width_B << Width_T << Width_W << Width_Z << Width_H << Coupl_Alpha_QED << Coupl_Alpha_QCD << _implementation << _fixedalphaS << _alphasfixedvalue << _includeWidths << ounit(_alphascale,GeV) << ounit(_basescale,GeV) << _fixedscale << _includeBquark << _scalemultiplier << _id << _idtri << _idbox << _idint;
}

void MEHiggsPairOL::persistentInput(PersistentIStream & is, int) {
  is >> _theModel >> _selfcoupling >> _process >> _higgs >> iunit(_topmass,GeV) >> iunit(_bottommass,GeV) >> iunit(_zmass,GeV) >> iunit(_m1,GeV) >> iunit(_m2,GeV) >> iunit(_mh,GeV) >> iunit(_wh,GeV) >> _hmass >> Mass_E >> Mass_M >> Mass_L >> Mass_T >> Mass_U >> Mass_C >> Mass_D >> Mass_S >> Mass_B >> Mass_Z >> Mass_W >> Mass_H >> Width_C >> Width_B >> Width_T >> Width_W >> Width_Z >> Width_H >> Coupl_Alpha_QED >> Coupl_Alpha_QCD >> _implementation >> _fixedalphaS >> _alphasfixedvalue >> _includeWidths >> iunit(_alphascale,GeV) >> iunit(_basescale,GeV) >> _fixedscale >> _includeBquark >> _scalemultiplier >> _id >> _idtri >> _idbox >> _idint;

}

Energy2 MEHiggsPairOL::scale() const {
  if(_fixedscale==0 || _fixedscale==5 ) { return sqr(_scalemultiplier)*sHat(); } 
  else if (_fixedscale==1 || _fixedscale==2 || _fixedscale==3) { return (sqr(_scalemultiplier)*sqr(_basescale)); }
}

// Definition of the static class description member.
ClassDescription<MEHiggsPairOL> MEHiggsPairOL::initMEHiggsPairOL;

void MEHiggsPairOL::Init() {
  static ClassDocumentation<MEHiggsPairOL> documentation
    ("The MEHiggsPairOL class implements the hh 2->2 processes in hadron-hadron"
     " collisions");
}


CrossSection MEHiggsPairOL::dSigHatDR() const {
  using Constants::pi;
  return me2()/(16.0*pi*sqr(sHat())/GeV/GeV)*sqr(hbarc);
}

Selector<MEBase::DiagramIndex>
MEHiggsPairOL::diagrams(const DiagramVector & diags) const {
  // select the diagram, this is easy for us as we have already done it
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    sel.insert(1.0, i);
  }
  return sel;
}

//diagrams with box or triangle
void MEHiggsPairOL::getDiagrams() const {
  // get the particle data objects
  PDPtr gluon = getParticleData(ParticleID::g);
  PDPtr boxon = getParticleData(99926);
  PDPtr triangon = getParticleData(99927);
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
MEHiggsPairOL::colourGeometries(tcDiagPtr diag) const {
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
double MEHiggsPairOL::me2() const {
  using Constants::pi;
  double Mc = _mh/GeV;
  double Md = _mh/GeV;
  double s = sHat()/GeV/GeV;
  //XLAM = beta*S
  double XLAM=sqrt(sqr(s+sqr(Mc)-sqr(Md))-4.*s*sqr(Mc));
  double ALPS(0.);
  if(_fixedalphaS == 0) { ALPS=SM().alphaS(scale()); }
  else if(_fixedalphaS == 1) { ALPS = SM().alphaS(sqr(_alphascale)); }
  else if(_fixedalphaS == 2) { ALPS = _alphasfixedvalue; }
  double mesq(0.);

  /* 
   * HPAIR starts here
   */ 
  if(_implementation == 1 || _implementation == 2) { 
    //get the masses and Mandestam variables
    Energy mh(_mh);

    double t = tHat()/GeV/GeV;
    double u = uHat()/GeV/GeV;
    
    //calculate the matrix element from the form factors
    mesq = MATRIX(s, t, u, Mc, Md);
    //QCD factors
    double QCD=sqr(ALPS);
  
    //the jacobian is jacobian()*XLAM
    mesq *= jacobian()*XLAM; 
    mesq *= QCD / (2048.*sqr(pi));
    if(_implementation == 1) { return mesq; }
  }
  
  /* 
   * OpenLoops starts here
   */
  double m2ol(0.);
  if(_implementation == 0 || _implementation == 2) { 
    double m2_tree, m2_loop[3], acc;
    double pp[5*ol_n_external(_id)];
    ol_setparameter_double("alpha_s", ALPS);
    double fermiConst = SM().fermiConstant()*GeV*GeV;
    double Coupl_Alpha_QED_ = 1/((pi / ( sqrt(2) *  fermiConst * sqr(Mass_W) )) * 1/( 1 - sqr(Mass_W)/sqr(Mass_Z) ));
    ol_setparameter_double("alpha",Coupl_Alpha_QED_);
    // Set parameter: renormalization scale
    ol_setparameter_double("mu", sqrt(scale()/GeV/GeV));

    //fill in particle momenta
    for(int kk = 0; kk < 4; kk++) { 
      pp[5*kk] = meMomenta()[kk].e()/GeV;
      pp[5*kk+1] = meMomenta()[kk].x()/GeV;
      pp[5*kk+2] = meMomenta()[kk].y()/GeV;
      pp[5*kk+3] = meMomenta()[kk].z()/GeV;
      pp[5*kk+4] = -1;
    }
    /*  for (int k = 0; k < 4; k++) {
       std::cout << "P[" << k+1 << "] = " << pp[5*k] << "  " << pp[5*k+1]
       << "  " << pp[5*k+2] << "  " << pp[5*k+3] << std::endl;
       }*/
    if(_process == 0) { ol_evaluate_loop2(_id, pp, &m2_tree, &acc); }
    if(_process == 1) { ol_evaluate_loop2(_idtri, pp, &m2_tree, &acc); }
    if(_process == 2) { ol_evaluate_loop2(_idbox, pp, &m2_tree, &acc); }
    m2ol = m2_tree;
    if(_implementation == 0) { /*cout << "  m2_tree*((XLAM)) = " <<  m2_tree*((XLAM)) << endl;*/ return m2ol*((XLAM)); }
  }
  cout << setprecision(15) << "r = " << mesq/(m2ol*(XLAM)) << " implementation = " <<  _implementation  << endl;	
  return mesq;
}

bool MEHiggsPairOL::generateKinematics(const double * r) { 
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
  double cth = sqrt(1 - sqr(pt)/sqr(q));
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

void MEHiggsPairOL::setKinematics() {
  HwMEBase::setKinematics();
}


vector<Complex> MEHiggsPairOL::iniscal(double AMQ, double S, double T,double U, double M1, double M2) const {
//INITIALIZATION OF SCALAR INTEGRALS
//taken from hpair

  Complex C0AB,C0AC,C0AD,C0BC,C0BD,C0CD,D0ABC,D0BAC,D0ACB;
  double EPM = 1.E-8;
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
double MEHiggsPairOL::MATRIX(double S, double T,double U, double M1, double M2) const {
  //C--LO MATRIX ELEMENT FOR GG -> HH
   using Constants::pi;

  Complex A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2;
  Complex F1,F2,PROH, PROHEAVY;
  
  //bottom and top masses
  double AMT = _topmass/GeV;
  double AMB = _bottommass/GeV;

  //higs width and mass
  double FACH=1.; 
  //  double GAMH= _wh/GeV; 
  double GAMH = 0.;
  double AMH = _mh/GeV;

  //vacuum expectation value, self-coupling, top and bottom Yukawas
  double V =  1./sqrt(sqrt(2.) * SM().fermiConstant()*GeV*GeV);
  double GHHH = _selfcoupling * 3. * sqr(AMH) / V;
  double GHT=AMT/V;
  double GHB=AMB/V;
  
  //Higgs propagator
  PROH = Complex(S-sqr(AMH),AMH*GAMH*FACH);
  
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
    F1 = F1 + form_.H1*(GHT*GHHH/PROH);
  }
  //box
  if(_process==0||_process==2||_process==3) { 
    F1 = F1 + form_.HH1*GHT*GHT;
    F2 = F2 + form_.HH2*GHT*GHT;    
  }

  //form factors for Bottom quark contributions
  if(_includeBquark == 1) {
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
    if(_process==0||_process==1) { 
      F1 = F1 + (form_.H1*GHB*GHHH/PROH);//comment out for no b-quark
    }
    //box
    if(_process==0||_process==2||_process==3) { 
    F1 = F1 + form_.HH1*GHB*GHB; //comment out for no b
    F2 = F2 + form_.HH2*GHB*GHB; //comment out for no b
    }
  }
  
  //square
  double DMAT = 2. * (norm(F1) + norm(F2));  
  return DMAT;
}


 
void MEHiggsPairOL::InitOpenLoops() { 
  //ol_setparameter_int("order_ew", 2);
  /*ol_setparameter_int("redlib1",1);
  ol_setparameter_int("redlib2",7);
  ol_setparameter_int("verbose",2);

  ol_setparameter_int("stability_mode",21);
  ol_setparameter_double("mass(23)", Mass_Z);
  ol_setparameter_double("mass(24)", Mass_W);
  ol_setparameter_double("mass(25)", Mass_H);
  ol_setparameter_double("mass(1)", Mass_D);
  ol_setparameter_double("mass(2)", Mass_U);
  ol_setparameter_double("mass(3)", Mass_S);
  ol_setparameter_double("mass(4)", Mass_C);
  ol_setparameter_double("mass(5)", Mass_B);
  ol_setparameter_double("mass(6)", Mass_T);
  
  ol_setparameter_double("width(23)", Width_Z);
  ol_setparameter_double("width(24)", Width_W);
  ol_setparameter_double("width(25)", Width_H);
  
  _id = ol_register_process("21 21 -> 25 25", 12);
  char only3h[] = "only3h";
  ol_setparameter_string("approx", only3h);
  _idtri = ol_register_process("21 21 -> 25 25", 12);
  char no3h[] = "no3h";
  ol_setparameter_string("approx", no3h);
  _idbox = ol_register_process("21 21 -> 25 25", 12);
  char interf3h[] = "interf3h";
  ol_setparameter_string("approx", interf3h);
  _idint = ol_register_process("21 21 -> 25 25", 12);

  // Initialize OpenLoops
  ol_start();*/
}
