// -*- C++ -*-
//
// MEHiggsPairJet.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2009-2011 The Herwig Collaboration
//
// Herwig is licenaaaced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEHiggsPairJet class.
//

#include "MEHiggsPairJet.h"
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

using namespace Herwig;


MEHiggsPairJet::MEHiggsPairJet()
  :_selfcoupling(1.0),_hhHcoupling(1.0), _process(0), _mh(), _wh(), _fixedalphaS(0), _alphasfixedvalue(0.123601), _alphascale(100.*GeV), _fixedscale(0),
    _implementation(0),
   _includeWidths(0), _includeBquark(1), _maxflavour(5), _basescale(125.*GeV), _subprocessjet(0), _alphasreweight(0), _scalemultiplier(1.0),
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
    Coupl_Alpha_QCD(0.1)
{
  massOption(vector<unsigned int>(2,0));
}

int MEHiggsPairJet::nDim() const {
  return 5;
}

void MEHiggsPairJet::doinitrun() {
  //  InitOpenLoops();
}

void MEHiggsPairJet::doinit() {
  _theModel = generator()->standardModel();
  tcHiggsPairPtr hwHiggsPair=dynamic_ptr_cast<tcHiggsPairPtr>(_theModel);
  if(hwHiggsPair){
    _selfcoupling=hwHiggsPair->selfcoupling();
    _process=hwHiggsPair->process();
    _includeBquark=hwHiggsPair->includeBquark();
    _includeWidths=hwHiggsPair->includeWidths();
    _implementation=hwHiggsPair->implementation();
    _maxflavour=hwHiggsPair->maxflavour();
    _fixedscale=hwHiggsPair->fixedscale();
    _alphascale=hwHiggsPair->alphascale();
    _fixedalphaS=hwHiggsPair->fixedalphaS();
    _basescale=hwHiggsPair->basescale();
    _subprocessjet=hwHiggsPair->subprocessjet();
    _alphasreweight=hwHiggsPair->alphasreweight();
    _scalemultiplier=hwHiggsPair->scalemultiplier();
    
    _id1 =  hwHiggsPair->id1();
    _id1tri = hwHiggsPair->id1tri();
    _id1box = hwHiggsPair->id1box();
    _id1int = hwHiggsPair->id1int();

    _id2 =  hwHiggsPair->id2();
    _id2tri = hwHiggsPair->id2tri();
    _id2box = hwHiggsPair->id2box();
    _id2int = hwHiggsPair->id2int();

    _id3 =  hwHiggsPair->id3();
    _id3tri = hwHiggsPair->id3tri();
    _id3box = hwHiggsPair->id3box();
    _id3int = hwHiggsPair->id3int();

    _id4 =  hwHiggsPair->id4();
    _id4tri = hwHiggsPair->id4tri();
    _id4box = hwHiggsPair->id4box();
    _id4int = hwHiggsPair->id4int();

  }
 
  //  cout << _selfcoupling << " " << _process << " " << _includeBquark << " " << _includeWidths << " " << _implementation << " " << _maxflavour << " " << _fixedscale << " " << _alphascale << " " << _fixedalphaS << endl;

  if(_implementation == 0 && _process != 0) { cerr << "The OpenLoops implementation currently only contains SM production" << endl; exit(1); }

  if(_process == 4 || _process == 5) { cerr << "HH and hH production not implemented yet, please choose an hh production subprocess" << endl; exit(1); } 

  higgs(getParticleData(ParticleID::h0));
  PDPtr top = getParticleData(ParticleID::t);
  _topmass = top->mass();
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
  if(_includeBquark == 0) { Mass_B = 0; }

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
  Coupl_Alpha_QCD = SM().alphaS(sqr(2*Mass_H*GeV)); //change to correct scale

  //  InitOpenLoops();

  _mh = _higgs->mass();
  _m1 = _higgs->mass();
  _m2 = _higgs->mass();
  _wh = _higgs->width();

  if(_higgs->massGenerator()) {
    _hmass=dynamic_ptr_cast<GenericMassGeneratorPtr>(_higgs->massGenerator());
  }
  HwMEBase::doinit();
}

void MEHiggsPairJet::rebind(const TranslationMap & trans) {
  HwMEBase::rebind(trans);
}

IVector MEHiggsPairJet::getReferences() {
  IVector ret = HwMEBase::getReferences();
  return ret;
}

IBPtr MEHiggsPairJet::clone() const {
  return new_ptr(*this);
}
  
IBPtr MEHiggsPairJet::fullclone() const {
    return new_ptr(*this);
  }

void MEHiggsPairJet::persistentOutput(PersistentOStream & os) const {
  os <<_theModel << _selfcoupling << _hhHcoupling << _process << _higgs << ounit(_topmass,GeV) << ounit(_bottommass,GeV) << ounit(_zmass,GeV) << ounit(_m1,GeV) << ounit(_m2,GeV) << ounit(_mh,GeV) << ounit(_wh,GeV) << _hmass << Mass_E << Mass_M << Mass_L << Mass_T << Mass_U << Mass_C << Mass_D << Mass_S << Mass_B << Mass_Z << Mass_W << Mass_H << Width_C << Width_B << Width_T << Width_W << Width_Z << Width_H << Coupl_Alpha_QED << Coupl_Alpha_QCD << _implementation << _fixedalphaS << _alphasfixedvalue << _includeWidths << ounit(_alphascale,GeV) << ounit(_basescale,GeV) << _fixedscale << _includeBquark << _scalemultiplier << _id1 << _id1tri << _id1box << _id1int << _id2 << _id2tri << _id2box << _id2int << _id3 << _id3tri << _id3box << _id3int << _id4 << _id4tri << _id4box << _id4int;

}

void MEHiggsPairJet::persistentInput(PersistentIStream & is, int) {
  is >> _theModel >> _selfcoupling >> _hhHcoupling >> _process >> _higgs >> iunit(_topmass,GeV) >> iunit(_bottommass,GeV) >> iunit(_zmass,GeV) >> iunit(_m1,GeV) >> iunit(_m2,GeV) >> iunit(_mh,GeV) >> iunit(_wh,GeV) >> _hmass >> Mass_E >> Mass_M >> Mass_L >> Mass_T >> Mass_U >> Mass_C >> Mass_D >> Mass_S >> Mass_B >> Mass_Z >> Mass_W >> Mass_H >> Width_C >> Width_B >> Width_T >> Width_W >> Width_Z >> Width_H >> Coupl_Alpha_QED >> Coupl_Alpha_QCD >> _implementation >> _fixedalphaS >> _alphasfixedvalue >> _includeWidths >> iunit(_alphascale,GeV) >> iunit(_basescale,GeV) >> _fixedscale >> _includeBquark >> _scalemultiplier >> _id1 >> _id1tri >> _id1box >> _id1int  >> _id2 >> _id2tri >> _id2box >> _id2int >> _id3 >> _id3tri >> _id3box >> _id3int  >> _id4 >> _id4tri >> _id4box >> _id4int;
}

Energy2 MEHiggsPairJet::scale() const {
  return _scale; 
}  


// Definition of the static class description member.
ClassDescription<MEHiggsPairJet> MEHiggsPairJet::initMEHiggsPairJet;

void MEHiggsPairJet::Init() {

  static ClassDocumentation<MEHiggsPairJet> documentation
    ("The MEHiggsPairJet class implements the Higgs boson pair + jet 2->3 processes in hadron-hadron"
     " collisions");  
}

CrossSection MEHiggsPairJet::dSigHatDR() const {
  using Constants::pi;
  double mesq = me2();
  return mesq*jacobian()/(16.0*sqr(Constants::pi)*sHat())*sqr(hbarc);
}



Selector<MEBase::DiagramIndex>
MEHiggsPairJet::diagrams(const DiagramVector & diags) const {
  // select the diagram, this is easy for us as we have already done it
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    sel.insert(1.0, i);
  }
  return sel;
}

//diagrams with box or triangle
void MEHiggsPairJet::getDiagrams() const {

  // get the particle data objects
  PDPtr gluon = getParticleData(ParticleID::g);

  PDPtr boxon = getParticleData(99926);
  PDPtr triangon = getParticleData(99927);
 
  if(_process==0||_process==1||_process==3) {  
    if(_subprocessjet == 0|| _subprocessjet == 1) { add(new_ptr((Tree2toNDiagram(3), gluon, gluon, gluon, 1, boxon, 2, gluon, 4, higgs(), 4, higgs(), -2))); }
    for (int i = 1; i <= _maxflavour; ++i ) {
      tcPDPtr q = getParticleData(i);
      tcPDPtr qb = q->CC();
      if(_subprocessjet == 0 || _subprocessjet == 3) { add(new_ptr((Tree2toNDiagram(3), q, gluon, gluon, 1, boxon, 2, q,  4, higgs(), 4, higgs(), -4))); }
      if(_subprocessjet == 0 || _subprocessjet == 2) { add(new_ptr((Tree2toNDiagram(3), qb, gluon, gluon, 1, boxon, 2, qb,  4, higgs(), 4, higgs(), -6))); }
      if(_subprocessjet == 0 || _subprocessjet == 4) { add(new_ptr((Tree2toNDiagram(2), q, qb, 1, gluon, 3, boxon, 3, gluon, 4, higgs(), 4, higgs(), -8))); }
    }
  }
  if(_process==0||_process==2||_process==3) {
    if(_subprocessjet == 0|| _subprocessjet == 1) { add(new_ptr((Tree2toNDiagram(3), gluon, gluon, gluon, 1, triangon, 2, gluon,  4, higgs(), 4, higgs(), -1))); }
    for (int i = 1; i <= _maxflavour; ++i ) {
      tcPDPtr q = getParticleData(i);
      tcPDPtr qb = q->CC();
      if(_subprocessjet == 0 || _subprocessjet == 3) { add(new_ptr((Tree2toNDiagram(3), q, gluon, gluon, 1, triangon, 2, q,  4, higgs(), 4, higgs(), -3))); }
      if(_subprocessjet == 0 || _subprocessjet == 2) { add(new_ptr((Tree2toNDiagram(3), qb, gluon, gluon, 1, triangon, 2, qb,  4, higgs(), 4, higgs(), -5))); }
      if(_subprocessjet == 0 || _subprocessjet == 4) { add(new_ptr((Tree2toNDiagram(2), q, qb, 1, gluon, 3, triangon, 3, gluon, 4, higgs(), 4, higgs(), -7))); }
    }
  }
}

//both diagrams possess the same colour structure
Selector<const ColourLines *>
MEHiggsPairJet::colourGeometries(tcDiagPtr diag) const {
  // colour lines for gg to hh

  static const ColourLines cgghhg1("1 2 5, -1 -2 3, -3 -5");  
  static const ColourLines cgghhg2("-1 -2 -5, 1 2 -3, 3 5");  
  static const ColourLines cqghhq("3 2 5, 1 -2 -3");  
  static const ColourLines cqbghhqb("-3 -2 -5, -1 2 3");  
  static const ColourLines cqqbhhg("1 3 5, -2 -3 -5");  
  static const ColourLines cqbqhhg("2 3 5, -1 -3 -5");  


  Selector<const ColourLines *> sel;
  switch(abs(diag->id())) {
  case 1: case 2: 
    sel.insert(1.0, &cgghhg1);
    sel.insert(1.0, &cgghhg2);
    break;
  case 3: case 4: 
    sel.insert(1.0, &cqghhq);
    break;
  case 5: case 6: 
    sel.insert(1.0, &cqbghhqb);
    break;
  case 7: case 8: 
    sel.insert(1.0, &cqqbhhg);
    break;
  case 9: case 10: 
    sel.insert(1.0, &cqbqhhg);
    break;
  } 

  return sel;
}

//the matrix element squared with the appropriate factors
double MEHiggsPairJet::me2() const {
  using Constants::pi;
  double mesq(0.);
  /* 
   * OpenLoops starts here
   */
  //InitOpenLoops();
  // parameters_write_();
  // loop_parameters_write_();

  /* 
   * OpenLoops ME
   */
  //define relevant arrays to pass by reference

  
  /* 
   * OpenLoops starts here
   */
  double m2_tree, m2_loop[3], acc;

  double ALPS(0.);
  if(_fixedalphaS == 0) { ALPS=SM().alphaS(scale()); }
  else if(_fixedalphaS == 1) { ALPS = SM().alphaS(sqr(_alphascale)); }
  else if(_fixedalphaS == 2) { ALPS = _alphasfixedvalue; }

  ol_setparameter_double("alpha_s", ALPS);
  double fermiConst = SM().fermiConstant()*GeV*GeV;
  double Coupl_Alpha_QED_ = 1/((pi / ( sqrt(2) *  fermiConst * sqr(Mass_W) )) * 1/( 1 - sqr(Mass_W)/sqr(Mass_Z) ));
  ol_setparameter_double("alpha",Coupl_Alpha_QED_);
  // Set parameter: renormalization scale
  ol_setparameter_double("mu", sqrt(scale()/GeV/GeV));

  int _idtoeval;
  // q g to q HH
  if(mePartonData()[0]->id()<=6&&mePartonData()[0]->id()>0&&
     mePartonData()[1]->id()==ParticleID::g) {
    _idtoeval = _id1;
    
  }
  // qbar g to qbar HH
  else if(mePartonData()[0]->id()>=-6&&mePartonData()[0]->id()<0&&
	  mePartonData()[1]->id()==ParticleID::g) {
    _idtoeval = _id2;
  }
  // g g to g HH
  else if(mePartonData()[0]->id()==ParticleID::g&&
	  mePartonData()[1]->id()==ParticleID::g) {
    _idtoeval = _id3;

  }
  // q qbar to g HH
  else if(mePartonData()[0]->id()<=6&&
	  mePartonData()[1]->id()<=6) {
    _idtoeval = _id4;	
  }

  double pp[5*ol_n_external(_idtoeval)];

  //  cout << "ALPS = " << ALPS << " mu = " <<  sqrt(scale()/GeV/GeV) << " Coupl_Alpha_QED_ = " << Coupl_Alpha_QED_ << endl;
    
  //fill in particle momenta
  for(int kk = 0; kk < 5; kk++) { 
    pp[5*kk] = meMomenta()[kk].e()/GeV;
    pp[5*kk+1] = meMomenta()[kk].x()/GeV;
    pp[5*kk+2] = meMomenta()[kk].y()/GeV;
    pp[5*kk+3] = meMomenta()[kk].z()/GeV;
    pp[5*kk+4] = -1;
  }

  // q g to q HH
  if(mePartonData()[0]->id()<=6&&mePartonData()[0]->id()>0&&
     mePartonData()[1]->id()==ParticleID::g) {
    if(_process == 0) { ol_evaluate_loop2(_id1, pp, &m2_tree, &acc); }
    if(_process == 1) { ol_evaluate_loop2(_id1tri, pp, &m2_tree, &acc); }
    if(_process == 2) { ol_evaluate_loop2(_id1box, pp, &m2_tree, &acc); }
    //   cout << "qg->qHH" << endl;

  }
  // qbar g to qbar HH
  else if(mePartonData()[0]->id()>=-6&&mePartonData()[0]->id()<0&&
	  mePartonData()[1]->id()==ParticleID::g) {
    if(_process == 0) { ol_evaluate_loop2(_id2, pp, &m2_tree, &acc); }
    if(_process == 1) { ol_evaluate_loop2(_id2tri, pp, &m2_tree, &acc); }
    if(_process == 2) { ol_evaluate_loop2(_id2box, pp, &m2_tree, &acc); }
    // cout << "qbg->qbHH" << endl;

  }
  // g g to g HH
  else if(mePartonData()[0]->id()==ParticleID::g&&
	  mePartonData()[1]->id()==ParticleID::g) {
    if(_process == 0) { ol_evaluate_loop2(_id3, pp, &m2_tree, &acc); }
    if(_process == 1) { ol_evaluate_loop2(_id3tri, pp, &m2_tree, &acc); }
    if(_process == 2) { ol_evaluate_loop2(_id3box, pp, &m2_tree, &acc); }
    // cout << "gg->gHH" << endl;
  }
  // q qbar to g HH
  else if(mePartonData()[0]->id()<=6&&
	  mePartonData()[1]->id()<=6) {
    if(_process == 0) { ol_evaluate_loop2(_id4, pp, &m2_tree, &acc); }
    if(_process == 1) { ol_evaluate_loop2(_id4tri, pp, &m2_tree, &acc); }
    if(_process == 2) { ol_evaluate_loop2(_id4box, pp, &m2_tree, &acc); }
    //cout << "qqb->gHH" << endl;
	
  }
  mesq = m2_tree;
  /*cout << "--" << endl;
  cout << "m2_tree = " << m2_tree << endl;
  cout << "mesq*sHat()/GeV/GeV = " << mesq*sHat()/GeV/GeV << endl;
  cout << "_process = " << _process << endl;
  cout << "printing ids" << endl;
  cout << _id1 << " " << _id1tri << " " << _id1box << " " << _id1int <<" " <<  _id2 << " " << _id2tri <<" " <<  _id2box << " " << _id2int << " " << _id3 <<" " <<  _id3tri <<" " <<  _id3box << " " << _id3int <<" " <<  _id4 << " " << _id4tri <<" " <<  _id4box << " " << _id4int << endl; */
  return mesq*sHat()/GeV/GeV;
}

bool MEHiggsPairJet::generateKinematics(const double * r) { 
  // initialize jacobian
  jacobian(1.);
  // cms energy
  Energy ecm=sqrt(sHat());
  // first generate the mass of the off-shell gauge boson
  // minimum mass of the 
  tcPDVector ptemp;
  ptemp.push_back(mePartonData()[3]);
  ptemp.push_back(mePartonData()[4]);
  Energy2 minMass2 = max(lastCuts().minSij(mePartonData()[3],mePartonData()[4]),
			 lastCuts().minS(ptemp));
  // minimum pt of the jet
  Energy ptmin = lastCuts().minKT(mePartonData()[2]);

  // maximum mass of the so pt is possible
  Energy2 maxMass2 = min(ecm*(ecm-2.*ptmin),lastCuts().maxS(ptemp));


  if(maxMass2<=ZERO||minMass2<ZERO) return false;
  
  // Energy mbt= sqrt(minMass2) + r[1] * (sqrt(maxMass2) - sqrt(minMass2));
  _mbt2=minMass2*maxMass2/(minMass2+r[1]*(maxMass2-minMass2));
  Energy mbt = sqrt(_mbt2);
  InvEnergy2 emjac1 = minMass2*maxMass2/(maxMass2-minMass2)/sqr(_mbt2);
  jacobian(jacobian()/sHat()/emjac1);
  // set the masses of the outgoing particles to 2-2 scattering
  meMomenta()[2].setMass(ZERO);
  Lorentz5Momentum pz(mbt);
  // generate the polar angle of the hard scattering
  double ctmin(-1.0), ctmax(1.0);
  Energy q(ZERO);
  try {
    q = SimplePhaseSpace::getMagnitude(sHat(), meMomenta()[2].mass(),mbt);
  } 
  catch ( ImpossibleKinematics ) {
    return false;
  }	    
  Energy2 pq = sqrt(sHat())*q;
  if ( ptmin > ZERO ) {
    double ctm = 1.0 - sqr(ptmin/q);
    if ( ctm <= 0.0 ) return false;
    ctmin = max(ctmin, -sqrt(ctm));
    ctmax = min(ctmax,  sqrt(ctm));
  }
  if ( ctmin >= ctmax ) return false;
  double cth = getCosTheta(ctmin, ctmax, r[0]);
  Energy pt  = q*sqrt(1.0-sqr(cth));
  double phi = 2.0*Constants::pi*r[2];
  meMomenta()[2].setVect(Momentum3( pt*sin(phi), pt*cos(phi), q*cth));
  pz.setVect(            Momentum3(-pt*sin(phi),-pt*cos(phi),-q*cth));
  meMomenta()[2].rescaleEnergy();
  pz.rescaleEnergy();
 


  //cout << LLSudakov(0.*GeV, sqrt(_scale), 0.*GeV, 22) << endl;
  //cout << "sqrt(_scale)/pt = " << sqrt(_scale)/GeV << " " << pt/GeV << endl;
  // generate the momenta of the fake resonance decay products
  meMomenta()[3].setMass(mePartonData()[3]->mass());
  meMomenta()[4].setMass(mePartonData()[4]->mass());
  Energy q2 = ZERO;
  try {
    q2 = SimplePhaseSpace::getMagnitude(_mbt2, meMomenta()[3].mass(),
					meMomenta()[4].mass());
  } catch ( ImpossibleKinematics ) {
    return false;
  }
  
  // set the scale
  if(_fixedscale==0) { _scale = sqr(_scalemultiplier)*sHat(); } 
  else if(_fixedscale==1) { _scale = sqr(_scalemultiplier)*(sqr(_basescale) + sqr(pt)); }
  else if(_fixedscale==2) { _scale = sqr(_scalemultiplier)*sqr(_basescale); }
  else if(_fixedscale==3) { _scale = sqr(_scalemultiplier)*sqr(_basescale + pt); }
  else if(_fixedscale==5) { _scale= sqr(_scalemultiplier)*_mbt2; }

  double cth2 =-1.+2.*r[3];
  double phi2=Constants::twopi*r[4];
  Energy pt2 =q2*sqrt(1.-sqr(cth2));
  Lorentz5Momentum pl[2]={Lorentz5Momentum( pt2*cos(phi2), pt2*sin(phi2), q2*cth2,ZERO,
					    meMomenta()[3].mass()),
			  Lorentz5Momentum(-pt2*cos(phi2),-pt2*sin(phi2),-q2*cth2,ZERO,
					   meMomenta()[4].mass())};
  pl[0].rescaleEnergy();
  pl[1].rescaleEnergy();
  Boost boostv(pz.boostVector());
  pl[0].boost(boostv);
  pl[1].boost(boostv);
  meMomenta()[3] = pl[0];
  meMomenta()[4] = pl[1];
  // check passes all the cuts
  vector<LorentzMomentum> out(3);
  tcPDVector tout(3);
  for(unsigned int ix=0;ix<3;++ix) {
    out[ ix] = meMomenta()[ix+2];
    tout[ix] = mePartonData()[ix+2];
  }
  if ( !lastCuts().passCuts(tout, out, mePartonData()[0], mePartonData()[1]) )
    return false;

  /* alphaS reweighting peformed here,
   * if chosen by AlphaSreweighting interface
   * via _alphasreweight.
   */
  if(_alphasreweight == 1) {					       
    //  cout << sqr(pt) << endl;
    double alphasFactor = sqrt( SM().alphaS(sqr(pt))/SM().alphaS(scale()) );
    //cout << pt << "\t" << sqrt(scale()) << "\t" << alphasFactor << endl;
    jacobian(jacobian()*alphasFactor);
  }

  // jacobian
  jacobian((pq/sHat())*Constants::pi*jacobian()/8./sqr(Constants::pi)*q2/mbt);
  return true;
}

void MEHiggsPairJet::setKinematics() {
  HwMEBase::setKinematics();
}
void MEHiggsPairJet::InitOpenLoops() { 
 //  ol_setparameter_int("order_ew", 2);
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


  _id1 = ol_register_process("2 21 -> 2 25 25", 12);
  _id2 = ol_register_process("-2 21 -> -2 25 25", 12);
  _id3 = ol_register_process("21 21 -> 21 25 25", 12);
  _id4 = ol_register_process("2 -2 -> 21 25 25", 12);
  char only3h[] = "only3h";
  ol_setparameter_string("approx", only3h);
  _id1tri = ol_register_process("2 21 -> 2 25 25", 12);
  _id2tri = ol_register_process("-2 21 -> -2 25 25", 12);
  _id3tri = ol_register_process("21 21 -> 21 25 25", 12);
  _id4tri = ol_register_process("2 -2 -> 21 25 25", 12);
  char no3h[] = "no3h";
  ol_setparameter_string("approx", no3h);
  _id1box = ol_register_process("2 21 -> 2 25 25", 12);
  _id2box = ol_register_process("-2 21 -> -2 25 25", 12);
  _id3box = ol_register_process("21 21 -> 21 25 25", 12);
  _id4box = ol_register_process("2 -2 -> 21 25 25", 12); 
  char interf3h[] = "interf3h";
  ol_setparameter_string("approx", interf3h);
  _id1int = ol_register_process("2 21 -> 2 25 25", 12);
  _id2int = ol_register_process("-2 21 -> -2 25 25", 12);
  _id3int = ol_register_process("21 21 -> 21 25 25", 12);
  _id4int = ol_register_process("2 -2 -> 21 25 25", 12); 
    
  
 // Initialize OpenLoops
 ol_start();*/
}
