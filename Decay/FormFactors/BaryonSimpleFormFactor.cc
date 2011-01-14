// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BaryonSimpleFormFactor class.
//

#include "BaryonSimpleFormFactor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

BaryonSimpleFormFactor::BaryonSimpleFormFactor() {
  // axial vector coupling in beta decay
  _gA=1.25;
  // SU(3) breaking paramters
  // D/(D+F) ratio
  _alphaD=0.6;
  // eta_V parameter
  _etaV=0.970;
  // eta_A parameter
  _etaA=1.080;
  // SU(3) breaking factors for the electric dipole moment
  _rhoE=0.094;
  // SU(3) breaking factors for the magnetic dipole moment
  _rhoM=0.860;
  // the various decay modes handled by the model
  addFormFactor(2112,2212,2,2,2,1,1,2);
  addFormFactor(3222,3122,2,2,3,2,1,2);
  addFormFactor(3112,3122,2,2,3,1,1,2);
  addFormFactor(3112,3212,2,2,3,1,1,2);
  addFormFactor(3212,3222,2,2,3,2,1,2);
  addFormFactor(3312,3322,2,2,3,3,1,2);
  addFormFactor(3122,2212,2,2,2,1,3,2);
  addFormFactor(3212,2212,2,2,2,1,3,2);
  addFormFactor(3112,2112,2,2,1,1,3,2);
  addFormFactor(3312,3122,2,2,3,1,3,2);
  addFormFactor(3312,3212,2,2,3,1,3,2);
  addFormFactor(3322,3222,2,2,3,2,3,2);
  // set the inital number of form factors
  initialModes(numberOfFactors());
}

void BaryonSimpleFormFactor::doinit() {
  BaryonFormFactor::doinit();
  _f1.clear();
  _f2.clear();
  _g1.clear();
  _g2.clear();
  _f1.resize(numberOfFactors());
  _f2.resize(numberOfFactors());
  _g1.resize(numberOfFactors());
  _g2.resize(numberOfFactors());
  // calculate the couplings for the different modes
  int id0,id1;
  double root23(sqrt(2./3.)),root2(sqrt(2.)),root32(sqrt(3./2.));
  for(unsigned int ix=0;ix<numberOfFactors();++ix) {
    // get the particle ids for the mode
    particleID(ix,id0,id1);
    // work out the couplings
    // neutron beta decay
    if(id0==2112&&id1==2212) {
      _f1[ix] = 1. ;
      _g1[ix] = _gA;
      _f2[ix] = 3.7*_gA-1.0;
      _g2[ix] = 0.;
    }
    // sigma+ to Lambda
    else if(id0==3222&&id1==3122) {
      _f1[ix] = 0.;
      _g1[ix] = _gA*root23*_alphaD;
      _f2[ix] = 4.55*_g1[ix]-1.0*_f1[ix];
      _g2[ix] = -0.03*_g1[ix];
    }
    // sigma
    else if(id0==3112&&id1==3122) {
      _f1[ix] = 0.;
      _g1[ix] = _gA*root23*_alphaD;
      _f2[ix] = 4.55*_g1[ix]-1.0*_f1[ix];
      _g2[ix] = -0.03*_g1[ix];
    }
    else if(id0==3112&&id1==3212) {
      _f1[ix] = root2;
      _g1[ix] = _gA*root2*(1.-_alphaD);
      _f2[ix] = 4.69*_g1[ix]-_f1[ix];
      _g2[ix] = 0.;
    }
    else if(id0==3212&&id1==3222) {
      _f1[ix] = -root2;
      _g1[ix] = -_gA*root2*(1.-_alphaD);
      _f2[ix] = 4.69*_g1[ix]-_f1[ix];
      _g2[ix] = 0.;
    }
    else if(id0==3312&&id1==3322) {
      _f1[ix] = -1.;
      _g1[ix] = -_gA*(1.-2.*_alphaD);
      _f2[ix] = 5.21*_g1[ix]-_f1[ix];
      _g2[ix] = 0.;
    }
    else if(id0==3122&&id1==2212) {
      _f1[ix] = -root32*_etaV;
      _g1[ix] = -_gA*root32*_etaA*(1.-2.*_alphaD/3.);
      _f2[ix] = 4.05*_rhoM*_g1[ix]-1.01*_f1[ix];
      _g2[ix] = (4.05*_rhoE-0.09)*_g1[ix];
    }
    else if(id0==3212&&id1==2212) {
      _f1[ix] = -_etaV/root2;
      _g1[ix] = -_gA/root2*_etaA*(1.-2.*_alphaD);
      _f2[ix] = 4.2*_rhoM*_g1[ix]-1.03*_f1[ix];
      _g2[ix] = (4.2*_rhoE-0.12)*_g1[ix];
    }
    else if(id0==3112&&id1==2112) {
      _f1[ix] = -_etaV;
      _g1[ix] = -_gA*_etaA*(1.-2.*_alphaD);
      _f2[ix] = 4.2*_rhoM*_g1[ix]-1.03*_f1[ix];
      _g2[ix] = (4.2*_rhoE-0.12)*_g1[ix];
    }
    else if(id0==3312&&id1==3122) {
      _f1[ix] = root32*_etaV;
      _g1[ix] = _gA*root32*_etaA*(1.-4./3.*_alphaD);
      _f2[ix] = 4.8*_rhoM*_g1[ix]-1.01*_f1[ix];
      _g2[ix] = (4.8*_rhoE-0.08)*_g1[ix];
    }
    else if(id0==3312&&id1==3212) {
      _f1[ix] = _etaV/root2;
      _g1[ix] = _gA*_etaA/root2;
      _f2[ix] = 4.95*_rhoM*_g1[ix]-1.*_f1[ix];
      _g2[ix] = (4.95*_rhoE-0.05)*_g1[ix];
    }
    else if(id0==3322&&id1==3222) {
      _f1[ix] = _etaV;
      _g1[ix] = _gA*_etaA;
      _f2[ix] = 4.95*_rhoM*_g1[ix]-1.*_f1[ix];
      _g2[ix] = (4.95*_rhoE-0.05)*_g1[ix];
    }
    else {
      throw InitException() << "Mode not recognised in BaryonSimpleFormFactor "
			    << id0 << "  " << id1 
			    << Exception::abortnow;
    }
  }
}

void BaryonSimpleFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _gA << _alphaD << _etaV << _etaA << _rhoE << _rhoM 
     << _f1 << _f2 << _g1 << _g2;}

void BaryonSimpleFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _gA >> _alphaD >> _etaV >> _etaA >> _rhoE >> _rhoM 
     >> _f1 >> _f2 >> _g1 >> _g2;}

ClassDescription<BaryonSimpleFormFactor> 
BaryonSimpleFormFactor::initBaryonSimpleFormFactor;
// Definition of the static class description member.

void BaryonSimpleFormFactor::Init() {

  static ClassDocumentation<BaryonSimpleFormFactor> documentation
    ("The BaryonSimpleFormFactor class implements the"
     " quark model calculation of the form-factors from PRD25, 206",
     "The BaryonSimpleFormFactor class which implements the results"
     " of \\cite{Donoghue:1981uk}was used for the weak decays of the light baryons.",
     "\\bibitem{Donoghue:1981uk}\n"
     "J.~F.~Donoghue and B.~R.~Holstein,\n"
     "Phys.\\ Rev.\\  D {\\bf 25} (1982) 206.\n"
     "%%CITATION = PHRVA,D25,206;%%\n");

  static Parameter<BaryonSimpleFormFactor,double> interfaceg_A
    ("g_A",
     "The axial-vector coupling in neutron beta decay.",
     &BaryonSimpleFormFactor::_gA, 1.25, 0.0, 10.0,
     false, false, true);

  static Parameter<BaryonSimpleFormFactor,double> interfacealpha_D
    ("alpha_D",
     "SU(3) breaking parameter which is the ratio D/(D+F). ",
     &BaryonSimpleFormFactor::_alphaD, 0.6, 0.0, 1.0,
     false, false, true);

  static Parameter<BaryonSimpleFormFactor,double> interfaceeta_V
    ("eta_V",
     "The eta_V SU(3) breaking parameter",
     &BaryonSimpleFormFactor::_etaV, .97, 0.0, 10.0,
     false, false, true);

  static Parameter<BaryonSimpleFormFactor,double> interfaceeta_A
    ("eta_A",
     "The eta_A SU(3) breaking parameter",
     &BaryonSimpleFormFactor::_etaA, 1.08, 0.0, 10.0,
     false, false, true);

  static Parameter<BaryonSimpleFormFactor,double> interfacerho_E
    ("rho_E",
     "The SU(3) breaking parameter for the electric dipole moment.",
     &BaryonSimpleFormFactor::_rhoE, 0.094, 0.0, 10.0,
     false, false, true);

  static Parameter<BaryonSimpleFormFactor,double> interfacerho_M
    ("rho_M",
     "The SU(3) breaking parameter for the magentic dipole moment.",
     &BaryonSimpleFormFactor::_rhoM, 0.860, 0.0, 10.0,
     false, false, true);
}

void BaryonSimpleFormFactor::
SpinHalfSpinHalfFormFactor(Energy2,int iloc, int,int,Energy,Energy,
			   Complex & f1v,Complex & f2v,Complex & f3v,
			   Complex & f1a,Complex & f2a,Complex & f3a) {
  useMe();
  f1v =  _f1[iloc];
  f1a = -_g1[iloc];
  f2v = -_f2[iloc];
  f2a =  _g2[iloc];
  f3v = 0.;
  f3a = 0.;
}

void BaryonSimpleFormFactor::dataBaseOutput(ofstream& output,bool header,
					    bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::BaryonSimpleFormFactor " 
		    << name() << " \n";
  output << "newdef " << name() << ":g_A " <<  _gA << " \n";
  output << "newdef " << name() << ":alpha_D " << _alphaD  << " \n";
  output << "newdef " << name() << ":eta_V " << _etaV  << " \n";
  output << "newdef " << name() << ":eta_A " <<  _etaA << " \n";
  output << "newdef " << name() << ":rho_E " << _rhoE  << " \n";
  output << "newdef " << name() << ":rho_M " << _rhoM  << " \n";
  BaryonFormFactor::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
