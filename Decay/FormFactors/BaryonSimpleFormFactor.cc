// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BaryonSimpleFormFactor class.
//

#include "BaryonSimpleFormFactor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "BaryonSimpleFormFactor.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

void BaryonSimpleFormFactor::doinit() throw(InitException) {
  BaryonFormFactor::doinit();
  _f1.resize(0);_f2.resize(0);_g1.resize(0);_g2.resize(0);
  // calculate the couplings for the different modes
  int id0,id1;
  double root23(sqrt(2./3.)),root2(sqrt(2)),root32(sqrt(3./2.));
  for(unsigned int ix=0;ix<numberOfFactors();++ix)
    {
      // get the particle ids for the mode
      particleID(ix,id0,id1);
      // work out the couplings
      // neutron beta decay
      if(id0==2112&&id1==2212)
	{
	  _f1.push_back(1.) ;_g1.push_back(_gA);
	  _f2.push_back(3.7*_gA-1.0);_g2.push_back(0.);
	}
      // sigma+ to Lambda
      else if(id0==3222&&id1==3122)
	{
	  _f1.push_back(0.);_g1.push_back(_gA*root23*_alphaD);
	  _f2.push_back(4.55*_g1.back()-1.0*_f1.back());
	  _g2.push_back(-0.03*_g1.back());
	}
      // sigma
      else if(id0==3112&&id1==3122)
	{
	  _f1.push_back(0.);_g1.push_back(_gA*root23*_alphaD);
	  _f2.push_back(4.55*_g1.back()-1.0*_f1.back());
	  _g2.push_back(-0.03*_g1.back());
	}
      else if(id0==3112&&id1==3212)
	{
	  _f1.push_back(root2);_g1.push_back(_gA*root2*(1.-_alphaD));
	  _f2.push_back(4.69*_g1.back()-_f1.back());_g2.push_back(0.);
	}
      else if(id0==3212&&id1==3222)
	{
	  _f1.push_back(-root2);_g1.push_back(-_gA*root2*(1.-_alphaD));
	  _f2.push_back(4.69*_g1.back()-_f1.back());_g2.push_back(0.);
	}
      else if(id0==3312&&id1==3322)
	{
	  _f1.push_back(-1.);_g1.push_back(-_gA*(1.-2.*_alphaD));
	  _f2.push_back(5.21*_g1.back()-_f1.back());_g2.push_back(0.);
	}
      else if(id0==3122&&id1==2212)
	{
	  _f1.push_back(-root32*_etaV);
	  _g1.push_back(-_gA*root32*_etaA*(1.-2.*_alphaD/3.));
	  _f2.push_back(4.05*_rhoM*_g1.back()-1.01*_f1.back());
	  _g2.push_back((4.05*_rhoE-0.09)*_g1.back());
	}
      else if(id0==3212&&id1==2212)
	{
	  _f1.push_back(-_etaV/root2);
	  _g1.push_back(-_gA/root2*_etaA*(1.-2.*_alphaD));
	  _f2.push_back(4.2*_rhoM*_g1.back()-1.03*_f1.back());
	  _g2.push_back((4.2*_rhoE-0.12)*_g1.back());
	}
      else if(id0==3112&&id1==2112)
	{
	  _f1.push_back(-_etaV);
	  _g1.push_back(-_gA*_etaA*(1.-2.*_alphaD));
	  _f2.push_back(4.2*_rhoM*_g1.back()-1.03*_f1.back());
	  _g2.push_back((4.2*_rhoE-0.12)*_g1.back());
	}
      else if(id0==3312&&id1==3122)
	{
	  _f1.push_back(root32*_etaV);
	  _g1.push_back(_gA*root32*_etaA*(1.-4./3.*_alphaD));
	  _f2.push_back(4.8*_rhoM*_g1.back()-1.01*_f1.back());
	  _g2.push_back((4.8*_rhoE-0.08)*_g1.back());
	}
      else if(id0==3312&&id1==3212)
	{
	  _f1.push_back(_etaV/root2);
	  _g1.push_back(_gA*_etaA/root2);
	  _f2.push_back(4.95*_rhoM*_g1.back()-1.*_f1.back());
	  _g2.push_back((4.95*_rhoE-0.05)*_g1.back());
	}
      else if(id0==3322&&id1==3222)
	{
	  _f1.push_back(_etaV);
	  _g1.push_back(_gA*_etaA);
	  _f2.push_back(4.95*_rhoM*_g1.back()-1.*_f1.back());
	  _g2.push_back((4.95*_rhoE-0.05)*_g1.back());
	}
      else
	{throw InitException() << "Mode not recognised in BaryonSimpleFormFactor "
			       << id0 << "  " << id1 
			       << Exception::abortnow;}
    }
}

BaryonSimpleFormFactor::~BaryonSimpleFormFactor() {}

void BaryonSimpleFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _gA << _alphaD << _etaV << _etaA << _rhoE << _rhoM 
     << _f1 << _f2 << _g1 << _g2;}

void BaryonSimpleFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _gA >> _alphaD >> _etaV >> _etaA >> _rhoE >> _rhoM 
     >> _f1 >> _f2 >> _g1 >> _g2;}

ClassDescription<BaryonSimpleFormFactor> BaryonSimpleFormFactor::initBaryonSimpleFormFactor;
// Definition of the static class description member.

void BaryonSimpleFormFactor::Init() {

  static ClassDocumentation<BaryonSimpleFormFactor> documentation
    ("The BaryonSimpleFormFactor class implements the"
     " quark model calculation of the form-factors from PRD25, 206");

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
			   Complex & f1a,Complex & f2a,Complex & f3a)
{
  f1v =  _f1[iloc];
  f1a = -_g1[iloc];
  f2v = -_f2[iloc];
  f2a =  _g2[iloc];
  f3v = 0.;
  f3a = 0.;
}

void BaryonSimpleFormFactor::dataBaseOutput(ofstream& output,bool header,
					    bool create) const 
{
  if(header){output << "update decayers set parameters=\"";}
  if(create)
    {output << "create Herwig++::BaryonSimpleFormFactor " << fullName() << " \n";}
  output << "set " << fullName() << ":g_A " <<  _gA << " \n";
  output << "set " << fullName() << ":alpha_D " << _alphaD  << " \n";
  output << "set " << fullName() << ":eta_V " << _etaV  << " \n";
  output << "set " << fullName() << ":eta_A " <<  _etaA << " \n";
  output << "set " << fullName() << ":rho_E " << _rhoE  << " \n";
  output << "set " << fullName() << ":rho_M " << _rhoM  << " \n";
  BaryonFormFactor::dataBaseOutput(output,false,false);
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}
}

