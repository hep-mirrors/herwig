// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RadiativeHeavyBaryonDecayer class.
//

#include "RadiativeHeavyBaryonDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void RadiativeHeavyBaryonDecayer::doinitrun() {
  Baryon1MesonDecayerBase::doinitrun();
  if(initialize()) {
    _maxweight.clear();
    for(unsigned int ix=0;ix<numberModes();++ix) {
      if(mode(ix)) _maxweight.push_back(mode(ix)->maxWeight());
      else         _maxweight.push_back(1.);
    }
  }
}

RadiativeHeavyBaryonDecayer::RadiativeHeavyBaryonDecayer() {
  // Xi_c' to Xi_c
  _incoming.push_back(4312);_outgoingB.push_back(4132);_maxweight.push_back(2.1);
  _modetype.push_back(1);_e1coupling.push_back(ZERO);_m1coupling.push_back(1.1004e-5/MeV);
  _incoming.push_back(4322);_outgoingB.push_back(4232);_maxweight.push_back(2.1);
  _modetype.push_back(1);_e1coupling.push_back(ZERO);_m1coupling.push_back(9.4102e-5/MeV);
  // Xi_b' to Xi_b
  _incoming.push_back(5312);_outgoingB.push_back(5132);_maxweight.push_back(2.1);
  _modetype.push_back(1);_e1coupling.push_back(ZERO);_m1coupling.push_back(1.1004e-5/MeV);
  _incoming.push_back(5322);_outgoingB.push_back(5232);_maxweight.push_back(2.1);
  _modetype.push_back(1);_e1coupling.push_back(ZERO);_m1coupling.push_back(9.4102e-5/MeV);
  // sigma_c*+ to Lambda_c+
  _incoming.push_back(4214);_outgoingB.push_back(4122);_maxweight.push_back(0.05);
  _modetype.push_back(1);_e1coupling.push_back(ZERO);_m1coupling.push_back(3.76e-4/MeV);
  // sigma_b*+ to Lambda_b+
  _incoming.push_back(5214);_outgoingB.push_back(5122);_maxweight.push_back(0.05);
  _modetype.push_back(1);_e1coupling.push_back(ZERO);_m1coupling.push_back(3.76e-4/MeV);
  // Xi_c* to Xi_c
  _incoming.push_back(4314);_outgoingB.push_back(4132);_maxweight.push_back(0.005);
  _modetype.push_back(1);_e1coupling.push_back(ZERO);_m1coupling.push_back(0.374e-4/MeV);
  _incoming.push_back(4324);_outgoingB.push_back(4232);_maxweight.push_back(0.05);
  _modetype.push_back(1);_e1coupling.push_back(ZERO);_m1coupling.push_back(3.2486e-4/MeV);
  // Xi_b* to Xi_b
  _incoming.push_back(5314);_outgoingB.push_back(5132);_maxweight.push_back(0.005);
  _modetype.push_back(1);_e1coupling.push_back(ZERO);_m1coupling.push_back(1.0132e-4/MeV);
  _incoming.push_back(5324);_outgoingB.push_back(5232);_maxweight.push_back(0.05);
  _modetype.push_back(1);_e1coupling.push_back(ZERO);_m1coupling.push_back(3.2486e-4/MeV);
  // sigma_c+ to Lambda_c+
  _incoming.push_back(4212);_outgoingB.push_back(4122);_maxweight.push_back(0.1);
  _modetype.push_back(1);_e1coupling.push_back(ZERO);_m1coupling.push_back(1.088e-4/MeV);
  // sigma_b+ to Lambda_b+
  _incoming.push_back(5212);_outgoingB.push_back(5122);_maxweight.push_back(0.1);
  _modetype.push_back(1);_e1coupling.push_back(ZERO);_m1coupling.push_back(1.088e-4/MeV);
  // sigma_c* to sigma_c
  _incoming.push_back(4224);_outgoingB.push_back(4222);_maxweight.push_back(0.008);
  _modetype.push_back(1);_e1coupling.push_back(ZERO);_m1coupling.push_back( 3.337e-4/MeV);
  _incoming.push_back(4214);_outgoingB.push_back(4212);_maxweight.push_back(0.003);
  _modetype.push_back(1);_e1coupling.push_back(ZERO);_m1coupling.push_back( 0.834e-4/MeV);
  _incoming.push_back(4114);_outgoingB.push_back(4112);_maxweight.push_back(0.0012);
  _modetype.push_back(1);_e1coupling.push_back(ZERO);_m1coupling.push_back(-1.688e-4/MeV);
  // Xi_c* to Xi'_c
  _incoming.push_back(4314);_outgoingB.push_back(4312);_maxweight.push_back(0.012);
  _modetype.push_back(1);_e1coupling.push_back(ZERO);_m1coupling.push_back(-6.110e-4/MeV);
  _incoming.push_back(4324);_outgoingB.push_back(4322);_maxweight.push_back(0.006);
  _modetype.push_back(1);_e1coupling.push_back(ZERO);_m1coupling.push_back(-3.607e-4/MeV);
  // Omega_c* to Omega_c
  _incoming.push_back(4334);_outgoingB.push_back(4332);_maxweight.push_back(2.4);
  _modetype.push_back(1);_e1coupling.push_back(ZERO);_m1coupling.push_back( 1.055e-3/MeV);
  // sigma_b* to sigma_b
  _incoming.push_back(5224);_outgoingB.push_back(5222);_maxweight.push_back(0.06);
  _modetype.push_back(1);_e1coupling.push_back(ZERO);_m1coupling.push_back( 9.945e-4/MeV);
  _incoming.push_back(5214);_outgoingB.push_back(5212);_maxweight.push_back(0.004);
  _modetype.push_back(1);_e1coupling.push_back(ZERO);_m1coupling.push_back( 2.486e-4/MeV);
  _incoming.push_back(5114);_outgoingB.push_back(5112);_maxweight.push_back(0.015);
  _modetype.push_back(1);_e1coupling.push_back(ZERO);_m1coupling.push_back(-4.973e-4/MeV);
  // Xi_b* to Xi'_b
  _incoming.push_back(5314);_outgoingB.push_back(5312);_maxweight.push_back(0.007);
  _modetype.push_back(1);_e1coupling.push_back(ZERO);_m1coupling.push_back(-4.371e-4/MeV);
  _incoming.push_back(5324);_outgoingB.push_back(5322);_maxweight.push_back(0.005);
  _modetype.push_back(1);_e1coupling.push_back(ZERO);_m1coupling.push_back( 3.088e-4/MeV);
  // Omega_b* to Omega_b
  _incoming.push_back(5334);_outgoingB.push_back(5332);_maxweight.push_back(2.1);
  _modetype.push_back(1);_e1coupling.push_back(ZERO);_m1coupling.push_back(-3.769e-4/MeV);
  // initial size of the arrays
  _initsize=_incoming.size();
  // intermediates
  generateIntermediates(false);
}

void RadiativeHeavyBaryonDecayer::doinit() {
  Baryon1MesonDecayerBase::doinit();
  // check the parameters are consistent
  unsigned int isize(_incoming.size());
  if(isize!=_outgoingB.size() ||isize!=_maxweight.size()||isize!=_m1coupling.size()||
     isize!=_e1coupling.size()||isize!=_modetype.size())
    throw InitException() << "Inconsistent parameters in "
			  << "RadiativeHeavyBaryonDecayer::doinit()" 
			  << Exception::abortnow;
  vector<double> wgt(0);
  tPDVector extpart(3);
  DecayPhaseSpaceModePtr mode;
  // the decay modes
  extpart[2]=getParticleData(ParticleID::gamma);
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    extpart[0]=getParticleData(_incoming[ix]);
    extpart[1]=getParticleData(_outgoingB[ix]);
    if(extpart[0]&&extpart[1]) {
      mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
      addMode(mode,_maxweight[ix],wgt);
    }
    else {
      addMode(DecayPhaseSpaceModePtr(),_maxweight[ix],wgt);
    }
  }
}

int RadiativeHeavyBaryonDecayer::modeNumber(bool & cc,tcPDPtr parent,
					    const tPDVector & children) const {
  int imode(-1);
  // must be two outgoing particles
  if(children.size()!=2){return imode;}
  // ids of the particles
  int id0(parent->id());
  int id1(children[0]->id());
  int id2(children[1]->id()),ibaryon;
  if(id1==ParticleID::gamma){ibaryon=id2;}
  else if(id2==ParticleID::gamma){ibaryon=id1;}
  else {return imode;}
  unsigned int ix(0);
  do {
    if(     id0== _incoming[ix]&&ibaryon== _outgoingB[ix]) {
      imode=ix;
      cc=false;
    }
    else if(id0==-_incoming[ix]&&ibaryon==-_outgoingB[ix]) {
      imode=ix;
      cc=true;
    }
    ++ix;
  }
  while(ix<_incoming.size()&&imode<0);
  return imode;
}

void RadiativeHeavyBaryonDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(_m1coupling,1./GeV) << ounit(_e1coupling,1./GeV2) 
     << _incoming << _outgoingB << _modetype << _maxweight;
}

void RadiativeHeavyBaryonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_m1coupling,1./GeV) >> iunit(_e1coupling,1./GeV2) 
     >> _incoming >> _outgoingB >> _modetype >> _maxweight;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<RadiativeHeavyBaryonDecayer,Baryon1MesonDecayerBase>
describeHerwigRadiativeHeavyBaryonDecayer("Herwig::RadiativeHeavyBaryonDecayer", "HwBaryonDecay.so");

void RadiativeHeavyBaryonDecayer::Init() {

  static ClassDocumentation<RadiativeHeavyBaryonDecayer> documentation
    ("The RadiativeHeavyBaryonDecayer class is designed for the radiative decays of"
     " heavy baryons.",
     "The radiative decays of the heavy baryons were simulated using the results of"
     "\\cite{Ivanov:1999bk,Ivanov:1998wj}.",
     "\\bibitem{Ivanov:1999bk}\n"
     "M.~A.~Ivanov, J.~G.~Korner, V.~E.~Lyubovitskij and A.~G.~Rusetsky,\n"
     "Phys.\\ Rev.\\  D {\\bf 60} (1999) 094002\n"
     "[arXiv:hep-ph/9904421].\n"
     "%%CITATION = PHRVA,D60,094002;%%\n"
     "\\bibitem{Ivanov:1998wj}\n"
     "M.~A.~Ivanov, J.~G.~Korner and V.~E.~Lyubovitskij,\n"
     "Phys.\\ Lett.\\  B {\\bf 448} (1999) 143 [arXiv:hep-ph/9811370].\n"
     "%%CITATION = PHLTA,B448,143;%%\n");

  static ParVector<RadiativeHeavyBaryonDecayer,InvEnergy> interfaceM1Coupling
    ("M1Coupling",
     "The coupling for the M1 modes",
     &RadiativeHeavyBaryonDecayer::_m1coupling, 
     1./MeV, -1, ZERO, ZERO, ZERO,
     false, false, false);

  static ParVector<RadiativeHeavyBaryonDecayer,InvEnergy2> interfaceE1Coupling
    ("E1Coupling",
     "The coupling for the E1 modes",
     &RadiativeHeavyBaryonDecayer::_e1coupling,
     1./MeV/MeV, -1,  ZERO, ZERO, ZERO,
     false, false, false);

  static ParVector<RadiativeHeavyBaryonDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code of the incoming baryon",
     &RadiativeHeavyBaryonDecayer::_incoming, -1, 0, 0, 1000000,
     false, false, true);

  static ParVector<RadiativeHeavyBaryonDecayer,int> interfaceOutgoingB
    ("OutgoingB",
     "The PDG code of the outgoing baryon",
     &RadiativeHeavyBaryonDecayer::_outgoingB, -1, 0, 0, 1000000,
     false, false, true);

  static ParVector<RadiativeHeavyBaryonDecayer,int> interfaceModeType
    ("ModeType",
     "The type of mode. 0 is E1, 1 is M1",
     &RadiativeHeavyBaryonDecayer::_modetype, -1, 0, 0, 2,
     false, false, true);

  static ParVector<RadiativeHeavyBaryonDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &RadiativeHeavyBaryonDecayer::_maxweight,
     0, 0, 0, 0., 100., false, false, true);
}

void RadiativeHeavyBaryonDecayer::dataBaseOutput(ofstream & output,bool header) const {
  if(header) output << "update decayers set parameters=\"";
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    if(ix<_initsize) {
      output << "newdef " << name() << ":M1Coupling " 
	     << ix << " " << _m1coupling[ix]*MeV << "\n";
      output << "newdef " << name() << ":E1Coupling " 
	     << ix << " " << _e1coupling[ix]*MeV2 << "\n";
      output << "newdef " << name() << ":Incoming " 
	     << ix << " " << _incoming[ix] << "\n";
      output << "newdef " << name() << ":OutgoingB " 
	     << ix << " " << _outgoingB[ix] << "\n";
      output << "newdef " << name() << ":ModeType " 
	     << ix << " " << _modetype[ix] << "\n";
      output << "newdef " << name() << ":MaxWeight " 
	     << ix << " " << _maxweight[ix] << "\n";
    }
    else {
      output << "insert " << name() << ":M1Coupling " 
	     << ix << " " << _m1coupling[ix]*MeV << "\n";
      output << "insert " << name() << ":E1Coupling " 
	     << ix << " " << _e1coupling[ix]*MeV2 << "\n";
      output << "insert " << name() << ":Incoming " 
	     << ix << " " << _incoming[ix] << "\n";
      output << "insert " << name() << ":OutgoingB " 
	     << ix << " " << _outgoingB[ix] << "\n";
      output << "insert " << name() << ":ModeType " 
	     << ix << " " << _modetype[ix] << "\n";
      output << "insert " << name() << ":MaxWeight " 
	     << ix << " " << _maxweight[ix] << "\n";
    }
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

void RadiativeHeavyBaryonDecayer::halfHalfVectorCoupling(int imode,Energy m0,Energy m1,
							 Energy,Complex& A1,
							 Complex& A2,Complex& B1,
							 Complex& B2) const {
  useMe();
  if(_modetype[imode]==0) {
    B1=-0.5*(sqr(m0) - sqr(m1))*_e1coupling[imode];
    B2=-    (m0+m1)*(m0+m1)*_e1coupling[imode];
    A1=0.;
    A2=0.;
  }
  else if(_modetype[imode]==1) {
    A1=-2.*(m0+m1)*_m1coupling[imode];
    A2= 4.*(m0+m1)*_m1coupling[imode];
    B1=0.;
    B2=0.;
  }
  else {
    throw DecayIntegratorError() << "Unknown type of mode " << _modetype[imode] 
				 << " in RadiativeHeavyBaryonDecayer::"
				 << "halfHalfVectorCoupling()" << Exception::abortnow;
  }
}

void RadiativeHeavyBaryonDecayer::threeHalfHalfVectorCoupling(int imode,Energy m0,
							      Energy m1,Energy,
							      Complex& A1,Complex& A2,
							      Complex& A3,Complex& B1,
							      Complex& B2,
							      Complex& B3) const {
  useMe();
  if(_modetype[imode]==0) {
    A1=-0.5*(m0*m0 - m1*m1)*_e1coupling[imode];
    A3=-    (m0+m1)*(m0+m1)*_e1coupling[imode];    
    A2=0.;
    B1=0.;
    B2=0.;
    B3=0.;
  }
  else if(_modetype[imode]==1) {
    B1=-(m0+m1)*_m1coupling[imode];
    B2= (m0+m1)*_m1coupling[imode];
    B3=0.;
    A1=0.;
    A2=0.;
    A3=0.;
  }
  else {
    throw DecayIntegratorError() << "Unknown type of mode " << _modetype[imode] 
				 << " in RadiativeHeavyBaryonDecayer::"
				 << "threeHalfHalfVectorCoupling()" 
				 << Exception::abortnow;
  }
}
