// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ChengHeavyBaryonFormFactor class.
//

#include "ChengHeavyBaryonFormFactor.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ClassDocumentation.h" 

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ChengHeavyBaryonFormFactor.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

ChengHeavyBaryonFormFactor::~ChengHeavyBaryonFormFactor() {}
  
void ChengHeavyBaryonFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _mu << _md << _ms << _mc << _mb << _Nfi << _eta << _f1 << _f2 << _f3 
     << _g1 << _g2 << _g3 << _mVbc << _mVbs << _mVcs << _mVbd << _mVcu << _mAbc 
     << _mAbs << _mAcs << _mAbd << _mAcu;
 }
  

void ChengHeavyBaryonFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _mu >> _md >> _ms >> _mc >> _mb >> _Nfi >> _eta >> _f1 >> _f2 >> _f3 
     >> _g1 >> _g2 >> _g3 >> _mVbc >> _mVbs >> _mVcs >> _mVbd >> _mVcu >> _mAbc 
     >> _mAbs >> _mAcs >> _mAbd >> _mAcu;
}

ClassDescription<ChengHeavyBaryonFormFactor> ChengHeavyBaryonFormFactor::initChengHeavyBaryonFormFactor;
// Definition of the static class description member.

void ChengHeavyBaryonFormFactor::Init() {

  static ClassDocumentation<ChengHeavyBaryonFormFactor> documentation
    ("The \\classname{ChengHeavyBaryonFormFactor} class is the implementation"
     " of the form-factors of PRD53, 1457 and PRD56, 2799 for the weak decay of"
     "baryons containing a heavy quark. This model can be used for either"
     "semi-leptonic decays, or with the factorization approximation for"
     " non-leptonic weak decays");

  static Parameter<ChengHeavyBaryonFormFactor,Energy> interfaceUpMass
    ("DownMass",
     "The consituent mass of the down quark",
     &ChengHeavyBaryonFormFactor::_md, GeV, 0.332*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ChengHeavyBaryonFormFactor,Energy> interfaceDownMass
    ("UpMass",
     "The consituent mass of the up quark",
     &ChengHeavyBaryonFormFactor::_mu, GeV, 0.338*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ChengHeavyBaryonFormFactor,Energy> interfacStrangeMass
    ("StrangeMass",
     "The consituent mass of the strange quark",
     &ChengHeavyBaryonFormFactor::_ms, GeV, 0.510*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ChengHeavyBaryonFormFactor,Energy> interfaceCharmMass
    ("CharmMass",
     "The consituent mass of the charm quark",
     &ChengHeavyBaryonFormFactor::_mc, GeV, 1.6*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ChengHeavyBaryonFormFactor,Energy> interfaceBottomMass
    ("BottomMass",
     "The consituent mass of the bottom quark",
     &ChengHeavyBaryonFormFactor::_mb, GeV, 5.0*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ChengHeavyBaryonFormFactor,Energy> interfaceVectorMassbc
    ("VectorMassbc",
     "The vector mass for the b->c transitions.",
     &ChengHeavyBaryonFormFactor::_mVbc, GeV, 6.34*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ChengHeavyBaryonFormFactor,Energy> interfaceAxialMassbc
    ("AxialMassbc",
     "The axial-vector mass for the b->c transitions.",
     &ChengHeavyBaryonFormFactor::_mAbc, GeV, 6.73*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ChengHeavyBaryonFormFactor,Energy> interfaceVectorMassbs
    ("VectorMassbs",
     "The vector mass for the b->s transitions.",
     &ChengHeavyBaryonFormFactor::_mVbs, GeV, 5.42*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ChengHeavyBaryonFormFactor,Energy> interfaceAxialMassbs
    ("AxialMassbs",
     "The axial-vector mass for the b->s transitions.",
     &ChengHeavyBaryonFormFactor::_mAbs, GeV, 5.86*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ChengHeavyBaryonFormFactor,Energy> interfaceVectorMassbd
    ("VectorMassbd",
     "The vector mass for the b->d transitions.",
     &ChengHeavyBaryonFormFactor::_mVbd, GeV, 5.32*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ChengHeavyBaryonFormFactor,Energy> interfaceAxialMassbd
    ("AxialMassbd",
     "The axial-vector mass for the b->d transitions.",
     &ChengHeavyBaryonFormFactor::_mAbd, GeV, 5.71*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ChengHeavyBaryonFormFactor,Energy> interfaceVectorMasscs
    ("VectorMasscs",
     "The vector mass for the c->s transitions.",
     &ChengHeavyBaryonFormFactor::_mVcs, GeV, 2.11*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ChengHeavyBaryonFormFactor,Energy> interfaceAxialMasscs
    ("AxialMasscs",
     "The axial-vector mass for the c->s transitions.",
     &ChengHeavyBaryonFormFactor::_mAcs, GeV, 2.54*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ChengHeavyBaryonFormFactor,Energy> interfaceVectorMasscu
    ("VectorMasscu",
     "The vector mass for the c->u transitions.",
     &ChengHeavyBaryonFormFactor::_mVcu, GeV, 2.01*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ChengHeavyBaryonFormFactor,Energy> interfaceAxialMasscu
    ("AxialMasscu",
     "The axial-vector mass for the c->u transitions.",
     &ChengHeavyBaryonFormFactor::_mAcu, GeV, 2.42*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

}

// form factor for spin-1/2 to spin-1/2
void ChengHeavyBaryonFormFactor::
SpinHalfSpinHalfFormFactor(Energy2 q2,int iloc,int id0,int id1, Energy m0,Energy m1,
			   Complex & f1v,Complex & f2v,Complex & f3v,
			   Complex & f1a,Complex & f2a,Complex & f3a)
{
  id0=abs(id0);
  id1=abs(id1);
  // masses for the energy dependence of the form-factors
  Energy mV,mA;
   if((id0==4122&&id1==3122)||(id0==4232&&id1==3322)||(id0==4132&&id1==3312)||
      (id0==4332&&id1==3334))
     {mA=_mAcs;mV=_mVcs;}
   else if((id0==4332&&id1==3322)||(id0==4332&&id1==3324))
     {mA=_mAcu;mV=_mVcu;}
   else if((id0==5122&&id1==4122)||(id0==5232&&id1==4232)||(id0==5132&&id1==4132)||
	   (id0==5332&&id1==4332)||(id0==5332&&id1==4334))
     {mA=_mAbc;mV=_mVbc;}
   else if((id0==5122&&id1==3122)||(id0==5132&&id1==3312)||(id0==5232&&id1==3322)||
	   (id0==5332&&id1==3334))
     {mA=_mAbs;mV=_mVbs;}
   else if((id0==5122&&id1==2112)||(id0==5132&&id1==3112)||(id0==5232&&id1==3212)||
	   (id0==5332&&id1==3312)||(id0==5232&&id1==3122)||(id0==5332&&id1==3314))
     {mA=_mAbd;mV=_mVbd;}
   Energy delta=m0-m1;
   double Vfact = (1.-delta*delta/mV/mV)/(1.-q2/mV/mV);
   Vfact *=Vfact;
   double Afact = (1.-delta*delta/mA/mA)/(1.-q2/mA/mA);
   Afact *=Afact;
   f1v = _f1[iloc]*Vfact;
   f2v = _f2[iloc]*Vfact;
   f3v = _f3[iloc]*Vfact;
   f1a =-_g1[iloc]*Afact;
   f2a =-_g2[iloc]*Afact;
   f3a =-_g3[iloc]*Afact;
}

// form factor for spin-1/2 to spin-3/2
void ChengHeavyBaryonFormFactor::
SpinHalfSpinThreeHalfFormFactor(Energy2 q2,int iloc, int id0, int id1,
				Energy m0, Energy m1,
				Complex & g1v,Complex & g2v,Complex & g3v,
				Complex & g4v,Complex & g1a,Complex & g2a,
				Complex & g3a,Complex & g4a)
{
  id0=abs(id0);
  id1=abs(id1);
  // masses for the energy dependence of the form-factors
  Energy mV,mA;
   if((id0==4122&&id1==3122)||(id0==4232&&id1==3322)||(id0==4132&&id1==3312)||
      (id0==4332&&id1==3334))
     {mA=_mAcs;mV=_mVcs;}
   else if((id0==4332&&id1==3322)||(id0==4332&&id1==3324))
     {mA=_mAcu;mV=_mVcu;}
   else if((id0==5122&&id1==4122)||(id0==5232&&id1==4232)||(id0==5132&&id1==4132)||
	   (id0==5332&&id1==4332)||(id0==5332&&id1==4334))
     {mA=_mAbc;mV=_mVbc;}
   else if((id0==5122&&id1==3122)||(id0==5132&&id1==3312)||(id0==5232&&id1==3322)||
	   (id0==5332&&id1==3334))
     {mA=_mAbs;mV=_mVbs;}
   else if((id0==5122&&id1==2112)||(id0==5132&&id1==3112)||(id0==5232&&id1==3212)||
	   (id0==5332&&id1==3312)||(id0==5232&&id1==3122)||(id0==5332&&id1==3314))
     {mA=_mAbd;mV=_mVbd;}
   Energy delta=m0-m1;
   double Vfact = (1.-delta*delta/mV/mV)/(1.-q2/mV/mV);
   Vfact *=Vfact;
   double Afact = (1.-delta*delta/mA/mA)/(1.-q2/mA/mA);
   Afact *=Afact;
   g1v = _f1[iloc]*Vfact;
   g2v = _f2[iloc]*Vfact;
   g3v = _f3[iloc]*Vfact;
   g4v = 0.;
   g1a =-_g1[iloc]*Afact;
   g2a =-_g2[iloc]*Afact;
   g3a =-_g3[iloc]*Afact;
   g4a = 0.;
}

// output the information for the database
void ChengHeavyBaryonFormFactor::dataBaseOutput(ofstream& output)
{
  output << "create /Herwig++/ChengHeavyBaryonFormFactor " << fullName() << " \n";
  output << "set " << fullName() << ":DownMass     " << _md/GeV << " \n";
  output << "set " << fullName() << ":UpMass       " << _mu/GeV << " \n";
  output << "set " << fullName() << ":StrangeMass  " << _ms/GeV << " \n";
  output << "set " << fullName() << ":CharmMass    " << _mc/GeV << " \n";
  output << "set " << fullName() << ":BottomMass   " << _mb/GeV << " \n";
  output << "set " << fullName() << ":VectorMassbc " << _mVbc/GeV << " \n";
  output << "set " << fullName() << ":AxialMassbc  " << _mAbc/GeV << " \n";
  output << "set " << fullName() << ":VectorMassbs " << _mVbs/GeV << " \n";
  output << "set " << fullName() << ":AxialMassbs  " << _mAbs/GeV << " \n";
  output << "set " << fullName() << ":VectorMassbd " << _mVbd/GeV << " \n";
  output << "set " << fullName() << ":AxialMassbd  " << _mAbd/GeV << " \n";
  output << "set " << fullName() << ":VectorMasscs " << _mVcs/GeV << " \n";
  output << "set " << fullName() << ":AxialMasscs  " << _mAcs/GeV << " \n";
  output << "set " << fullName() << ":VectorMasscu " << _mVcu/GeV << " \n";
  output << "set " << fullName() << ":AxialMasscu  " << _mAcu/GeV << " \n";
}

}
