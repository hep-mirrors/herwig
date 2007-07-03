// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ChengHeavyBaryonFormFactor class.
//

#include "ChengHeavyBaryonFormFactor.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/ClassDocumentation.h" 

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ChengHeavyBaryonFormFactor.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;
ChengHeavyBaryonFormFactor::ChengHeavyBaryonFormFactor() 
{
  // consituent quark masses
  _mu = 338*MeV;
  _md = 322*MeV;
  _ms = 510*MeV;
  _mc = 1.6*GeV;
  _mb = 5.0*GeV;
  // masses for the q^2 dependence
  _mVbc = 6.34*GeV;
  _mVbs = 5.42*GeV;
  _mVcs = 2.11*GeV;
  _mVbd = 5.32*GeV;
  _mVcu = 2.01*GeV;
  _mAbc = 6.73*GeV;
  _mAbs = 5.86*GeV;
  _mAcs = 2.54*GeV;
  _mAbd = 5.71*GeV;
  _mAcu = 2.42*GeV;
  double one3(1./sqrt(3.)),one2(1./sqrt(2.));
  // lambda_b to lambda_c
  addFormFactor(5122,4122,2,2,1,2,5,4);_Nfi.push_back(1.     );_eta.push_back(1.);
  // lambda_b to lambda
  addFormFactor(5122,3122,2,2,1,2,5,3);_Nfi.push_back(one3   );_eta.push_back(1.);
  // lambda_b to n
  addFormFactor(5122,2112,2,2,1,2,5,1);_Nfi.push_back(one2   );_eta.push_back(1.);
  // xi_b to xi_c
  addFormFactor(5232,4232,2,2,2,3,5,4);_Nfi.push_back(1.     );_eta.push_back(1.);
  addFormFactor(5132,4132,2,2,1,3,5,4);_Nfi.push_back(1.     );_eta.push_back(1.);
  // xi_b to xi
  addFormFactor(5232,3322,2,2,2,3,5,3);_Nfi.push_back(one2   );_eta.push_back(1.);
  addFormFactor(5132,3312,2,2,1,3,5,3);_Nfi.push_back(one2   );_eta.push_back(1.);
  // xi_b to sigma
  addFormFactor(5232,3212,2,2,2,3,5,1);_Nfi.push_back(0.5    );_eta.push_back(1.);
  addFormFactor(5132,3112,2,2,1,3,5,1);_Nfi.push_back(0.5    );_eta.push_back(1.);
  // xi_b to lambda
  addFormFactor(5232,3122,2,2,2,3,5,1);_Nfi.push_back(one3/2.);_eta.push_back(1.);
  // omega_b to omega_c
  addFormFactor(5332,4332,2,2,3,3,5,4);_Nfi.push_back(1.     );_eta.push_back(-1./3.);
  // omega_b to xi
  addFormFactor(5332,3312,2,2,3,3,5,1);_Nfi.push_back(one3   );_eta.push_back(-1./3.);
  // omega_b to omega_c*
  addFormFactor(5332,4334,2,4,3,3,5,4);_Nfi.push_back(1.     );_eta.push_back(0.);
  // omega_b to omega
  addFormFactor(5332,3334,2,4,3,3,5,3);_Nfi.push_back(1.     );_eta.push_back(0.);
  // omega_b to xi*
  addFormFactor(5332,3314,2,4,3,3,5,1);_Nfi.push_back(one3   );_eta.push_back(0.);
  // omega_c to omega
  addFormFactor(4332,3334,2,4,3,3,4,3);_Nfi.push_back(1.     );_eta.push_back(0.);
  // omega_c to xi*
  addFormFactor(4332,3324,2,4,3,3,4,2);_Nfi.push_back(one3   );_eta.push_back(0.);
  // lambda_c to lambda_0
  addFormFactor(4122,3122,2,2,1,2,4,3);_Nfi.push_back(1./sqrt(3.));_eta.push_back(1.);
  // xi_c to xi
  addFormFactor(4232,3322,2,2,2,3,4,3);_Nfi.push_back(1./sqrt(3.));_eta.push_back(1.);
  addFormFactor(4132,3312,2,2,1,3,4,3);_Nfi.push_back(1./sqrt(3.));_eta.push_back(1.);
  // initial number of form factors
  initialModes(numberOfFactors());
}

void ChengHeavyBaryonFormFactor::doinit() throw(InitException) {
  BaryonFormFactor::doinit();
  // check the parameters are consistent
  unsigned int isize(numberOfFactors());
  if(isize!=_eta.size()||isize!=_Nfi.size())
    {throw InitException() << "Inconsistent paramters in ChengHeavyBaryon" 
			   << "FormFactor::doinit() " << Exception::abortnow;}
  Energy mi,mf,mq,mQ,lambda,delta,msum;
  int id0,id1,inspin,outspin,isp1,isp2,inq,outq;
  for(unsigned int ix=0;ix<numberOfFactors();++ix)
    {
      // ids of the external particles
      particleID(ix,id0,id1);
      formFactorInfo(ix,inspin,outspin,isp1,isp2,inq,outq);
      id0=abs(id0);id1=abs(id1);
      mi=getParticleData(id0)->mass();
      mf=getParticleData(id1)->mass();
      msum=mi+mf;
      // masses of the incoming and outgoing quarks
      if((id0==4122&&id1==3122)||(id0==4232&&id1==3322)||(id0==4132&&id1==3312)||
	 (id0==4332&&id1==3334))
	{mq=_ms;mQ=_mc;}
      else if((id0==4332&&id1==3322)||(id0==4332&&id1==3324))
	{mq=_mu;mQ=_mc;}
      else if((id0==5122&&id1==4122)||(id0==5232&&id1==4232)||(id0==5132&&id1==4132)||
	      (id0==5332&&id1==4332)||(id0==5332&&id1==4334))
	{mq=_mc;mQ=_mb;}
      else if((id0==5122&&id1==3122)||(id0==5132&&id1==3312)||(id0==5232&&id1==3322)||
	      (id0==5332&&id1==3334))
	{mq=_ms;mQ=_mb;}
      else if((id0==5122&&id1==2112)||(id0==5132&&id1==3112)||(id0==5232&&id1==3212)||
	      (id0==5332&&id1==3312)||(id0==5232&&id1==3122)||(id0==5332&&id1==3314))
	{mq=_md;mQ=_mb;}
      else
	{throw InitException() << "Unknown decay in ChengHeavyBaryon" 
			       << "FormFactor::doinit() " << Exception::abortnow;}
      // parameters
      lambda = mf-mq;
      delta  = mi-mf;
      // compute the form-factors
      if(inspin==2&&outspin==2)
	{
	  _f1.push_back(_Nfi[ix]*(1.-0.5*delta/mi
				  +0.25*delta/mi/mq*(1.-0.5*lambda/mf)*
				  (mi+mf-_eta[ix]*delta)
				  -0.125*delta/mi/mf/mQ*lambda*(mi+mf+_eta[ix]*delta)));
	  _f2.push_back(_Nfi[ix]*msum*(0.5/mi+0.25/mi/mq*(1.-0.5*lambda/mf)*
				       (delta-(mi+mf)*_eta[ix])
				       -0.125*lambda/mi/mf/mQ*(delta+(mi+mf)*_eta[ix])));
	  _f3.push_back(_Nfi[ix]*msum*(0.5/mi-0.25/mi/mq*(1.-0.5*lambda/mf)*
				       (mi+mf-_eta[ix]*delta)
				       +0.125*lambda/mi/mf/mQ*(mi+mf+_eta[ix]*delta)));
	  _g1.push_back(_Nfi[ix]*_eta[ix]*(1.+0.25*delta*lambda*(1./mi/mq-1./mf/mQ)));
	  _g2.push_back(-0.25*msum*_Nfi[ix]*_eta[ix]*lambda*(1./mi/mq-1./mf/mQ));
	  _g3.push_back(-0.25*msum*_Nfi[ix]*_eta[ix]*lambda*(1./mi/mq+1./mf/mQ));
	}
      else if(inspin==2&&outspin==4)
	{
	  _f1.push_back(2.*_Nfi[ix]/sqrt(3.)*(1.+0.5*lambda*(1./mq+1./mQ)));
	  _f2.push_back(_Nfi[ix]*msum/sqrt(3.)/mi*(1.+0.5*lambda*(1./mq+1./mQ)));
	  _f3.push_back(-_Nfi[ix]*msum*msum/mi/mf/sqrt(3.)*
			(1.+0.5*lambda*(1./mq+1./mQ)));
	  _g1.push_back(-2./sqrt(3.)*_Nfi[ix]);
	  _g2.push_back(-_Nfi[ix]*msum/sqrt(3.)*lambda/mq/mi);
	  _g3.push_back(-_f3.back());
	}
      else
	{throw InitException() << "Unknown spin combination in ChengHeavyBaryon" 
			       << "FormFactor::doinit() "  << Exception::abortnow;}
    }
  for(unsigned int ix=0;ix<numberOfFactors();++ix)
    {
      int id0,id1;
      particleID(ix,id0,id1);
      tcPDPtr part0=getParticleData(id0);Energy m0=part0->mass();
      tcPDPtr part1=getParticleData(id1);Energy m1=part1->mass();
      Complex f1v,f2v,f3v,f4v,f1a,f2a,f3a,f4a;
      if(part1->iSpin()==2)
	{SpinHalfSpinHalfFormFactor(0.*GeV2,ix,id0,id1,m0,m1,f1v,f2v,f3v,f1a,f2a,f3a);}
      else
	{SpinHalfSpinThreeHalfFormFactor(0.*GeV2,ix,id0,id1,m0,m1,f1v,f2v,f3v,
					 f4v,f1a,f2a,f3a,f4a);}
    }
}

ChengHeavyBaryonFormFactor::~ChengHeavyBaryonFormFactor() {}
  
void ChengHeavyBaryonFormFactor::persistentOutput(PersistentOStream & os) const {
  os << ounit(_mu,MeV) << ounit(_md,MeV) << ounit(_ms,MeV) << ounit(_mc,MeV) << ounit(_mb,MeV) 
     << _Nfi << _eta << _f1 << _f2 << _f3 
     << _g1 << _g2 << _g3 << ounit(_mVbc,MeV) << ounit(_mVbs,MeV) << ounit(_mVcs,MeV) 
     << ounit(_mVbd,MeV) << ounit(_mVcu,MeV) << ounit(_mAbc,MeV) 
     << ounit(_mAbs,MeV) << ounit(_mAcs,MeV) << ounit(_mAbd,MeV) << ounit(_mAcu,MeV);
 }
  

void ChengHeavyBaryonFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_mu,MeV) >> iunit(_md,MeV) >> iunit(_ms,MeV) >> iunit(_mc,MeV) >> iunit(_mb,MeV) 
     >> _Nfi >> _eta >> _f1 >> _f2 >> _f3 
     >> _g1 >> _g2 >> _g3 >> iunit(_mVbc,MeV) >> iunit(_mVbs,MeV) >> iunit(_mVcs,MeV) 
     >> iunit(_mVbd,MeV) >> iunit(_mVcu,MeV) >> iunit(_mAbc,MeV) 
     >> iunit(_mAbs,MeV) >> iunit(_mAcs,MeV) >> iunit(_mAbd,MeV) >> iunit(_mAcu,MeV);
}

ClassDescription<ChengHeavyBaryonFormFactor> ChengHeavyBaryonFormFactor::initChengHeavyBaryonFormFactor;
// Definition of the static class description member.

void ChengHeavyBaryonFormFactor::Init() {

  static ClassDocumentation<ChengHeavyBaryonFormFactor> documentation
    ("The ChengHeavyBaryonFormFactor class is the implementation"
     " of the form-factors of PRD53, 1457 and PRD56, 2799 for the weak decay of"
     "baryons containing a heavy quark. This model can be used for either"
     "semi-leptonic decays, or with the factorization approximation for"
     " non-leptonic weak decays");

  static Parameter<ChengHeavyBaryonFormFactor,Energy> interfaceUpMass
    ("DownMass",
     "The consituent mass of the down quark",
     &ChengHeavyBaryonFormFactor::_md, GeV, 0.322*GeV, 0.0*GeV, 10.0*GeV,
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

  static ParVector<ChengHeavyBaryonFormFactor,double> interfaceNfi
    ("Nfi",
     "The prefactor for a given form factor",
     &ChengHeavyBaryonFormFactor::_Nfi, -1, 1.0, -10.0, 10.0,
     false, false, true);

  static ParVector<ChengHeavyBaryonFormFactor,double> interfaceEta
    ("Eta",
     "The eta parameter for the form factor",
     &ChengHeavyBaryonFormFactor::_eta, -1, 0.0, -10.0, 10.0,
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
  Energy mV(0.*GeV),mA(0.*GeV);
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
  Energy mV(0.*GeV),mA(0.*GeV);
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
void ChengHeavyBaryonFormFactor::dataBaseOutput(ofstream& output,bool header,
						bool create) const
{
  if(header){output << "update decayers set parameters=\"";}
  if(create)
    {output << "create /Herwig++/ChengHeavyBaryonFormFactor " << fullName() << " \n";}
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
  for(unsigned int ix=0;ix<numberOfFactors();++ix)
    {
      if(ix<initialModes())
	{
	  output << "set " << fullName() << ":Nfi " << ix << "  " 
		<< _Nfi[ix] << endl;
	  output << "set " << fullName() << ":Eta " << ix << "  " 
		<< _eta[ix] << endl;
	}
      else
	{
	  output << "insert " << fullName() << ":Nfi " << ix << "  " 
		<< _Nfi[ix] << endl;
	  output << "insert " << fullName() << ":Eta " << ix << "  " 
		<< _eta[ix] << endl;
	}
    }
  BaryonFormFactor::dataBaseOutput(output,false,false);
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}

}
