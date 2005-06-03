// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BaryonWidthGenerator class.
//

#include "BaryonWidthGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/RefVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "BaryonWidthGenerator.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Decay/Baryon/Baryon1MesonDecayerBase.h"

namespace Herwig {
using namespace ThePEG;


BaryonWidthGenerator::~BaryonWidthGenerator() {}

void BaryonWidthGenerator::persistentOutput(PersistentOStream & os) const {
  os << _baryondecayers << _modeloc
     << _Afact1 << _Afact2 << _Afact3 << _Afact4 << _Afact5 << _Afact6
     << _Bfact1 << _Bfact2 << _Bfact3 << _Bfact4 << _Bfact5 << _Bfact6;
}

void BaryonWidthGenerator::persistentInput(PersistentIStream & is, int) {
  is >> _baryondecayers >> _modeloc
     >> _Afact1 >> _Afact2 >> _Afact3 >> _Afact4 >> _Afact5 >> _Afact6
     >> _Bfact1 >> _Bfact2 >> _Bfact3 >> _Bfact4 >> _Bfact5 >> _Bfact6;
}

ClassDescription<BaryonWidthGenerator> BaryonWidthGenerator::initBaryonWidthGenerator;
// Definition of the static class description member.

void BaryonWidthGenerator::Init() {

  static ClassDocumentation<BaryonWidthGenerator> documentation
    ("The BaryonWidthGenerator class is designed for the calculation of the"
     " running width for baryons.");

  static RefVector<BaryonWidthGenerator,Baryon1MesonDecayerBase> interfaceBaryonDecayers
    ("BaryonDecayers",
     "Pointers to the baryon decayers to get the couplings",
     &BaryonWidthGenerator::_baryondecayers, -1, false, false, true, true, false);

  static ParVector<BaryonWidthGenerator,int> interfaceModeLocation
    ("ModeLocation",
     "The location of the mode in the decayer",
     &BaryonWidthGenerator::_modeloc, 0, -1, 0, 0,
     false, false, false);

  static ParVector<BaryonWidthGenerator,double> interfaceA1factor
    ("A1factor",
     "The first factor for the partial width for the non-gamma_5 terms.",
     &BaryonWidthGenerator::_Afact1, -1, 0.0, 0, 0,
     false, false, false);

  static ParVector<BaryonWidthGenerator,double> interfaceA2factor
    ("A2factor",
     "The second factor for the partial width for the non-gamma_5 terms.",
     &BaryonWidthGenerator::_Afact2, -1, 0.0, 0, 0,
     false, false, false);

  static ParVector<BaryonWidthGenerator,double> interfaceA3factor
    ("A3factor",
     "The third factor for the partial width for the non-gamma_5 terms.",
     &BaryonWidthGenerator::_Afact3, -1, 0.0, 0, 0,
     false, false, false);

  static ParVector<BaryonWidthGenerator,double> interfaceA4factor
    ("A4factor",
     "The fourth factor for the partial width for the non-gamma_5 terms.",
     &BaryonWidthGenerator::_Afact4, -1, 0.0, 0, 0,
     false, false, false);

  static ParVector<BaryonWidthGenerator,double> interfaceA5factor
    ("A5factor",
     "The fifth factor for the partial width for the non-gamma_5 terms.",
     &BaryonWidthGenerator::_Afact5, -1, 0.0, 0, 0,
     false, false, false);

  static ParVector<BaryonWidthGenerator,double> interfaceA6factor
    ("A6factor",
     "The sixth factor for the partial width for the non-gamma_5 terms.",
     &BaryonWidthGenerator::_Afact6, -1, 0.0, 0, 0,
     false, false, false);

  static ParVector<BaryonWidthGenerator,double> interfaceB1factor
    ("B1factor",
     "The first factor for the partial width for the gamma_5 terms.",
     &BaryonWidthGenerator::_Bfact1, -1, 0.0, 0, 0,
     false, false, false);

  static ParVector<BaryonWidthGenerator,double> interfaceB2factor
    ("B2factor",
     "The second factor for the partial width for the gamma_5 terms.",
     &BaryonWidthGenerator::_Bfact2, -1, 0.0, 0, 0,
     false, false, false);

  static ParVector<BaryonWidthGenerator,double> interfaceB3factor
    ("B3factor",
     "The third factor for the partial width for the gamma_5 terms.",
     &BaryonWidthGenerator::_Bfact3, -1, 0.0, 0, 0,
     false, false, false);

  static ParVector<BaryonWidthGenerator,double> interfaceB4factor
    ("B4factor",
     "The fourth factor for the partial width for the gamma_5 terms.",
     &BaryonWidthGenerator::_Bfact4, -1, 0.0, 0, 0,
     false, false, false);

  static ParVector<BaryonWidthGenerator,double> interfaceB5factor
    ("B5factor",
     "The fifth factor for the partial width for the gamma_5 terms.",
     &BaryonWidthGenerator::_Bfact5, -1, 0.0, 0, 0,
     false, false, false);

  static ParVector<BaryonWidthGenerator,double> interfaceB6factor
    ("B6factor",
     "The sixth factor for the partial width for the gamma_5 terms.",
     &BaryonWidthGenerator::_Bfact6, -1, 0.0, 0, 0,
     false, false, false);
}
 
void BaryonWidthGenerator::setupMode(tcDMPtr mode, tDecayIntegratorPtr decayer,
				     unsigned int imode)
{
  // cast the decayer
  tBaryon1MesonDecayerBasePtr 
    baryon(dynamic_ptr_cast<tBaryon1MesonDecayerBasePtr>(decayer));
  if(baryon)
    {
      ParticleMSet::const_iterator pit = mode->products().begin();
      tcPDPtr part0(mode->parent());
      tcPDPtr part1(*pit);++pit;
      tcPDPtr part2(*pit);
      Energy m0(part0->mass()),m1(part1->mass()),m2(part2->mass());
      int dmode(baryon->findMode(*mode));
      Complex A1(0.),A2(0.),A3(0.),B1(0.),B2(0.),B3(0.);
      if(dmode<0)
	{
	  _Afact1.push_back(0.);_Afact3.push_back(0.);_Afact5.push_back(0.);
	  _Afact2.push_back(0.);_Afact4.push_back(0.);_Afact6.push_back(0.);
	  _Bfact1.push_back(0.);_Bfact3.push_back(0.);_Bfact5.push_back(0.);
	  _Bfact2.push_back(0.);_Bfact4.push_back(0.);_Bfact6.push_back(0.);
	  _baryondecayers.push_back(Baryon1MesonDecayerBasePtr());
	  _modeloc.push_back(-1);
	  return;
	}
      // 1/2 -> 1/2 0
      if(part0->iSpin()==2&&((part1->iSpin()==2&&part2->iSpin()==1)||
			     (part1->iSpin()==1&&part2->iSpin()==2)))
	{baryon->halfHalfScalarCoupling(dmode,m0,m1,m2,A1,B1);}
      // 1/2 -> 1/2 1
      else if(part0->iSpin()==2&&((part1->iSpin()==2&&part2->iSpin()==3)||
				  (part1->iSpin()==3&&part2->iSpin()==2)))
	{baryon->halfHalfVectorCoupling(dmode,m0,m1,m2,A1,A2,B1,B2);}
      // 1/2 -> 3/2 0
      else if(part0->iSpin()==2&&((part1->iSpin()==4&&part2->iSpin()==1)||
				  (part1->iSpin()==1&&part2->iSpin()==4)))
	{baryon->halfThreeHalfScalarCoupling(dmode,m0,m1,m2,A1,B1);}
      // 1/2 -> 3/2 1
      else if(part0->iSpin()==2&&((part1->iSpin()==4&&part2->iSpin()==3)||
				  (part1->iSpin()==3&&part2->iSpin()==4)))
	{baryon->halfThreeHalfVectorCoupling(dmode,m0,m1,m2,A1,A2,A3,B1,B2,B3);}
      // 3/2 -> 1/2 0
      else if(part0->iSpin()==4&&((part1->iSpin()==2&&part2->iSpin()==1)||
				  (part1->iSpin()==1&&part2->iSpin()==2)))
	{baryon->threeHalfHalfScalarCoupling(dmode,m0,m1,m2,A1,B1);}
      // 3/2 -> 1/2 1
      else if(part0->iSpin()==4&&((part1->iSpin()==2&&part2->iSpin()==3)||
				  (part1->iSpin()==3&&part2->iSpin()==2)))
	{baryon->threeHalfHalfVectorCoupling(dmode,m0,m1,m2,A1,A2,A3,B1,B2,B3);}
      // 3/2 -> 3/2 0
      else if(part0->iSpin()==4&&((part1->iSpin()==4&&part2->iSpin()==1)||
				  (part1->iSpin()==1&&part2->iSpin()==4)))
	{baryon->threeHalfThreeHalfScalarCoupling(dmode,m0,m1,m2,A1,A2,B1,B2);}
      else
	{cout << "unimplemented spin " 
	      << part0->iSpin() << " " << part1->iSpin() << " " 
	      << part2->iSpin() << endl;}
      _baryondecayers.push_back(baryon);
      _modeloc.push_back(dmode);
      _Afact1.push_back((A1*conj(A1)).real());
      _Afact2.push_back((A2*conj(A2)).real());
      _Afact3.push_back((A3*conj(A3)).real());
      _Afact4.push_back((A1*conj(A2)+conj(A1)*A2).real());
      _Afact5.push_back((A1*conj(A3)+conj(A1)*A3).real());
      _Afact6.push_back((A2*conj(A3)+conj(A2)*A3).real());
      _Bfact1.push_back((B1*conj(B1)).real());
      _Bfact2.push_back((B2*conj(B2)).real());
      _Bfact3.push_back((B3*conj(B3)).real());
      _Bfact4.push_back((B1*conj(B2)+conj(B1)*B2).real());
      _Bfact5.push_back((B1*conj(B3)+conj(B1)*B3).real());
      _Bfact6.push_back((B2*conj(B3)+conj(B2)*B3).real());
    }
  else
    {
      _baryondecayers.push_back(Baryon1MesonDecayerBasePtr());
      _modeloc.push_back(-1);
      _Afact1.push_back(0.);_Afact3.push_back(0.);_Afact5.push_back(0.);
      _Afact2.push_back(0.);_Afact4.push_back(0.);_Afact6.push_back(0.);
      _Bfact1.push_back(0.);_Bfact3.push_back(0.);_Bfact5.push_back(0.);
      _Bfact2.push_back(0.);_Bfact4.push_back(0.);_Bfact6.push_back(0.);
    }
}

void BaryonWidthGenerator::dataBaseOutput(ofstream & output, bool header)
{
  cout << "testing " << this << "  " << fullName() << endl;
  if(header){output << "update Width_Generators set parameters=\"";}
  // info from the base class
  GenericWidthGenerator::dataBaseOutput(output,false);
  // info from this class
  cout << "testing the size " << _Afact1.size() << " " << _baryondecayers.size() << endl;
  for(unsigned int ix=0;ix<_Afact1.size();++ix)
    {
      if(_baryondecayers[ix])
	{
	  cout << "testing A " << ix << " " << _baryondecayers[ix] << endl;
	  output << "insert " << fullName() << ":BaryonDecayers " << ix 
		<< " " << _baryondecayers[ix]->fullName() << "\n";
	}
      else
	{output << "insert " << fullName() << ":BaryonDecayers " << ix 
		<< " NULL \n";}
      output << "insert " << fullName() << ":ModeLocation " << ix 
	     << " " << _modeloc[ix] << "\n";
      output << "insert " << fullName() << ":A1factor " << ix 
	     << " " << _Afact1[ix] << "\n";
      output << "insert " << fullName() << ":A2factor " << ix 
	     << " " << _Afact2[ix] << "\n";
      output << "insert " << fullName() << ":A3factor " << ix 
	     << " " << _Afact3[ix] << "\n";
      output << "insert " << fullName() << ":A4factor " << ix 
	     << " " << _Afact4[ix] << "\n";
      output << "insert " << fullName() << ":A5factor " << ix 
	     << " " << _Afact5[ix] << "\n";
      output << "insert " << fullName() << ":A6factor " << ix 
	     << " " << _Afact6[ix] << "\n";
      output << "insert " << fullName() << ":B1factor " << ix 
	     << " " << _Bfact1[ix] << "\n";
      output << "insert " << fullName() << ":B2factor " << ix 
	     << " " << _Bfact2[ix] << "\n";
      output << "insert " << fullName() << ":B3factor " << ix 
	     << " " << _Bfact3[ix] << "\n";
      output << "insert " << fullName() << ":B4factor " << ix 
	     << " " << _Bfact4[ix] << "\n";
      output << "insert " << fullName() << ":B5factor " << ix 
	     << " " << _Bfact5[ix] << "\n";
      output << "insert " << fullName() << ":B6factor " << ix 
	     << " " << _Bfact6[ix] << "\n";
    }
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}

Energy BaryonWidthGenerator::partial2BodyWidth(int imode, Energy q,Energy m1,
					       Energy m2) const
{
  Complex A1,A2,A3,B1,B2,B3;
  if(q<m1+m2){return 0.;}
  // mode from the base class
  int mecode(MEcode(imode));
  if(mecode<=100){return GenericWidthGenerator::partial2BodyWidth(imode,q,m1,m2);}
  // calcluate the decay momentum
  Energy2 q2(q*q),m12(m1*m1),m22(m2*m2),
    pcm2(0.25*(q2*(q2-2.*m12-2.*m22)+(m12-m22)*(m12-m22))/q2);
  Energy pcm(sqrt(pcm2)),gam(0.),msum(q+m1);
  // Energy m0(mass());
  Energy2 fact1((q+m1)*(q+m1)-m22),fact2((q-m1)*(q-m1)-m22),fact3(q2+m12-m22);
  // 1/2 -> 1/2 0
  if(mecode==101)
    {
      _baryondecayers[imode]->halfHalfScalarCoupling(_modeloc[imode],q,m1,m2,A1,B1);
      double Afact1((A1*conj(A1)).real()),Bfact1((B1*conj(B1)).real());
      gam = 0.125/pi/q2*pcm*(Afact1*fact1+Bfact1*fact2);
    }
  // 1/2 -> 1/2 1
  else if(mecode==102)
    {
      _baryondecayers[imode]->halfHalfVectorCoupling(_modeloc[imode],q,m1,m2,
						     A1,A2,B1,B2);
      double 
	Afact1((A1*conj(A1)).real()),Afact2((A2*conj(A2)).real()),
	Bfact1((B1*conj(B1)).real()),Bfact2((B2*conj(B2)).real()),
	Afact4((A1*conj(A2)+A2*conj(A1)).real()),
	Bfact4((B1*conj(B2)+B2*conj(B1)).real());
      double me(2.*(fact1*Bfact1+fact2*Afact1));
      gam = pcm/8./pi/q2*me;
    }
  // 1/2 -> 3/2 0
  else if(mecode==103)
    {
      _baryondecayers[imode]->halfThreeHalfScalarCoupling(_modeloc[imode],q,m1,m2,A1,B1);
      double Afact1((A1*conj(A1)).real()),Bfact1((B1*conj(B1)).real());
      gam = 0.25/pi/3./msum/msum/m12*pcm*pcm2*(Afact1*fact1+Bfact1*fact2);
    }
  // 3/2 -> 1/2 0
  else if(mecode==105)
    {
      _baryondecayers[imode]->threeHalfHalfScalarCoupling(_modeloc[imode],q,m1,m2,A1,B1);
      double Afact1((A1*conj(A1)).real()),Bfact1((B1*conj(B1)).real());
      gam = 0.125/3./pi/msum/msum/q2*pcm*pcm2*(Afact1*fact1+Bfact1*fact2);
    }
  // 3/2 -> 1/2 1
  else if(mecode==106)
    {
      _baryondecayers[imode]->threeHalfHalfVectorCoupling(_modeloc[imode],q,m1,m2,
							  A1,A2,A3,B1,B2,B3);
      A2 /=msum;A3 /=(msum*msum);
      B2 /=msum;B3 /=(msum*msum);
      double 
	Afact1((A1*conj(A1)).real()),Afact2((A2*conj(A2)).real()),
	Bfact1((B1*conj(B1)).real()),Bfact2((B2*conj(B2)).real()),
	Afact4((A1*conj(A2)+A2*conj(A1)).real()),
	Bfact4((B1*conj(B2)+B2*conj(B1)).real());
      double me((fact1*(4.*Afact1-fact2/q*Afact4+fact2*fact2/q2*Afact2)+
		 fact2*(4.*Bfact1+fact1/q*Bfact4+fact1*fact1/q2*Bfact2))/6.);
      gam = pcm/8./pi/q2*me;
    }
  // 3/2 -> 3/2 0
  else if(mecode==107)
    {
      _baryondecayers[imode]->threeHalfThreeHalfScalarCoupling(_modeloc[imode],q,m1,m2,
							       A1,A2,B1,B2);
      A2/=(msum*msum);
      B2/=(msum*msum);
      double 
	Afact1((A1*conj(A1)).real()),Afact2((A2*conj(A2)).real()),
	Bfact1((B1*conj(B1)).real()),Bfact2((B2*conj(B2)).real()),
	Afact4((A1*conj(A2)+A2*conj(A1)).real()),
	Bfact4((B1*conj(B2)+B2*conj(B1)).real());
      gam = pcm/36./pi/q2/q2/m12*
	( fact1*(Afact2*q2*q2*pcm2*pcm2
		 +0.25*Afact1*(fact3*fact1+10.*q2*m12)
		 +0.5*Afact4*q2*pcm*pcm*(fact3+q*m1))+
	  fact2*(Bfact2*q2*q2*pcm2*pcm2
		 +0.25*Bfact1*(fact3*fact2+10.*q2*m12)
		 +0.5*Bfact4*q2*pcm*pcm*(fact3-q*m1)));
      /*
      cout << "testing piecesA " <<  pcm/36./pi/q2/q2/m12 << " " 
	   << fact1 << " " << Afact2*q2*q2*pcm2*pcm2 << " " 
	   << +0.25*Afact1 << " " << fact3*fact1 << " " << 10.*q2*m12 << " " 
	   << +0.5*Afact4*q2*pcm*pcm*(fact3+q*m1) << " " << fact2 << " " 
	   << Bfact2*q2*q2*pcm2*pcm2 << " " 
	   << +0.25*Bfact1*(fact3*fact2+10.*q2*m12) << " "
	   << +0.5*Bfact4*q2*pcm*pcm*(fact3-q*m1) << endl;
      */
    }
  else
    {throw Exception() << "Unknown type of mode " << mecode 
		       << " in BaryonWidthGenerator::partial2BodyWidth() " 
		       << Exception::abortnow;}
  return gam*MEcoupling(imode)*MEcoupling(imode);
}

}
