// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BaryonWidthGenerator class.
//

#include "BaryonWidthGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Decay/Baryon/Baryon1MesonDecayerBase.h"

using namespace Herwig;
using namespace ThePEG;

void BaryonWidthGenerator::persistentOutput(PersistentOStream & os) const {
  os << _baryondecayers << _modeloc;
}

void BaryonWidthGenerator::persistentInput(PersistentIStream & is, int) {
  is >> _baryondecayers >> _modeloc;
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
}
 
void BaryonWidthGenerator::setupMode(tcDMPtr mode, tDecayIntegratorPtr decayer,
				     unsigned int) {
  // cast the decayer
  tBaryon1MesonDecayerBasePtr 
    baryon(dynamic_ptr_cast<tBaryon1MesonDecayerBasePtr>(decayer));
  if(baryon) {
    int dmode(baryon->findMode(*mode));
    if(dmode<0) {
      _baryondecayers.push_back(Baryon1MesonDecayerBasePtr());
      _modeloc.push_back(-1);
      return;
    }
    else {
      _baryondecayers.push_back(baryon);
      _modeloc.push_back(dmode);
    }
  }
  else {
    _baryondecayers.push_back(Baryon1MesonDecayerBasePtr());
    _modeloc.push_back(-1);
  }
}

void BaryonWidthGenerator::dataBaseOutput(ofstream & output, bool header) {
  if(header) output << "update Width_Generators set parameters=\"";
  // info from the base class
  GenericWidthGenerator::dataBaseOutput(output,false);
  // info from this class
  for(unsigned int ix=0;ix<_baryondecayers.size();++ix) {
    if(_baryondecayers[ix]) {
      output << "insert " << name() << ":BaryonDecayers " << ix 
	     << " " << _baryondecayers[ix]->fullName() << "\n";
    }
    else {
      output << "insert " << name() << ":BaryonDecayers " << ix 
	     << " NULL \n";
    }
    output << "insert " << name() << ":ModeLocation " << ix 
	   << " " << _modeloc[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" << name() << "\";" 
		    << endl;
}

Energy BaryonWidthGenerator::partial2BodyWidth(int imode, Energy q,Energy m1,
					       Energy m2) const {
  using Constants::pi;
  Complex A1,A2,A3,B1,B2,B3;
  if(q<m1+m2) 
    return ZERO;
  // mode from the base class
  int mecode(MEcode(imode));
  if(mecode<=100){ 
    return GenericWidthGenerator::partial2BodyWidth(imode,q,m1,m2);
  }
  // calcluate the decay momentum
  Energy2 q2(q*q),m12(m1*m1),m22(m2*m2),
    pcm2(0.25*(q2*(q2-2.*m12-2.*m22)+(m12-m22)*(m12-m22))/q2);
  Energy pcm(sqrt(pcm2)),gam(ZERO),msum(q+m1);
  // Energy m0(mass());
  Energy2 fact1((q+m1)*(q+m1)-m22),fact2((q-m1)*(q-m1)-m22),fact3(q2+m12-m22);
  // 1/2 -> 1/2 0
  if(mecode==101) {
    _baryondecayers[imode]->halfHalfScalarCoupling(_modeloc[imode],q,m1,m2,A1,B1);
    double Afact1((A1*conj(A1)).real()),Bfact1((B1*conj(B1)).real());
    gam = 0.125/pi/q2*pcm*(Afact1*fact1+Bfact1*fact2);
    /*
      Energy Qp(sqrt(pow(q+m1,2)-pow(m2,2))),Qm(sqrt(pow(q-m1,2)-pow(m2,2)));
      Complex h1(2.*Qp*A1),h2(-2.*Qm*B1);
      cout << "testing 1/2->1/2 0 " 
      << gam << "   " 
      << real(h1*conj(h1)+h2*conj(h2))/32./pi*pcm/q2     << "   " 
      << real(h1*conj(h1)+h2*conj(h2))/32./pi*pcm/q2/gam << endl;
    */
  }
  // 1/2 -> 1/2 1
  else if(mecode==102) {
    _baryondecayers[imode]->halfHalfVectorCoupling(_modeloc[imode],q,m1,m2,
						   A1,A2,B1,B2);
    double 
      Afact1((A1*conj(A1)).real()),Afact2((A2*conj(A2)).real()),
      Bfact1((B1*conj(B1)).real()),Bfact2((B2*conj(B2)).real()),
      Afact4((A1*conj(A2)+A2*conj(A1)).real()),
      Bfact4((B1*conj(B2)+B2*conj(B1)).real());
    Energy2 me(2.*(fact1*Bfact1+fact2*Afact1));
    if(m2>1e-10*GeV) {
      me+=1./m22*(fact1*Bfact1*(q-m1)*(q-m1)-2.*q*pcm*Bfact4*q*pcm*(q-m1)/msum
		  +fact2*q2*Bfact2*pcm2/msum/msum
		  +fact2*Afact1*msum  *msum  +2.*q*pcm*Afact4*q*pcm             
		  +fact1*q2*Afact2*pcm2/msum/msum);
    }
    gam = pcm/8./pi/q2*me;
    // test of the matrix element
    /*
      Energy2 Qp(sqrt(pow(q+m1,2)-pow(m2,2))),Qm(sqrt(pow(q-m1,2)-pow(m2,2)));
      double r2(sqrt(2.));
      Complex h1(2.*r2*Qp*B1),h2(-2.*r2*Qm*A1),h3(0.),h4(0.);
      if(m2>1e-10*GeV)
      {
      h3=2./m2*(Qp*(q-m1)*B1-Qm*q*B2*pcm/(q+m1));
      h4=2./m2*(Qm*(q+m1)*A1+Qp*q*A2*pcm/(q+m1));
      }
      cout << "testing 1/2->1/2 0 " 
      << gam << "   " 
      << real(h1*conj(h1)+h2*conj(h2)+h3*conj(h3)+h4*conj(h4))/32./pi*pcm/q2     
      << "   " 
      << real(h1*conj(h1)+h2*conj(h2)+h3*conj(h3)+h4*conj(h4))/32./pi*pcm/q2/gam 
      << endl;
    */
  }
  // 1/2 -> 3/2 0
  else if(mecode==103) {
    _baryondecayers[imode]->halfThreeHalfScalarCoupling(_modeloc[imode],q,m1,m2,A1,B1);
    double Afact1((A1*conj(A1)).real()),Bfact1((B1*conj(B1)).real());
    gam = 0.25/pi/3./msum/msum/m12*pcm*pcm2*(Afact1*fact1+Bfact1*fact2);
    /*
      Energy2 Qp(sqrt(pow(q+m1,2)-pow(m2,2))),Qm(sqrt(pow(q-m1,2)-pow(m2,2)));
      double r23(sqrt(2./3.));
      Complex h1(-2.*r23*pcm*q/m1*Qm*B1/(q+m1)),h2( 2.*r23*pcm*q/m1*Qp*A1/(q+m1));
      cout << "testing 1/2->3/2 0 "
      << gam << "   " 
      << real(h1*conj(h1)+h2*conj(h2))/32./pi*pcm/q2     << "   " 
      << real(h1*conj(h1)+h2*conj(h2))/32./pi*pcm/q2/gam << endl;
    */
  }
  // 1/2 -> 3/2 1
  else if(mecode==104) {
    Energy Qp(sqrt(fact1)),Qm(sqrt(fact2));
    double r2(sqrt(2.)),r3(sqrt(3.));
    complex<Energy> h1(-2.*Qp*A1),h2(2.*Qm*B1),h5(ZERO),h6(ZERO);
    complex<Energy> h3(-2./r3*Qp*(A1-fact2/m1*A2/msum));
    complex<Energy> h4( 2./r3*Qm*(B1-fact1/m1*B2/msum));
    if(m2>1e-10*GeV) {
      h5=-2.*r2/r3/m1/m2*Qp*(0.5*(q2-m12-m22)*A1+0.5*fact2*(q+m1)*A2/msum
			     +q2*pcm*pcm*A3/sqr(msum));
      h6= 2.*r2/r3/m1/m2*Qm*(0.5*(q2-m12-m22)*B1-0.5*fact1*(q-m1)*B2/msum
			     +q2*pcm*pcm*B3/sqr(msum));
    }
    gam=real(+h1*conj(h1)+h2*conj(h2)+h3*conj(h3)
	     +h4*conj(h4)+h5*conj(h5)+h6*conj(h6))/32./pi*pcm/q2;
  }
  // 3/2 -> 1/2 0
  else if(mecode==105) {
    _baryondecayers[imode]->threeHalfHalfScalarCoupling(_modeloc[imode],q,m1,m2,A1,B1);
    double Afact1((A1*conj(A1)).real()),Bfact1((B1*conj(B1)).real());
    gam = 0.125/3./pi/msum/msum/q2*pcm*pcm2*(Afact1*fact1+Bfact1*fact2);
    /*
      Energy2 Qp(sqrt(pow(q+m1,2)-pow(m2,2))),Qm(sqrt(pow(q-m1,2)-pow(m2,2)));
      double r23(sqrt(2./3.));
      Complex h1(-2.*r23*pcm*Qm*B1/(q+m1)),
      h2( 2.*r23*pcm*Qp*A1/(q+m1));
      cout << "testing 3/2->1/2 0 "
      << gam << "   " 
      << real(h1*conj(h1)+h2*conj(h2))/64./pi*pcm/q2     << "   " 
      << real(h1*conj(h1)+h2*conj(h2))/64./pi*pcm/q2/gam << endl;
    */
  }
  // 3/2 -> 1/2 1
  else if(mecode==106) {
    _baryondecayers[imode]->threeHalfHalfVectorCoupling(_modeloc[imode],q,m1,m2,
							A1,A2,A3,B1,B2,B3);
    Energy Qp(sqrt(fact1)),Qm(sqrt(fact2));
    double r2(sqrt(2.)),r3(sqrt(3.));
    complex<Energy> h1(-2.*Qp*A1),h2(2.*Qm*B1),h5(ZERO),h6(ZERO);
    complex<Energy> h3(-2./r3*Qp*(A1-fact2/q*(A2/msum)));
    complex<Energy> h4( 2./r3*Qm*(B1-fact1/q*(B2/msum)));
    if(m2>1e-10*GeV) {
      h5=-2.*r2/r3/q/m2*Qp*(0.5*(m12-q2-m22)*A1+0.5*fact2*(m1+q)*A2/msum
			    +q2*pcm*pcm*(A3/sqr(msum)));
      h6= 2.*r2/r3/q/m2*Qm*(0.5*(m12-q2-m22)*B1-0.5*fact1*(m1-q)*(B2/msum)
			    +q2*pcm*pcm*(B3/sqr(msum)));
    }
    gam=real(+h1*conj(h1)+h2*conj(h2)+h3*conj(h3)
	     +h4*conj(h4)+h5*conj(h5)+h6*conj(h6))/64./pi*pcm/q2;
  }
  // 3/2 -> 3/2 0
  else if(mecode==107) {
    _baryondecayers[imode]->threeHalfThreeHalfScalarCoupling(_modeloc[imode],q,m1,m2,
							     A1,A2,B1,B2);
    complex<InvEnergy2> A2byM2 = A2/sqr(msum);
    complex<InvEnergy2> B2byM2 = B2/sqr(msum);
    double Afact1((A1*conj(A1)).real());
    InvEnergy4 Afact2((A2byM2*conj(A2byM2)).real());
    double Bfact1((B1*conj(B1)).real());
    InvEnergy4 Bfact2((B2byM2*conj(B2byM2)).real());
    InvEnergy2 Afact4((A1*conj(A2byM2)+A2byM2*conj(A1)).real());
    InvEnergy2 Bfact4((B1*conj(B2byM2)+B2byM2*conj(B1)).real());
    gam = pcm/36./pi/q2/q2/m12*
      ( fact1*(Afact2*q2*q2*pcm2*pcm2
	       +0.25*Afact1*(fact3*fact1+10.*q2*m12)
	       +0.5*Afact4*q2*pcm*pcm*(fact3+q*m1))+
	fact2*(Bfact2*q2*q2*pcm2*pcm2
	       +0.25*Bfact1*(fact3*fact2+10.*q2*m12)
	       +0.5*Bfact4*q2*pcm*pcm*(fact3-q*m1)));
  }
  else {
    throw Exception() << "Unknown type of mode " << mecode 
		      << " in BaryonWidthGenerator::partial2BodyWidth() " 
		      << Exception::abortnow;
  }
  return gam*MEcoupling(imode)*MEcoupling(imode);
}

void BaryonWidthGenerator::doinit() {
  if(initialize()) { 
    _baryondecayers.clear();
    _modeloc.clear();
  }
  GenericWidthGenerator::doinit();
}
