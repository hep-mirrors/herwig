// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the KMatrix class.
//

#include "KMatrix.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/lu.hpp>

using namespace Herwig;

KMatrix::KMatrix(FlavourInfo flavour, vector<Channels> channels,
		 vector<Energy2> poles, vector<vector<Energy> > g)
  : flavour_(flavour), channels_(channels), poles_(poles), g_(g)
{}

void KMatrix::persistentOutput(PersistentOStream & os) const {
  os << ounit(poles_,GeV2) << ounit(g_,GeV)
     << ounit(mPiPlus_,GeV) << ounit(mPi0_,GeV)
     << ounit(mKPlus_,GeV) << ounit(mK0_,GeV) << ounit(mEta_,GeV)
     << ounit(mEtaPrime_,GeV);
}

void KMatrix::persistentInput(PersistentIStream & is, int) {
  is >> iunit(poles_,GeV2) >> iunit(g_,GeV)
     >> iunit(mPiPlus_,GeV) >> iunit(mPi0_,GeV)
     >> iunit(mKPlus_,GeV) >> iunit(mK0_,GeV) >> iunit(mEta_,GeV)
     >> iunit(mEtaPrime_,GeV);
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<KMatrix,Interfaced>
describeHerwigKMatrix("Herwig::KMatrix", "Herwig.so");

void KMatrix::Init() {

  static ClassDocumentation<KMatrix> documentation
    ("The KMatrix class provides a base class for the implementation of "
     "K-matrix parameterizations in Herwig");

}

void KMatrix::doinit() {
  Interfaced::doinit();
  // The charged pion mass
  mPiPlus_=getParticleData(ParticleID::piplus)->mass();
  // The neutral pion mass
  mPi0_=getParticleData(ParticleID::pi0)->mass();
  // The charged kaon mass
  mKPlus_=getParticleData(ParticleID::Kplus)->mass();
  // The neutral kaon mass
  mK0_=getParticleData(ParticleID::K0)->mass();
  // The eta mass
  mEta_=getParticleData(ParticleID::eta)->mass();
  // The eta' mass
  mEtaPrime_=getParticleData(ParticleID::etaprime)->mass();
}

namespace {

  double kallen(const Energy2 &s, const Energy &m1, const Energy & m2) {
    return (1.-sqr(m1+m2)/s)*(1.-sqr(m1-m2)/s);
  }
}

ublas::matrix<Complex> KMatrix::rho(Energy2 s) const {
  size_t msize = channels_.size();
  ublas::diagonal_matrix<Complex> rho(msize,msize);
  for(unsigned int iChan=0;iChan<msize;++iChan) {
    double val(0);
    switch (channels_[iChan]) {
    case PiPi:
      val=kallen(s,mPiPlus_,mPiPlus_);
      break;
    case KPi:
      val=kallen(s,mKPlus_,mPiPlus_);
      break;
    case KEta:
      val=kallen(s,mKPlus_,mEta_);
      break;
    case KEtaPrime:
      val=kallen(s,mKPlus_,mEtaPrime_);
      break;
    default:
      assert(false);
    }
    if(val>=0)
      rho(iChan,iChan) = sqrt(val);
    else
      rho(iChan,iChan) = Complex(0.,1.)*sqrt(-val);
  }
  return rho;
}

ublas::vector<Complex> KMatrix::
amplitudes(Energy2 s, ublas::vector<Complex> pVector, bool multiplyByPoles) const {
  static const Complex ii(0.,1.);
  const ublas::identity_matrix<Complex> I(channels_.size());
  double fact(1.);
  if(multiplyByPoles) {
    for (Energy2 pole : poles_ ) fact *= 1.-s/pole;
  }
  // matrix for which we need the inverse
  ublas::matrix<Complex> m = fact*I-ii*prod(K(s,multiplyByPoles),rho(s));
  // inverse matrix
  ublas::matrix<Complex> inverse = ublas::identity_matrix<Complex>(m.size1());
  // 1x1 just a number
  if(m.size1()==1) {
    inverse(0,0) = 1./m(0,0);
  }
  // compute the inverse
  else {
    // create a permutation matrix for the LU-factorization
    ublas::permutation_matrix<std::size_t>  pm(m.size1());
    // perform LU-factorization
    int res = ublas::lu_factorize(m,pm);
    if( res != 0 ) {
      cerr << "problem with factorization\n";
      exit(1);
    }
    // back substitute to get the inverse
    ublas::lu_substitute(m, pm, inverse);
  }
  // compute the amplitudes
  return prod(inverse,pVector);
}
