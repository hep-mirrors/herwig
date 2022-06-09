// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the KPiIHalfFOCUSKMatrix class.
//

#include "KPiIHalfFOCUSKMatrix.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

KPiIHalfFOCUSKMatrix::KPiIHalfFOCUSKMatrix()
  : KMatrix(FlavourInfo(IsoSpin::IHalf,IsoSpin::I3Unknown,
			Strangeness::PlusOne,Charm::Zero,
			Beauty::Zero),
	    vector<Channels>({KMatrix::KPi,KMatrix::KEtaPrime}),
	    vector<Energy2>({1.7919*GeV2}),
	    vector<vector<Energy> >(1,vector<Energy>({0.31072*GeV,-0.02323*GeV}))),
    C11_({0.79299,-0.15099,0.00811}),
    C22_({0.15040,-0.038266,0.0022596}),
    C12_({0.17054,-0.0219,0.00085655}),
	    sHalf_(0.23*GeV2)
{}

IBPtr KPiIHalfFOCUSKMatrix::clone() const {
  return new_ptr(*this);
}

IBPtr KPiIHalfFOCUSKMatrix::fullclone() const {
  return new_ptr(*this);
}

void KPiIHalfFOCUSKMatrix::persistentOutput(PersistentOStream & os) const {
  os << C11_ << C22_ << C12_ << ounit(sHalf_,GeV2) << ounit(sNorm_,GeV2);
}

void KPiIHalfFOCUSKMatrix::persistentInput(PersistentIStream & is, int) {
  is >> C11_ >> C22_ >> C12_ >> iunit(sHalf_,GeV2) >> iunit(sNorm_,GeV2);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<KPiIHalfFOCUSKMatrix,KMatrix>
describeHerwigKPiIHalfFOCUSKMatrix("Herwig::KPiIHalfFOCUSKMatrix", "HwFormFactors.so");

void KPiIHalfFOCUSKMatrix::Init() {

  static ClassDocumentation<KPiIHalfFOCUSKMatrix> documentation
    ("The KPiIHalfFOCUSKMatrix class implements the K-matrix fit of "
     "the FOCUS collaboration (Phys.Lett. B653 (2007) 1-11) for the I=1/2 "
     "component of the Kpi K-matrix.",
     "The KPiIHalfFOCUSKMatrix class implements the K-matrix fit of "
     "the FOCUS collaboration \\cite{Pennington:2007se} for the $I=1/2$ "
     "component of the $K\\pi$ K-matrix.",
     "\\bibitem{Pennington:2007se}"
     "J.~M.~Link {\\it et al.} [FOCUS Collaboration],"
     "%``Dalitz plot analysis of the $D^{+} \\to K^{-} \\pi^{+} \\pi^{+}$ decay in the FOCUS experiment,''"
     "Phys.\\ Lett.\\ B {\\bf 653} (2007) 1"
     "doi:10.1016/j.physletb.2007.06.070"
     "[arXiv:0705.2248 [hep-ex]]."
     "%%CITATION = doi:10.1016/j.physletb.2007.06.070;%%"
     "%79 citations counted in INSPIRE as of 14 Jan 2020");

}

void KPiIHalfFOCUSKMatrix::doinit() {
  KMatrix::doinit();
  Energy mK = getParticleData(ParticleID::Kplus)->mass();
  Energy mpi= getParticleData(ParticleID::piplus)->mass();
  sNorm_ = sqr(mK)+sqr(mpi);
}

boost::numeric::ublas::matrix<double> KPiIHalfFOCUSKMatrix::K(Energy2 s, bool multiplyByPoles) const {
  double st = s/sNorm_-1.;
  double pre = (s-sHalf_)/sNorm_;
  Energy2 denom = !multiplyByPoles ? poles()[0]-s : poles()[0];
  boost::numeric::ublas::matrix<double> output =
    boost::numeric::ublas::zero_matrix<double>(2,2);
  output(0,0) = poleCouplings()[0][0]*poleCouplings()[0][0]/denom;
  output(0,1) = poleCouplings()[0][0]*poleCouplings()[0][1]/denom;
  output(1,1) = poleCouplings()[0][1]*poleCouplings()[0][1]/denom;
  double param = !multiplyByPoles ? 1. : (1.-s/poles()[0]);
  for(unsigned int ix=0;ix<C11_.size();++ix) {
    output(0,0) += C11_[ix]*param;
    output(1,1) += C22_[ix]*param;
    output(0,1) += C12_[ix]*param;
    param *= st;
  }
  output(1,0) = output(0,1);
  output *= pre;
  return output;
}
