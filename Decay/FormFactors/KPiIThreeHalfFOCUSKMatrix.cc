// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the KPiIThreeHalfFOCUSKMatrix class.
//

#include "KPiIThreeHalfFOCUSKMatrix.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

KPiIThreeHalfFOCUSKMatrix::KPiIThreeHalfFOCUSKMatrix()
  : KMatrix(FlavourInfo(IsoSpin::IThreeHalf,IsoSpin::I3Unknown,
			Strangeness::PlusOne,Charm::Zero,
			Beauty::Zero),
	    vector<Channels>({KMatrix::KPi})),
    D_({-0.22147,0.026637,-0.00092057}),
    sThreeHalf_(0.27*GeV2)
{}

IBPtr KPiIThreeHalfFOCUSKMatrix::clone() const {
  return new_ptr(*this);
}

IBPtr KPiIThreeHalfFOCUSKMatrix::fullclone() const {
  return new_ptr(*this);
}

void KPiIThreeHalfFOCUSKMatrix::persistentOutput(PersistentOStream & os) const {
  os << D_ << ounit(sThreeHalf_,GeV2) << ounit(sNorm_,GeV2);
}

void KPiIThreeHalfFOCUSKMatrix::persistentInput(PersistentIStream & is, int) {
  is >> D_ >> iunit(sThreeHalf_,GeV2) >> iunit(sNorm_,GeV2);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<KPiIThreeHalfFOCUSKMatrix,KMatrix>
describeHerwigKPiIThreeHalfFOCUSKMatrix("Herwig::KPiIThreeHalfFOCUSKMatrix", "HwFormFactors.so");

void KPiIThreeHalfFOCUSKMatrix::Init() {

  static ClassDocumentation<KPiIThreeHalfFOCUSKMatrix> documentation
    ("The KPiIThreeHalfFOCUSKMatrix class implements the K-matrix fit of "
     "the FOCUS collaboration (Phys.Lett. B653 (2007) 1-11) for the I=3/2 "
     "component of the Kpi K-matrix.",
     "The KPiIThreeHalfFOCUSKMatrix class implements the K-matrix fit of "
     "the FOCUS collaboration \\cite{Pennington:2007se} for the $I=3/2$ "
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

void KPiIThreeHalfFOCUSKMatrix::doinit() {
  KMatrix::doinit();
  Energy mK = getParticleData(ParticleID::Kplus)->mass();
  Energy mpi= getParticleData(ParticleID::piplus)->mass();
  sNorm_ = sqr(mK)+sqr(mpi);
}

boost::numeric::ublas::matrix<double> KPiIThreeHalfFOCUSKMatrix::K(Energy2 s,bool) {
  double st = s/sNorm_-1.;
  double param=1.;
  boost::numeric::ublas::matrix<double> output =
    boost::numeric::ublas::zero_matrix<double>(1,1);
  for(unsigned int ix=0;ix<D_.size();++ix) {
    output(0,0) += D_[ix]*param;
    param *= st;
  }
  output *=(s-sThreeHalf_)/sNorm_;
  return output;
}
