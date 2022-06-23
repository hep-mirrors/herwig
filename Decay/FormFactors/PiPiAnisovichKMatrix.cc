// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PiPiAnisovichKMatrix class.
//

#include "PiPiAnisovichKMatrix.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

PiPiAnisovichKMatrix::PiPiAnisovichKMatrix()
  : KMatrix(FlavourInfo(IsoSpin::IZero,IsoSpin::I3Unknown,
			Strangeness::Zero,Charm::Zero,
			Beauty::Zero),
	    vector<Channels>({KMatrix::PiPi,KMatrix::KK,KMatrix::FourPi,KMatrix::EtaEta,KMatrix::EtaEtaPrime}),
	    vector<Energy2>({sqr(0.65100*GeV), sqr(1.20720*GeV), sqr(1.56122*GeV), sqr(1.21257*GeV), sqr(1.81746*GeV) }),
	    vector<vector<Energy> >({{0.24844*GeV,-0.52523*GeV,0.00000*GeV,-0.38878*GeV,-0.36397*GeV},
				     {0.91779*GeV,0.55427*GeV,0.00000*GeV,0.38705*GeV,0.29448*GeV},
				     {0.37024*GeV,0.23591*GeV,0.62605*GeV,0.18409*GeV,0.18923*GeV},
				     {0.34501*GeV,0.39642*GeV,0.97644*GeV,0.19746*GeV,0.00357*GeV},
				     {0.15770*GeV,-0.17915*GeV,-0.90100*GeV,-0.00931*GeV,0.20689*GeV}})),
    s0Scatt_(-3.30564*GeV2), f1a_({0.26681,0.16583,-0.19840,0.32808,0.31193}),sA_(1.),sA0_(-0.2)
{}

IBPtr PiPiAnisovichKMatrix::clone() const {
  return new_ptr(*this);
}

IBPtr PiPiAnisovichKMatrix::fullclone() const {
  return new_ptr(*this);
}

void PiPiAnisovichKMatrix::persistentOutput(PersistentOStream & os) const {
  os << ounit(s0Scatt_,GeV2) << f1a_ << sA_ << sA0_ << ounit(mPi_,GeV);
}

void PiPiAnisovichKMatrix::persistentInput(PersistentIStream & is, int) {
  is >> iunit(s0Scatt_,GeV2) >> f1a_ >> sA_ >> sA0_ >> iunit(mPi_,GeV);
}

void PiPiAnisovichKMatrix::doinit() {
  KMatrix::doinit();
  mPi_ = getParticleData(211)->mass();
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<PiPiAnisovichKMatrix,KMatrix>
describeHerwigPiPiAnisovichKMatrix("Herwig::PiPiAnisovichKMatrix", "HwFormFactors.so");

void PiPiAnisovichKMatrix::Init() {

  static ClassDocumentation<PiPiAnisovichKMatrix> documentation
    ("The PiPiAnisovichKMatrix class implements the K matrix pararmeterization from Eur.Phys.J.A 16 (2003) 229-258",
     "The K-matrix parameterization of Anisovich and Sarantsev \\cite{Anisovich:2002ij} was used.",
     "\bibitem{Anisovich:2002ij}\n"
     "V.~V.~Anisovich and A.~V.~Sarantsev,"
     "%``K matrix analysis of the (I J**(PC) = 00++)-wave in the mass region below 1900 MeV,''\n"
     "Eur. Phys. J. A \\textbf{16} (2003), 229-258\n"
     "doi:10.1140/epja/i2002-10068-x\n"
     "[arXiv:hep-ph/0204328 [hep-ph]].\n"
     "%170 citations counted in INSPIRE as of 20 Jun 2022\n");
  
  static Parameter<PiPiAnisovichKMatrix,Energy2> interfaces0Scatt
    ("s0Scatt",
     "The s0 parameters for scattering",
     &PiPiAnisovichKMatrix::s0Scatt_, GeV2, -3.30564*GeV2, -100.0*GeV2, 100.0*GeV2,
     false, false, Interface::limited);

  static ParVector<PiPiAnisovichKMatrix,double> interfaceF1A
    ("F1A",
     "The non-pole coefficients",
     &PiPiAnisovichKMatrix::f1a_, 5, 1.0, -100, 100.,
     false, false, Interface::limited);

  static Parameter<PiPiAnisovichKMatrix,double> interfacesA
    ("sA",
     "The sA parameter",
     &PiPiAnisovichKMatrix::sA_, 1., -10., 10.0,
     false, false, Interface::limited);

  static Parameter<PiPiAnisovichKMatrix,double> interfacesA0
    ("sA0",
     "The sA0 parameter",
     &PiPiAnisovichKMatrix::sA0_, -0.2, -10., 10.0,
     false, false, Interface::limited);
}

boost::numeric::ublas::matrix<double> PiPiAnisovichKMatrix::K(Energy2 s, bool multiplyByPoles) const {
  double pre = (s-0.5*sA_*sqr(mPi_))*(1.-sA0_)/(s-sA0_*GeV2);
  double coeff = (GeV2-s0Scatt_)/(s-s0Scatt_);
  if(multiplyByPoles)
    for (Energy2 pole : poles() ) coeff *= 1.-s/pole;
  boost::numeric::ublas::matrix<double> output =
    boost::numeric::ublas::zero_matrix<double>(5,5);
  for(unsigned int im=0;im<poles().size();++im) {
    InvEnergy2 term;
    if(multiplyByPoles) {
      term = 1./poles()[im];
      for(unsigned int iz=0;iz<poles().size();++iz) {
	if(iz==im) continue;
	term *= 1. - s/poles()[iz];
      }
    }
    else
      term = 1./(poles()[im]-s);
    for(unsigned int ix=0;ix<5;++ix)
      for(unsigned int iy=ix;iy<5;++iy)
	output(ix,iy) += term*poleCouplings()[im][ix]*poleCouplings()[im][iy];
  }
  for(unsigned int iy=0;iy<5;++iy) output(0,iy) += coeff*f1a_[iy];
  for(unsigned int ix=0;ix<5;++ix)
    for(unsigned int iy=ix+1;iy<5;++iy)
      output(iy,ix)=output(ix,iy);
  output *= pre;
  return output;
}
