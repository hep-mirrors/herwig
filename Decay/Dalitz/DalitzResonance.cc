// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DalitzResonance class.
//

#include "DalitzResonance.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/EnumIO.h"
#include "FlatteResonance.h"
#include "FlatteResonance2.h"
#include "MIPWA.h"
#include "PiPiI2.h"
#include "DalitzKMatrix.h"
#include "DalitzLASS.h"
#include "DalitzGS.h"
#include "DalitzSigma.h"

using namespace Herwig;

void DalitzResonance::persistentOutput(PersistentOStream & os) const {
  os << id << oenum(type) << ounit(mass,GeV) << ounit(width,GeV)
     << daughter1 << daughter2 << spectator
     << amp << ounit(R,1./GeV);
}

void DalitzResonance::persistentInput(PersistentIStream & is, int) {
  is >> id >> ienum(type) >> iunit(mass,GeV) >> iunit(width,GeV)
     >> daughter1 >> daughter2 >> spectator
     >> amp >> iunit(R,1./GeV);
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<DalitzResonance,Base>
  describeHerwigDalitzResonance("Herwig::DalitzResonance", "HwDalitzDecay.so");

void DalitzResonance::Init() {

  static ClassDocumentation<DalitzResonance> documentation
    ("The DalitzResonance class provides a container class for"
     " information on resonances in multi-body dalitz decays.");

}

Complex DalitzResonance::BreitWigner(const Energy & mAB, const Energy & mA, const Energy & mB) const {
  static const Complex ii = Complex(0.,1.);
  // non-resonant pieces
  if(abs(type)/10==10) return 1.;
  // momenta for the resonance decay
  // off-shell
  Energy pAB=sqrt(0.25*sqr(sqr(mAB) -sqr(mA)-sqr(mB)) - sqr(mA*mB))/mAB;
  if(type==ResonanceType::BABARf0) {
    double rho = 2.*pAB/mAB;
    return GeV2/(sqr(mass)-sqr(mAB)-ii*mass*width*rho);
  }
  else if (type==ResonanceType::Spin0Complex) {
    complex<Energy> sR(mass,width);
    return GeV2/(sqr(sR)-sqr(mAB));
  }
  else if (type==ResonanceType::Flattef0  ||
	   type==ResonanceType::Flatte2a0 ||
	   type==ResonanceType::Flatte2Kstar0) {
    assert(false);
  }
  //  on-shell
  Energy  pR=sqrt(0.25*sqr( mass*mass - sqr(mA) - sqr(mB)) - sqr(mA*mB))/mass;
  // Blatt-Weisskopf factors
  double fR=1;
  unsigned int power(1);
  if(type!=ResonanceType::Spin0 &&
     type!=ResonanceType::Spin0E691) {
    double r1A(R*pR),r1B(R*pAB);
    // Blatt-Weisskopf factors and spin piece
    switch (type) {
    case ResonanceType::Spin0Gauss:
      fR = exp(-(r1B-r1A)/12.);
      // cerr << "testing scalar B " <<   exp(+r1A/12.) << "\n";
      break;
    case ResonanceType::Spin1: case ResonanceType::Spin1E691 :
      fR=sqrt( (1. + sqr(r1A)) / (1. + sqr(r1B)) );
      // cerr << "testing vector B " << sqrt(1. + sqr(r1A))   << "\n";
      power=3;
      break;
    case ResonanceType::Spin2: case ResonanceType::Spin2E691:
      fR = sqrt( (9. + sqr(r1A)*(3.+sqr(r1A))) / (9. + sqr(r1B)*(3.+sqr(r1B))));
      // cerr << "testing tensor B " <<  sqrt( (9. + sqr(r1A)*(3.+sqr(r1A))))  << "\n";
      power=5;
      break;
    default :
      assert(false);
    }
  }
  // multiply by Breit-Wigner piece and return
  if (type/10 == 1 ) {
    return fR*sqrt(0.5*width/GeV/Constants::pi)*GeV/(mAB-mass-complex<Energy>(ZERO,0.5*width));
  }
  else {
    Energy gam = width*pow(pAB/pR,power)*(mass/mAB)*fR*fR;
    return fR*GeV2/(sqr(mass)-sqr(mAB)-mass*gam*ii);
  }
}

void DalitzResonance::dataBaseOutput(ofstream & output) {
  output << id << " " << oenum(type) << " "
	 << mass/GeV << " " << width/GeV << " "
	 << daughter1 << " " << daughter2 << " "
	 << spectator << " " 
	 << abs(amp) << " " << arg(amp) << " "
	 << R*GeV; 
}

DalitzResonancePtr DalitzResonance::readResonance(string arg, string & error) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long id = stoi(stype);
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  ResonanceType::Type type = static_cast<ResonanceType::Type>(stoi(stype));
  // mass and width
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  Energy mass = stof(stype)*GeV;
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  Energy width = stof(stype)*GeV;
  // children
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  int d1 = stoi(stype);
  // children
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  int d2 = stoi(stype);
  // children
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  int sp = stoi(stype);
  if (sp==d1 || sp ==d2 || d1 == d2) {
    error =  "Daughters and spectator must all be different not " + std::to_string(d1) + ", " + std::to_string(d2) + ", " + std::to_string(sp);
    return DalitzResonancePtr();
  }
  // magnitude and phase
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double mag = stof(stype);
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double phi = stof(stype);
  // radius
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  InvEnergy r = stof(stype)/GeV;
  // special for flate
  if (type==ResonanceType::Flattef0) {
    // Flatte parameters
    // magnitude and phase
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    double fpi = stof(stype);
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    double fK  = stof(stype);
    // add to list
    return new_ptr(FlatteResonance(id,type,mass,width,d1,d2,sp,mag,phi,r,fpi,fK));
  }
  else if( type==ResonanceType::Flatte2a0 ||
	   type==ResonanceType::Flatte2Kstar0) {
    // Flatte parameters
    // magnitude and phase
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    Energy fpi = stof(stype)*GeV;
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    Energy fK  = stof(stype)*GeV;
    // add to list
    return new_ptr(FlatteResonance2(id,type,mass,width,d1,d2,sp,mag,phi,r,fpi,fK));
  }
  // MIPWA
  else if(type==ResonanceType::Spin0MIPWA) {
    // no of entries in table
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    int nn = stoi(stype);
    vector<Energy> en; en.reserve(nn);
    vector<double> mag2,phase2; mag2.reserve(nn); phase2.reserve(nn);
    for(int ix=0;ix<nn;++ix) {
      stype = StringUtils::car(arg);
      arg   = StringUtils::cdr(arg);
      Energy ee = stof(stype)*GeV;
      en.push_back(ee);
      stype = StringUtils::car(arg);
      arg   = StringUtils::cdr(arg);
      double mm=stof(stype);
      mag2.push_back(mm);
      stype = StringUtils::car(arg);
      arg   = StringUtils::cdr(arg);
      double pp=stof(stype);
      phase2.push_back(pp);
    }
    return new_ptr(MIPWA(id,type,mass,width,d1,d2,sp,mag,phi,r,en,mag2,phase2));
  }
  // I=2 pipi
  else if(type==ResonanceType::PiPiI2) {
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    InvEnergy a = stof(stype)/GeV;
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    InvEnergy2 b = stof(stype)/GeV2;
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    InvEnergy4 c = stof(stype)/GeV2/GeV2;
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    InvEnergy6 d = stof(stype)/GeV2/GeV2/GeV2;
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    Energy mmin = stof(stype)*GeV;
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    Energy mmax = stof(stype)*GeV;
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    double dEta=stof(stype);
    return new_ptr(PiPiI2(id,type,mass,width,d1,d2,sp,mag,phi,r,
			  a,b,c,d,mmin,mmax,dEta));
  }
  // K-matrix
  else if(type==ResonanceType::KMatrix) {
    // no of poles
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    unsigned int npole= stoi(stype);
    // no of channels
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    unsigned int nchannels= stoi(stype);
    // matrix location
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    unsigned int imat = stoi(stype);
    // this channel
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    unsigned int chan = stoi(stype);
    // expansion point for the constants terms
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    Energy2 sc = GeV2*stof(stype);
    // type of expansion
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    unsigned int itype= stoi(stype);
    vector<pair<double,double> > beta;
    // first loop over the coefficients of the poles
    for(unsigned int ix=0;ix<npole;++ix) {
      stype = StringUtils::car(arg);
      arg   = StringUtils::cdr(arg);
      double b = stof(stype);
      stype = StringUtils::car(arg);
      arg   = StringUtils::cdr(arg);
      beta.push_back(make_pair(b,stof(stype)));
    }
    // now over the power series for the different channels
    vector<pair<double,vector<double > > > coeffs(nchannels);
    for(unsigned int ix=0;ix<nchannels;++ix) {
      stype = StringUtils::car(arg);
      arg   = StringUtils::cdr(arg);
      unsigned int nterms = stoi(stype);
      for(unsigned int iy=0;iy<nterms;++iy) {
    	stype = StringUtils::car(arg);
    	arg   = StringUtils::cdr(arg);
    	coeffs[ix].second.push_back(stof(stype));
      }
      stype = StringUtils::car(arg);
      arg   = StringUtils::cdr(arg);
      coeffs[ix].first = stof(stype);
    }
    // finally make the channel
    return new_ptr(DalitzKMatrix(id,type,mass,width,d1,d2,sp,mag,phi,r,imat,chan,sc,itype,beta,coeffs));
  }
  // LASS
  else if(type==ResonanceType::LASS) {
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    unsigned int iopt = stoi(stype);
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    double FNR = stof(stype);
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    double phiNR = stof(stype);
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    double FRes = stof(stype);
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    double phiRes = stof(stype);
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    InvEnergy ascat = stof(stype)/GeV;
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    InvEnergy reff = stof(stype)/GeV;
    // finally make the channel
    return new_ptr(DalitzLASS(id,type,mass,width,d1,d2,sp,mag,phi,r,iopt,
			      FNR,phiNR,FRes,phiRes,ascat,reff));
  }
  // Bugg sigma form
  else if(type==ResonanceType::Sigma) {
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    Energy2 a = stof(stype)*GeV2;
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    Energy b1 = stof(stype)*GeV;
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    InvEnergy b2 = stof(stype)/GeV;
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    Energy g4pi = stof(stype)*GeV;
    return new_ptr(DalitzSigma(id,type,mass,width,d1,d2,sp,mag,phi,r,a,b1,b2,g4pi));
  }
  // GS form
  else if(type==ResonanceType::Spin1GS) {
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    Energy mpi = stof(stype)*GeV;
    return new_ptr(DalitzGS(id,type,mass,width,d1,d2,sp,mag,phi,r,mpi));
  }
  // otherwise add to list
  else {
    return new_ptr(DalitzResonance(id,type,mass,width,d1,d2,sp,mag,phi,r));
  }
}
