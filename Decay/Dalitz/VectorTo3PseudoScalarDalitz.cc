// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorTo3PseudoScalarDalitz class.
//

#include "VectorTo3PseudoScalarDalitz.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/epsilon.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;

IBPtr VectorTo3PseudoScalarDalitz::clone() const {
  return new_ptr(*this);
}

IBPtr VectorTo3PseudoScalarDalitz::fullclone() const {
  return new_ptr(*this);
}

void VectorTo3PseudoScalarDalitz::persistentOutput(PersistentOStream & os) const {
  os << ounit(coupling_,1./GeV);
}

void VectorTo3PseudoScalarDalitz::persistentInput(PersistentIStream & is, int) {
  is >> iunit(coupling_,1./GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<VectorTo3PseudoScalarDalitz,DalitzBase>
describeHerwigVectorTo3PseudoScalarDalitz("Herwig::VectorTo3PseudoScalarDalitz", "HwDalitzDecay.so");

void VectorTo3PseudoScalarDalitz::Init() {

  static ClassDocumentation<VectorTo3PseudoScalarDalitz> documentation
    ("The VectorTo3PseudoScalarDalitz class provides a base class "
     "for the decay of vector mesons to 3 pseudoscalar mesons");

  static Parameter<VectorTo3PseudoScalarDalitz,InvEnergy> interfaceCouopling
    ("Coupling",
     "The coupling for the normalisation of the mode",
     &VectorTo3PseudoScalarDalitz::coupling_, 1./GeV, 1./GeV, 0./GeV, 1000./GeV,
     false, false, Interface::limited);

}

void VectorTo3PseudoScalarDalitz::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  VectorWaveFunction::constructSpinInfo(vectors_,const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  // set up the spin information for the decay products
  for(unsigned int ix=0;ix<3;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
}

double VectorTo3PseudoScalarDalitz::me2(const int ichan, const Particle & part,
				    const tPDVector & ,
				    const vector<Lorentz5Momentum> & momenta,
				    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin0,PDT::Spin0,PDT::Spin0)));
  useMe();
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_,rho_,
  						const_ptr_cast<tPPtr>(&part),
  						incoming,false);
  }
  // set the kinematics
  mD_ = part.mass();
  for(unsigned int ix=0;ix<momenta.size();++ix) {
    mOut_[ix]=momenta[ix].mass();
    for(unsigned int iy=ix+1;iy<momenta.size();++iy) {
      m2_[ix][iy]=(momenta[ix]+momenta[iy]).m();
      m2_[iy][ix]=m2_[ix][iy];
    }
  }
  // now compute the matrix element
  complex<InvEnergy2> amp = amplitude(ichan);
  // polarization vector piece
  LorentzVector<complex<Energy3> > scalar = epsilon(momenta[0],momenta[1],momenta[2]);
  // compute the matrix element
  for(unsigned int ix=0;ix<3;++ix) {
    (*ME())(ix,0,0,0) = Complex(coupling_*amp*scalar.dot(vectors_[ix]));
  }
  // return the answer
  return (ME()->contract(rho_)).real();
}

double VectorTo3PseudoScalarDalitz::
threeBodyMatrixElement(const int , const Energy2 q2,
		       const  Energy2 s3, const Energy2 s2, const Energy2 s1, const 
		       Energy m1, const Energy m2, const Energy m3) const {
  mD_ = sqrt(q2);
  mOut_[0] = m1;
  mOut_[1] = m2;
  mOut_[2] = m3;
  m2_[0][1]=m2_[1][0]=sqrt(s3);
  m2_[0][2]=m2_[2][0]=sqrt(s2);
  m2_[1][2]=m2_[2][1]=sqrt(s1);
  // now compute the matrix element
  // amplitide
  complex<InvEnergy2> amp = amplitude(-1);
  // epsilon piece
  Energy6 kin = (pow<4,1>(m1)*(-2*(sqr(m2) + sqr(m3)) + s1) + pow<4,1>(m2)*(-2*sqr(m3) + s2) +
		 s3*(pow<4,1>(m3) + s1*s2 - sqr(m3)*(s1 + s2 + s3)) - 
     sqr(m1)*(2*pow<4,1>(m2) + 2*pow<4,1>(m3) + sqr(m2)*(4*sqr(m3) - 3*(s1 + s2) - s3) + s1*(s1 + s2 + s3) - 
        sqr(m3)*(3*s1 + s2 + 3*s3)) - sqr(m2)*(2*pow<4,1>(m3) + s2*(s1 + s2 + s3) - sqr(m3)*(s1 + 3*(s2 + s3))))/12.;
  return norm(amp*coupling_*GeV*GeV2)*kin/GeV2/GeV2/GeV2;
} 

WidthCalculatorBasePtr 
VectorTo3PseudoScalarDalitz::threeBodyMEIntegrator(const DecayMode & ) const {
  int imode=0;
  // construct the integrator
  vector<double> inweights;
  vector<int> intype;
  vector<Energy> inmass,inwidth;
  inweights.reserve(resonances().size());
  intype.reserve(resonances().size());
  inmass.reserve(resonances().size());
  inwidth.reserve(resonances().size());
  int iloc=-1;
  vector<double> inpow(2,0.0);
  for(unsigned int ix=0;ix<resonances().size();++ix) {
    tPDPtr resonance = getParticleData(resonances()[ix]->id);
    if(resonance) {
      ++iloc;
      inweights.push_back(weights()[iloc]);
      inmass.push_back(resonances()[ix]->mass);
      inwidth.push_back(abs(resonances()[ix]->width));
      intype.push_back(resonances()[ix]->spectator+1);
    }
  }
  tcDecayIntegratorPtr decayer(this);
  WidthCalculatorBasePtr output(
    new_ptr(ThreeBodyAllOnCalculator<VectorTo3PseudoScalarDalitz>
  	    (inweights,intype,inmass,inwidth,inpow,
  	     *this,imode,
	     mode(0)->outgoing()[0]->mass(),
	     mode(0)->outgoing()[1]->mass(),
	     mode(0)->outgoing()[2]->mass())));
  return output;
}

void VectorTo3PseudoScalarDalitz::dataBaseOutput(ofstream & output,
					     bool header) const {
  if(header){output << "update decayers set parameters=\"";}
  output << "newdef " << name() << ":Coupling " << coupling_*GeV << "\n";
  // parameters for the DalitzBase base class
  DalitzBase::dataBaseOutput(output,false);
  if(header){output << "\n\" where BINARY ThePEGName=\"" 
   		    << fullName() << "\";" << endl;}
}

complex<InvEnergy2> VectorTo3PseudoScalarDalitz::resAmp(unsigned int i) const {
  // can't have a scalar here on spin/parity grounds
  assert(resonances()[i]->type%10!=1);
  // shouldn't have E691 stuff either
  assert(resonances()[i]->type%10!=1);
  // amplitude
  Complex output = resonances()[i]->amp;
  if (resonances()[i]->type==ResonanceType::NonResonant) return output/GeV2;
  // mass of the resonance
  const Energy & mR = resonances()[i]->mass ;
  // locations of the outgoing particles
  const int &d1 = resonances()[i]->daughter1;
  const int &d2 = resonances()[i]->daughter2;
  const int &sp = resonances()[i]->spectator;
  // epsilon piece  =eps(d1,d2,sp)
  double sign = (sp-d1)*(d1-d2)*(d2-sp)/2.;
  // compute the Breit-Wigner times resonance form factor piece
  output *= sign*resonances()[i]->BreitWigner(m2_[d1][d2],mOut_[d1],mOut_[d2]);
  // Blatt-Weisskopf factors
  // for the D decay
  Energy pD  = sqrt(max(ZERO,(0.25*sqr(sqr(mD_)-sqr(mR)-sqr(mOut_[sp])) - sqr(mR*mOut_[sp]))/sqr(mD_)));
  Energy pDAB= sqrt( 0.25*sqr(sqr(mD_)-sqr(m2_[d1][d2])-sqr(mOut_[sp])) - sqr(m2_[d1][d2]*mOut_[sp]))/mD_;
  double r2A(parentRadius()   *pD),r2B(parentRadius()   *pDAB);
  // Blatt-Weisskopf factors and spin piece
  switch (resonances()[i]->type) {
  case ResonanceType::Spin1: case ResonanceType::Spin1GS : 
    output *= sqrt( (1. + sqr(r2A)) / (1. + sqr(r2B)) );
    break;
  case ResonanceType::Spin2:
    output *= sqrt( (9. + sqr(r2A)*(3.+sqr(r2A))) / (9. + sqr(r2B)*(3.+sqr(r2B))));
    // spin piece
    output *= (sqr(mD_) - sqr(mOut_[sp]) + sqr(m2_[d1][d2]))/(sqr(mD_)*sqr(m2_[d1][d2]))/GeV2*
	       ((-sqr(m2_[sp][d1]) + sqr(m2_[sp][d2]))*sqr(m2_[d1][d2]) +
		(mD_ - mOut_[sp])*(mD_ + mOut_[sp])*(mOut_[d1] - mOut_[d2])*(mOut_[d1] + mOut_[d2]));
    break;
  default :
    assert(false);
  }
  return output/GeV2;
}
