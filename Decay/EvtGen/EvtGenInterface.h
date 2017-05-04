// -*- C++ -*-
#ifndef Herwig_EvtGenInterface_H
#define Herwig_EvtGenInterface_H
//
// This is the declaration of the EvtGenInterface class.
//

#include "EvtGenInterface.fh"
#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/Vectors/Lorentz5Vector.h"
#include "ThePEG/Helicity/ScalarSpinInfo.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "ThePEG/Helicity/RSFermionSpinInfo.h"
#include "ThePEG/Helicity/TensorSpinInfo.h"
#include "ThePEG/EventRecord/Particle.h"

#include "EvtGenRandom.h"
#include "EvtGen/EvtGen.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtSpinDensity.hh"
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtRaritaSchwinger.hh"
#include "EvtGenBase/EvtDecayAmp.hh"

namespace Herwig {

using namespace ThePEG;

/**
 * The EvtGenInterface class is the main class for the use of the EvtGen decay
 * package with Herwig.
 *
 * @see \ref EvtGenInterfaceInterfaces "The interfaces"
 * defined for EvtGenInterface.
 */
class EvtGenInterface: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  EvtGenInterface();

  /**
   * The copy constructor (explicit as cannot copy streams)
   */
  EvtGenInterface(const EvtGenInterface &);

  /**
   * The destructor.
   */
  virtual ~EvtGenInterface();
  //@}

public:

  /**
   * Use EvtGen to perform a decay
   * @param parent The decaying particle
   * @param dm The decaymode
   * @return The decay products
   */
  ParticleVector decay(const Particle &parent,
		       bool recursive, const DecayMode & dm) const;

public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /** @name Functions to convert between EvtGen and Herwig classes */
  //@{
  /**
   * Convert a particle to an EvtGen particle.
   * @param part The particle to be converted.
   */
  EvtParticle *EvtGenParticle(const Particle & part) const;

  /**
   *  Return the decay products of an EvtGen particle in as ThePEG particles
   * @param evtpart The EvtGen particle
   */
  ParticleVector decayProducts(EvtParticle* evtpart, bool boost=true) const;

  /**
   * Convert a Lorentz5Momentum to a real EvtGen 4-vector
   * @param mom The momentum to be converted
   */
  EvtVector4R EvtGenMomentum(const Lorentz5Momentum & mom) const {
    return EvtVector4R(mom.t()/GeV,mom.x()/GeV,mom.y()/GeV,mom.z()/GeV);
  }

  /**
   * Convert a PDG code from ThePEG into an EvtGen particle id
   * @param id The PDG code
   * @param exception Whether or not to throw an Exception if fails
   */
  EvtId EvtGenID(int id,bool exception=true) const;

  /**
   * Convert a LorentzSpinor to an EvtGen one. The spinor is converted to the 
   * EvtGen Dirac representation/
   * @param sp The LorentzSpinor
   */
  EvtDiracSpinor EvtGenSpinor(const LorentzSpinor<SqrtEnergy> & sp) const {
    InvSqrtEnergy norm(sqrt(0.5)/sqrt(GeV));
    EvtDiracSpinor output;
    output.set(EvtGenComplex(-norm*( sp.s1()+sp.s3())),
	       EvtGenComplex(-norm*( sp.s2()+sp.s4())),
	       EvtGenComplex(-norm*(-sp.s1()+sp.s3())),
	       EvtGenComplex(-norm*(-sp.s2()+sp.s4())));
    return output;
  }

  /**
   * Convert a LorentzPolarizationVector to a complex EvtGen 4-vector
   * @param eps The polarization vector to be converted
   */
  EvtVector4C EvtGenPolarization(const LorentzPolarizationVector & eps) const {
    return EvtVector4C(EvtGenComplex(eps.t()),EvtGenComplex(eps.x()),
		       EvtGenComplex(eps.y()),EvtGenComplex(eps.z()));
  }
  
  /**
   * Convert our Rarita-Schwinger spinor to the EvtGen one
   * @param sp Our  RS Spinor
   */
  EvtRaritaSchwinger EvtGenRSSpinor(const LorentzRSSpinor<SqrtEnergy> & sp) const {
    InvSqrtEnergy norm(sqrt(0.5)/sqrt(GeV));
    complex<double> out[4][4];
    for(unsigned int ix=0;ix<4;++ix) {
      out[ix][0] =-norm*( sp(ix,0)+sp(ix,2));
      out[ix][1] =-norm*( sp(ix,1)+sp(ix,3));
      out[ix][2] =-norm*(-sp(ix,0)+sp(ix,2));
      out[ix][3] =-norm*(-sp(ix,1)+sp(ix,3));
    }
    EvtRaritaSchwinger output;
    unsigned int ix,iy;
    // remember we have vec,spin and evtgen spin,vec
    for(ix=0;ix<4;++ix) {
      for(iy=0;iy<4;++iy) output.set(ix,iy,EvtGenComplex(out[iy][ix]));
    }
    return output;
  }
  
  /**
   * Convert our tensor to the EvtGen one.
   * @param ten Our tensor
   */
  EvtTensor4C EvtGenTensor(const LorentzTensor<double> & ten) const {
    EvtTensor4C output;
    unsigned int ix,iy;
    for(ix=0;ix<4;++ix){
      for(iy=0;iy<4;++iy) output.set(ix,iy,EvtGenComplex(ten(ix,iy)));
    }
    return output;
  }

  /**
   * Convert a spin density matrix to an EvtGen spin density matrix.
   * @param rho The spin density matrix to be converted.
   */
  EvtSpinDensity EvtGenSpinDensity(const RhoDMatrix & rho) const {
    EvtSpinDensity rhoout;
    unsigned int ix,iy,ispin(rho.iSpin());
    rhoout.setDim(ispin);
    for(ix=0;ix<ispin;++ix) {
      for(iy=0;iy<ispin;++iy)
	rhoout.set(ix,iy,EvtGenComplex(rho(ix,iy)));
    }
    return rhoout;
  }

  /**
   * Convert from our complex to the EvtGen one
   */
  EvtComplex EvtGenComplex(Complex c) const {
    return EvtComplex(c.real(),c.imag());
  }
  //@}

  /**
   *  Functions to convert between EvtGen and Herwig classes
   */
  //@{
  /**
   * Convert a particle from an EvtGen one to ThePEG one.
   * @param part The EvtGen particle.
   * @param pd Pointer to the particle data object of ThePEG for the particle.
   * @param spin Convert the spin information as well
   */
  PPtr ThePEGParticle(EvtParticle *part, tcPDPtr pd,bool spin=true) const {
    PPtr output(new_ptr(Particle(pd)));
    output->set5Momentum(ThePEGMomentum(part->getP4(),part->mass()));
    if(spin) ThePEGSpin(output,part);
    // make the daughters 
    ParticleVector daug(decayProducts(part,spin));
    for(unsigned int ix=0;ix<daug.size();++ix) output->addChild(daug[ix]);
    return output;
  }

  /**
   * Set the SpinInfo for a ThePEG particle using an EvtGen particle
   * @param pegpart ThePEG particle.
   * @param evtpart The EvtGen particle.
   */
  void ThePEGSpin(PPtr pegpart,EvtParticle *evtpart) const;

  /**
   * Convert an EvtGen EvtId to a PDG code in our conventions
   * @param id The EvtGen ID.
   * @param exception Whether or not to throw an Exception if fails
   */
  int ThePEGID(EvtId id,bool exception=true) const;

  /**
   * Convert from EvtGen momentum to Lorentz5Momentum
   * @param mom The EvtGen 4-momentum
   * @param mass The mass
   */
  Lorentz5Momentum ThePEGMomentum(const EvtVector4R & mom,double mass) const  {
    return Lorentz5Momentum(mom.get(1)*GeV,mom.get(2)*GeV,
			    mom.get(3)*GeV,mom.get(0)*GeV,mass*GeV);
  }
  /**
   * Convert from EvtGen complex to ours
   */
  Complex ThePEGComplex(EvtComplex c) const {
    return Complex(real(c),imag(c));
  }

  /**
   * Convert a spin density to a ThePEG one from an EvtGen one
   * @param rho The spin density matrix to be converted
   * @param id The PDG code of the particle to get special cases right.
   */
  RhoDMatrix ThePEGSpinDensity(const EvtSpinDensity & rho, int id) const;
  
  /**
   * Convert an EvtDiracSpinor a LorentzSpinor. This spinor is converted to 
   * the default Dirac matrix representation used by ThePEG.
   * @param sp The EvtDiracSpinor
   */
  LorentzSpinor<SqrtEnergy> ThePEGSpinor(const EvtDiracSpinor & sp) const {
    SqrtEnergy norm(sqrt(0.5)*sqrt(GeV));
    vector<complex<SqrtEnergy> > evtSpin(4);
    for(unsigned int ix=0;ix<4;++ix) evtSpin[ix] = -norm*ThePEGComplex(sp.get_spinor(ix));
    return LorentzSpinor<SqrtEnergy>(evtSpin[0]-evtSpin[2],evtSpin[1]-evtSpin[3],
				     evtSpin[0]+evtSpin[2],evtSpin[1]+evtSpin[3]);
  }

  /**
   * Convert an EvtGen complex 4-vector to a LorentzPolarizationVector
   * @param eps The complex 4-vector to be converted.
   */
  LorentzPolarizationVector ThePEGPolarization(const EvtVector4C & eps) const {
    return LorentzPolarizationVector(conj(ThePEGComplex(eps.get(1))),
				     conj(ThePEGComplex(eps.get(2))),
				     conj(ThePEGComplex(eps.get(3))),
				     conj(ThePEGComplex(eps.get(0))));
  }
  
  /**
   * Convert an EvtGen Rarita-Schwinger spinor to ours
   * @param sp The EvtGen RS spinor.
   */
  LorentzRSSpinor<SqrtEnergy> ThePEGRSSpinor(const EvtRaritaSchwinger & sp) const {
    complex<SqrtEnergy> evtSpin[4][4];
    SqrtEnergy norm(sqrt(0.5)*sqrt(GeV));
    // normalisation and swap vec,spin order
    for(unsigned int ix=0;ix<4;++ix) {
      for(unsigned int iy=0;iy<4;++iy) evtSpin[ix][iy]=-norm*ThePEGComplex(sp.get(iy,ix));
    }
    LorentzRSSpinor<SqrtEnergy> output;
    for(unsigned int ix=0;ix<4;++ix) {
      output(ix,0) = evtSpin[ix][0] - evtSpin[ix][2];
      output(ix,1) = evtSpin[ix][1] - evtSpin[ix][3];
      output(ix,2) = evtSpin[ix][0] + evtSpin[ix][2];
      output(ix,3) = evtSpin[ix][1] + evtSpin[ix][3];
    }
    // output.changeRep(Helicity::defaultDRep);
    return output;
  }
  
  /**
   * Convert an EvtGen tensor to ThePEG
   * @param ten The EvtGen tensor
   */
  LorentzTensor<double> ThePEGTensor(const EvtTensor4C & ten) const {
    LorentzTensor<double> output;
    unsigned int ix,iy;
    for(ix=0;ix<4;++ix) {
      for(iy=0;iy<4;++iy)output(ix,iy)=conj(ThePEGComplex(ten.get(ix,iy)));
    }
    return output;
  }
  //@}

  /**
   *  Check the conversion of particles between Herwig and EvtGen
   */
  void checkConversion() const;

  /**
   * Output the EvtGen decay modes for a given particle
   * @param id The PDG code of the particle to output
   */
  void outputEvtGenDecays(long id) const;

  /**
   *  Find the location in the EvtGen list of decay channels for
   *  a given decay mode.
   */
  int EvtGenChannel(const DecayMode &dm) const;

  /**
   * Check the particle has SpinInfo and if not create it
   * @param part The particle
   */
  tSpinPtr getSpinInfo(const Particle &part) const {
    // return spin info if exists
    if(part.spinInfo()) {
      return dynamic_ptr_cast<tSpinPtr>(const_ptr_cast<tPPtr>(&part)->spinInfo());
    }
    // otherwise make it
    tPPtr ptemp(const_ptr_cast<tPPtr>(&part));
    PDT::Spin spin(part.dataPtr()->iSpin());
    SpinPtr pspin;
    if(spin==PDT::Spin0)          pspin=new_ptr(ScalarSpinInfo(part.momentum(),true));
    else if(spin==PDT::Spin1Half) pspin=new_ptr(FermionSpinInfo(part.momentum(),true));
    else if(spin==PDT::Spin1)     pspin=new_ptr(VectorSpinInfo(part.momentum(),true));
    else if(spin==PDT::Spin3Half) pspin=new_ptr(RSFermionSpinInfo(part.momentum(),true));
    else if(spin==PDT::Spin2)     pspin=new_ptr(TensorSpinInfo(part.momentum(),true));
    else throw Exception() << "Can't create spinInfo for decaying particle in "
			   << "EvtGen::checkSpinInfo for spin " << spin << "particle " 
			   << Exception::eventerror;
    ptemp->spinInfo(pspin);
    return pspin;
  }

  /**
   *  Construct the DecayVertex for Herwig using the information from
   *  EvtGen
   * @param parent The decaying particle
   * @param products The outgoing particles
   * @param damp Pointer to the EvtGen decayer
   */
  void constructVertex(const Particle & parent,ParticleVector products,
		       EvtDecayAmp* damp) const;

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  EvtGenInterface & operator=(const EvtGenInterface &);

private:

  /**
   *    Names of the various EvtGen parameter files
   */
  //@{
  /**
   *  The name of the file containing the decays
   */
  string decayName_;

  /**
   *  The name of the file containing the particle data
   */
  string pdtName_;

  /**
   *  Names of addition user specified decays
   */
  vector<string> userDecays_; 
  //@}

  /**
   *  Whether or not to redirect cout and cerr when EvtGen is running
   */
  bool reDirect_;

  /**
   *  Check the conversion of the particles
   */
  bool checkConv_;

  /**
   *  Particles for which to output the EvtGen decays
   */
  vector<long> convID_;

  /**
   *  Location of the PYTHIA8 data directory
   */
  string p8Data_;

private:

  /**
   * Pointer to the random number generator for EvtGen
   */
  EvtRandomEngine * evtrnd_;

  /** 
   * Main EvtGen object
   */
  EvtGen * evtgen_;

  /**
   *  File to output the log info to
   */
  mutable ofstream logFile_;

};

}

#endif /* Herwig_EvtGenInterface_H */
