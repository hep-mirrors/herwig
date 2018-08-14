// -*- C++ -*-
#ifndef Herwig_PhaseSpaceChannel_H
#define Herwig_PhaseSpaceChannel_H
//
// This is the declaration of the PhaseSpaceChannel class.
//

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/EnumIO.h"
#include "PhaseSpaceMode.fh"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"

namespace Herwig {

using namespace ThePEG;

  /** \ingroup Decay
   *
   * This class is designed to store the information needed for a given
   * phase-space channel for use by the multi-channel phase space decayer
   * and perform the generation of the phase space for that channel.
   *
   * The decay channel is specified as a series of \f$1\to2\f$ decays to either
   * external particles or other intermediates. For each intermediate
   * the Jacobian to be used can be either a Breit-Wigner(0) or a power-law
   * (1).
   *
   * The class is then capable of generating a phase-space point using this
   * channel and computing the weight of a given point for use in a multi-channel
   * phase space integration using the <code>PhaseSpaceMode</code> class.
   *
   * The class is designed so that the phase-space channels can either by specified
   * using the addIntermediate method directly or via the repository.
   * (In practice at the moment all the channels are constructed by the relevant decayers
   *  using the former method at the moment.)
   *
   * @see PhaseSpaceMode
   * @see DecayIntegrator2
   *
   * @author  Peter Richardson
   */
class PhaseSpaceChannel  {

public:

   friend PersistentOStream;
   friend PersistentIStream;
  
/**
 *  Struct for the intermediates in the phase-space channel
 */
struct PhaseSpaceResonance {

  /**
   *   Enum for the jacobians
   */
  enum Jacobian {BreitWigner,Power,OnShell};
  
   /**
    *   The constructor
    */
PhaseSpaceResonance() {};
   /**
    *   The constructor
    */
PhaseSpaceResonance(cPDPtr part) :particle(part), mass2(sqr(part->mass())), mWidth(part->mass()*part->width()),
				  jacobian(BreitWigner), power(0.), children(make_pair(0,0))
{};
  /**
   *  The particle
   */
  cPDPtr particle;

  /**
   *  Mass squared
   */
  Energy2 mass2;

  /**
   *  Mass times width
   */
  Energy2 mWidth;

  /**
   *   Type of jacobian
   */
  Jacobian jacobian;

  /**
   *   The power for a power law
   */
  double power;
  
  /**
   *  The children
   */
  pair<int,int> children;
  
  /**
   *   The final descendents
   */
  vector<int> descendents;
  
};
  
public:
  
  /**
   *  Default constructior
   */
  PhaseSpaceChannel() : weight_(1.)  {};
  
  /** 
   *  Constructor with incoming particles
   */
  PhaseSpaceChannel(tPhaseSpaceModePtr inm);
  
  /**
   * If less than zero indicate that this channel is competed. Otherwise
   * signal the parent of the next added parton.
   */
  PhaseSpaceChannel & operator , (tPDPtr res) {
    intermediates_.push_back(PhaseSpaceResonance(res));
    if(iAdd_<0) return *this;
    if(intermediates_[iAdd_].children.first==0)
      intermediates_[iAdd_].children.first  = 1-int(intermediates_.size());
    else
      intermediates_[iAdd_].children.second = 1-int(intermediates_.size());
    iAdd_=-1;
    return *this;
  }
  
  /**
   * If less than zero indicate that this channel is competed. Otherwise
   * signal the parent of the next added parton.
   */
  PhaseSpaceChannel & operator , (int o) {
    if(iAdd_<0) iAdd_ = o;
    else if(o>=0) {
      if(intermediates_[iAdd_].children.first==0)
	intermediates_[iAdd_].children.first  = o;
      else
	intermediates_[iAdd_].children.second = o;
      iAdd_=-1;
    }
    else if(o<0) {
      assert(false);
    }
    return *this;
  }
  
  /**
   *  Set the jacobian for a given resonance
   */
  void setJacobian(unsigned int ires, PhaseSpaceResonance::Jacobian jac, double power) {
    intermediates_[ires].jacobian = jac;
    intermediates_[ires].power    = power;
  } 

public:

  /**
   *  Initialize the channel
   */
  void init(tPhaseSpaceModePtr mode);
  
  /**
   *  Initialize the channel
   */
  void initrun(tPhaseSpaceModePtr mode);

  /**
   *  The weight
   */
  const double & weight() const {return weight_;}
  
  /**
   *  Set the weight
   */
  void  weight(double in) {weight_=in;}
  
  /**
   * Reset the properties of an intermediate particle. This member is used
   * when a Decayer is used a different value for either the mass or width
   * of a resonace to that in the ParticleData object. This improves the 
   * efficiency of the integration.
   * @param part A pointer to the particle data object for the intermediate.
   * @param mass The new mass of the intermediate
   * @param width The new width of the intermediate.
   */
  void resetIntermediate(tcPDPtr part,Energy mass,Energy width) {
    if(!part) return;
    for(PhaseSpaceResonance & res : intermediates_) {
      if(res.particle!=part) continue;
      res.mass2 = sqr(mass);
      res.mWidth = mass*width;
    }
  }
  
  /**
   * Generate the momenta of the external particles. This member uses the masses
   * of the external particles generated by the DecayPhaseMode class and the
   * intermediates for the channel to generate the momenta of the decay products.
   * @param pin The momenta of the decay products. This is outputed by the member.
   * @param massext The masses of the particles. This is to allow inclusion of
   * off-shell effects for the external particles.
   */
  vector<Lorentz5Momentum> generateMomenta(const Lorentz5Momentum & pin,
					   const vector<Energy> & massext) const;
  
  
  /**
   * Generate the weight for this channel given a phase space configuration.
   * This member generates the weight for a given phase space configuration
   * and is used by the DecayPhaseSpaceMode class to compute the denominator
   * of the weight in the multi-channel integration.
   * @param output The momenta of the outgoing particles.
   */
  double generateWeight(const vector<Lorentz5Momentum> & output) const;

  /**
   * Generate the final-state particles including the intermediate resonances.
   * This method takes the outgoing particles and adds the intermediate particles
   * specified by this phase-space channel. This is to allow a given set of 
   * intermediates to be specified even when there is interference between different
   * intermediate states.
   * @param cc Whether the particles are the mode specified or its charge conjugate.
   * @param in The incoming particles.
   * @param out The outgoing particles.
   * 
   */
  void generateIntermediates(bool cc,const Particle & in, ParticleVector & out);

  /**
   *   Create a \f$2\to n\f$ diagrams for the channel
   */
  ThePEG::Ptr<ThePEG::Tree2toNDiagram>::pointer createDiagram() const;
  
public:

  /** 
   * Output operator to allow the structure to be persistently written
   * @param os The output stream
   * @param x The intermediate
   */
  inline friend PersistentOStream & operator<<(PersistentOStream & os, 
					       const PhaseSpaceChannel  & x) {
    os << x.weight_ << x.intermediates_;
    return os;
  }

  /** 
   * Input operator to allow persistently written data to be read in
   * @param is The input stream
   * @param x The NBVertex 
   */
  inline friend PersistentIStream & operator>>(PersistentIStream & is,
					       PhaseSpaceChannel & x) {
    is >> x.weight_ >> x.intermediates_;
    return is;
  }
  
  /**
   *  A friend output operator to allow the channel to be outputted for
   * debugging purposes
   */
  friend ostream & operator<<(ostream & os, const PhaseSpaceChannel & channel);

private:

  /**
   *  Find the external particles which are the children of a given resonance
   */
  void findChildren(const PhaseSpaceResonance & res,
		    vector<int> & children) {
    if(res.children.first>0)
      children.push_back(res.children.first);
    else
      findChildren(intermediates_[abs(res.children.first)],children);
    if(!res.particle) return; 
    if(res.children.second>0)
      children.push_back(res.children.second);
    else
      findChildren(intermediates_[abs(res.children.second)],children);
  }

  /**
   * Calculate the momenta for a two body decay
   * The return value indicates success or failure.
   * @param p The momentum of the decaying particle
   * @param m1 The mass of the first decay product
   * @param m2 The mass of the second decay product
   * @param p1 The momentum of the first decay product
   * @param p2 The momentum of the second decay product
   */
  void twoBodyDecay(const Lorentz5Momentum & p, 
		    const Energy m1, const Energy m2,
		    Lorentz5Momentum & p1, Lorentz5Momentum & p2) const;
    
  /** @name Mass Generation Members */
  //@{
  /**
   * Generate the mass of a resonance.
   * @param ires The resonance to be generated.
   * @param lower The lower limit on the particle's mass.
   * @param upper The upper limit on the particle's mass. 
   */
  Energy generateMass(const PhaseSpaceResonance & res,
		      Energy lower,Energy upper,
		      const double & rnd) const;
  
  /**
   * Return the weight for a given resonance.
   * @param ires The resonance to be generated.
   * @param moff The mass of the resonance.
   * @param lower The lower limit on the particle's mass.
   * @param upper The upper limit on the particle's mass. 
   */
  InvEnergy2 massWeight(const PhaseSpaceResonance & res,
			Energy moff,Energy lower,Energy upper) const;
  //@}
  
  /**
   * Helper function for the weight calculation.
   * @param ires The resonance to be generated.
   * @param limit The limit on the particle's mass. 
   */
  double atanhelper(const PhaseSpaceResonance & res, Energy limit) const;
  
private:

  /**
   *  Pointer to the phase-space mode
   */
  tPhaseSpaceModePtr mode_;
  
  /**
   *  The intermediates
   */
  vector<PhaseSpaceResonance> intermediates_;
  
  /**
   *   Integer to keep track of what we are adding
   */
  int iAdd_ = -1;
  
  /**
   *  The weight
   */
  double weight_;

};


/** 
 * Output operator to allow the structure to be persistently written
 * @param os The output stream
 * @param x The intermediate
 */
inline PersistentOStream & operator<<(PersistentOStream & os, 
				      const PhaseSpaceChannel::PhaseSpaceResonance  & x) {
  os << x.particle << ounit(x.mass2,GeV2) << ounit(x.mWidth,GeV2)
     << oenum(x.jacobian) << x.power << x.children << x.descendents;
  return os;
}
  
/** 
 * Input operator to allow persistently written data to be read in
 * @param is The input stream
 * @param x The NBVertex 
 */
inline PersistentIStream & operator>>(PersistentIStream & is,
				      PhaseSpaceChannel::PhaseSpaceResonance & x) {
  is >> x.particle >> iunit(x.mass2,GeV2) >> iunit(x.mWidth,GeV2)
     >> ienum(x.jacobian) >> x.power >> x.children >> x.descendents;
  return is;
}
  
/**
 *  A friend output operator to allow the channel to be outputted for
 * debugging purposes
 */
ostream & operator<<(ostream & os, const PhaseSpaceChannel & channel);
  
/**
 * exception for this class and those inheriting from it
 */
class PhaseSpaceError: public Exception {};
  
}

#endif /* Herwig_PhaseSpaceChannel_H */
