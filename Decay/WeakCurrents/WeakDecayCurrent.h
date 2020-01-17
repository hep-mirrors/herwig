// -*- C++ -*-
//
// WeakDecayCurrent.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_WeakDecayCurrent_H
#define HERWIG_WeakDecayCurrent_H
//
// This is the declaration of the WeakDecayCurrent class.
//
#include "ThePEG/Interface/Interfaced.h"
#include "WeakDecayCurrent.fh"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"
#include "Herwig/Decay/DecayPhaseSpaceChannel.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"

namespace Herwig {
using namespace ThePEG;
using ThePEG::Helicity::LorentzPolarizationVector;
using ThePEG::Helicity::LorentzPolarizationVectorE;

/** \ingroup Decay
 *
 *  The <code>WeakDecayCurrent</code> class is the base class for the hadronic
 *  currents produced in weak decays. This is designed so it can be used in
 *  any weak decay. In general the currents which are implemented are either
 *  simple one meson production or taken from tau decay.
 *
 *  In classes inheriting from this one a number of member functions must be implemented
 *
 * - createMode which takes a vector of partially completed DecayPhaseSpaceChannel 
 *   and adds the extra information required for the current. This method should
 *   assume that the particles from the current are the last ones specified in the
 *   DecayPhaseSpaceMode. This method will then construct the DecayPhaseSpaceMode
 *   for the decay.
 *
 * - particles() which returns the external particles produced by the current.
 *
 * - current() which given the decay products calculates the decay current
 *
 * - accept() which uses the PDG codes for the particles in the current to 
 *   decide if a given mode is allowed.
 *
 * - decayMode() which uses the PDG codes for the particles in the current to 
 *   workout which modes is being performed.
 *
 * - dataBaseOutput() which should output the information on all the Interfaces so
 *   that the WeakDecayCurrent can be reconstructed by the Herwig particle
 *   properties database.
 *
 * @see Interfaced.
 * 
 */
class WeakDecayCurrent: public Interfaced {

public:

  /**
   * Default constructor
   */
  WeakDecayCurrent() : _numbermodes(0) {}

public:

  /** @name Methods for the construction of the phase space integrator. */
  //@{
  /**
   * Complete the construction of the decay mode for integration.classes inheriting
   * from this one.
   * This method is purely virtual and must be implemented in the classes inheriting
   * from WeakDecayCurrent.
   * @param icharge The total charge of the outgoing particles in the current.
   * @param imode   The mode in the current being asked for.
   * @param mode    The phase space mode for the integration
   * @param iloc    The location of the of the first particle from the current in
   *                the list of outgoing particles.
   * @param ires    The location of the first intermediate for the current.
   * @param phase   The prototype phase space channel for the integration.
   * @param upp     The maximum possible mass the particles in the current are
   *                allowed to have.
   * @return Whether the current was sucessfully constructed.
   */
  virtual bool createMode(int icharge,unsigned int imode,DecayPhaseSpaceModePtr mode,
			  unsigned int iloc,unsigned int ires,
			  DecayPhaseSpaceChannelPtr phase,Energy upp)=0;

  /**
   * The particles produced by the current. This method is purely virtual and
   * must be implemented in classes inheriting from this one.
   * @param icharge The total charge of the particles in the current.
   * @param imode The mode for which the particles are being requested
   * @param iq The PDG code for the quark
   * @param ia The PDG code for the antiquark
   * @return The external particles for the current.
   */
  virtual tPDVector particles(int icharge, unsigned int imode, int iq, int ia)=0;
  //@}

  /**
   * Return the number of modes handled by this current
   */
  unsigned int numberOfModes() const {return _quark.size();}

  /**
   * Hadronic current. This method is purely virtual and must be implemented in
   * all classes inheriting from this one.
   * @param imode The mode
   * @param ichan The phase-space channel the current is needed for.
   * @param scale The invariant mass of the particles in the current.
   * @param decay The decay products
   * @param meopt Option for the calculation of the matrix element
   * @return The current. 
   */
  virtual vector<LorentzPolarizationVectorE> 
  current(const int imode, const int ichan,Energy & scale, 
	  const ParticleVector & decay, DecayIntegrator::MEOption meopt) const=0;

  /**
   * Accept the decay. This method is purely virtual and must be implemented in any class
   * inheriting from this one.
   * @param id The id's of the particles in the current.
   * @return Can this current have the external particles specified.
   */
  virtual bool accept(vector<int> id)=0;

  /**
   * Return the decay mode number for a given set of particles in the current. 
   * This method is purely virtual and must be implemented in any class
   * inheriting from this one.
   * @param id The id's of the particles in the current.
   * @return The number of the mode
   */
  virtual unsigned int decayMode(vector<int> id)=0;

  /**
   *  Information on a decay mode
   * @param imode The mode
   * @param iq The PDG code of the quark.
   * @param ia The PDG code of the antiquark.
   */
  void decayModeInfo(unsigned int imode, int& iq, int& ia) const {
    if(imode<_quark.size()) {
      iq=_quark[imode];
      ia=_antiquark[imode];
    }
    else {
      iq=0;
      ia=0;
    }
  }

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   * @param create Whether or not to add a statement creating the object
   */
  virtual void dataBaseOutput(ofstream & os,bool header,bool create) const;

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
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

protected:

  /**
   *  Add a decay mode to the list.
   * @param iq The PDG code for the quark.
   * @param ia The PDG code for the antiquark.
   */
  void addDecayMode(int iq,int ia) {
    _quark.push_back(iq);
    _antiquark.push_back(ia);
  }

  /**
   *  Set initial number of modes
   * @param nmodes The number of modes.
   */
  void setInitialModes(unsigned int nmodes) {_numbermodes=nmodes;}

private:

  /**
   * Private and non-existent assignment operator.
   */
  WeakDecayCurrent & operator=(const WeakDecayCurrent &) = delete;

private:

  /**
   * The PDG codes for the quarks contained in the current.
   */
  vector<int> _quark;

  /**
   * The PDG codes for the antiquarks contained in the current.
   */
  vector<int> _antiquark;

  /**
   * The initial number of modes
   */
  unsigned int _numbermodes;

};

}


#endif /* HERWIG_WeakDecayCurrent_H */
