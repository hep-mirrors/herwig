// -*- C++ -*-
#ifndef HERWIG_DecayIntegrator_H
#define HERWIG_DecayIntegrator_H
//
// This is the declaration of the DecayIntegrator class.
//
#include <ThePEG/PDT/Decayer.h>
#include "DecayPhaseSpaceChannel.h"
#include <ThePEG/PDT/EnumParticles.h>
#include <Herwig++/Helicity/Correlations/DecayVertex.h>
#include "ThePEG/Utilities/Timer.h"
#include <ThePEG/Helicity/SpinInfo.h>
#include "DecayPhaseSpaceMode.fh"
#include "DecayIntegrator.fh"
#include "Herwig++/PDT/WidthCalculatorBase.h"
#include <iostream>

namespace Herwig {
using namespace ThePEG;
using Herwig::Helicity::DecayMatrixElement;

  /** \ingroup Decay
   *
   *  This class is designed to be the base class for Herwig++ decays including
   *  the implementation of a multichannel decayer or n-body phase space decays.
   *
   *  The <code>DecayIntegrator</code> class inherits from ThePEG's Decayer class
   *  and makes use of the <code>DecayPhaseSpaceMode</code> class to specify a number
   *  of decay modes.
   *
   *  Additional modes can be added using the addMode method. In practice the 
   *  phase space channels for a particular mode are usually constructed in the 
   *  doinit member of a Decayer and then the modes added to the Decayer.
   *
   *  For the majority of the decays currently implemented the 
   *  phase-space integration has been optimised and the maximum weight set.
   *  If the parameters of the decay model are changed the Initialize interface 
   *  can be used to optimise the integration and calculate the maximum weight.
   *
   *  In classes inheriting from this the me2() member which gives the matrix element
   *  squared must be implemented. This should be combined with the setting of the
   *  phase space channels, and the setting of which channels to use and their
   *  initial weights in the doinit() member. The different decay modes should then
   *  be initialized in the initrun() member if needed. The generate member can then
   *  be called from the decay() member to generate a phase-space configuration for a 
   *  decay.
   *   
   *
   * @see DecayPhaseSpaceMode
   * @see DecayPhaseSpaceChannel
   * 
   * \author Peter Richardson
   */

class DecayIntegrator: public Decayer {

  /**
   *  The output operator is a friend, this is mainly for debugging
   */    
  friend ostream & operator<<(ostream &, const DecayIntegrator &);

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline DecayIntegrator();

  /**
   * Copy-constructor.
   */
  DecayIntegrator(const DecayIntegrator &);

  /**
   * Destructor.
   */
  virtual ~DecayIntegrator();
  //@}
      
public:
  
  /**
   * Accept member which is called at initialization to see if this Decayer can
   * handle a given decay mode. As this is the base class it returns false and
   * should be overridden in class implementing the decays.
   * @param dm The DecayMode
   * @return Whether the mode can be handled.
   *
   */
  virtual bool accept(const DecayMode & dm) const;
  
  /**
   * For a given decay mode and a given particle instance, perform the
   * decay and return the decay products. As this is the base class this
   * is not implemented.
   * @param dm The DecayMode
   * @param part The Particle instant being decayed.
   * @return The vector of particles produced in the decay.
   */
  virtual ParticleVector decay(const DecayMode & dm, const Particle & part) const;
  
  /**
   * Add a phase-space mode to the list
   * @param mode The mode being added.
   * @param maxwgt The maximum weight for the phase space integration.
   * @param wgts The weights of the different channels in the multichannel approach.
   */
  void addMode(DecayPhaseSpaceModePtr mode,double maxwgt,
	       const vector<double> wgts) const;

  /**
   * Return the matrix element squared for a given mode and phase-space channel.
   * This function is purely virtual and must be implemented in classes inheriting
   * from DecayIntegrator.
   * @param vertex Output the information on the vertex for spin correlations
   * @param ichan The channel we are calculating the matrix element for. 
   * @param part The decaying Particle.
   * @param decay The particles produced in the decay.
   * @return The matrix element squared for the phase-space configuration.
   */
  virtual double me2(bool vertex, const int ichan, const Particle & part,
		      const ParticleVector & decay) const=0;
  
  /**
   * The helicity amplitude matrix element for spin correlations.
   */
  inline const DecayMatrixElement & ME() const;

  /**
   * Specify the \f$1\to2\f$ matrix element to be used in the running width calculation.
   * @param dm The DecayMode
   * @param mecode The code for the matrix element as described
   *               in the GenericWidthGenerator class.
   * @param coupling The coupling for the matrix element.
   * @return True or False if this mode can be handled.
   */
  virtual bool twoBodyMEcode(const DecayMode & dm, int & mecode,
			     double & coupling) const;
  
  /**
   * Method to return an object to calculate the 3 (or higher body) partial width
   * @param dm The DecayMode
   * @return A pointer to a WidthCalculatorBase object capable of calculating the width
   */
  virtual WidthCalculatorBasePtr threeBodyMEIntegrator(const DecayMode & dm) const;
  
  /**
   * The matrix element to be integrated for the three-body decays as a function
   * of the invariant masses of pairs of the outgoing particles.
   * @param imode The mode for which the matrix element is needed.
   * @param q2 The scale, \e i.e. the mass squared of the decaying particle.
   * @param s3 The invariant mass squared of particles 1 and 2, \f$s_3=m^2_{12}\f$.
   * @param s2 The invariant mass squared of particles 1 and 3, \f$s_2=m^2_{13}\f$.
   * @param s1 The invariant mass squared of particles 2 and 3, \f$s_1=m^2_{23}\f$.
   * @param m1 The mass of the first  outgoing particle.
   * @param m2 The mass of the second outgoing particle.
   * @param m3 The mass of the third  outgoing particle.
   * @return The matrix element
   */
  virtual double threeBodyMatrixElement(int imode,Energy2 q2, Energy2 s3,Energy2 s2,
					Energy2 s1,Energy m1,Energy m2,Energy m3);
  
  /**
   * The differential three body decay rate with one integral performed.
   * @param imode The mode for which the matrix element is needed.
   * @param q2 The scale, \e i.e. the mass squared of the decaying particle.
   * @param s  The invariant mass which still needs to be integrate over.
   * @param m1 The mass of the first  outgoing particle.
   * @param m2 The mass of the second outgoing particle.
   * @param m3 The mass of the third  outgoing particle.
   * @return The differential rate \f$\frac{d\Gamma}{ds}\f$
   */
  virtual double threeBodydGammads(int imode,Energy q2, Energy2 s,Energy m1,Energy m2,
				   Energy m3);
  
  /**
   * Set the code for the partial width. Finds the partial width in the
   * GenericWidthGenerator class which corresponds to the decay mode.
   * @param dm The DecayMode
   * @param imode The mode. 
   */
  void setPartialWidth(const DecayMode & dm, int imode);

  /**
   * Output the setup information for the particle database
   */
  virtual void dataBaseOutput(ofstream &);
  
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
   * Generate the momenta for the decay
   * @param inter Generate the intermediates produced in the decay as well as the
   * final particles.
   * @param cc Is this the mode defined or its charge conjugate.
   * @param imode The mode being generated.
   * @param inpart The decaying particle.
   * @return The particles produced inthe decay.
   */
  ParticleVector generate(bool inter,bool cc, const unsigned int & imode,
			  const Particle & inpart) const;  

  /**
   * The mode being used for this decay
   */
  inline int imode() const;

  /**
   * Set the mode being use for this decay.
   */
  inline void imode(int);
  
  /**
   * Set the helicity matrix element for the decay.
   */
  inline void ME(const DecayMatrixElement &) const;
   
  /**
   * Reset the properities of all intermediates.
   * @param part The intermediate particle being reset.
   * @param mass The mass of the particle.
   * @param width The width of the particle.
   */
  void resetIntermediate(tcPDPtr part, Energy mass, Energy width);

  /**
   * Number of decay modes
   */
  inline unsigned int numberModes() const;

  /**
   * Pointer to a mode
   */
  tDecayPhaseSpaceModePtr mode(unsigned int);
  
protected:
  
  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Check sanity of the object during the setup phase.
   */
  inline virtual void doupdate() throw(UpdateException);

  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);

  /**
   * Initialize this object to the begining of the run phase.
   */
  virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  inline virtual void dofinish();

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given pointer.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in
   * this object.
   * @return a vector of pointers.
   */
  inline virtual IVector getReferences();
  //@}

private:
  
  /**
   * Describe an abstract base class with persistent data.
   */
  static AbstractClassDescription<DecayIntegrator> initDecayIntegrator;
  
  /**
   * Private and non-existent assignment operator.
   */
  DecayIntegrator & operator=(const DecayIntegrator &);
  
protected:

  /**
   * Number of iterations for th initialization.
   */
  int _niter;

  /**
   * Number of points for initialisation
   */
  int _npoint;

  /**
   * number of attempts to generate the decay
   */
  int _ntry;

private:
  
  /**
   * List of the decay modes
   */
  mutable vector<DecayPhaseSpaceModePtr> _modes;

  /**
   * perform initialisation
   */
  bool _Initialize;


private:

  /**
   * mode currently being generated  
   */
  mutable int _imode;

  /**
   * The helicity matrix element for the current decay
   */
  DecayMatrixElement _matrixelement;
  
};
  /**
   * Output information on the DecayIntegrator for debugging purposes
   */
  ostream & operator<<(ostream &, const DecayIntegrator &);

  /**
   * Exception for this class and those inheriting from it
   */
  class DecayIntegratorError: public Exception {};

}


namespace ThePEG {
  
  /**
   * The following template specialization informs ThePEG about the
   * base class of DecayIntegrator.
   */
  template <>
  struct BaseClassTrait<Herwig::DecayIntegrator,1> {
    /** Typedef of the base class of DecayIntegrator. */
    typedef Decayer NthBase;
  };
  
  /**
   * The following template specialization informs ThePEG about the
   * name of this class and the shared object where it is defined.
   */
  template <>
  struct ClassTraits<Herwig::DecayIntegrator>
    : public ClassTraitsBase<Herwig::DecayIntegrator> {
    /** Return the class name. */
    static string className() { return "Herwig++::DecayIntegrator"; }
    /**
     * Return the name of the shared library to be loaded to get
     * access to this class and every other class it uses
     * (except the base class).
     */
    static string library() { return "libHwDecay.so"; }

  };

}

#include "DecayIntegrator.icc"

#endif /* HERWIG_DecayIntegrator_H */
