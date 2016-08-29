// -*- C++ -*-

// couplings/lo_alpha_s.h is part of matchbox
// (C) 2008 Simon Platzer -- sp@particle.uni-karlsruhe.de

#ifndef matchbox_couplings_lo_alpha_s_h
#define matchbox_couplings_lo_alpha_s_h

#include "alpha_s.h"

namespace matchbox {

  using namespace ThePEG;

  /**
   * LO running alpha_s
   *
   * @see \ref lo_alpha_sInterfaces "The interfaces"
   * defined for lo_alpha_s.
   */
  class lo_alpha_s
    : public alpha_s {

  public:

    /** @name Standard constructors and destructors. */
    //@{
    /**
     * The default constructor.
     */
    lo_alpha_s();

    /**
     * The destructor.
     */
    virtual ~lo_alpha_s();
    //@}

  public:

    /// return alpha_s as function of scale, QCD scale
    /// and number of active flavours
    virtual double operator () (Energy2 scale,
				Energy2 lambda2,
				unsigned int nf) const;

    /// return the number of loops which determine this running
    virtual unsigned int nloops () const { return 1; }

  public:

    /** @name Functions used by the persistent I/O system. */
    //@{
    /**
     * Function used to write out object persistently.
     * @name os the persistent output stream written to.
     */
    void persistentOutput(PersistentOStream & os) const;

    /**
     * Function used to read in object persistently.
     * @name is the persistent input stream read from.
     * @name version the version number of the object when written.
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

    /** @name Standard Interfaced functions. */
    //@{

    /**
     * Initialize this object after the setup phase before saving an
     * EventGenerator to disk.
     * @throws InitException if object could not be initialized properly.
     */
    virtual inline void doinit() throw(InitException) {
      freezing_scale_ *= scale_factor();
      alpha_s::doinit();
    }

    //@}

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

  private:

    /**
     * The static object used to initialize the description of this class.
     * Indicates that this is an abstract class with persistent data.
     */
    static ClassDescription<lo_alpha_s> initlo_alpha_s;

    /**
     * The assignment operator is private and must never be called.
     * In fact, it should not even be implemented.
     */
    lo_alpha_s & operator=(const lo_alpha_s &);

  private:

    Energy freezing_scale_;

  };

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

  /** @cond TRAITSPECIALIZATIONS */

  /** This template specialization informs ThePEG about the
   *  base classes of lo_alpha_s. */
  template <>
  struct BaseClassTrait<matchbox::lo_alpha_s,1> {
    /** Typedef of the first base class of lo_alpha_s. */
    typedef matchbox::alpha_s NthBase;
  };

  /** This template specialization informs ThePEG about the name of
   *  the lo_alpha_s class and the shared object where it is defined. */
  template <>
  struct ClassTraits<matchbox::lo_alpha_s>
    : public ClassTraitsBase<matchbox::lo_alpha_s> {
    /** Return a platform-independent class name */
    static string className() { return "matchbox::lo_alpha_s"; }
    /**
     * The name of a file containing the dynamic library where the class
     * lo_alpha_s is implemented. It may also include several, space-separated,
     * libraries if the class lo_alpha_s depends on other classes (base classes
     * excepted). In this case the listed libraries will be dynamically
     * linked in the order they are specified.
     */
    static string library() { return "HwDipoleShowerAlphaS.so"; }
  };

  /** @endcond */

}

#endif /* matchbox_couplings_lo_alpha_s_h */
