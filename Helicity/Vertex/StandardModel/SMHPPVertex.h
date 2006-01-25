// -*- C++ -*-
#ifndef HERWIG_SMHPPVertex_H
#define HERWIG_SMHPPVertex_H
//
// This is the declaration of the SMHPPVertex class.
//

#include "Herwig++/Helicity/Vertex/Scalar/SVVLoopVertex.h"
#include "SMHPPVertex.fh"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"

namespace Herwig {
  namespace Helicity {
    using namespace ThePEG;
    
    /**
     * The <code>SMHPPVertex</code> class implements the
     * setCoupling member for the Standard Model Higgs to  
     * gamma,gamma decay mode.
     */
    class SMHPPVertex: public SVVLoopVertex {
    
    public:
    
      /** @name Standard constructors and destructors. */
      //@{
      /**
       * The default constructor.
       */
      inline SMHPPVertex();
    
      /**
       * The copy constructor.
       */
      inline SMHPPVertex(const SMHPPVertex &);

      /**
       * The destructor.
       */
      virtual ~SMHPPVertex();
      //@}

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

      /**
       * Calculate couplings
       *@param q2 Scale at which to evaluate coupling
       *@param part1 ParticleData pointer to first particle
       *@param part2 ParticleData pointer to first particle
       *@param part3 ParticleData pointer to first particle
       */
      virtual void setCoupling(Energy q2, tcPDPtr part1, tcPDPtr part2,
			       tcPDPtr part3);
    
    protected:

      /** @name Clone Methods. */
      //@{
      /**
       * Make a simple clone of this object.
       * @return a pointer to the new object.
       */
      inline virtual IBPtr clone() const;

      /** Make a clone of this object, possibly modifying the cloned object
       * to make it sane.
       * @return a pointer to the new object.
       */
      inline virtual IBPtr fullclone() const;
      //@}

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
       * Initialize this object. Called in the run phase just before
       * a run begins.
       */
      inline virtual void doinitrun();

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
       * @throws RebindException if no cloned object was found for a given
       * pointer.
       */
      inline virtual void rebind(const TranslationMap & trans)
	throw(RebindException);

      /**
       * Return a vector of all pointers to Interfaced objects used in this
       * object.
       * @return a vector of pointers.
       */
      inline virtual IVector getReferences();
      //@}

    private:

      /**
       * The static object used to initialize the description of this class.
       * Indicates that this is a concrete class with persistent data.
       */
      static ClassDescription<SMHPPVertex> initSMHPPVertex;

      /**
       * The assignment operator is private and must never be called.
       * In fact, it should not even be implemented.
       */
      SMHPPVertex & operator=(const SMHPPVertex &);

      /**
       *Storage of couplings
       */
      //@{
      /**
       *Last value of the coupling calculated
       */
      Complex _couplast;
    
      /**
       
      *The scale \f$q^2\f$ at which coupling was last evaluated
      */
      Energy2 _q2last;
      //@}
    
      /**
       *Pointer to Standard Model object
       */
      Ptr<Herwig::StandardModel>::pointer _theSM;

      /**
       *Mass of W boson for higgs coupling
       */
      Energy _mw;
    
      /**
       *Storage of \f$\sin\theta_W\f$
       */
      double _sw;

    };
  }
}

// CLASSDOC OFF

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** This template specialization informs ThePEG about the
 *  base classes of SMHPPVertex. */
template <>
struct BaseClassTrait<Herwig::Helicity::SMHPPVertex,1> {
  /** Typedef of the first base class of SMHPPVertex. */
  typedef Herwig::Helicity::SVVLoopVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SMHPPVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::Helicity::SMHPPVertex>
  : public ClassTraitsBase<Herwig::Helicity::SMHPPVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::SMHPPVertex"; }
  /** Return the name of the shared library be loaded to get
   *  access to the SMHPPVertex class and every other class it uses
   *  (except the base class). */
  static string library() { return "HwSMVertex.so"; }
};

}

#include "SMHPPVertex.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "SMHPPVertex.tcc"
#endif

#endif /* HERWIG_SMHPPVertex_H */
