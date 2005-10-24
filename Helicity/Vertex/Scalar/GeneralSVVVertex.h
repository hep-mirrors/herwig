// -*- C++ -*-
#ifndef HERWIG_GeneralSVVVertex_H
#define HERWIG_GeneralSVVVertex_H
//
// This is the declaration of the GeneralSVVVertex class.
//

#include "Herwig++/Helicity/Vertex/VertexBase.h"
#include "GeneralSVVVertex.fh"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"

namespace Herwig {
  namespace Helicity {
    using namespace ThePEG;
    
    /**
     * The <code>GeneralSVVVertex<\code> class implements a
     * general Scalar-Vector-Vector vertex allowing for decay modes 
     * that only enter at the one-loop level
     * 
     * The loop integral is calculated by Passarino-Veltman reduction 
     * and the coefficients are stored here. They must be calculated
     * in the inheriting class along with implementation of the
     * setCoupling member.
     */
    class GeneralSVVVertex: public VertexBase {
      
    public:
      
      /** @name Standard constructors and destructors. */
      //@{
      /**
       * The default constructor.
       */
      inline GeneralSVVVertex();
      
      /**
       * The copy constructor.
       */
      inline GeneralSVVVertex(const GeneralSVVVertex &);
      
      /**
       * The destructor.
       */
      virtual ~GeneralSVVVertex();
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
       * Member functions to calculate helicity amplitudes for vertices
       */
      //@{
      
      /**
       * Evaluate the vertex
       * @param q2 Scale at which to evaluate the coupling
       * @param sca Scalar wavefunction 
       * @param vec1 Wavefunction of first vector particle
       * @param vec2 Wavefunction of second vector particle
       */
      Complex evaluate(Energy2 q2, const ScalarWaveFunction & sca,
		       const VectorWaveFunction & vec1,
		       const VectorWaveFunction & vec2);
      
      /**
       * Calculate coupling.
       *@param q2 Scale at which to evaluate couplings
       *@param part1 ParticleDataPointer to first particle 
       *@param part2 ParticleDataPointer to second particle
       *@param part3 ParticleDataPointer to third particle 
      */
      virtual void setCoupling(Energy2 q2,tcPDPtr part1, tcPDPtr part2,
			       tcPDPtr part3)=0;
    
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
      
       /**@name Tensor Coefficients access and setting functions. */
       //@{
      /**
       * Set tensor coefficients 
       */
      inline void a00(const Complex & val);
      
      inline void a11(const Complex & val);

      inline void a12(const Complex & val);
      
      inline void a21(const Complex & val);
      
      inline void a22(const Complex & val);

      inline void aEp(const Complex & val);

      /**
       *Access to tensor coefficients
       */
      inline Complex a00() const;
      
      inline Complex a11() const;
      
      inline Complex a12() const;
      
      inline Complex a21() const;
      
      inline Complex a22() const;
      
      inline Complex aEp() const;
      //@}

    private:
      
      /**
       * The static object used to initialize the description of this class.
       * Indicates that this is an abstract class with persistent data.
       */
      static AbstractClassDescription<GeneralSVVVertex> initGeneralSVVVertex;

      /**
       * The assignment operator is private and must never be called.
       * In fact, it should not even be implemented.
       */
      GeneralSVVVertex & operator=(const GeneralSVVVertex &);
      
      /**
       *@name Store tensor coefficients. 
       */

      Complex _a00,_a11,_a12,_a21,_a22,_aEp;
      
    };
  }
}

#include "GeneralSVVVertex.icc"

// CLASSDOC OFF

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** This template specialization informs ThePEG about the
 *  base classes of GeneralSVVVertex. */
template <>
struct BaseClassTrait<Herwig::Helicity::GeneralSVVVertex,1> {
  /** Typedef of the first base class of GeneralSVVVertex. */
  typedef Herwig::Helicity::VertexBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the GeneralSVVVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::Helicity::GeneralSVVVertex>
  : public ClassTraitsBase<Herwig::Helicity::GeneralSVVVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::Helicity::GeneralSVVVertex"; }
  /** Return the name of the shared library be loaded to get
   *  access to the GeneralSVVVertex class and every other class it uses
   *  (except the base class). */
  static string library() { return "libHwSVertex.so"; }
};

}


#endif /* HERWIG_GeneralSVVVertex_H */
