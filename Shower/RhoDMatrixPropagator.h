// -*- C++ -*-
#ifndef HERWIG_RhoDMatrixPropagator_H
#define HERWIG_RhoDMatrixPropagator_H
//
// This is the declaration of the <!id>RhoDMatrixPropagator<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This class is reponsible for the computation of the spin density matrix (rho) <BR>
// and decay matrix (D), in terms of the rho and D matrices of the particles <BR>
// (parent, and/or siblings, and/or children) connected to the same "vertex" <BR>
// (hard subprocess, or decay vertex, or splitting vertex) and the amplitude of the <BR>
// "vertex" in terms of the helicities of its connected particles. <BR>
// This class has also a switch that allows to turn OFF the rhoD matrix <BR>
// propagation completely. <BR>
// 
// Notice that:
// <UL>
//  <LI> there is an unique method <!id>matrixElement(...)<!!id> to calculate <BR> 
//       the matrix element, and an unique method <!id>computeRhoD(...)<!!id> <BR>
//       to update the rhoD matrix, independently from the kind of vertex, <BR>
//       whether hard subprocess or decay or splitting. <BR>
//       Furthermore, the computation of rhoD matrix is also independent <BR>
//       from whether the evolution is forward or backward. <BR>
//       These two not evident properties are based on the fact that <BR>
//       in <!class>ShowerParticle<!!class> we use a single rhoD matrix to store, <BR>
//       in different time, either the spin density matrix rho, or the decay matrix D. <BR>
//       This allows a very symmetric and compact writing of the various formulas.
//  <LI> the codes of both <!id>matrixElement(...)<!!id> and <!id>computeRhoD(...)<!!id> <BR>
//       do not have to separate between the various vertex multiplicity cases: <BR>
//       for hard subprocesses:    <I> 2 -&GT; 2 , 2 -&GT; 3 , ... , 2 -&GT; N </I> <BR>
//       for decays or splittings: <I> 1 -&GT; 2 , 1 -&GT; 3 , ... , 1 -&GT; N </I> <BR>
//       The general case of generic <I> N &GT;= 2 </I> is treated directly. <BR>
//       This is possible by avoiding to use explicitly nested for loops <BR> 
//       (which, of course, would require to know exactly how many of them we have <BR>
//        to write down) and generating instead automatically all the needed <BR>
//       helicity configurations. However, only in the case of splitting processes 
//       <I> 1-&GT;N </I>, <BR>
//       you have to distinguish between the various <I>N</I> in the implementation <BR>
//       (but not in the interface) of the method <!id>evaluateAmplitudes(...)<!!id>, <BR>
//       using downcasting. 
// </UL>
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:ShowerParticle.html">ShowerParticle.h</a>.
// 

#include "Pythia7/Handlers/HandlerBase.h"
#include "ShowerConfig.h"
#include "Herwig++/Config/GlobalParameters.h"
#include "Pythia7/MatrixElement/MEBase.h"


namespace Herwig {

using namespace Pythia7;

class Pythia7::PartialCollisionHandler;  // forward declaration
class Pythia7::Decayer;                  // forward declaration


class RhoDMatrixPropagator: public Pythia7::HandlerBase {

public:

  inline RhoDMatrixPropagator();
  inline RhoDMatrixPropagator(const RhoDMatrixPropagator &);
  virtual ~RhoDMatrixPropagator();
  // Standard ctors and dtor.

  inline bool isRhoDPropagationON() const;
  // This method returns true (false) if the rhoD propagation 
  // is switched ON (OFF).

  double matrixElement( const tMEPtr hardSubME, const CollecShoParPtr & collecShoPar,
		        const tShoParPtr aParticle ); 
  // Given <!id>aParticle<!!id> coming from a given vertex <BR>
  // --- this must be the decaying particle in the case of a decay vertex; 
  // or the emitting particle in the case of a splitting vertex; or any of the 
  // particles entering the hard subprocess, in the case of the hard subprocess vertex. 
  // Notice that only in the latter case the <!id>hardSubME<!!id> is set 
  // properly (not null) and used to access the hard subprocess matrix element, 
  // and only in this case the collections of shower particles is needed in order
  // to find the ones connected with the hard process --- <BR>
  // it computes the real matrix element to be used to generate such vertex. 
  // If something goes wrong, the method returns 0.
  // In the case the <!id>onoffSwitchMode<!!id> is <I>0 (OFF)</I>, the method 
  // does nothing and returns <I>1.0</I>.
    
  bool computeRhoD( const tMEPtr hardSubME, const CollecShoParPtr & collecShoPar,
		    const tShoParPtr theParticle ); 
  // It computes the rhoD matrix of <!id>theParticle<!!id> in terms of the rhoD 
  // matrices of the other particles entering the same vertex, and the
  // amplitude (and its conjugate) of the vertex. The "vertex" can be
  // a hard subprocess one --- in which case, and only in this one, the 
  // pointer to the hard subprocess matrix element <!id>hardSubME<!!id>, 
  // is properly set (not null) and also only in this case the collections 
  // of shower particles is needed in order to find the ones connected with 
  // the hard subprocess --- 
  // or a decay one, or a splitting one.
  // The method returns true if it succeed, false otherwise.
  // In the case the <!id>onoffSwitchMode<!!id> is <I>0 (OFF)</I>, 
  // the method does nothing and returns true.

public:

  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.

  static void Init();
  // Standard Init function used to initialize the interfaces.

protected:

  inline virtual IBPtr clone() const;
  inline virtual IBPtr fullclone() const;
  // Standard clone methods.

protected:

  inline virtual void doupdate() throw(UpdateException);
  inline virtual void doinit() throw(InitException);
  inline virtual void dofinish();
  // Standard Interfaced virtual functions.

  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);
  // Change all pointers to Interfaced objects to corresponding clones.

  inline virtual IVector getReferences();
  // Return pointers to all Interfaced objects refered to by this.

private:

  static ClassDescription<RhoDMatrixPropagator> initRhoDMatrixPropagator;
  // Describe a concrete class with persistent data.

  RhoDMatrixPropagator & operator=(const RhoDMatrixPropagator &);
  // Private and non-existent assignment operator.

  Complex trace( const ComplexMatrix & matrix ) const;
  // Returns the trace of the matrix.

  bool normalize( ComplexMatrix & matrix );
  // If the trace of the matrix is not zero then it normalizes the matrix
  // (that is the new matrix has trace equal <I>1</I>) and returns true; 
  // otherwise it returns false.

  void evaluateAmplitude( const tShoParPtr particle, 
			  const tcPDVector & dataParticles, 
			  const vector<Lorentz5Momentum> & momenta, 
			  const vector<int> & helicities, const vector<int> & helicitiesPrime, 
			  Complex & amplitudeValue, Complex & amplitudePrimeValue ); 
  // We can treat similarly hard <I>2-&GT;N</I> processes
  // and decay <I>1-&GT;N</I> processes, because both depend on momenta 
  // and helicities, but the case of splitting functions must instead be dealt 
  // differently, because the dependency is, in the common <I>1-&GT;2</I> case, 
  // on <I>(z,phi)</I> and helicities rather than in momenta and helicities.
  // This method receives in input the following parameters:
  // --- the pointer to the emitting particle, in the case of
  //     a splitting process (whereas it is ignored in the other
  //     cases, hard processes and decay processes);
  // --- a vector of data particles for all particles attached to the "vertex" 
  //     (whatever it is: hard process, decay process, splitting process);
  // --- momenta of all the particles attached to the "vertex"
  // --- the two helicities configurations, helicities and helicitiesPrime;
  // and it returns the values of the amplitudes associated to the vertex
  // (hard process, decay process, splitting process) in correspondence
  // of the two helicities configurations.
  // Notice that, whereas for hard <I>2-&GT;N</I> processes and 
  // decay <I>1-&GT;N</I> processes the procedure holds for any 
  // final state multiplicity <I>N</I>, in the case of 
  // splitting <I>1-&GT;N</I> process the interface of this method remains
  // unchanged, but in the implementation of this method you have
  // to esplicitly distinguish between difference cases, using downcasting
  // (see comment in the implementation code).

  bool nextIndexConfiguration( const vector<int> & sizes,
			       vector<int> & indeces, vector<int> & indecesPrime );
  // This method determine the next index configuration, given in input
  // the vector (sizes) of max values for each index, and the current 
  // index configuration (vectors: indeces and indecesPrime). The returned
  // new configuration is overwritten to the current one.
  // The method returns false if there is not anymore other index 
  // configurations to be considered; true otherwise.

  double matrixElement( const CollecShoParPtr & particles ); 
  // Given all of the <!id>particles<!!id> coming from the same vertex, it computes
  // the real matrix element to be used to generate such vertex.
  // If something goes wrong, the method returns <I>0</I>.

  bool computeRhoD( const tShoParPtr theParticle, const CollecShoParPtr & particles ); 
  // Given all of the <!id>particles<!!id> coming from the same vertex, it computes
  // the rhoD matrix of one of them, <!id>theParticle<!!id> . 
  // It returns true if it succeed, false otherwise.

  bool setVertexPointer( const tMEPtr hardSubME, const tShoParPtr particlePtr ); 
  // It sets the hard subprocess matrix element pointer, or the decayer pointer, 
  // or the splitFun pointer depending if the vertex to which the <!id>particlePtr<!!id> 
  // is connected to is respectively a hard subprocess vertex, a decay vertex, 
  // a splitting vertex. Notice that in the former case, the input particle, 
  // <!id>particlePtr<!!id>, can be any of the incoming or outgoing particles 
  // in the hard subprocess, whereas in the latter two cases the particle must be 
  // the decaying one or the emitting (splitting) one. 
  // Notice that in the case of an outgoing particle from the hard
  // subprocess, such particle could also be, at the same time, a decaying 
  // or emitting (splitting) particle, therefore there could be in 
  // principle an ambiguity about which vertex to consider for it.
  // We use the convention that <!id>hardSubME<!!id> is passed not null
  // only when we want to consider the hard process as vertex.
  // The method returns false if it does not succeed, true otherwise.

  void findVertexParticles( const CollecShoParPtr & collecShoPar , const tShoParPtr aParticle,
			    CollecShoParPtr & particles );
  // This method puts in the vector <!id>particles<!!id> all the particles 
  // (indeed pointers to ShowerParticles objects) that belong to
  // the same "vertex" (hard subprocess, or decay, or splitting) as
  // the the input particle <!id>aParticle<!!id>. It needs the collection of
  // all shower particles only in the case that the vertex is the
  // hard subprocess (because the subprocess object would give only
  // the incoming and outgoing Pythia7 particles that enter the
  // hard subprocess, but not the ShowerParticle objects we want.

  // In order to access the amplitude of the "vertex", as a function of 
  // the helicities of the particles connected to that vertex, we need
  // a pointer to either the hard subprocess matrix element, or the decayer, 
  // or the splitFun. Notice that only one of them should be not null.
  tMEPtr       _hardSubME;
  tDecayerPtr  _decayer;
  tSplitFunPtr _splitFun;

  int _onoffSwitchMode;

};

}

// CLASSDOC OFF

namespace Pythia7 {

// The following template specialization informs Pythia7 about the
// base class of RhoDMatrixPropagator.
template <>
struct BaseClassTrait<Herwig::RhoDMatrixPropagator,1> {
  typedef Pythia7::HandlerBase NthBase;
};

// The following template specialization informs Pythia7 about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::RhoDMatrixPropagator>: public ClassTraitsBase<Herwig::RhoDMatrixPropagator> {
  static string className() { return "/Herwig++/RhoDMatrixPropagator"; }
  // Return the class name.
  static string library() { return "libHwShower.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "RhoDMatrixPropagator.icc"

#endif /* HERWIG_RhoDMatrixPropagator_H */
