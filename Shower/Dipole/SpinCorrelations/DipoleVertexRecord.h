// -*- C++ -*-
#ifndef Herwig_DipoleVertexRecord_H
#define Herwig_DipoleVertexRecord_H
//
// This is the declaration of the DipoleVertexRecord class.
//

#include "ThePEG/Config/ThePEG.h"
#include "Herwig/Shower/Dipole/Base/DipoleSplittingInfo.h"
#include "Herwig/Shower/Dipole/Base/Dipole.h"
#include "ThePEG/EventRecord/RhoDMatrix.h"
#include "DipoleShowerParticle.h"

namespace Herwig {

  using namespace ThePEG;

  /**
   * Here is the documentation of the DipoleVertexRecord class.
   */
  class DipoleVertexRecord: public Base {

  public:

    /** @name Standard constructors and destructors. */
    //@{
    /**
     * The default constructor.
     */
    DipoleVertexRecord() {}

    /**
     * The destructor.
     */
    virtual ~DipoleVertexRecord() { clear(); }
    //@}

  public:

    /**
     * Prepare the emitter and spectator
     * for the spin correlations computations.
     **/
    void prepareSplitting( const DipoleSplittingInfo& dInfo, const Dipole& dip);

    /**
     * Correctly initialise the decay matrix 
     * to a delta matrix for an external particle.
     */
    void initDecayMatrix(PPtr& particle, Helicity::Direction dir);

    /**
     * Compute the spin density matrix for the given emitter.
     * This tracks the path between the given emitter and 
     * the previous emitter, calculating a rho/decay matrix
     * at each vertex as appropriate.
     */
    RhoDMatrix emitterDensityMatrix(PPtr emitter);

    /**
     * Generate the spin-correlated azimuthal angle for a splitting.
     */
    void generatePhi(DipoleSplittingInfo& dInfo, Dipole& dip);


    /**
     * Identify the type of particle and use the appropriate function
     * to set up the spin info.
     * Required for e.g. MPI
     */
    void createSpinInfo(PPtr& part,
                        const Helicity::Direction& dir);
    
    /**
     * Create and set up fermion spin info.
     * Required for e.g. MPI
     */
    void createFermionSpinInfo(PPtr& part,
                               const Helicity::Direction& dir);
      
    /**
     * Create and set up vector spin info. 
     * Required for e.g. MPI
     */
    void createVectorSpinInfo(PPtr& part,
                              const Helicity::Direction& dir);
    
    /**
     * Update the vertex record following a splitting.
     */
    void update(const DipoleSplittingInfo& dInfo);

    /**
     * For spectators. Set new particle spin info the that of the 
     * old particle. Update the spin info to include any momentum changes.
     */
    void updateSpinInfo( PPtr& oldPart,
                         PPtr& newPart );
    
    /**
     * Set the stopUpdate flag in the spin info of a particle
     * incoming to the current decay.
     */
    void prepareParticleDecay( const PPtr& parent );

    /**
     * Update the spin info of the incoming to the decay
     * following showering of the decay.
     */
    void updateParticleDecay();

    /**
     * SW 06/02/2019: Required for NearestNeighbourDipoleAnalysis tests.
     * Access the emitter info record.
     */
    //map<PPtr,DipoleSplittingInfo> emitterInfoRecord() const {
    //return theEmitterInfoRecord;
    //}    
    
    /**
     * SW 06/02/2019: Required for NearestNeighbourDipoleAnalysis tests.
     * Add a splitting to the emitter info record.
     */
    // void addToRecord(const DipoleSplittingInfo& dInfo) {
    //   assert(dInfo.emitter());
    //   theEmitterInfoRecord[dInfo.emitter()] = dInfo;
    // }
    
    /**
     * Clear the vertex record: Give up ownership
     * on any object involved in the evolution.
     */
    virtual void clear();
     
    /**
     * The standard Init function used to initialize the interfaces.
     * Called exactly once for each class by the class description system
     * before the main function starts or
     * when this class is dynamically loaded.
     */
    static void Init();
    
  private:

    /**
     * The current emitter.
     */
    DipoleShowerParticle theCurrentEmitter;

    /**
     * SW 06/02/2019: Required for NearestNeighbourDipoleAnalysis tests.
     * Record of the splittings as 
     * required for the testing analysis.
     */
    //map<PPtr, DipoleSplittingInfo> theEmitterInfoRecord;

    /**
     * The spin info of a particle incoming to the decay
     * under consideration.
     */
    tcSpinPtr theDecayParentSpinInfo;
    
    /**
     * The assignment operator is private and must never be called.
     * In fact, it should not even be implemented.
     */
    DipoleVertexRecord & operator=(const DipoleVertexRecord &) = delete;

  };

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

  /** @cond TRAITSPECIALIZATIONS */

  /** This template specialization informs ThePEG about the
   *  base classes of DipoleVertexRecord. */
  template <>
  struct BaseClassTrait<Herwig::DipoleVertexRecord,1> {
    /** Typedef of the first base class of DipoleVertexRecord. */
    typedef Base NthBase;
  };

  /** This template specialization informs ThePEG about the name of
   *  the DipoleVertexRecord class and the shared object where it is defined. */
  template <>
  struct ClassTraits<Herwig::DipoleVertexRecord>
    : public ClassTraitsBase<Herwig::DipoleVertexRecord> {
    /** Return a platform-independent class name */
    static string className() { return "Herwig::DipoleVertexRecord"; }
    /**
     * The name of a file containing the dynamic library where the class
     * DipoleVertexRecord is implemented. It may also include several, space-separated,
     * libraries if the class DipoleVertexRecord depends on other classes (base classes
     * excepted). In this case the listed libraries will be dynamically
     * linked in the order they are specified.
     */
    static string library() { return "HwDipoleShower.so"; }
  };

  /** @endcond */

}
#endif /* Herwig_DipoleVertexRecord_H */
