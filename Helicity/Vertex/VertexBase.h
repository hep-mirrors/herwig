// -*- C++ -*-
#ifndef HERWIG_VertexBase_H
#define HERWIG_VertexBase_H
//
// This is the declaration of the VertexBase class.

#include <ThePEG/Config/Complex.h>
#include "ThePEG/Interface/Interfaced.h"
#include <ThePEG/PDT/ParticleData.h>
#include <ThePEG/PDT/WidthGenerator.h>
#include <ThePEG/Helicity/HelicityDefinitions.h>

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

/** \ingroup Helicity
 * 
 *  The VertexBase class is the base class for all helicity amplitude
 *  vertices inside Herwig++. In implements the storage of the particles 
 *  which are allowed to interact at the vertex and some simple functions 
 *  which are often needed by the classes which implement the specific 
 *  vertices.
 *
 *  In practice little use is made of this information and it is mainly
 *  included for future extensions. It can also be used at the development
 *  and debugging stage.
 *
 */
class VertexBase  : public Interfaced {
  
/**
 *  The output operator is a friend to avoid the data being public.
 */
friend ostream & operator<<(ostream &, const VertexBase &);

public:
  
  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline VertexBase();

  /**
   * Copy-constructor.
   */
  inline VertexBase(const VertexBase &);

  /**
   * Constructor for three point vertices.
   * @param ispin1 \f$2S+1\f$ for the first  particle.
   * @param ispin2 \f$2S+1\f$ for the second particle.
   * @param ispin3 \f$2S+1\f$ for the third  particle.
   * @param id1 The PDG codes for the first  set of particles.
   * @param id2 The PDG codes for the second set of particles.
   * @param id3 The PDG codes for the third  set of particles.
   * @param kine Whether the kinematic invariants should be calculated.
   */
  inline VertexBase(int ispin1,int ispin2,int ispin3,
		    vector<int> id1,vector<int> id2,vector<int> id3,bool kine=false);

  /**
   * Constructor for four point vertices.
   * @param ispin1 \f$2S+1\f$ for the first  particle.
   * @param ispin2 \f$2S+1\f$ for the second particle.
   * @param ispin3 \f$2S+1\f$ for the third  particle.
   * @param ispin4 \f$2S+1\f$ for the fourth  particle.
   * @param id1 The PDG codes for the first  set of particles.
   * @param id2 The PDG codes for the second set of particles.
   * @param id3 The PDG codes for the third  set of particles.
   * @param id4 The PDG codes for the fourth  set of particles.
   * @param kine Whether the kinematic invariants should be calculated.
   */
  inline VertexBase(int ispin1,int ispin2,int ispin3,int ispin4,
		    vector<int> id1,vector<int> id2,vector<int> id3,vector<int> id4,
		    bool kine=false);

  /**
   * Constructor for five point vertices.
   * @param ispin1 \f$2S+1\f$ for the first  particle.
   * @param ispin2 \f$2S+1\f$ for the second particle.
   * @param ispin3 \f$2S+1\f$ for the third  particle.
   * @param ispin4 \f$2S+1\f$ for the fourth particle.
   * @param ispin5 \f$2S+1\f$ for the fifth  particle.
   * @param id1 The PDG codes for the first  set of particles.
   * @param id2 The PDG codes for the second set of particles.
   * @param id3 The PDG codes for the third  set of particles.
   * @param id4 The PDG codes for the fourth set of particles.
   * @param id5 The PDG codes for the fifth  set of particles.
   * @param kine Whether the kinematic invariants should be calculated.
   */
  inline VertexBase(int ispin1,int ispin2,int ispin3,int ispin4,int ispin5,
		    vector<int> id1,vector<int> id2,vector<int> id3,vector<int> id4,
		    vector<int> id5,bool kine=false);

  /**
   * Constructor for \f$n\f$-point vertices.
   * @param npoint The number of external particles.
   * @param kine Whether the kinematic invariants should be calculated.
   */
  inline VertexBase(int npoint,bool kine=false);

  /**
   * Destructor.
   */
  virtual ~VertexBase();
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
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();
  
public:

  /**
   * Add particles to the lists.
   */
  //@{ 
  /**
   * Add item to three point list.
   * @param id1 PDG code of the first particle.
   * @param id2 PDG code of the second particle.
   * @param id3 PDG code of the third particle.
   */
  void add(int id1,int id2,int id3);

  /**
   * Add item to four point list.
   * @param id1 PDG code of the first particle.
   * @param id2 PDG code of the second particle.
   * @param id3 PDG code of the third particle.
   * @param id4 PDG code of the fourth particle.
   */
  void add(int id1,int id2,int id3,int id4);

  /**
   * Add item to five  point list.
   * @param id1 PDG code of the first particle.
   * @param id2 PDG code of the second particle.
   * @param id3 PDG code of the third particle.
   * @param id4 PDG code of the fourth particle.
   * @param id5 PDG code of the fifth particle.
   */
  void add(int id1,int id2,int id3,int id4,int id5);
  //@}

  /**
   *  Access to the particle information
   */
  //@{
  /**
   * Number of different particle combinations allowed.
   */
  inline unsigned int size();

  /**
   * Get the number of external particles.
   */
  inline int getNpoint();

  /**
   * Is a particle allowed as an incoming particle?
   * @param id The PDG code
   */
  inline bool incoming(int id);

  /**
   * Is a particle allowed as an outgoing particle?
   * @param id The PDG code
   */
  inline bool outgoing(int id);

  /**
   * Get the list of incoming particles.
   */
  inline vector<PDPtr> getIncoming();

  /**
   * Get the list of outgoing particles.
   */
  inline vector<PDPtr> getOutgoing();

  /**
   * Get the coupling.
   */
  inline const Complex & getNorm();

  /**
   * Function to search the list.
   * @param ilist Which list to search
   * @param id The PDG code to look for.
   */
  vector<PDPtr> search(int ilist,int id);

  /**
   * Is a given combination allowed.
   * @param id1 PDG code of the first particle.
   * @param id2 PDG code of the second particle.
   * @param id3 PDG code of the third particle.
   */
  bool allowed(int id1,int id2,int id3);

  /**
   * Is a given combination allowed.
   * @param id1 PDG code of the first particle.
   * @param id2 PDG code of the second particle.
   * @param id3 PDG code of the third particle.
   * @param id4 PDG code of the fourth particle.
   */
  bool allowed(int id1,int id2,int id3,int id4);

  /**
   * Is a given combination allowed.
   * @param id1 PDG code of the first particle.
   * @param id2 PDG code of the second particle.
   * @param id3 PDG code of the third particle.
   * @param id4 PDG code of the fourth particle.
   * @param id5 PDG code of the fifth particle.
   */
  bool allowed(int id1,int id2,int id3,int id4,int id5);
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
   * Initialize this object to the begining of the run phase.
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

protected:

  /**
   *  Members to set-up the particles
   */
  //@{
  /**
   * Setup the spins of the particles for a three point vertex.
   * @param ispin1 \f$2S+1\f$ for the first  particle.
   * @param ispin2 \f$2S+1\f$ for the second particle.
   * @param ispin3 \f$2S+1\f$ for the third  particle.
   */
  inline void setSpin(int ispin1,int ispin2,int ispin3);

  /**
   * Setup the spins of the particles for a four point vertex.
   * @param ispin1 \f$2S+1\f$ for the first  particle.
   * @param ispin2 \f$2S+1\f$ for the second particle.
   * @param ispin3 \f$2S+1\f$ for the third  particle.
   * @param ispin4 \f$2S+1\f$ for the fourth particle.
   */
  inline void setSpin(int ispin1,int ispin2,int ispin3,int ispin4);

  /**
   * Setup the spins of the particles for a five point vertex.
   * @param ispin1 \f$2S+1\f$ for the first  particle.
   * @param ispin2 \f$2S+1\f$ for the second particle.
   * @param ispin3 \f$2S+1\f$ for the third  particle.
   * @param ispin4 \f$2S+1\f$ for the fourth particle.
   * @param ispin5 \f$2S+1\f$ for the fifth  particle.
   */
  inline void setSpin(int ispin1,int ispin2,int ispin3,int ispin4,int ispin5);

  /**
   * Set up the lists of particles for the three point vertex.
   * @param id1 The PDG codes for the first  set of particles.
   * @param id2 The PDG codes for the second set of particles.
   * @param id3 The PDG codes for the third  set of particles.
   */
  inline void setList(vector<int> id1,vector<int> id2,vector<int> id3);

  /**
   * Set up the lists of particles for the three point vertex.
   * @param id1 The PDG codes for the first  set of particles.
   * @param id2 The PDG codes for the second set of particles.
   * @param id3 The PDG codes for the third  set of particles.
   * @param id4 The PDG codes for the fourth set of particles.
   */
  inline void setList(vector<int> id1,vector<int> id2,vector<int> id3,vector<int> id4);

  /**
   * Set up the lists of particles for the three point vertex.
   * @param id1 The PDG codes for the first  set of particles.
   * @param id2 The PDG codes for the second set of particles.
   * @param id3 The PDG codes for the third  set of particles.
   * @param id4 The PDG codes for the fourth set of particles.
   * @param id5 The PDG codes for the fifth  set of particles.
   */
  inline void setList(vector<int> id1,vector<int> id2,vector<int> id3,vector<int> id4,
		      vector<int> id5);

  /**
   * Set the list of incoming particles.
   */
  void setIncoming();

  /**
   * Set the list of outgoing particles.
   */
  void setOutgoing();

  /**
   * Set the number of external particles.
   * @param npoint The number of external particles.
   */
  inline void setNpoint(int npoint);
  //@}

  /**
   *  Members for the amplitude calculations
   */
  //@{
  /**
   * Set the coupling.
   * @param coup The coupling.
   */
  inline void setNorm(const Complex & coup);

  /**
   * Calculate the propagator for a diagram.
   * @param iopt The option for the Breit-Wigner shape
   * @param q2 The scale
   * @param part The ParticleData pointer for the off-shell particle.
   */
  inline Complex propagator(int iopt, Energy2 q2,tcPDPtr part);

  /**
   * Calculate propagator multiplied by coupling.
   * @param iopt The option for the Breit-Wigner shape
   * @param q2 The scale
   * @param part The ParticleData pointer for the off-shell particle.
   */
  inline Complex normPropagator(int iopt, Energy2 q2,tcPDPtr part);
  //@}    

  /** @name Kinematic invariants for loop diagrams */
  //@{

  /**
   * Whether or not to calculate the kinematics invariants
   */
  inline bool kinematics();

  /**
   * Set whether or not to calculate the kinematics invariants
   */
  inline void kinematics(bool );

  /**
   *  Calculate the kinematics for a 3-point vertex
   */
  inline void calculateKinematics(const Lorentz5Momentum &,const Lorentz5Momentum &,
				  const Lorentz5Momentum &);
  /**
   *  Calculate the kinematics for a 4-point vertex
   */
  inline void calculateKinematics(const Lorentz5Momentum &,const Lorentz5Momentum &,
				  const Lorentz5Momentum &,const Lorentz5Momentum &);
  /**
   *  Calculate the kinematics for a 5-point vertex
   */
  inline void calculateKinematics(const Lorentz5Momentum &,const Lorentz5Momentum &,
				  const Lorentz5Momentum &,const Lorentz5Momentum &,
				  const Lorentz5Momentum &);
  /**
   *  Calculate the kinematics for a n-point vertex
   */
  inline void calculateKinematics(const vector<Lorentz5Momentum> &);

  /**
   * Get one of the kinematic invariants
   */
  inline Energy2 invariant(unsigned int,unsigned int);
  //@}
private:
  
  /**
   * Describe a concrete class with persistent data.
   */
  static AbstractClassDescription<Herwig::Helicity::VertexBase> initVertexBase;
  
  /**
   * Private and non-existent assignment operator.
   */
  VertexBase & operator=(const VertexBase &);
  
private:

  /**
   * Storage of the particles.   
   */
  //@{
  /**
   * PDG codes for the first set of  particles.
   */
  vector<int> _iparticlea;

  /**
   * PDG codes for the second set of  particles.
   */
  vector<int> _iparticleb;

  /**
   * PDG codes for the third set of  particles.
   */
  vector<int> _iparticlec;

  /**
   * PDG codes for the fourth set of  particles.
   */
  vector<int> _iparticled;

  /**
   * PDG codes for the fifth set of  particles.
   */
  vector<int> _iparticlee;

  /**
   * Particle data pointers for the first set of particles.
   */
  vector<PDPtr> _particlea;

  /**
   * Particle data pointers for the second set of particles.
   */
  vector<PDPtr> _particleb;

  /**
   * Particle data pointers for the third set of particles.
   */
  vector<PDPtr> _particlec;

  /**
   * Particle data pointers for the fourth set of particles.
   */
  vector<PDPtr> _particled;

  /**
   * Particle data pointers for the fifth set of particles.
   */
  vector<PDPtr> _particlee;

  /**
   * Spin.
   */
  vector<int> _ispin;

  /**
   *  Number of particles at the vertex
   */
  int _npoint;

  /**
   *  Number of particle combinations at the vertex
   */
  unsigned int _nsize;

  /**
   * ParticleData pointers for the allowed incoming particles.
   */
  vector<PDPtr> _inpart;

  /**
   * PDG codes for the allowed incoming particles.
   */
  vector <int> _iinpart;

  /**
   * ParticleData pointers for the allowed outgoing particles.
   */
  vector<PDPtr> _outpart;

  /**
   * PDG codes for the allowed outgoing particles.
   */
  vector <int> _ioutpart;
  //@}

  /**
   * The overall coupling.
   */
  Complex _norm;

  /**
   * Whether or not to calculate the kinematic invariants for the vertex
   */
  bool _calckinematics;

  /**
   * Kinematica quantities needed for loop vertices
   */
  Energy2 _kine[5][5];
};
  
/**
 * Output the information on the vertex.
 */
ostream & operator<<(ostream &, const VertexBase &);

}
}


namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of VertexBase.
 */
template <>
struct BaseClassTrait<Herwig::Helicity::VertexBase,1> {
  /** Typedef of the base class of VertexBase. */
  typedef Interfaced NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::Helicity::VertexBase>
  : public ClassTraitsBase<Herwig::Helicity::VertexBase> {
  /**
   * Return the class name.
   */
  static string className() { return "Herwig++::Helicity::VertexBase"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwVertex.so"; }

};

}

#include "VertexBase.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "VertexBase.tcc"
#endif

#endif /* HERWIG_VertexBase_H */
