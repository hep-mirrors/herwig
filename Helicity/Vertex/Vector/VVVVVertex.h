// -*- C++ -*-
#ifndef HERWIG_VVVVVertex_H
#define HERWIG_VVVVVertex_H
//
// This is the declaration of the VVVVVertex class.

#include "Herwig++/Helicity/Vertex/VertexBase.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "VVVVVertex.fh"

namespace Herwig {
namespace Helicity{
using namespace ThePEG;

/** \ingroup Helicity
 *
 * This is the implementation of the four vector vertex. 
 * It is based on the VertexBase class for the storage of particles 
 * which are allowed to interact at the vertex.
 * Classes implementation a specific vertex should inherit from this 
 * one and implement the virtual setCoupling member.
 *
 * The form of the vertex is
 * \f[ic^2\left[
 *   2\epsilon_1\cdot\epsilon_2\epsilon_3\cdot\epsilon_4-
 *   \epsilon_1\cdot\epsilon_3\epsilon_2\cdot\epsilon_4-
 *   \epsilon_1\cdot\epsilon_4\epsilon_2\cdot\epsilon_3
 * \right]\f]
 *  optional the additional diagrams from the three point vertices can be included.
 *
 * @see VertexBase
 */
class VVVVVertex: public VertexBase {
  
public:

  /**
   * Default constructor.
   */
  inline VVVVVertex();

  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();
  
public:
  
  /**
   * Members to calculate the helicity amplitude expressions for vertices
   * and off-shell particles.
   */
  //@{
  /**
   * Evaluate the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Evaluation option, 0 just evaluate the four point vertex, 1
   * include all the three point diagrams as well.
   * @param vec1 The wavefunction for the first  vector.
   * @param vec2 The wavefunction for the second vector.
   * @param vec3 The wavefunction for the third  vector.
   * @param vec4 The wavefunction for the fourth vector.
   */
  Complex evaluate(Energy2 q2, int iopt,
		   const VectorWaveFunction & vec1, const VectorWaveFunction & vec2,
		   const VectorWaveFunction & vec3, const VectorWaveFunction & vec4);
  //@}

  /**
   * Calculate the couplings. This method is virtual and must be implemented in 
   * classes inheriting from this.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   * @param part4 The ParticleData pointer for the fourth particle.
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,tcPDPtr part3,
			   tcPDPtr part4)=0;
  
protected:
  
  /**
   * Set the order of the particles.
   * @param id1 The PDG code of the first  particle.
   * @param id2 The PDG code of the second particle.
   * @param id3 The PDG code of the third  particle.
   * @param id4 The PDG code of the fourth particle.
   */
  inline void setOrder(int id1,int id2,int id3,int id4);

  /**
   * Set the type of the vertex.
   * @param itype The type of vertex (QCD=1 or electroweak=2).
   */
  inline void setType(int itype);

  /**
   * Set the intermediate particles if including s/u/t channel terms.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param c1 The coupling for the first  particle.
   * @param c2 The coupling for the second particle.
   */
  inline void setIntermediate(tcPDPtr part1,tcPDPtr part2,Complex c1,Complex c2);
  
private:
  
  /**
   * Describe an abstract base class with persistent data.
   */
  static AbstractNoPIOClassDescription<VVVVVertex> initVVVVVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  VVVVVertex & operator=(const VVVVVertex &);
  
private:

  /**
   * Type of vertex 1=QCD 2=EW.
   */
  int _itype;

  /**  
   * Order of the particles.
   */
  int _iorder[4];

  /**
   *  Intermediate particles
   */
  tcPDPtr _inter[2];

  /**
   * Couplings of the intermediate particles.
   */
  Complex _coup[2];

};
}
}

#include "VVVVVertex.icc"

namespace ThePEG {
  
  /**
   * The following template specialization informs ThePEG about the
   * base class of VVVVVertex.
   */
  template <>
  struct BaseClassTrait<Herwig::Helicity::VVVVVertex,1> {
  /** Typedef of the base class of VVVVVertex. */
    typedef Herwig::Helicity::VertexBase NthBase;
  };
  
  /**
   * The following template specialization informs ThePEG about the
   * name of this class and the shared object where it is defined.
   */
  template <>
  struct ClassTraits<Herwig::Helicity::VVVVVertex>
    : public ClassTraitsBase<Herwig::Helicity::VVVVVertex> {

    /**
     * Return the class name.
     */
    static string className() { return "Herwig++::VVVVVertex"; }
  };
  
}

#endif /* HERWIG_VVVVVertex_H */
