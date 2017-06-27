/* KrkNLO Event Reweight (Dipole) */

#ifndef Herwig_KrknloEventReweight_H
#define Herwig_KrknloEventReweight_H

#include "Herwig/Shower/Dipole/Base/DipoleEventReweight.h"

namespace Herwig {

using namespace ThePEG;

class KrknloEventReweight : public DipoleEventReweight
{
public:
  double weight(const PPair& in, const PList& out, const PList& hard,  Ptr<AlphaSBase>::tptr as) const;
  double weightCascade(const PPair& in, const PList& out, const PList& hard,  Ptr<AlphaSBase>::tptr as) const;

  bool firstInteraction() const ;
  bool secondaryInteractions() const;
  
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
    void persistentInput(PersistentIStream & is, int);
    //@}

    /**
     * The standard Init function used to initialize the interfaces.
     * Called exactly once for each class by the class description system
     * before the main function starts or
     * when this class is dynamically loaded.
     */
    static void Init();

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


  // If needed, insert declarations of virtual function defined in the
  // InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

private:

  double scaleFactor_         = 1.0;

  unsigned int alphaSScale_R_ = 0;
  unsigned int alphaSScale_V_ = 0;
  
  unsigned int mode_          = 0;
  unsigned int mcPDF_         = 1;
  

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  KrknloEventReweight & operator=(const KrknloEventReweight &) = delete;

};

}

#endif /* Herwig_KrknloEventReweight_H */
