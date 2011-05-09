// -*- C++ -*-
#ifndef Herwig_MEPP2GauginoGauginoPowheg_H
#define Herwig_MEPP2GauginoGauginoPowheg_H
//
// This is the declaration of the MEPP2GauginoGauginoPowheg class.
//

#include "NLODrellYanBase.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the MEPP2GauginoGauginoPowheg class.
 *
 * @see \ref MEPP2GauginoGauginoPowhegInterfaces "The interfaces"
 * defined for MEPP2GauginoGauginoPowheg.
 */
class MEPP2GauginoGauginoPowheg: public NLODrellYanBase {

public:

  /**
   * The default constructor.
   */
  MEPP2GauginoGauginoPowheg();

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /**
   *  Finite part of the virtual correction
   */
  double finiteVirtual(Energy ms, Energy2 mb2,
		       vector<Complex> Cl, vector<Complex> Cr,
		       vector<double > Cs, vector<Complex> Ct,
		       vector<Complex> Cv) const;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEPP2GauginoGauginoPowheg & operator=(const MEPP2GauginoGauginoPowheg &);

};

}

#endif /* Herwig_MEPP2GauginoGauginoPowheg_H */
