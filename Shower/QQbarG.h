#ifndef _QQBARG_H_
#define _QQBARG_H_
#include "ShowerConfig.h"
#include "ShowerParticle.h"
//#include "ThePEG/Repository/UseRandom.h"

using namespace ThePEG;

namespace Herwig {

/** \ingroup Shower
 *  
 *  A comment should be added!
 */
class QQbarG {

public:

  /**
   * Constructors
   * given Q and m, symmetric kappa_ini by default: 
   */
  QQbarG(Energy Q, Energy m);

  /**
   * Give kta explicitly
   */
  QQbarG(Energy Q, Energy m, double k);

  ~QQbarG();

  /**
   * Set initial conditions in various ways.
   * kappa~_(antiquark) is always adjusted.
   */
  void setKtilde(double);
  void setKtildeSymm();
  void setKtildeLargest();
  void setKtildeSmallest();

  /**
   * The ratio of ME/PS for Vector exchange.
   */
  double getRatioV(double, double);
  double qWeight(double, double); 
  double qbarWeight(double, double);
  double qWeightX(Energy qtilde, double z);
  double qbarWeightX(Energy qtilde, double z);
  double getHard(double &, double &);
  //  vector<Lorentz5Momentum> applyHard(const ShowerParticleVector &p);
  vector<Lorentz5Momentum> applyHard(const PVector &p);
  // (x, xbar) -> (ktilde, z)
  double getZfromX(double, double);
  double getKfromX(double, double);
  inline Energy getQ() {return d_Q;};
  inline Energy getM() {return d_m;};
  
private:

  /**
   * CM energy and mass.
   */
  Energy d_Q, d_m;

  /**
   * The parameters rho and v, to be calculated in the constructor. 
   */
  double d_rho, d_v; 

  /**
   * The initial kappa-tilde values for radiation from quark and
   * antiquark resp. 
   */
  double d_kt1, d_kt2;

private:

  void setRho(double);

  /**
   * Set second value, given the first.
   */
  void setKtilde2();

  /**
   * Get x, xbar from kappa~, z.
   */
  void getXXbar(double, double, double &, double &);

  /**
   * CM energy and quark mass, to be given in constructor, only to
   * determine rho/v. kappa...
   * helper.
   */
  double u(double);

  /**
   * The ME, given x1, x2.
   */
  double PS(double, double);

public:

  double MEV(double, double);
  double MEA(double, double);

};

}
#endif
