// -*- C++ -*-
//
// Colorea.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_Colorea_H
#define HERWIG_Colorea_H
//
// This is the declaration of the Colorea class.
//
#include "Herwig/Shower/Dipole/Base/DipoleChain.fh"
#include "Herwig/Shower/Dipole/Base/Dipole.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup DipoleShower
 * \author Johannes Bellm
 *
 * This class implemets the method described in arXiv:1801.06113.
 * and it's extension described in arXiv:180X.XXXXX.
 * It allows to calculate the probability to rearrange
 * dipole chains with simple matrix elements.
 *
 * The implementation currently allows:
 * - rearranging pure FF dipole chains with
 *   processes: eeuugg eeuuggg and eeuugggg
 * - IF or FI dipoles can be rearranged with eeuugg. 
 * - II should not be rearranged.
 * - masses of the q qbar particles are set.
 *
 */
class Colorea {

public:

  /**
   * Default constructor
   */
  Colorea();
  
  /**
   * main function to rearrange the dipole chain.
   * - dipmax: chains with dipmax are rearranged with matix elements including
   *           dipmax dipoles:
   *                  dipmax=3 -> eeuugg
   *                  dipmax=4 -> eeuuggg
   *                  dipmax=5 -> eeuugggg
   *           if a chain contains more than dipmax dipoles, we treat the chain
   *           as long.
   * - diplong: this parameter allows to change the behaviour for Coloreas
   *            that are longer than dipmax (long chains).
   *            If diplong is 3 the chain is rearranged with eeuugg
   *            matrix elements.
   *            If diplong is 4 the chain is rearranged with eeuuggg
   *            matrix elements. diplong > 4 is currently not implemented.
   *            Note: We dont observe a difference between the diplong=3 or
   *                  diplong=4 treatment.
   */
  void rearrange(int dipmax,int diplong);
  
  /**
   * Set the chain to use in the rearrangement.
   */
  void setChain(Ptr<DipoleChain>::tptr ch){theChain=ch;};
  
   // Everything below is private. The rearranging is acting only on the current chain.
  private:
  /**
   * Main function for rearranging tripple dipoles.
   */
  void rearrange3(list<Dipole>::iterator dipi,
                  list<Dipole>::iterator dipj,
                  list<Dipole>::iterator dipk);
  
  /**
   * Function to rearrange triple dipoles in the all FF dipole case.
   */
  void rearrange3_FF_FF_FF(list<Dipole>::iterator dipi,
                           list<Dipole>::iterator dipj,
                           list<Dipole>::iterator dipk);
  
  /**
   * Function to rearrange triple dipoles if dipj is FI dipole.
   * This function is also used for IF.
   */
  void rearrange3_FF_FI_IF(list<Dipole>::iterator dipi,
                           list<Dipole>::iterator dipj,
                           list<Dipole>::iterator dipk);
  /**
   * Main function for rearranging 4 dipoles.
   * This funtion currently only implements the case of 4 FF dipoles.
   */
  void rearrange4(list<Dipole>::iterator dipi,
                  list<Dipole>::iterator dipj,
                  list<Dipole>::iterator dipk,
                  list<Dipole>::iterator dipl);
  /**
   * Main function for rearranging 5 dipoles.
   * This funtion currently only implements the case of 5 FF dipoles.
   */
  void rearrange5(list<Dipole>::iterator dipi,
                  list<Dipole>::iterator dipj,
                  list<Dipole>::iterator dipk,
                  list<Dipole>::iterator dipl,
                  list<Dipole>::iterator dipm);
  /**
   * Rearrange long chains.
   */
  void rearrangeLong(int diplong);
  
  /**
   * Produce the possible swapping or the dipj gluons.
   */
  bool produceSwapping(list<Dipole>::iterator dipi,
                       list<Dipole>::iterator dipj,
                       list<Dipole>::iterator dipk,bool FF);

    // helper function to remove particles from color line if color line connects them.
  void rmcol(tColinePtr A,tColinePtr B, list<Dipole>::iterator & dip);

    // remove connecting color lines from these dipoles
  void removeColors(list<Dipole>::iterator dipi,
                    list<Dipole>::iterator dipj,
                    list<Dipole>::iterator dipk);

    // Access the dipoles of the current chain.
  list<Dipole>& dipoles(); 

  private:
    // the current chain to rearrange.
  Ptr<DipoleChain>::tptr theChain;

};

inline ostream& operator << (ostream& os, const Colorea& ) {
  return os;
}

}

#endif /* HERWIG_Colorea_H */
