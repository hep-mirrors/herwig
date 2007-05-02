// file HwDebug.h  

#ifndef HERWIG_HwDebug_H
#define HERWIG_HwDebug_H

namespace Herwig {

/** \ingroup Utilities
 *
 *  Definition of the class HwDebug , for Herwig specific debugging.
 *
 *  The meaning of the four debugging levels are the following:
 *  1)  nodebug = 0    this should be used after the successfull
 *                     debugging phase.
 *  2)  minimal = 1    it just signal a wrong situations, but the
 *                      clutter of the log files is minimal: you can
 *                      use it for any number of events, and the only
 *                      price you have to pay is in terms of CPU performance,
 *                      because it goes through a huge number of checks
 *                      which slows down considerably the performance. 
 *   3)  full    = 10   it provides full, detailed information; although
 *                      this is very useful for debugging, it clutters
 *                      a lot the log file, increasing its size enormously.
 *                      Therefore you can use it only for a limited 
 *                      number of events (let's say <= 1000).
 *   4)  extreme = 100  it provides an extreme quantities of information,
 *                      but you should use it only for few events (let's
 *                      say <= 10).
 *  The typical use of these different levels is the following:
 *   --- when the code is finished, use  * extreme *  for about 10 events.
 *       After that the output makes sense, move to the next step;
 *   --- use  * full *  for about 1000 events.
 *       After that the output makes sense, move to the next step;
 *   --- use  * minimal *  for the largest number of events as possible
 *       (typically >= 10,000).
 *       After that the output does not show any error or serious
 *       warning message,  move to the next, last step;
 *   --- use  * nodebug *  for real use.
 * 
 *  During the development phase, it is useful to debug a single package
 *  (like: HardSubprocess, Shower, Hadronization, Decay) separately, one after
 *  the other. For this purpose, package-specific debugging levels have
 *  been introduced. An example of their use is the following. Suppose
 *  that we want to debug the package Shower, after having already
 *  debugged Hadronization and HardSubprocess. Then one should modify 
 *  HwDebug below as follows:
 *       minimal_HardSubprocess    = extreme+1,
 *       full_HardSubprocess       = extreme+1,
 *       extreme_HardSubprocess    = extreme+1,
 *       minimal_Shower            = minimal,
 *       full_Shower               = full,
 *       extreme_Shower            = extreme,
 *       minimal_Hadronization     = extreme+1,
 *       full_Hadronization        = extreme+1,
 *       extreme_Hadronization     = extreme+1,
 *       minimal_Decay             = minimal+1,
 *       full_Decay                = full+1,
 *       extreme_Decay             = extreme+1,
 */
  class HwDebug {

  public:

    /**
     * The different levels.
     */
    enum Levels {

      /**
       *  No debugging
       */
      noDebug = 0,

      /**
       * General.
       */
      minimal = 1,  
      full    = 10,
      extreme = 100,

      /**
       * Specific for the HardSubprocess part.
       */
      minimal_HardSubprocess = extreme + 1,
      full_HardSubprocess    = extreme + 1,
      extreme_HardSubprocess = extreme + 1,

      /**
       * Specific for the Shower part.
       */
      minimal_Shower = minimal,
      full_Shower    = full,
      extreme_Shower = extreme,

      /**
       * Specific for the Hadronization part
       *   minimal_Hadronization = extreme+1,
       *   full_Hadronization    = extreme+1,
       *   extreme_Hadronization = extreme+1,
       */
      minimal_Hadronization = minimal,
      full_Hadronization    = full,
      extreme_Hadronization = extreme,

      /**
       * Specific for the Decay part.
       */
      minimal_Decay = minimal+1,
      full_Decay    = full+1,
      extreme_Decay = extreme+1

    };

    /**
     * The current level.
     */
    static int level;
    
  };

}


#endif /* HERWIG_Debug_H */
