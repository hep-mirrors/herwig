// -*- C++ -*-
//
// EWProcess.h is a part of Herwig - A multi-purpose Monte Carlo event generator
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
//
#ifndef HERWIG_EWProcess_H
#define HERWIG_EWProcess_H

namespace Herwig {

namespace EWProcess {

  /**
   * Enumerates the processes for which we have SCET Wilson Coefficients
   */
  enum Process { QQQQ,
		 QQQQiden,
		 QtQtQQ,
		 QQUU,
		 QtQtUU,
		 QQtRtR,
		 QQDD,
		 QtQtDD,
		 QQLL,
		 QQEE,
		 UUUU,
		 UUUUiden,
		 tRtRUU,
		 UUDD,
		 tRtRDD,
		 UULL,
		 UUEE,
		 DDDD,
		 DDDDiden,
		 DDLL,
		 DDEE,
		 LLLL,
		 LLLLiden,
		 LLEE,
		 EEEE,
		 EEEEiden,
		 QQWW,
		 QQPhiPhi,
		 QQWG,
		 QQBG,
		 QQGG,
		 QtQtGG,
		 LLWW,
		 LLPhiPhi,
		 UUBB,
		 UUPhiPhi,
		 UUBG,
		 UUGG,
		 tRtRGG,
		 DDBB,
		 DDPhiPhi,
		 DDBG,
		 DDGG,
		 EEBB,
		 EEPhiPhi,
		 LAST};
}
}

#endif // HERWIG_EWProcess_H 
