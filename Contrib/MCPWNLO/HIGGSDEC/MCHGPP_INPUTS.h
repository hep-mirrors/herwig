/*  MC@NLO input file*/
#include<iostream>
#include <string>

class Input

{
 public:

  //  Input();
  // ~Input();

/* Center of mass energy/GeV */

  double cme() {return 1800.;}

  /* Vector Boson ID
   Options : 23  for Z boson
           : 24  for W boson
  */

  int VBID() {return 24;}

  /* Type of collider. Set to true if proton-proton or false if ppbar */

  bool pp() {return false;}

  /*  Mass of Higgs/GeV */
 
  double mH() {return 114.;}

 /* Pole mass of Z/GeV */

 double Mz() {return 91.2;}

/* Pole mass of W/GeV */

 double Mw() {return 80.3;}

 /* Mass of b quark/GeV */

 double mb() {return 5.0;}

 /* AlphaS (MH) */

 double alphasmH() {return 0.114138;}


 /* Number of parton flavours */

 int nf() {return 1;}

/* Maximum number of iterations for Newton-Raphson boost. About 4 is recommended  */

 int it() {return 4;}

/* Implement resolution cut for xb xbbar phase space? */

bool resolve() {return true;}

/* PDFset e.g cteq5m.LHgrid, MRST2004nlo.LHgrid etc Make sure the chosen set is compatible with rscheme defined above! */

 char* PDFset() { return "cteq5m.LHgrid"; }

/* Number of events for integration */

 int nevint() {return 10000;}

/* Number of events to generate */

 int nevgen() {return 10000;}

 
/* Random number seed */

 int rseed() {return 1;}

};



