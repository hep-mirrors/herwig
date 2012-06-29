/*  MC@NLO input file*/
#include<iostream>
#include <string>

class Input

{
 public:

  //  Input();
  // ~Input();

/* POWHEG or MCNLO ? True if MCNLO and false if POWHEG*/

bool MCNLO() {return true;}


 /* Implement truncated shower for POWHEG? */

bool trunc() {return true;}

/* Vector Boson ID 
   Options : 23  for Z boson
           : 24  for W boson  
*/

 int VBID() {return 24;}

/* Center of mass energy/GeV */

 double cme() {return 1800.;}

/*pp accelerator? True if pp, False if ppbar  */

bool acc()  {return false;}

/* Pole mass of Z/GeV */

 double Mz() {return 91.2;}

/* Pole mass of W/GeV */

 double Mw() {return 80.3;}

/* Zero width or Breit-Wigner resonance */

bool zerowidth() {return false;}

/* Width of W/GeV */

 double widthw() {return 2.032;}

/* Width of Z/GeV */

double widthz() {return 2.486;};

/* Number of half-widths from resonance Z/W peak/GeV */

double halfwidth() {return 20.;}

/* Renormalization scheme 
   Options : "DIS"
           : "MSbar"
*/

char* rscheme() {return "MSbar";}


/* PDFset e.g cteq5m.LHgrid, MRST2004nlo.LHgrid etc Make sure the chosen set is compatible with rscheme defined above! */

char* PDFset() { return "cteq5m.LHgrid"; }

/* Number of events for integration */

 int nevint() {return 10000;}

/* Number of events to generate */

 int nevgen() {return 1;}

/* Random number seed */

 int rseed() {return 5;}

};



