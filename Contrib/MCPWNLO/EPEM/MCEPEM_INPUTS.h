/*  MC@NLO input file*/
#include<iostream>
#include <string>

class Input

{
 public:

  //  Input();
  // ~Input();

/* Center of mass energy/GeV */

 double cme() {return 91.2;}

/* Pole mass of Z/GeV */

 double Mz() {return 91.2;}

 /* AlphaS (MZ) */

 double alphasmz() {return 0.118;}


 /* Number of parton flavours */

 int nf() {return 5;}

 /* Use massless matrix element? Options: true or false */

bool massiveME() {return true;}

 /*Boost masses to true quark masses? (if using massless ME)Options: true or false */

bool bst() {return true;}

/* Maximum number of iterations for Newton-Raphson boost. 4 is recommended for sqrt(s)=91.2GeV */

 int it() {return 4;}

/* Implement resolution cut? */

bool resolve() {return true;}

/* Number of events for integration */

 int nevint() {return 100000;}

/* Number of events to generate */

 int nevgen() {return 100000;}

 
/* Random number seed */

 int rseed() {return 1;}

};



