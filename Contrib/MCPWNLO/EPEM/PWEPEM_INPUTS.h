/*  POWHEG input file*/
#include<iostream>
#include <string>

class Input

{
 public:


/* Center of mass energy/GeV */

 double cme() {return 91.2;}

/* Pole mass of Z/GeV */

 double Mz() {return 91.2;}

 /* AlphaS (MZ) */

 double alphasmz() {return 0.118;}

 /* LambdaQCD/GeV */

 double Lambda() { return 0.2;} 

 /* Number of parton flavours */

 int nf() {return 5;}

 /* Use massless matrix element? Options: true or false */

bool massiveME() {return true;}

 /*Boost masses to true quark masses? (if using massless ME)Options: true or false */

bool bst() {return false;}

/* Maximum number of iterations for Newton-Raphson boost. 4 is recommended for sqrt(s)=91.2GeV */

 int it() {return 4;}


/* Number of events to generate */

 int nevgen() {return 100000;}

 
/* Random number seed */

 int rseed() {return 3;}

 /* Implement truncated shower? */

bool trunc() {return true;}

};



