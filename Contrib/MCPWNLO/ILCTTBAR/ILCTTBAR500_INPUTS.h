/*  MC@NLO input file*/
#include<iostream>
#include <string>

class Input

{
 public:

  //  Input();
  // ~Input();

/* Center of mass energy/GeV */

  double cme() {return 500.;} /*This code is hard wired for 500 GeV centre of mass energy */

/* Polarization of e- (Pem) and e+ (Pep). +1 for RH and -1 for LH */

  int Pem() {return 1;}

  int Pep()  {return -1;}

/* Pole mass of Z/GeV */

 double Mz() {return 91.2;}

/* Number of events to generate */

 int nevgen() {return 100000;};

 
/* Random number seed */

 int rseed() {return 3;}

/* Implement POWHEG for ttbar production? i.e want NLO production?*/

bool POWHEGprod() {return true;}

/*Implement POWHEG for ttbar decays? i.e want NLO decays?*/

bool POWHEGdecay() {return true;}

 /* Implement truncated shower for production? */

bool truncpro() {return true;}

 /* Implement truncated shower for decays? */

bool truncdec() {return true;}

};



