/*  MC@NLO input file*/
#include<iostream>
#include <string>

class Input

{
 public:

  //  Input();
  // ~Input();

  /* Do you want to interface this program to another which generates Higgs bosons independently? 
   If YES and you want to generate associated ZH events with this code, choose true. 
  If NO and want to use the internal associated Higgs production code, choose false. */

  bool iface() {return true;}


  /*********************************************************************************************************/
  /*           INFORMATION REQUIRED FOR INTERFACING ANOTHER PROCESS TO NLO HIGGS  DECAY                    */
  /*                                                                                                       */ 
  /*********************************************************************************************************/
 /*If interfacing, do you require the b, \bar{b} and g momenta in the Higgs rest frame or the lab frame?
    If rest frame, set to true. If lab frame, set to false. Please provide fileof lab frame Higgs' 5-momenta*/

  bool rest() {return false;}

  /* Have you got a properly formatted file containing Higgs 5-momenta (px, py, pz, E, mass) to interface? You
     must do if rest() == false! */

  bool ifile() {return true;}

  /* Address of the file */
  
  char* filename() { return "/usera/seyi/MCPWNLO/HIGGSDEC/Higgs.dat"; }

  /* Number of events in ifile */

  int nfile() {return 1; } 

  /*********************************************************************************************************/

/* Center of mass energy/GeV */

  double cme() {return 500.0;}

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

/* Maximum number of iterations for Newton-Raphson boost. About 4 is recommended  */

 int it() {return 4;}

/* Implement resolution cut for xb,xbbar phase space? */

bool resolve() {return true;}

/* Number of events for integration */

 int nevint() {return 10000;}

/* Number of events to generate */

 int nevgen() {return 100;}

 
/* Random number seed */

 int rseed() {return 1;}

};



