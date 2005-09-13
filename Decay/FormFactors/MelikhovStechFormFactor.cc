// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MelikhovStechFormFactor class.
//

#include "MelikhovStechFormFactor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "MelikhovStechFormFactor.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/ParticleData.h"

namespace Herwig {
using namespace ThePEG;

MelikhovStechFormFactor::MelikhovStechFormFactor() 
{
  // form factors for D to K
  addFormFactor(421,-321,0,-2,4,3);
  addFormFactor(411,-311,0,-1,4,3);
  _fplus0.push_back(0.78);_sigma1fp.push_back(0.24);_sigma2fp.push_back(0.00);
  _f00.push_back(0.78)   ;_sigma1f0.push_back(0.38);_sigma2f0.push_back(0.46);
  _fT0.push_back(0.75)   ;_sigma1fT.push_back(0.27);_sigma2fT.push_back(0.00);
  _V0.push_back(0.00)    ;_sigma1V0.push_back(0.00);_sigma2V0.push_back(0.00);
  _A00.push_back(0.00)   ;_sigma1A0.push_back(0.00);_sigma2A0.push_back(0.00);
  _A10.push_back(0.00)   ;_sigma1A1.push_back(0.00);_sigma2A1.push_back(0.00);
  _A20.push_back(0.00)   ;_sigma1A2.push_back(0.00);_sigma2A2.push_back(0.00);
  _T10.push_back(0.00)   ;_sigma1T1.push_back(0.00);_sigma2T1.push_back(0.00);
  _T20.push_back(0.00)   ;_sigma1T2.push_back(0.00);_sigma2T2.push_back(0.00);
  _T30.push_back(0.00)   ;_sigma1T3.push_back(0.00);_sigma2T3.push_back(0.00);
  _massP.push_back(1.97*GeV);_massV.push_back(2.11*GeV);
  _fplus0.push_back(0.78);_sigma1fp.push_back(0.24);_sigma2fp.push_back(0.00);
  _f00.push_back(0.78)   ;_sigma1f0.push_back(0.38);_sigma2f0.push_back(0.46);
  _fT0.push_back(0.75)   ;_sigma1fT.push_back(0.27);_sigma2fT.push_back(0.00);
  _V0.push_back(0.00)    ;_sigma1V0.push_back(0.00);_sigma2V0.push_back(0.00);
  _A00.push_back(0.00)   ;_sigma1A0.push_back(0.00);_sigma2A0.push_back(0.00);
  _A10.push_back(0.00)   ;_sigma1A1.push_back(0.00);_sigma2A1.push_back(0.00);
  _A20.push_back(0.00)   ;_sigma1A2.push_back(0.00);_sigma2A2.push_back(0.00);
  _T10.push_back(0.00)   ;_sigma1T1.push_back(0.00);_sigma2T1.push_back(0.00);
  _T20.push_back(0.00)   ;_sigma1T2.push_back(0.00);_sigma2T2.push_back(0.00);
  _T30.push_back(0.00)   ;_sigma1T3.push_back(0.00);_sigma2T3.push_back(0.00);
  _massP.push_back(1.97*GeV);_massV.push_back(2.11*GeV);
  // form factors for D to K*
  addFormFactor(421,-323,1,-2,4,3);
  addFormFactor(411,-313,1,-1,4,3);
  _fplus0.push_back(0.00);_sigma1fp.push_back(0.00);_sigma2fp.push_back(0.00);
  _f00.push_back(0.00)   ;_sigma1f0.push_back(0.00);_sigma2f0.push_back(0.00);
  _fT0.push_back(0.00)   ;_sigma1fT.push_back(0.00);_sigma2fT.push_back(0.00);
  _V0.push_back(1.03)    ;_sigma1V0.push_back(0.27);_sigma2V0.push_back(0.00);
  _A00.push_back(0.76)   ;_sigma1A0.push_back(0.17);_sigma2A0.push_back(0.00);
  _A10.push_back(0.66)   ;_sigma1A1.push_back(0.30);_sigma2A1.push_back(0.20);
  _A20.push_back(0.49)   ;_sigma1A2.push_back(0.67);_sigma2A2.push_back(0.16);
  _T10.push_back(0.78)   ;_sigma1T1.push_back(0.25);_sigma2T1.push_back(0.00);
  _T20.push_back(0.78)   ;_sigma1T2.push_back(0.02);_sigma2T2.push_back(1.80);
  _T30.push_back(0.45)   ;_sigma1T3.push_back(1.23);_sigma2T3.push_back(0.34);
  _massP.push_back(1.97*GeV);_massV.push_back(2.11*GeV);
  _fplus0.push_back(0.00);_sigma1fp.push_back(0.00);_sigma2fp.push_back(0.00);
  _f00.push_back(0.00)   ;_sigma1f0.push_back(0.00);_sigma2f0.push_back(0.00);
  _fT0.push_back(0.00)   ;_sigma1fT.push_back(0.00);_sigma2fT.push_back(0.00);
  _V0.push_back(1.03)    ;_sigma1V0.push_back(0.27);_sigma2V0.push_back(0.00);
  _A00.push_back(0.76)   ;_sigma1A0.push_back(0.17);_sigma2A0.push_back(0.00);
  _A10.push_back(0.66)   ;_sigma1A1.push_back(0.30);_sigma2A1.push_back(0.20);
  _A20.push_back(0.49)   ;_sigma1A2.push_back(0.67);_sigma2A2.push_back(0.16);
  _T10.push_back(0.78)   ;_sigma1T1.push_back(0.25);_sigma2T1.push_back(0.00);
  _T20.push_back(0.78)   ;_sigma1T2.push_back(0.02);_sigma2T2.push_back(1.80);
  _T30.push_back(0.45)   ;_sigma1T3.push_back(1.23);_sigma2T3.push_back(0.34);
  _massP.push_back(1.97*GeV);_massV.push_back(2.11*GeV);
  // form factors for D to pi
  addFormFactor(421,-211,0,-2,4,1);
  addFormFactor(421, 111,0,-2,4,2);
  addFormFactor(411, 111,0,-1,4,1);
  addFormFactor(411, 211,0,-1,4,2);
  _fplus0.push_back(0.69);_sigma1fp.push_back(0.50);_sigma2fp.push_back(0.00);
  _f00.push_back(0.69)   ;_sigma1f0.push_back(0.54);_sigma2f0.push_back(0.32);
  _fT0.push_back(0.60)   ;_sigma1fT.push_back(0.34);_sigma2fT.push_back(0.00);
  _V0.push_back(0.00)    ;_sigma1V0.push_back(0.00);_sigma2V0.push_back(0.00);
  _A00.push_back(0.00)   ;_sigma1A0.push_back(0.00);_sigma2A0.push_back(0.00);
  _A10.push_back(0.00)   ;_sigma1A1.push_back(0.00);_sigma2A1.push_back(0.00);
  _A20.push_back(0.00)   ;_sigma1A2.push_back(0.00);_sigma2A2.push_back(0.00);
  _T10.push_back(0.00)   ;_sigma1T1.push_back(0.00);_sigma2T1.push_back(0.00);
  _T20.push_back(0.00)   ;_sigma1T2.push_back(0.00);_sigma2T2.push_back(0.00);
  _T30.push_back(0.00)   ;_sigma1T3.push_back(0.00);_sigma2T3.push_back(0.00);
  _massP.push_back(1.87*GeV);_massV.push_back(2.01*GeV);
  _fplus0.push_back(0.69);_sigma1fp.push_back(0.50);_sigma2fp.push_back(0.00);
  _f00.push_back(0.69)   ;_sigma1f0.push_back(0.54);_sigma2f0.push_back(0.32);
  _fT0.push_back(0.60)   ;_sigma1fT.push_back(0.34);_sigma2fT.push_back(0.00);
  _V0.push_back(0.00)    ;_sigma1V0.push_back(0.00);_sigma2V0.push_back(0.00);
  _A00.push_back(0.00)   ;_sigma1A0.push_back(0.00);_sigma2A0.push_back(0.00);
  _A10.push_back(0.00)   ;_sigma1A1.push_back(0.00);_sigma2A1.push_back(0.00);
  _A20.push_back(0.00)   ;_sigma1A2.push_back(0.00);_sigma2A2.push_back(0.00);
  _T10.push_back(0.00)   ;_sigma1T1.push_back(0.00);_sigma2T1.push_back(0.00);
  _T20.push_back(0.00)   ;_sigma1T2.push_back(0.00);_sigma2T2.push_back(0.00);
  _T30.push_back(0.00)   ;_sigma1T3.push_back(0.00);_sigma2T3.push_back(0.00);
  _massP.push_back(1.87*GeV);_massV.push_back(2.01*GeV);
  _fplus0.push_back(0.69);_sigma1fp.push_back(0.50);_sigma2fp.push_back(0.00);
  _f00.push_back(0.69)   ;_sigma1f0.push_back(0.54);_sigma2f0.push_back(0.32);
  _fT0.push_back(0.60)   ;_sigma1fT.push_back(0.34);_sigma2fT.push_back(0.00);
  _V0.push_back(0.00)    ;_sigma1V0.push_back(0.00);_sigma2V0.push_back(0.00);
  _A00.push_back(0.00)   ;_sigma1A0.push_back(0.00);_sigma2A0.push_back(0.00);
  _A10.push_back(0.00)   ;_sigma1A1.push_back(0.00);_sigma2A1.push_back(0.00);
  _A20.push_back(0.00)   ;_sigma1A2.push_back(0.00);_sigma2A2.push_back(0.00);
  _T10.push_back(0.00)   ;_sigma1T1.push_back(0.00);_sigma2T1.push_back(0.00);
  _T20.push_back(0.00)   ;_sigma1T2.push_back(0.00);_sigma2T2.push_back(0.00);
  _T30.push_back(0.00)   ;_sigma1T3.push_back(0.00);_sigma2T3.push_back(0.00);
  _massP.push_back(1.87*GeV);_massV.push_back(2.01*GeV);
  _fplus0.push_back(0.69);_sigma1fp.push_back(0.50);_sigma2fp.push_back(0.00);
  _f00.push_back(0.69)   ;_sigma1f0.push_back(0.54);_sigma2f0.push_back(0.32);
  _fT0.push_back(0.60)   ;_sigma1fT.push_back(0.34);_sigma2fT.push_back(0.00);
  _V0.push_back(0.00)    ;_sigma1V0.push_back(0.00);_sigma2V0.push_back(0.00);
  _A00.push_back(0.00)   ;_sigma1A0.push_back(0.00);_sigma2A0.push_back(0.00);
  _A10.push_back(0.00)   ;_sigma1A1.push_back(0.00);_sigma2A1.push_back(0.00);
  _A20.push_back(0.00)   ;_sigma1A2.push_back(0.00);_sigma2A2.push_back(0.00);
  _T10.push_back(0.00)   ;_sigma1T1.push_back(0.00);_sigma2T1.push_back(0.00);
  _T20.push_back(0.00)   ;_sigma1T2.push_back(0.00);_sigma2T2.push_back(0.00);
  _T30.push_back(0.00)   ;_sigma1T3.push_back(0.00);_sigma2T3.push_back(0.00);
  _massP.push_back(1.87*GeV);_massV.push_back(2.01*GeV);
  // form factors for D to rho
  addFormFactor(421,-213,1,-2,4,1);
  addFormFactor(421, 113,1,-2,4,2);
  addFormFactor(411, 113,1,-1,4,1);
  addFormFactor(411, 213,1,-1,4,2);
  _fplus0.push_back(0.00);_sigma1fp.push_back(0.00);_sigma2fp.push_back(0.00);
  _f00.push_back(0.00)   ;_sigma1f0.push_back(0.00);_sigma2f0.push_back(0.00);
  _fT0.push_back(0.00)   ;_sigma1fT.push_back(0.00);_sigma2fT.push_back(0.00);
  _V0.push_back(0.90)    ;_sigma1V0.push_back(0.46);_sigma2V0.push_back(0.00);
  _A00.push_back(0.66)   ;_sigma1A0.push_back(0.36);_sigma2A0.push_back(0.00);
  _A10.push_back(0.59)   ;_sigma1A1.push_back(0.50);_sigma2A1.push_back(0.00);
  _A20.push_back(0.49)   ;_sigma1A2.push_back(0.89);_sigma2A2.push_back(0.00);
  _T10.push_back(0.66)   ;_sigma1T1.push_back(0.44);_sigma2T1.push_back(0.00);
  _T20.push_back(0.66)   ;_sigma1T2.push_back(0.38);_sigma2T2.push_back(0.50);
  _T30.push_back(0.31)   ;_sigma1T3.push_back(1.10);_sigma2T3.push_back(0.17);
  _massP.push_back(1.87*GeV);_massV.push_back(2.01*GeV);
  _fplus0.push_back(0.00);_sigma1fp.push_back(0.00);_sigma2fp.push_back(0.00);
  _f00.push_back(0.00)   ;_sigma1f0.push_back(0.00);_sigma2f0.push_back(0.00);
  _fT0.push_back(0.00)   ;_sigma1fT.push_back(0.00);_sigma2fT.push_back(0.00);
  _V0.push_back(0.90)    ;_sigma1V0.push_back(0.46);_sigma2V0.push_back(0.00);
  _A00.push_back(0.66)   ;_sigma1A0.push_back(0.36);_sigma2A0.push_back(0.00);
  _A10.push_back(0.59)   ;_sigma1A1.push_back(0.50);_sigma2A1.push_back(0.00);
  _A20.push_back(0.49)   ;_sigma1A2.push_back(0.89);_sigma2A2.push_back(0.00);
  _T10.push_back(0.66)   ;_sigma1T1.push_back(0.44);_sigma2T1.push_back(0.00);
  _T20.push_back(0.66)   ;_sigma1T2.push_back(0.38);_sigma2T2.push_back(0.50);
  _T30.push_back(0.31)   ;_sigma1T3.push_back(1.10);_sigma2T3.push_back(0.17);
  _massP.push_back(1.87*GeV);_massV.push_back(2.01*GeV);
  _fplus0.push_back(0.00);_sigma1fp.push_back(0.00);_sigma2fp.push_back(0.00);
  _f00.push_back(0.00)   ;_sigma1f0.push_back(0.00);_sigma2f0.push_back(0.00);
  _fT0.push_back(0.00)   ;_sigma1fT.push_back(0.00);_sigma2fT.push_back(0.00);
  _V0.push_back(0.90)    ;_sigma1V0.push_back(0.46);_sigma2V0.push_back(0.00);
  _A00.push_back(0.66)   ;_sigma1A0.push_back(0.36);_sigma2A0.push_back(0.00);
  _A10.push_back(0.59)   ;_sigma1A1.push_back(0.50);_sigma2A1.push_back(0.00);
  _A20.push_back(0.49)   ;_sigma1A2.push_back(0.89);_sigma2A2.push_back(0.00);
  _T10.push_back(0.66)   ;_sigma1T1.push_back(0.44);_sigma2T1.push_back(0.00);
  _T20.push_back(0.66)   ;_sigma1T2.push_back(0.38);_sigma2T2.push_back(0.50);
  _T30.push_back(0.31)   ;_sigma1T3.push_back(1.10);_sigma2T3.push_back(0.17);
  _massP.push_back(1.87*GeV);_massV.push_back(2.01*GeV);
  _fplus0.push_back(0.00);_sigma1fp.push_back(0.00);_sigma2fp.push_back(0.00);
  _f00.push_back(0.00)   ;_sigma1f0.push_back(0.00);_sigma2f0.push_back(0.00);
  _fT0.push_back(0.00)   ;_sigma1fT.push_back(0.00);_sigma2fT.push_back(0.00);
  _V0.push_back(0.90)    ;_sigma1V0.push_back(0.46);_sigma2V0.push_back(0.00);
  _A00.push_back(0.66)   ;_sigma1A0.push_back(0.36);_sigma2A0.push_back(0.00);
  _A10.push_back(0.59)   ;_sigma1A1.push_back(0.50);_sigma2A1.push_back(0.00);
  _A20.push_back(0.49)   ;_sigma1A2.push_back(0.89);_sigma2A2.push_back(0.00);
  _T10.push_back(0.66)   ;_sigma1T1.push_back(0.44);_sigma2T1.push_back(0.00);
  _T20.push_back(0.66)   ;_sigma1T2.push_back(0.38);_sigma2T2.push_back(0.50);
  _T30.push_back(0.31)   ;_sigma1T3.push_back(1.10);_sigma2T3.push_back(0.17);
  _massP.push_back(1.87*GeV);_massV.push_back(2.01*GeV);
  // B to D
  addFormFactor(521,-421,0,2,-5,-4);
  addFormFactor(511,-411,0,2,-5,-4);
  _fplus0.push_back(0.67);_sigma1fp.push_back(0.57);_sigma2fp.push_back(0.00);
  _f00.push_back(0.67)   ;_sigma1f0.push_back(0.78);_sigma2f0.push_back(0.00);
  _fT0.push_back(0.69)   ;_sigma1fT.push_back(0.56);_sigma2fT.push_back(0.00);
  _V0.push_back(0.00)    ;_sigma1V0.push_back(0.00);_sigma2V0.push_back(0.00);
  _A00.push_back(0.00)   ;_sigma1A0.push_back(0.00);_sigma2A0.push_back(0.00);
  _A10.push_back(0.00)   ;_sigma1A1.push_back(0.00);_sigma2A1.push_back(0.00);
  _A20.push_back(0.00)   ;_sigma1A2.push_back(0.00);_sigma2A2.push_back(0.00);
  _T10.push_back(0.00)   ;_sigma1T1.push_back(0.00);_sigma2T1.push_back(0.00);
  _T20.push_back(0.00)   ;_sigma1T2.push_back(0.00);_sigma2T2.push_back(0.00);
  _T30.push_back(0.00)   ;_sigma1T3.push_back(0.00);_sigma2T3.push_back(0.00);
  _massP.push_back(6.4*GeV);_massV.push_back(6.4*GeV);
  _fplus0.push_back(0.67);_sigma1fp.push_back(0.57);_sigma2fp.push_back(0.00);
  _f00.push_back(0.67)   ;_sigma1f0.push_back(0.78);_sigma2f0.push_back(0.00);
  _fT0.push_back(0.69)   ;_sigma1fT.push_back(0.56);_sigma2fT.push_back(0.00);
  _V0.push_back(0.00)    ;_sigma1V0.push_back(0.00);_sigma2V0.push_back(0.00);
  _A00.push_back(0.00)   ;_sigma1A0.push_back(0.00);_sigma2A0.push_back(0.00);
  _A10.push_back(0.00)   ;_sigma1A1.push_back(0.00);_sigma2A1.push_back(0.00);
  _A20.push_back(0.00)   ;_sigma1A2.push_back(0.00);_sigma2A2.push_back(0.00);
  _T10.push_back(0.00)   ;_sigma1T1.push_back(0.00);_sigma2T1.push_back(0.00);
  _T20.push_back(0.00)   ;_sigma1T2.push_back(0.00);_sigma2T2.push_back(0.00);
  _T30.push_back(0.00)   ;_sigma1T3.push_back(0.00);_sigma2T3.push_back(0.00);
  _massP.push_back(6.4*GeV);_massV.push_back(6.4*GeV);
  // B to D*
  addFormFactor(521,-423,1,2,-5,-4);
  addFormFactor(511,-413,1,1,-5,-4);
  _fplus0.push_back(0.00);_sigma1fp.push_back(0.00);_sigma2fp.push_back(0.00);
  _f00.push_back(0.00)   ;_sigma1f0.push_back(0.00);_sigma2f0.push_back(0.00);
  _fT0.push_back(0.00)   ;_sigma1fT.push_back(0.00);_sigma2fT.push_back(0.00);
  _V0.push_back(0.76)    ;_sigma1V0.push_back(0.57);_sigma2V0.push_back(0.00);
  _A00.push_back(0.69)   ;_sigma1A0.push_back(0.58);_sigma2A0.push_back(0.00);
  _A10.push_back(0.66)   ;_sigma1A1.push_back(0.78);_sigma2A1.push_back(0.00);
  _A20.push_back(0.62)   ;_sigma1A2.push_back(1.40);_sigma2A2.push_back(0.41);
  _T10.push_back(0.68)   ;_sigma1T1.push_back(0.57);_sigma2T1.push_back(0.00);
  _T20.push_back(0.68)   ;_sigma1T2.push_back(0.64);_sigma2T2.push_back(0.00);
  _T30.push_back(0.33)   ;_sigma1T3.push_back(1.46);_sigma2T3.push_back(0.50);
  _massP.push_back(6.4*GeV);_massV.push_back(6.4*GeV);
  _fplus0.push_back(0.00);_sigma1fp.push_back(0.00);_sigma2fp.push_back(0.00);
  _f00.push_back(0.00)   ;_sigma1f0.push_back(0.00);_sigma2f0.push_back(0.00);
  _fT0.push_back(0.00)   ;_sigma1fT.push_back(0.00);_sigma2fT.push_back(0.00);
  _V0.push_back(0.76)    ;_sigma1V0.push_back(0.57);_sigma2V0.push_back(0.00);
  _A00.push_back(0.69)   ;_sigma1A0.push_back(0.58);_sigma2A0.push_back(0.00);
  _A10.push_back(0.66)   ;_sigma1A1.push_back(0.78);_sigma2A1.push_back(0.00);
  _A20.push_back(0.62)   ;_sigma1A2.push_back(1.40);_sigma2A2.push_back(0.41);
  _T10.push_back(0.68)   ;_sigma1T1.push_back(0.57);_sigma2T1.push_back(0.00);
  _T20.push_back(0.68)   ;_sigma1T2.push_back(0.64);_sigma2T2.push_back(0.00);
  _T30.push_back(0.33)   ;_sigma1T3.push_back(1.46);_sigma2T3.push_back(0.50);
  _massP.push_back(6.4*GeV);_massV.push_back(6.4*GeV);
  // B to K
  addFormFactor(521,321,0,2,-5,-3);
  addFormFactor(511,311,0,1,-5,-3);
  _fplus0.push_back(0.36);_sigma1fp.push_back(0.43);_sigma2fp.push_back(0.00);
  _f00.push_back(0.36)   ;_sigma1f0.push_back(0.70);_sigma2f0.push_back(0.27);
  _fT0.push_back(0.35)   ;_sigma1fT.push_back(0.43);_sigma2fT.push_back(0.00);
  _V0.push_back(0.00)    ;_sigma1V0.push_back(0.00);_sigma2V0.push_back(0.00);
  _A00.push_back(0.00)   ;_sigma1A0.push_back(0.00);_sigma2A0.push_back(0.00);
  _A10.push_back(0.00)   ;_sigma1A1.push_back(0.00);_sigma2A1.push_back(0.00);
  _A20.push_back(0.00)   ;_sigma1A2.push_back(0.00);_sigma2A2.push_back(0.00);
  _T10.push_back(0.00)   ;_sigma1T1.push_back(0.00);_sigma2T1.push_back(0.00);
  _T20.push_back(0.00)   ;_sigma1T2.push_back(0.00);_sigma2T2.push_back(0.00);
  _T30.push_back(0.00)   ;_sigma1T3.push_back(0.00);_sigma2T3.push_back(0.00);
  _massP.push_back(5.37*GeV);_massV.push_back(5.42*GeV);
  _fplus0.push_back(0.36);_sigma1fp.push_back(0.43);_sigma2fp.push_back(0.00);
  _f00.push_back(0.36)   ;_sigma1f0.push_back(0.70);_sigma2f0.push_back(0.27);
  _fT0.push_back(0.35)   ;_sigma1fT.push_back(0.43);_sigma2fT.push_back(0.00);
  _V0.push_back(0.00)    ;_sigma1V0.push_back(0.00);_sigma2V0.push_back(0.00);
  _A00.push_back(0.00)   ;_sigma1A0.push_back(0.00);_sigma2A0.push_back(0.00);
  _A10.push_back(0.00)   ;_sigma1A1.push_back(0.00);_sigma2A1.push_back(0.00);
  _A20.push_back(0.00)   ;_sigma1A2.push_back(0.00);_sigma2A2.push_back(0.00);
  _T10.push_back(0.00)   ;_sigma1T1.push_back(0.00);_sigma2T1.push_back(0.00);
  _T20.push_back(0.00)   ;_sigma1T2.push_back(0.00);_sigma2T2.push_back(0.00);
  _T30.push_back(0.00)   ;_sigma1T3.push_back(0.00);_sigma2T3.push_back(0.00);
  _massP.push_back(5.37*GeV);_massV.push_back(5.42*GeV);
  // B to K*
  addFormFactor(521,323,1,2,-5,-3);
  addFormFactor(511,313,1,1,-5,-3);
  _fplus0.push_back(0.00);_sigma1fp.push_back(0.00);_sigma2fp.push_back(0.00);
  _f00.push_back(0.00)   ;_sigma1f0.push_back(0.00);_sigma2f0.push_back(0.00);
  _fT0.push_back(0.00)   ;_sigma1fT.push_back(0.00);_sigma2fT.push_back(0.00);
  _V0.push_back(0.44)    ;_sigma1V0.push_back(0.45);_sigma2V0.push_back(0.00);
  _A00.push_back(0.45)   ;_sigma1A0.push_back(0.46);_sigma2A0.push_back(0.00);
  _A10.push_back(0.36)   ;_sigma1A1.push_back(0.64);_sigma2A1.push_back(0.36);
  _A20.push_back(0.32)   ;_sigma1A2.push_back(1.23);_sigma2A2.push_back(0.38);
  _T10.push_back(0.39)   ;_sigma1T1.push_back(0.45);_sigma2T1.push_back(0.00);
  _T20.push_back(0.39)   ;_sigma1T2.push_back(0.72);_sigma2T2.push_back(0.62);
  _T30.push_back(0.27)   ;_sigma1T3.push_back(1.31);_sigma2T3.push_back(0.41);
  _massP.push_back(5.37*GeV);_massV.push_back(5.42*GeV);
  _fplus0.push_back(0.00);_sigma1fp.push_back(0.00);_sigma2fp.push_back(0.00);
  _f00.push_back(0.00)   ;_sigma1f0.push_back(0.00);_sigma2f0.push_back(0.00);
  _fT0.push_back(0.00)   ;_sigma1fT.push_back(0.00);_sigma2fT.push_back(0.00);
  _V0.push_back(0.44)    ;_sigma1V0.push_back(0.45);_sigma2V0.push_back(0.00);
  _A00.push_back(0.45)   ;_sigma1A0.push_back(0.46);_sigma2A0.push_back(0.00);
  _A10.push_back(0.36)   ;_sigma1A1.push_back(0.64);_sigma2A1.push_back(0.36);
  _A20.push_back(0.32)   ;_sigma1A2.push_back(1.23);_sigma2A2.push_back(0.38);
  _T10.push_back(0.39)   ;_sigma1T1.push_back(0.45);_sigma2T1.push_back(0.00);
  _T20.push_back(0.39)   ;_sigma1T2.push_back(0.72);_sigma2T2.push_back(0.62);
  _T30.push_back(0.27)   ;_sigma1T3.push_back(1.31);_sigma2T3.push_back(0.41);
  _massP.push_back(5.37*GeV);_massV.push_back(5.42*GeV);
  // B to pi
  addFormFactor(521, 111,0,2,-5,-2);
  addFormFactor(511,-211,0,1,-5,-2);
  addFormFactor(521, 211,0,2,-5,-1);
  addFormFactor(511, 111,0,1,-5,-1);
  _fplus0.push_back(0.29);_sigma1fp.push_back(0.48);_sigma2fp.push_back(0.00);
  _f00.push_back(0.29)   ;_sigma1f0.push_back(0.76);_sigma2f0.push_back(0.28);
  _fT0.push_back(0.28)   ;_sigma1fT.push_back(0.48);_sigma2fT.push_back(0.00);
  _V0.push_back(0.00)    ;_sigma1V0.push_back(0.00);_sigma2V0.push_back(0.00);
  _A00.push_back(0.00)   ;_sigma1A0.push_back(0.00);_sigma2A0.push_back(0.00);
  _A10.push_back(0.00)   ;_sigma1A1.push_back(0.00);_sigma2A1.push_back(0.00);
  _A20.push_back(0.00)   ;_sigma1A2.push_back(0.00);_sigma2A2.push_back(0.00);
  _T10.push_back(0.00)   ;_sigma1T1.push_back(0.00);_sigma2T1.push_back(0.00);
  _T20.push_back(0.00)   ;_sigma1T2.push_back(0.00);_sigma2T2.push_back(0.00);
  _T30.push_back(0.00)   ;_sigma1T3.push_back(0.00);_sigma2T3.push_back(0.00);
  _massP.push_back(5.27*GeV);_massV.push_back(5.32*GeV);
  _fplus0.push_back(0.29);_sigma1fp.push_back(0.48);_sigma2fp.push_back(0.00);
  _f00.push_back(0.29)   ;_sigma1f0.push_back(0.76);_sigma2f0.push_back(0.28);
  _fT0.push_back(0.28)   ;_sigma1fT.push_back(0.48);_sigma2fT.push_back(0.00);
  _V0.push_back(0.00)    ;_sigma1V0.push_back(0.00);_sigma2V0.push_back(0.00);
  _A00.push_back(0.00)   ;_sigma1A0.push_back(0.00);_sigma2A0.push_back(0.00);
  _A10.push_back(0.00)   ;_sigma1A1.push_back(0.00);_sigma2A1.push_back(0.00);
  _A20.push_back(0.00)   ;_sigma1A2.push_back(0.00);_sigma2A2.push_back(0.00);
  _T10.push_back(0.00)   ;_sigma1T1.push_back(0.00);_sigma2T1.push_back(0.00);
  _T20.push_back(0.00)   ;_sigma1T2.push_back(0.00);_sigma2T2.push_back(0.00);
  _T30.push_back(0.00)   ;_sigma1T3.push_back(0.00);_sigma2T3.push_back(0.00);
  _massP.push_back(5.27*GeV);_massV.push_back(5.32*GeV);
  _fplus0.push_back(0.29);_sigma1fp.push_back(0.48);_sigma2fp.push_back(0.00);
  _f00.push_back(0.29)   ;_sigma1f0.push_back(0.76);_sigma2f0.push_back(0.28);
  _fT0.push_back(0.28)   ;_sigma1fT.push_back(0.48);_sigma2fT.push_back(0.00);
  _V0.push_back(0.00)    ;_sigma1V0.push_back(0.00);_sigma2V0.push_back(0.00);
  _A00.push_back(0.00)   ;_sigma1A0.push_back(0.00);_sigma2A0.push_back(0.00);
  _A10.push_back(0.00)   ;_sigma1A1.push_back(0.00);_sigma2A1.push_back(0.00);
  _A20.push_back(0.00)   ;_sigma1A2.push_back(0.00);_sigma2A2.push_back(0.00);
  _T10.push_back(0.00)   ;_sigma1T1.push_back(0.00);_sigma2T1.push_back(0.00);
  _T20.push_back(0.00)   ;_sigma1T2.push_back(0.00);_sigma2T2.push_back(0.00);
  _T30.push_back(0.00)   ;_sigma1T3.push_back(0.00);_sigma2T3.push_back(0.00);
  _massP.push_back(5.27*GeV);_massV.push_back(5.32*GeV);
  _fplus0.push_back(0.29);_sigma1fp.push_back(0.48);_sigma2fp.push_back(0.00);
  _f00.push_back(0.29)   ;_sigma1f0.push_back(0.76);_sigma2f0.push_back(0.28);
  _fT0.push_back(0.28)   ;_sigma1fT.push_back(0.48);_sigma2fT.push_back(0.00);
  _V0.push_back(0.00)    ;_sigma1V0.push_back(0.00);_sigma2V0.push_back(0.00);
  _A00.push_back(0.00)   ;_sigma1A0.push_back(0.00);_sigma2A0.push_back(0.00);
  _A10.push_back(0.00)   ;_sigma1A1.push_back(0.00);_sigma2A1.push_back(0.00);
  _A20.push_back(0.00)   ;_sigma1A2.push_back(0.00);_sigma2A2.push_back(0.00);
  _T10.push_back(0.00)   ;_sigma1T1.push_back(0.00);_sigma2T1.push_back(0.00);
  _T20.push_back(0.00)   ;_sigma1T2.push_back(0.00);_sigma2T2.push_back(0.00);
  _T30.push_back(0.00)   ;_sigma1T3.push_back(0.00);_sigma2T3.push_back(0.00);
  _massP.push_back(5.27*GeV);_massV.push_back(5.32*GeV);
  // B to rho
  addFormFactor(521, 113,1,2,-5,-2);
  addFormFactor(511,-213,1,1,-5,-2);
  addFormFactor(521, 213,1,2,-5,-1);
  addFormFactor(511, 113,1,1,-5,-1);
  _fplus0.push_back(0.00);_sigma1fp.push_back(0.00);_sigma2fp.push_back(0.00);
  _f00.push_back(0.00)   ;_sigma1f0.push_back(0.00);_sigma2f0.push_back(0.00);
  _fT0.push_back(0.00)   ;_sigma1fT.push_back(0.00);_sigma2fT.push_back(0.00);
  _V0.push_back(0.31)    ;_sigma1V0.push_back(0.59);_sigma2V0.push_back(0.00);
  _A00.push_back(0.30)   ;_sigma1A0.push_back(0.54);_sigma2A0.push_back(0.00);
  _A10.push_back(0.26)   ;_sigma1A1.push_back(0.73);_sigma2A1.push_back(0.10);
  _A20.push_back(0.24)   ;_sigma1A2.push_back(1.40);_sigma2A2.push_back(0.50);
  _T10.push_back(0.27)   ;_sigma1T1.push_back(0.60);_sigma2T1.push_back(0.00);
  _T20.push_back(0.27)   ;_sigma1T2.push_back(0.74);_sigma2T2.push_back(0.19);
  _T30.push_back(0.19)   ;_sigma1T3.push_back(1.42);_sigma2T3.push_back(0.51);
  _massP.push_back(5.27*GeV);_massV.push_back(5.32*GeV);
  _fplus0.push_back(0.00);_sigma1fp.push_back(0.00);_sigma2fp.push_back(0.00);
  _f00.push_back(0.00)   ;_sigma1f0.push_back(0.00);_sigma2f0.push_back(0.00);
  _fT0.push_back(0.00)   ;_sigma1fT.push_back(0.00);_sigma2fT.push_back(0.00);
  _V0.push_back(0.31)    ;_sigma1V0.push_back(0.59);_sigma2V0.push_back(0.00);
  _A00.push_back(0.30)   ;_sigma1A0.push_back(0.54);_sigma2A0.push_back(0.00);
  _A10.push_back(0.26)   ;_sigma1A1.push_back(0.73);_sigma2A1.push_back(0.10);
  _A20.push_back(0.24)   ;_sigma1A2.push_back(1.40);_sigma2A2.push_back(0.50);
  _T10.push_back(0.27)   ;_sigma1T1.push_back(0.60);_sigma2T1.push_back(0.00);
  _T20.push_back(0.27)   ;_sigma1T2.push_back(0.74);_sigma2T2.push_back(0.19);
  _T30.push_back(0.19)   ;_sigma1T3.push_back(1.42);_sigma2T3.push_back(0.51);
  _massP.push_back(5.27*GeV);_massV.push_back(5.32*GeV);
  _fplus0.push_back(0.00);_sigma1fp.push_back(0.00);_sigma2fp.push_back(0.00);
  _f00.push_back(0.00)   ;_sigma1f0.push_back(0.00);_sigma2f0.push_back(0.00);
  _fT0.push_back(0.00)   ;_sigma1fT.push_back(0.00);_sigma2fT.push_back(0.00);
  _V0.push_back(0.31)    ;_sigma1V0.push_back(0.59);_sigma2V0.push_back(0.00);
  _A00.push_back(0.30)   ;_sigma1A0.push_back(0.54);_sigma2A0.push_back(0.00);
  _A10.push_back(0.26)   ;_sigma1A1.push_back(0.73);_sigma2A1.push_back(0.10);
  _A20.push_back(0.24)   ;_sigma1A2.push_back(1.40);_sigma2A2.push_back(0.50);
  _T10.push_back(0.27)   ;_sigma1T1.push_back(0.60);_sigma2T1.push_back(0.00);
  _T20.push_back(0.27)   ;_sigma1T2.push_back(0.74);_sigma2T2.push_back(0.19);
  _T30.push_back(0.19)   ;_sigma1T3.push_back(1.42);_sigma2T3.push_back(0.51);
  _massP.push_back(5.27*GeV);_massV.push_back(5.32*GeV);
  _fplus0.push_back(0.00);_sigma1fp.push_back(0.00);_sigma2fp.push_back(0.00);
  _f00.push_back(0.00)   ;_sigma1f0.push_back(0.00);_sigma2f0.push_back(0.00);
  _fT0.push_back(0.00)   ;_sigma1fT.push_back(0.00);_sigma2fT.push_back(0.00);
  _V0.push_back(0.31)    ;_sigma1V0.push_back(0.59);_sigma2V0.push_back(0.00);
  _A00.push_back(0.30)   ;_sigma1A0.push_back(0.54);_sigma2A0.push_back(0.00);
  _A10.push_back(0.26)   ;_sigma1A1.push_back(0.73);_sigma2A1.push_back(0.10);
  _A20.push_back(0.24)   ;_sigma1A2.push_back(1.40);_sigma2A2.push_back(0.50);
  _T10.push_back(0.27)   ;_sigma1T1.push_back(0.60);_sigma2T1.push_back(0.00);
  _T20.push_back(0.27)   ;_sigma1T2.push_back(0.74);_sigma2T2.push_back(0.19);
  _T30.push_back(0.19)   ;_sigma1T3.push_back(1.42);_sigma2T3.push_back(0.51);
  _massP.push_back(5.27*GeV);_massV.push_back(5.32*GeV);
  // D_s to K
  addFormFactor(431,311,0,-3,4,1);
  addFormFactor(431,321,0,-3,4,2);
  _fplus0.push_back(0.72);_sigma1fp.push_back(0.20);_sigma2fp.push_back(0.00);
  _f00.push_back(0.72)   ;_sigma1f0.push_back(0.41);_sigma2f0.push_back(0.70);
  _fT0.push_back(0.77)   ;_sigma1fT.push_back(0.24);_sigma2fT.push_back(0.00);
  _V0.push_back(0.00)    ;_sigma1V0.push_back(0.00);_sigma2V0.push_back(0.00);
  _A00.push_back(0.00)   ;_sigma1A0.push_back(0.00);_sigma2A0.push_back(0.00);
  _A10.push_back(0.00)   ;_sigma1A1.push_back(0.00);_sigma2A1.push_back(0.00);
  _A20.push_back(0.00)   ;_sigma1A2.push_back(0.00);_sigma2A2.push_back(0.00);
  _T10.push_back(0.00)   ;_sigma1T1.push_back(0.00);_sigma2T1.push_back(0.00);
  _T20.push_back(0.00)   ;_sigma1T2.push_back(0.00);_sigma2T2.push_back(0.00);
  _T30.push_back(0.00)   ;_sigma1T3.push_back(0.00);_sigma2T3.push_back(0.00);
  _massP.push_back(1.87*GeV);_massV.push_back(2.01*GeV);
  _fplus0.push_back(0.72);_sigma1fp.push_back(0.20);_sigma2fp.push_back(0.00);
  _f00.push_back(0.72)   ;_sigma1f0.push_back(0.41);_sigma2f0.push_back(0.70);
  _fT0.push_back(0.77)   ;_sigma1fT.push_back(0.24);_sigma2fT.push_back(0.00);
  _V0.push_back(0.00)    ;_sigma1V0.push_back(0.00);_sigma2V0.push_back(0.00);
  _A00.push_back(0.00)   ;_sigma1A0.push_back(0.00);_sigma2A0.push_back(0.00);
  _A10.push_back(0.00)   ;_sigma1A1.push_back(0.00);_sigma2A1.push_back(0.00);
  _A20.push_back(0.00)   ;_sigma1A2.push_back(0.00);_sigma2A2.push_back(0.00);
  _T10.push_back(0.00)   ;_sigma1T1.push_back(0.00);_sigma2T1.push_back(0.00);
  _T20.push_back(0.00)   ;_sigma1T2.push_back(0.00);_sigma2T2.push_back(0.00);
  _T30.push_back(0.00)   ;_sigma1T3.push_back(0.00);_sigma2T3.push_back(0.00);
  _massP.push_back(1.87*GeV);_massV.push_back(2.01*GeV);
  // D_s to K*
  addFormFactor(431,313,1,-3,4,1);
  addFormFactor(431,323,1,-3,4,2);
  _fplus0.push_back(0.00);_sigma1fp.push_back( 0.00);_sigma2fp.push_back(0.00);
  _f00.push_back(0.00)   ;_sigma1f0.push_back( 0.00);_sigma2f0.push_back(0.00);
  _fT0.push_back(0.00)   ;_sigma1fT.push_back( 0.00);_sigma2fT.push_back(0.00);
  _V0.push_back(1.04)    ;_sigma1V0.push_back( 0.24);_sigma2V0.push_back(0.00);
  _A00.push_back(0.67)   ;_sigma1A0.push_back( 0.20);_sigma2A0.push_back(0.00);
  _A10.push_back(0.57)   ;_sigma1A1.push_back( 0.29);_sigma2A1.push_back(0.42);
  _A20.push_back(0.42)   ;_sigma1A2.push_back( 0.58);_sigma2A2.push_back(0.00);
  _T10.push_back(0.71)   ;_sigma1T1.push_back( 0.22);_sigma2T1.push_back(0.00);
  _T20.push_back(0.71)   ;_sigma1T2.push_back(-0.06);_sigma2T2.push_back(0.44);
  _T30.push_back(0.45)   ;_sigma1T3.push_back( 1.08);_sigma2T3.push_back(0.68);
  _massP.push_back(1.87*GeV);_massV.push_back(2.01*GeV);
  _fplus0.push_back(0.00);_sigma1fp.push_back( 0.00);_sigma2fp.push_back(0.00);
  _f00.push_back(0.00)   ;_sigma1f0.push_back( 0.00);_sigma2f0.push_back(0.00);
  _fT0.push_back(0.00)   ;_sigma1fT.push_back( 0.00);_sigma2fT.push_back(0.00);
  _V0.push_back(1.04)    ;_sigma1V0.push_back( 0.24);_sigma2V0.push_back(0.00);
  _A00.push_back(0.67)   ;_sigma1A0.push_back( 0.20);_sigma2A0.push_back(0.00);
  _A10.push_back(0.57)   ;_sigma1A1.push_back( 0.29);_sigma2A1.push_back(0.42);
  _A20.push_back(0.42)   ;_sigma1A2.push_back( 0.58);_sigma2A2.push_back(0.00);
  _T10.push_back(0.71)   ;_sigma1T1.push_back( 0.22);_sigma2T1.push_back(0.00);
  _T20.push_back(0.71)   ;_sigma1T2.push_back(-0.06);_sigma2T2.push_back(0.44);
  _T30.push_back(0.45)   ;_sigma1T3.push_back( 1.08);_sigma2T3.push_back(0.68);
  _massP.push_back(1.87*GeV);_massV.push_back(2.01*GeV);
  // D_s to eta
  addFormFactor(431,221,0,-3,4,3);
  _fplus0.push_back(0.78);_sigma1fp.push_back(0.23);_sigma2fp.push_back(0.00);
  _f00.push_back(0.78)   ;_sigma1f0.push_back(0.33);_sigma2f0.push_back(0.38);
  _fT0.push_back(0.80)   ;_sigma1fT.push_back(0.24);_sigma2fT.push_back(0.00);
  _V0.push_back(0.00)    ;_sigma1V0.push_back(0.00);_sigma2V0.push_back(0.00);
  _A00.push_back(0.00)   ;_sigma1A0.push_back(0.00);_sigma2A0.push_back(0.00);
  _A10.push_back(0.00)   ;_sigma1A1.push_back(0.00);_sigma2A1.push_back(0.00);
  _A20.push_back(0.00)   ;_sigma1A2.push_back(0.00);_sigma2A2.push_back(0.00);
  _T10.push_back(0.00)   ;_sigma1T1.push_back(0.00);_sigma2T1.push_back(0.00);
  _T20.push_back(0.00)   ;_sigma1T2.push_back(0.00);_sigma2T2.push_back(0.00);
  _T30.push_back(0.00)   ;_sigma1T3.push_back(0.00);_sigma2T3.push_back(0.00);
  _massP.push_back(1.97*GeV);_massV.push_back(2.11*GeV);
  // D_s to eta'
  addFormFactor(431,331,0,-3,4,3);
  _fplus0.push_back(0.78);_sigma1fp.push_back(0.23);_sigma2fp.push_back(0.00);
  _f00.push_back(0.78)   ;_sigma1f0.push_back(0.21);_sigma2f0.push_back(0.76);
  _fT0.push_back(0.94)   ;_sigma1fT.push_back(0.24);_sigma2fT.push_back(0.00);
  _V0.push_back(0.00)    ;_sigma1V0.push_back(0.00);_sigma2V0.push_back(0.00);
  _A00.push_back(0.00)   ;_sigma1A0.push_back(0.00);_sigma2A0.push_back(0.00);
  _A10.push_back(0.00)   ;_sigma1A1.push_back(0.00);_sigma2A1.push_back(0.00);
  _A20.push_back(0.00)   ;_sigma1A2.push_back(0.00);_sigma2A2.push_back(0.00);
  _T10.push_back(0.00)   ;_sigma1T1.push_back(0.00);_sigma2T1.push_back(0.00);
  _T20.push_back(0.00)   ;_sigma1T2.push_back(0.00);_sigma2T2.push_back(0.00);
  _T30.push_back(0.00)   ;_sigma1T3.push_back(0.00);_sigma2T3.push_back(0.00);
  _massP.push_back(1.97*GeV);_massV.push_back(2.11*GeV);
  // D_s to phi
  addFormFactor(431,333,1,-3,4,3);
  _fplus0.push_back(0.00);_sigma1fp.push_back(0.00);_sigma2fp.push_back(0.00);
  _f00.push_back(0.00)   ;_sigma1f0.push_back(0.00);_sigma2f0.push_back(0.00);
  _fT0.push_back(0.00)   ;_sigma1fT.push_back(0.00);_sigma2fT.push_back(0.00);
  _V0.push_back(1.10)    ;_sigma1V0.push_back(0.26);_sigma2V0.push_back(0.00);
  _A00.push_back(0.73)   ;_sigma1A0.push_back(0.10);_sigma2A0.push_back(0.00);
  _A10.push_back(0.64)   ;_sigma1A1.push_back(0.29);_sigma2A1.push_back(0.00);
  _A20.push_back(0.47)   ;_sigma1A2.push_back(0.63);_sigma2A2.push_back(0.00);
  _T10.push_back(0.77)   ;_sigma1T1.push_back(0.25);_sigma2T1.push_back(0.00);
  _T20.push_back(0.77)   ;_sigma1T2.push_back(0.02);_sigma2T2.push_back(2.01);
  _T30.push_back(0.46)   ;_sigma1T3.push_back(1.34);_sigma2T3.push_back(0.45);
  _massP.push_back(1.97*GeV);_massV.push_back(2.11*GeV);
  // B_s to K
  addFormFactor( 531,-311  ,0, 3,-5,-1);
  addFormFactor( 531,-321  ,0, 3,-5,-2);
  _fplus0.push_back(0.31);_sigma1fp.push_back(0.63);_sigma2fp.push_back(0.33);
  _f00.push_back(0.31)   ;_sigma1f0.push_back(0.93);_sigma2f0.push_back(0.70);
  _fT0.push_back(0.31)   ;_sigma1fT.push_back(0.61);_sigma2fT.push_back(0.30);
  _V0.push_back(0.00)    ;_sigma1V0.push_back(0.00);_sigma2V0.push_back(0.00);
  _A00.push_back(0.00)   ;_sigma1A0.push_back(0.00);_sigma2A0.push_back(0.00);
  _A10.push_back(0.00)   ;_sigma1A1.push_back(0.00);_sigma2A1.push_back(0.00);
  _A20.push_back(0.00)   ;_sigma1A2.push_back(0.00);_sigma2A2.push_back(0.00);
  _T10.push_back(0.00)   ;_sigma1T1.push_back(0.00);_sigma2T1.push_back(0.00);
  _T20.push_back(0.00)   ;_sigma1T2.push_back(0.00);_sigma2T2.push_back(0.00);
  _T30.push_back(0.00)   ;_sigma1T3.push_back(0.00);_sigma2T3.push_back(0.00);
  _massP.push_back(5.27*GeV);_massV.push_back(5.32*GeV);
  _fplus0.push_back(0.31);_sigma1fp.push_back(0.63);_sigma2fp.push_back(0.33);
  _f00.push_back(0.31)   ;_sigma1f0.push_back(0.93);_sigma2f0.push_back(0.70);
  _fT0.push_back(0.31)   ;_sigma1fT.push_back(0.61);_sigma2fT.push_back(0.30);
  _V0.push_back(0.00)    ;_sigma1V0.push_back(0.00);_sigma2V0.push_back(0.00);
  _A00.push_back(0.00)   ;_sigma1A0.push_back(0.00);_sigma2A0.push_back(0.00);
  _A10.push_back(0.00)   ;_sigma1A1.push_back(0.00);_sigma2A1.push_back(0.00);
  _A20.push_back(0.00)   ;_sigma1A2.push_back(0.00);_sigma2A2.push_back(0.00);
  _T10.push_back(0.00)   ;_sigma1T1.push_back(0.00);_sigma2T1.push_back(0.00);
  _T20.push_back(0.00)   ;_sigma1T2.push_back(0.00);_sigma2T2.push_back(0.00);
  _T30.push_back(0.00)   ;_sigma1T3.push_back(0.00);_sigma2T3.push_back(0.00);
  _massP.push_back(5.27*GeV);_massV.push_back(5.32*GeV);
  // B_s to K*
  addFormFactor( 531,-313  ,1, 3,-5,-1);
  addFormFactor( 531,-323  ,1, 3,-5,-2);
  _fplus0.push_back(0.00);_sigma1fp.push_back(0.00);_sigma2fp.push_back(0.00);
  _f00.push_back(0.00)   ;_sigma1f0.push_back(0.00);_sigma2f0.push_back(0.00);
  _fT0.push_back(0.00)   ;_sigma1fT.push_back(0.00);_sigma2fT.push_back(0.00);
  _V0.push_back(0.38)    ;_sigma1V0.push_back(0.66);_sigma2V0.push_back(0.30);
  _A00.push_back(0.37)   ;_sigma1A0.push_back(0.60);_sigma2A0.push_back(0.16);
  _A10.push_back(0.29)   ;_sigma1A1.push_back(0.86);_sigma2A1.push_back(0.60);
  _A20.push_back(0.26)   ;_sigma1A2.push_back(1.32);_sigma2A2.push_back(0.54);
  _T10.push_back(0.32)   ;_sigma1T1.push_back(0.66);_sigma2T1.push_back(0.31);
  _T20.push_back(0.32)   ;_sigma1T2.push_back(0.98);_sigma2T2.push_back(0.90);
  _T30.push_back(0.23)   ;_sigma1T3.push_back(1.42);_sigma2T3.push_back(0.62);
  _massP.push_back(5.27*GeV);_massV.push_back(5.32*GeV);
  _fplus0.push_back(0.00);_sigma1fp.push_back(0.00);_sigma2fp.push_back(0.00);
  _f00.push_back(0.00)   ;_sigma1f0.push_back(0.00);_sigma2f0.push_back(0.00);
  _fT0.push_back(0.00)   ;_sigma1fT.push_back(0.00);_sigma2fT.push_back(0.00);
  _V0.push_back(0.38)    ;_sigma1V0.push_back(0.66);_sigma2V0.push_back(0.30);
  _A00.push_back(0.37)   ;_sigma1A0.push_back(0.60);_sigma2A0.push_back(0.16);
  _A10.push_back(0.29)   ;_sigma1A1.push_back(0.86);_sigma2A1.push_back(0.60);
  _A20.push_back(0.26)   ;_sigma1A2.push_back(1.32);_sigma2A2.push_back(0.54);
  _T10.push_back(0.32)   ;_sigma1T1.push_back(0.66);_sigma2T1.push_back(0.31);
  _T20.push_back(0.32)   ;_sigma1T2.push_back(0.98);_sigma2T2.push_back(0.90);
  _T30.push_back(0.23)   ;_sigma1T3.push_back(1.42);_sigma2T3.push_back(0.62);
  _massP.push_back(5.27*GeV);_massV.push_back(5.32*GeV);
  // B_s to eta
  addFormFactor( 531, 221  ,0, 3,-5,-3);
  _fplus0.push_back(0.36);_sigma1fp.push_back(0.60);_sigma2fp.push_back(0.20);
  _f00.push_back(0.36)   ;_sigma1f0.push_back(0.80);_sigma2f0.push_back(0.40);
  _fT0.push_back(0.36)   ;_sigma1fT.push_back(0.58);_sigma2fT.push_back(0.18);
  _V0.push_back(0.00)    ;_sigma1V0.push_back(0.00);_sigma2V0.push_back(0.00);
  _A00.push_back(0.00)   ;_sigma1A0.push_back(0.00);_sigma2A0.push_back(0.00);
  _A10.push_back(0.00)   ;_sigma1A1.push_back(0.00);_sigma2A1.push_back(0.00);
  _A20.push_back(0.00)   ;_sigma1A2.push_back(0.00);_sigma2A2.push_back(0.00);
  _T10.push_back(0.00)   ;_sigma1T1.push_back(0.00);_sigma2T1.push_back(0.00);
  _T20.push_back(0.00)   ;_sigma1T2.push_back(0.00);_sigma2T2.push_back(0.00);
  _T30.push_back(0.00)   ;_sigma1T3.push_back(0.00);_sigma2T3.push_back(0.00);
  _massP.push_back(5.37*GeV);_massV.push_back(5.42*GeV);
  // B_s to eta'
  addFormFactor( 531, 331  ,0, 3,-5,-3);
  _fplus0.push_back(0.36);_sigma1fp.push_back(0.60);_sigma2fp.push_back(0.20);
  _f00.push_back(0.36)   ;_sigma1f0.push_back(0.80);_sigma2f0.push_back(0.45);
  _fT0.push_back(0.39)   ;_sigma1fT.push_back(0.58);_sigma2fT.push_back(0.18);
  _V0.push_back(0.00)    ;_sigma1V0.push_back(0.00);_sigma2V0.push_back(0.00);
  _A00.push_back(0.00)   ;_sigma1A0.push_back(0.00);_sigma2A0.push_back(0.00);
  _A10.push_back(0.00)   ;_sigma1A1.push_back(0.00);_sigma2A1.push_back(0.00);
  _A20.push_back(0.00)   ;_sigma1A2.push_back(0.00);_sigma2A2.push_back(0.00);
  _T10.push_back(0.00)   ;_sigma1T1.push_back(0.00);_sigma2T1.push_back(0.00);
  _T20.push_back(0.00)   ;_sigma1T2.push_back(0.00);_sigma2T2.push_back(0.00);
  _T30.push_back(0.00)   ;_sigma1T3.push_back(0.00);_sigma2T3.push_back(0.00);
  _massP.push_back(5.37*GeV);_massV.push_back(5.42*GeV);
  // B_s to phi
  addFormFactor( 531, 333  ,1, 3,-5,-3);
  _fplus0.push_back(0.00);_sigma1fp.push_back(0.00);_sigma2fp.push_back(0.00);
  _f00.push_back(0.00)   ;_sigma1f0.push_back(0.00);_sigma2f0.push_back(0.00);
  _fT0.push_back(0.00)   ;_sigma1fT.push_back(0.00);_sigma2fT.push_back(0.00);
  _V0.push_back(0.44)    ;_sigma1V0.push_back(0.62);_sigma2V0.push_back(0.20);
  _A00.push_back(0.42)   ;_sigma1A0.push_back(0.55);_sigma2A0.push_back(0.12);
  _A10.push_back(0.34)   ;_sigma1A1.push_back(0.73);_sigma2A1.push_back(0.42);
  _A20.push_back(0.31)   ;_sigma1A2.push_back(1.30);_sigma2A2.push_back(0.52);
  _T10.push_back(0.38)   ;_sigma1T1.push_back(0.62);_sigma2T1.push_back(0.20);
  _T20.push_back(0.38)   ;_sigma1T2.push_back(0.83);_sigma2T2.push_back(0.71);
  _T30.push_back(0.26)   ;_sigma1T3.push_back(1.41);_sigma2T3.push_back(0.57);
  _massP.push_back(5.37*GeV);_massV.push_back(5.42*GeV);
  // set the initial number of modes
  initialModes(numberOfFactors());
  // eta-eta' mixing angle
  _thetaeta = 2./9.*pi;
}

MelikhovStechFormFactor::~MelikhovStechFormFactor() {}

void MelikhovStechFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _fplus0 << _sigma1fp << _sigma2fp << _f00 << _sigma1f0 << _sigma2f0 << _fT0 
     << _sigma1fT << _sigma2fT << _V0 << _sigma1V0 << _sigma2V0 << _A00 << _sigma1A0 
     << _sigma2A0 << _A10 << _sigma1A1 << _sigma2A1 << _A20 << _sigma1A2 << _sigma2A2 
     << _T10 << _sigma1T1 << _sigma2T1 << _T20 << _sigma1T2 << _sigma2T2 << _T30 
     << _sigma1T3 << _sigma2T3 << _massP << _massV << _thetaeta;
}

void MelikhovStechFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _fplus0 >> _sigma1fp >> _sigma2fp >> _f00 >> _sigma1f0 >> _sigma2f0 >> _fT0 
     >> _sigma1fT >> _sigma2fT >> _V0 >> _sigma1V0 >> _sigma2V0 >> _A00 >> _sigma1A0 
     >> _sigma2A0 >> _A10 >> _sigma1A1 >> _sigma2A1 >> _A20 >> _sigma1A2 >> _sigma2A2 
     >> _T10 >> _sigma1T1 >> _sigma2T1 >> _T20 >> _sigma1T2 >> _sigma2T2 >> _T30 
     >> _sigma1T3 >> _sigma2T3 >> _massP >> _massV >> _thetaeta;
}

ClassDescription<MelikhovStechFormFactor> MelikhovStechFormFactor::initMelikhovStechFormFactor;
// Definition of the static class description member.

void MelikhovStechFormFactor::Init() {

  static ClassDocumentation<MelikhovStechFormFactor> documentation
    ("The MelikhovStechFormFactor class is the implementation of the"
     " form factors from hep-ph/0001113");

  static ParVector<MelikhovStechFormFactor,double> interfaceFPlus0
    ("FPlus0",
     "The value of the F_+ form factor at q^2=0",
     &MelikhovStechFormFactor::_fplus0, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceFpsigma_1
    ("F+sigma_1",
     "The value of sigma_1 for the F_+ form factor",
     &MelikhovStechFormFactor::_sigma1fp, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceFpsigma_2
    ("F+sigma_2",
     "The value of sigma_2 for the F_+ form factor",
     &MelikhovStechFormFactor::_sigma2fp, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceF00
    ("F00",
     "The value of the F_0 form factor at q^2=0",
     &MelikhovStechFormFactor::_f00, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceF0sigma_1
    ("F0sigma_1",
     "The value of sigma_1 for the F_0 form factor",
     &MelikhovStechFormFactor::_sigma1f0, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceF0sigma_2
    ("F0sigma_2",
     "The value of sigma_2 for the F_0 form factor",
     &MelikhovStechFormFactor::_sigma2f0, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceFT0
    ("FT0",
     "The value of the F_T form factor at q^2=0",
     &MelikhovStechFormFactor::_fT0, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceFTsigma_1
    ("FTsigma_1",
     "The value of sigma_1 for the F_T form factor",
     &MelikhovStechFormFactor::_sigma1fT, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceFTsigma_2
    ("FTsigma_2",
     "The value of sigma_2 for the F_T form factor",
     &MelikhovStechFormFactor::_sigma2fT, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceV00
    ("V00",
     "The value of the V_0 form factor at q^2=0",
     &MelikhovStechFormFactor::_V0, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceV0sigma_1
    ("V0sigma_1",
     "The value of sigma_1 for the V_0 form factor",
     &MelikhovStechFormFactor::_sigma1V0, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceV0sigma_2
    ("V0sigma_2",
     "The value of sigma_2 for the V_0 form factor",
     &MelikhovStechFormFactor::_sigma2V0, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceA00
    ("A00",
     "The value of the A_0 form factor at q^2=0",
     &MelikhovStechFormFactor::_A00, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceA0sigma_1
    ("A0sigma_1",
     "The value of sigma_1 for the A_0 form factor",
     &MelikhovStechFormFactor::_sigma1A0, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceA0sigma_2
    ("A0sigma_2",
     "The value of sigma_2 for the A_0 form factor",
     &MelikhovStechFormFactor::_sigma2A0, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceA10
    ("A10",
     "The value of the A_1 form factor at q^2=0",
     &MelikhovStechFormFactor::_A10, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceA1sigma_1
    ("A1sigma_1",
     "The value of sigma_1 for the A_1 form factor",
     &MelikhovStechFormFactor::_sigma1A1, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceA1sigma_2
    ("A1sigma_2",
     "The value of sigma_2 for the A_1 form factor",
     &MelikhovStechFormFactor::_sigma2A1, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceA20
    ("A20",
     "The value of the A_2 form factor at q^2=0",
     &MelikhovStechFormFactor::_A20, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceA2sigma_1
    ("A2sigma_1",
     "The value of sigma_1 for the A_2 form factor",
     &MelikhovStechFormFactor::_sigma1A2, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceA2sigma_2
    ("A2sigma_2",
     "The value of sigma_2 for the A_2 form factor",
     &MelikhovStechFormFactor::_sigma2A2, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceT10
    ("T10",
     "The value of the T_1 form factor at q^2=0",
     &MelikhovStechFormFactor::_T10, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceT1sigma_1
    ("T1sigma_1",
     "The value of sigma_1 for the T_1 form factor",
     &MelikhovStechFormFactor::_sigma1T1, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceT1sigma_2
    ("T1sigma_2",
     "The value of sigma_2 for the T_1 form factor",
     &MelikhovStechFormFactor::_sigma2T1, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceT20
    ("T20",
     "The value of the T_2 form factor at q^2=0",
     &MelikhovStechFormFactor::_T20, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceT2sigma_1
    ("T2sigma_1",
     "The value of sigma_1 for the T_2 form factor",
     &MelikhovStechFormFactor::_sigma1T2, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceT2sigma_2
    ("T2sigma_2",
     "The value of sigma_2 for the T_2 form factor",
     &MelikhovStechFormFactor::_sigma2T2, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceT30
    ("T30",
     "The value of the T_3 form factor at q^2=0",
     &MelikhovStechFormFactor::_T30, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceT3sigma_1
    ("T3sigma_1",
     "The value of sigma_1 for the T_3 form factor",
     &MelikhovStechFormFactor::_sigma1T3, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceT3sigma_2
    ("T3sigma_2",
     "The value of sigma_2 for the T_3 form factor",
     &MelikhovStechFormFactor::_sigma2T3, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,Energy> interfaceMassP
    ("MassP",
     "The mass of the pseudoscalar for the q^2 dependence of the form factors.",
     &MelikhovStechFormFactor::_massP, GeV, -1, 0.0*GeV, 0*GeV, 10.*GeV,
     false, false, false);

  static ParVector<MelikhovStechFormFactor,Energy> interfaceMassV
    ("MassV",
     "The mass of the vector for the q^2 dependence of the form factors.",
     &MelikhovStechFormFactor::_massV, GeV, -1, 0.0*GeV, 0*GeV, 10.*GeV,
     false, false, false);

  static Parameter< MelikhovStechFormFactor,double> interfaceThetaEtaEtaPrime
    ("ThetaEtaEtaPrime",
     "The eta-eta' mixing angle",
     & MelikhovStechFormFactor::_thetaeta, 2.*pi/9., -pi, pi,
     false, false, true);
}

// form-factor for scalar to scalar
void MelikhovStechFormFactor::ScalarScalarFormFactor(Energy2 q2,unsigned int mode,
						     int id0,int id1,
						     Energy m0, Energy m1,Complex & f0,
						     Complex & fp) const
{
  double ratio(q2/_massV[mode]/_massV[mode]);
  fp = _fplus0[mode]/(1.-ratio)/(1.-ratio*(_sigma1fp[mode]-_sigma2fp[mode]*ratio));
  f0 = _f00[mode]              /(1.-ratio*(_sigma1f0[mode]-_sigma2f0[mode]*ratio));
  int jspin,spect,inquark,outquark;
  formFactorInfo(mode,jspin,spect,inquark,outquark);
  if(abs(outquark)==abs(spect))
    {
      double fact;
      if(id1==ParticleID::eta)
	{
	  if(abs(outquark)==3){fact=-sin(_thetaeta);}
	  else{fact=cos(_thetaeta)*sqrt(0.5);}
	}
      else if(id1==ParticleID::etaprime)
	{
	  if(abs(outquark)==3){fact=cos(_thetaeta);}
	  else{fact=sin(_thetaeta);}
	}
      else if(id1==ParticleID::pi0&&abs(outquark)==1){fact=-sqrt(0.5);}
      else{fact= sqrt(0.5);}
      f0*=fact;fp*=fact;
    }
}
  
void MelikhovStechFormFactor::ScalarScalarSigmaFormFactor(Energy2 q2,unsigned int mode,
							  int id0,int id1,
							  Energy m0, Energy m1,
							  Complex & fT) const
{
  double ratio(q2/_massV[mode]/_massV[mode]);
  fT = _fT0[mode]              /(1.-ratio*(_sigma1fT[mode]-_sigma2fT[mode]*ratio));
  int jspin,spect,inquark,outquark;
  formFactorInfo(mode,jspin,spect,inquark,outquark);
  if(abs(outquark)==abs(spect))
    {
      double fact;
      if(id1==ParticleID::eta)
	{
	  if(abs(outquark)==3){fact=-sin(_thetaeta);}
	  else{fact=cos(_thetaeta)*sqrt(0.5);}
	}
      else if(id1==ParticleID::etaprime)
	{
	  if(abs(outquark)==3){fact=cos(_thetaeta);}
	  else{fact=sin(_thetaeta);}
	}
      else if(id1==ParticleID::pi0&&abs(outquark)==1){fact=-sqrt(0.5);}
      else{fact= sqrt(0.5);}
      fT*=fact;
    }
}

void MelikhovStechFormFactor::ScalarVectorFormFactor(Energy2 q2, unsigned int mode,
						     int id0, int id1,
						     Energy m0, Energy m1,Complex & A0,
						     Complex & A1,Complex & A2,
						     Complex & V) const
{
  double ratio(q2/_massV[mode]/_massV[mode]),ratioP(q2/_massP[mode]/_massP[mode]);
  A0= _A00[mode]/(1.-ratioP)/(1.-ratioP*(_sigma1A0[mode]-_sigma2A0[mode]*ratioP));
  A1= _A10[mode]            /(1.-ratio *(_sigma1A1[mode]-_sigma2A1[mode]*ratio));
  A2= _A20[mode]            /(1.-ratio *(_sigma1A2[mode]-_sigma2A2[mode]*ratio));
  V = _V0[mode] /(1.-ratio )/(1.-ratio *(_sigma1V0[mode]-_sigma2V0[mode]*ratio ));
  int jspin,spect,inquark,outquark;
  formFactorInfo(mode,jspin,spect,inquark,outquark);
  if(abs(outquark)==abs(spect)&&abs(spect)<3)
    {
      double fact(sqrt(0.5));
      if(id1==ParticleID::rho0&&abs(outquark)==1){fact=-fact;}
      A0*=fact;A1*=fact;A2*=fact;V*=fact;
    }
}

void MelikhovStechFormFactor::ScalarVectorSigmaFormFactor(Energy2 q2,unsigned int mode,
							  int id0,int id1,
							  Energy m0, Energy m1,
							  Complex & T1,Complex & T2,
							  Complex & T3) const
{
  double ratio(q2/_massV[mode]/_massV[mode]);
  T1= _T10[mode]/(1.-ratio)/(1.-ratio*(_sigma1T1[mode]-_sigma2T1[mode]*ratio));
  T2= _T20[mode]/(1.-ratio)/(1.-ratio*(_sigma1T2[mode]-_sigma2T2[mode]*ratio));
  T3= _T30[mode]           /(1.-ratio*(_sigma1T3[mode]-_sigma2T3[mode]*ratio));
  int jspin,spect,inquark,outquark;
  formFactorInfo(mode,jspin,spect,inquark,outquark);
  if(abs(outquark)==abs(spect)&&abs(spect)<3)
    {
      double fact(sqrt(0.5));
      if(id1==ParticleID::rho0&&abs(outquark)==1){fact=-fact;}
      T1*=fact;T2*=fact;T3*=fact;
    }
}

void MelikhovStechFormFactor::dataBaseOutput(ofstream & output,bool header,
					     bool create) const
{
  if(header){output << "update decayers set parameters=\"";}
  if(create)
    {output << "create Herwig++::MelikhovStechFormFactor " << fullName() << " \n";}
  output << "set " << fullName() << ":ThetaEtaEtaPrime " << _thetaeta  << "\n";
  for(unsigned int ix=0;ix<numberOfFactors();++ix)
    {
      if(ix<initialModes())
	{
	  output << "set " << fullName() << ":FPlus0 " 
		 << ix << "  " << _fplus0[ix] << "\n";
	  output << "set " << fullName() << ":F+sigma_1 " 
		 << ix << "  " << _sigma1fp[ix] << "\n";
	  output << "set " << fullName() << ":F+sigma_2 " 
		 << ix << "  " << _sigma2fp[ix] << "\n";
	  output << "set " << fullName() << ":F00 " 
		 << ix << "  " << _f00[ix] << "\n";
	  output << "set " << fullName() << ":F0sigma_1 " 
		 << ix << "  " << _sigma1f0[ix] << "\n";
	  output << "set " << fullName() << ":F0sigma_2 " 
		 << ix << "  " << _sigma2f0[ix] << "\n";
	  output << "set " << fullName() << ":FT0 " 
		 << ix << "  " << _fT0[ix] << "\n";
	  output << "set " << fullName() << ":FTsigma_1 " 
		 << ix << "  " << _sigma1fT[ix] << "\n";
	  output << "set " << fullName() << ":FTsigma_2 " 
		 << ix << "  " << _sigma2fT[ix] << "\n";
	  output << "set " << fullName() << ":V00 " 
		 << ix << "  " << _V0[ix] << "\n";
	  output << "set " << fullName() << ":V0sigma_1 " 
		 << ix << "  " << _sigma1V0[ix] << "\n";
	  output << "set " << fullName() << ":V0sigma_2 " 
		 << ix << "  " << _sigma2V0[ix] << "\n";
	  output << "set " << fullName() << ":A00 " 
		 << ix << "  " << _A00[ix] << "\n";
	  output << "set " << fullName() << ":A0sigma_1 " 
		 << ix << "  " << _sigma1A0[ix] << "\n";
	  output << "set " << fullName() << ":A0sigma_2 " 
		 << ix << "  " << _sigma2A0[ix] << "\n";
	  output << "set " << fullName() << ":A10 " 
		 << ix << "  " << _A10[ix] << "\n";
	  output << "set " << fullName() << ":A1sigma_1 " 
		 << ix << "  " << _sigma1A1[ix] << "\n";
	  output << "set " << fullName() << ":A1sigma_2 " 
		 << ix << "  " << _sigma2A1[ix] << "\n";
	  output << "set " << fullName() << ":A20 " 
		 << ix << "  " << _A20[ix] << "\n";
	  output << "set " << fullName() << ":A2sigma_1 " 
		 << ix << "  " << _sigma1A2[ix] << "\n";
	  output << "set " << fullName() << ":A2sigma_2 " 
		 << ix << "  " << _sigma2A2[ix] << "\n";
	  output << "set " << fullName() << ":T10 " 
		 << ix << "  " << _T10[ix] << "\n";
	  output << "set " << fullName() << ":T1sigma_1 " 
		 << ix << "  " << _sigma1T1[ix] << "\n";
	  output << "set " << fullName() << ":T1sigma_2 " 
		 << ix << "  " << _sigma2T1[ix] << "\n";
	  output << "set " << fullName() << ":T20 " 
		 << ix << "  " << _T20[ix] << "\n";
	  output << "set " << fullName() << ":T2sigma_1 " 
		 << ix << "  " << _sigma1T2[ix] << "\n";
	  output << "set " << fullName() << ":T2sigma_2 " 
		 << ix << "  " << _sigma2T2[ix] << "\n";
	  output << "set " << fullName() << ":T30 " 
		 << ix << "  " << _T30[ix] << "\n";
	  output << "set " << fullName() << ":T3sigma_1 " 
		 << ix << "  " << _sigma1T3[ix] << "\n";
	  output << "set " << fullName() << ":T3sigma_2 " 
		 << ix << "  " << _sigma2T3[ix] << "\n";
	  output << "set " << fullName() << ":MassP " 
		 << ix << "  " << _massP[ix]/GeV << "\n";
	  output << "set " << fullName() << ":MassV " 
		 << ix << "  " << _massV[ix]/GeV << "\n";
	}
      else
	{
	  output << "insert " << fullName() << ":FPlus0 " 
		 << ix << "  " << _fplus0[ix] << "\n";
	  output << "insert " << fullName() << ":F+sigma_1 " 
		 << ix << "  " << _sigma1fp[ix] << "\n";
	  output << "insert " << fullName() << ":F+sigma_2 " 
		 << ix << "  " << _sigma2fp[ix] << "\n";
	  output << "insert " << fullName() << ":F00 " 
		 << ix << "  " << _f00[ix] << "\n";
	  output << "insert " << fullName() << ":F0sigma_1 " 
		 << ix << "  " << _sigma1f0[ix] << "\n";
	  output << "insert " << fullName() << ":F0sigma_2 " 
		 << ix << "  " << _sigma2f0[ix] << "\n";
	  output << "insert " << fullName() << ":FT0 " 
		 << ix << "  " << _fT0[ix] << "\n";
	  output << "insert " << fullName() << ":FTsigma_1 " 
		 << ix << "  " << _sigma1fT[ix] << "\n";
	  output << "insert " << fullName() << ":FTsigma_2 " 
		 << ix << "  " << _sigma2fT[ix] << "\n";
	  output << "insert " << fullName() << ":V00 " 
		 << ix << "  " << _V0[ix] << "\n";
	  output << "insert " << fullName() << ":V0sigma_1 " 
		 << ix << "  " << _sigma1V0[ix] << "\n";
	  output << "insert " << fullName() << ":V0sigma_2 " 
		 << ix << "  " << _sigma2V0[ix] << "\n";
	  output << "insert " << fullName() << ":A00 " 
		 << ix << "  " << _A00[ix] << "\n";
	  output << "insert " << fullName() << ":A0sigma_1 " 
		 << ix << "  " << _sigma1A0[ix] << "\n";
	  output << "insert " << fullName() << ":A0sigma_2 " 
		 << ix << "  " << _sigma2A0[ix] << "\n";
	  output << "insert " << fullName() << ":A10 " 
		 << ix << "  " << _A10[ix] << "\n";
	  output << "insert " << fullName() << ":A1sigma_1 " 
		 << ix << "  " << _sigma1A1[ix] << "\n";
	  output << "insert " << fullName() << ":A1sigma_2 " 
		 << ix << "  " << _sigma2A1[ix] << "\n";
	  output << "insert " << fullName() << ":A20 " 
		 << ix << "  " << _A20[ix] << "\n";
	  output << "insert " << fullName() << ":A2sigma_1 " 
		 << ix << "  " << _sigma1A2[ix] << "\n";
	  output << "insert " << fullName() << ":A2sigma_2 " 
		 << ix << "  " << _sigma2A2[ix] << "\n";
	  output << "insert " << fullName() << ":T10 " 
		 << ix << "  " << _T10[ix] << "\n";
	  output << "insert " << fullName() << ":T1sigma_1 " 
		 << ix << "  " << _sigma1T1[ix] << "\n";
	  output << "insert " << fullName() << ":T1sigma_2 " 
		 << ix << "  " << _sigma2T1[ix] << "\n";
	  output << "insert " << fullName() << ":T20 " 
		 << ix << "  " << _T20[ix] << "\n";
	  output << "insert " << fullName() << ":T2sigma_1 " 
		 << ix << "  " << _sigma1T2[ix] << "\n";
	  output << "insert " << fullName() << ":T2sigma_2 " 
		 << ix << "  " << _sigma2T2[ix] << "\n";
	  output << "insert " << fullName() << ":T30 " 
		 << ix << "  " << _T30[ix] << "\n";
	  output << "insert " << fullName() << ":T3sigma_1 " 
		 << ix << "  " << _sigma1T3[ix] << "\n";
	  output << "insert " << fullName() << ":T3sigma_2 " 
		 << ix << "  " << _sigma2T3[ix] << "\n";
	  output << "insert " << fullName() << ":MassP " 
		 << ix << "  " << _massP[ix]/GeV << "\n";
	  output << "insert " << fullName() << ":MassV " 
		 << ix << "  " << _massV[ix]/GeV << "\n";
	}
    }
  ScalarFormFactor::dataBaseOutput(output,false,false);
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}
}
