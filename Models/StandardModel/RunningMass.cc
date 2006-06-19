// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RunningMass class.
//

#include "RunningMass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Parameter.h"

namespace Herwig {
using namespace ThePEG;

void RunningMass::persistentOutput(PersistentOStream & os) const {
  os << _theQCDOrder << _thePower << _theCoefficient << _theMaxFlav
     << _theStandardModel;
}

void RunningMass::persistentInput(PersistentIStream & is, int) {
  is >> _theQCDOrder >> _thePower >> _theCoefficient >> _theMaxFlav 
     >> _theStandardModel;
}

ClassDescription<RunningMass> RunningMass::initRunningMass;
// Definition of the static class description member.

void RunningMass::Init() {
  
static Reference<RunningMass,StandardModelBase> interfaceStandardModel
  ("StandardModel",
   "Reference to the Standard Model.",
   &RunningMass::_theStandardModel, false, false, true, false); 

  static Parameter<RunningMass,unsigned int> interfaceQCDOrder
    ("QCDOrder",
     "The order in alpha_S",
     &RunningMass::_theQCDOrder, 1, 1, 2,
     false, false, true);
  
  
  static Parameter<RunningMass,unsigned int> interfaceMaxFlav
    ("MaxFlav",
     "maximum number of flavours",
     &RunningMass::_theMaxFlav, 6, 3, 6,
     false, false, true);
  
  static ClassDocumentation<RunningMass> documentation
    ("The RunningMass class is the implementation of the"
     " QCD running mass to one or two loop in alpha_S");
  
}
// Return the masses used.
vector<Energy> RunningMass::mass() const
{
  vector<Energy> masses;
  Energy massf;
  double coeff, as, pi=acos(-1.0);
  for ( unsigned long f = 1; f <= _theMaxFlav; ++f ) {
    PDPtr p = getParticleData(f);
    if(_theQCDOrder==2)
      {coeff=_theCoefficient[f-1]+4./3./pi;}
    else
      {coeff=0.;}
    if ( p ){massf=p->mass();}
    else{massf=0.;}
    as = _theStandardModel->alphaS(massf*massf);
    if(as>0)
      {massf = massf/(1.+coeff*as)/pow(as,_thePower[f-1]);}
    else
      {massf = 0.;}
    masses.push_back(massf);
  }
  return masses;
}
// Return the running mass for a given scale and particle type.
Energy RunningMass::value(Energy2 scale, tcPDPtr part) const
{
  Energy output;
  unsigned int id=abs(part->id());
  // calculate the running mass
  if(id<=_theMaxFlav)
    {
      // calculate the value of alphaS nad number of flavours
      unsigned int nf=  _theStandardModel->Nf(scale);
      double as = _theStandardModel->alphaS(scale);
      id=id-1;
      output = massElement(id)*pow(as,_thePower[nf]);
      if(_theQCDOrder==2){output*=(1.+as*_theCoefficient[nf]);}
    }
  // 
  else
    {output= part->mass();}
  return output;
}
}
