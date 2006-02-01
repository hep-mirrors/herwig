// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the YODA class.
//

#include "YODA.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "YODA.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "FFDipole.h"
#include "IFDipole.h"
#include "GeneralDipole.h"
using namespace Herwig;

YODA::YODA() {}

YODA::YODA(const YODA & x) : DecayRadiationGenerator(x), _ffdipole(x._ffdipole),
			     _ifdipole(x._ifdipole), _generaldipole(x._generaldipole)
{}

YODA::~YODA() {}

void YODA::persistentOutput(PersistentOStream & os) const {
  os << _ffdipole << _ifdipole << _generaldipole;
}

void YODA::persistentInput(PersistentIStream & is, int) {
  is >> _ffdipole >> _ifdipole >> _generaldipole;
}

ClassDescription<YODA> YODA::initYODA;
// Definition of the static class description member.

void YODA::Init() {

  static ClassDocumentation<YODA> documentation
    ("There is no documentation for the YODA class");

  static Reference<YODA,FFDipole> interfaceFFDipole
    ("FFDipole",
     "The final-final dipole",
     &YODA::_ffdipole, false, false, true, false, false);

  static Reference<YODA,IFDipole> interfaceIFDipole
    ("IFDipole",
     "_ifdipole",
     &YODA::_ifdipole, false, false, true, false, false);

  static Reference<YODA,GeneralDipole> interfaceGeneralDipole
    ("GeneralDipole",
     "The general dipole",
     &YODA::_generaldipole, false, false, true, false, false);

}

ParticleVector YODA::generatePhotons(const Particle & p,ParticleVector children)
{
  if(children.size()==2)
    {
      // final-final dipole
      if(p.dataPtr()->iCharge()==0)
	{
	  if(children[0]->dataPtr()->iCharge()!=0&&
	     children[1]->dataPtr()->iCharge()!=0)
	    {
	      //return _generaldipole->generatePhotons(p,children);
	      return _ffdipole->generatePhotons(p,children);
	    }
	  else
	    {return children;}
	}
      // initial final dipole
      else
	{
	  if((children[0]->dataPtr()->iCharge()==0&&
	      children[1]->dataPtr()->iCharge()!=0)||
             (children[0]->dataPtr()->iCharge()!=0&&
	      children[1]->dataPtr()->iCharge()==0)
            )
	    {return _ifdipole->generatePhotons(p,children);}
	  else
	    {return children;}
	}
    }
  else
    {
      return _generaldipole->generatePhotons(p,children);
    }
}
