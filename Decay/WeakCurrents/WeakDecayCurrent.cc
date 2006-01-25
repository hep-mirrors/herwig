// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the WeakDecayCurrent class.
//

#include "WeakDecayCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "WeakDecayCurrent.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

WeakDecayCurrent::~WeakDecayCurrent() {}

void WeakDecayCurrent::persistentOutput(PersistentOStream & os) const {
  os << _quark << _antiquark << _numbermodes;
}

void WeakDecayCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _quark >> _antiquark >> _numbermodes;
}

AbstractClassDescription<WeakDecayCurrent> WeakDecayCurrent::initWeakDecayCurrent;
// Definition of the static class description member.

void WeakDecayCurrent::Init() {

  static ClassDocumentation<WeakDecayCurrent> documentation
    ("The WeakDecayCurrent class is the basse class for the"
     " implementation of hadronic currents in weak decays.");

  static ParVector<WeakDecayCurrent,int> interfaceQuark
    ("Quark",
     "The PDG code for the quark.",
     &WeakDecayCurrent::_quark,
     0, 0, 0, 0, 16, false, false, true);

  static ParVector<WeakDecayCurrent,int> interfaceAntiQuark
    ("AntiQuark",
     "The PDG code for the antiquark.",
     &WeakDecayCurrent::_antiquark,
     0, 0, 0, -16, 0, false, false, true);
}

void WeakDecayCurrent::dataBaseOutput(ofstream & output,bool header,bool create) const
{
  if(header){output << "update decayers set parameters=\"";}
  if(create)
    {output << "create Herwig++::WeakDecayCurrent " << fullName() << " \n";}
  for(unsigned int ix=0;ix<_quark.size();++ix)
    {
      if(ix<_numbermodes)
	{
	  output << "set " << fullName() << ":Quark "     
		 << ix << "  " << _quark[ix]     << endl;
	  output << "set " << fullName() << ":AntiQuark " 
		 << ix << "  " << _antiquark[ix] << endl;
	}
      else
	{
	  output << "insert "  << fullName() << ":Quark "     
		 << ix << "  " << _quark[ix]     << endl;
	  output << "insert "  << fullName() << ":AntiQuark " 
		 << ix << "  " << _antiquark[ix] << endl;
	}
    }
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}
}
