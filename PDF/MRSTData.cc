#include "MRSTData.h"

using namespace ThePEG;
using namespace Herwig;

MRSTData::MRSTData() : Interfaced() {}
MRSTData::MRSTData(const MRSTData &x) : Interfaced(x) {}
MRSTData::~MRSTData() {}

AbstractNoPIOClassDescription<MRSTData> MRSTData::initMRSTData;

void MRSTData::persistentOutput(PersistentOStream &) const {}
void MRSTData::peristentInput(PersistentIStream &, int) {}
void MRSTData::Init() {}
void MRSTData::doinit() throw(InitException) { Interfaced::doinit(); }
void MRSTData::dofinish() { Interfaced::dofinish(); }
void MRSTData::doupdate() throw(UpdateException) { Interfaced::doupdate(); }
void MRSTData::rebind(const TranslationMap &t) throw(RebindException) { 
  Interfaced::rebind(t); 
}
IVector MRSTData::getReferences() { return Interfaced::getReferences(); }
