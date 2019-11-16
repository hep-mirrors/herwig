// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEMultiChannel class.
//

#include "MEMultiChannel.h"
#include "Herwig/Decay/DecayIntegrator.fh"
#include "Herwig/Decay/PhaseSpaceMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Cuts/Cuts.h"

using namespace Herwig;

MEMultiChannel::~MEMultiChannel() {}

void MEMultiChannel::getDiagrams() const {
  int ndiag=0;
  channelMap_.clear();
  for(PhaseSpaceModePtr mode : modes_) {
    channelMap_.push_back(map<int,int>());
    for(unsigned int ix=0;ix<mode->channels().size();++ix) {
      ThePEG::Ptr<ThePEG::Tree2toNDiagram>::pointer diag = mode->channels()[ix].createDiagram();
       ndiag+=1;
      diag = new_ptr((*diag,-ndiag));
      channelMap_.back()[ndiag]= ix;
      add(diag);
    }
  }
}

Energy2 MEMultiChannel::scale() const {
  return sHat();
}

int MEMultiChannel::nDim() const {
  return 0;
  // return modes_[0]->nRand();
}

bool MEMultiChannel::generateKinematics(const double * r) {
  // first find the mode
  int imode = -1;
  for(unsigned int ix=0;ix<modes_.size();++ix) {
    bool matched=true;
    for(unsigned int iy=0;iy<modes_[ix]->outgoing().size();++iy) {
      if(mePartonData()[iy+2]!=modes_[ix]->outgoing()[iy]) {
	matched=false;
	break;
      }
    }
    if(matched) {
      imode=ix;
      break;
    }
  }
  assert(imode>=0);
  // fill the stack of random numbers
  modes_[imode]->fillStack();
  //modes_[imode]->fillStack(r);
  // generate the momenta
  int ichan;
  vector<Lorentz5Momentum> momenta(meMomenta().size()-2);
  Lorentz5Momentum pcm=meMomenta()[0]+meMomenta()[1];
  Energy wgt = modes_[imode]->weight(ichan,pcm,momenta);
  // set the jacobian
  jacobian(wgt/sqrt(sHat()));
  // and momenta
  for(unsigned int ix=0;ix<momenta.size();++ix)
    meMomenta()[ix+2]=momenta[ix];
  // check the cuts
  tcPDVector tout(momenta.size());
  vector<LorentzMomentum> out(meMomenta().size()-2);
  for(unsigned int ix=2;ix<meMomenta().size();++ix) {
    tout[ix-2] = mePartonData()[ix];
    out[ix-2]=meMomenta()[ix];
  }
  if ( !lastCuts().passCuts(tout, out, mePartonData()[0], mePartonData()[1]) )
    return false;
  // passed cuts
  return true;
}

CrossSection MEMultiChannel::dSigHatDR() const {
  return sqr(hbarc)*me2()*jacobian()/sHat();
}

Selector<MEBase::DiagramIndex>
MEMultiChannel::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    double wgt = me2(channelMap_[iMode_][-diags[i]->id()] );
    sel.insert(wgt, i);
  }
  return sel;
}

Selector<const ColourLines *>
MEMultiChannel::colourGeometries(tcDiagPtr ) const {
  static ColourLines none("");
  Selector<const ColourLines *> sel;
  sel.insert(1.0, &none);
  return sel;
}

void MEMultiChannel::persistentOutput(PersistentOStream & os) const {
  os << modes_ << channelMap_;
}

void MEMultiChannel::persistentInput(PersistentIStream & is, int) {
  is >> modes_ >> channelMap_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<MEMultiChannel,MEBase>
  describeHerwigMEMultiChannel("Herwig::MEMultiChannel", "Herwig.so");

void MEMultiChannel::Init() {

  static ClassDocumentation<MEMultiChannel> documentation
    ("There is no documentation for the MEMultiChannel class");

}


void MEMultiChannel::doinit() {
  MEBase::doinit();
for(tPhaseSpaceModePtr mode : modes_)
  mode->init();
}

void MEMultiChannel::doinitrun() {
  MEBase::doinitrun();
for(tPhaseSpaceModePtr mode : modes_)
  mode->initrun();
}
