// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BFragmentationAnalysisHandler class.
//

#include "BFragmentationAnalysisHandler.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Event.h"
#include "Herwig++/Utilities/StandardSelectors.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "BFragmentationAnalysisHandler.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

BFragmentationAnalysisHandler::~BFragmentationAnalysisHandler() {}

void BFragmentationAnalysisHandler::analyze(tEventPtr event, long,
					    int loop, int state) 
{
  if ( loop > 0 || state != 0 || !event ) return;
  transform(event);

  ///////////////////////////
  // Hadron Level Analysis //
  ///////////////////////////
  // extract the weakly decaying B hadrons using set to avoid double counting
  set<PPtr> allParticles;
  event->select(inserter(allParticles),WeakBHadronSelector());
  // convert to vector
  tPVector particles(allParticles.begin(),allParticles.end());
  // numerator
  _emax = 0.5*generator()->maximumCMEnergy();   
  analyze(particles); 

  ////////////////////////////////////////
  // Parton Level Analysis (e+e-->bbar) //
  ////////////////////////////////////////
  // Get all the particles from the perturbative bit of the events. Find the
  // b and bbar coming straight out of the Z/photon (b_orig, bbar_orig and 
  // ZGamma respectively). Then go and find the b and bbar that have their 
  // maximal off-shellness i.e. just before the first gluon is emitted. 
  // Finally go off down the b quark lines to find the b's just before they
  // hadronize (b_end and bbar_end respectively).
  ParticleSet pert=event->primaryCollision()->step(1)->all();
  analyze_bquarks(pert);
}

LorentzRotation BFragmentationAnalysisHandler::transform(tEventPtr) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void BFragmentationAnalysisHandler::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

void BFragmentationAnalysisHandler::analyze(tPPtr part) 
{
  *_fragBxE  += part->momentum().e()/_emax;
  *_fragBxEa += part->momentum().e()/_emax;
}

void BFragmentationAnalysisHandler::analyze_bquarks(ParticleSet pert)
{
  ParticleSet::const_iterator pit;
  PPtr b_orig,bbar_orig,ZGamma,b_start,bbar_start,b_end,bbar_end;
  // First go through all the particles looking for b's coming out of Z/gamma's:
  for(pit=pert.begin();pit!=pert.end();++pit) {
      PPtr bline = *pit;
      if(abs((*pit)->id())==5) {
	  while(abs(bline->parents()[0]->id())==5) bline=bline->parents()[0];  
      }
      if(bline->id()== 5&&(bline->parents()[0]->id()!=21)) b_orig    = bline;
      if(bline->id()==-5&&(bline->parents()[0]->id()!=21)) bbar_orig = bline;
  }
  if(!b_orig) return;
  if(!bbar_orig) return;
  // Note down the Z/Photon that decays to the b & bbar:
  if(b_orig->parents()[0]==bbar_orig->parents()[0]) 
      ZGamma = b_orig->parents()[0];
  PPtr root_b[] = {b_orig,bbar_orig};
  // Now go and look for the b & bbar just before the first gluon is emitted:
  for(int ix=0;ix<=1;ix++) {
      while(root_b[ix]->momentum().m()<=
	    getParticleData(ParticleID::b)->mass()+1.e-8*GeV) {
	  for(unsigned int jx=0;jx<root_b[ix]->children().size();jx++) 
	      if(root_b[ix]->id()==root_b[ix]->children()[jx]->id()) 
		  root_b[ix]=root_b[ix]->children()[jx]; 
      }
  }
  b_start    = root_b[0];
  bbar_start = root_b[1];
  // Now go and find the b and bbar quarks at the end of the shower before 
  // they turn into hadrons.
  for(unsigned int ix=0;ix<=1;ix++) {
      while(root_b[ix]->momentum().m()>=
	    getParticleData(ParticleID::b)->constituentMass()+1.e-8*GeV) {
	  for(unsigned int jx=0;jx<root_b[ix]->children().size();jx++) 
	      if(root_b[ix]->id()==root_b[ix]->children()[jx]->id()) 
		  root_b[ix]=root_b[ix]->children()[jx]; 
      }
  }
  b_end    = root_b[0];
  bbar_end = root_b[1];

  // Fill the energy fraction histograms with that of the b quarks.
  *_fragbquarkxE  += b_end->momentum().e()/_emax;
  *_fragbquarkxE  += bbar_end->momentum().e()/_emax;
  *_fragbquarkjetmass += b_start->momentum().m()/GeV;
  *_fragbquarkjetmass += bbar_start->momentum().m()/GeV;
}

void BFragmentationAnalysisHandler::persistentOutput(PersistentOStream &) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void BFragmentationAnalysisHandler::persistentInput(PersistentIStream &, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<BFragmentationAnalysisHandler> BFragmentationAnalysisHandler::initBFragmentationAnalysisHandler;
// Definition of the static class description member.

void BFragmentationAnalysisHandler::Init() {

  static ClassDocumentation<BFragmentationAnalysisHandler> documentation
    ("The BFragmentationAnalysisHandler class performs analysis"
     " of the B fragmentation function");
  
}
