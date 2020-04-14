// -*- C++ -*-
//
// DipoleEventRecord.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DipoleEventRecord class.
//

#include "DipoleEventRecord.h"
#include "Herwig/Shower/Dipole/DipoleShowerHandler.h"
#include "Herwig/Shower/Dipole/Utility/DipolePartonSplitter.h"
#include "Herwig/Shower/ShowerHandler.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig/Decay/HwDecayerBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDF/PartonExtractor.h"
#include "Herwig/Shower/RealEmissionProcess.h"
#include "Herwig/MatrixElement/Matchbox/Base/MatchboxMEBase.h"
#include "Herwig/MatrixElement/Matchbox/Base/MatchboxAmplitude.h"

#include <boost/utility.hpp>
#include <algorithm>
#include <iterator>

using namespace Herwig;

PList DipoleEventRecord::colourOrdered(PPair & in,
                                       PList & out) {
  
  PList colour_ordered;
  size_t done_size = out.size();
  
  if (in.first->coloured())
    ++done_size;
  if (in.second && in.second->coloured())
    ++done_size;
  
  while (colour_ordered.size() != done_size) {
    
    PPtr current;
    
    // start with singlets, as long as we have some
    
    if (find(colour_ordered.begin(),colour_ordered.end(),in.first) ==
        colour_ordered.end() && in.first->coloured()) {
      if (!in.first->hasColour() || !in.first->hasAntiColour())
	current = in.first;
    }
    
    if (!current) {
      for (PList::iterator p = out.begin();
           p != out.end(); ++p) {
        if (find(colour_ordered.begin(),colour_ordered.end(),*p) ==
            colour_ordered.end() && (**p).coloured()) {
          if (!(**p).hasColour() || !(**p).hasAntiColour()) {
            current = *p;
            break;
          }
        }
      }
    }
    
    if (!current) {
      if (in.second && find(colour_ordered.begin(),colour_ordered.end(),in.second) ==
          colour_ordered.end() && in.second->coloured()) {
        if (!in.second->hasColour() || !in.second->hasAntiColour())
	  current = in.second;
      }
    }
    
    // then go on with anything else
    
    if (!current) {
      if (find(colour_ordered.begin(),colour_ordered.end(),in.first) ==
          colour_ordered.end() && in.first->coloured()) {
        current = in.first;
      }
    }
    
    if (!current) {
      for (PList::iterator p = out.begin();
           p != out.end(); ++p) {
        if (find(colour_ordered.begin(),colour_ordered.end(),*p) ==
            colour_ordered.end() && (**p).coloured()) {
          current = *p;
          break;
        }
      }
    }
    
    if (!current) {
      if (in.second && find(colour_ordered.begin(),colour_ordered.end(),in.second) ==
          colour_ordered.end() && in.second->coloured()) {
        current = in.second;
      }
    }
    
    assert(current);
    
    PPtr next;
    Ptr<ColourLine>::ptr walk_the_line;
    
    while (true) {
      
      if (!walk_the_line) {
        if (current->hasColour()) {
          walk_the_line = current->colourLine();
        }
        else if (current->hasAntiColour()) {
          walk_the_line = current->antiColourLine();
        }
      }
      
      if (!next)
	for (tPVector::const_iterator p = walk_the_line->coloured().begin();
	     p != walk_the_line->coloured().end(); ++p) {
	  if (*p == current)
	    continue;
	  if (find(out.begin(),out.end(),*p) != out.end() ||
	      *p == in.first ||
	      (in.second && *p == in.second)) {
	    next = *p;
	    if (next->hasColour() && next->hasAntiColour()) {
	      walk_the_line = walk_the_line == next->colourLine() ? next->antiColourLine() : next->colourLine();
	    }
	    break;
	  }
	}
      
      if (!next)
	for (tPVector::const_iterator p = walk_the_line->antiColoured().begin();
	     p != walk_the_line->antiColoured().end(); ++p) {
	  if (*p == current)
	    continue;
	  if (find(out.begin(),out.end(),*p) != out.end() ||
	      *p == in.first ||
	      (in.second && *p == in.second)) {
	    next = *p;
	    if (next->hasColour() && next->hasAntiColour()) {
	      walk_the_line = walk_the_line == next->colourLine() ? next->antiColourLine() : next->colourLine();
	    }
	    break;
	  }
	}

      assert(next);

      colour_ordered.push_back(current);
      current = next;

      // done if next is not a gluon or next is already in colour_ordered

      if ((current->hasColour() && !current->hasAntiColour()) ||
	  (!current->hasColour() && current->hasAntiColour())) {
	colour_ordered.push_back(current);
	break;
      }

      if (next->hasColour() && next->hasAntiColour()) {
	if (find(colour_ordered.begin(),colour_ordered.end(),next) != colour_ordered.end())
	  break;
      }

      next = PPtr();

    }

  }  

  return colour_ordered;

}

void DipoleEventRecord::popChain() { 
  assert(!theChains.empty());
  theDoneChains.push_back(DipoleChain());
  theDoneChains.back().dipoles().splice(theDoneChains.back().dipoles().begin(),theChains.front().dipoles());
  theChains.pop_front();
}

void DipoleEventRecord::popChain(list<DipoleChain>::iterator ch) {
  assert(!theChains.empty());
  theDoneChains.push_back(DipoleChain());
  theDoneChains.back().dipoles().splice(theDoneChains.back().dipoles().begin(),ch->dipoles());
  theChains.erase(ch);
}

void DipoleEventRecord::popChains(const list<list<DipoleChain>::iterator>& chs) {

  assert(!theChains.empty());

  for ( list<list<DipoleChain>::iterator>::const_iterator ch =
	  chs.begin(); ch != chs.end(); ++ch ) {
    theDoneChains.push_back(DipoleChain());
    theDoneChains.back().dipoles().splice(theDoneChains.back().dipoles().begin(),(*ch)->dipoles());
  }

  for ( list<list<DipoleChain>::iterator>::const_iterator ch =
	  chs.begin(); ch != chs.end(); ++ch )
    theChains.erase(*ch);

}

DipoleIndex 
DipoleEventRecord::mergeIndex(list<Dipole>::iterator firstDipole, const pair<bool,bool>& whichFirst,
			      list<Dipole>::iterator secondDipole, const pair<bool,bool>& whichSecond) const {
  tcPDPtr emitterData = 
    whichFirst.first ? firstDipole->leftParticle()->dataPtr() : firstDipole->rightParticle()->dataPtr();
  tcPDPtr spectatorData = 
    whichSecond.first ? secondDipole->leftParticle()->dataPtr() : secondDipole->rightParticle()->dataPtr();
  const PDF& emitterPDF =
    whichFirst.first ? firstDipole->leftPDF() : firstDipole->rightPDF();
  const PDF& spectatorPDF =
    whichSecond.first ? secondDipole->leftPDF() : secondDipole->rightPDF();
  return DipoleIndex(emitterData,spectatorData,emitterPDF,spectatorPDF);
}


SubleadingSplittingInfo 
DipoleEventRecord::mergeSplittingInfo(list<DipoleChain>::iterator firstChain, list<Dipole>::iterator firstDipole, 
				      const pair<bool,bool>& whichFirst,
				      list<DipoleChain>::iterator secondChain, list<Dipole>::iterator secondDipole, 
				      const pair<bool,bool>& whichSecond) const {
  SubleadingSplittingInfo res;
  res.index(mergeIndex(firstDipole,whichFirst,secondDipole,whichSecond));
  res.emitter(whichFirst.first ? firstDipole->leftParticle() : firstDipole->rightParticle());
  res.spectator(whichSecond.first ? secondDipole->leftParticle() : secondDipole->rightParticle());
  res.emitterX(whichFirst.first ? firstDipole->leftFraction() : firstDipole->rightFraction());
  res.spectatorX(whichSecond.first ? secondDipole->leftFraction() : secondDipole->rightFraction());
  res.configuration(whichFirst);
  res.spectatorConfiguration(whichSecond);
  res.emitterChain(firstChain);
  res.emitterDipole(firstDipole);
  res.spectatorChain(secondChain);
  res.spectatorDipole(secondDipole);
  return res;
}

void DipoleEventRecord::getSubleadingSplittings(list<SubleadingSplittingInfo>& res) {
  static pair<bool,bool> left(true,false);
  static pair<bool,bool> right(false,true);
  res.clear();
  for ( list<DipoleChain>::iterator cit = theChains.begin();
	cit != theChains.end(); ++cit ) {
    for ( list<Dipole>::iterator dit = cit->dipoles().begin();
	  dit != cit->dipoles().end(); ++dit ) {
      for ( list<Dipole>::iterator djt = dit;
	    djt != cit->dipoles().end(); ++djt ) {
	res.push_back(mergeSplittingInfo(cit,dit,left,cit,djt,left));
	res.push_back(mergeSplittingInfo(cit,dit,right,cit,djt,right));
	if ( dit != djt ) {
	  res.push_back(mergeSplittingInfo(cit,dit,left,cit,djt,right));
	  res.push_back(mergeSplittingInfo(cit,dit,right,cit,djt,left));
	}
      }
    }
    list<DipoleChain>::iterator cjt = cit; ++cjt;
    for ( ; cjt != theChains.end(); ++cjt ) {
      for ( list<Dipole>::iterator dit = cit->dipoles().begin();
	    dit != cit->dipoles().end(); ++dit ) {
	for ( list<Dipole>::iterator djt = cjt->dipoles().begin();
	      djt != cjt->dipoles().end(); ++djt ) {
	  res.push_back(mergeSplittingInfo(cit,dit,left,cjt,djt,left));
	  res.push_back(mergeSplittingInfo(cit,dit,right,cjt,djt,right));
	  res.push_back(mergeSplittingInfo(cit,dit,left,cjt,djt,right));
	  res.push_back(mergeSplittingInfo(cit,dit,right,cjt,djt,left));
	}
      }
    }
  }
}

void DipoleEventRecord::splitSubleading(SubleadingSplittingInfo& dsplit,
					pair<list<Dipole>::iterator,list<Dipole>::iterator>& childIterators,
					DipoleChain*& firstChain, DipoleChain*& secondChain) {
  if ( dsplit.emitterDipole() == dsplit.spectatorDipole() ) {
    assert(dsplit.emitterChain() == dsplit.spectatorChain());
    split(dsplit.emitterDipole(),dsplit.emitterChain(),dsplit,
	  childIterators,firstChain,secondChain,false);
  } else {
    // first need to recoil, then split
    recoil(dsplit.spectatorDipole(),dsplit.spectatorChain(),dsplit);
    split(dsplit.emitterDipole(),dsplit.emitterChain(),dsplit,
	  childIterators,firstChain,secondChain,true);
  }
}

void DipoleEventRecord::findChains(const PList& ordered, 
				   const set<long>& offShellPartons,
				   const bool decay) {
  
  // All uses of findChains should guarantee
  // a non-empty list of particles
  assert( !ordered.empty() );
  
  theChains.clear();
  theDoneChains.clear();

  DipoleChain current_chain;

  // this whole thing needs to have a more elegant implementation at some point

  bool startIsTriplet =
    (ordered.front()->hasColour() && !ordered.front()->hasAntiColour()) ||
    (!ordered.front()->hasColour() && ordered.front()->hasAntiColour());
  bool endIsTriplet =
    (ordered.back()->hasColour() && !ordered.back()->hasAntiColour()) ||
    (!ordered.back()->hasColour() && ordered.back()->hasAntiColour());


  if (!( ordered.size() == 2 && startIsTriplet && endIsTriplet)) {
    
    PList::const_iterator theStart = ordered.begin();
    bool onceMore = false;

    for (PList::const_iterator p = ordered.begin();
	 p != ordered.end(); ++p) {

      PList::const_iterator next_it =
	p != --ordered.end() ? std::next(p) : ordered.begin();

      if (!DipolePartonSplitter::colourConnected(*p,*next_it)) {
	// it may have happened that we need to close the chain due to another
	// chain starting right now; see the above global comment for this fix
	bool startIsOctet =
	  (**theStart).hasColour() && (**theStart).hasAntiColour();
	bool endIsOctet =
	  (**p).hasColour() && (**p).hasAntiColour();
	if ( DipolePartonSplitter::colourConnected(*p,*theStart) &&
	     startIsOctet && endIsOctet ) {
	  swap(next_it,theStart);
	  onceMore = true;
	} else {
          theStart = next_it;
          current_chain.check();
          // Randomize the chains agains biasing of directions.
          if(UseRandom::rndbool()) theChains.push_back(current_chain);
	  else theChains.insert(theChains.begin(),current_chain);
          current_chain.dipoles().clear();
          continue;
	}
      }

      pair<bool,bool> initial_state (false,false);
      initial_state.first = (*p == incoming().first || *p == incoming().second);
      initial_state.second = (*next_it == incoming().first || *next_it == incoming().second);

      pair<int,int> which_in (-1,-1);
      if (initial_state.first)
	which_in.first = *p == incoming().first ? 0 : 1;
      if (initial_state.second)
	which_in.second = *next_it == incoming().first ? 0 : 1;

      pair<double,double> xs (1.,1.);
      if (initial_state.first)
	xs.first = *p == incoming().first ? fractions().first : fractions().second;
      if (initial_state.second)
	xs.second = *next_it == incoming().first ? fractions().first : fractions().second;

      pair<PDF,PDF> pdf;

      if ( which_in.first == 0 )
	pdf.first = pdfs().first;
      else if ( which_in.first == 1 )
	pdf.first = pdfs().second;

      if ( which_in.second == 0 )
	pdf.second = pdfs().first;
      else if ( which_in.second == 1 )
	pdf.second = pdfs().second;
      
      // In the case of a decay process register which
      // parton is incoming to the decay
      pair<bool,bool> decayed_parton (false,false);
      if (decay) {
        decayed_parton.first = (*p == currentDecay()->incoming()[0].first);
        decayed_parton.second = (*next_it == currentDecay()->incoming()[0].first);
      }
     
      // Identify if either parton can have an off-shell mass
      // The first test for partons with zero nominal mass should
      // avoid issues of e.g. non-zero mass gluons
      pair<bool,bool> off_shell (false,false);

      // Note we could do away with the offShellPartons set but,
      // to be safe in the case of an off-shell parton with a mass
      // *very* close to its on-shell mass, we would need to include tests on the
      // offShell indicators in the DipoleIndex == and < operators AND
      // in canHandle and canHandleEquivalent in each massive kernel.
      // Testing these in every splitting will probably be more expensive
      // than doing the following checks for each hard process and decay process

      // Only do off-shell check if the nominal mass is non-zero
      if ( (*p)->nominalMass() != ZERO ) {
	if ( offShellPartons.find(abs((*p)->id())) != offShellPartons.end() )
	  off_shell.first = true;
	else 
	  assert( abs((*p)->mass() - (*p)->nominalMass()) < (*p)->nominalMass()*1.e-5
		  && "There is an off-shell coloured particle in the hard process or a decay"
		  "which  needs to be added to DipoleShowerHandler:OffShellInShower." );
      }
     
      if ( (*next_it)->nominalMass() != ZERO ) {
	if ( offShellPartons.find(abs((*next_it)->id())) != offShellPartons.end() )
	  off_shell.second = true;
	else
	  assert( abs((*next_it)->mass() - (*next_it)->nominalMass())
		  < (*next_it)->nominalMass()*1.e-5
		  && "There is an off-shell coloured particle in the hard process or a decay"
		  "which  needs to be added to DipoleShowerHandler:OffShellInShower." );
      }

      
      current_chain.dipoles().push_back(Dipole({*p,*next_it},pdf,xs,
					       decayed_parton, off_shell));
      
      if ( onceMore ) {
        next_it = theStart;
        current_chain.check();
        // Randomize the chains agains biasing of directions.
        if(UseRandom::rndbool()) theChains.push_back(current_chain);
        else theChains.insert(theChains.begin(),current_chain);
        current_chain.dipoles().clear();
        onceMore = false;
      }
      
    }
  } else {
    
    // treat 2 -> singlet, singlet -> 2 and 1 + singlet -> 1 + singlet special
    // to prevent duplicate dipole

    assert(DipolePartonSplitter::colourConnected(ordered.front(),ordered.back()));

    pair<bool,bool> initial_state (false,false);
    initial_state.first = (ordered.front() == incoming().first || ordered.front() == incoming().second);
    initial_state.second = (ordered.back() == incoming().first || ordered.back() == incoming().second);

    pair<int,int> which_in (-1,-1);
    if (initial_state.first)
      which_in.first = ordered.front() == incoming().first ? 0 : 1;
    if (initial_state.second)
      which_in.second = ordered.back() == incoming().first ? 0 : 1;

    pair<double,double> xs (1.,1.);
    if (initial_state.first)
      xs.first = ordered.front() == incoming().first ? fractions().first : fractions().second;
    if (initial_state.second)
      xs.second = ordered.back() == incoming().first ? fractions().first : fractions().second;

    pair<PDF,PDF> pdf;

    if ( which_in.first == 0 )
      pdf.first = pdfs().first;
    else if ( which_in.first == 1 )
      pdf.first = pdfs().second;

    if ( which_in.second == 0 )
      pdf.second = pdfs().first;
    else if ( which_in.second == 1 )
      pdf.second = pdfs().second;
    
    
    // In the case of a decay process register which
    // parton is incoming to the decay
    pair<bool,bool> decayed_parton (false,false);
    if (decay) {
      decayed_parton.first = (ordered.front() == currentDecay()->incoming()[0].first);
      decayed_parton.second = (ordered.back() == currentDecay()->incoming()[0].first);
    }

    
    // Identify if either parton can have an off-shell mass
    // The first test for partons with zero nominal mass should
    // avoid issues of e.g. non-zero mass gluons
    pair<bool,bool> off_shell (false,false);

    // Only do off-shell check if the nominal mass is non-zero
    if ( ordered.front()->nominalMass() != ZERO ) {
      if ( offShellPartons.find(abs(ordered.front()->id())) != offShellPartons.end() )
	off_shell.first = true;
      else
	assert( abs(ordered.front()->mass() - ordered.front()->nominalMass())
		< ordered.front()->nominalMass()*1.e-5
		&& "There is an off-shell coloured particle in the hard process or a decay"
		"which  needs to be added to DipoleShowerHandler:OffShellInShower." );
    }
   
    if ( ordered.back()->nominalMass() != ZERO ) {
      if ( offShellPartons.find(abs(ordered.back()->id())) != offShellPartons.end() )
	off_shell.second = true;
      else
	assert( abs(ordered.back()->mass() - ordered.back()->nominalMass())
		< ordered.back()->nominalMass()*1.e-5
		&& "There is an off-shell coloured particle in the hard process or a decay"
		"which  needs to be added to DipoleShowerHandler:OffShellInShower." );
    }
        
    current_chain.dipoles().push_back(Dipole({ordered.front(),ordered.back()},
					     pdf,xs, decayed_parton, off_shell));
    
  }

  if (!current_chain.dipoles().empty()) {
    current_chain.check();
    // Randomize the chains agains biasing of directions.
    if(UseRandom::rndbool()) theChains.push_back(current_chain);
    else theChains.insert(theChains.begin(),current_chain);
  }

}

const map<PPtr,PPtr>&
DipoleEventRecord::prepare(tSubProPtr subpro,
                           tStdXCombPtr xc,
			   StepPtr step,
                           const pair<PDF,PDF>& pdf,tPPair beam,
			   bool firstInteraction,
			   const set<long>& offShellPartons,
			   bool dipoles) {
  // set the subprocess
  subProcess(subpro);
  // clear the event record
  outgoing().clear();
  theHard.clear();
  theOriginals.clear();
  theDecays.clear();
  theCurrentDecay = PerturbativeProcessPtr();

  subEmDone = 0;
  if ( doSubleadingNc && firstInteraction ) {

    theDensityOperator.clear();
    continueSubleadingNc = true;
   
    const Ptr<MatchboxXComb>::tptr MBXCombPtr
      = dynamic_ptr_cast<Ptr<MatchboxXComb>::tptr >(xc);
    if ( !MBXCombPtr ) {
      throw Exception() << "Cannot cast StandardXComb as MatchboxXComb. "
			<< "Matchbox is required for "
			<< "colour matrix element corrections."
			<< Exception::runerror; 
    }
    // Set the colour basis if it has not been set
    if ( !theDensityOperator.colourBasis() ) {
      theDensityOperator.colourBasis(MBXCombPtr->matchboxME()->matchboxAmplitude()->colourBasis());
    } else if ( theDensityOperator.colourBasis() != 
	 MBXCombPtr->matchboxME()->matchboxAmplitude()->colourBasis() ) {
      throw Exception() << "The colour basis used in the colour matrix "
			<< "element corrections should not change between events. "
			<< Exception::runerror; 
    }

  } else {
    continueSubleadingNc = false;
  }

  // extract incoming particles
  PPair in = subpro->incoming();
  // get the incoming momentum fractions
  // don't take these from the XComb as it may be null
  pair<double,double> xs;
  ThePEG::Direction<0> dir(true);
  xs.first = in.first->momentum().dirPlus()/beam.first->momentum().dirPlus();
  dir.reverse();
  xs.second = in.second->momentum().dirPlus()/beam.second->momentum().dirPlus();

  xcombPtr(xc);
  pdfs() = pdf;
  fractions() = xs;
  // use ShowerHandler to split up the hard process
  PerturbativeProcessPtr hard;
  DecayProcessMap decay;

  // Special handling for the first interaction:
  // If a post subprocess handler (e.g. QED radiation)
  // is applied, there may be particles in the step object not
  // present in the subprocess object (other than any remnants).
  // These need to be included in any transformations due to
  // II splittings in ::update.
  if ( firstInteraction ) {

    // Initialise a PVector for the outgoing
    tPVector hardProcOutgoing;

    // Include all outgoing particles that are not remnants
    for ( auto & part : step->particles() )
      if ( part->id() != 82 ) {
	hardProcOutgoing.push_back(part);
      }
    ShowerHandler::currentHandler()->splitHardProcess(hardProcOutgoing,
						      hard, decay);
  }

  // For secondary collisions we must use the
  // subProcess object and not the step as the
  // step stores all outgoing from the entire collision
  else
    ShowerHandler::currentHandler()->splitHardProcess(tPVector(subpro->outgoing().begin(),
							       subpro->outgoing().end()),
						      hard,decay);
  
  // vectors for originals and copies of the particles
  vector<PPtr> original;
  vector<PPtr> copies;
  // fill originals
  for(unsigned int ix=0;ix<2;++ix)
    original.push_back(hard->incoming()[ix].first);
  for(unsigned int ix=0;ix<hard->outgoing().size();++ix)
    original.push_back(hard->outgoing()[ix].first);
  for(DecayProcessMap::const_iterator it=decay.begin();it!=decay.end();++it) {
    fillFromDecays(it->second, original);
  }
  // and make copies
  for ( vector<PPtr>::const_iterator p = original.begin();
	p != original.end(); ++p ) {
    PPtr copy = new_ptr(Particle(**p));
    copies.push_back(copy);
    theOriginals[*p] = copy;
  }
  // isolate the colour of the copies from the originals
  colourIsolate(original,copies);
  
  // set the incoming particles
  incoming().first = copies[0];
  ParticleVector children = incoming().first->children();
  for ( ParticleVector::const_iterator c = children.begin();
	c != children.end(); ++c )
    incoming().first->abandonChild(*c);
  incoming().second = copies[1];
  children = incoming().second->children();
  for ( ParticleVector::const_iterator c = children.begin();
	c != children.end(); ++c )
    incoming().second->abandonChild(*c);
  
  // set the outgoing particles for the hard process
  for(unsigned int ix=0;ix<hard->outgoing().size();++ix) {
    if(hard->outgoing()[ix].first->coloured())
      outgoing().push_back(theOriginals[hard->outgoing()[ix].first]);
    else
      theHard.push_back(theOriginals[hard->outgoing()[ix].first]);
  }

  if ( dipoles ) {
    PList cordered = colourOrdered(incoming(),outgoing());
    if ( !cordered.empty() )
      findChains(cordered, offShellPartons, false);
  }
  
  
  
  // sort out the decays
  for(auto const & dec : decay) {
    
    // If the decay particle is in original it needs
    // to be added to the decays and the decay needs to be
    // changed to the copied particles.
    if ( theOriginals.find(dec.second->incoming()[0].first) != theOriginals.end() ) {
      theDecays[theOriginals[dec.second->incoming()[0].first]] = dec.second;
      PerturbativeProcessPtr decayProc = theDecays[theOriginals[dec.second->incoming()[0].first]];
      separateDecay(decayProc);
    }
    
    else {
      assert( find( copies.begin(), copies.end(), dec.second->incoming()[0].first ) != copies.end() );
      theDecays[dec.second->incoming()[0].first] = dec.second;
    }
  }
  
  
  PList::const_iterator XFirst, XLast;

  if ( !theHard.empty() ) {
    XFirst = theHard.begin();
    XLast = theHard.end();
  } else {
    XFirst = outgoing().begin();
    XLast = outgoing().end();
  }

  thePX = (**XFirst).momentum();
  ++XFirst;
  for ( ; XFirst != XLast; ++XFirst )
    thePX += (**XFirst).momentum();
  identifyEventType();


  if ( doSubleadingNc ) {
    theParticlesBefore.clear();
    theParticlesAfter.clear();
    theMomentaAfter.clear();
    theParticleIndices.clear();
    
    // Set the particles and fill the dictionary
    theParticleIndices[incoming().first] = 0;
    theParticleIndices[incoming().second] = 1;
    size_t i = 2;
    theParticlesAfter.reserve(2 + outgoing().size()  + theHard.size());
    theParticlesAfter.push_back(incoming().first->dataPtr());
    theParticlesAfter.push_back(incoming().second->dataPtr());
    theMomentaAfter.push_back(incoming().first->momentum());
    theMomentaAfter.push_back(incoming().second->momentum());
    for ( PList::const_iterator it = outgoing().begin(); it != outgoing().end(); it++ ) {
      theParticlesAfter.push_back((*it)->dataPtr());
      theMomentaAfter.push_back((*it)->momentum());
      theParticleIndices[*it] = i;
      i++;
    }
    // theHard is not added to theParticleIndices, as they aren't needed there
    for ( PList::const_iterator it = theHard.begin(); it != theHard.end(); it++ ) {
      theParticlesAfter.push_back((*it)->dataPtr());
      theMomentaAfter.push_back((*it)->momentum());
    }

    // theParticlesAfter is required for fill
    const Ptr<MatchboxXComb>::tptr MBXCombPtr
      = dynamic_ptr_cast<Ptr<MatchboxXComb>::tptr >(xc);
    if ( !MBXCombPtr ) {
      throw Exception() << "Cannot cast StandardXComb as MatchboxXComb. "
			<< "Matchbox is required for "
			<< "colour matrix element corrections."
			<< Exception::runerror; 
    }
    theDensityOperator.fill(MBXCombPtr,theParticlesAfter,theMomentaAfter);
  }
  
  return theOriginals;
  
}

void DipoleEventRecord::slimprepare(tSubProPtr subpro,
                                    tStdXCombPtr xc,
                                    const pair<PDF,PDF>& pdf,tPPair beam,			 
				    const set<long>& offShellPartons,
                                    bool dipoles) {
  // set the subprocess
  subProcess(subpro);
  // clear the event record
  outgoing().clear();
  theHard.clear();
  theOriginals.clear();
  theDecays.clear();
  theCurrentDecay = PerturbativeProcessPtr();
  // extract incoming particles
  PPair in = subpro->incoming();
  // get the beam
  // get the incoming momentum fractions
  // don't take these from the XComb as it may be null
  pair<double,double> xs;
  ThePEG::Direction<0> dir(true);
  xs.first = in.first->momentum().dirPlus()/beam.first->momentum().dirPlus();
  dir.reverse();
  xs.second = in.second->momentum().dirPlus()/beam.second->momentum().dirPlus();
  xcombPtr(xc);
  
  pdfs() = pdf;
  fractions() = xs;
  incoming()  = in;
  
  for(unsigned int ix=0;ix<subpro->outgoing().size();++ix) {
    if(subpro->outgoing()[ix]->coloured())
      outgoing().push_back(subpro->outgoing()[ix]);
  }
  
  
  if ( dipoles ) {
    PList cordered = colourOrdered(incoming(),outgoing());
    if ( !cordered.empty() )
      findChains(cordered, offShellPartons, false);
  }
  
}


void DipoleEventRecord::fillFromDecays(PerturbativeProcessPtr decayProc, vector<PPtr>& original) {
  
  // Loop over the outgoing of the given perturbative process
  for ( auto const & outIt : decayProc->outgoing() ) {
    
    // Add the outgoing particle to the vector of original particles
    original.push_back(outIt.first);
    
    // Iterate through the outgoing
    if ( outIt.second )
      fillFromDecays( outIt.second, original);
  }
}


void DipoleEventRecord::separateDecay(PerturbativeProcessPtr decayProc) {
  
  // Iteratively replace all entries in the incoming
  // with their copies.
  for ( auto  & inIt : decayProc->incoming() ) {
    
    if ( theOriginals.find( inIt.first ) != theOriginals.end() )
      inIt.first = theOriginals[inIt.first];
  }
  
  // Iteratively replace all entries in the outgoing
  // with their copies.
  for ( auto  & outIt : decayProc->outgoing()) {
    
    if ( theOriginals.count( outIt.first ) )
      outIt.first = theOriginals[outIt.first];
    
    if ( outIt.second )
      separateDecay(outIt.second);
  }
}


void DipoleEventRecord::clear() {
  ShowerEventRecord::clear();
  theDecays.clear();
  theHard.clear();
  theChains.clear();
  theDoneChains.clear();
  theOriginals.clear();
  theDensityOperator.clear();
  theParticlesBefore.clear();
  theParticlesAfter.clear();
  theMomentaAfter.clear();
  theNextDecays.clear();
}



pair<PVector,PVector> DipoleEventRecord::tmpupdate(DipoleSplittingInfo& dsplit) {
  
  PVector inc;
  PVector out;
  
  tcPPtr IF = incoming().first;
  tcPPtr IS = incoming().second;
  
  tcPPtr DE = dsplit.emitter();
  tcPPtr DS = dsplit.spectator();
  
  if ( IF != DE && IF != DS ) {
    PPtr p = IF->data().produceParticle(IF->momentum());
    inc.push_back(p);
  }
  else if ( IF == DE ) inc.push_back( dsplit.splitEmitter() );
  else if ( IF == DS ) inc.push_back( dsplit.splitSpectator() );
  
  if ( IS != DE && IS != DS ) {
    PPtr p = IS->data().produceParticle(IS->momentum());
    inc.push_back(p);
  }
  else if ( IS == DE ) inc.push_back( dsplit.splitEmitter() );
  else if ( IS == DS ) inc.push_back( dsplit.splitSpectator() );
  
  
  if ( IF != DE && IS != DE)
    out.push_back( dsplit.splitEmitter());
  
  if ( IF != DS && IS != DS)
    out.push_back( dsplit.splitSpectator());
  
  out.push_back( dsplit.emission());
  
  for ( tcPPtr h : theHard ){
    PPtr p = h->data().produceParticle(h->momentum());
    if ( dsplit.splittingKinematics()->doesTransform() ) {
      dsplit.splittingKinematics()->transform(p);
    }
    out.push_back(p);
  }
  
  for ( tcPPtr p : outgoing() )
    if ( p != DE &&
         p != DS &&
         p != dsplit.emission() ){
    
      PPtr ou = p->data().produceParticle(p->momentum());;
      if ( dsplit.splittingKinematics()->doesTransform() ){
        dsplit.splittingKinematics()->transform(ou);
      }
      out.push_back(ou);
    }
  return {inc,out};
}



void DipoleEventRecord::update(DipoleSplittingInfo& dsplit) {

  if ( continueSubleadingNc ) {
    subEmDone++;
    theParticlesBefore = theParticlesAfter;
  }
  if ( incoming().first == dsplit.emitter() ) {
    intermediates().push_back(dsplit.emitter());
    incoming().first = dsplit.splitEmitter();
    fractions().first /= dsplit.lastEmitterZ();
    if ( continueSubleadingNc ) {
      theParticleIndices[dsplit.splitEmitter()] = 0;
      theParticlesAfter[0] = dsplit.splitEmitter()->dataPtr();
      theEmitterEmissionIndices.first = 0;
      theEmitterEmissionIndices.second.first = 0;
    }
  } else if ( incoming().first == dsplit.spectator() ) {
    intermediates().push_back(dsplit.spectator());
    incoming().first = dsplit.splitSpectator();
    fractions().first /= dsplit.lastSpectatorZ();
    if ( continueSubleadingNc ) {
      theParticleIndices[dsplit.splitSpectator()] = 0;
      theParticlesAfter[0] = dsplit.splitSpectator()->dataPtr();
      theSpectatorIndices.first = 0;
      theSpectatorIndices.second = 0;
    }
  }

  if ( incoming().second == dsplit.emitter() ) {
    intermediates().push_back(dsplit.emitter());
    incoming().second = dsplit.splitEmitter();
    fractions().second /= dsplit.lastEmitterZ();
    if ( continueSubleadingNc ) {
      theParticleIndices[dsplit.splitEmitter()] = 1;
      theParticlesAfter[1] = dsplit.splitEmitter()->dataPtr();
      theEmitterEmissionIndices.first = 1;
      theEmitterEmissionIndices.second.first = 1;
    }
  } else if ( incoming().second == dsplit.spectator() ) {
    intermediates().push_back(dsplit.spectator());
    incoming().second = dsplit.splitSpectator();
    fractions().second /= dsplit.lastSpectatorZ();    
    if ( continueSubleadingNc ) {
      theParticleIndices[dsplit.splitSpectator()] = 1;
      theParticlesAfter[1] = dsplit.splitSpectator()->dataPtr();
      theSpectatorIndices.first = 1;
      theSpectatorIndices.second = 1;
    }
  }

  PList::iterator pos;

  pos = find(outgoing().begin(), outgoing().end(), dsplit.emitter());
  if (pos != outgoing().end()) {
    intermediates().push_back(*pos);
    *pos = dsplit.splitEmitter();
    if ( continueSubleadingNc ) {
      // The two first elements in theParticlesBefore/After are the incoming
      theEmitterEmissionIndices.first = 2 + distance(outgoing().begin(), pos);
      theEmitterEmissionIndices.second.first = theEmitterEmissionIndices.first;
      theParticlesAfter[theEmitterEmissionIndices.second.first] = dsplit.splitEmitter()->dataPtr();
      theParticleIndices[dsplit.splitEmitter()] = theEmitterEmissionIndices.second.first;
    }
  }

  pos = find(outgoing().begin(), outgoing().end(), dsplit.spectator());
  if (pos != outgoing().end()) {
    intermediates().push_back(*pos);
    *pos = dsplit.splitSpectator();
    if ( continueSubleadingNc ) {
      // The two first elements in theParticlesBefore/After are the incoming
      theSpectatorIndices.first = 2 + distance(outgoing().begin(), pos);
      theSpectatorIndices.second = theSpectatorIndices.first;
      theParticlesAfter[theSpectatorIndices.second] = dsplit.splitSpectator()->dataPtr();
      theParticleIndices[dsplit.splitSpectator()] = theSpectatorIndices.second;
    }
  }

  if ( continueSubleadingNc ) {
    theEmitterEmissionIndices.second.second = 2 + outgoing().size();
    theParticlesAfter.insert(theParticlesAfter.begin()+theEmitterEmissionIndices.second.second,
			     dsplit.emission()->dataPtr());
    theMomentaAfter.insert(theMomentaAfter.begin()+theEmitterEmissionIndices.second.second,
			   dsplit.emission()->momentum());
    theParticleIndices[dsplit.emission()] = theEmitterEmissionIndices.second.second;
  }

  outgoing().push_back(dsplit.emission());

  if (dsplit.splittingKinematics()->doesTransform()) {

    for (PList::iterator h = theHard.begin();
         h != theHard.end(); ++h)
      dsplit.splittingKinematics()->transform(*h);

    for (PList::iterator p = intermediates().begin();
         p != intermediates().end(); ++p) 
      dsplit.splittingKinematics()->transform(*p);

    for (PList::iterator p = outgoing().begin();
         p != outgoing().end(); ++p) {
      if ((*p) != dsplit.splitEmitter() &&
	  (*p) != dsplit.splitSpectator() &&
	  (*p) != dsplit.emission())
        dsplit.splittingKinematics()->transform(*p);	
    }
    
    if ( continueSubleadingNc ) {
      theMomentaAfter[0] = incoming().first->momentum();
      theMomentaAfter[1] = incoming().second->momentum();
      size_t i = 2;
      for (PList::iterator p = outgoing().begin();
	   p != outgoing().end(); p++) {
	theMomentaAfter[i] = (*p)->momentum();
	i++;
      }
      for (PList::iterator p = theHard.begin();
       	   p != theHard.end(); p++) {
       	theMomentaAfter[i] = (*p)->momentum();
       	i++;
      }
    }
    
  } else if ( continueSubleadingNc ) {
    theMomentaAfter[theEmitterEmissionIndices.second.first] = dsplit.splitEmitter()->momentum();//
    theMomentaAfter[theSpectatorIndices.second] = dsplit.splitSpectator()->momentum();//
  }

  // Stop with subleading emissions if the limit has been reached
  if ( doSubleadingNc ) 
    if ( subEmDone == subleadingNcEmissionsLimit ) 
      continueSubleadingNc = false;

  // Handle updates related to decays
  // Showering of decay processes
  // Treat the evolution of the incoming
  // decayed particle as in backward evolution

  if ( dsplit.isDecayProc() ) {
    
    // Create a pointer to the decay process
    PerturbativeProcessPtr decayProc = currentDecay();
    
    // Add the emission to the outgoing of the decay process
    decayProc->outgoing().push_back( {dsplit.emission(), PerturbativeProcessPtr() });
    // Bools to be used throughout
    const bool decayedEmtr = dsplit.index().incomingDecayEmitter();
    const bool decayedSpec = dsplit.index().incomingDecaySpectator();
    
    
    /*
      In the current implementation, **following the hard process**
      all particles in theDecays evolve independently
      e.g. if we have W -> XYZ where all X, Y and Z need to be
      showered and decayed, we only identify them as needing decaying
      (and hence put them in theDecays) AFTER showering the decay of W.
      Hence, XYZ are not even in theDecays until W has been fully
      showered and then they are decayed and showered completely independently
      KEY POINT - Never need to update other entries of theDecays
    
      Note: The PPtr in theDecays should remain unchanged and all changes
      should be made to the relative PerturbativeProcess.
    */
    
    // Splittings from dipoles in the decay process which
    // do not have the decayed parton as emitter or spectator.
    // Update the decay process in theDecays
    if ( !decayedEmtr && !decayedSpec ) {
      
      // Find and replace the old spectator and
      // emitter in the outgoing of the decay process
      bool decayProcEm = false;
      bool decayProcSp = false;
      
      for ( auto & outIt : decayProc->outgoing() ) {
        if ( !decayProcEm && outIt.first == dsplit.emitter() ) {
          outIt = {dsplit.splitEmitter(), PerturbativeProcessPtr()};
          decayProcEm = true;
        }
        
        if ( !decayProcSp && outIt.first == dsplit.spectator() ) {
          outIt = {dsplit.splitSpectator(), PerturbativeProcessPtr() };
          decayProcSp = true;
        }
	
        if ( decayProcEm && decayProcSp )
	  break;
      }
      
      // Test that nothing strange is happening
      assert( (decayProcEm && decayProcSp) );
      
      return;
    }
    
    
    // The spectator is the decayed particle
    else if ( decayedSpec ) {
      
      // Update the dipole event record intermediates
      intermediates().push_back(dsplit.splitSpectator());
      
      // Update the the decayProcess incoming
      decayProc->incoming().clear();
      decayProc->incoming().push_back({dsplit.splitSpectator(),decayProc});
      
      // Update the decay process outgoing
      // Replace the old emitter with the new emitter
      for ( auto & outEmtrIt : decayProc->outgoing() ) {
        if ( outEmtrIt.first == dsplit.emitter() ){
          outEmtrIt = {dsplit.splitEmitter(), PerturbativeProcessPtr() };
          break;
        }
      }
      
      // Perform the recoil transformation
      // Find all particles in the recoil system
      PList recoilSystem;
      for ( auto const & outIt : decayProc->outgoing() ) {
        if ( outIt.first != dsplit.splitEmitter() && outIt.first != dsplit.emission() ) {
          recoilSystem.push_back(outIt.first);
        }
      }
      dsplit.splittingKinematics()->decayRecoil( recoilSystem );
      
      return;
    }
    
    
    // The emitter is the decayed particle
    else  {
      throw Exception()
	<< "DipoleEventRecord: The emitter as a decayed particle is currently not implemented."
	<< Exception::runerror;
      
      assert( currentDecay()->incoming()[0].first == dsplit.emitter() && decayedEmtr && !decayedSpec );
      
      // Update the dipole event record intermediates
      intermediates().push_back(dsplit.splitEmitter());
      
      // Update the the decayProcess incoming
      decayProc->incoming().clear();
      decayProc->incoming().push_back({dsplit.splitEmitter(),decayProc});
      
      // Update the decay process outgoing
      // Replace the old spectator with the new spectator
      for (auto & outSpecIt : decayProc->outgoing() ) {
        if ( outSpecIt.first == dsplit.spectator() ){
          outSpecIt = { dsplit.splitSpectator(), PerturbativeProcessPtr() };
          break;
        }
      }
      
      // Perform the recoil transformation
      assert(dsplit.splittingKinematics()->isDecay());
      // Find all particles in the recoil system
      PList recoilSystem;
      for ( auto const & outIt : decayProc->outgoing() ) {
        if ( outIt.first != dsplit.splitSpectator() && outIt.first != dsplit.emission() ) {
          recoilSystem.push_back(outIt.first);
        }
      }
      dsplit.splittingKinematics()->decayRecoil( recoilSystem );
      
      return;
    }
    
  }

  if ( continueSubleadingNc ) {
    // Fixed alphaS
    double alphaS = 0.118;

    map<pair<size_t,size_t>,Complex> Vijk;
    double Vtemp;
    const Lorentz5Momentum pEmission = dsplit.emission()->momentum();

    // Special cases for the density operator evolution
    // g->qqbar splitting
    bool splitAGluon = (dsplit.emitter()->id() == ParticleID::g) && 
      (dsplit.emission()->id() != ParticleID::g);
    // initial state g->qqbar splitting 
    bool initialGluonSplitting = (dsplit.splitEmitter()->id() == ParticleID::g) && 
      (dsplit.emission()->id() != ParticleID::g);
    if ( initialGluonSplitting )
      assert(dsplit.splitEmitter() == incoming().first
	     || dsplit.splitEmitter() == incoming().second);

    // Set up the dictionary
    std::tuple<size_t,size_t,size_t> tmpTuple;
    map<size_t,size_t> tmpMap;
    size_t n = theEmitterEmissionIndices.second.second;
    theEmissionsMap.clear();
    if ( splitAGluon || initialGluonSplitting ) {
      tmpTuple = std::make_tuple(theEmitterEmissionIndices.first,
				 theEmitterEmissionIndices.second.first,
				 theEmitterEmissionIndices.second.second);
      tmpMap.clear();
      for ( size_t j = 0; j < theParticlesBefore.size(); j++ ) {
	if ( j != theEmitterEmissionIndices.first )
	  tmpMap[j] = j; 
      }
      theEmissionsMap[tmpTuple] = tmpMap;
    } else {
      for ( size_t i = 0; i < theParticlesBefore.size(); i++ ) {
	if ( theParticlesBefore[i]->coloured() ) {
	  tmpTuple = std::make_tuple(i,i,n);
	  tmpMap.clear();
	  for ( size_t j = 0; j < theParticlesBefore.size(); j++ ) {
	    if ( j != i )
	      tmpMap[j] = j; 
	  }
	  theEmissionsMap[tmpTuple] = tmpMap;
	}
      }
    }

    Energy2 pEmitpEmis;
    Energy2 pEmispSpec;

    Lorentz5Momentum pEmitter;
    Lorentz5Momentum pSpectator;

    // Calculate all required dipole factors
    int i,k;
    typedef map<std::tuple<size_t,size_t,size_t>,map<size_t,size_t> > dictMap;
    for(dictMap::const_iterator ijit = theEmissionsMap.begin();
	ijit != theEmissionsMap.end(); ijit++) {
      i = std::get<1>(ijit->first);
      pEmitter = theMomentaAfter[i];
      pEmitpEmis = pEmitter*pEmission;
      for(dictMap::const_iterator kit = theEmissionsMap.begin();
	  kit != theEmissionsMap.end(); kit++) {
	// For gluon splitting ijit == kit
	if ( ijit != kit ) {
	  k = std::get<1>(kit->first);
	  pSpectator = theMomentaAfter[k];
	  pEmispSpec = pEmission*pSpectator;
	  Vtemp = 4*Constants::pi*alphaS*dipoleKernelForEvolution(i, k,
								  pEmitter*pSpectator, pEmitpEmis, 
								  pEmispSpec);
	  Vijk.insert(make_pair(make_pair(i,k),Complex(Vtemp,0.0)));
	} else if ( splitAGluon || initialGluonSplitting ) {
	  k = std::get<1>(kit->first);
	  Vijk.insert(make_pair(make_pair(i,k),Complex(1.0,0.0)));
	}
      }
    }

    theDensityOperator.evolve(Vijk,theParticlesBefore,theParticlesAfter,
			      theEmissionsMap,splitAGluon,initialGluonSplitting);
  }
}

double 
DipoleEventRecord::dipoleKernelForEvolution(size_t em, size_t spec,
					    Energy2 pEmitpSpec, Energy2 pEmitpEmis, 
					    Energy2 pEmispSpec) {
  double Vijk;

  if ( densityOperatorEvolution == 3 ) {
    if ( em == theEmitterEmissionIndices.second.first && 
	 spec == theSpectatorIndices.second ) { 
      Vijk = 1.0;
    } else {
      Vijk = 0.0;
    }
  } else if ( densityOperatorEvolution == 2 ) {
    Vijk = 1.0;
  } else {
    if ( densityOperatorEvolution == 0 ) {
      if ( pEmitpEmis < densityOperatorCutoff ) pEmitpEmis = densityOperatorCutoff; 
      if ( pEmispSpec < densityOperatorCutoff ) pEmispSpec = densityOperatorCutoff;
    }
    Vijk = ((pEmitpSpec)/GeV2)/((pEmitpEmis/GeV2)*
				(pEmispSpec/GeV2));
  }
  
  return Vijk;
}

void
DipoleEventRecord::split(list<Dipole>::iterator dip,
                         list<DipoleChain>::iterator ch,
                         DipoleSplittingInfo& dsplit,
                         pair<list<Dipole>::iterator,list<Dipole>::iterator>& childIterators,
                         DipoleChain*& firstChain, DipoleChain*& secondChain,
                         bool colourSpectator) {
  
  static DipoleChain empty;

  pair<Dipole,Dipole> children = dip->split(dsplit,colourSpectator,
					    continueSubleadingNc);

  list<Dipole>::iterator breakup =
    ch->insertSplitting(dip,children,childIterators);

  if ( breakup == ch->dipoles().end() ) {
    firstChain = &(*ch);
    secondChain = &empty;
  } else {

    DipoleChain other;
    other.dipoles().splice(other.dipoles().end(),ch->dipoles(),breakup,ch->dipoles().end());

    chains().push_back(other);
    firstChain = &(*ch);
    secondChain = &(chains().back());

    // explicitly fix iterators in case the splice implementation
    // at hand does invalidate iterators (the SGI docu says, it doesn't,
    // but it seems that this behaviour is not part of the standard)
    childIterators.first = --firstChain->dipoles().end();
    childIterators.second = secondChain->dipoles().begin();  }
  
  if ( !colourSpectator ) {
    update(dsplit); // otherwise done by recoil(...)
    
  }
}

 

pair<PVector,PVector> DipoleEventRecord::tmpsplit(list<Dipole>::iterator dip,
                                                  list<DipoleChain>::iterator ,
                                                  DipoleSplittingInfo& dsplit,
                                                  pair<list<Dipole>::iterator,list<Dipole>::iterator>& ,
                                                  DipoleChain*& , DipoleChain*& ,
                                                  bool colourSpectator) {
  
  
  dip->tmpsplit(dsplit,colourSpectator);
  return tmpupdate(dsplit); // otherwise done by recoil(...)
  
}

void DipoleEventRecord::recoil(list<Dipole>::iterator dip,
                               list<DipoleChain>::iterator ch,
                               DipoleSplittingInfo& dsplit) {
  
  dip->recoil(dsplit);
  ch->updateDipole(dip);

  update(dsplit);

}

list<pair<list<Dipole>::iterator,list<DipoleChain>::iterator> >
DipoleEventRecord::inDipoles() {

  list<pair<list<Dipole>::iterator,list<DipoleChain>::iterator> > res;

  for ( list<DipoleChain>::iterator chit = theDoneChains.begin();
	chit != theDoneChains.end(); ++chit ) {

    bool haveOne = false;

    for ( list<Dipole>::iterator dit = chit->dipoles().begin();
	  dit != chit->dipoles().end(); ++dit ) {
      if ( dit->leftPDF().pdf() || dit->rightPDF().pdf() ) {
	haveOne = true;
	break;
      }
    }

    if ( haveOne ) {
      theChains.splice(theChains.begin(),theDoneChains,chit);
      for ( list<Dipole>::iterator dit = theChains.front().dipoles().begin();
	    dit != theChains.front().dipoles().end(); ++dit ) {
        if ( dit->leftPDF().pdf() || dit->rightPDF().pdf() ) {
          res.push_back({dit,theChains.begin()});
        }
      }
    }

  }

  return res;

}

void DipoleEventRecord::transform(const LorentzRotation& rot) {


  Lorentz5Momentum tmp;

  for (PList::iterator p = intermediates().begin();
       p != intermediates().end(); ++p) {
    tmp = (**p).momentum();
    if ( (*p)->spinInfo() )
      (*p)->spinInfo()->transform(tmp, rot);
    tmp = rot * tmp;
    (**p).set5Momentum(tmp);
  }

  for (PList::iterator h = theHard.begin();
       h != theHard.end(); ++h) {
    tmp = (**h).momentum();
    if ( (*h)->spinInfo() )
      (*h)->spinInfo()->transform(tmp, rot);
    tmp = rot * tmp;
    (**h).set5Momentum(tmp);
  }

  for (PList::iterator p = outgoing().begin();
       p != outgoing().end(); ++p) {
    tmp = (**p).momentum();
    if ( (*p)->spinInfo() )
      (*p)->spinInfo()->transform(tmp, rot);
    tmp = rot * tmp;
    (**p).set5Momentum(tmp);
  }

}

tPPair DipoleEventRecord::fillEventRecord(StepPtr step, bool firstInteraction, bool) {

  PPtr inSubPro = subProcess()->incoming().first;
  PPtr inParticle;
  
  if ( !(inSubPro->parents().empty()) )
    inParticle = inSubPro->parents()[0];
  else
    inParticle = inSubPro;
  
  PPtr inParton = theOriginals[inSubPro];
  theOriginals.erase(inSubPro);
  updateColour(incoming().first,true);
  
  if ( inParticle != inSubPro )
    inParticle->abandonChild(inSubPro);
  inParton->addChild(inSubPro);
  if ( inParticle != inSubPro )
    inParticle->addChild(incoming().first);
  intermediates().push_back(inSubPro);
  intermediates().push_back(inParton);
  
  // Repeat all the above for the second incoming particle
  inSubPro = subProcess()->incoming().second;
  if ( !(inSubPro->parents().empty()) )
    inParticle = inSubPro->parents()[0];
  else
    inParticle = inSubPro;
  inParton = theOriginals[inSubPro];
  theOriginals.erase(inSubPro);
  updateColour(incoming().second,true);
  if ( inParticle != inSubPro )
    inParticle->abandonChild(inSubPro);
  inParton->addChild(inSubPro);
  if ( inParticle != inSubPro )
    inParticle->addChild(incoming().second);
  intermediates().push_back(inSubPro);
  intermediates().push_back(inParton);
  
  // theOriginals is populated in ::prepare and contains all of the incoming and outgoing particles of the original hard process
  // Here outgoing particles from theOriginals are added into the intermediates()
  while ( !theOriginals.empty() ) {
    PPtr outSubPro = theOriginals.begin()->first;
    PPtr outParton = theOriginals.begin()->second;
    outSubPro->setLifeLength(Lorentz5Distance());
    outSubPro->setVertex(LorentzPoint());
    outParton->setLifeLength(Lorentz5Distance());
    outParton->setVertex(LorentzPoint());
    // fix for displacemet of onstable particles
    Energy width = outParton->dataPtr()->generateWidth(outParton->mass());
    if ( width > ZERO ) {
      Time lifetime = outParton->dataPtr()->generateLifeTime(outParton->mass(), width);
      Lorentz5Distance lLength;
      lLength.setTau(lifetime);
      lLength.setVect(outParton->momentum().vect()*(lifetime / max(outParton->mass(), Constants::epsilon*GeV)));
      lLength.rescaleEnergy();
      outParton->setLifeLength(lLength);
    }
    // workaround for OS X Mavericks LLVM libc++
#ifdef _LIBCPP_VERSION
    map<PPtr,PPtr>::const_iterator beg = theOriginals.begin();
#else
    map<PPtr,PPtr>::iterator beg = theOriginals.begin();
#endif
    theOriginals.erase(beg);
    updateColour(outParton,true);
    outSubPro->addChild(outParton);
    intermediates().push_back(outSubPro);
  }
  
  // Update the intermediates of the step
  step->addIntermediates(intermediates().begin(),intermediates().end());
  
  for (auto const & p : outgoing())
    step->addDecayProduct( p );
  
  for (auto const & p : theHard)
    step->addDecayProduct( p );
  
  if ( firstInteraction &&
       (incoming().first->coloured() ||
	incoming().second->coloured() ) ) {
    ShowerHandler::currentHandler()->lastExtractor()
      ->newRemnants(subProcess()->incoming(),incoming(),step);
  }
  
  step->addIntermediate(incoming().first);
  step->addIntermediate(incoming().second);
  
  return incoming();
  
}


 bool DipoleEventRecord::prepareDecay( PerturbativeProcessPtr decayProc,
				       const set<long>& offShellPartons ) {
  
  // Create objects containing the incoming and outgoing partons,
  // required as inputs for colourOrdered.
  PList out;
  for( auto const & dec : decayProc->outgoing()) {
    if(dec.first->coloured()) {
      out.push_back(dec.first);
    }
  }
  
  // Only need to shower if we have coloured outgoing particles
  if ( out.empty() )
    return false;
  
  else {
    // For the incoming, use a PPair containing the incoming and a null pointer
    PPair in;
    in.first = decayProc->incoming()[0].first;
    

    // Chains are found later if the subleading shower is used
    if ( !doSubleadingNc ) {
      // Create an ordered list of particles
      PList cordered;
      cordered = colourOrdered(in,out);
    
      // Find the dipole chains for this decay
      findChains(cordered,offShellPartons,true);
    }    
    return true;
  }
}

Energy DipoleEventRecord::decay(PPtr incoming, bool& powhegEmission) {
  // get the process
  PerturbativeProcessPtr process = theDecays[incoming];
  assert(process);
  //tDMPtr decayMode = new_ptr(DecayMode());
  tDMPtr decayMode = DMPtr();
  // Do not decay particles that have already been decayed
  // Note the herwig decayer deals with colour connections
  if ( process->outgoing().empty() ) {
    process->incoming()[0].first = incoming;
    DecayProcessMap decay;
    // Decay the particle, returning a pointer to the decay mode
    decayMode = ShowerHandler::currentHandler()->decay(process,decay,true);
  }
  
  
  // Sort out the colour connections of particles already decayed
  else {
    // sort out the colour of the incoming
    map<tColinePtr,tColinePtr> cmap;
    if(incoming->colourLine())
      cmap[process->incoming()[0].first->colourLine()] = incoming->colourLine();
    if(incoming->antiColourLine())
      cmap[process->incoming()[0].first->antiColourLine()] = incoming->antiColourLine();
    // fix colours of outgoing
    for(auto const & outg : process->outgoing()) {
      map<tColinePtr,tColinePtr>::iterator it =
	cmap.find(outg.first->colourLine());
      if(it!=cmap.end()) {
        ColinePtr c1=outg.first->colourLine();
        c1->removeColoured(outg.first);
        it->second->addColoured(outg.first);
      }
      it = cmap.find(outg.first->antiColourLine());
      if(it!=cmap.end()) {
        ColinePtr c1=outg.first->antiColourLine();
        c1->removeAntiColoured(outg.first);
        it->second->addAntiColoured(outg.first);
      }
    }
    // swap the incoming
    process->incoming()[0].first = incoming;
  }
  
  // Set the scale of all particles involved in the decay process to the
  // mass of the decaying particle
  
  // Initialise the scale for the evolution of
  // the parton shower following the decay
  Energy showerScale = ZERO;
  
  // Set the scale for the evolution of the shower
  showerScale = process->incoming()[0].first->momentum().m();
  
  Energy2 decayScaleSqr = sqr( showerScale );
  process->incoming()[0].first->scale( decayScaleSqr );
  
  for(auto & outg : process->outgoing()) {
    outg.first->scale( decayScaleSqr );
  }
  
  // Update the decaying particle in the process and the event
  PList::iterator posOut = find(outgoing().begin(), outgoing().end(), incoming);
  PList::iterator posHard = find(hard().begin(), hard().end(), incoming);
  assert((posOut!=outgoing().end() && posHard==hard().end()) ||
         (posOut==outgoing().end() && posHard!=hard().end()) );
  
  
  if ( posOut!=outgoing().end() ) {
    outgoing().erase(posOut);
  }
  
  else {
    hard().erase(posHard);
  }
  intermediates().push_back(process->incoming()[0].first);
  
  // Populate the children of the incoming
  for(auto const & outg : process->outgoing()) {
    PPtr outgoing = outg.first;
    process->incoming()[0].first->addChild(outgoing);
  }
  
  
  // If a decayed particle is not decayed above,
  // e.g. a W in a 3-body top decay, find its decaymode.
  if ( powhegEmission && !decayMode ) {
    
    string tag = incoming->dataPtr()->name() + "->";
    
    // Must use OrderedParticles for a tag search
    ShowerHandler::OrderedParticles decayOut;
    for(auto const & outg : process->outgoing()) {
      decayOut.insert(outg.first->dataPtr());
    }
    
    // Construct the tag
    for(auto const & dec : decayOut) {
      if( dec!=*decayOut.begin() ) tag += ",";
      tag +=dec->name();
    }
    tag += ";";
    
    // Find the decay mode
    decayMode = ShowerHandler::currentHandler()->findDecayMode(tag);
  }
  
  
  
  // Perform the powheg emission
  if ( powhegEmission ) {

    if ( decayMode ) {

      HwDecayerBasePtr decayer;
      decayer = dynamic_ptr_cast<HwDecayerBasePtr>(decayMode->decayer());

      if ( decayer->hasPOWHEGCorrection() ) {

	// Construct a real emission process and populate its
	// incoming and outcoming prior to any powheg emission
	RealEmissionProcessPtr born = new_ptr( RealEmissionProcess() );
	born->bornIncoming().push_back( incoming );

	for(auto const & outg : process->outgoing()) {
	  born->bornOutgoing().push_back(outg.first);
	}
        
	// Generate any powheg emission, returning 'real'
        RealEmissionProcessPtr real = decayer->generateHardest( born );
        
	// If an emission has been attempted
	// (Note if the emission fails, a null ptr is returned)
        if ( real ) {
	  
          showerScale = real->pT()[ShowerInteraction::QCD];

	  // If an emission is generated sort out the particles
	  if ( !real->outgoing().empty() ) {
	    
	    // Update the decay process
	    // Note: Do not use the new incoming particle
	    PPtr oldEmitter;
	    PPtr newEmitter;
	    
	    // Use the name recoiler to avoid confusion with
	    // the spectator in the POWHEGDecayer
	    // i.e. the recoiler can be coloured or non-coloured
          PPtr oldRecoiler;
          PPtr newRecoiler;
          
          if ( real->emitter() == 1 ) {
            oldEmitter = real->bornOutgoing()[0];
            oldRecoiler = real->bornOutgoing()[1];
            newEmitter = real->outgoing()[0];
            newRecoiler = real->outgoing()[1];
          }
          else if ( real->emitter() == 2) {
            oldEmitter = real->bornOutgoing()[1];
            oldRecoiler = real->bornOutgoing()[0];
            newEmitter = real->outgoing()[1];
            newRecoiler = real->outgoing()[0];
          }
          
          PPtr emitted = real->outgoing()[ real->emitted()-1];
          
	  // Update the scales
          newRecoiler->scale(oldRecoiler->scale());
          
          newEmitter->scale(sqr(showerScale));
          emitted->scale(sqr(showerScale));
          
	  // Update the colour flow of the new outgoing particles
	  // Note the emitted and newEmitter are already colour
	  // connected by the powheg emission function
          emitted->incomingColour(oldEmitter, oldEmitter->id()<0);
          
          if ( newRecoiler->coloured() )
	    newRecoiler->incomingColour(oldRecoiler, oldRecoiler->id()<0);
          
	  // Update the children of the outgoing
          oldRecoiler->addChild( newRecoiler );
          oldEmitter->addChild( newEmitter );
          oldEmitter->addChild( emitted );
          
          
	  // Note: The particles in the pert proc outgoing and both outgoing
	  // vectors of the real emission proc are in the same order
          for(unsigned int ix=0;ix<real->bornOutgoing().size();++ix) {
            
	    // Update the decay process
            assert(process->outgoing()[ix].first == real->bornOutgoing()[ix]);
            process->outgoing()[ix].first = real->outgoing()[ix];
            
	    // Add the outgoing from the born
	    // decay to the event intermediates
            intermediates().push_back(real->bornOutgoing()[ix]);
          }
          
	  // Add the emitted to the outgoing of the decay process
          process->outgoing().push_back( { emitted, PerturbativeProcessPtr() } );
	  }


	  // Else, if no emission above pTmin, set particle scales
	  else {
	    for(auto & outg : process->outgoing()) {
	      outg.first->scale( sqr(showerScale) );
	    }
	    powhegEmission = false;
	  }

	}
	
	// No powheg emission occurred:
        else
	  powhegEmission = false;
        
      }
      
      // No powheg emission occurred:
      else
	powhegEmission = false;
    }
    
    // No powheg emission occurred:
    else
      powhegEmission = false;
  }
  
  // Copy the outgoing from the decay
  // process to the event record
  for(auto const & outg : process->outgoing()) {
    if ( outg.first->coloured() )
      outgoing().push_back(outg.first);
    else
      hard().push_back(outg.first);
  }
  
  return showerScale;
}

void DipoleEventRecord::updateDecayMom( PPtr decayParent, PerturbativeProcessPtr decayProc ) {
  
  // Only particles that have already been decayed
  // should be passed to this function
  assert( !(decayProc->outgoing().empty()) );
  
  // Create a list of the children to update their momenta
  PList children;
  for ( auto const & outg : decayProc->outgoing() ) {
    children.push_back( outg.first );
  }
  
  // Boost the children
  PList::iterator beginChildren = children.begin();
  PList::iterator endChildren = children.end();

  const Momentum3 transformMom = decayParent->momentum().vect();
  Lorentz5Momentum sum = ThePEG::UtilityBase::sumMomentum(beginChildren, endChildren);
  LorentzRotation rot = ThePEG::UtilityBase::transformToCMS(sum);
  rot = ThePEG::UtilityBase::transformFromCMS
    (Lorentz5Momentum(transformMom, sqrt(transformMom.mag2() + sum.m2()))) * rot;
      
  // Must transform the spinInfo using the momentum prior to transforming
  for ( const auto& p : children ) {
    if ( p->spinInfo() )
      p->spinInfo()->transform(p->momentum(),rot);
  }
      
  ThePEG::UtilityBase::transform(beginChildren, endChildren, rot );

}


void DipoleEventRecord::updateDecayChainMom( PPtr decayParent, PerturbativeProcessPtr decayProc ) {
  
  // Note - this updates the momenta of the
  // outgoing of the given decay process
  
  // Update the momenta of the outgoing from this decay
  updateDecayMom( decayParent, decayProc );
  
  // Iteratively update the momenta of the rest of the decay chain
  for ( auto & outg : decayProc->outgoing() ) {

    // If a child has a corresponding pert proc
    // then it has decay products
    if ( outg.second ) {

      for ( auto & dec : theDecays ) {
        if ( dec.second == outg.second ) {

          // If the particle has spininfo
          if ( dec.first->spinInfo() ) {
            
            // Copied from DipoleVertexRecord::updateSpinInfo,
            // would be better to use a common function
            // Update any spin information
            const Lorentz5Momentum& oldMom = dec.first->momentum();
            const Lorentz5Momentum& newMom = outg.first->momentum();
            
            // Rotation from old momentum to +ve z-axis
            LorentzRotation oldToZAxis;
            Axis axisOld(oldMom.vect().unit());
            if( axisOld.perp2() > 1e-12 ) {
              double sinth(sqrt(1.-sqr(axisOld.z())));
              oldToZAxis.rotate( -acos(axisOld.z()),Axis(-axisOld.y()/sinth,axisOld.x()/sinth,0.));
            }
            
            // Rotation from new momentum to +ve z-axis
            LorentzRotation newToZAxis;
            Axis axisNew(newMom.vect().unit());
            if( axisNew.perp2() > 1e-12 ) {
              double sinth(sqrt(1.-sqr(axisNew.z())));
              newToZAxis.rotate( -acos(axisNew.z()),Axis(-axisNew.y()/sinth,axisNew.x()/sinth,0.));
            }
            
            // Boost from old momentum to new momentum along z-axis
            Lorentz5Momentum momOldRotated = oldToZAxis*Lorentz5Momentum(oldMom);
            Lorentz5Momentum momNewRotated = newToZAxis*Lorentz5Momentum(newMom);
            
            Energy2 a = sqr(momOldRotated.z()) + sqr(momNewRotated.t());
            Energy2 b = 2.*momOldRotated.t()*momOldRotated.z();
            Energy2 c = sqr(momOldRotated.t()) - sqr(momNewRotated.t());
            double beta;
            
            // The rotated momentum should always lie along the +ve z-axis
            if ( momOldRotated.z() > ZERO )
              beta = (-b + sqrt(sqr(b)-4.*a*c)) / 2. / a;
            else
              beta = (-b - sqrt(sqr(b)-4.*a*c)) / 2. / a;
            
            LorentzRotation boostOldToNew(0., 0., beta);
            
            // Total transform
            LorentzRotation transform = (newToZAxis.inverse())*boostOldToNew*oldToZAxis;

            // Transform spin info and mom
            dec.first->spinInfo()->transform(oldMom, transform);
          }

          dec.first->setMomentum(outg.first->momentum());          
          break;
        }
      }
      
      // Iteratively update any decay products
      if ( !outg.second->outgoing().empty() )
	updateDecayChainMom( outg.first, outg.second );
    }
  }
}


void DipoleEventRecord::updateDecays(PerturbativeProcessPtr decayProc, bool iterate) {
  
  // Note - This does not update the momenta of the outgoing
  // of decayProc. 
  // i.e. it is for use following the (non-)showering 
  // of a decay when the daughter momentum are correct.
  // With iterate = true, this updates the rest of the decay chain.

  // Update the list of next decays
  if ( decayProc == theCurrentDecay && !theNextDecays.empty() ) {
    assert( theNextDecays.back() == decayProc->incoming()[0].first );
    theNextDecays.pop_back();
  }
  
  // Loop over the outgoing from this decay
  for ( auto & outg : decayProc->outgoing() ) {
        if ( outg.second && !outg.second->outgoing().empty() ) {
      // Outgoing particles which have already been decayed
      PPtr newDecayed = outg.first;
      PerturbativeProcessPtr newDecayProc = outg.second;
      
      // Update the outgoing momenta from this decay
      updateDecayMom( newDecayed, newDecayProc);
      
      // If this decay is already in theDecays then erase it
      for ( auto const & dec : theDecays ) {
        if(dec.second==newDecayProc) {
          theDecays.erase(dec.first);
          break;
        }
      }
      // Add to theDecays
      theDecays[newDecayed] = newDecayProc;

      // Update the list of next decays
      if ( decayProc == theCurrentDecay )
        theNextDecays.push_back(newDecayed);
	
      // Iteratively update theDecays from the decay chain
      if ( iterate ) 
	updateDecays( newDecayProc );
      
	}
    
    // Deal with any outgoing which need to be decayed
    else if ( ShowerHandler::currentHandler()->decaysInShower(outg.first->id()) ) {
      PerturbativeProcessPtr newDecay=new_ptr(PerturbativeProcess());
      newDecay->incoming().push_back({ outg.first , decayProc } );
      theDecays[outg.first] = newDecay;
      
      // Update the list of next decays
      if ( decayProc )
        theNextDecays.push_back(outg.first);
      
    }
  }
  
}


void DipoleEventRecord::debugLastEvent(ostream& os) const {
  
  bool first = ShowerHandler::currentHandler()->firstInteraction();
  
  os << "--- DipoleEventRecord ----------------------------------------------------------\n";
  
  os << " the " << (first ? "hard" : "secondary") << " subprocess is:\n"
     << (*subProcess());

  os << " using PDF's " << pdfs().first.pdf() << " and " 
     << pdfs().second.pdf() << "\n";

  os << " chains showering currently:\n";

  for ( list<DipoleChain>::const_iterator chit = theChains.begin();
	chit != theChains.end(); ++chit )
    os << (*chit);

  os << " chains which finished showering:\n";

  for ( list<DipoleChain>::const_iterator chit = theDoneChains.begin();
	chit != theDoneChains.end(); ++chit )
    os << (*chit);

  os << "--------------------------------------------------------------------------------\n";

  os << flush;

}
