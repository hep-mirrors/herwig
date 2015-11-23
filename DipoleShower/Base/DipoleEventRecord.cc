// -*- C++ -*-
//
// DipoleEventRecord.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DipoleEventRecord class.
//

#include "DipoleEventRecord.h"
#include "Herwig/DipoleShower/Utility/DipolePartonSplitter.h"
#include "Herwig/Shower/ShowerHandler.h"

#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDF/PartonExtractor.h"

#include <boost/utility.hpp>

#include <algorithm>

using namespace Herwig;

PList DipoleEventRecord::colourOrdered() {
  
  PList colour_ordered;

  size_t done_size = outgoing().size();
  if (incoming().first->coloured())
    ++done_size;
  if (incoming().second->coloured())
    ++done_size;

  while (colour_ordered.size() != done_size) {

    PPtr current;

    // start with singlets, as long as we have some

    if (find(colour_ordered.begin(),colour_ordered.end(),incoming().first) ==
	colour_ordered.end() && incoming().first->coloured()) {
      if (!incoming().first->hasColour() || !incoming().first->hasAntiColour())
	current = incoming().first;
    }

    if (!current) {
      for (PList::iterator p = outgoing().begin();
	   p != outgoing().end(); ++p) {
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
      if (find(colour_ordered.begin(),colour_ordered.end(),incoming().second) ==
	  colour_ordered.end() && incoming().second->coloured()) {
	if (!incoming().second->hasColour() || !incoming().second->hasAntiColour())
	  current = incoming().second;
      }
    }

    // then go on with anything else

    if (!current) {
      if (find(colour_ordered.begin(),colour_ordered.end(),incoming().first) ==
	  colour_ordered.end() && incoming().first->coloured()) {
	current = incoming().first;
      }
    }

    if (!current) {
      for (PList::iterator p = outgoing().begin();
	   p != outgoing().end(); ++p) {
	if (find(colour_ordered.begin(),colour_ordered.end(),*p) ==
	    colour_ordered.end() && (**p).coloured()) {
	  current = *p;
	  break;
	}
      }
    }

    if (!current) {
      if (find(colour_ordered.begin(),colour_ordered.end(),incoming().second) ==
	  colour_ordered.end() && incoming().second->coloured()) {
	current = incoming().second;
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
	  if (find(outgoing().begin(),outgoing().end(),*p) != outgoing().end() ||
	      *p == incoming().first ||
	      *p == incoming().second) {
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
	  if (find(outgoing().begin(),outgoing().end(),*p) != outgoing().end() ||
	      *p == incoming().first ||
	      *p == incoming().second) {
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

void DipoleEventRecord::findChains(const PList& ordered) {

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

  bool is33bar =
    ordered.size() == 2 && startIsTriplet && endIsTriplet;

  if (!is33bar) {

    PList::const_iterator theStart = ordered.begin();
    bool onceMore = false;

    for (PList::const_iterator p = ordered.begin();
	 p != ordered.end(); ++p) {

      PList::const_iterator next_it =
	p != --ordered.end() ? boost::next(p) : ordered.begin();

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
	  theChains.push_back(current_chain);
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

      current_chain.dipoles().push_back(Dipole(make_pair(*p,*next_it),pdf,xs));

      if ( onceMore ) {
	next_it = theStart;
	current_chain.check();
	theChains.push_back(current_chain);
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

    current_chain.dipoles().push_back(Dipole(make_pair(ordered.front(),ordered.back()),pdf,xs));

  }

  if (!current_chain.dipoles().empty()) {
    current_chain.check();
    theChains.push_back(current_chain);
  }

}

void DipoleEventRecord::getAll(const ParticleVector& childs,
			       set<PPtr>& hardSet,
			       set<PPtr>& outgoingSet) {

  for ( ParticleVector::const_iterator p = childs.begin();
	p != childs.end(); ++p ) {
    if ( ShowerHandler::currentHandler()->eventHandler()->currentCollision()->isRemnant(*p) )
      continue;
    if ( (**p).children().empty() ) {
      if ( (**p).coloured() &&
	   outgoingSet.find(*p) == outgoingSet.end() )
	outgoingSet.insert(*p);
      else if ( !(**p).coloured() &&
		hardSet.find(*p) == hardSet.end() )
	hardSet.insert(*p);
    } else {
      getAll((**p).children(),hardSet,outgoingSet);
    }
  }

}

// shamelessly stolen from ShowerTree
void DipoleEventRecord::colourIsolate(const vector<PPtr> & original,
				      const vector<PPtr> & copy) {
  // vectors must have same size
  assert(original.size()==copy.size());
  // create a temporary map with all the particles to make looping easier
  vector<PPair> particles;
  particles.reserve(original.size());
  for(unsigned int ix=0;ix<original.size();++ix)
    particles.push_back(make_pair(copy[ix],original[ix]));
  // reset the colour of the copies
  vector<PPair>::const_iterator cit,cjt;
  for(cit=particles.begin();cit!=particles.end();++cit)
    if((*cit).first->colourInfo()) (*cit).first->colourInfo(new_ptr(ColourBase()));
  map<tColinePtr,tColinePtr> cmap;
  // make the colour connections of the copies
  for(cit=particles.begin();cit!=particles.end();++cit) {
    ColinePtr c1,newline;
    // if particle has a colour line
    if((*cit).second->colourLine()&&!(*cit).first->colourLine()) {
      c1=(*cit).second->colourLine();
      newline=ColourLine::create((*cit).first);
      cmap[c1]=newline;
      for(cjt=particles.begin();cjt!=particles.end();++cjt) {
	if(cjt==cit) continue;
	if((*cjt).second->colourLine()==c1)
	  newline->addColoured((*cjt).first);
	else if((*cjt).second->antiColourLine()==c1)
	  newline->addColoured((*cjt).first,true);
      }
    }
    // if anticolour line
    if((*cit).second->antiColourLine()&&!(*cit).first->antiColourLine()) {
      c1=(*cit).second->antiColourLine();
      newline=ColourLine::create((*cit).first,true);
      cmap[c1]=newline;
      for(cjt=particles.begin();cjt!=particles.end();++cjt) {
	if(cjt==cit) continue;
	if((*cjt).second->colourLine()==c1)
	  newline->addColoured((*cjt).first);
	else if((*cjt).second->antiColourLine()==c1)
	  newline->addColoured((*cjt).first,true);
      }
    }
  }
  for ( map<tColinePtr,tColinePtr>::const_iterator c = cmap.begin();
	c != cmap.end(); ++c ) {
    theColourLines[c->second] = c->first;
  }
  // sort out sinks and sources
  for(cit=particles.begin();cit!=particles.end();++cit) {
    tColinePtr cline[2];
    tColinePair cpair;
    for(unsigned int ix=0;ix<4;++ix) {
      cline[0] = ix<2 ? cit->second->colourLine() : cit->second->antiColourLine();
      cline[1] = ix<2 ? cit->first ->colourLine() : cit->first ->antiColourLine();
      if(cline[0]) {
	switch (ix) {
	case 0: case 2:
 	  cpair = cline[0]->sinkNeighbours();
	  break;
	case 1: case 3:
	  cpair = cline[0]->sourceNeighbours();
	  break;
	};
      }
      else {
	cpair = make_pair(tColinePtr(),tColinePtr());
      }
      if(cline[0]&&cpair.first) {
 	map<tColinePtr,tColinePtr>::const_iterator 
	  mit[2] = {cmap.find(cpair.first),cmap.find(cpair.second)};
	if(mit[0]!=cmap.end()&&mit[1]!=cmap.end()) {
	  if(ix==0||ix==2) {
	    cline[1]->setSinkNeighbours(mit[0]->second,mit[1]->second);
	  }
	  else {
	    cline[1]->setSourceNeighbours(mit[0]->second,mit[1]->second);
	  }
	}
      }
    }
  }
}

// shamelessly stolen from ShowerTree
void DipoleEventRecord::updateColour(PPtr particle) {
  // if attached to a colour line
  if(particle->colourLine()) {
    bool reset=false;
    // if colour line from hard process reconnect
    if(theColourLines.find(particle->colourLine())!=theColourLines.end()) {
      ColinePtr c1=particle->colourLine();
      c1->removeColoured(particle);
      theColourLines[c1]->addColoured(particle);
      reset=true;
    }
    // ensure properly connected to the line
    if(!reset) {
      ColinePtr c1=particle->colourLine();
      c1->removeColoured(particle);
      c1->addColoured(particle);
    }
  }
  // if attached to an anticolour line
  if(particle->antiColourLine()) {
    bool reset=false;
    // if anti colour line from hard process reconnect
    if(theColourLines.find(particle->antiColourLine())!=theColourLines.end()) {
      ColinePtr c1=particle->antiColourLine();
      c1->removeColoured(particle,true);
      theColourLines[c1]->addColoured(particle,true);
      reset=true;
    }
    if(!reset) {
      ColinePtr c1=particle->antiColourLine();
      c1->removeColoured(particle,true);
      c1->addColoured(particle,true);
    }
  }
  for ( ParticleVector::const_iterator c = particle->children().begin();
	c != particle->children().end(); ++c ) {
    updateColour(*c);
  }
}


const map<PPtr,PPtr>& 
DipoleEventRecord::prepare(tSubProPtr subpro,
			   tStdXCombPtr xc,
			   const pair<PDF,PDF>& pdf,
			   bool dipoles) {

  theSubProcess = subpro;

  theOutgoing.clear();
  theHard.clear();
  theColourLines.clear();
  theOriginals.clear();

  PPair in = subpro->incoming();

  assert(ShowerHandler::currentHandler());
  tPPair beam = ShowerHandler::currentHandler()->generator()->currentEvent()->incoming();

  // don't take these from the XComb as it may be null
  pair<double,double> xs;
  Direction<0> dir(true);
  xs.first = in.first->momentum().dirPlus()/beam.first->momentum().dirPlus();
  dir.reverse();
  xs.second = in.second->momentum().dirPlus()/beam.second->momentum().dirPlus();

  theXComb = xc;
  thePDFs = pdf;
  theFractions = xs;

  set<PPtr> allHard;
  set<PPtr> allOutgoing;

  getAll(in.first->children(),allHard,allOutgoing);
  getAll(in.second->children(),allHard,allOutgoing);

  vector<PPtr> original;
  vector<PPtr> copies;
  original.push_back(in.first);
  original.push_back(in.second);
  copy(allOutgoing.begin(),allOutgoing.end(),back_inserter(original));
  for ( vector<PPtr>::const_iterator p = original.begin();
	p != original.end(); ++p ) {
    PPtr copy = new_ptr(Particle(**p));
    copies.push_back(copy);
    theOriginals[*p] = copy;
  }

  colourIsolate(original,copies);

  theIncoming.first = copies[0];
  ParticleVector children = theIncoming.first->children();
  for ( ParticleVector::const_iterator c = children.begin();
	c != children.end(); ++c )
    theIncoming.first->abandonChild(*c);
  theIncoming.second = copies[1];
  children = theIncoming.second->children();
  for ( ParticleVector::const_iterator c = children.begin();
	c != children.end(); ++c )
    theIncoming.second->abandonChild(*c);
  copy(copies.begin()+2,copies.end(),back_inserter(theOutgoing));

  for ( set<PPtr>::const_iterator p = allHard.begin(); p != allHard.end(); ++p ) {
    PPtr copy = new_ptr(Particle(**p));
    theHard.push_back(copy);
    theOriginals[*p] = copy;
  }

  if ( dipoles ) {
    PList cordered = colourOrdered();
    findChains(cordered);
  }

  PList::const_iterator XFirst, XLast;

  if ( !theHard.empty() ) {
    XFirst = theHard.begin();
    XLast = theHard.end();
  } else {
    XFirst = theOutgoing.begin();
    XLast = theOutgoing.end();
  }

  thePX = (**XFirst).momentum();
  ++XFirst;
  for ( ; XFirst != XLast; ++XFirst )
    thePX += (**XFirst).momentum();

  return theOriginals;

}

void DipoleEventRecord::clear() {
  theSubProcess = SubProPtr();
  theXComb = StdXCombPtr();
  thePDFs = pair<PDF,PDF>();
  theIncoming = PPair();
  theOutgoing.clear();
  theIntermediates.clear();
  theHard.clear();
  theChains.clear();
  theDoneChains.clear();
  theOriginals.clear();
  theColourLines.clear();
}

void DipoleEventRecord::update(DipoleSplittingInfo& dsplit) {

  if ( incoming().first == dsplit.emitter() ) {
    theIntermediates.push_back(dsplit.emitter());
    theIncoming.first = dsplit.splitEmitter();
    theFractions.first /= dsplit.lastEmitterZ();
  } else if ( incoming().first == dsplit.spectator() ) {
    theIntermediates.push_back(dsplit.spectator());
    theIncoming.first = dsplit.splitSpectator();
    theFractions.first /= dsplit.lastSpectatorZ();    
  }

  if ( incoming().second == dsplit.emitter() ) {
    theIntermediates.push_back(dsplit.emitter());
    theIncoming.second = dsplit.splitEmitter();
    theFractions.second /= dsplit.lastEmitterZ();
  } else if ( incoming().second == dsplit.spectator() ) {
    theIntermediates.push_back(dsplit.spectator());
    theIncoming.second = dsplit.splitSpectator();
    theFractions.second /= dsplit.lastSpectatorZ();    
  }

  PList::iterator pos;

  pos = find(theOutgoing.begin(), theOutgoing.end(), dsplit.emitter());
  if (pos != theOutgoing.end()) {
    theIntermediates.push_back(*pos);
    *pos = dsplit.splitEmitter();
  }

  pos = find(theOutgoing.begin(), theOutgoing.end(), dsplit.spectator());
  if (pos != theOutgoing.end()) {
    theIntermediates.push_back(*pos);
    *pos = dsplit.splitSpectator();
  }

  theOutgoing.push_back(dsplit.emission());

  if (dsplit.splittingKinematics()->doesTransform()) {

    for (PList::iterator p = theIntermediates.begin();
	 p != theIntermediates.end(); ++p) {
      (**p).set5Momentum(dsplit.splittingKinematics()->transform((**p).momentum()));
    }

    for (PList::iterator h = theHard.begin();
	 h != theHard.end(); ++h) {
      (**h).set5Momentum(dsplit.splittingKinematics()->transform((**h).momentum()));
    }

    for (PList::iterator p = theOutgoing.begin();
	 p != theOutgoing.end(); ++p)
      if ((*p) != dsplit.splitEmitter() &&
	  (*p) != dsplit.splitSpectator() &&
	  (*p) != dsplit.emission())
	(**p).set5Momentum(dsplit.splittingKinematics()->transform((**p).momentum()));

  }

}

void
DipoleEventRecord::split(list<Dipole>::iterator dip,
			 list<DipoleChain>::iterator ch,
			 DipoleSplittingInfo& dsplit,
			 pair<list<Dipole>::iterator,list<Dipole>::iterator>& childIterators,
			 DipoleChain*& firstChain, DipoleChain*& secondChain,
			 bool colourSpectator) {

  static DipoleChain empty;

  pair<Dipole,Dipole> children = dip->split(dsplit,colourSpectator);

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
    childIterators.second = secondChain->dipoles().begin();

  }

  if ( !colourSpectator )
    update(dsplit); // otherwise done by recoil(...)

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
	  res.push_back(make_pair(dit,theChains.begin()));
	}
      }
    }

  }

  return res;

}

void DipoleEventRecord::transform(const SpinOneLorentzRotation& rot) {


  Lorentz5Momentum tmp;

  for (PList::iterator p = theIntermediates.begin();
       p != theIntermediates.end(); ++p) {
    tmp = (**p).momentum(); tmp = rot * tmp;
    (**p).set5Momentum(tmp);
  }

  for (PList::iterator h = theHard.begin();
       h != theHard.end(); ++h) {
    tmp = (**h).momentum(); tmp = rot * tmp;
    (**h).set5Momentum(tmp);
  }

  for (PList::iterator p = theOutgoing.begin();
       p != theOutgoing.end(); ++p) {
    tmp = (**p).momentum(); tmp = rot * tmp;
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
  updateColour(incoming().first);
  if ( inParticle != inSubPro )
    inParticle->abandonChild(inSubPro);
  inParton->addChild(inSubPro);
  if ( inParticle != inSubPro )
    inParticle->addChild(incoming().first);
  theIntermediates.push_back(inSubPro);
  theIntermediates.push_back(inParton);

  inSubPro = subProcess()->incoming().second;
  if ( !(inSubPro->parents().empty()) )
    inParticle = inSubPro->parents()[0];
  else
    inParticle = inSubPro;
  inParton = theOriginals[inSubPro];
  theOriginals.erase(inSubPro);
  updateColour(incoming().second);
  if ( inParticle != inSubPro )
    inParticle->abandonChild(inSubPro);
  inParton->addChild(inSubPro);
  if ( inParticle != inSubPro )
    inParticle->addChild(incoming().second);
  theIntermediates.push_back(inSubPro);
  theIntermediates.push_back(inParton);

  while ( !theOriginals.empty() ) {
    PPtr outSubPro = theOriginals.begin()->first;
    PPtr outParton = theOriginals.begin()->second;
    // workaround for OS X Mavericks LLVM libc++
#ifdef _LIBCPP_VERSION
    map<PPtr,PPtr>::const_iterator beg = theOriginals.begin();
#else
    map<PPtr,PPtr>::iterator beg = theOriginals.begin();
#endif
    theOriginals.erase(beg);
    updateColour(outParton);
    outSubPro->addChild(outParton);
    theIntermediates.push_back(outSubPro);
  }

  step->addIntermediates(theIntermediates.begin(),theIntermediates.end());

  for (PList::const_iterator p = theOutgoing.begin();
       p != theOutgoing.end(); ++p)
    step->addDecayProduct(*p);

  for (PList::const_iterator p = theHard.begin();
       p != theHard.end(); ++p)
    step->addDecayProduct(*p);

  if ( firstInteraction &&
       (theIncoming.first->coloured() ||
	theIncoming.second->coloured() ) ) {
      ShowerHandler::currentHandler()->lastExtractor()->newRemnants(theSubProcess->incoming(),theIncoming,step);
  }

  step->addIntermediate(theIncoming.first);
  step->addIntermediate(theIncoming.second);

  return theIncoming;

}

void DipoleEventRecord::debugLastEvent(ostream& os) const {

  bool first = ShowerHandler::currentHandler()->firstInteraction();

  os << "--- DipoleEventRecord ----------------------------------------------------------\n";

  os << " the " << (first ? "hard" : "secondary") << " subprocess is:\n"
     << (*theSubProcess);

  os << " using PDF's " << thePDFs.first.pdf() << " and " 
     << thePDFs.second.pdf() << "\n";

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
