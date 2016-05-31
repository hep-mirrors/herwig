// -*- C++ -*-
//
// QTildeReconstructor.tcc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined templated member
// functions of the QTildeReconstructor class.
//
using namespace Herwig;
using namespace ThePEG;

namespace {
/**
 *  find showering particle for hard branchings
 */
tShowerParticlePtr SHOWERINGPARTICLE(HardBranchingPtr a) {
  return a->branchingParticle();
}

/**
 *  find showering particle for progenitors
 */
tShowerParticlePtr SHOWERINGPARTICLE(ShowerProgenitorPtr a) {
  return a->progenitor();
}


/**
 * Return colour line progenitor pointer for ShowerProgenitor
 */
template<typename Value>
Ptr<ThePEG::ColourLine>::transient_pointer
CL(Value a, unsigned int index=0) {
  return const_ptr_cast<ThePEG::tColinePtr>(SHOWERINGPARTICLE(a)->colourInfo()->colourLines()[index]);
}

/**
 * Return progenitor colour line size for ShowerProgenitor
 */
template<typename Value>
unsigned int CLSIZE(Value a) {
  return SHOWERINGPARTICLE(a)->colourInfo()->colourLines().size();
}

/**
 * Return anti-colour line progenitor pointer for ShowerProgenitor
 */
template<typename Value>
Ptr<ThePEG::ColourLine>::transient_pointer 
ACL(Value a, unsigned int index=0) {
  return const_ptr_cast<ThePEG::tColinePtr>(SHOWERINGPARTICLE(a)->colourInfo()->antiColourLines()[index]);
}

/**
 * Return progenitor anti-colour line size for ShowerProgenitor
 */
template<typename Value>
unsigned int ACLSIZE(Value a) {
  return SHOWERINGPARTICLE(a)->colourInfo()->antiColourLines().size();
}
}

template<typename Value> void QTildeReconstructor::
findPartners(Value jet,set<Value> & done,
             const set<Value> & jets,
             vector<Value> & system) const {
  tShowerParticlePtr part=SHOWERINGPARTICLE(jet);
  unsigned int partNumColourLines  = part->colourInfo()->    colourLines().size();
  unsigned int partNumAColourLines = part->colourInfo()->antiColourLines().size();
  for(typename set<Value>::const_iterator cit=jets.begin();cit!=jets.end();++cit) {
    if(done.find(*cit)!=done.end()||!SHOWERINGPARTICLE(*cit)->coloured())
      continue;
    bool isPartner = false;
    // one initial one final
    if(part->isFinalState()!=SHOWERINGPARTICLE(*cit)->isFinalState()) {
      //loop over all the colours of both
      for(unsigned int ix=0; ix<partNumColourLines; ++ix) {
	for(unsigned int jx=0; jx<CLSIZE(*cit); ++jx) {
	  if(CL(jet,ix) && CL(jet,ix)==CL(*cit,jx)) {
	    isPartner = true;
	    break;
	  }
	}
	if(isPartner) break;
      }
      if(!isPartner) {
	//loop over anti colours of both
	for(unsigned int ix=0; ix<partNumAColourLines; ++ix) {
	  for(unsigned int jx=0; jx<ACLSIZE(*cit); ++jx) {
	    if(ACL(jet,ix) && ACL(jet,ix)==ACL(*cit,jx)) {
	      isPartner = true;
	      break;
	    }
	  }
	  if(isPartner) break;
	}
      }
    }
    // both in either initial or final state
    else {
      // loop over the colours of the first and the anti-colours of the other
      if(part->colourLine()) {
	for(unsigned int ix=0; ix<partNumColourLines; ++ix) {
	  for(unsigned int jx=0; jx<ACLSIZE(*cit); ++jx) {
	    if(CL(jet,ix) && CL(jet,ix)==ACL(*cit,jx)) {
	      isPartner = true;
	      break;
	    }
	  }
	  if(isPartner) break;
	}
      }
      //loop over the anti-colours of the first and the colours of the other
      if(part->antiColourLine()&&!isPartner) {
	for(unsigned int ix=0; ix<partNumAColourLines; ++ix) {
	  for(unsigned int jx=0; jx<CLSIZE(*cit); jx++) {
	    if(ACL(jet,ix) && ACL(jet,ix)==CL(*cit,jx)) {
	      isPartner = true;
	    }
	  }
	  if(isPartner) break;
	}
      }
    }
    if(isPartner) {
      system.push_back(*cit);
      done.insert(*cit);
      findPartners(*cit,done,jets,system);
      continue;
    }
    // special for sources/sinks
    if(part->colourLine()) {
      if(part->colourLine()->sourceNeighbours().first) {
        tColinePair lines = part->colourLine()->sourceNeighbours();
        if(lines.first ==  CL(*cit)  || lines.first ==  ACL(*cit) ||
           lines.second == CL(*cit)  || lines.second == ACL(*cit) )
          isPartner = true;
      }
      if(part->colourLine()->sinkNeighbours().first) {
        tColinePair lines = part->colourLine()->sinkNeighbours();
        if(lines.first  == CL(*cit)  || lines.first  == ACL(*cit) ||
           lines.second == CL(*cit)  || lines.second == ACL(*cit) )
          isPartner = true;
      }
    }
    if(part->antiColourLine()) {
      if(part->antiColourLine()->sourceNeighbours().first) {
        tColinePair lines = part->antiColourLine()->sourceNeighbours();
        if(lines.first  == CL(*cit)  || lines.first  == ACL(*cit) ||
           lines.second == CL(*cit)  || lines.second == ACL(*cit) )
          isPartner = true;
      }
      if(part->antiColourLine()->sinkNeighbours().first) {
        tColinePair lines = part->antiColourLine()->sinkNeighbours();
        if(lines.first  == CL(*cit)  || lines.first  == ACL(*cit) ||
           lines.second == CL(*cit)  || lines.second == ACL(*cit) )
          isPartner = true;
      }
    }
    if(isPartner) {
      system.push_back(*cit);
      done.insert(*cit);
      findPartners(*cit,done,jets,system);
    }
  }
}

template<typename Value >
typename Herwig::ColourSinglet<Value>::VecType QTildeReconstructor::
identifySystems(set<Value> jets,
		unsigned int & nnun,unsigned int & nnii,unsigned int & nnif,
		unsigned int & nnf ,unsigned int & nni ) const {
  vector<ColourSinglet<Value> > systems;
  set<Value> done;
  for(typename set<Value>::const_iterator it=jets.begin();it!=jets.end();++it) {
    // if not treated create new system
    if(done.find(*it)!=done.end()) continue;
    done.insert(*it);
    systems.push_back(ColourSinglet<Value> (UNDEFINED,*it));
    if(!SHOWERINGPARTICLE(*it)->coloured()) continue;
    findPartners(*it,done,jets,systems.back().jets);
  }
  for(unsigned int ix=0;ix<systems.size();++ix) {
    unsigned int ni(0),nf(0);
    for(unsigned int iy=0;iy<systems[ix].jets.size();++iy) {
      if(SHOWERINGPARTICLE(systems[ix].jets[iy])->isFinalState()) ++nf;
      else                                                        ++ni;
    }
    // type
    // initial-initial
    if(ni==2&&nf==0) {
      systems[ix].type = II;
      ++nnii;
    }
    // initial only
    else if(ni==1&&nf==0) {
      systems[ix].type = I;
      ++nni;
    }
    // initial-final
    else if(ni==1&&nf>0) {
      systems[ix].type = IF;
      ++nnif;
    }
    // final only
    else if(ni==0&&nf>0) {
      systems[ix].type = F;
      ++nnf;
    }
    // otherwise unknown
    else {
      systems[ix].type = UNDEFINED;
      ++nnun;
    }
  }
  return systems;
}

template<typename Value >
void QTildeReconstructor::combineFinalState(vector<ColourSinglet<Value> > & systems) const {
  // check that 1 particle final-state systems which can be combine
  bool canCombine(true);
  for(unsigned int ix=0;ix<systems.size();++ix) {
    if(systems[ix].type!=F) continue;
    if(systems[ix].jets.size()!=1) canCombine = false;
  }
  // return if can't combine
  if(!canCombine) return;
  // otherwise combine them
  vector<ColourSinglet<Value> > oldsystems=systems;
  systems.clear();
  ColourSinglet<Value> finalState;
  finalState.type = F;
  for(unsigned int ix=0;ix<oldsystems.size();++ix) {
    if(oldsystems[ix].type==F) {
      for(unsigned int iy=0;iy<oldsystems[ix].jets.size();++iy)
	finalState.jets.push_back(oldsystems[ix].jets[iy]);
    }
    else
      systems.push_back(oldsystems[ix]);
  }
  systems.push_back(finalState);
}