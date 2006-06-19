// -*- C++ -*-

#include "ShowerTree.h"
#include "Herwig++/Shower/Kinematics/ShowerParticle.h"
#include "ThePEG/Utilities/Timer.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Handlers/XComb.h"
#include "ShowerProgenitor.h"

#include <cassert>

namespace Herwig {
  
using namespace ThePEG;

// constructor from hard process
ShowerTree::ShowerTree(tEHPtr eh, 
		       const ParticleVector & out,ShowerVarsPtr vars,
		       multimap<Energy,ShowerTreePtr>& decay) 
  : _wasHard(true),_parent(ShowerTreePtr()), _showerVariables(vars), _hasShowered(false)
{
  Timer<1400> timer("ShowerTree::ShowerTree::hard");
  assert(eh);
  PPtr in1  =eh->lastPartons().first;
  PPtr in2  =eh->lastPartons().second;
  double x1 =eh->lastX1();
  double x2 =eh->lastX2();
  // must have two incoming particles
  assert(in1 && in2);
  // set the parent tree
  _parent=ShowerTreePtr();
  // tempory vectors to contain all the particles before insertion into
  // the data structure
  vector<PPtr> original,copy;
  vector<ShowerParticlePtr> shower;
  // create copies of ThePEG particles for the incoming particles
  original.push_back(in1);copy.push_back(new_ptr(Particle(*in1)));
  original.push_back(in2);copy.push_back(new_ptr(Particle(*in2)));
  // and same for outgoing
  map<PPtr,ShowerTreePtr> trees;
  for (ParticleVector::const_iterator it = out.begin();it != out.end(); ++it)
    {
      // if decayed or should be decayed in shower make the tree
      PPtr orig=*it;
      if(!orig->children().empty()||
	 _showerVariables->decayInShower(orig->id()))
	{
	  ShowerTreePtr newtree=new_ptr(ShowerTree(orig,_showerVariables,
						   decay,eh));
	  newtree->setParents();
	  trees.insert(make_pair(orig,newtree));
 	  Energy width=orig->dataPtr()->generateWidth(orig->mass());
 	  decay.insert(make_pair(width,newtree));
	}
      original.push_back(orig);
      copy.push_back(new_ptr(Particle(*orig)));
    }
  // colour isolate the hard process
  colourIsolate(original,copy);
  // now create the Shower particles
  // create ShowerParticles for the incoming particles
  assert(original.size() == copy.size());
  for(unsigned int ix=0;ix<original.size();++ix)
    {
      ShowerParticlePtr temp=new_ptr(ShowerParticle(*copy[ix],1));
      temp->setFinalState(ix>=2);
      // incoming
      if(ix<2)
	{
	  if(ix==0)      temp->x(x1);
	  else if(ix==1) temp->x(x2);
	  _incomingLines.insert(make_pair(new_ptr(ShowerProgenitor(original[ix],
								   copy[ix],temp)),
					  temp));
	}
      // outgoing
      else
	{
	  _outgoingLines.insert(make_pair(new_ptr(ShowerProgenitor(original[ix],
								   copy[ix],temp)),
					  temp));
	}
    }
  // set up the map of daughter trees
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator mit;
  for(mit=_outgoingLines.begin();mit!=_outgoingLines.end();++mit)
    {
      map<PPtr,ShowerTreePtr>::const_iterator tit=trees.find(mit->first->original());
      if(tit!=trees.end())
	_treelinks.insert(make_pair(tit->second,mit->first));
    }
}

ShowerTree::ShowerTree(PPtr in,ShowerVarsPtr vars,
		       multimap<Energy,ShowerTreePtr>& decay,
		       tEHPtr ch)
    : _wasHard(false), _showerVariables(vars), _hasShowered(false)
{
  Timer<1401> timer("ShowerTree::ShowerTree::decay");
  // there must be an incoming particle
  assert(in);
  // tempory vectors to contain all the particles before insertion into
  // the data structure
  vector<PPtr> original,copy;
  // insert place holder for incoming particle
  original.push_back(in);copy.push_back(PPtr());
  // we need to deal with the decay products if decayed
  map<PPtr,ShowerTreePtr> trees;
  if(!in->children().empty())
    {
      ParticleVector children=in->children();
      for(unsigned int ix=0;ix<children.size();++ix)
 	{
	  // if decayed or should be decayed in shower make the tree
	  PPtr orig=children[ix];
 	  in->abandonChild(orig);
	  // does not work if child is in hard process
	  ch->currentEvent()->removeParticle(orig);
	  if(!orig->children().empty()||
	     _showerVariables->decayInShower(orig->id()))
	    {
	      ShowerTreePtr newtree=new_ptr(ShowerTree(orig,_showerVariables,decay,ch));
	      trees.insert(make_pair(orig,newtree));
	      Energy width=orig->dataPtr()->generateWidth(orig->mass());
	      decay.insert(make_pair(width,newtree));
	    }
	  original.push_back(orig);
	  copy.push_back(new_ptr(Particle(*orig)));
	}
    }
  // create the incoming particle
  copy[0]     = new_ptr(Particle(*in));
  // create the parent
  ShowerParticlePtr sparent(new_ptr(ShowerParticle(*copy[0],2)));
  // make the new children if needed
//   if(copy.size()>1)
//     {
//       for(unsigned int ix=1;ix<copy.size();++ix)
// 	{copy[0]->addChild(copy[ix]);}
//     }
  // isolate the colour
  colourIsolate(original,copy);
  sparent->setFinalState(false);
  _incomingLines.insert(make_pair(new_ptr(ShowerProgenitor(original[0],copy[0],sparent))
 				  ,sparent));
  // return if not decayed
  if(original.size()==1) return;
  // create the children
  assert(copy.size() == original.size());
  for (unsigned int ix=1;ix<original.size();++ix)
    {
      ShowerParticlePtr stemp= new_ptr(ShowerParticle(*copy[ix],2));
      stemp->setFinalState(true);
      _outgoingLines.insert(make_pair(new_ptr(ShowerProgenitor(original[ix],copy[ix],
 							       stemp)),
 				      stemp));
    } 
  // set up the map of daughter trees
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator mit;
  for(mit=_outgoingLines.begin();mit!=_outgoingLines.end();++mit)
    {
      map<PPtr,ShowerTreePtr>::const_iterator tit=trees.find(mit->first->original());
      if(tit!=trees.end())
	_treelinks.insert(make_pair(tit->second,mit->first));
    }
}


void ShowerTree::updateFinalStateShowerProduct(ShowerProgenitorPtr progenitor,
					       ShowerParticlePtr parent,
					       const ShowerParticleVector & children)
{
  assert(children.size()==2);
  bool matches[2];
  for(unsigned int ix=0;ix<2;++ix)
    {matches[ix]=children[ix]->id()==progenitor->id();}
  ShowerParticlePtr newpart;
  if(matches[0]&&matches[1])
    {
      double z[2]={children[0]->sudAlpha(),children[1]->sudAlpha()};
      if(z[0]>z[1]) newpart=children[0];
      else          newpart=children[1];
    }
  else if(matches[0]) newpart=children[0];
  else if(matches[1]) newpart=children[1];
  else                newpart=ShowerParticlePtr();
  _outgoingLines[progenitor]=newpart;
}

void ShowerTree::updateInitialStateShowerProduct(ShowerProgenitorPtr progenitor,
						 ShowerParticlePtr oldParent,
						 ShowerParticlePtr newParent)
{
  if(oldParent->id()==newParent->id())
    _incomingLines[progenitor]=newParent;
  else
    _incomingLines[progenitor]=ShowerParticlePtr();
}

void ShowerTree::colourIsolate(vector<PPtr> original,vector<PPtr> copy)
{
  // vectors must have same size
  assert(original.size()==copy.size());
  // create a temporary map with all the particles to make looping easier
  map<PPtr,PPtr> particles;
  for(unsigned int ix=0;ix<original.size();++ix)
    {particles.insert(make_pair(copy[ix],original[ix]));}
  // reset the colour of the copies
  map<PPtr,PPtr>::const_iterator cit,cjt;
  for(cit=particles.begin();cit!=particles.end();++cit)
    if((*cit).first->colourInfo()) (*cit).first->colourInfo(new_ptr(ColourBase()));
  // make the colour connections of the copies
  for(cit=particles.begin();cit!=particles.end();++cit)
    {
      // if particle has a colour line
      if((*cit).second->colourLine()&&!(*cit).first->colourLine())
	{
	  ColinePtr c1=(*cit).second->colourLine();
	  ColinePtr newline=ColourLine::create((*cit).first);
 	  for(cjt=particles.begin();cjt!=particles.end();++cjt)
 	    {
	      if(cjt==cit) continue;
 	      if((*cjt).second->colourLine()==c1)
		newline->addColoured((*cjt).first);
 	      else if((*cjt).second->antiColourLine()==c1)
		newline->addColoured((*cjt).first,true);
	    }
	}
      // if anticolour line
      if((*cit).second->antiColourLine()&&!(*cit).first->antiColourLine())
	{
	  ColinePtr c1=(*cit).second->antiColourLine();
	  ColinePtr newline=ColourLine::create((*cit).first,true);
 	  for(cjt=particles.begin();cjt!=particles.end();++cjt)
 	    {
	      if(cjt==cit) continue;
 	      if((*cjt).second->colourLine()==c1)
		newline->addColoured((*cjt).first);
 	      else if((*cjt).second->antiColourLine()==c1)
		newline->addColoured((*cjt).first,true);
 	    }
	}
    }
}

void ShowerTree::insertHard(StepPtr pstep,bool ISR, bool FSR)
{
  assert(_incomingLines.size()==2);
  _colour.clear();
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cit;
  // construct the map of colour lines for hard process
  for(cit=_incomingLines.begin();cit!=_incomingLines.end();++cit)
    {
      if(!cit->first->perturbative()) continue; 
      if((*cit).first->copy()->colourLine()) 
	_colour.insert(make_pair((*cit).first->copy()->colourLine(),
				 (*cit).first->original()->colourLine()));
      if((*cit).first->copy()->antiColourLine())
	_colour.insert(make_pair((*cit).first->copy()->antiColourLine(),
				 (*cit).first->original()->antiColourLine()));
    }
  for(cit=_outgoingLines.begin();cit!=_outgoingLines.end();++cit)
    {
      if(!cit->first->perturbative()) continue;
      if((*cit).first->copy()->colourLine()) 
	_colour.insert(make_pair((*cit).first->copy()->colourLine(),
				 (*cit).first->original()->colourLine()));
      if((*cit).first->copy()->antiColourLine())
	_colour.insert(make_pair((*cit).first->copy()->antiColourLine(),
				 (*cit).first->original()->antiColourLine()));
    }
  // initial-state radiation
  if(ISR)
    {
      for(cit=incomingLines().begin();cit!=incomingLines().end();++cit)
 	{
 	  ShowerParticlePtr init=(*cit).first->progenitor();
 	  if(!init->getThePEGBase()) 
 	    throw Exception() << "Initial-state particle must have a ThePEGBase"
 			      << " in ShowerTree::fillEventRecord()" 
 			      << Exception::runerror;
	  // if not from a matrix element correction
	  if(cit->first->perturbative())
	    {
	      PPtr original = (*cit).first->original();
	      PPtr hadron = original->parents()[0];
	      PPtr intermediate = original->children()[0];
	      // break mother/daugther relations
	      init->abandonChild(intermediate);
	      init->addChild(original);
	      hadron->abandonChild(original);
	      // if particle showers add shower
	      if(hadron->children().size() > 1) 
		addInitialStateShower(init,pstep,false);
	      // no showering for this particle
	      else 
		{
		  updateColour(init);
		  hadron->abandonChild(init);
		  hadron->addChild(init);
		  pstep->addIntermediate(init);
		}
	    }
	  else
	    {
 	      PPtr original = (*cit).first->original();
 	      PPtr hadron = original->parents()[0];
 	      PPtr intermediate = original->children()[0];
	      PPtr copy=cit->first->copy();
 	      // break mother/daugther relations
 	      init->abandonChild(intermediate);
	      copy->abandonChild(intermediate);
	      copy->addChild(original);
	      updateColour(copy);
 	      hadron->abandonChild(original);
	      hadron->abandonChild(copy);
 	      init->addChild(copy);
	      pstep->addIntermediate(copy);
 	      // if particle showers add shower
 	      if(hadron->children().size() > 1) 
 		addInitialStateShower(init,pstep,false);
	      // no showering for this particle
	      else 
		{
		  updateColour(init);
		  hadron->abandonChild(init);
		  hadron->addChild(init);
		  pstep->addIntermediate(init);
		}
	    }
	}
     }
   // final-state radiation
  PPair incoming=pstep->collision()->primarySubProcess()->incoming();
   if(FSR)
     {
       for(cit=outgoingLines().begin();cit!=outgoingLines().end();++cit)
 	{
 	  ShowerParticlePtr init=(*cit).first->progenitor();
 	  if(!init->getThePEGBase()) 
 	    throw Exception() << "Final-state particle must have a ThePEGBase"
 			      << " in ShowerTree::fillEventRecord()" 
 			      << Exception::runerror;
	  // if not from a matrix element correction
	  if(cit->first->perturbative())
	    {
	      // register the shower particle as a 
	      // copy of the one from the hard process
	      tParticleVector parents=init->parents();
	      for(unsigned int ix=0;ix<parents.size();++ix)
		{parents[ix]->abandonChild(init);}
	      (*cit).first->original()->addChild(init);
	      pstep->addDecayProduct(init);
	    }
	  // from a matrix element correction
	  else
	    {
	      if(cit->first->original()==incoming.first||
		 cit->first->original()==incoming.second)
		{
		  updateColour((*cit).first->copy());
 		  (*cit).first->original()->parents()[0]->
		    addChild((*cit).first->copy());
 		  pstep->addDecayProduct((*cit).first->copy());
 		  (*cit).first->copy()->addChild(init);
 		  pstep->addDecayProduct(init);
		}
	      else
		{
		  updateColour((*cit).first->copy());
		  (*cit).first->original()->addChild((*cit).first->copy());
		  pstep->addDecayProduct((*cit).first->copy());
		  (*cit).first->copy()->addChild(init);
		  pstep->addDecayProduct(init);
		}
	    }
	  updateColour(init);
	  // insert shower products
	  addFinalStateShower(init,pstep);
	}
     }
  _colour.clear();
}

void ShowerTree::addFinalStateShower(PPtr p, StepPtr s) {
  if(!p->children().empty()) 
    {
      ParticleVector::const_iterator child;
      for(child=p->children().begin(); child != p->children().end(); ++child) 
	{
	  updateColour(*child);
	  s->addDecayProduct(*child);
	  addFinalStateShower(*child,s);
	}
    }
}

void ShowerTree::updateColour(PPtr particle)
{
  // if attached to a colour line
  if(particle->colourLine())
    {
      // if colour line from hard process reconnect
      if(_colour.find(particle->colourLine())!=_colour.end())
	{
	  ColinePtr c1=particle->colourLine();
	  c1->removeColoured(particle);
	  _colour[c1]->addColoured(particle);
	}
    }
  // if attached to an anticolour line
  if(particle->antiColourLine())
    {
      // if anti colour line from hard process reconnect
      if(_colour.find(particle->antiColourLine())!=_colour.end())
	{
	  ColinePtr c1=particle->antiColourLine();
	  c1->removeColoured(particle,true);
	  _colour[c1]->addColoured(particle,true);
	}
    }
}

void ShowerTree::addInitialStateShower(PPtr p, StepPtr s, bool addchildren) {
  // Each parton here should only have one parent
  if(!p->parents().empty()) 
    {
      if(p->parents().size()!=1) 
	throw Exception() << "Particle must only have one parent in ShowerTree"
			  << "::addInitialStateShower" << Exception::runerror;
      addInitialStateShower(p->parents()[0],s);
    }
  updateColour(p);
  ParticleVector::const_iterator child;
  // if not adding children return 
  if(!addchildren) return;
  // add children
  for(child = p->children().begin(); child != p->children().end(); ++child) 
    {
      // if a final-state particle update the colour
      ShowerParticlePtr schild = 
	dynamic_ptr_cast<ShowerParticlePtr>(*child);
      if(schild && schild->isFinalState()) updateColour(*child);
      // if there are grandchildren of p
      if(!(*child)->children().empty()) 
	{
	  // Add child as intermediate
	  s->addIntermediate(*child);
	  // If child is shower particle and final-state, add children
	  if(schild && schild->isFinalState()) addFinalStateShower(schild,s);
	} 
      else 
	s->addDecayProduct(*child);
    }
}

void ShowerTree::decay(multimap<Energy,ShowerTreePtr> & decay,
		       tEHPtr ch)
{
  // must be one incoming particle
  assert(_incomingLines.size()==1);
  // if not already decayed decay it
  if(_outgoingLines.empty())
    {
      // now we need to replace the particle with a new copy after the shower
      // find particle after the shower
      ShowerProgenitorPtr pthis=_parent->_treelinks[this];
      assert(pthis);
      ShowerParticlePtr newparent=_parent->_outgoingLines[pthis];
      // now make the new progenitor
      vector<PPtr> original,copy;
      original.push_back(newparent);
      copy.push_back(new_ptr(Particle(*newparent)));
      // reisolate the colour
      colourIsolate(original,copy);
      // make the new progenitor
      ShowerParticlePtr stemp=new_ptr(ShowerParticle(*copy[0],2));
      stemp->setFinalState(false);
      ShowerProgenitorPtr newprog=new_ptr(ShowerProgenitor(original[0],copy[0],stemp));
      _incomingLines.clear();
      _incomingLines.insert(make_pair(newprog,stemp));
      // now we need to decay the copy
      PPtr parent=copy[0];
      unsigned int ntry = 0;
       while (true)
 	{
 	  // exit if fails
 	  if (++ntry>=200)
 	    throw Exception() << "Failed to perform decay in ShowerTree::decay()"
 			      << " after " << 200
 			      << " attempts for " << parent->PDGName() 
 			      << Exception::eventerror;
 	  // select decay mode
 	  tDMPtr dm(parent->data().selectMode(*parent));
 	  if(!dm) 
 	    throw Exception() << "Failed to select decay  mode in ShowerTree::decay()"
 			      << Exception::eventerror;
 	  if(!dm->decayer()) 
 	    throw Exception() << "No Decayer for selected decay mode "
 			      << " in ShowerTree::decay()"
 			      << Exception::runerror;
 	  // start of try block
 	  try 
 	    {
 	      ParticleVector children = dm->decayer()->decay(*dm, *parent);
 	      // if no children have another go
 	      if(children.empty()) continue;
 	      // set up parent
 	      parent->decayMode(dm);
 	      // add children
	      for (unsigned int i = 0, N = children.size(); i < N; ++i )
 		{
 		  children[i]->setLabVertex(parent->labDecayVertex());
 		  parent->addChild(children[i]);
 		  parent->scale(0.0*GeV2);
 		}
 	      // if succeeded break out of loop
 	      break;
 	    }
 	  catch(Veto) {}
	}
       // insert the trees from the children
       ParticleVector children=parent->children();
       map<PPtr,ShowerTreePtr> trees;
       for(unsigned int ix=0;ix<children.size();++ix)
  	{
	  PPtr orig=children[ix];
  	  parent->abandonChild(orig);
	  // if particle has children or decays in shower
	  if(!orig->children().empty()||
	     _showerVariables->decayInShower(orig->id()))
	    {
	      ShowerTreePtr newtree=new_ptr(ShowerTree(orig,_showerVariables,decay,ch));
	      trees.insert(make_pair(orig,newtree));
	      Energy width=orig->dataPtr()->generateWidth(orig->mass());
	      decay.insert(make_pair(width,newtree));
	    }
  	  // now create the shower progenitors
  	  PPtr ncopy=new_ptr(Particle(*orig));
  	  copy[0]->addChild(ncopy);
  	  ShowerParticlePtr nshow=new_ptr(ShowerParticle(*ncopy,2));
  	  nshow->setFinalState(true);
  	  ShowerProgenitorPtr prog=new_ptr(ShowerProgenitor(children[ix],
  							    ncopy,nshow));
  	  _outgoingLines.insert(make_pair(prog,nshow));
	}
       // set up the map of daughter trees
       map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator mit;
       for(mit=_outgoingLines.begin();mit!=_outgoingLines.end();++mit)
	 {
	   map<PPtr,ShowerTreePtr>::const_iterator tit=trees.find(mit->first->original());
	   if(tit!=trees.end())
	     {
	       _treelinks.insert(make_pair(tit->second,mit->first));
	       tit->second->_parent=this;
	     }
	 }
    }
  // all ready decayed
   else
     {
       // need to boost the system to conserve momentum
       // find parent tree and particle
       ShowerTreePtr ptree=ShowerTreePtr(this);
       ShowerProgenitorPtr pthis=_parent->_treelinks[ptree];
       ShowerParticlePtr newparent=_parent->_outgoingLines[pthis];
       // workout the lorentz boost
       LorentzRotation boost(_incomingLines.begin()->first->original()->momentum().
			     findBoostToCM());
       boost.boost(-newparent->momentum().findBoostToCM());
       // now boost all the particles
       map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cit;
       for(cit=_incomingLines.begin();cit!=_incomingLines.end();++cit)
	 {
	   cit->first->progenitor()->deepTransform(boost);
	   cit->first->original()->deepTransform(boost);
	   cit->first->copy()->deepTransform(boost);
	 }
       for(cit=_outgoingLines.begin();cit!=_outgoingLines.end();++cit)
	 {
	   cit->first->progenitor()->deepTransform(boost);
	   cit->first->original()->deepTransform(boost);
	   cit->first->copy()->deepTransform(boost);
	 }
     }
}

void ShowerTree::insertDecay(StepPtr pstep,bool ISR, bool FSR)
{
  assert(_incomingLines.size()==1);
  _colour.clear();
  // find final particle from previous tree
  ShowerParticlePtr final=
    _parent->getFinalStateShowerProduct(_parent->_treelinks[this]);
  // construct the map of colour lines
  PPtr copy=_incomingLines.begin()->first->copy();
  if(copy->colourLine()) 
    _colour.insert(make_pair(copy->colourLine(),final->colourLine()));
  if(copy->antiColourLine())
    _colour.insert(make_pair(copy->antiColourLine(),final->antiColourLine()));
  // initial-state radiation
  if(ISR&&!_incomingLines.begin()->first->progenitor()->children().empty())
    {
      ShowerParticlePtr init=_incomingLines.begin()->first->progenitor();
      // initial particle is copy of final from previous shower
      tParticleVector parents=init->parents();
      for(unsigned int ix=0;ix<parents.size();++ix)
	{parents[ix]->abandonChild(init);}
      updateColour(init);
      final->addChild(init);
      pstep->addDecayProduct(init);
      // insert shower products
      addFinalStateShower(init,pstep);
      // sort out colour
      final=_incomingLines.begin()->second;
      _colour.clear();
      if(copy->colourLine()) 
	_colour.insert(make_pair(copy->colourLine(),final->colourLine()));
      if(copy->antiColourLine())
	_colour.insert(make_pair(copy->antiColourLine(),final->antiColourLine()));
    }
  // get the decaying particles
  // make the copy
  updateColour(copy);
  // copy of the one from the hard process
  tParticleVector dpar=copy->parents();
  for(unsigned int ix=0;ix<dpar.size();++ix)
    {dpar[ix]->abandonChild(copy);}
  final->addChild(copy);
  pstep->addDecayProduct(copy);
  ParticleVector copyc=copy->children();
  for(unsigned int ix=0;ix<copyc.size();++ix)
      copy->abandonChild(copyc[ix]);
  // final-state radiation
  if(FSR)
    {
      map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cit;
      for(cit=outgoingLines().begin();cit!=outgoingLines().end();++cit)
 	{
 	  ShowerParticlePtr init=cit->first->progenitor();
 	  if(!init->getThePEGBase()) 
 	    throw Exception() << "Final-state particle must have a ThePEGBase"
 			      << " in ShowerTree::fillEventRecord()" 
 			      << Exception::runerror;
// if not from matrix element correction
	  if(cit->first->perturbative())
	    {
	      // add the child
	      updateColour(cit->first->copy());
	      PPtr orig=cit->first->original();
	      copy->addChild(orig);
	      pstep->addDecayProduct(orig);
	      orig->addChild(cit->first->copy());
	      pstep->addDecayProduct(cit->first->copy());
	      // register the shower particle as a 
	      // copy of the one from the hard process
	      tParticleVector parents=init->parents();
	      for(unsigned int ix=0;ix<parents.size();++ix)
		{parents[ix]->abandonChild(init);}
	      (*cit).first->copy()->addChild(init);
	      pstep->addDecayProduct(init);
	      updateColour(init);
	    }
// from a matrix element correction
	  else
	  {
	      if(copy->children().end()==
		 find(copy->children().begin(),copy->children().end(),
		      cit->first->original()))
	      {
		  updateColour(cit->first->original());
		  copy->addChild(cit->first->original());
		  pstep->addDecayProduct(cit->first->original());
	      }
	      updateColour(cit->first->copy());
	      cit->first->original()->addChild(cit->first->copy());
	      pstep->addDecayProduct(cit->first->copy());
	      // register the shower particle as a 
	      // copy of the one from the hard process
	      tParticleVector parents=init->parents();
	      for(unsigned int ix=0;ix<parents.size();++ix)
	      {parents[ix]->abandonChild(init);}
	      (*cit).first->copy()->addChild(init);
	      pstep->addDecayProduct(init);
	      updateColour(init);
	  }
	  // insert shower products
	  addFinalStateShower(init,pstep);
	}
    }
  _colour.clear();
}

}



