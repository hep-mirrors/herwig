// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PartonSplitter class.
//

#include "PartonSplitter.h"
#include <ThePEG/Interface/ClassDocumentation.h>
#include <ThePEG/Interface/Reference.h>
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/PDT/EnumParticles.h>
#include <ThePEG/EventRecord/Step.h>
#include <ThePEG/Repository/EventGenerator.h>
#include "Herwig++/Utilities/Kinematics.h"
#include "Herwig++/Utilities/HwDebug.h"

using namespace Herwig;

void PartonSplitter::persistentOutput(PersistentOStream &) const {
}


void PartonSplitter::persistentInput(PersistentIStream &, int) {
}


ClassDescription<PartonSplitter> PartonSplitter::initPartonSplitter;
// Definition of the static class description member.


void PartonSplitter::Init() {

  static ClassDocumentation<PartonSplitter> documentation
    ("This class is reponsible of the nonperturbative splitting of partons"
     "\n (mainly time-like gluons, but also space-like sea partons");
 
}


tPVector PartonSplitter::split(const tPVector & tagged, tStepPtr pstep) {

  tParticleVector newPartons;  // only for debugging

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
    //debuggingInfo(tagged,newPartons);
    generator()->log() << "######## Particles before PartonSplitter #######\n";
    generator()->log() << *pstep;
    generator()->log() << "################################################\n";
  }
  tPVector newtag;
  // Loop over all of the particles in the event.
  for(tPVector::const_iterator pit = tagged.begin(); pit!=tagged.end(); ++pit) 
    {
      // split the gluons
      if ( (**pit).data().id() == ParticleID::g ) 
	{
	  // time like gluon
	  if ( (**pit).momentum().m2() > 0.0 )
	    {
	      PPtr ptrQ = PPtr();
	      PPtr ptrQbar = PPtr();
	      splitTimeLikeGluon(*pit,ptrQ,ptrQbar);
	      Energy Q0 = getParticleData(ParticleID::g)->constituentMass();
	      ptrQ->scale(Q0*Q0);
	      ptrQbar->scale(Q0*Q0);
	      pstep->addDecayProduct(*pit,ptrQ);
	      pstep->addDecayProduct(*pit,ptrQbar);
	      newtag.push_back(ptrQ);
	      newtag.push_back(ptrQbar);
	      pstep->fixColourFlow();   
	      // Debugging
	      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
		newPartons.push_back(ptrQ);
		newPartons.push_back(ptrQbar);
	      }
	    }
	  // space-like gluons
	  else if( (**pit).momentum().m2() < 0.0 ) 
	    {
	      // ... write the code
	      // splitSpaceLikeGluon(*pit,ptrQ,ptrQbar);      
	      // ... write the code
	      cerr << "\nSpacelike gluon in PartonSplitter::split()\n";
	      throw Exception() << "Spacelike gluon in PartonSplitter::split()"
			       << Exception::eventerror;
	    }
	  // q^2=0 gluon 
	  else 
	    {
	      generator()->logWarning( Exception("PartonSplitter::split "
						 "***Gluon on the mass shell (m=0)***", 
						 Exception::warning) );
	      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
		generator()->log() << "         ===>" << (**pit).momentum() << '\n' << '\n';
	      }
	    }
	}
      else 
	{
	  newtag.push_back(*pit);
	  // write a method in Herwig++/Utilities/CheckId (or use/add 
	  //  ThePEG/PDT/StandardMatchers.h) for checking if it is a 
	  // sea quark or antiquark
	  //if ( ( (**pit).momentum().m2() < 0.0 ) &&               // space like and
	  //     ( SeaQuarkMatcher::Check( (**pit).data() ) ||         // (sea quark or
	  //       SeaAntiQuarkMatcher::Check( (**pit).data() ) ) ) {  //  sea anti-quark)
	  //   ... write the code
	  //   splitSpaceLikeSeaQuark(*pit,ptrGluon,ptrSeaQ1); 
	  //   splitSpaceLikeGluon(ptrGluon,ptrQ,ptrQbar);      
	  //   ... write the code
	  //}
	}
    } // end for loop over tagged 
  
  // Debugging
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
    //debuggingInfo(tagged,newPartons);
    generator()->log() << "######## Particles after PartonSplitter ########\n";
    generator()->log() << *pstep;
    generator()->log() << "################################################\n";
  }
  return newtag;
}


void PartonSplitter::splitTimeLikeGluon(tcPPtr ptrGluon, PPtr & ptrQ, PPtr & ptrQbar){

  // Choose the flavour of the quark (u or d with equal 50% probability)
  long newId = 0;
  if ( UseRandom::rndbool() ){
    newId = ParticleID::u;
  } else {
    newId = ParticleID::d;
  } 

  // Solve the kinematics of the two body decay  G --> Q + Qbar
  Lorentz5Momentum momentumQ = Lorentz5Momentum();
  Lorentz5Momentum momentumQbar = Lorentz5Momentum();
  double cosThetaStar = UseRandom::rnd( -1.0 , 1.0 );
  double phiStar = UseRandom::rnd( -pi , pi );
  Energy constituentQmass = getParticleData(newId)->constituentMass();

  if (ptrGluon->momentum().m() < 2.0*constituentQmass) {
    throw Exception() << "Impossible Kinematics in PartonSplitter::splitTimeLikeGluon()" 
		      << Exception::eventerror;
  }
  Kinematics::twoBodyDecay(ptrGluon->momentum(), constituentQmass, 
			   constituentQmass, cosThetaStar, phiStar, momentumQ, 
			   momentumQbar ); 

  // Create quark and anti-quark particles of the chosen flavour 
  // and set they 5-momentum (the mass is the constituent one).
  ptrQ    = getParticle(newId);
  ptrQbar = getParticle(-newId);
  ptrQ->set5Momentum( momentumQ );
  ptrQbar->set5Momentum( momentumQbar );

  // Sanity check (normally skipped) to see if the energy-momentum is conserved.
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Hadronization ) {    
    Lorentz5Momentum diff = ptrGluon->momentum() - 
                            ( ptrQ->momentum() + ptrQbar->momentum() );
    Energy ediff = fabs( diff.m() );
    if ( ediff > 1e-3*GeV ) {
      generator()->logWarning( Exception("PartonSplitter::splitTimeLikeGluon " 
					 "***Violation of energy-momentum conservation***", 
					 Exception::warning) );    
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
	generator()->log() << "         ===> g -> " << newId << " " << -newId << '\n'
			   << "   diff " << diff << '\n'
			   << "   " << ptrGluon->momentum() << " ---> "
			   << ( ptrQ->momentum() + ptrQbar->momentum() ) << '\n' 
			   << " = " << ptrQ->momentum() << " + " << ptrQbar->momentum() 
			   << '\n' << '\n';
      }      
    }
  }
}


/*void PartonSplitter::splitSpaceLikeGluon(tcPPtr ptrGluon, PPtr & ptrQ, PPtr & ptrQbar){
  // write the code
}


void PartonSplitter::splitSpaceLikeSeaQuark(tcPPtr ptrSeaQ0, PPtr & ptrGluon, PPtr & ptrSeaQ1){
  // write the code
}

*/
void PartonSplitter::debuggingInfo(const tPVector & tagged, const set<tPPtr> & newPartons) {

  // Print information about all coloured particles present before the
  // parton splitting, and then print the ones created by the parton
  // splitting. Before doing that, associate to each of them an integer 
  // number (starting with 1), in order to translate in a readable way 
  // the cross references between particles. 
  int count = 0;
  map<tPPtr,int> orderingMap;
  for ( tPVector::const_iterator it = tagged.begin(); it != tagged.end(); ++it ) {
    if ( (*it)->coloured() ) {
      orderingMap.insert( orderingMap.end(), pair<tPPtr,int>(*it,++count) ); 
    }
  }
  int save_count = count;
  for (set<tPPtr>::const_iterator it = newPartons.begin(); it != newPartons.end(); ++it) {
    orderingMap.insert( orderingMap.end(), pair<tPPtr,int>(*it,++count) ); 
    generator()->log() << "Inserting " << *(*it) << "\ninto orderingMap\n";
  }

  // Do the same for colour and anti-colour lines.
  int countColines = 0;
  map<tColinePtr,int> colineMap;
  for ( tPVector::const_iterator it = tagged.begin(); it != tagged.end(); ++it ) {
    if ( (*it)->colourLine() ) {
      colineMap.insert( colineMap.end(), pair<tColinePtr,int>((*it)->colourLine(),++countColines) ); 
    }
    if ( (*it)->antiColourLine() ) {
      colineMap.insert( colineMap.end(), pair<tColinePtr,int>((*it)->antiColourLine(),++countColines) ); 
    }
  }

  generator()->log() << "PartonSplitter::debuggingInfo ===> START DEBUGGING <=== " 
		     << "   EventNumber=" << generator()->currentEventNumber() << '\n'
		     << "  INITIAL PARTONS : num = " << save_count << '\n'
		     << "  NEW     PARTONS : num = " << count - save_count << '\n';

  /*  for ( map<tPPtr,int>::const_iterator it = orderingMap.begin(); 
	it != orderingMap.end(); ++it ) {
    tPPtr pptr = it->first;  
    int i = it->second;  
    generator()->log() << "  --- Parton --- " << i << '\n'
		       << "\t id = " << pptr->id() << "     " << pptr->PDGName() << '\n'
                       << "\t masses:  constituent=" << pptr->data().constituentMass()
		       << "  current=" << pptr->mass() 
		       << "  invariant=" << pptr->momentum().m() << '\n'
		       << "\t momentum= " << pptr->momentum() << '\n'
                       << "\t scale=" 
                       << ( pptr->scale() > 0 ? sqrt( fabs( pptr->scale() ) ) : 
			    - sqrt( fabs( pptr->scale() ) ) ) 
		       <<  '\n';

    if ( pptr->final()->children() != ParticleVector() ) {
      generator()->log() << "\t Childrens :  "; 
      for ( ParticleVector::const_iterator jt = pptr->final()->children().begin();
	    jt != pptr->final()->children().end(); ++jt ) {	
	int child = 0;
	if ( orderingMap.find( *jt ) != orderingMap.end() ) {
	  child = orderingMap.find( *jt )->second;
	} else {
	  child = -999;  // Error: it shouldn't happen!
	}
	generator()->log() << child << "   ";	
      }
      generator()->log() << '\n';
    }

    if ( i > save_count ){            
      if ( pptr->parents() != tParticleVector() ) {
	generator()->log() << "\t Parents :  "; 
	for ( tParticleVector::const_iterator jt = pptr->parents().begin();
	      jt != pptr->parents().end(); ++jt ) {	
          int parent = -999;
	  tPPtr parentPtr = *jt;
          while ( parentPtr  &&  orderingMap.find( parentPtr ) == orderingMap.end() ) {
	    parentPtr = parentPtr->previous();
	  }
	  if ( parentPtr ) parent = orderingMap.find( parentPtr )->second;
	  generator()->log() << parent << "   ";	
	}
	generator()->log() << '\n';
      }      
    }

    // Print info about colour lines, and eventual colour sink/source.
    for (int k = 0; k < 4; ++k) {
      bool condition = false;
      string message;
      int numLine = 0;
      tColinePair colinePair = tColinePair();
      switch (k) {
      case 0: {               // colourLine and eventual sink colourLine
	if ( pptr->colourLine() ) {
	  generator()->log() << "\t Colour line :  "; 
	  if ( colineMap.find( pptr->colourLine() ) != colineMap.end() ) {
	    numLine = colineMap.find( pptr->colourLine() )->second;
	  } else {
	    numLine= -999;  // Error: it shouldn't happen!
	  }
	  generator()->log() << numLine << '\n';
	  condition = ( pptr->colourLine()->sinkNeighbours() != tColinePair() );     
	  message = "\t *** Sink colour line :  "; 
	  colinePair = pptr->colourLine()->sinkNeighbours();
	}	
	break;
      }
      case 1: {               // colourLine and eventual source colourLine
	if ( pptr->colourLine() ) {
	  if ( colineMap.find( pptr->colourLine() ) != colineMap.end() ) {
	    numLine = colineMap.find( pptr->colourLine() )->second;
	  } else {
	    numLine= -999;  // Error: it shouldn't happen!
	  }
	  condition = ( pptr->colourLine()->sourceNeighbours() != tColinePair() );     
	  message = "\t *** Source colour line :  "; 
	  colinePair = pptr->colourLine()->sourceNeighbours();
	}
	break;
      }
      case 2: {               // antiColourLine and eventual sink antiColourLine
	if ( pptr->antiColourLine() ) {
	  generator()->log() << "\t AntiColour line :  "; 
	  if ( colineMap.find( pptr->antiColourLine() ) != colineMap.end() ) {
	    numLine = colineMap.find( pptr->antiColourLine() )->second;
	  } else {
	    numLine= -999;  // Error: it shouldn't happen!
	  }
	  generator()->log() << numLine << '\n';
	  condition = ( pptr->antiColourLine()->sinkNeighbours() != tColinePair() );     
	  message = "\t *** Sink antiColour line :  "; 
	  colinePair = pptr->antiColourLine()->sinkNeighbours();
	}
	break;
      }
      case 3: {               // antiColourLine and eventual source antiColourLine
	if ( pptr->antiColourLine() ) {
	  if ( colineMap.find( pptr->antiColourLine() ) != colineMap.end() ) {
	    numLine = colineMap.find( pptr->antiColourLine() )->second;
	  } else {
	    numLine= -999;  // Error: it shouldn't happen!
	  }
	  condition = ( pptr->antiColourLine()->sourceNeighbours() != tColinePair() );     
	  message = "\t *** Source antiColour line :  "; 
	  colinePair = pptr->antiColourLine()->sourceNeighbours();
	}
	break;
      }
      }
      if ( condition ) {
	generator()->log() << message << numLine; 
	for ( int kk = 0; kk < 2; ++kk ) {
	  tColinePtr colinePtr = colinePair.first;
	  if (kk == 1) colinePtr = colinePair.second;
	  if ( colinePtr != tColinePtr() ) {
	    int colineNum = 0;
	    if ( colineMap.find( colinePtr ) != colineMap.end() ) {
	      colineNum = colineMap.find( colinePtr )->second;
	    } else {
	      colineNum = -999;  // Error: it shouldn't happen!
	    }
	    generator()->log() << "  " << colineNum;
	  }
	}
	generator()->log() << '\n';
      }
    }
  } // end main for loop over orderingMap;
  */
  for(tPVector::const_iterator it = tagged.begin(); it != tagged.end(); ++it) {
    generator()->log() << *(*it);
  }
  for (set<tPPtr>::const_iterator it = newPartons.begin(); 
       it != newPartons.end(); ++it) {
    generator()->log() << *(*it);
  }
  if ( save_count == count ) {
    generator()->log() << " \t NO New Partons \n";
  }
  
  generator()->log() << "PartonSplitter::debuggingInfo ===> END DEBUGGING <===" << endl;

}


