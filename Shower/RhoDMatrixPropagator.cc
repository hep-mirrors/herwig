// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RhoDMatrixPropagator class.
//

#include "RhoDMatrixPropagator.h"
#include "Pythia7/Interface/ClassDocumentation.h"
#include "Pythia7/Persistency/PersistentOStream.h"
#include "Pythia7/Persistency/PersistentIStream.h"
#include "ShowerParticle.h"
#include "Pythia7/Handlers/PartialCollisionHandler.h"
#include "Pythia7/PDT/Decayer.h"
#include "SplitFun.h"
#include "Pythia7/Repository/EventGenerator.h"
#include "Herwig++/Utilities/HwDebug.h"
#include "Pythia7/Interface/Switch.h"
#include "QtildaShowerKinematics1to2.h"
#include "SplitFun1to2.h"
// #include "Pythia7/MatrixElement/Amplitude.h"

using namespace Herwig;


RhoDMatrixPropagator::~RhoDMatrixPropagator() {}


void RhoDMatrixPropagator::persistentOutput(PersistentOStream & os) const {
  os << _onoffSwitchMode; 
}


void RhoDMatrixPropagator::persistentInput(PersistentIStream & is, int) {
  is >> _onoffSwitchMode; 
}


ClassDescription<RhoDMatrixPropagator> RhoDMatrixPropagator::initRhoDMatrixPropagator;
// Definition of the static class description member.


void RhoDMatrixPropagator::Init() {

  static ClassDocumentation<RhoDMatrixPropagator> documentation
    ("This class is responsible for the propagation of spin density and decay matrices");

  static Switch<RhoDMatrixPropagator, int> interfaceOnOffSwitchMode
    ("OnOffRhoDPropagationMode",
     "Choice of the on-off rhoD matrix propagation switch mode",
     &RhoDMatrixPropagator::_onoffSwitchMode, 1, false, false);
  static SwitchOption interfaceOnOffSwitchMode0
    (interfaceOnOffSwitchMode,"RhoDPropagation-OFF","rhoD propagation is OFF", 0);
  static SwitchOption interfaceOnOffSwitchMode1
    (interfaceOnOffSwitchMode,"RhoDPropagation-ON","rhoD propagation is ON", 1);

}


// ------------------------------ PUBLIC ---------------------------------------

double RhoDMatrixPropagator::matrixElement( const tMEPtr hardSubME,  
					    const CollecShoParPtr & collecShoPar,
					    const tShoParPtr aParticle ) {
  double result = 1.0;
  if ( isRhoDPropagationON() ) {
    result = 0.0;
    if ( setVertexPointer(hardSubME,aParticle) ) {
      CollecShoParPtr particles;
      findVertexParticles(collecShoPar, aParticle, particles);
      result = matrixElement(particles);
    }
  }
  return result;
}


bool RhoDMatrixPropagator::computeRhoD( const tMEPtr hardSubME, 
				        const CollecShoParPtr & collecShoPar,
				        const tShoParPtr theParticle ) {
  bool isOK = true;
  if ( isRhoDPropagationON() ) {
    isOK = setVertexPointer(hardSubME,theParticle);
    if ( isOK ) {
      CollecShoParPtr particles;
      findVertexParticles(collecShoPar, theParticle, particles);
      if ( computeRhoD(theParticle,particles) ) {
	normalize( theParticle->rhoD() );
      } else {
	isOK = false;
      }
    }
  }
  return isOK;
}


// ------------------------------ PRIVATE ---------------------------------------

bool RhoDMatrixPropagator::setVertexPointer( const tMEPtr hardSubME, 
					     const tShoParPtr particlePtr ) {
  // First set to null the three possible "vertex" pointers
  _hardSubME = hardSubME;
  _decayer = tDecayerPtr();
  _splitFun = tSplitFunPtr(); 
  // Now set at most one of the three vertices. If  hardSubME  is not null
  // then the hard subprocess is considered; otherwise it should be either 
  // the decay vertex or the splitting one.
  // Notice that in the latter two cases the particle in input (particlePtr).
  // must be the decaying or the splitting one, not one of the decay or
  // splitting products. This is also true for intial state radiation
  // (backward evolution), because we define "parent" in ShowerParticle
  // the particle which is closer to the hard subprocess, although the
  // physically it is a child of the incoming particle after an initial
  // state radiation. 
  bool foundVertex = false;
  if ( hardSubME ) {
    foundVertex = true;
  } else {
    if ( particlePtr ) {
      if ( particlePtr->splitFun() ) {
	_splitFun = particlePtr->splitFun();
	foundVertex = true;
      } else if ( particlePtr->decayer() ) {
	_decayer = particlePtr->decayer();
	foundVertex = true;
      }
    }
  }
  return foundVertex;
} 


void RhoDMatrixPropagator::findVertexParticles( const CollecShoParPtr & collecShoPar,
					        const tShoParPtr aParticle,
					        CollecShoParPtr & particles ) {
  // Given a pointer  aParticle  to a shower particle, the output of this
  // method is the vector  particles  filled with all the pointers 
  // (including  aParticle  itself) to the shower particles entering
  // the same vertex. In the case of decay or splitting vertex, this
  // is easy done by considering the children of  aParticle  , being
  // for assumption the parent of that vertex. In the case of hard
  // subprocess vertex, we need the full collection of shower particles
  //  collecShoPar  in order to find the ones that belong to the hard
  // subprocess: in fact, the subprocess methods
  //  _subprocess->incoming()   and   _subprocess->outgoing()
  // cannot be used because they return Pythia7 particle objects, not
  // ShowerParticle ones as we want.
  if ( ! aParticle ) return; 
  if ( _hardSubME ) {
    for ( CollecShoParPtr::const_iterator cit = collecShoPar.begin();
	  cit != collecShoPar.end(); ++cit ) {
      if ( (*cit)->isFromHardSubprocess() ) {
	particles.push_back( *cit );
      }
    }
  } else if ( _decayer || _splitFun ) {
    particles.push_back( aParticle );
    for ( CollecShoParPtr::const_iterator cit = aParticle->children().begin();
	  cit != aParticle->children().end(); ++cit ) {
      particles.push_back( *cit );
    }
  }
}


Complex RhoDMatrixPropagator::trace( const ComplexMatrix & matrix ) const {
  Complex sumDiag(0.0,0.0);
  int numRows = matrix.size();
  for (int i=0; i < numRows; ++i) {
    sumDiag += matrix[i][i];
  }
  return sumDiag;
}


bool RhoDMatrixPropagator::normalize( ComplexMatrix & matrix ) {
  // If the trace of the matrix is zero then return false;
  // If the trace is 1.0 don't do anything and return true;
  // otherwise, divide all elements of the matrix by the trace value
  // and return true;
  Complex traceVal = trace( matrix );
  if ( norm(traceVal) < 1.0e-20 ) return false;
  if ( norm( traceVal - Complex(1.0,0.0) ) > 1.0e-20 ) {
    int numRows = matrix.size();
    for ( int i=0; i < numRows; i++ ) {
      for ( int j=0; j < numRows; j++ ) {
	matrix[i][j] /= traceVal;
      }
    }
  }
  return true;
}


double RhoDMatrixPropagator::matrixElement( const CollecShoParPtr & particles ) {
  // This method returns the real value obtained by the product of all the
  // rho-D matrices entering the vertex and the amplitude and its
  // conjugate, summing over all helicity indeces.
  // The vector of ShowerParticle pointers, particles, contains the
  // pointers to all the particles entering the considered vertex.
  // It is important to notice that we are assuming that all the rhoD
  // matrices involved are already normalised (to trace equal 1). That
  // should be automatic, because the initialized value assign to such
  // matrices (in ShowerParticle constructor) is the normalise 2x2 
  // identity, and then each time a rhoD matrix is computed, calling
  // the method  RhoDMatrixPropagator::computeRhoD(...)  
  // it automatically normalises it. 

  Complex sum(0.0,0.0); // at the end we check that the result should be real.
  int numParticles = particles.size();
  // Keep track of the helicity range for each particle, which is equal
  // to the size of the rhoD matrix associated to it (if the particle
  // is massive and has spin S>0, the helicity range is between 0 and
  // 2*S, and the rhoD matrix is a  (2*S+1)x(2*S+1) complex matrix).
  tcPDVector dataParticles( numParticles );
  vector<int> sizes( numParticles );
  vector<int> indeces( numParticles );
  vector<int> indecesPrime( numParticles );
  for ( int i=0; i < numParticles; i++ ) {
    dataParticles[i] = particles[i]->dataPtr();
    sizes[i] = particles[i]->rhoD().size();
    indeces[i] = indecesPrime[i] = 0;
  }
  bool newIndexConfiguration = true;
  while ( newIndexConfiguration ) {
    // Calculate the factor caming from the product of all the rho-D matrices,
    // and initialize the vectors of pairs of momenta and indeces, needed for
    // the amplitude evaluation.
    Complex matrixFactors(1.0);
    vector<Lorentz5Momentum> momenta;     
    for (int i=0; i < numParticles; i++) {
      matrixFactors *= particles[i]->rhoD()[ indeces[i] ][ indecesPrime[i] ];
      momenta[i] = particles[i]->momentum();
    }
    Complex amplitudeValue = 0.0, amplitudePrimeValue = 0.0;
    evaluateAmplitude( particles[0], dataParticles, momenta, indeces, indecesPrime,
		       amplitudeValue, amplitudePrimeValue );
    // Calculate the full contribution, for this index configuration,
    // by taking into account also the amplitude and its conjugate.
    sum += matrixFactors * amplitudeValue * conj( amplitudePrimeValue );
    // Determine the next index configuration. 
    newIndexConfiguration = nextIndexConfiguration( sizes, indeces, indecesPrime );
  };
  // Check if the result is real, as expected (rho and D matrices are hermitian)
  if ( fabs( sum.imag() ) > 1.0e-10 ) {
    generator()->logWarning( Exception("RhoDMatrixPropagator::matrixElement "
                                       "***matrixElement NOT Real*** ", 
                                       Exception::warning) );
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
      generator()->log() << "         ===>" << " matrixElement = " << sum 
                         << endl << endl;
    }
  }
  return sum.real();
}


bool RhoDMatrixPropagator::computeRhoD( const tShoParPtr theParticle, 
				        const CollecShoParPtr & particles ) {
  //  particles  is a vector of pointers to all the particles (objects
  // of class ShowerParticle) entering the same vertex. This method 
  // computes the rhoD matrix of one of them, pointed by  theParticle  ,
  // in terms of the rhoD matrices of the other particles and the 
  // amplitude (and its conjugate) of the vertex as a function of the 
  // helicities of the particles.
  // The method returns true if it succeed, false otherwise.

  if ( ! theParticle  ) return false;
  // If the particle is a scalar, or a massless spin 1/2 particles,
  // then there is nothing to compute. 
  if ( theParticle->rhoD().size() == 1 ) return true;
  int thePosition = -1;
  int numParticles = particles.size();
  // Keep track of the helicity range for each particle, which is equal
  // to the size of the rhoD matrix associated to it (if the particle
  // is massive and has spin S>0, the helicity range is between 0 and
  // 2*S, and the rhoD matrix is a  (2*S+1)x(2*S+1) complex matrix).
  tcPDVector dataParticles( numParticles );
  vector<int> sizes( numParticles );
  vector<int> indeces( numParticles );
  vector<int> indecesPrime( numParticles );
  for ( int i=0; i < numParticles; i++ ) {
    if ( particles[i] == theParticle ) {
      thePosition = i;
    }
    dataParticles[i] = particles[i]->dataPtr();
    sizes[i] = particles[i]->rhoD().size();
    indeces[i] = indecesPrime[i] = 0;
  }
  //  thePosition  is the index, in the vector particles, of 
  // the particle whose rhoD we have to determine.
  if ( thePosition < 0 ) return false;
  // Reset the previous rhoD matrix.
  for ( int i=0; i < sizes[thePosition]; i++ ) {
    for ( int j=0; j < sizes[thePosition]; j++ ) {
      theParticle->rhoD()[i][j] = 0.0;
    }
  }
  bool newIndexConfiguration = true;
  while ( newIndexConfiguration ) {
    // Calculate the factor caming from the product of all the rhoD matrices
    // excluding only the one we want to calculate its rhoD matrix;
    // initialize also the vectors of pairs of momenta and indeces, needed for
    // the amplitude evaluation.
    Complex matrixFactors(1.0);
    vector<Lorentz5Momentum> momenta;     
    for (int i=0; i < numParticles; i++) {
      if ( i != thePosition ) {
	matrixFactors *= particles[i]->rhoD()[ indeces[i] ][ indecesPrime[i] ];
      }
      momenta[i] = particles[i]->momentum();
    }
    Complex amplitudeValue = 0.0, amplitudePrimeValue = 0.0;
    evaluateAmplitude( particles[0], dataParticles, momenta, indeces, indecesPrime,
		       amplitudeValue, amplitudePrimeValue );
    // Calculate the full contribution, for this index configuration,
    // by taking into account also the amplitude and its conjugate.
    // Notice that for different values of  indeces[thePosition] and/or
    // indecesPrime[thePosition]  we actually compute different elements
    // of the rhoD matrix.  
    theParticle->rhoD()[ indeces[thePosition] ][ indecesPrime[thePosition] ] += 
      matrixFactors * amplitudeValue * conj( amplitudePrimeValue );
    // Determine the next index configuration. 
    newIndexConfiguration = nextIndexConfiguration(sizes, indeces, indecesPrime);
  }
  return true;
}


void RhoDMatrixPropagator::
evaluateAmplitude( const tShoParPtr particle,                                           // in
		   const tcPDVector & dataParticles,                                    // in
		   const vector<Lorentz5Momentum> & momenta,                            // in
		   const vector<int> & helicities, const vector<int> & helicitiesPrime, // in
		   Complex & amplitudeValue, Complex & amplitudePrimeValue ) {          // out
  amplitudeValue = amplitudePrimeValue = 1.0;
  if ( _hardSubME  &&  _hardSubME->amplitude() ) {
    amplitudeValue = _hardSubME->amplitude()->value( dataParticles, momenta, helicities );
    amplitudePrimeValue =  _hardSubME->amplitude()->value( dataParticles, momenta, helicitiesPrime );
  } else if ( _decayer  &&  _decayer->amplitude() ) {
    amplitudeValue = _decayer->amplitude()->value( dataParticles, momenta, helicities );
    amplitudePrimeValue = _decayer->amplitude()->value( dataParticles, momenta, helicitiesPrime );     
  } else if ( _splitFun  &&  particle  &&  particle->showerKinematics() ) {
    //***LOOKHERE*** Only in this implementation part of the all
    // RhoDMatrixPropagator class we have to esplicitly distinguish
    // between the various multiplicities n of final states in the
    // splitting process 1->n. The reason is that for the common
    // 1->2 case, treated below, the dependence of the splitting
    // function in on (z,phi) and helicities, whereas for 1->3 case
    // the (z,phi) dependence is replaced by another, with more
    // variables, one (and also on helicities, but the latter are
    // in the form of a generic vector, therefore they do not need
    // different treatment as a function of n).
    // As you can see below, the right, safe way to distinguish
    // between the various n is by using downcasting.
    // If you want to treat the case of 1->3 splitting processes
    // you should add some code very similar to the one below,
    // replacing "2" with "3", and (z,phi) with other variables.
    typedef Ptr<SplitFun1to2>::transient_pointer tSplitFun1to2Ptr;
    tSplitFun1to2Ptr splitFun1to2 = dynamic_ptr_cast< tSplitFun1to2Ptr >( _splitFun );
    if ( splitFun1to2 ) {
      typedef Ptr<QtildaShowerKinematics1to2>::transient_pointer tQtildaShoKin1to2Ptr;
      tQtildaShoKin1to2Ptr shoKin = 
	dynamic_ptr_cast< tQtildaShoKin1to2Ptr >( particle->showerKinematics() );
      if ( shoKin ) {
	double z = shoKin->z();
	Energy2 qtilde2 = sqr( shoKin->qtilda() );
	double phi = shoKin->phi();
	amplitudeValue = splitFun1to2->
	  fullFunWithHelicities( z, qtilde2, phi, helicities[0], helicities[1], helicities[2] );
	amplitudePrimeValue = splitFun1to2->
	  fullFunWithHelicities( z, qtilde2, phi, helicitiesPrime[0], helicitiesPrime[1], helicitiesPrime[2] );
      }
    }
  }
}


bool RhoDMatrixPropagator::nextIndexConfiguration( const vector<int> & sizes,
						   vector<int> & indeces, 
						   vector<int> & indecesPrime ) {
  // If we knew in advance the number of particles which enter a
  // "vertex" (either hard subprocess, or decay, or splitting) then 
  // we would get the summation over all helicities configurations
  // by simply nesting a certain number (precisely, 2 * number of
  // particles) of for loops (as done in Fortran Herwig).
  // However, to avoid code duplication, we prefer to keep the
  // code general enough to treat any number of particles, 
  // therefore we have to generate ourselves all the possible
  // helicity configurations, without neither missing any of them
  // nor overcounting them. The way we generate such configurations
  // is better to be explained by showing explicitly the order in 
  // which they are produces:
  //
  //  indeces[0]  indecesPrime[0]  ... indeces[n]  indecesPrime[n]
  //     0               0         ...      0             0
  //     0               0         ...      0             1
  //    ...             ...        ...     ...           ...
  //     0               0         ...      0           sizes[n]-1
  //     0               0         ...      1             0
  //     0               0         ...      1             1
  //    ...             ...        ...     ...           ...
  //     0               0         ...    sizes[n]-1    sizes[n]-1
  //    ...             ...        ...     ...           ...
  //     0               1         ...      0             0
  //    ...             ...        ...     ...           ...
  //     0             sizes[0]-1  ...    sizes[n]-1    sizes[n]-1
  //     1               0         ...      0             0
  //    ...             ...        ...     ...           ...
  //   sizes[0]-1      sizes[0]-1  ...    sizes[n]-1    sizes[n]-1                      
  // 
  // The vector  sizes  contains the max values for all the indeces
  // and it is not modify by this method; the other two vectors,
  //  indeces  and  indecesPrime  contain the current (old, that is
  // already used) index configuration, and at the end of the method,
  // if there is still new index configurations to be considered
  // (in this case the method returns true; otherwise false),
  // the new index configuration is overwritten in the same
  //  indeces  and  indecesPrime  vectors.
 
  int currentPos = indeces.size() - 1;
  bool incremented = false;
  do {
    if ( indecesPrime[currentPos]  <  sizes[currentPos] - 1 ) {
      indecesPrime[currentPos]++;
      incremented = true;
    } else {
      indecesPrime[currentPos] = 0;
      if ( indeces[currentPos]  <  sizes[currentPos] - 1 ) {
	indeces[currentPos]++;
	incremented = true;
      } else {
	indeces[currentPos] = 0;
	currentPos--;
      }
    }
  } while ( currentPos >= 0   &&  ! incremented );
  return incremented;
}

