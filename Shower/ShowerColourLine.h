// -*- C++ -*-
#ifndef HERWIG_ShowerColourLine_H
#define HERWIG_ShowerColourLine_H
//
// This is the declaration of the <!id>ShowerColourLine<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// The <!id>ShowerColourLine<!!id> class represents colour lines <BR> 
// connecting <!class>ShowerParticle<!!class>s. <BR>
// A <!id>ShowerColourLine<!!id> keeps track of the shower particles <BR> 
// connected to it. To connect a shower particle to a colour line, <BR>
// the <!id>addColoured<!!id> and <!id>addAntiColoured<!!id> <BR>
// methods should be used - these will automatically set up the <BR>
// <!class>ShowerParticle<!!class> correctly. <BR>
//
// If a colour line stems from a colour source or ends in a colour <BR>
// sink, it is possible to obtain the neighbouring colour lines. <BR>
// This is also the way sinks and sources are implemented.
//
// Notice that it has been necessary to introduce this new class <BR>
// <!id>ShowerColourLine<!!id> because the similar Pythia7 class <BR>
// <!id>ColourLine<!!id> can be used only with Pythia7 Particle. <BR>
// However, <!id>ShowerColourLine<!!id> class is almost identical to <BR>
// Pythia7 <!id>ColourLine<!!id> class.

// CLASSDOC SUBSECTION See also:
//
// <a href="http:ColourLine.html">ColourLine.h</a>,
// <a href="http:ShowerParticle.html">ShowerParticle.h</a>,
// 

#include "ShowerConfig.h"
#include "Pythia7/Pointer/Ptr.h"
#include "Pythia7/Pointer/ReferenceCounted.h"
#include "Pythia7/Pointer/PtrTraits.h"
#include "Pythia7/Pointer/RCPtr.h"
// #include "Pythia7/EventRecord/ColourLine.h"
#include "ShowerParticle.h"


namespace Herwig {

  using namespace Pythia7;

  class ShowerColourLine: public ReferenceCounted {

  public:

    inline ShowerColourLine();
    inline ShowerColourLine(const ShowerColourLine &);
    virtual ~ShowerColourLine();
    // Standard ctors and dtor.

    void addAntiColoured(tShoParPtr);
    // Add an anti-coloured shower particle to this line.
    
    void addColoured(tShoParPtr, bool anti = false);
    // Add a coloured shower particle to this line.
    
    void removeAntiColoured(tShoParPtr);
    // Remove an anti-coloured shower particle to this line.
    
    void removeColoured(tShoParPtr, bool anti = false);
    // Remove a coloured shower particle to this line.

    inline const tCollecShoParPtr & coloured() const;
    inline const tCollecShoParPtr & antiColoured() const;
    // Return the vectors of shower particles connected to this line 
    // with their coloures or anti colours.

    inline tShoParPtr startParticle() const;
    // Return the first particle on this colour line. Returns null if
    // this line stems from a colour source. If the particle is
    // outgoing, its anti colour is connected, otherwise its colour is
    // connected.

    inline tShoParPtr endParticle() const;
    // Return the last particle on this colour line. Returns null if
    // this line ends in a colour sink. If the particle is outgoing, its
    // colour is connected, otherwise its anti colour is connected.

    inline tShoColinePair sinkNeighbours() const;
    // If this colour line ends in a colour sink, these two colour lines
    // ends in the same.

    inline tShoColinePair sourceNeighbours() const;
    // If this colour line stems from a colour source, these two colour
    // lines stems from the same.

    inline void setSinkNeighbours(tShoColinePtr, tShoColinePtr);
    inline void setSourceNeighbours(tShoColinePtr, tShoColinePtr);
    // Add two colour lines as neighbour to this lines sink or
    // source. Also the neighbors are set up correspondingly.
    
  private:

    tCollecShoParPtr theColoured;
    tCollecShoParPtr theAntiColoured;
    // The particles connecting to this colour line, either following
    // the incoming colour flow or the outgoing one.

    tShoColinePair theSourceNeighbours;
    // If this colour line stems from a colour source, these two colour
    // lines stems from the same.
    
    tShoColinePair theSinkNeighbours;
    // If this colour line ends in a colour sink, these two colour lines
    // ends in the same.
    
  private:
    
    ShowerColourLine & operator=(const ShowerColourLine &);
    //  Private and non-existent assignment operator.
    
  };

}

#include "ShowerColourLine.icc"

#endif /* HERWIG_ShowerColourLine_H */
