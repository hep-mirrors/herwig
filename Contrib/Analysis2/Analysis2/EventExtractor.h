// -*- C++ -*-

// (C) 2007-2009 Simon Plaetzer -- sp@particle.uni-karlsruhe.de

#ifndef Analysis2_EventExtractor_H
#define Analysis2_EventExtractor_H
//
// This is the declaration of the EventExtractor class.
//

#include "ThePEG/Interface/Interfaced.h"

#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Vectors/Lorentz5Vector.h"

namespace Analysis2 {

using namespace ThePEG;

/**\ingroup Analysis2
 * 
 * @author Simon Plaetzer
 *
 * @see \ref EventExtractorInterfaces "The interfaces"
 * defined for EventExtractor.
 */
class EventExtractor: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline EventExtractor()
    : _lastEvent(), _ok(false), _doneEvent(false), _lastMomenta(), _lastWeight(0.) {}

  /**
   * The destructor.
   */
  virtual ~EventExtractor();
  //@}

public:

  virtual void prepare ();

  inline void use (tcEventPtr ev) { _lastEvent = ev; prepare(); }

  inline tcEventPtr lastEvent() const { return _lastEvent; }

public:

  virtual inline EventExtractor& operator >> (pair<vector<Lorentz5Momentum>,double> & ev) {
    if (_doneEvent) {
      _ok = false;
      return *this;
    }
    ev = make_pair(_lastMomenta,_lastWeight);
    _doneEvent = true;
    return *this;
  }

public:

  inline operator bool () const { return _ok; }

  inline bool operator ! () const { return !_ok; }

protected:

  inline bool& ok () { return _ok; }

public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr clone() const { return new_ptr(*this); }

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const { return new_ptr(*this); }
  //@}

private:

  /// the event to be analysed
  tcEventPtr _lastEvent;

  /// stream state
  bool _ok;

  bool _doneEvent;

  /// last final state momenta
  vector<Lorentz5Momentum> _lastMomenta;

  /// last weight
  double _lastWeight;

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static ClassDescription<EventExtractor> initEventExtractor;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  EventExtractor & operator=(const EventExtractor &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"
#include "ThePEG/Config/Pointers.h"

namespace ThePEG {

ThePEG_DECLARE_POINTERS(Analysis2::EventExtractor,EventExtractorPtr);

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of EventExtractor. */
template <>
struct BaseClassTrait<Analysis2::EventExtractor,1> {
  /** Typedef of the first base class of EventExtractor. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the EventExtractor class and the shared object where it is defined. */
template <>
struct ClassTraits<Analysis2::EventExtractor>
  : public ClassTraitsBase<Analysis2::EventExtractor> {
  /** Return a platform-independent class name */
  static string className() { return "Analysis2::EventExtractor"; }
  /**
   * The name of a file containing the dynamic library where the class
   * EventExtractor is implemented. It may also include several, space-separated,
   * libraries if the class EventExtractor depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "Analysis2.so"; }
};

/** @endcond */

}

#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "EventExtractor.tcc"
#endif

#endif /* Analysis2_EventExtractor_H */
