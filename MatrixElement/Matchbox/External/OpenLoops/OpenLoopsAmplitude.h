// -*- C++ -*-
//
// OpenLoopsAmplitude.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_OpenLoopsAmplitude_H
#define Herwig_OpenLoopsAmplitude_H
//
// This is the declaration of the OpenLoopsAmplitude class.
//

#include "Herwig++/MatrixElement/Matchbox/Base/MatchboxOLPME.h"
#include "ThePEG/Utilities/DynamicLoader.h"




namespace Herwig {

using namespace ThePEG;


class openloopsprocinfo{

	public:
		openloopsprocinfo(){};
		openloopsprocinfo(int HID,int GID, string procstr,string typestr):
			theHOlpId(HID),theGOlpId(GID),theProcstr(procstr),theTypestr(typestr){
		}
		~openloopsprocinfo(){}
		int HID() const {return theHOlpId;}
		int GID() const {return theGOlpId;}
		string Pstr() const {return theProcstr;}
		string Tstr() const {return theTypestr;}
		void setGID(int g){theGOlpId=g;}
		void setOAs(int i){ orderAlphas=i;}
		int orderAs(){return orderAlphas;}
	private:
		int theHOlpId;
		int theGOlpId;
		string theProcstr;
		string theTypestr;
		int orderAlphas;
	public:
		void persistentOutput(PersistentOStream & os) const{os<<theHOlpId<<theGOlpId<<theProcstr<<theTypestr<<orderAlphas;}
		void persistentInput(PersistentIStream &is) {is>>theHOlpId>>theGOlpId>>theProcstr>>theTypestr>>orderAlphas;}
};
/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief OpenLoopsAmplitude implements an interface to NJets
 */
class OpenLoopsAmplitude: public MatchboxOLPME {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  OpenLoopsAmplitude();

  /**
   * The destructor.
   */
  virtual ~OpenLoopsAmplitude();
  //@}

public:

  virtual void fillOrderFile(const map<pair<Process,int>,int>& procs);


  virtual bool isCS() const { return false; }
  virtual bool isExpanded() const { return true; }
  virtual bool isBDK() const { return false; }
  //virtual bool isDR() const { return true; }
  /**
   * Start the one loop provider, if appropriate, giving order and
   * contract files
   */

  virtual bool checkOLPContract();

  /**
   * Start the one loop provider, if appropriate
   */
  virtual void startOLP(const string&, int& status);
  virtual void getids() const ;

  virtual Energy2 mu2() const { return lastSHat(); }
  /**
   * Start the one loop provider, if appropriate. This default
   * implementation writes an BLHA 2.0 order file and starts the OLP
   */
  virtual bool startOLP(const map<pair<Process,int>,int>& procs);


    /**
   * Return true, if this amplitude already includes averaging over
   * incoming parton's quantum numbers.
   */
  virtual bool hasInitialAverage() const { return true; }

  /**
   * Return true, if this amplitude already includes symmetry factors
   * for identical outgoing particles.
   */
  virtual bool hasFinalStateSymmetry() const { return true; }

  /**
   * Call OLP_EvalSubProcess and fill in the results
   */
   void evalSubProcess() const;

  /**
   * Fill in results for the given colour correlator
   */
  virtual void evalColourCorrelator(pair<int,int> ij) const;


  /**
   * Fill in results for the given colour/spin correlator
   */
  virtual void evalSpinColourCorrelator(pair<int,int> ij) const;


  /**
   * Return the colour and spin correlated matrix element.
   */
  virtual double spinColourCorrelatedME2(pair<int,int> ij,
					 const SpinCorrelationTensor& c) const;


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

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  OpenLoopsAmplitude & operator=(const OpenLoopsAmplitude &);

  /**
   * Store colour correlator results
   */
  mutable vector<double> colourCorrelatorResults;

  /**
   * Store spin colour correlator results
   */
  mutable vector<double> spinColourCorrelatorResults;


  /**
   * first is the olp id from herwig, second the answer from openloops
   */
  mutable map< int , int > idpair;

  map<int , openloopsprocinfo > processmap;

  mutable string openloopsInstallPath;

  bool theCodeExists;

  mutable string extraOpenLoopsPath;


};
//inline PersistentOStream& operator<<(PersistentOStream& os,
//				     const openloopsprocinfo& h) {
//  h.persistentOutput(os);
//  return os;
//}

//inline PersistentIStream& operator>>(PersistentIStream& is,
//		openloopsprocinfo& h) {
//  h.persistentInput(is);
//
//return is;
//}
}

#endif /* Herwig_OpenLoopsAmplitude_H */
