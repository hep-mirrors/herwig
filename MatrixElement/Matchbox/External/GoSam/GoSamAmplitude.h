// -*- C++ -*-
//
// GoSamAmplitude.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_GoSamAmplitude_H
#define Herwig_GoSamAmplitude_H
//
// This is the declaration of the GoSamAmplitude class.
//

#include "Herwig/MatrixElement/Matchbox/Base/MatchboxOLPME.h"
#include "ThePEG/Utilities/DynamicLoader.h"


namespace Herwig {

using namespace ThePEG;


class gosamprocinfo{

	public:
		gosamprocinfo(){};
		gosamprocinfo(int HID,int GID, string procstr,string typestr):
			theHOlpId(HID),theGOlpId(GID),theProcstr(procstr),theTypestr(typestr){
		}
		~gosamprocinfo(){}
		int HID() const {return theHOlpId;}
		int GID() const {return theGOlpId;}
		string Pstr() const {return theProcstr;}
		string Tstr() const {return theTypestr;}
		void setGID(int g){theGOlpId=g;}
		void setOAs(int i){ orderAlphas=i;}
		int orderAs(){return orderAlphas;}
		void setOAew(int j){ orderAlphaew=j;}
		int orderAew(){return orderAlphaew;}
	private:
		int theHOlpId;
		int theGOlpId;
		string theProcstr;
		string theTypestr;
		int orderAlphas;
		int orderAlphaew;
	public:
		void persistentOutput(PersistentOStream & os) const{os<<theHOlpId<<theGOlpId<<theProcstr<<theTypestr<<orderAlphas<<orderAlphaew;}
		void persistentInput(PersistentIStream &is) {is>>theHOlpId>>theGOlpId>>theProcstr>>theTypestr>>orderAlphas>>orderAlphaew;}
};


/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief GoSamAmplitude implements an interface to GoSam
 */
class GoSamAmplitude: public MatchboxOLPME {

public:

  /** @name Standard constructors and destructors. */
  //@{

  /**
   * The default constructor.
   */
  GoSamAmplitude();

  /**
   * The destructor.
   */
  virtual ~GoSamAmplitude();

  //@}


public:

  virtual void fillOrderFile(const map<pair<Process,int>,int>& procs, string OrderFileName);
  
  virtual bool isCS() const { return false; }
  virtual bool isExpanded() const { return true; }
  virtual bool isBDK() const { return false; }
  virtual bool isDR() const { return isitDR; }
  virtual bool isDRbar() const {return false;}

  /**
   * Start the one loop provider, if appropriate, giving order and
   * contract files
   */
  virtual void signOLP(const string&, const string&);

  virtual bool checkOLPContract(string contractFileName);

  /**
   * Return true, if this amplitude already includes symmetry factors
   * for identical outgoing particles.
   */
  // virtual bool hasFinalStateSymmetry() const { return true; }
  virtual bool hasFinalStateSymmetry() const { return false; }
  
  virtual bool buildGoSam();

  /**
   * Start the one loop provider, if appropriate
   */
  virtual void startOLP(const string&, int& status);

  /**
   * Return the value of the dimensional regularization
   * parameter. Note that renormalization scale dependence is fully
   * restored in DipoleIOperator.
   */
  // virtual Energy2 mu2() const { return lastSHat(); }

  /**
   * Start the one loop provider, if appropriate. This default
   * implementation writes an BLHA 2.0 order file and starts the OLP
   */
  virtual bool startOLP(const map<pair<Process,int>,int>& procs);

  /**
   * Call OLP_EvalSubProcess and fill in the results
   */
  void evalSubProcess() const;

  /**
   * Fill in results for the given colour correlator
   */
  virtual void evalColourCorrelator(pair<int,int> ij) const;

  /**
   * Return a positive helicity polarization vector for a gluon of
   * momentum p (with reference vector n) to be used when evaluating
   * spin correlations.
   */
  virtual LorentzVector<Complex> plusPolarization(const Lorentz5Momentum& p,
						  const Lorentz5Momentum& n, int inc) const;

  /**
   * Fill in results for the given colour/spin correlator
   */
  virtual void evalSpinColourCorrelator(pair<int,int> ij) const;

  void getids() const;

  int accuracyTargetNegExp() const { return theAccuracyTarget; }

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

  /** @name Standard Interfaced functions. */
  //@{

  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}


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
  GoSamAmplitude & operator=(const GoSamAmplitude &) = delete;

  /**
   * Store colour correlator results
   */
  mutable vector<double> colourCorrelatorResults;

  /**
   * Store spin colour correlator results
   */
  mutable vector<double> spinColourCorrelatorResults;

  /**
   * first is the olp id from herwig, second the answer from gosam
   */
  mutable vector<int> idpair;

  /**
   * first is the olp id from herwig, second the amplitude type
   */
  mutable vector<string> idtypepair;

   /**
   * Map to store all processes handled by this Amplitude
   */
  map<int , gosamprocinfo > processmap;

  mutable string gosamPathInterface;
  mutable string gosamSetupInFileNameInterface;
  mutable string gosamBuildScript;

  mutable string gosamPath;
  mutable string gosamSourcePath;
  mutable string gosamInstallPath;

  mutable string gosamSetupInFileName;

  mutable string orderFileTitle;
  mutable string contractFileTitle;
  mutable string contractFileName;
  mutable string orderFileName;

  mutable string accuracyFileTitle;
  mutable string accuracyFile;
  int theAccuracyTarget;

  bool theCodeExists;
  bool theFormOpt;
  bool theNinja;
  bool theHiggsEff;
  bool theMassiveLeptons;

  int theLoopInducedOption;

  bool isitDR;

  mutable bool doneGoSamInit;
  mutable bool doneGoSamInitRun;
  
  /**
   * The PDG codes of those quarks with mass
   */
  vector<int> massiveParticles; //theMassiveParticles;

  /**
   * Method to create the setup.in file for GoSam
   */
  void setupGoSamIn(string setupGoSamInFile);

protected:

  /**
   *   Location of the installed executables
   */
  string bindir_;

  /**
   *   Location of the data files
   */
  string pkgdatadir_;

  /**
   *  Location of GOSAM
   */
  string GoSamPrefix_;

}; // end "class GoSamAmplitude: public MatchboxOLPME"

  inline PersistentOStream& operator<<(PersistentOStream& os, const gosamprocinfo& p) {
    p.persistentOutput(os); return os;
  }

  inline PersistentIStream& operator>>(PersistentIStream& is, gosamprocinfo& p) {
    p.persistentInput(is); return is;
  }

} // end "namespace Herwig"

#endif /* Herwig_GoSamAmplitude_H */
