namespace Herwig {
using namespace ThePEG;

/** \ingroup Hadronization
   *  \class HadronInfo
   *  \brief Class used to store all the hadron information for easy access.
   *  \author Philip Stephens
   *
   *  Note that:
   *  - the hadrons in _table can be filled in any ordered
   *    w.r.t. the mass value, and flavours for different
   *    groups (for instance, (u,s) hadrons don't need to
   *    be placed after (d,s) or any other flavour), but
   *    all hadrons with the same flavours must be consecutive
   *    ( for instance you cannot alternate hadrons of type
   *    (d,s) with those of flavour (u,s) ).
   *    Furthermore, it is assumed that particle and antiparticle
   *    have the same weights, and therefore only one of them
   *    must be entered in the table: we have chosen to refer
   *    to the particle, defined as PDG id > 0, although if
   *    an anti-particle is provided in input it is automatically
   *    transform to its particle, simply by taking the modulus
   *    of its id.
   */
  class HadronInfo {

  public:

    /**
     *  Constructor
     * @param idin The PDG code of the hadron
     * @param datain The pointer to the ParticleData object
     * @param swtin  The singlet/decuplet/orbital factor
     * @param massin The mass of the hadron
     */
    HadronInfo(long idin=0, tPDPtr datain=tPDPtr(),
	       double swtin=1., Energy massin=ZERO)
      : id(idin), ptrData(datain), swtef(swtin), wt(1.0), overallWeight(0.0),
	mass(massin)
    {}

    /**
     *  Comparision operator on mass
     */
     bool operator<(const HadronInfo &x) const {
      if(mass!=x.mass) return mass < x.mass;
      else return id < x.id;
    }

    /**
     * The hadrons id.
     */
    long  id;

    /**
     * pointer to ParticleData, to get the spin, etc...
     */
    tPDPtr ptrData;

    /**
     * singlet/decuplet/orbital factor
     */
    double swtef;

    /**
     * mixing factor
     */
    double wt;

    /**
     * (2*J+1)*wt*swtef
     */
    double overallWeight;

    /**
     * The hadrons mass
     */
    Energy mass;

    /**
     *  Rescale the weight for a given hadron
     */
    void rescale(double x) const {
      const_cast<HadronInfo*>(this)->overallWeight *= x;
    }

    /**
     * Friend method used to print the value of a table element.
     */
    friend PersistentOStream & operator<< (PersistentOStream & os,
					   const HadronInfo & hi ) {
      os << hi.id << hi.ptrData << hi.swtef << hi.wt
	 << hi.overallWeight << ounit(hi.mass,GeV);
      return os;
    }

    /**
     * debug output
     */
    friend ostream & operator<< (ostream & os, const HadronInfo & hi ) {
      os << std::scientific << std::showpoint
	 << std::setprecision(4)
	 << setw(2)
	 << hi.id << '\t'
	 << hi.swtef << '\t'
	 << hi.wt << '\t'
	 << hi.overallWeight << '\t'
	 << ounit(hi.mass,GeV);
      return os;
    }

    /**
     * Friend method used to read in the value of a table element.
     */
    friend PersistentIStream & operator>> (PersistentIStream & is,
					   HadronInfo & hi ) {
      is >> hi.id >> hi.ptrData >> hi.swtef >> hi.wt
	 >> hi.overallWeight >> iunit(hi.mass,GeV);
      return is;
    }
  };
}
