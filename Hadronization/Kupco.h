namespace Herwig {
using namespace ThePEG;

  /** \ingroup Hadronization
   *  \class Kupco
   *  \brief Class designed to make STL routines easy to use.
   *  \author Philip Stephens
   *
   *  This class is used to generate a list of the hadron pairs which can
   *  be produced that allows easy traversal and quick access.
   */
  class Kupco {

  public:

    /**
     *  Constructor
     * @param inidQ PDG code of the quark drawn from the vacuum.
     * @param inhad1 ParticleData for the first hadron produced.
     * @param inhad2 ParticleData for the second hadron produced.
     * @param inwgt  The weight for the hadron pair
     */
     Kupco(tcPDPtr inidQ,tcPDPtr inhad1,tcPDPtr inhad2, Energy inwgt)
      : idQ(inidQ),hadron1(inhad1),hadron2(inhad2),weight(inwgt)
    {}

    /**
     * id of the quark drawn from the vacuum.
     */
    tcPDPtr idQ;

    /**
     * The ParticleData object for the first hadron produced.
     */
    tcPDPtr hadron1;

    /**
     * The ParticleData object for the second hadron produced.
     */
    tcPDPtr hadron2;

    /**
     * Weight factor of this componation.
     */
    Energy weight;
  };

  /**
   *  Type defs 
   */
  //@{
  /**
   * The type is used to contain all the hadrons info of a given flavour.
   */
  typedef set<HadronInfo> KupcoData;
  //@}

  /**
   * The hadron table type.
   */
  typedef map<pair<long,long>,KupcoData> HadronTable;
}
