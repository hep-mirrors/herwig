#include "HerwigVersion.h"

std::string Herwig::HerwigVersion::versionstring = "";

#ifdef HERWIG_DATADIR
std::string Herwig::HerwigVersion::pkgdatadir = HERWIG_DATADIR;
#else
#error "HERWIG_DATADIR must be defined."
#endif
