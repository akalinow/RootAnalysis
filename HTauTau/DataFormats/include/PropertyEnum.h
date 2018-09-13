//! File containing an enum type of particle properties and an interface to get them using std::string names.
/*!
  \author Rafal Maselek && Artur Kalinowski
  \date May 2018
  
  This file contains an enum type declaration, which enumerates 29 different properties of particles defined
  by HTTParticle object type. In order to access them via parsing an external file, class PropertyEnumString
  was added with a static member std::map -- an interface to get PropertyEnum by using std::string name of
  property. 
  See also PropertyEnum.cc
*/


#include <map>

enum class PropertyEnum { PDGId = 0, 
charge = 1, 
decayMode = 2, 
discriminator = 3, 
muonID = 4, 
typeOfMuon = 5, 
byCombinedIsolationDeltaBetaCorrRaw3Hits = 6, 
byIsolationMVArun2v1DBoldDMwLTraw = 7, 
againstElectronMVA5category = 8, 
dxy = 9, 
dz = 10, 
SIP = 11, 
tauID = 12, 
combreliso = 13, 
leadChargedParticlePt = 14, 
isGoodTriggerType = 15, 
FilterFired = 16, 
L3FilterFired = 17, 
L3FilterFiredLast = 18, 
mc_match = 19, 
rawPt = 20, 
area = 21, 
PUJetID = 22, 
jecUnc = 23, 
Flavour = 24, 
bDiscriminator = 25, 
bCSVscore = 26, 
PFjetID = 27, 
NONE = 28
};

// Defined in PropertyEnum.cc
class PropertyEnumString
{
	public:
    	static const std::map<std::string, PropertyEnum> enumMap;
};

