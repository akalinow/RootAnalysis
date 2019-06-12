enum class PropertyEnum {
PDGId = 0, 
charge = 1, 
decayMode = 2, 
discriminator = 3, 
muonID = 4, 
typeOfMuon = 5, 
byCombinedIsolationDeltaBetaCorrRaw3Hits = 6, 
byIsolationMVArun2v1DBoldDMwLTraw = 7, 
byIsolationMVArun2v1DBoldDMwLTraw2017v2 = 8, 
byIsolationMVArun2v1DBnewDMwLTraw2017v2 = 9, 
dxy = 10, 
dz = 11, 
SIP = 12, 
tauID = 13, 
combreliso = 14, 
leadChargedParticlePt = 15, 
isGoodTriggerType = 16, 
FilterFired = 17, 
L3FilterFired = 18, 
L3FilterFiredLast = 19, 
mc_match = 20, 
rawPt = 21, 
area = 22, 
PUJetID = 23, 
jecUnc = 24, 
Flavour = 25, 
bDiscriminator = 26, 
bCSVscore = 27, 
PFjetID = 28, 
NONE = 29
};

#include <map>
class PropertyEnumString { 
public:
static const std::map<std::string, PropertyEnum> enumMap;
};
