#include "HSCPEvent.h"

#include <algorithm>
////////////////////////////////////////////////
////////////////////////////////////////////////
HSCPEvent::HSCPEvent(){

}
////////////////////////////////////////////////
////////////////////////////////////////////////
void HSCPEvent::clear(){

  runId = 0;
  eventId = 0;

}
////////////////////////////////////////////////
////////////////////////////////////////////////
HSCPParticle::HSCPParticle(const HSCPEvent *aEvent){

  clear();
  myEvent = aEvent;

}
////////////////////////////////////////////////                                                                                                                                            ////////////////////////////////////////////////                                                                                                                                           
void HSCPParticle::init(int index){

  pt = myEvent->pt[index];
  eta = myEvent->eta[index];
  eta = myEvent->phi[index];

}
////////////////////////////////////////////////                                                                                                                                            ////////////////////////////////////////////////
HSCPParticle HSCPEvent::getCandidate(unsigned int index) const{

  HSCPParticle aCandidate(this);
  aCandidate.clear();
  if((int)index>=numberOfCandidates) return aCandidate;

  aCandidate.init(index);

  return aCandidate;
}
////////////////////////////////////////////////
////////////////////////////////////////////////
void HSCPParticle::clear(){

  pt = -999;
  eta = -999;
  phi = -999;
 
}
////////////////////////////////////////////////
////////////////////////////////////////////////
