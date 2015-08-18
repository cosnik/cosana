#include "track.h"

track::track(TVector3 x0, TVector3 x1, Int_t i){
   // track classification
   itype = i;

   // position on the track
   pos = x0;

   // direction of the track
   TVector3 delta = x1 - x0;
   dir = delta.Unit();

}
