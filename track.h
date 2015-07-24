#ifndef __TRACK_H__
#define __TRACK_H__

#include "TVector3.h"

#define TOP    0
#define BOTTOM 1

class track{
public:
   track(TVector3 x0, TVector3 x1, Int_t i);
   TVector3 getPosition(){return pos;};
   TVector3 getDirection(){return dir;};
   Int_t getType(){return itype;};
private:
   TVector3 pos;
   TVector3 dir;
   Int_t    itype;
};

#endif
