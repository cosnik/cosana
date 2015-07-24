#ifndef __INTERSECT_H__
#define __INTERSECT_H__

#include "track.h"

class intersect{
public:
   intersect(track t0, track t1);
   TVector3 getPosition(){return xint;};
   Double_t getDist(){return dist;};
   Double_t getCos(){return cos_psi;};
private:
   TVector3 xint;
   Double_t dist;
   Double_t cos_psi;
};

#endif
