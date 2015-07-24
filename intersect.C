#include "intersect.h"

#include "TMatrixD.h"
#include "TVectorD.h"

intersect::intersect(track track0, track track1){
  // calculate the intersect of two tracks

  TVector3 delta = track0.getPosition() - track1.getPosition();
  cos_psi = track0.getDirection() * track1.getDirection();
  
  TMatrixD mm(2,2);
  mm(0,0) = 1.;
  mm(0,1) = -cos_psi;
  mm(1,0) = -cos_psi;
  mm(1,1) = 1.;   

  TVectorD gg(2);
  gg(0) = -delta*track0.getDirection();
  gg(1) =  delta*track1.getDirection();

  TVectorD solution(2);

  TDecompLU lu(mm);
  Bool_t ok;
  solution = lu.Solve(gg,ok);

  TVector3 s0 = track0.getPosition()+solution(0)*track0.getDirection();
  TVector3 s1 = track1.getPosition()+solution(1)*track1.getDirection();
  xint = (s0+s1)*0.5;
  dist = (s0-s1).Mag();
  
}
