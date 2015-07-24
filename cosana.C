#define cosana_cxx
#include "cosana.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

/*--------------------------------------------------------------------------------*/
void cosana::make_tracks(){
   Int_t nhit = (Int_t)xp->size();
   // nested loop over hits to find pairs of hits in layer 0+1 or in layer 2+3
   for(int ih = 0; ih < nhit-1; ih++){
      Int_t ilayer = (Int_t)collid->at(ih);
      for(int jh = ih+1; jh < nhit; jh++){
        Int_t jlayer = (Int_t)collid->at(jh);
        Int_t itype  = -1;
        if( ( ilayer==0 && jlayer==1) || (ilayer==1 && jlayer==0) ) itype = BOTTOM;
        if( ( ilayer==2 && jlayer==3) || (ilayer==3 && jlayer==2) ) itype = TOP;
        // did we find a good track pair?
        if(itype != -1) {
           track t(TVector3(xp->at(ih),yp->at(ih),zp->at(ih)),TVector3(xp->at(jh),yp->at(jh),zp->at(jh)),itype);
           tracks.push_back(t);
        }
      } 
   }
}
/*--------------------------------------------------------------------------------*/
void cosana::make_intersects(){ 
   Int_t ntrack = (Int_t)tracks.size();
   for(int it = 0; it<ntrack-1; it++){
     Int_t itype = tracks.at(it).getType();
     for(int jt=it+1; jt<ntrack; jt++){
        Int_t jtype = tracks.at(jt).getType();
        if((itype==BOTTOM && jtype==TOP) || (itype==TOP && jtype==BOTTOM)){
           intersect xing(tracks.at(it),tracks.at(jt));
           intersects.push_back(xing);
        }
     }
   }
}
/*--------------------------------------------------------------------------------*/
void cosana::fill_histograms(){
   for(Int_t ix=0; ix<(Int_t)intersects.size(); ix++){
      intersect xing = intersects.at(ix);
      Double_t xx = xing.getPosition().x();
      Double_t yy = xing.getPosition().y();
      Double_t zz = xing.getPosition().z();
      _ip0_yz->Fill(yy,zz);
      if(xing.getDist()<1.0 && xing.getCos()<0.999) _ip0_yz_cut->Fill(yy,zz);
   }
}
/*--------------------------------------------------------------------------------*/
void cosana::book_histograms(){

   _f = new TFile(outFile.c_str(),"RECREATE");
//   TH1F *_nhit = new TH1F("nhit","nhit",10,-0.5,9.5);
//   TH1F *_cost = new TH1F("cost","cost",1000,-1,1);
//   TH1F *_dist = new TH1F("dist","dist",1000,0.,10.);
//   TH2F *_xy0 = new TH2F("xy0","xy0",500,-500,500,500,-500,500);
//   TH2F *_xy1 = new TH2F("xy1","xy1",500,-500,500,500,-500,500);
//   TH2F *_xy2 = new TH2F("xy2","xy2",500,-500,500,500,-500,500);
//   TH2F *_xy3 = new TH2F("xy3","xy3",500,-500,500,500,-500,500);

   _ip0_yz     = new TH2F("yz0","yz0",100,-200,200,100,-200,200);
   _ip0_yz_cut = new TH2F("yz0_cut","yz0_cut",100,-200,200,100,-200,200);

}
/*--------------------------------------------------------------------------------*/
void cosana::reset_variables(){
   tracks.clear();
   intersects.clear();
}
/*--------------------------------------------------------------------------------*/
void cosana::Loop()
{
   //
   // book histograms
   //
   book_histograms();
   //
   // event loop
   //
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      //
      // read an event from the event TTree/TChain
      //
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if(jentry%1000000 == 0) cout << "processed "<<jentry<<" events"<<endl;
      //
      // reset variables
      // 
      reset_variables();
      //
      // make the track elements
      //
      make_tracks();
      //
      // make intersects
      //
      make_intersects();
      //
      // fill histograms
      //
      fill_histograms();
   } // end loop over events
   _f->Write();
}
