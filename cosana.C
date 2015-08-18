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
           int ihit=0;
           int jhit=0;
           // always make sure that the start position of the track is measured in the same silicon layer
           // this may be needed when making a sinogram
           if(ilayer > jlayer) {
             ihit = ih;
             jhit = jh;
           } else {
             ihit = jh;
             jhit = ih;
           }
           track t(TVector3(xp->at(ihit),yp->at(ihit),zp->at(ihit)),TVector3(xp->at(jhit),yp->at(jhit),zp->at(jhit)),itype);
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
void cosana::make_sinogram(){ 
   // find events with one track in the top detector and one in the bottom detector
   Int_t ntrack = (Int_t)tracks.size();
   if (ntrack != 2 ) return;

   for(int it = 0; it<ntrack-1; it++){
     Int_t itype = tracks.at(it).getType();
     for(int jt=it+1; jt<ntrack; jt++){
        Int_t jtype = tracks.at(jt).getType();
        // one bottom track & one top track
        if((itype==BOTTOM && jtype==TOP) || (itype==TOP && jtype==BOTTOM)){
          // get the angle between ingoing and outgoing track 
          Double_t cos_psi = fabs(tracks.at(it).getDirection()*tracks.at(jt).getDirection());
          // get the index of the top track
          int itop=-1;
          if(itype==TOP) {
             itop = it;
          } else {
             itop = jt;
          }
          // get the angle wrt +y-axis
          TVector3 ydir(0,1,0);
          TVector3 yzdir(0,tracks.at(itop).getDirection().y(), tracks.at(itop).getDirection().z());
          yzdir = yzdir.Unit();
          Double_t cos_alfa = yzdir*ydir;
          Double_t alfa = acos(cos_alfa);
          // get the start position of the top track
          Double_t y0 = tracks.at(itop).getPosition().y();
          //
          if(cos_psi < 0.9995) {
              _sinogram->Fill(alfa,y0,0.);
          } else {
              _sinogram->Fill(alfa,y0,1.);
          }

           
        }
     }
   }
}
/*--------------------------------------------------------------------------------*/
void cosana::fill_histograms(){
   _nhit->Fill((Double_t)xp->size());

   for(Int_t ix=0; ix<(Int_t)intersects.size(); ix++){
      
      intersect xing = intersects.at(ix);
      _dist->Fill(xing.getDist());
      _cost->Fill(xing.getCos());

      Double_t xx = xing.getPosition().x();
      Double_t yy = xing.getPosition().y();
      Double_t zz = xing.getPosition().z();
      _ip0_yz->Fill(yy,zz);
//      if(xing.getDist()<0.2 && xing.getCos()<0.99995) {
      if(xing.getDist()<1.0 && xing.getCos()<0.999) {
        _ip0_yz_cut->Fill(yy,zz);
      }
   }
}
/*--------------------------------------------------------------------------------*/
void cosana::book_histograms(){

   _f = new TFile(outFile.c_str(),"RECREATE");

   _nhit = new TH1F("nhit","nhit",10,-0.5,9.5);
   _cost = new TH1F("cost","cost",1000,-1,1);
   _dist = new TH1F("dist","dist",1000,0.,10.);
   _xy0 = new TH2F("xy0","xy0",500,-500,500,500,-500,500);
   _xy1 = new TH2F("xy1","xy1",500,-500,500,500,-500,500);
   _xy2 = new TH2F("xy2","xy2",500,-500,500,500,-500,500);
   _xy3 = new TH2F("xy3","xy3",500,-500,500,500,-500,500);

   _testit   = new TH2F("testit","testit",500,0.,3.1415,500,-500,500);
   _sinogram = new TProfile2D("sino","sino",500,0.,3.1415,500,-500,500,0.,2);

   _ip0_yz     = new TH2F("yz0","yz0",100,-200,200,100,-200,200);
   _ip0_yz_cut = new TH2F("yz0_cut","yz0_cut",100,-200,200,100,-200,200);

//   _out_tree = new TTree("ana");
   //_out_tree->Branch("x0",&_t_x0,"x0/D");
   //_out_tree->Branch("y0",&_t_y0,"y0/D");
   //_out_tree->Branch("z0",&_t_z0,"z0/D");
   //_out_tree->Branch("tx0",&_t_tx0,"tx/D");
   //_out_tree->Branch("ty0",&_t_ty0,"ty/D");
   //_out_tree->Branch("tz0",&_t_tz0,"tz/D");
   //_out_tree->Branch("tx1",&_t_tx1,"tx/D");
   //_out_tree->Branch("ty1",&_t_ty1,"ty/D");
   //_out_tree->Branch("tz1",&_t_tz1,"tz/D");
   //_out_tree->Branch("dth",&_t_dth,"dth/D");

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
      // make a sinogram
      //
      make_sinogram();
      //
      // fill histograms
      //
      fill_histograms();
      //
      //
      //
   } // end loop over events
   _f->Write();
}
