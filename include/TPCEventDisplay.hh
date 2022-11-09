#include "TPCManager.hh"
class TPCMainFrame : public TPCManager {
   RQ_OBJECT("TPCMainFrame")

public:
  TPCMainFrame(const TGWindow *p,UInt_t w,UInt_t h);
  virtual ~TPCMainFrame();
  void DoDraw();
  Bool_t InitializeTPC();
  Bool_t loadFile(){return true;};
  Bool_t loadEvent();

  void GoForward();
  void GoBackward();
  void GoTo();

private:
  TGMainFrame         *fMain;
  TRootEmbeddedCanvas *fEcanvas;
  TCanvas *fCanvas;
  TGeoVolume *fTPC;
  TH2Poly             *fTPC2dPoly;
  TH2D             *fTPCxz;
  TH2D             *fTPCxz_raw;
  TH2D             *fTPCxz_cluster;
  TGTextEntry      *fStatusBar;
  TGTextEntry      *fEvtHandler;
  Int_t eventID;
};


TPCMainFrame::TPCMainFrame(const TGWindow *p,UInt_t w,UInt_t h) {
   // Create a main frame
   fMain = new TGMainFrame(p,w,h);

   // Create canvas widget
   fEcanvas = new TRootEmbeddedCanvas("Ecanvas",fMain,1600,800);
   fMain->AddFrame(fEcanvas, new TGLayoutHints(kLHintsExpandX |
                   kLHintsExpandY, 10,10,10,1));
   
   fTPC2dPoly = new TH2Poly("fTPC2dPoly","fTPC2dPoly",-300,300,-300,300);
   fTPCxz = new TH2D("fTPCxz","fTPCxz",300,-300,300,300,-300,300);
   fTPCxz_cluster = new TH2D("fTPCxz_cluster","fTPCxz_cluster",300,-300,300,300,-300,300);
   fTPCxz_raw = new TH2D("fTPCxz_raw","fTPCxz_raw",300,-300,300,300,-300,300);

   fTPC = TPCGeometry();

   // Create a horizontal frame widget with buttons
   TGHorizontalFrame *hframe1 = new TGHorizontalFrame(fMain,1600,40);
   TString iconDir(TString::Format("%s/icons",gSystem->Getenv("ROOTSYS")));
   TGPictureButton *bwd = new TGPictureButton(hframe1, gClient->GetPicture(iconDir +"GoBack.gif"));
   bwd->Connect("Clicked()","TPCMainFrame",this, "GoBackward()");
   hframe1->AddFrame(bwd,new TGLayoutHints(kLHintsCenterX | kLHintsExpandX,
					  2,2,2,2));

   TGPictureButton *fwd = new TGPictureButton(hframe1, gClient->GetPicture(iconDir +"GoFoward.gif"));
   fwd->Connect("Clicked()","TPCMainFrame",this, "GoForward()");
   hframe1->AddFrame(fwd, new TGLayoutHints(kLHintsCenterX | kLHintsExpandX,
					  2,2,2,2));

   TGTextButton *shortcut = new TGTextButton(hframe1);
   shortcut->Connect("Clicked()","TPCMainFrame",this, "GoTo()");
   hframe1->AddFrame(shortcut, new TGLayoutHints(kLHintsLeft | kLHintsExpandX | kLHintsExpandY,
                                            2,2,2,2));

   fEvtHandler = new TGTextEntry(hframe1);
   fEvtHandler->SetEnabled(true);
   fEvtHandler->SetInsertMode();
   hframe1->AddFrame(fEvtHandler,new TGLayoutHints(kLHintsCenterX | kLHintsExpandX | kLHintsExpandY,
						     2,2,2,2));
   
   TGTextButton *exit = new TGTextButton(hframe1,"&Exit",
                                "gApplication->Terminate(0)");
   hframe1->AddFrame(exit, new TGLayoutHints(kLHintsLeft | kLHintsExpandX | kLHintsExpandY,
                                            2,2,2,2));
   TGHorizontalFrame *hframe2 = new TGHorizontalFrame(fMain,1600,40);

   fStatusBar = new TGTextEntry(hframe2);
   fStatusBar->SetEnabled(kFALSE);
   hframe2->AddFrame(fStatusBar,new TGLayoutHints(kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2)); 

   fMain->AddFrame(hframe1, new TGLayoutHints(kLHintsCenterX | kLHintsExpandX,
                                             2,2,2,2));
   fMain->AddFrame(hframe2, new TGLayoutHints(kLHintsCenterX| kLHintsExpandX ,
                                             2,2,2,2));
   // Set a name to the main frame
   fMain->SetWindowName("Simple Example");

   // Map all subwindows of main frame
   fMain->MapSubwindows();

   // Initialize the layout algorithm
   fMain->Resize(fMain->GetDefaultSize());

   // Map main frame
   fMain->MapWindow();

   fCanvas = fEcanvas->GetCanvas();
   fCanvas->Divide(2);
//   loadEvent();	
		DoDraw();
/*
   if (loadFile())
*/
}


Bool_t TPCMainFrame::InitializeTPC() {

   Double_t X[5]; Double_t Y[5];
   
   for (int l = 0 ; l < 32; l++) {
     Double_t pLength = tpc::padParameter[l][5];
     Double_t st = (180. - ( 360./tpc::padParameter[l][3]) * tpc::padParameter[l][1]/2.);
     Double_t sTheta = (-1 + st/180.) * TMath::Pi();
     Double_t dTheta = (360./tpc::padParameter[l][3])/180.*TMath::Pi();
     Double_t cRad = tpc::padParameter[l][2];
     Double_t nPad = tpc::padParameter[l][1];

     for (int j = 0 ; j < nPad ; j++) {

       X[1] = (cRad+(pLength/2.))*TMath::Cos(j*dTheta+sTheta);
       X[2] = (cRad+(pLength/2.))*TMath::Cos((j+1)*dTheta+sTheta);
       X[3] = (cRad-(pLength/2.))*TMath::Cos((j+1)*dTheta+sTheta);
       X[4] = (cRad-(pLength/2.))*TMath::Cos(j*dTheta+sTheta);
       X[0] = X[4];
       Y[1] = (cRad+(pLength/2.))*TMath::Sin(j*dTheta+sTheta);
       Y[2] = (cRad+(pLength/2.))*TMath::Sin((j+1)*dTheta+sTheta);
       Y[3] = (cRad-(pLength/2.))*TMath::Sin((j+1)*dTheta+sTheta);
       Y[4] = (cRad-(pLength/2.))*TMath::Sin(j*dTheta+sTheta);
       Y[0] = Y[4];

       for (int k = 0 ; k < 5 ; k++)
	 X[k] += tpc::ZTarget;

       fTPC2dPoly->AddBin(5,X,Y);
     }
   }

   fTPC2dPoly->SetStats(0);

   return true;
}

TPCMainFrame::~TPCMainFrame() {
   // Clean up used widgets: frames, buttons, layout hints  
   fMain->Cleanup();
   delete fMain;
}



void TPCMainFrame::DoDraw() {
   // Draws function graphics in randomly chosen interval
/*
	fCanvas = fEcanvas->GetCanvas();
  fCanvas->Divide(2);
  fCanvas->cd(1);
  
  fTPC2dPoly->Draw("col");
  gPad->SetLogz(1);
  fCanvas->Update();
  fCanvas->Modified();

  fCanvas->cd(2);
  fTPC->Draw("");
  fCanvas->Update();
  fCanvas->Modified();
*/
	loadEvent();
}

Bool_t TPCMainFrame::loadEvent(){
 
	fTPC2dPoly->Reset("");
  fTPCxz->Reset("");
  fTPCxz_raw->Reset("");
  fTPCxz_cluster->Reset("");
  fCanvas->cd(1);
	
  fTPC2dPoly->Draw("col");
  fTPCxz_raw->Draw("SAME");
  fTPCxz_raw->SetMarkerSize(1);
  fTPCxz_raw->SetMarkerStyle(2);
  fTPCxz_raw->SetMarkerColor(kOrange-2);

  fTPCxz->Draw("SAME");
  fTPCxz->SetMarkerStyle(8);
  fTPCxz->SetMarkerColor(kRed);

  fTPCxz_cluster->Draw("SAME");
  fTPCxz_cluster->SetMarkerStyle(8);
  fTPCxz_cluster->SetMarkerColor(kSpring);
  fCanvas->Modified();
  fCanvas->Update();
  
	fCanvas->cd(2);
  TView3D *view = (TView3D*) TView::CreateView(1);
  fTPC->Draw();
//  view->ZoomView(gPad,1.5);

  fCanvas->Update();
  fCanvas->Modified();

	return true;
}

void TPCMainFrame::GoForward() {

  fCanvas->Clear("D");
  /*
  fTPC2dPoly->Reset("");
  fTPCxz->Reset("");
  fTPCxz_raw->Reset("");
  fTPCxz_cluster->Reset("");
  */
  //  fTPC->CleanAll();
 	cout<<"GoingForward"<<endl; 
 	loadEvent();
	if (eventID < 100) {
    eventID++;
//    loadEvent();
  }
  else {
    fStatusBar->SetTextColor(0xff0000);
    fStatusBar->SetText("Already at the last event");
    printf("Already at the last event.\n");
  }
}

void TPCMainFrame::GoBackward() {

  fCanvas->Clear("D");

  //  fTPC->CleanAll();
 	cout<<"GoingBackward"<<endl; 
 	loadEvent();
  
  if (eventID > 0 ) {
    eventID--;
 //   loadEvent();
  }
  else {
    fStatusBar->SetTextColor(0xff0000);
    fStatusBar->SetText("Already at the first event");
    printf("Already at the first event.\n");
  }
}

void TPCMainFrame::GoTo() {

  fCanvas->Clear("D");

  //  fTPC->CleanAll();
  
  TString evnum = fEvtHandler->GetMarkedText();
	int evi = evnum.Atoi();
  if (evi >=0 ) {
    eventID = evi;
		cout<<"Event : "<<evi<<endl;
		//    loadEvent();
  }
  else {
    fStatusBar->SetTextColor(0xff0000);
    fStatusBar->SetText("Out of range");
    printf("Out of range\n");
  }
}





