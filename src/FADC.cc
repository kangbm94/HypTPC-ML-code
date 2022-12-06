#include "../include/FADC.hh"
#include "TPCManager.cc"
vector<int>fpads = {60 ,72 ,73 ,74 ,111 ,135 ,136 ,186 ,223 ,333 ,468 ,627 ,810 ,922 ,1017 ,1141 ,1364 ,1380 ,1381 ,1382 ,1383 ,1384 ,1385 ,1386 ,1387 ,1389 ,1390 ,1443 ,1449 ,1451 ,1453 ,1454 ,1504 ,1507 ,1508 ,1511 ,1512 ,1513 ,1514 ,1515 ,1516 ,1527 ,1528 ,1529 ,1530 ,1565 ,1566 ,1585 ,1586 ,1587 ,1588 ,1589 ,1590 ,1591 ,1592 ,1661 ,1692 ,1704 ,1705 ,1714 ,1729 ,1730 ,1731 ,1752 ,1756 ,1758 ,1760 ,1762 ,1768 ,1779 ,1819 ,1820 ,1875 ,1876 ,1921 ,2037 ,2038 ,2097 ,2098 ,2099 ,2246 ,2247 ,2248 ,2308 ,2348 ,2353 ,2355 ,2356 ,2358 ,2359 ,2362 ,2368 ,2456 ,2457 ,2487 ,2488 ,2580 ,2582 ,2583 ,2584 ,2669 ,2670 ,2795 ,2796 ,2797 ,2798 ,2799 ,2800 ,2801 ,2802 ,2803 ,2804 ,2805 ,2806 ,2807 ,2808 ,2809 ,2810 ,2888 ,3003 ,3004 ,3005 ,3006 ,3007 ,3008 ,3009 ,3010 ,3011 ,3012 ,3013 ,3014 ,3015 ,3016 ,3017 ,3018 ,3019 ,3027 ,3038 ,3039 ,3040 ,3041 ,3042 ,3043 ,3044 ,3045 ,3046 ,3047 ,3048 ,3049 ,3050 ,3051 ,3052 ,3053 ,3054 ,3055 ,3112 ,3113 ,3184 ,3186 ,3187 ,3201 ,3203 ,3205 ,3207 ,3209 ,3210 ,3211 ,3212 ,3213 ,3215 ,3217 ,3228 ,3229 ,3230 ,3231 ,3232 ,3240 ,3254 ,3290 ,3291 ,3292 ,3293 ,3294 ,3295 ,3296 ,3297 ,3298 ,3343 ,3344 ,3414 ,3415 ,3416 ,3417 ,3418 ,3419 ,3420 ,3424 ,3427 ,3429 ,3442 ,3443 ,3444 ,3446 ,3447 ,3448 ,3449 ,3450 ,3451 ,3452 ,3454 ,3455 ,3456 ,3474 ,3475 ,3478 ,3479 ,3480 ,3481 ,3482 ,3484 ,3486 ,3487 ,3488 ,3490 ,3491 ,3492 ,3494 ,3495 ,3498 ,3499 ,3502 ,3503 ,3504 ,3510 ,3541 ,3542 ,3543 ,3544 ,3545 ,3546 ,3547 ,3581 ,3660 ,3661 ,3662 ,3663 ,3664 ,3665 ,3666 ,3693 ,3698 ,3815 ,3891 ,3904 ,3905 ,3906 ,3907 ,3908 ,3909 ,3933 ,3934 ,3935 ,3937 ,3939 ,4036 ,4067 ,4107 ,4135 ,4136 ,4137 ,4138 ,4158 ,4159 ,4160 ,4161 ,4162 ,4163 ,4164 ,4165 ,4166 ,4167 ,4168 ,4169 ,4170 ,4171 ,4172 ,4354 ,4355 ,4356 ,4357 ,4376 ,4378 ,4379 ,4380 ,4381 ,4382 ,4383 ,4384 ,4385 ,4386 ,4387 ,4388 ,4389 ,4411 ,4412 ,4413 ,4414 ,4415 ,4416 ,4417 ,4418 ,4419 ,4420 ,4421 ,4422 ,4423 ,4424 ,4425 ,4426 ,4427 ,4428 ,4429 ,4430 ,4431 ,4432 ,4433 ,4435 ,4436 ,4437 ,4438 ,4439 ,4440 ,4444 ,4448 ,4536 ,4566 ,4567 ,4568 ,4569 ,4587 ,4588 ,4589 ,4590 ,4591 ,4592 ,4593 ,4594 ,4595 ,4596 ,4597 ,4598 ,4613 ,4617 ,4619 ,4721 ,4722 ,4723 ,4724 ,4725 ,4726 ,4727 ,4728 ,4729 ,4736 ,4740 ,4741 ,4743 ,4747 ,4775 ,4776 ,4777 ,4790 ,4794 ,4795 ,4796 ,4798 ,4800 ,4804 ,4880 ,4930 ,4931 ,4932 ,4933 ,4934 ,4935 ,4936 ,4937 ,4941 ,4942 ,4982 ,4983 ,4998 ,5001 ,5003 ,5004 ,5006 ,5008 ,5009 ,5074 ,5113 ,5114 ,5135 ,5136 ,5137 ,5138 ,5139 ,5140 ,5141 ,5184 ,5185 ,5202 ,5205 ,5207 ,5209 ,5215 ,5216 ,5298 ,5328 ,5329 ,5330 ,5331 ,5332 ,5333 ,5334 ,5374 ,5375 ,5405 ,5488 ,5489 ,5490 ,5491 ,5492 ,5493 ,5531 ,5532 ,5533 ,5534 ,5561 ,5563 ,5565 ,5613 ,5614 ,5615 ,5616 ,5617 ,5655 ,5656 ,5681 ,5685 ,5690 ,5717 ,5718 ,5719 ,5720 ,5721 ,5758 ,5759};
vector<int> flayers={25,25,25,23,23,12,13,14,14,24,26};
vector<int> frows  ={56,57,60,163,175,103,94,91,92,75,68};
int n_fadc = frows.size();
FADCManager FMan;
void FADC(){
	int runnum=5814;
	TString dir = base_dir+"MayRun/rootfiles/FADC/";
	TString filename = Form("run0%d_TPCWaveform_.root",runnum);
	FMan.LoadFile(dir+filename);
	cout<<"FADC(int dum)"<<endl;
	cout<<"FADCBase()"<<endl;
}
void FADC(int dum){
	for(int i=0;i<n_fadc;++i){
		flayers.push_back(24);
	}
		flayers.push_back(16);
		frows.push_back(29);
	
	TCanvas* c1 = new TCanvas("c2","c2",1500,1200);
//	c1->Divide(2,1);
	TPad* pad1 = new TPad("pad1","pad1",0.,0.,1,1);
	TPad* pad2 = new TPad("pad2","pad2",0.,0.,1,1);
	pad2->Divide(6,3);
//	int ent2 = gTPCManager.GetEntries();

	int ent = FMan.GetEntries();
#if 1	//WaveformDisplay	
	int start = 00;
	for(int i=start;i<start+ent;++i){
		FMan.SetEvent(i);
		cout<<i<<endl;
		cout<<"Loading WF"<<endl;
		FMan.LoadWaveform();
		cout<<"Loading WF"<<endl;
		int wn = FMan.GetWaveNum();
		bool flag = false;
		double par[10];
		TH1D* wh[31];TH1D* rh[31];
		int cnt = 0;
		pad1->cd();
		auto base = FMan.GetBaseline();
		//base->Draw();
		for(int j=0;j<18;++j){
			wh[j] = FMan.GetWaveformHist(flayers[j],frows[j]);
			rh[j] = FMan.GetRawWaveformHist(flayers[j],frows[j]);
			if(cnt>11)continue;
			if(wh[j]!=nullptr){
				cnt++;
				pad2->cd(cnt);
				FMan.DoFit(flayers[j],frows[j],par);	
				wh[j]->Draw();
				rh[j]->SetLineColor(kGreen);
				rh[j]->Draw("same");
				rh[j]->SetTitle(Form("LR=(%d,%d),ev=%d",flayers[j],frows[j],i));
			}
		}
		for(int j=cnt;j<12;++j)	
		{
			pad2->cd(j+1);
			rh[j]= new TH1D(Form("rdum%d_%d",j,i),Form("rdum%d_%d",j,i),10,0,10);
			rh[j]->Draw("col");
		}
		c1->cd();
		pad2->Draw();

		pad1->Modified();
		pad1->Update();
		pad2->Modified();
		pad2->Update();
		pad2->cd(18);
		base->Draw();
		c1->Modified();
		c1->Update();
		gSystem->ProcessEvents();
		cin.ignore();
		FMan.Clear();
	}
#endif
#if 0//Draw BC
	TString tpcdir = base_dir+"MayRun/rootfiles/Defocus/";
	TString tpcfilename = Form("run0%d_DstTPCBcOutOld.root",runnum);
	gTPCManager.LoadFile(tpcdir+tpcfilename);
	gTPCManager.LoadBcOut();
	gTPCManager.LoadClusterChain();
	gTPCManager.InitializeHistograms();

	TGraph* BcGr ;
	TGraph* ClGr ;
	TF1* Bcf = new TF1("Bcf","pol1",-300,300);	
	TF1* Clf = new TF1("Clf","pol1",-300,300);	
	Clf->SetLineColor(kGreen);
	if(ent!=ent2) cout<<"Warning! event number does not match!"<<endl;
	auto* h = gTPCManager.GetPadHistogram();

	for(int i=599;i<ent;++i){	
		gTPCManager.ClearHistogram();
		//TPCDisplay	
		gTPCManager.SetEvent(i);
		FMan.SetEvent(i);
		int nh = gTPCManager.GetNhits(0);
		int ncl = gTPCManager.GetNhits(1);
		if(gTPCManager.GetBCnt()!=1||ncl>35){
			continue;
		}
		bool go = true;
		for(int irow = 60;irow<70;++irow){
			if(FMan.CheckLayerRowHit(31,irow)) go = false;
		}
		if(go) continue;
		cout<<"Go!"<<endl;
		Track Bc = gTPCManager.GetTrack();
		double bcx=Bc.GetPosition(-K18HS).X();
		cout<<"FillingHist"<<endl;
		for(int j=0;j<nh;++j){
			gTPCManager.FillHist(j);
		}
		double bcpx[100]={0},bcpz[100]={0};
		double clpx[100]={0},clpz[100]={0};
		cout<<ncl<<endl;
		int hnp = 0;
		for(int j=0;j<ncl;++j){
			TVector3 ClPos = gTPCManager.GetClusterPosition(j);
			double clx=ClPos.X(),clz=ClPos.Z();
			double bcx=Bc.GetPosition(clz).X();
			int l,r;
			int bcpad = tpc::findPadID(clz,bcx);
			GetLayerRow(bcpad,l,r);
			if(l>28) continue;
			bcpx[hnp]=bcx;bcpz[hnp]=clz;
			clpx[hnp]=clx;clpz[hnp]=clz;
			hnp++;
			//			cout<<Form("BcPos: (%f,%f)",bcpz[j],bcpx[j])<<endl;
		}
		BcGr = new TGraph(hnp,bcpz,bcpx);
		ClGr = new TGraph(hnp,clpz,clpx);
		pad1->cd();	
		cout<<"Hist"<<endl;
		h->SetTitle(Form("Run05754,Event%d",i));
		//	h->Draw("colz");
		h->Draw("");
		h->GetXaxis()->SetRangeUser(-250,250);
		h->GetYaxis()->SetRangeUser(-50,120);
		//		gTPCManager.SetPadContent(31,59,2);
		//		gTPCManager.SetPadContent(31,70,2);
		BcGr->Draw("Psame");
		BcGr->SetMarkerStyle(1);
		BcGr->Fit("Bcf");
		ClGr->Draw("Psame");
		ClGr->SetMarkerStyle(20);
		ClGr->SetMarkerColor(kBlue);
		ClGr->Fit("Clf");
		c1->cd(1);
		pad1->Draw();

		cout<<"BCTrack"<<endl;
	}
#endif
}


void FADCBase(){

	int ent = FMan.GetEntries();
	TH1D* h = new TH1D("h","h",1000,0,100);
	for(int i=0;i<ent;++i){
		FMan.SetEvent(i);
		FMan.LoadWaveform();
		auto base = FMan.GetBaseline();
		double rms = base->GetRMS();	
		h->Fill(rms);

	}
	h->Draw();


}
