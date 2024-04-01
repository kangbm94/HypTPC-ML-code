
double mXi = 1.321;
double mL = 1.115;
double mP = 0.938;
double mK = 0.493;
double mPi = 0.139;
double ResP = 0.1;
double ResPi1 = 0.1;
double ResPi2 = 0.1;
double ResTh = 0.018; 
double ResPh = 0.038; 
double ResThP = 0.026; 
double ResPhP = 0.029; 
double ResThV = 0.05; 
double ResPhV = 0.05; 

double ResPLd = 0.03;
double ResThLd = 0.01;
double ResPhLd = 0.01;

double PLd,ThLd,PhLd; 
double PP,ThP,PhP; 
double PPi1,ThPi1,PhPi1; 
double PPi2,ThPi2,PhPi2; 
double PXi,ThXi,PhXi; 
double PLdMeas,ThLdMeas,PhLdMeas; 
double PLdCorMeas,ThLdCorMeas,PhLdCorMeas; 
double PXiMeas,ThXiMeas,PhXiMeas; 
double PPMeas,ThPMeas,PhPMeas; 
double PPi1Meas,ThPi1Meas,PhPi1Meas; 
double PPi2Meas,ThPi2Meas,PhPi2Meas; 
double PLdCor,ThLdCor,PhLdCor; 
double PLdFit,ThLdFit,PhLdFit; 
double PLdCorCor,ThLdCorCor,PhLdCorCor; 
double PXiCor,ThXiCor,PhXiCor; 
double PPCor,ThPCor,PhPCor; 
double PPi1Cor,ThPi1Cor,PhPi1Cor; 
double PPi2Cor,ThPi2Cor,PhPi2Cor; 
double InvMLd,InvMXi,InvMLdCor,InvMXiCor;
double Chi2,Chi2Xi;
double PullPLd,PullThLd,PullPhLd,PullPP,PullThP,PullPhP,PullPPi1,PullThPi1,PullPhPi1;
double PullPLdCor,PullThLdCor,PullPhLdCor,PullPPi2,PullThPi2,PullPhPi2;
int NStep,BestStep;
int Step[200];
double StepChi2[200];
double StepMassDiff[200];
int iev;
void SetBranches(TTree* tree){
	tree->Branch("iev",&iev);
	tree->Branch("PLd",&PLd);
	tree->Branch("ThLd",&ThLd);
	tree->Branch("PhLd",&PhLd);
	tree->Branch("PXi",&PXi);
	tree->Branch("ThXi",&ThXi);
	tree->Branch("PhXi",&PhXi);
	tree->Branch("PP",&PP);
	tree->Branch("ThP",&ThP);
	tree->Branch("PhP",&PhP);
	tree->Branch("PPi1",&PPi1);
	tree->Branch("ThPi1",&ThPi1);
	tree->Branch("PhPi1",&PhPi1);
	tree->Branch("PPi2",&PPi2);
	tree->Branch("ThPi2",&ThPi2);
	tree->Branch("PhPi2",&PhPi2);
	
	tree->Branch("InvMLd",&InvMLd);	
	tree->Branch("InvMXi",&InvMXi);	
	tree->Branch("PLdMeas",&PLdMeas);
	tree->Branch("ThLdMeas",&ThLdMeas);
	tree->Branch("PhLdMeas",&PhLdMeas);
	tree->Branch("PPMeas",&PPMeas);
	tree->Branch("ThPMeas",&ThPMeas);
	tree->Branch("PhPMeas",&PhPMeas);
	tree->Branch("PPi1Meas",&PPi1Meas);
	tree->Branch("ThPi1Meas",&ThPi1Meas);
	tree->Branch("PhPi1Meas",&PhPi1Meas);
	tree->Branch("PPi2Meas",&PPi2Meas);
	tree->Branch("ThPi2Meas",&ThPi2Meas);
	tree->Branch("PhPi2Meas",&PhPi2Meas);
	
	tree->Branch("InvMLdCor",&InvMLdCor);	
	tree->Branch("InvMXiCor",&InvMXiCor);	
	tree->Branch("PLdCor",&PLdCor);
	tree->Branch("ThLdCor",&ThLdCor);
	tree->Branch("PhLdCor",&PhLdCor);
	tree->Branch("PLdCorCor",&PLdCorCor);
	tree->Branch("ThLdCorCor",&ThLdCorCor);
	tree->Branch("PhLdCorCor",&PhLdCorCor);
	tree->Branch("PXiCor",&PXiCor);
	tree->Branch("ThXiCor",&ThXiCor);
	tree->Branch("PhXiCor",&PhXiCor);
	tree->Branch("PPCor",&PPCor);
	tree->Branch("ThPCor",&ThPCor);
	tree->Branch("PhPCor",&PhPCor);
	tree->Branch("PPi1Cor",&PPi1Cor);
	tree->Branch("ThPi1Cor",&ThPi1Cor);
	tree->Branch("PhPi1Cor",&PhPi1Cor);
	tree->Branch("PPi2Cor",&PPi2Cor);
	tree->Branch("ThPi2Cor",&ThPi2Cor);
	tree->Branch("PhPi2Cor",&PhPi2Cor);

	tree->Branch("NStep",&NStep);
	tree->Branch("BestStep",&BestStep);
	tree->Branch("Chi2",&Chi2);
	tree->Branch("Step",Step,"Step[NStep]/I");
	tree->Branch("StepChi2",StepChi2,"StepChi2[NStep]/D");
	tree->Branch("StepMassDiff",StepMassDiff,"StepMassDiff[NStep]/D");

}
	TH1D* HistLd = new TH1D("InvMLd","InvMLd",100,1.08,1.14);
	TH1D* HistXi = new TH1D("InvMXi","InvMXi",100,1.25,1.4);
	TH1D* HistLdFit = new TH1D("InvMLdFit","InvMLdFit",1000,1.08,1.14);
	TH1D* HistXiFit = new TH1D("InvMXiFit","InvMXiFit",100,1.25,1.4);
	TH1D* HistLdP = new TH1D("HistLdP","HistLdP",100,-0.3,0.3);
	TH1D* HistLdTh = new TH1D("HistLdTh","HistLdTh",300,-0.3,0.3);
	TH1D* HistLdPh = new TH1D("HistLdPh","HistLdPh",300,-0.3,0.3);
	TH1D* HistPP = new TH1D("HistPP","HistPP",100,-0.3,0.3);
	TH1D* HistPTh = new TH1D("HistPTh","HistPTh",300,-0.3,0.3);
	TH1D* HistPPh = new TH1D("HistPPh","HistPPh",300,-0.3,0.3);
	TH1D* HistPi1P = new TH1D("HistPi1P","HistPi1P",100,-0.05,0.05);
	TH1D* HistPi1Th = new TH1D("HistPi1Th","HistPi1Th",300,-0.3,0.3);
	TH1D* HistPi1Ph = new TH1D("HistPi1Ph","HistPi1Ph",300,-0.3,0.3);
	TH1D* HistXiP = new TH1D("HistXiP","HistXiP",100,-0.3,0.3);
	TH1D* HistXiTh = new TH1D("HistXiTh","HistXiTh",100,-0.3,0.3);
	TH1D* HistXiPh = new TH1D("HistXiPh","HistXiPh",100,-0.3,0.3);
	TH1D* HistLdCorP = new TH1D("HistLdCorP","HistLdCorP",100,-0.3,0.3);
	TH1D* HistLdCorTh = new TH1D("HistLdCorTh","HistLdCorTh",100,-0.3,0.3);
	TH1D* HistLdCorPh = new TH1D("HistLdCorPh","HistLdCorPh",100,-0.3,0.3);
	
	TH1D* HistPi2P = new TH1D("HistPi2P","HistPi2P",100,-0.1,0.1);
	TH1D* HistPi2Th = new TH1D("HistPi2Th","HistPi2Th",100,-0.3,0.3);
	TH1D* HistPi2Ph = new TH1D("HistPi2Ph","HistPi2Ph",100,-0.3,0.3);
	
	TH1D* HistDLdP = new TH1D("HistDLdP","HistDLdP",100,-0.3,0.3);
	TH1D* HistDLdTh = new TH1D("HistDLdTh","HistDLdTh",300,-0.3,0.3);
	TH1D* HistDLdPh = new TH1D("HistDLdPh","HistDLdPh",300,-0.3,0.3);
	TH1D* HistDPP = new TH1D("HistDPP","HistDPP",100,-0.3,0.3);
	TH1D* HistDPTh = new TH1D("HistDPTh","HistDPTh",300,-0.3,0.3);
	TH1D* HistDPPh = new TH1D("HistDPPh","HistDPPh",300,-0.3,0.3);
	TH1D* HistDPi1P = new TH1D("HistDPi1P","HistDPi1P",100,-0.05,0.05);
	TH1D* HistDPi1Th = new TH1D("HistDPi1Th","HistDPi1Th",300,-0.3,0.3);
	TH1D* HistDPi1Ph = new TH1D("HistDPi1Ph","HistDPi1Ph",300,-0.3,0.3);

	TH1D* HistDXiP = new TH1D("HistDXiP","HistDXiP",100,-0.3,0.3);
	TH1D* HistDXiTh = new TH1D("HistDXiTh","HistDXiTh",100,-0.3,0.3);
	TH1D* HistDXiPh = new TH1D("HistDXiPh","HistDXiPh",100,-0.3,0.3);
	TH1D* HistDLdCorP = new TH1D("HistDLdCorP","HistDLdCorP",100,-0.3,0.3);
	TH1D* HistDLdCorTh = new TH1D("HistDLdCorTh","HistDLdCorTh",100,-0.3,0.3);
	TH1D* HistDLdCorPh = new TH1D("HistDLdCorPh","HistDLdCorPh",100,-0.3,0.3);
	TH1D* HistDPi2P = new TH1D("HistDPi2P","HistDPi2P",100,-0.1,0.1);
	TH1D* HistDPi2Th = new TH1D("HistDPi2Th","HistDPi2Th",100,-0.3,0.3);
	TH1D* HistDPi2Ph = new TH1D("HistDPi2Ph","HistDPi2Ph",100,-0.3,0.3);



	TH1D* HistDifLdP = new TH1D("HistDifLdP","HistDifLdP",100,-0.1,0.1);
	TH1D* HistDifLdTh = new TH1D("HistDifLdTh","HistDifLdTh",100,-0.3,0.3);
	TH1D* HistDifLdPh = new TH1D("HistDifLdPh","HistDifLdPh",100,-0.3,0.3);
	TH1D* HistDifPP = new TH1D("HistDifPP","HistDifPP",100,-0.3,0.3);
	TH1D* HistDifPTh = new TH1D("HistDifPTh","HistDifPTh",100,-0.3,0.3);
	TH1D* HistDifPPh = new TH1D("HistDifPPh","HistDifPPh",100,-0.3,0.3);
	TH1D* HistDifPi1P = new TH1D("HistDifPi1P","HistDifPi1P",100,-0.1,0.1);
	TH1D* HistDifPi1Th = new TH1D("HistDifPi1Th","HistDifPi1Th",100,-0.3,0.3);
	TH1D* HistDifPi1Ph = new TH1D("HistDifPi1Ph","HistDifPi1Ph",100,-0.3,0.3);

	TH1D* HistChi2 = new TH1D("Chi2","Chi2",500,0,50);
	TH2D* HistChi22D = new TH2D("Chi2:InvM","Chi2:InvM",100,0,20,100,1.1,1.14);
	TH1D* HistChi2Xi = new TH1D("Chi2Xi","Chi2Xi",500,0,50);
	TH2D* HistChi2Xi2D = new TH2D("Chi2Xi:InvM","Chi2Xi:InvM",100,0,20,100,1.25,1.40);
	TH1D* HistProb = new TH1D("Prob","Prob",100,0,1.1);
	TH2D* HistProb2D = new TH2D("Prob:InvM","Prob:InvM",100,0,1.1,100,1.1,1.14);
	TH1D* HistProbXi = new TH1D("ProbXi","ProbXi",100,0,1.1);
	TH2D* HistProbXi2D = new TH2D("ProbXi:InvM","ProbXi:InvM",100,0,1.1,100,1.25,1.40);
	TH1D* HistMM = new TH1D("MM","MM",100,-0.3,0.3);
	TH1D* HistMMCor = new TH1D("MMCor","MMCor",100,-0.3,0.3);

	TH1D* HistPullPLd = new TH1D("HistPullPLd","HistPullPLd",1000,-5,5);
	TH1D* HistPullThLd = new TH1D("HistPullThLd","HistPullThLd",1000,-5,5);
	TH1D* HistPullPhLd = new TH1D("HistPullPhLd","HistPullPhLd",1000,-5,5);
	TH1D* HistPullPP = new TH1D("HistPullPP","HistPullPP",100,-5,5);
	TH1D* HistPullThP = new TH1D("HistPullThP","HistPullThP",100,-5,5);
	TH1D* HistPullPhP = new TH1D("HistPullPhP","HistPullPhP",100,-5,5);
	TH1D* HistPullPPi1 = new TH1D("HistPullPPi1","HistPullPPi1",100,-5,5);
	TH1D* HistPullThPi1 = new TH1D("HistPullThPi1","HistPullThPi1",100,-5,5);
	TH1D* HistPullPhPi1 = new TH1D("HistPullPhPi1","HistPullPhPi1",100,-5,5);
	

	TH1D* HistPullPLdCor = new TH1D("HistPullPLdCor","HistPullPLdCor",100,-5,5);
	TH1D* HistPullThLdCor = new TH1D("HistPullThLdCor","HistPullThLdCor",100,-5,5);
	TH1D* HistPullPhLdCor = new TH1D("HistPullPhLdCor","HistPullPhLdCor",100,-5,5);
	TH1D* HistPullPPi2 = new TH1D("HistPullPPi2","HistPullPPi2",100,-5,5);
	TH1D* HistPullThPi2 = new TH1D("HistPullThPi2","HistPullThPi2",100,-5,5);
	TH1D* HistPullPhPi2 = new TH1D("HistPullPhPi2","HistPullPhPi2",100,-5,5);

	TH1D* HistPx = new TH1D("HistPx","HistPx",1000,-0.2,0.2);
	TH1D* HistPy = new TH1D("HistPy","HistPy",1000,-0.2,0.2);
	TH1D* HistPz = new TH1D("HistPz","HistPz",1000,-0.1,0.1);
	TH1D* HistE = new TH1D("HistE","HistE",1000,-0.01,0.01);
	
	TH1D* HistDPx = new TH1D("HistDPx","HistDPx",1000,-0.2,0.2);
	TH1D* HistDPy = new TH1D("HistDPy","HistDPy",1000,-0.2,0.2);
	TH1D* HistDPz = new TH1D("HistDPz","HistDPz",1000,-0.1,0.1);
	TH1D* HistDE = new TH1D("HistDE","HistDE",1000,-0.01,0.01);
