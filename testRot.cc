void testRot(){
	TLorentzVector P(2,8,1,10);
	TLorentzVector Q(-1,-1,-3,5);

	auto Zaxis =(P + Q).Vect();
	auto vP = P.Vect();
	auto vQ = Q.Vect();
	auto Yaxis = vP.Cross(vQ);
	double Th_F = Zaxis.Theta();
	double Ph_F = Zaxis.Phi();
	cout<<Form("F = %g,%g",Th_F,Ph_F)<<endl;
	double RotZ[9] ={
		cos(-Ph_F),	-sin(-Ph_F),	0,	
		sin(-Ph_F),	cos(-Ph_F),		0,
		0,					0,						1
	};
	TMatrixD R_Z(3,3,RotZ);
	double Th_FY = Zaxis.Theta();
	double Ph_FY = Zaxis.Phi();
	double hp = 0.5 * acos(-1);
	double RotY[9] ={
		cos(Th_F),		0,				-sin(Th_F),
		0,						-1,				0,
		sin(Th_F),	0,				cos(Th_F)
	};
	TMatrixD R_Y(3,3,RotY);
	TMatrixD R_F = R_Y * R_Z;
	cout<<Form("Y(%g,%g,%g)",Yaxis.x(),Yaxis.y(),Yaxis.z())<<endl;
	cout<<Form("Z(%g,%g,%g)",Zaxis.x(),Zaxis.y(),Zaxis.z())<<endl;
	Yaxis = R_F * Yaxis;
	Zaxis = R_F * Zaxis;
	double Th_Z = Zaxis.Theta();
	double Ph_Z = Zaxis.Phi();
	double Th_Y = Yaxis.Theta();
	double Ph_Y = Yaxis.Phi();
	cout<<Form("Y(%g,%g,%g)",Yaxis.x(),Yaxis.y(),Yaxis.z())<<endl;
	cout<<Form("Z(%g,%g,%g)",Zaxis.x(),Zaxis.y(),Zaxis.z())<<endl;
	double RotX[9] ={
		sin(-Ph_Y),	cos(-Ph_Y),	0,	
		cos(-Ph_Y),	-sin(-Ph_Y),		0,
		0,					0,						1
	};
	TMatrixD RX(3,3,RotX);
	Yaxis = RX * Yaxis;
	Zaxis = RX * Zaxis;
	cout<<Form("Y(%g,%g,%g)",Yaxis.x(),Yaxis.y(),Yaxis.z())<<endl;
	cout<<Form("Z(%g,%g,%g)",Zaxis.x(),Zaxis.y(),Zaxis.z())<<endl;
	

}
