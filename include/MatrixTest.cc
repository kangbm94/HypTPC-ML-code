#include"Matrix.hh"
void MatrixTest(){
	Matrix Mat(3,3);
	double b[5][5] = { {1,2,3,4,0},
		{5,6,7,8,0},
		{9,10,12,13,0},
		{13,13,15,18,0}
	};
	double constants[3]={7,6,4};
	Mat.SetConstants(constants);
	double variable[3];
	Mat.Initialize(b);
	cout<<"Before Gaussian Elimination"<<endl;
//	Mat.Show();
	cout<<"After Gaussian Elimination"<<endl;
	Mat.GaussianElimination(variable);
//	Mat.Show();
	cout<<"Solution : "<<endl;
	for(int i=0;i<3;++i){
		cout<<Form("%f /12, ",12*variable[i])<<endl;
	}
}
