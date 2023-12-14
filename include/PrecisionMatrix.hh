#ifndef PrecisionMatrix_h
#define PrecisionMatrix
#include <gmpxx.h>
mpf_set_default_prec(256);
void ReduceDim(double* Min, double* Mout, int dim, int row, int col){
	int cnt=0;
	for(int i=0;i<dim;i++){
		if(i%dim==col){
			continue;
		}
		else{
			for(int j=0;j<dim;j++){
				if(j%dim==row){
					continue;
				}
				else{
					Mout[cnt]=Min[dim*i+j];
					cnt++;
				
				}
			}
		}
	}
}
class Mat4D{
	private:
		double M_[16];
		int dim=4;
	public:
		Mat4D(){};
		Mat4D(double* M){
			for(int i=0;i<16;i++){
				M_[i]=M[i];
			}
		}
		void List(){
			for(int i=0;i<dim*dim;i++){
				cout<<M_[i]<<endl;
			}
		}
		double GetElement(int i, int j){
			if(dim*i+j<dim*dim){
				return M_[dim*i+j];
			}
			else{
				return -99999;
			}
		}
		void GetElements(double* M){
			for(int i=0;i<16;i++){
				M[i]=M_[i];
		}
		}
		void Sum(Mat4D A){
			for(int i=0;i<16;i++){
				M_[i]=M_[i]+A.GetElement((dim*i)%dim,i%dim);
			}
		}
		void Sub(Mat4D A){
			for(int i=0;i<16;i++){
				M_[i]=M_[i]-A.GetElement((dim*i)%dim,i%dim);
			}
		}
		
		Mat3D Minor(int col, int row){
			double Mred[9];
			ReduceDim(M_,Mred,dim,col,row);
			return Mat3D(Mred);
		}
		double Determinant(){
			double val=0;
			for(int i=0;i<dim;i++){
				val+=M_[i]*(Minor(i,0).Determinant())*PM(i);
			}
			return val;
		}
		Mat4D Transpose(){
			double MOut[dim*dim];
			for(int i=0;i<dim;i++){
				for(int j=0;j<dim;j++){
					MOut[dim*i+j]=M_[dim*j+i];
				}
			}
			Mat4D MT(MOut);
			return MT;
		}
		Mat4D Invert(){
			double MOut[dim*dim];
			for(int i=0;i<dim;i++){
				for(int j=0;j<dim;j++){
					MOut[j*dim+i]=PM(i+j)*(Minor(j,i).Determinant()) /Determinant();
				}
			}
			Mat4D InvMat(MOut);
			return InvMat;
		}
		Mat4D Product( Mat4D B){
			double MOut[dim*dim];
				for(int i=0;i<dim;i++){
					for(int j=0;j<dim;j++){
							MOut[dim*i+j]=0;
						for(int k=0;k<dim;k++){
							MOut[dim*i+j]+=(GetElement(i,k))*(B.GetElement(k,j));
						}
					}
				}
			Mat4D C(MOut);
			return C;
		}
};

#endif
