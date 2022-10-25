#include "Math.hh"
#ifndef Matrix_h
#define Matrix_h
bool IsZero(double a){
	if(abs(a)<1e-6) return true;
	else return false;
}
class Matrix{
	private:
		static const int max_dim = 5;
		static const int max_row = max_dim;
		static const int max_col = max_dim;
		 double element[max_row][max_col];
		int nrow,ncol;
		 double constants[max_row]={0};
		 double variable[max_row]={0};
	public:
		Matrix(int row,int col){
			nrow=row;ncol=col;
		}
		Matrix(int row,int col,double el[][max_dim]){
			nrow=row;ncol=col;
			Initialize(el);
		}
		int GetNRow(){return nrow;}
		int GetNCol(){return ncol;}
		void Initialize( double el[][max_dim]);
		 double GetElement(int row,int col){
			return element[row][col];
		}

		void Show();

		void AssignRow(int row, double* cols){
			for(int i=0;i<max_col;++i) cols[i]=element[row][i];
		}
		void AssignCol(int col, double* rows){
			for(int i=0;i<max_row;++i) rows[i]=element[i][col];
		}



		Matrix operator=(Matrix B);
		Matrix operator+(Matrix B);
		Matrix operator-(Matrix B);
		Matrix operator*(Matrix B);


		void SetConstants( double* sol){
			for(int i=0;i<max_row;++i) constants[i]=sol[i];
		}
		void GetSolution( double* var){
			for(int i=0;i<max_row;++i) var[i] = variable[i];
		}
		void MultiplyRow(int row, double val);
		void SubtractRow(int row, double* val);
		void RowSubtraction(int src,int dst, double multi);



		void GaussianElimination( double* var);

};
void Matrix::Show(){
	cout<<Form("Matrix of (%d,%d)",nrow,ncol)<<endl;
	for(int row=0;row<nrow;++row){
		cout<<"[	";
		for(int col=0;col<ncol;++col){
			cout<<Form(" %f , ",GetElement(row,col));
		}
		cout<<" | "<<constants[row]<<" ]"<<endl;
	}
}
void Matrix::Initialize( double el[][max_dim]){
	for(int row=0;row<nrow;++row){
		for(int col=0;col<ncol;++col){
			element[row][col]=el[row][col];
		}
	}
}
void Matrix::MultiplyRow(int row, double val){
	for(int col=0;col<ncol;++col){
		element[row][col]*=val;
	}
	constants[row]*=val;
}
void Matrix::SubtractRow(int row, double* val){
	for(int col=0;col<ncol;++col){
		element[row][col]-=val[col];
	}
	constants[row]-=val[ncol];
}
void Matrix::RowSubtraction(int src,int dst, double multi=1){
	for(int col=0;col<ncol;++col){
		element[dst][col]-=multi*element[src][col];
	}
	constants[dst]-=multi*constants[src];
}
Matrix Matrix::operator=(Matrix B){
	 double temp[max_dim][max_dim];
	for(int row=0;row<nrow;++row){
		for(int col=0;col<ncol;++col){
			temp[row][col]=B.GetElement(row,col);
		}
	}
	Matrix A(nrow,ncol);
	A.Initialize(temp);
	return A;
}
Matrix Matrix::operator+(Matrix B){
	 double temp[max_dim][max_dim];
	for(int row=0;row<nrow;++row){
		for(int col=0;col<ncol;++col){
			temp[row][col]=GetElement(row,col)+B.GetElement(row,col);
		}
	}
	Matrix A(nrow,ncol);
	A.Initialize(temp);
	return A;
}
Matrix Matrix::operator-(Matrix B){
	 double temp[max_dim][max_dim];
	for(int row=0;row<nrow;++row){
		for(int col=0;col<ncol;++col){
			temp[row][col]=GetElement(row,col)-B.GetElement(row,col);
		}
	}
	Matrix A(nrow,ncol);
	A.Initialize(temp);
	return A;
}
Matrix Matrix::operator*(Matrix B){
	 double temp[max_dim][max_dim];
	int nrow_ = GetNRow();
	int ncol_ = B.GetNCol();

	for(int row=0;row<nrow_;++row){
		for(int col=0;col<ncol_;++col){
			temp[row][col]=0;
			for(int k=0;k<ncol;++k){
				temp[row][col]+=GetElement(row,k)*B.GetElement(k,col);
			}
		}
	}
	Matrix A(nrow_,ncol_);
	A.Initialize(temp);
	return A;
}


void Matrix::GaussianElimination( double* var){
	cout<<"Before: "<<endl;
	Show();
	for(int row=0;row<nrow;++row){
		 double head = GetElement(row,row);
		MultiplyRow(row,1/head);
		for(int i = 0;i<nrow;++i){
			 double subhead  = GetElement(i,row);
			if(i!=row) RowSubtraction(row,i,subhead);
		}
	}
	for(int row = 0;row<nrow;++row){
		variable[row]=constants[row];
	}
	GetSolution(var);
	cout<<"After: "<<endl;
	Show();
}
#endif
