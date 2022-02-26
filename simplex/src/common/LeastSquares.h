//////////////////////////////////////////////////////////////////////////
// Least squares
// Copyright (c) (2018-), Bo Zhu, Xiangxin Kong, Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __LeastSquares_h__
#define __LeastSquares_h__
#include "Common.h"
#include "AuxFunc.h"
#include <iostream>

namespace LeastSquares
{
	template<int d,int p> class LS {};
	template<int d,int p> class WLS{};
	template<int d,int p> class MLS{};
	template<int d> class MLSPointSet{};

	//////////////////////////////////////////////////////////////////////////
	////LS

	////f(x)=c0+c1*x
	template<> class LS<1,1>
	{public:
		VectorX c;
	
		////[x0,f0,x1,f1,x2,f2,...]
		void Fit(const real* data,const int m)
		{
			MatrixX B(m,2);
			VectorX f(m);
			for(int i=0;i<m;i++){
				B.coeffRef(i,0)=(real)1;
				B.coeffRef(i,1)=data[i*2];
				f[i]=data[i*2+1];}
			c=B.colPivHouseholderQr().solve(f);
		}

		inline real operator() (const real& x) const {return c[0]+c[1]*x;}
		inline real operator() (const Vector1& x) const {return c[0]+c[1]*x[0];}
		inline Vector1 Grad(const Vector1& x) const {return Vector1(c[1]);}
	};

	////f(x)=c0+c1*x+c2*x^2
	template<> class LS<1,2>
	{public:
		VectorX c;
	
		////[x0,f0, x1,f1, x2,f2, ...]
		void Fit(const real* data,const int m)
		{
			MatrixX B(m,3);
			VectorX f(m);
			for(int i=0;i<m;i++){
				B.coeffRef(i,0)=(real)1;
				B.coeffRef(i,1)=data[i*2];
				B.coeffRef(i,2)=pow(data[i*2],2);
				f[i]=data[i*2+1];}
			c=B.colPivHouseholderQr().solve(f);
		}

		inline real operator() (const real& x) const {return c[0]+c[1]*x+c[2]*pow(x,2);}
		inline real operator() (const Vector1& x) const {return c[0]+c[1]*x[0]+c[2]*pow(x[0],2);}
		inline Vector1 Grad(const Vector1& x) const {return Vector1(c[1]+(real)2*c[2]*x[0]);}
	};

	////f(x)=c0+c1*x+c2*x^2+c3*x^3
	template<> class LS<1,3>
	{public:
		VectorX c;
	
		////[x0,f0, x1,f1, x2,f2, ...]
		void Fit(const real* data,const int m)
		{
			MatrixX B(m,4);
			VectorX f(m);
			for(int i=0;i<m;i++){
				B.coeffRef(i,0)=(real)1;
				B.coeffRef(i,1)=data[i*2];
				B.coeffRef(i,2)=pow(data[i*2],2);
				B.coeffRef(i,3)=pow(data[i*2],3);
				f[i]=data[i*2+1];}
			c=B.colPivHouseholderQr().solve(f);
		}

		inline real operator() (const real& x) const {return c[0]+c[1]*x+c[2]*pow(x,2)+c[3]*pow(x,3);}
		inline real operator() (const Vector1& x) const {return c[0]+c[1]*x[0]+c[2]*pow(x[0],2)+c[3]*pow(x[0],3);}
		inline Vector1 Grad(const Vector1& x) const {return Vector1(c[1]+(real)2*c[2]*x[0]+(real)3*c[3]*pow(x[0],2));}
	};


	////f(x,y)=c0+c1*x+c2*y
	template<> class LS<2,1>
	{public:
		VectorX c;
	
		////[x0,y0,f0, x1,y1,f1, x2,y2,f2, ...]
		void Fit(const real* data,const int m)
		{
			MatrixX B(m,3);
			VectorX f(m);
			for(int i=0;i<m;i++){
				B.coeffRef(i,0)=(real)1;
				B.coeffRef(i,1)=data[i*3];
				B.coeffRef(i,2)=data[i*3+1];
				f[i]=data[i*3+2];}
			c=B.colPivHouseholderQr().solve(f);
		}

		inline real operator() (const real& x,const real& y) const {return c[0]+c[1]*x+c[2]*y;}
		inline real operator() (const Vector2& x) const {return c[0]+c[1]*x[0]+c[2]*x[1];}
		inline Vector2 Grad(const Vector2& x) const {return Vector2(c[1],c[2]);}
	};

	////f(x,y)=c0+c1*x+c2*y+c3*x^2+c4*y^2+c5*xy
	template<> class LS<2, 2>
	{
	public:
		Vector<double, 6> c;

		////[x0,y0,f0, x1,y1,f1, x2,y2,f2, ...]
		void Fit(const real* data, const int m)
		{
			if (m == 0) { AuxFunc::Crash_With_Info("LS::Fit passed m=0"); }
			MatrixX B(m, 6);
			VectorX f(m);
			for (int i = 0; i < m; i++) {
				B.coeffRef(i, 0) = (real)1;
				B.coeffRef(i, 1) = data[i * 3];
				B.coeffRef(i, 2) = data[i * 3 + 1];
				B.coeffRef(i, 3) = pow(data[i * 3], 2);
				B.coeffRef(i, 4) = pow(data[i * 3 + 1], 2);
				B.coeffRef(i, 5) = data[i * 3] * data[i * 3 + 1];
				f[i] = data[i * 3 + 2];
			}
			c = B.colPivHouseholderQr().solve(f);
		}

		inline real operator() (const real& x, const real& y) const
		{
			return c[0] + c[1] * x + c[2] * y + c[3] * pow(x, 2) + c[4] * pow(y, 2) + c[5] * x * y;
		}

		inline real operator() (const Vector2& x) const
		{
			return c[0] + c[1] * x[0] + c[2] * x[1] + c[3] * pow(x[0], 2) + c[4] * pow(x[1], 2) + c[5] * x[0] * x[1];
		}

		inline Vector2 Grad(const Vector2& x) const
		{
			return Vector2(c[1] + (real)2 * c[3] * x[0] + c[5] * x[1], c[2] + (real)2 * c[4] * x[1] + c[5] * x[0]);
		}
	};

	//////////////////////////////////////////////////////////////////////////
	////WLS

	////f(x,y)=c0+c1*x+c2*y
	template<> class WLS<2,1>
	{public:
		VectorX c;
	
		////[x0,y0,f0,t0, x1,y1,f1,t1, x2,y2,f2,t2, ...]
		void Fit(const real* data,const int m)
		{
			MatrixX B(m,3);
			VectorX f(m);
			for(int i=0;i<m;i++){
				real w=sqrt(data[i*4+3]);
				B.coeffRef(i,0)=w*(real)1;
				B.coeffRef(i,1)=w*data[i*4];
				B.coeffRef(i,2)=w*data[i*4+1];
				f[i]=w*data[i*4+2];}
			c=B.colPivHouseholderQr().solve(f);
		}

		inline real operator() (const real& x,const real& y) const {return c[0]+c[1]*x+c[2]*y;}
		inline real operator() (const Vector2& x) const {return c[0]+c[1]*x[0]+c[2]*x[1];}
		inline Vector2 Grad(const Vector2& x) const {return Vector2(c[1],c[2]);}
	};

	//////////////////////////////////////////////////////////////////////////
	////MLS
	
	////f(x)=c0+c1*x+c2*x^2
	template<> class MLS<1, 2>
	{
	public:
		VectorX c;
		void Fit(const real* data, const int m, const real& x, const real& y) {
			MatrixX B(m, 3);
			VectorX f(m);
			for (int i = 0; i < m; i++) {
				real w = sqrt(Weight(data[i * 2], x));
				B.coeffRef(i, 0) = w * (real)1;
				B.coeffRef(i, 1) = w * data[i * 2];
				B.coeffRef(i, 2) = w * pow(data[i * 2], 2);
				f[i] = w * data[i * 2 + 1];
			}
			c = B.colPivHouseholderQr().solve(f);
		}
		inline real operator() (const Vector1& x) const {
			return c[0] + c[1] * x[0] + c[2] * pow(x[0], 2);
		}
	protected:
		inline real Weight(const real& x0, const real& x1) const
		{
			real dis2 = pow(x0 - x1, 2);
			return exp(-dis2);
		}
	};

	////f(x,y)=c0+c1*x+c2*y+c3*x^2+c4*y^2+c5*xy
	template<> class MLS<2, 2>
	{
	public:
		VectorX c;

		////[x0,y0,f0, x1,y1,f1, x2,y2,f2, ...]
		void Fit(const real* data, const int m, const real& x, const real& y)
		{
			MatrixX B(m, 6);
			VectorX f(m);
			for (int i = 0; i < m; i++) {
				real w = sqrt(Weight(data[i * 3], data[i * 3 + 1], x, y));
				B.coeffRef(i, 0) = w * (real)1;
				B.coeffRef(i, 1) = w * data[i * 3];
				B.coeffRef(i, 2) = w * data[i * 3 + 1];
				B.coeffRef(i, 3) = w * pow(data[i * 3], 2);
				B.coeffRef(i, 4) = w * pow(data[i * 3 + 1], 2);
				B.coeffRef(i, 5) = w * data[i * 3] * data[i * 3 + 1];
				f[i] = w * data[i * 3 + 2];
			}
			if (m <= 6) {
				std::cout << "MLS::Fit error: only " << m << " equations for 6 parameters. Fitting that may cause severe artifact\n";
				std::cout << "(x,y)=" << x << " , " << y << "\n";
				std::cout << "data: "; for (int i = 0; i < m * 3; i++) { std::cout << data[i] << " "; }std::cout << "\n";
				std::cout << "B:\n" << B << "\n";
				std::cout << "f: " << f.transpose() << "\n";
				exit(1);
			}
			c = B.colPivHouseholderQr().solve(f);
		}

		inline real operator() (const real& x, const real& y) const
		{
			return c[0] + c[1] * x + c[2] * y + c[3] * pow(x, 2) + c[4] * pow(y, 2) + c[5] * x * y;
		}

		inline real operator() (const Vector2& x) const
		{
			return c[0] + c[1] * x[0] + c[2] * x[1] + c[3] * pow(x[0], 2) + c[4] * pow(x[1], 2) + c[5] * x[0] * x[1];
		}

		inline Vector2 Grad(const Vector2& x) const
		{
			return Vector2(c[1] + (real)2 * c[3] * x[0] + c[5] * x[1], c[2] + (real)2 * c[4] * x[1] + c[5] * x[0]);
		}

	protected:
		inline real Weight(const real& x0, const real& y0, const real& x1, const real& y1) const
		{
			real dis2 = pow(x0 - x1, 2) + pow(y0 - y1, 2);
			return exp(-dis2);
		}
	};




	////f(x)=c0+c1*x+c2*x^2
	template<> class MLSPointSet<1>
	{public:
		VectorX c;
		int point_n=0;
		int p_idx=-1;

		////1D: [x0,f0, x1,f1, x2,f2, ...]
		////p is one particle among the particles with their positions specified in data
		void Fit(const real* data,const int m,const int p)
		{
			point_n=m;
			p_idx=p;

			MatrixX B(m,3);
			VectorX f(m);
			for(int i=0;i<m;i++){
				real w=sqrt(Weight(i));
				B.coeffRef(i,0)=w*(real)1;
				B.coeffRef(i,1)=w*data[i*2];
				B.coeffRef(i,2)=w*pow(data[i*2],2);
				f[i]=w*data[i*2+1];}
			c=B.colPivHouseholderQr().solve(f);
		}

		void Fit(const real* data,const real* ws,const int m)
		{
			point_n=m;

			MatrixX B(m,3);
			VectorX f(m);
			for(int i=0;i<m;i++){
				real w=sqrt(ws[i]);
				B.coeffRef(i,0)=w*(real)1;
				B.coeffRef(i,1)=w*data[i*2];
				B.coeffRef(i,2)=w*pow(data[i*2],2);
				f[i]=w*data[i*2+1];}
			c=B.colPivHouseholderQr().solve(f);
		}

		inline real operator() (const real& x) const 
		{return c[0]+c[1]*x+c[2]*pow(x,2);}	
		inline Vector1 Grad(const Vector1& x) const {return Vector1(c[1]+(real)2*c[2]*x[0]);}

	protected:
		inline real Weight(const int p) const
		{
			if(p==p_idx)return (real)1;
			else return (real)1/(real)point_n;
		}
	};

	////f(x,y)=c0+c1*x+c2*y+c3*x^2+c4*y^2+c5*xy
	template<> class MLSPointSet<2>
	{public:
		VectorX c;
		int point_n=0;
		int p_idx=-1;

		////2D: [x0,y0,f0, x1,y1,f1, x2,y2,f2, ...]
		void Fit(const real* data,const int m,const int p)
		{
			point_n=m;
			p_idx=p;

			MatrixX B(m,6);
			VectorX f(m);
			for(int i=0;i<m;i++){
				real w=sqrt(Weight(i));
				B.coeffRef(i,0)=w*(real)1;
				B.coeffRef(i,1)=w*data[i*3];
				B.coeffRef(i,2)=w*data[i*3+1];
				B.coeffRef(i,3)=w*pow(data[i*3],2);
				B.coeffRef(i,4)=w*pow(data[i*3+1],2);
				B.coeffRef(i,5)=w*data[i*3]*data[i*3+1];
				f[i]=w*data[i*3+2];}
			c=B.colPivHouseholderQr().solve(f);
		}

		////2D: [x0,y0,f0, x1,y1,f1, x2,y2,f2, ...]
		//// Fit with weights
		void Fit(const real* data,const real* ws,const int m)
		{
			point_n=m;

			MatrixX B(m,6);
			VectorX f(m);
			for(int i=0;i<m;i++){
				real w=sqrt(ws[i]);
				B.coeffRef(i,0)=w*(real)1;
				B.coeffRef(i,1)=w*data[i*3];
				B.coeffRef(i,2)=w*data[i*3+1];
				B.coeffRef(i,3)=w*pow(data[i*3],2);
				B.coeffRef(i,4)=w*pow(data[i*3+1],2);
				B.coeffRef(i,5)=w*data[i*3]*data[i*3+1];
				f[i]=w*data[i*3+2];}
			c=B.colPivHouseholderQr().solve(f);
		}

		inline real operator() (const real& x,const real& y) const 
		{return c[0]+c[1]*x+c[2]*y+c[3]*pow(x,2)+c[4]*pow(y,2)+c[5]*x*y;}	
		
		inline real operator() (const Vector2& x) const 
		{return c[0]+c[1]*x[0]+c[2]*x[1]+c[3]*pow(x[0],2)+c[4]*pow(x[1],2)+c[5]*x[0]*x[1];}	

		inline Vector2 Grad(const Vector2& x) const
		{return Vector2(c[1]+(real)2*c[3]*x[0]+c[5]*x[1],c[2]+(real)2*c[4]*x[1]+c[5]*x[0]);}

	protected:
		inline real Weight(const int p) const
		{
			if(p==p_idx)return (real)1;
			else return (real)1/(real)point_n;
		}
	};
};

#endif
