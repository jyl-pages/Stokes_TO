

#pragma warning( disable : 4786 ) 
//#define for if(0);else for

#include <vector>
#include <stack>
#include <string>
#include <math.h>
#include <assert.h>

#include "triangulation.h"

using namespace std;

const double MIN_TRI_AREA = 1.0e-10;

const unsigned int nNoEd  = 2;

const unsigned int nNoTri = 3;
const unsigned int nEdTri = 3;
const unsigned int noelTriEdge[nEdTri][nNoEd] = {
	{ 1, 2 },
	{ 2, 0 },
	{ 0, 1 },
};

static const char relTriTri[3][3] = {
	{ 0, 2, 1 }, //  0
	{ 2, 1, 0 }, //  1 
	{ 1, 0, 2 }, //  2
};

static const char indexRot3[3][3] = { 
	{ 0, 1, 2 },
	{ 1, 2, 0 },
	{ 2, 0, 1 },
};

static const char invRelTriTri[3] = {
	0, 1, 2
};

// (0に相当するノード番号)*3+(1に相当するノード番号)  →→　関係番号
static const char noel2RelTriTri[9] = {
	-1,	// 0 00
	-1,	// 1 01
	 0,	// 2 02
	 2, // 3 10
	-1, // 4 11
	-1,	// 5 12
	-1,	// 6 20
	 1,	// 7 21
	-1, // 8 22
};


//! ２次元ベクトルクラス
class CVector2D{
public:
	CVector2D(){}
	CVector2D( const CVector2D& rhs ){
		this->x = rhs.x;
		this->y = rhs.y;
	}
	CVector2D(double x, double y){
		this->x = x;
		this->y = y;
	}
	
	friend double Dot(const CVector2D&, const CVector2D&);

	// オペレータ定義
	inline CVector2D& operator+=(const CVector2D& rhs){
		x += rhs.x;
		y += rhs.y;
		return *this;
	}
	inline CVector2D& operator-=(const CVector2D& rhs){
		x -= rhs.x;
		y -= rhs.y;
		return *this;
	}
	inline CVector2D& operator*=(double scale){
		x *= scale;
		y *= scale;
		return *this;
	}
	inline CVector2D operator*(double scale) const {
		CVector2D v = *this;
		return v *= scale;
	}
	inline CVector2D operator+(const CVector2D& rhs) const {
		CVector2D v = *this;
		return v += rhs;
	}
	inline CVector2D operator-(const CVector2D& rhs) const {
		CVector2D v = *this;
		return v -= rhs;
	}

	inline void Normalize(){
		const double mag = Length();
		x /= mag;
		y /= mag;
	}
	inline void SetZero(){
		x = 0.0;
		y = 0.0;
	}

	//! ベクトルの長さを計算する
	double Length() const{
		return sqrt( x*x+y*y );
	}
	//! ベクトルの長さの２乗を計算する
	double SqLength() const{
		return x*x+y*y;
	}
public:
	double x;
	double y;
};

struct STri2D{
	unsigned int v[3];
	////////////////
	int g2[3];			// 隣接する要素配列ID(-1:隣接要素なし、-2:自分の要素配列に隣接)
	unsigned int s2[3];
	char r2[3];
};
	
class CPoint2D{
public:
	CPoint2D(){}
	CPoint2D( const CPoint2D& rhs )
		: p(rhs.p), e(rhs.e), d(rhs.d){}
	CPoint2D(double x, double y, int ielem, char idir)
		: p(x,y), e(ielem), d(idir){}
public:
	CVector2D p;
	int e;
	char d;
};

inline double TriArea(const CVector2D& v1, const CVector2D& v2, const CVector2D& v3){
	return 0.5*( (v2.x-v1.x)*(v3.y-v1.y) - (v3.x-v1.x)*(v2.y-v1.y) );
}

inline double SquareLength(const CVector2D& ipo0, const CVector2D& ipo1){
	return	( ipo1.x - ipo0.x )*( ipo1.x - ipo0.x ) + ( ipo1.y - ipo0.y )*( ipo1.y - ipo0.y );
}

inline double Length(const CVector2D& ipo0, const CVector2D& ipo1){
	return	sqrt( ( ipo1.x - ipo0.x )*( ipo1.x - ipo0.x ) + ( ipo1.y - ipo0.y )*( ipo1.y - ipo0.y ) );
}

inline int DetDelaunay(
		const CVector2D& p0, 
		const CVector2D& p1, 
		const CVector2D& p2, 
		const CVector2D& p3)
{
	const double area = TriArea(p0,p1,p2);
	if( fabs(area) < 1.0e-10 ){
		return 3;
	}
	const double tmp_val = 1.0/(area*area*16.0);

	const double dtmp0 = SquareLength(p1,p2);
	const double dtmp1 = SquareLength(p0,p2);
	const double dtmp2 = SquareLength(p0,p1);

	const double etmp0 = tmp_val*dtmp0*(dtmp1+dtmp2-dtmp0);
	const double etmp1 = tmp_val*dtmp1*(dtmp0+dtmp2-dtmp1);
	const double etmp2 = tmp_val*dtmp2*(dtmp0+dtmp1-dtmp2);

	const CVector2D out_center(
		etmp0*p0.x + etmp1*p1.x + etmp2*p2.x,
		etmp0*p0.y + etmp1*p1.y + etmp2*p2.y );

	const double qradius = SquareLength(out_center,p0);
	const double qdistance = SquareLength(out_center,p3);

//	assert( fabs( qradius - SquareLength(out_center,p1) ) < 1.0e-10*qradius );
//	assert( fabs( qradius - SquareLength(out_center,p2) ) < 1.0e-10*qradius );

	const double tol = 1.0e-20;
	if( qdistance > qradius*(1.0+tol) ){ return 2; }	// 外接円の外
	else{
		if( qdistance < qradius*(1.0-tol) ){ return 0; }	// 外接円の中
		else{ return 1;	}	// 外接円上
	}
	return 0;
}

////////////////////////////////


bool FlipEdge( unsigned int itri0, unsigned int ied0,
			  std::vector<CPoint2D>& aPo, std::vector<STri2D>& aTri )
{
	assert( itri0 < aTri.size() );
	assert( ied0 < 3 );
	assert( aTri[itri0].g2[ied0] == -2 );

	const unsigned int itri1 = aTri[itri0].s2[ied0];
	const unsigned int ied1  = relTriTri[ aTri[itri0].r2[ied0] ][ied0];
	assert( itri1 < aTri.size() );
	assert( ied1 < 3 );
	assert( aTri[itri1].g2[ied1] == -2 );

//	std::cout << itri0 << "-" << ied0 << "    " << itri1 << "-" << ied1 << std::endl;

	STri2D old0 = aTri[itri0];
	STri2D old1 = aTri[itri1];

	const unsigned int no0_0 = ied0;
	const unsigned int no1_0 = noelTriEdge[ied0][0];
	const unsigned int no2_0 = noelTriEdge[ied0][1];

	const unsigned int no0_1 = ied1;
	const unsigned int no1_1 = noelTriEdge[ied1][0];
	const unsigned int no2_1 = noelTriEdge[ied1][1];
	
	assert( old0.v[no1_0] == old1.v[no2_1] );
	assert( old0.v[no2_0] == old1.v[no1_1] );

	aPo[ old0.v[no1_0] ].e = itri0;	aPo[ old0.v[no1_0] ].d = 0;
	aPo[ old0.v[no0_0] ].e = itri0;	aPo[ old0.v[no0_0] ].d = 2;
	aPo[ old1.v[no1_1] ].e = itri1;	aPo[ old1.v[no1_1] ].d = 0;
	aPo[ old1.v[no0_1] ].e = itri1;	aPo[ old1.v[no0_1] ].d = 2;

	{
		STri2D& ref_tri = aTri[itri0];
		////////////////
		ref_tri.v[0]  = old0.v[no1_0];	ref_tri.v[1]  = old1.v[no0_1];	ref_tri.v[2]  = old0.v[no0_0];
		ref_tri.g2[0] = -2;				ref_tri.g2[1] = old0.g2[no2_0];	ref_tri.g2[2] = old1.g2[no1_1];
		ref_tri.s2[0] = itri1;			ref_tri.s2[1] = old0.s2[no2_0];	ref_tri.s2[2] = old1.s2[no1_1];
		////////////////
		ref_tri.r2[0] = 0;
		if( old0.g2[no2_0] == -2 || old0.g2[no2_0] == -3 ){
			assert( old0.r2[no2_0] < 3 );
			const char* rel = relTriTri[ old0.r2[no2_0] ];
			assert( old0.s2[no2_0] < aTri.size() );
			assert( old0.s2[no2_0] != itri0 );
			assert( old0.s2[no2_0] != itri1 );
			ref_tri.r2[1] = noel2RelTriTri[ rel[no1_0]*3 + rel[no2_0] ];
			assert( ref_tri.r2[1] >= 0 && ref_tri.r2[1] < 3 );
			aTri[ old0.s2[no2_0] ].s2[ rel[no2_0] ] = itri0;
			aTri[ old0.s2[no2_0] ].r2[ rel[no2_0] ] = invRelTriTri[ ref_tri.r2[1] ];
		}
		if( old1.g2[no1_1] == -2 || old1.g2[no1_1] == -3 ){
			assert( old1.r2[no1_1] < 3 );
			const char* rel = relTriTri[ old1.r2[no1_1] ];
			assert( old1.s2[no1_1] < aTri.size() );
			ref_tri.r2[2] = noel2RelTriTri[ rel[no2_1]*3 + rel[no0_1] ];
			assert( ref_tri.r2[2] >= 0 && ref_tri.r2[2] < 3 );
			aTri[ old1.s2[no1_1] ].s2[ rel[no1_1] ] = itri0;
			aTri[ old1.s2[no1_1] ].r2[ rel[no1_1] ] = invRelTriTri[ ref_tri.r2[2] ];			
		}
	}

	{
		STri2D& ref_tri = aTri[itri1];
		////////////////
		ref_tri.v[0] = old1.v[no1_1];	ref_tri.v[1]  = old0.v[no0_0];	ref_tri.v[2]  = old1.v[no0_1];
		ref_tri.g2[0] = -2;				ref_tri.g2[1] = old1.g2[no2_1];	ref_tri.g2[2] = old0.g2[no1_0];
		ref_tri.s2[0] = itri0;			ref_tri.s2[1] = old1.s2[no2_1];	ref_tri.s2[2] = old0.s2[no1_0];
		////////////////
		ref_tri.r2[0] = 0;
		if( old1.g2[no2_1] == -2 || old1.g2[no2_1] == -3 ){
			assert( old1.r2[no2_1] < 3 );
			const char* rel = relTriTri[ old1.r2[no2_1] ];
			assert( old1.s2[no2_1] < aTri.size() );
			ref_tri.r2[1] = noel2RelTriTri[ rel[no1_1]*3 + rel[no2_1] ];
			assert( ref_tri.r2[1] >= 0 && ref_tri.r2[1] < 3 );
			aTri[ old1.s2[no2_1] ].s2[ rel[no2_1] ] = itri1;
			aTri[ old1.s2[no2_1] ].r2[ rel[no2_1] ] = invRelTriTri[ ref_tri.r2[1] ];
		}
		if( old0.g2[no1_0] == -2 || old0.g2[no1_0] == -3 ){
			assert( old0.r2[no1_0] < 3 );
			const char* rel = relTriTri[ old0.r2[no1_0] ];
			assert( old0.s2[no1_0] < aTri.size() );
			ref_tri.r2[2] = noel2RelTriTri[ rel[no2_0]*3 + rel[no0_0] ];
			assert( ref_tri.r2[2] >= 0 && ref_tri.r2[2] < 3 );
			aTri[ old0.s2[no1_0] ].s2[ rel[no1_0] ] = itri1;
			aTri[ old0.s2[no1_0] ].r2[ rel[no1_0] ] = invRelTriTri[ ref_tri.r2[2] ];
		}
	}
	
	return true;
}

bool InsertPoint_ElemEdge( const unsigned int& ipo_ins, 
	const unsigned int& itri_ins, const unsigned int& ied_ins,
	std::vector<CPoint2D>& po, std::vector<STri2D>& tri )
{
	assert( itri_ins < tri.size() );
	assert( ipo_ins < po.size() );

	if( tri[itri_ins].g2[ied_ins] != -2 ){
//		std::cout << "未実装" << std::endl;
		assert(0);
	}

	const unsigned int itri_adj = tri[itri_ins].s2[ied_ins];
	const unsigned int ied_adj  = relTriTri[ tri[itri_ins].r2[ied_ins] ][ied_ins];
	assert( itri_adj < tri.size() );
	assert( ied_ins < 3 );

	const int itri0 = itri_ins;
	const int itri1 = itri_adj;
	const int itri2 = tri.size();
	const int itri3 = tri.size()+1;

	tri.resize( tri.size()+2 );

	STri2D old0 = tri[itri_ins];
	STri2D old1 = tri[itri_adj];

	const unsigned int ino0_0 = ied_ins;
	const unsigned int ino1_0 = noelTriEdge[ied_ins][0];
	const unsigned int ino2_0 = noelTriEdge[ied_ins][1];

	const unsigned int ino0_1 = ied_adj;
	const unsigned int ino1_1 = noelTriEdge[ied_adj][0];
	const unsigned int ino2_1 = noelTriEdge[ied_adj][1];
	
	assert( old0.v[ino1_0] == old1.v[ino2_1] );
	assert( old0.v[ino2_0] == old1.v[ino1_1] );
	assert( old0.s2[ino0_0 ] == itri1 );
	assert( old1.s2[ino0_1 ] == itri0 );

	po[ipo_ins].e = itri0;			po[ipo_ins].d = 0;
	po[ old0.v[ino2_0] ].e = itri0;	po[ old0.v[ino2_0] ].d = 1;
	po[ old0.v[ino0_0] ].e = itri1;	po[ old0.v[ino0_0] ].d = 1;
	po[ old1.v[ino2_1] ].e = itri2;	po[ old1.v[ino2_1] ].d = 1;
	po[ old1.v[ino0_1] ].e = itri3;	po[ old1.v[ino0_1] ].d = 1;

	{
		STri2D& ref_tri = tri[itri0];
		////////////////
		ref_tri.v[0]  = ipo_ins;			ref_tri.v[1]  = old0.v[ino2_0];	ref_tri.v[2]  = old0.v[ino0_0];
		ref_tri.g2[0] = old0.g2[ino1_0];	ref_tri.g2[1] = -2;				ref_tri.g2[2] = -2;
		ref_tri.s2[0] = old0.s2[ino1_0];	ref_tri.s2[1] = itri1;			ref_tri.s2[2] = itri3;
		////////////////
		if( old0.g2[ino1_0] == -2 || old0.g2[ino1_0] == -3 ){
			assert( old0.r2[ino1_0] < 3 );
			const char* rel = relTriTri[ old0.r2[ino1_0] ];
			ref_tri.r2[0] = noel2RelTriTri[ rel[ino1_0]*3 + rel[ino2_0] ];
			assert( ref_tri.r2[0] >= 0 && ref_tri.r2[0] < 3 );
			assert( old0.s2[ino1_0] < tri.size() );
			tri[ old0.s2[ino1_0] ].s2[ rel[ino1_0] ] = itri0;
			tri[ old0.s2[ino1_0] ].r2[ rel[ino1_0] ] = invRelTriTri[ ref_tri.r2[0] ];
		}
		ref_tri.r2[1] = 0;
		ref_tri.r2[2] = 0;
	}
	{
		STri2D& ref_tri = tri[itri1];
		////////////////
		ref_tri.v[0]  = ipo_ins;			ref_tri.v[1]  = old0.v[ino0_0];	ref_tri.v[2]  = old0.v[ino1_0];
		ref_tri.g2[0] = old0.g2[ino2_0];	ref_tri.g2[1] = -2;				ref_tri.g2[2] = -2;
		ref_tri.s2[0] = old0.s2[ino2_0];	ref_tri.s2[1] = itri2;			ref_tri.s2[2] = itri0;
		////////////////
		if( old0.g2[ino2_0] == -2 || old0.g2[ino2_0] == -3 ){
			assert( old0.r2[ino2_0] < 3 );
			const char* rel = relTriTri[ old0.r2[ino2_0] ];
			ref_tri.r2[0] = noel2RelTriTri[ rel[ino2_0]*3 + rel[ino0_0] ];
			assert( ref_tri.r2[0] >= 0 && ref_tri.r2[0] < 3 );
			assert( old0.s2[ino2_0] < tri.size() );
			tri[ old0.s2[ino2_0] ].s2[ rel[ino2_0] ] = itri1;
			tri[ old0.s2[ino2_0] ].r2[ rel[ino2_0] ] = invRelTriTri[ ref_tri.r2[0] ];
		}
		ref_tri.r2[1] = 0;
		ref_tri.r2[2] = 0;
	}
	{
		STri2D& ref_tri = tri[itri2];
		////////////////
		ref_tri.v[0]  = ipo_ins;			ref_tri.v[1]  = old1.v[ino2_1];	ref_tri.v[2]  = old1.v[ino0_1];
		ref_tri.g2[0] = old1.g2[ino1_1];	ref_tri.g2[1] = -2;				ref_tri.g2[2] = -2;
		ref_tri.s2[0] = old1.s2[ino1_1];	ref_tri.s2[1] = itri3;			ref_tri.s2[2] = itri1;
		////////////////
		if( old1.g2[ino1_1] == -2 || old0.g2[ino2_0] == -3 ){
			assert( old1.r2[ino1_1] < 3 );
			const char* rel = relTriTri[ old1.r2[ino1_1] ];
			ref_tri.r2[0] = noel2RelTriTri[ rel[ino1_1]*3 + rel[ino2_1] ];
			assert( ref_tri.r2[0] >= 0 && ref_tri.r2[0] < 3 );
			assert( old1.s2[ino1_1] < tri.size() );
			tri[ old1.s2[ino1_1] ].s2[ rel[ino1_1] ] = itri2;
			tri[ old1.s2[ino1_1] ].r2[ rel[ino1_1] ] = invRelTriTri[ ref_tri.r2[0] ];
		}
		ref_tri.r2[1] = 0;
		ref_tri.r2[2] = 0;
	}
	{
		STri2D& ref_tri = tri[itri3];
		ref_tri.v[0]  = ipo_ins;			ref_tri.v[1]  = old1.v[ino0_1];	ref_tri.v[2]  = old1.v[ino1_1];
		ref_tri.g2[0] = old1.g2[ino2_1];	ref_tri.g2[1] = -2;				ref_tri.g2[2] = -2;
		ref_tri.s2[0] = old1.s2[ino2_1];	ref_tri.s2[1] = itri0;			ref_tri.s2[2] = itri2;
		if( old1.g2[ino2_1] == -2 || old1.g2[ino2_1] == -3 ){
			assert( old1.r2[ino2_1] < 3 );
			const char* rel = relTriTri[ old1.r2[ino2_1] ];
			ref_tri.r2[0] = noel2RelTriTri[ rel[ino2_1]*3 + rel[ino0_1] ];
			assert( ref_tri.r2[0] >= 0 && ref_tri.r2[0] < 3 );
			assert( old1.s2[ino2_1] < tri.size() );
			tri[ old1.s2[ino2_1] ].s2[ rel[ino2_1] ] = itri3;
			tri[ old1.s2[ino2_1] ].r2[ rel[ino2_1] ] = invRelTriTri[ ref_tri.r2[0] ];
		}
		ref_tri.r2[1] = 0;
		ref_tri.r2[2] = 0;
	}
	return true;
}

bool InsertPoint_Elem( const unsigned int& ipo_ins, const unsigned int& itri_ins, 
					  std::vector<CPoint2D>& po, std::vector<STri2D>& tri )
{
	assert( itri_ins < tri.size() );
	assert( ipo_ins < po.size() );

	const int itri0 = itri_ins;
	const int itri1 = tri.size();
	const int itri2 = tri.size()+1;

	tri.resize( tri.size()+2 );

	const STri2D old_tri = tri[itri_ins];

	po[ipo_ins].e = itri0;			po[ipo_ins].d = 0;
	po[ old_tri.v[0] ].e = itri1;	po[ old_tri.v[0] ].d = 2;
	po[ old_tri.v[1] ].e = itri2;	po[ old_tri.v[1] ].d = 2;
	po[ old_tri.v[2] ].e = itri0;	po[ old_tri.v[2] ].d = 2;

	{
		STri2D& ref_tri = tri[itri0];

		ref_tri.v[0] = ipo_ins;			ref_tri.v[1] = old_tri.v[1];	ref_tri.v[2] = old_tri.v[2];
		ref_tri.g2[0] = old_tri.g2[0];	ref_tri.g2[1] = -2;				ref_tri.g2[2] = -2;
		ref_tri.s2[0] = old_tri.s2[0];	ref_tri.s2[1] = itri1;			ref_tri.s2[2] = itri2;

		if( old_tri.g2[0] == -2 || old_tri.g2[0] == -3 ){
			assert( old_tri.r2[0] < 3 );
			const char* rel = relTriTri[ old_tri.r2[0] ];
			ref_tri.r2[0] = noel2RelTriTri[ rel[0]*3 + rel[1] ];
			assert( ref_tri.r2[0] >= 0 && ref_tri.r2[0] < 3 );
			assert( old_tri.s2[0] < tri.size() );
			tri[ old_tri.s2[0] ].s2[ rel[0] ] = itri0;
			tri[ old_tri.s2[0] ].r2[ rel[0] ] = invRelTriTri[ ref_tri.r2[0] ];
		}
		ref_tri.r2[1] = 0;
		ref_tri.r2[2] = 0;
	}
	{
		STri2D& ref_tri = tri[itri1];

		ref_tri.v[0] = ipo_ins;			ref_tri.v[1] = old_tri.v[2];	ref_tri.v[2] = old_tri.v[0];
		ref_tri.g2[0] = old_tri.g2[1];	ref_tri.g2[1] = -2;				ref_tri.g2[2] = -2;
		ref_tri.s2[0] = old_tri.s2[1];	ref_tri.s2[1] = itri2;			ref_tri.s2[2] = itri0;

		if( old_tri.g2[1] == -2 || old_tri.g2[1] == -3 ){
			assert( old_tri.r2[1] < 3 );
			const char* rel = relTriTri[ old_tri.r2[1] ];
			ref_tri.r2[0] = noel2RelTriTri[ rel[1]*3 + rel[2] ];
			assert( ref_tri.r2[0] >= 0 && ref_tri.r2[0] < 3 );
			assert( old_tri.s2[1] < tri.size() );
			tri[ old_tri.s2[1] ].s2[ rel[1] ] = itri1;
			tri[ old_tri.s2[1] ].r2[ rel[1] ] = invRelTriTri[ ref_tri.r2[0] ];
		}
		ref_tri.r2[1] = 0;
		ref_tri.r2[2] = 0;
	}
	{
		STri2D& ref_tri = tri[itri2];

		ref_tri.v[0] = ipo_ins; 		ref_tri.v[1] = old_tri.v[0];	ref_tri.v[2] = old_tri.v[1];
		ref_tri.g2[0] = old_tri.g2[2];  ref_tri.g2[1] = -2;				ref_tri.g2[2] = -2;
		ref_tri.s2[0] = old_tri.s2[2];	ref_tri.s2[1] = itri0;			ref_tri.s2[2] = itri1;

		if( old_tri.g2[2] == -2 || old_tri.g2[2] == -3 ){
			assert( old_tri.r2[2] < 3 );
			const char* rel = relTriTri[ old_tri.r2[2] ];
			ref_tri.r2[0] = noel2RelTriTri[ rel[2]*3 + rel[0] ];
			assert( ref_tri.r2[0] >= 0 && ref_tri.r2[0] < 3 );
			assert( old_tri.s2[2] < tri.size() );
			tri[ old_tri.s2[2] ].s2[ rel[2] ] = itri2;
			tri[ old_tri.s2[2] ].r2[ rel[2] ] = invRelTriTri[ ref_tri.r2[0] ];
		}
		ref_tri.r2[1] = 0;
		ref_tri.r2[2] = 0;
	}

	return true;
}


// 辺[ipo0-ipo1]の左側の３角形itri0を探索する
// 三角形がなければ->falseを返す。
// 三角形があれば  ->true を返す。
// 但し、その場合
// tri[itri0].v[inotri0]==ipo0
// tri[itri0].v[inotri1]==ipo1
// を満たす
bool FindEdge( const unsigned int& ipo0, const unsigned int& ipo1,
	unsigned int& itri0, unsigned int& inotri0, unsigned int& inotri1,
	const std::vector<CPoint2D>& po, const std::vector<STri2D>& tri )
{
	const unsigned int itri_ini = po[ipo0].e;
	const unsigned int inotri_ini = po[ipo0].d;
	unsigned int inotri_cur = inotri_ini;
	unsigned int itri_cur = itri_ini;
	for(;;){	//　時計周りに検索する。
		assert( tri[itri_cur].v[inotri_cur] == ipo0 );
		{	// この要素がOKか調べる
			const unsigned int inotri2 = indexRot3[1][inotri_cur];
			if( tri[itri_cur].v[inotri2] == ipo1 ){
				itri0 = itri_cur;
				inotri0 = inotri_cur;
				inotri1 = inotri2;
				assert( tri[itri0].v[ inotri0 ] == ipo0 );
				assert( tri[itri0].v[ inotri1 ] == ipo1 );
				return true;
			}
		}
		{	// 次の要素へ進める
			const unsigned int inotri2 = indexRot3[2][inotri_cur];
			if( tri[itri_cur].g2[inotri2] != -2 && tri[itri_cur].g2[inotri2] != -3 ){ break; }
			const unsigned int itri_nex = tri[itri_cur].s2[inotri2];
			const char* rel = relTriTri[ tri[itri_cur].r2[inotri2] ];
			const unsigned int inotri3 = rel[inotri_cur];
			assert( tri[itri_nex].v[inotri3] == ipo0 );
			if( itri_nex == itri_ini ) return false;
			itri_cur = itri_nex;
			inotri_cur = inotri3;
		}
	}

	inotri_cur = inotri_ini;
	itri_cur = itri_ini;
	for(;;){	//　反時計周りの検索
		assert( tri[itri_cur].v[inotri_cur] == ipo0 );
		{	// 次の要素へ進める
			const unsigned int inotri2 = indexRot3[1][inotri_cur];
			if( tri[itri_cur].g2[inotri2] != -2 || tri[itri_cur].g2[inotri2] != -3 ){ break; }
			const unsigned int itri_nex = tri[itri_cur].s2[inotri2];
			const char* rel = relTriTri[ tri[itri_cur].r2[inotri2] ];
			const unsigned int inotri3 = rel[inotri_cur];
			assert( tri[itri_nex].v[inotri3] == ipo0 );
			if( itri_nex == itri_ini ){	// 一周したら終わり
				itri0 = 0;
				inotri0 = 0; inotri1 = 0;
				return false;
			}
			itri_cur = itri_nex;
			inotri_cur = inotri3;
		}
		{	// 要素の向きを調べる
			const unsigned int inotri2 = indexRot3[1][inotri_cur];
			if( tri[itri_cur].v[inotri2] == ipo1 ){
				itri0 = itri_cur;
				inotri0 = inotri_cur;
				inotri1 = inotri2;
				assert( tri[itri0].v[ inotri0 ] == ipo0 );
				assert( tri[itri0].v[ inotri1 ] == ipo1 );
				return true;
			}
		}
	}

	return false;
}

bool FindEdgePoint_AcrossEdge( const unsigned int& ipo0, const unsigned int& ipo1,
	unsigned int& itri0, unsigned int& inotri0, unsigned int& inotri1, double& ratio,
	std::vector<CPoint2D>& po, std::vector<STri2D>& tri )
{
	const unsigned int itri_ini = po[ipo0].e;
	const unsigned int inotri_ini = po[ipo0].d;
	unsigned int inotri_cur = inotri_ini;
	unsigned int itri_cur = itri_ini;
	for(;;){	//　反時計周りの検索
		assert( tri[itri_cur].v[inotri_cur] == ipo0 );
		{
			const unsigned int inotri2 = indexRot3[1][inotri_cur];
			const unsigned int inotri3 = indexRot3[2][inotri_cur];
			double area0 = TriArea( po[ipo0].p, po[ tri[itri_cur].v[inotri2] ].p, po[ipo1].p );
			if( area0 > -1.0e-20 ){
				double area1 =  TriArea( po[ipo0].p, po[ipo1].p, po[ tri[itri_cur].v[inotri3] ].p );
				if( area1 > -1.0e-20 ){
					assert( area0 + area1 > 1.0e-20 );
					ratio = area0 / ( area0 + area1 );
					itri0 = itri_cur;
					inotri0 = inotri2;
					inotri1 = inotri3;
					return true;
				}
			}
		}
		{	// 次の要素へ進める
			const unsigned int inotri2 = indexRot3[1][inotri_cur];
			if( tri[itri_cur].g2[inotri2] != -2 && tri[itri_cur].g2[inotri2] != -3 ){
				break;
			}
			unsigned int itri_nex = tri[itri_cur].s2[inotri2];
			const char* rel = relTriTri[ tri[itri_cur].r2[inotri2] ];
			const unsigned int inotri3 = rel[inotri_cur];
			assert( tri[itri_nex].v[inotri3] == ipo0 );
			if( itri_nex == itri_ini ){	// 一周したら終わり
				itri0 = 0;
				inotri0 = 0; inotri1 = 0;
				ratio = 0.0;
				return false;
			}
			itri_cur = itri_nex;
			inotri_cur = inotri3;
		}
	}

	// itri_iniを２回計算しているので少し無駄、いつか直す

	inotri_cur = inotri_ini;
	itri_cur = itri_ini;
	for(;;){	//　時計周りに検索する。
		assert( tri[itri_cur].v[inotri_cur] == ipo0 );
		{
			const unsigned int inotri2 = indexRot3[1][inotri_cur];
			const unsigned int inotri3 = indexRot3[2][inotri_cur];
			double area0 = TriArea( po[ipo0].p, po[ tri[itri_cur].v[inotri2] ].p, po[ipo1].p );
			if( area0 > -1.0e-20 ){
				double area1 =  TriArea( po[ipo0].p, po[ipo1].p, po[ tri[itri_cur].v[inotri3] ].p );
				if( area1 > -1.0e-20 ){
					assert( area0 + area1 > 1.0e-20 );
					ratio = area0 / ( area0 + area1 );
					itri0 = itri_cur;
					inotri0 = inotri2;
					inotri1 = inotri3;
					return true;
				}
			}
		}
		{	// 次の要素へ進める
			const unsigned int inotri2 = indexRot3[2][inotri_cur];
			if( tri[itri_cur].g2[inotri2] != -2 && tri[itri_cur].g2[inotri2] != -3 ){
				break;
			}
			unsigned int itri_nex = tri[itri_cur].s2[inotri2];
			const char* rel = relTriTri[ tri[itri_cur].r2[inotri2] ];
			const unsigned int inotri3 = rel[inotri_cur];
			assert( tri[itri_nex].v[inotri3] == ipo0 );
			if( itri_nex == itri_ini ){
				assert(0);	// 一周しないはず
			}
			itri_cur = itri_nex;
			inotri_cur = inotri3;
		}
	}

	// 失敗したときの値を入れる
	itri0 = 0;
	inotri0 = 0; inotri1 = 0;
	ratio = 0.0;

	return false;
}

////////////////////////////////

bool DelaunayAroundPoint( unsigned int ipo0,
	std::vector<CPoint2D>& aPo, std::vector<STri2D>& aTri )
{
	assert( ipo0 < aPo.size() );
	if( aPo[ipo0].e == -1 ) return true;

	assert( aPo[ipo0].e >= 0 && (unsigned int)aPo[ipo0].e < aTri.size() );
	assert( aTri[ aPo[ipo0].e ].v[ aPo[ipo0].d ] == ipo0 );

	const unsigned int itri0 = aPo[ipo0].e;
	unsigned int inotri0 = aPo[ipo0].d;

	unsigned int itri_cur = itri0;
	unsigned int inotri_cur = aPo[ipo0].d;
	bool flag_is_wall = false;
	for(;;){
		assert( aTri[itri_cur].v[inotri_cur] == ipo0 );

		if( aTri[itri_cur].g2[inotri_cur] == -2 ){
			// 向かいの要素を調べる
			const unsigned int itri_dia = aTri[itri_cur].s2[inotri_cur];
			const char* rel_dia = relTriTri[ aTri[itri_cur].r2[inotri_cur] ];
			const unsigned int inotri_dia = rel_dia[inotri_cur];
			assert( aTri[itri_dia].g2[inotri_dia] == -2 );
			assert( aTri[itri_dia].s2[inotri_dia] == itri_cur );
			const unsigned int ipo_dia = aTri[itri_dia].v[inotri_dia];
			if( DetDelaunay(
				aPo[ aTri[itri_cur].v[0] ].p,
				aPo[ aTri[itri_cur].v[1] ].p,
				aPo[ aTri[itri_cur].v[2] ].p,
				aPo[ ipo_dia ].p ) == 0 )	// Delaunay条件が満たされない場合
			{
				FlipEdge(itri_cur,inotri_cur,aPo,aTri);	// 辺を切り替える
				// FlipEdgeによってitri_curは時計回り側の３角形に切り替わる
				inotri_cur = 2;
				assert( aTri[itri_cur].v[inotri_cur] == ipo0 );
				// FlipによってaTri[itri0].v[inotri0] != ipo0 でなくなってしまうのを防ぐため
				if( itri_cur == itri0 ) inotri0 = inotri_cur;
				continue;	// ループの始めに戻る
			}
		}

		{	// 次の要素へ進める
			const unsigned int inotri1 = indexRot3[1][inotri_cur];
			if( aTri[itri_cur].g2[inotri1] != -2 && aTri[itri_cur].g2[inotri1] != -3 ){
				flag_is_wall = true;
				break;
			}
			const unsigned int itri_nex = aTri[itri_cur].s2[inotri1];
			const char* rel_nex = relTriTri[ aTri[itri_cur].r2[inotri1] ];
			const unsigned int inotri_nex = rel_nex[inotri_cur];
			assert( aTri[itri_nex].v[inotri_nex] == ipo0 );
			if( itri_nex == itri0 ) break;	// 一周したら終わり
			itri_cur = itri_nex;
			inotri_cur = inotri_nex;
		}
	}
	if( !flag_is_wall ) return true;

	////////////////////////////////
	// 逆向きへの回転

	itri_cur = itri0;
	inotri_cur = inotri0;
	for(;;){
		assert( aTri[itri_cur].v[inotri_cur] == ipo0 );

		if( aTri[itri_cur].g2[inotri_cur] == -2 ){
			// 向かいの要素を調べる
			const unsigned int itri_dia = aTri[itri_cur].s2[inotri_cur];
			const char* rel_dia = relTriTri[ aTri[itri_cur].r2[inotri_cur] ];
			const unsigned int inotri_dia = rel_dia[inotri_cur];
			assert( aTri[itri_dia].g2[inotri_dia] == -2 );
			assert( aTri[itri_dia].s2[inotri_dia] == itri_cur );
			const unsigned int ipo_dia = aTri[itri_dia].v[inotri_dia];
			if( DetDelaunay(
				aPo[ aTri[itri_cur].v[0] ].p,
				aPo[ aTri[itri_cur].v[1] ].p,
				aPo[ aTri[itri_cur].v[2] ].p,
				aPo[ ipo_dia ].p ) == 0 )	// Delaunay条件が満たされない場合
			{
				FlipEdge(itri_cur,inotri_cur,aPo,aTri);	// 辺を切り替える
				itri_cur = itri_dia;
				inotri_cur = 1;
				assert( aTri[itri_cur].v[inotri_cur] == ipo0 );
				continue;	// ループの始めに戻る
			}
		}

		{	// 次の要素へ進める
			const unsigned int inotri2 = indexRot3[2][inotri_cur];
			if( aTri[itri_cur].g2[inotri2] != -2 && aTri[itri_cur].g2[inotri2] != -3 ){
				return true;
			}
			const unsigned int itri_nex = aTri[itri_cur].s2[inotri2];
			const char* rel_nex = relTriTri[ aTri[itri_cur].r2[inotri2] ];
			const unsigned int inotri_nex = rel_nex[inotri_cur];
			assert( aTri[itri_nex].v[inotri_nex] == ipo0 );
			assert( itri_nex != itri0 );	// 一周したら終わり
			itri_cur = itri_nex;
			inotri_cur = inotri_nex;
		}
	}

	return true;
}


bool TriangulateOuterLoop(std::vector<CPoint2D>& aPo2D, std::vector<STri2D>& aTri_in, 
						  const std::vector<unsigned int>& aPtrVtxInd, const std::vector<unsigned int>& aVtxInd)
{
	std::vector<STri2D> aTri;	// 内側の三角形
	std::vector<unsigned int> aPoDel;	// 最終的に消す頂点

	{	// 与えられた点群を内部に持つ、大きな三角形を作る
		double max_len;
		double center[2];
		{
			double bound_2d[4];
			bound_2d[0] = aPo2D[0].p.x;
			bound_2d[1] = aPo2D[0].p.x;
			bound_2d[2] = aPo2D[0].p.y;
			bound_2d[3] = aPo2D[0].p.y;
			for(unsigned int ipoin=1;ipoin<aPo2D.size();ipoin++){
				if( aPo2D[ipoin].p.x < bound_2d[0] ){ bound_2d[0] = aPo2D[ipoin].p.x; }
				if( aPo2D[ipoin].p.x > bound_2d[1] ){ bound_2d[1] = aPo2D[ipoin].p.x; }
				if( aPo2D[ipoin].p.y < bound_2d[2] ){ bound_2d[2] = aPo2D[ipoin].p.y; }
				if( aPo2D[ipoin].p.y > bound_2d[3] ){ bound_2d[3] = aPo2D[ipoin].p.y; }
			}
			max_len = (bound_2d[1]-bound_2d[0]>bound_2d[3]-bound_2d[2]) ? bound_2d[1]-bound_2d[0] : bound_2d[3]-bound_2d[2];
			center[0] = (bound_2d[1]+bound_2d[0])*0.5;
			center[1] = (bound_2d[3]+bound_2d[2])*0.5;
		}

//		std::cout << "center " << center[0] << " " << center[1] << std::endl;
//		std::cout << "max_len " << max_len << std::endl;

		const double tri_len = max_len * 4.0;
		const double tmp_len = tri_len * sqrt(3.0) / 6.0;

		const int npo = aPo2D.size();
		aPoDel.push_back( aPo2D.size()+0 );
		aPoDel.push_back( aPo2D.size()+1 );
		aPoDel.push_back( aPo2D.size()+2 );
		aPo2D.resize(npo+3);
		aPo2D[npo+0].p.x = center[0];				aPo2D[npo+0].p.y = center[1]+2.0*tmp_len;	aPo2D[npo+0].e = 0;	aPo2D[npo+0].d = 0;
		aPo2D[npo+1].p.x = center[0]-0.5*tri_len;	aPo2D[npo+1].p.y = center[1]-tmp_len;		aPo2D[npo+1].e = 0;	aPo2D[npo+1].d = 1;
		aPo2D[npo+2].p.x = center[0]+0.5*tri_len;	aPo2D[npo+2].p.y = center[1]-tmp_len;		aPo2D[npo+2].e = 0;	aPo2D[npo+2].d = 2;

		aTri.resize(1);
		aTri[0].v[0] = npo+0;	aTri[0].v[1] = npo+1;	aTri[0].v[2] = npo+2;
		aTri[0].g2[0] = -1; aTri[0].g2[1] = -1; aTri[0].g2[2] = -1;
		aTri[0].s2[0] =  0; aTri[0].s2[1] =  0; aTri[0].s2[2] =  0;
		aTri[0].r2[0] =  0; aTri[0].r2[1] =  0; aTri[0].r2[2] =  0;
	}

	// Make Delaunay Division
	for(unsigned int ipoin=0;ipoin<aPo2D.size();ipoin++){
		if( aPo2D[ipoin].e >= 0 ) continue;	// 既にメッシュの一部である。
		const CVector2D& po_add = aPo2D[ipoin].p;
		int itri_in = -1;
		int iedge = -1;
		unsigned int iflg1 = 0, iflg2 = 0;
		for(unsigned int itri=0;itri<aTri.size();itri++){
			iflg1 = 0; iflg2 = 0;
			const STri2D& ref_tri = aTri[itri];
			if( TriArea(po_add, aPo2D[ref_tri.v[1]].p, aPo2D[ref_tri.v[2]].p ) > MIN_TRI_AREA ){
				iflg1++; iflg2 += 0;
			}
			if( TriArea(po_add, aPo2D[ref_tri.v[2]].p, aPo2D[ref_tri.v[0]].p ) > MIN_TRI_AREA ){
				iflg1++; iflg2 += 1;
			}
			if( TriArea(po_add, aPo2D[ref_tri.v[0]].p, aPo2D[ref_tri.v[1]].p ) > MIN_TRI_AREA ){
				iflg1++; iflg2 += 2;
			}
			if( iflg1 == 3 ){
				itri_in = itri;
				break;
			}
			else if( iflg1 == 2 ){
				const unsigned int ied0 = 3-iflg2;
				const unsigned int ipo_e0 = ref_tri.v[ noelTriEdge[ied0][0] ];
				const unsigned int ipo_e1 = ref_tri.v[ noelTriEdge[ied0][1] ];
				const char* rel = relTriTri[ ref_tri.r2[ied0] ];
				const unsigned int itri_s = ref_tri.s2[ied0];
				assert( aTri[itri_s].v[ rel[ noelTriEdge[ied0][0] ] ] == ipo_e0 );
				assert( aTri[itri_s].v[ rel[ noelTriEdge[ied0][1] ] ] == ipo_e1 );
				const unsigned int inoel_d = rel[ied0];
				assert( aTri[itri_s].s2[inoel_d] == itri );
				const unsigned int ipo_d = aTri[itri_s].v[inoel_d];
				assert( TriArea( po_add, aPo2D[ipo_e1].p, aPo2D[ aTri[itri].v[ied0] ].p ) > MIN_TRI_AREA );
				assert( TriArea( po_add, aPo2D[ aTri[itri].v[ied0] ].p, aPo2D[ipo_e0].p ) > MIN_TRI_AREA );
				if( TriArea( po_add, aPo2D[ipo_e0].p, aPo2D[ipo_d ].p ) < MIN_TRI_AREA ){ continue;	}
				if( TriArea( po_add, aPo2D[ipo_d ].p, aPo2D[ipo_e1].p ) < MIN_TRI_AREA ){ continue; }
				const unsigned int det_d =  DetDelaunay(po_add,aPo2D[ipo_e0].p,aPo2D[ipo_e1].p,aPo2D[ipo_d].p);
				if( det_d == 2 || det_d == 1 ) continue;
				itri_in = itri;
				iedge = ied0;
				break;
			}
		}
		if( itri_in == -1 ){
//			std::cout << "Super Triangle Failure " << ipoin << " " << po_add.x << " " << po_add.y << std::endl;
//			std::cout << aTri.size() << std::endl;
			assert(0);
			return false;
		}
		if( iedge == -1 ){
			InsertPoint_Elem(ipoin,itri_in,aPo2D,aTri);
		}
		else{
			InsertPoint_ElemEdge(ipoin,itri_in,iedge,aPo2D,aTri);
		}
		DelaunayAroundPoint(ipoin,aPo2D,aTri);
	}
	
	{	// エッジを回復する
		const unsigned int nloop = aPtrVtxInd.size()-1;
		for(unsigned int iloop=0;iloop<nloop;iloop++){
			const unsigned int nbar = aPtrVtxInd[iloop+1]-aPtrVtxInd[iloop];
			for(unsigned int ibar=0;ibar<nbar;ibar++){
				for(;;){ // EdgeをFlipしたら同じ辺について繰り返す				
					unsigned int ipoi0, ipoi1;
					{
						unsigned int ind0 = aPtrVtxInd[iloop];
						ipoi0 = aVtxInd[ind0+ibar];
						if( ibar != nbar-1 ){ ipoi1 = aVtxInd[ind0+ibar+1]; }
						else{ ipoi1 = aVtxInd[ind0]; }
					}
					assert( ipoi0 < aPo2D.size() ); assert( ipoi1 < aPo2D.size() );
					////////////////
					unsigned int itri0;
					unsigned int inotri0,inotri1;
					if( FindEdge(ipoi0,ipoi1,itri0,inotri0,inotri1,aPo2D,aTri) ){	// ループの内側に接する要素を見つける
						// Split Triangle
						assert( inotri0 != inotri1 );
						assert( inotri0 < 3 );
						assert( inotri1 < 3 );
						assert( aTri[itri0].v[ inotri0 ] == ipoi0 );
						assert( aTri[itri0].v[ inotri1 ] == ipoi1 );
						const unsigned int ied0 = 3 - inotri0 - inotri1;
						{
							const unsigned int itri1 = aTri[itri0].s2[ied0];
							const unsigned int ied1 = relTriTri[ aTri[itri0].r2[ied0] ][ied0];
							assert( aTri[itri1].s2[ied1] == itri0 );
							aTri[itri1].g2[ied1] = -3;	// -3は辺で区切られてるけど、隣接している３角形を指しているということ
							aTri[itri0].g2[ied0] = -3;	// -3は辺で区切られてるけど、隣接している３角形を指しているということ
						}
						break;	// 次のBarへ　for(;;)を抜ける
					}
					else{
						double ratio;
						if( !FindEdgePoint_AcrossEdge(ipoi0,ipoi1,
							itri0,inotri0,inotri1,ratio,
							aPo2D,aTri) ){ assert(0); }
						assert( ratio > -1.0e-20 && ratio < 1.0+1.0e-20 );
						assert( TriArea( aPo2D[ipoi0].p, aPo2D[ aTri[itri0].v[inotri0] ].p, aPo2D[ipoi1].p ) > 1.0e-20 );
						assert( TriArea( aPo2D[ipoi0].p, aPo2D[ipoi1].p, aPo2D[ aTri[itri0].v[inotri1] ].p ) > 1.0e-20 );
//						std::cout << ratio << std::endl;
						if( ratio < 1.0e-20 ){
	//						std::cout << "未実装 辺上に点がある場合" << std::endl;
							assert(0);
							return false;
						}
						else if( ratio > 1.0 - 1.0e-10 ){
	//						std::cout << "未実装 辺上に点がある場合" << std::endl;
							assert(0);
							return false;
						}
						else{
							const unsigned int ied0 = 3 - inotri0 - inotri1;
//							std::cout << aTri[itri0].g2[ied0] << std::endl;
							assert( aTri[itri0].g2[ied0] == -2 );
							const unsigned int itri1 = aTri[itri0].s2[ied0];
							const unsigned int ied1 = relTriTri[ aTri[itri0].r2[ied0] ][ied0];
							assert( aTri[itri1].s2[ied1] == itri0 );
							assert( aTri[itri1].g2[ied1] == -2 );
							FlipEdge(itri0,ied0,aPo2D,aTri);
							continue;
						}
					}
				}
			}
		}
	}

	{	// 外側の三角形の除去して結果をaTri_inに代入
		aTri_in.clear();	// 初期化する
		// 内側にある三角形をひとつ(itri0_ker)見つける
		unsigned int itri0_ker = aTri.size();
		{
			const unsigned int ipo0=aVtxInd[0], ipo1=aVtxInd[1];
			unsigned int itri0,inotri0,inotri1;
			bool res = FindEdge(ipo0, ipo1,
				itri0, inotri0, inotri1,
				aPo2D, aTri );
			assert( res == true );
			assert( aTri[itri0].v[inotri0] == aVtxInd[0] );
			assert( aTri[itri0].v[inotri1] == aVtxInd[1] );
			itri0_ker = itri0;
		}
		assert( itri0_ker < aTri.size() );

		// 領域の外の要素ならフラグが-1、そうでなければフラグは昇順の要素番号が入った配列inout_flgを作る
		unsigned int ntri_in;
		std::vector<int> inout_flg;	// フラグ配列
		{	// 上で見つけた内側の三角形を核として内側の三角形を周囲に拡大していく
			inout_flg.resize(aTri.size(),-1);
			inout_flg[itri0_ker] = 0;
			ntri_in = 1;
			std::stack<unsigned int> ind_stack;	// 周囲が探索されていない三角形
			ind_stack.push(itri0_ker);
			for(;;){
				if( ind_stack.empty() ) break;
				const unsigned int itri_cur = ind_stack.top();
				ind_stack.pop();
				for(unsigned int inotri=0;inotri<3;inotri++){
					if( aTri[itri_cur].g2[inotri] != -2 ) continue;
					const unsigned int itri_s = aTri[itri_cur].s2[inotri];
					if( inout_flg[itri_s] == -1 ){
						inout_flg[itri_s] = ntri_in;
						ntri_in++;
						ind_stack.push(itri_s);
					}
				}
			}
		}

		// フラグ配列に沿って内側の三角形を集めた配列aTri_inを作る
		aTri_in.resize( ntri_in );
		for(unsigned int itri=0;itri<aTri.size();itri++){
			if( inout_flg[itri] != -1 ){
				int itri_in = inout_flg[itri];
				assert( itri_in >= 0 && (unsigned int)itri_in < ntri_in );
				aTri_in[itri_in] = aTri[itri];
			}
		}
		// 内側の三角形配列のの隣接情報を作る
		for(unsigned int itri=0;itri<aTri_in.size();itri++){
			for(unsigned int ifatri=0;ifatri<3;ifatri++){
				if( aTri_in[itri].g2[ifatri] == -3 ){
					aTri_in[itri].g2[ifatri] = -1;
				}
				if( aTri_in[itri].g2[ifatri] != -2 ) continue;
				int itri_s0 = aTri_in[itri].s2[ifatri];
				assert( itri_s0 >= 0 && (unsigned int)itri_s0 < aTri.size() );
				int itri_in_s0 = inout_flg[itri_s0];
				assert( itri_in_s0 >= 0 && (unsigned int)itri_in_s0 < aTri_in.size() );
				aTri_in[itri].s2[ifatri] = itri_in_s0;
			}
		}
	}
	////////////////
	std::vector<int> map_po_del;
	unsigned int npo_pos;
	{
		map_po_del.resize( aPo2D.size(), -1 );
		for(unsigned int ipo=0;ipo<aPoDel.size();ipo++){
			map_po_del[ aPoDel[ipo] ] = -2;
		}
		npo_pos = 0;
		for(unsigned int ipo=0;ipo<aPo2D.size();ipo++){
			if( map_po_del[ipo] == -2 ) continue;
			map_po_del[ipo] = npo_pos;
			npo_pos++;
		}
	}
	{
		std::vector<CPoint2D> aPo_tmp = aPo2D;
		aPo2D.resize( npo_pos );
		for(unsigned int ipo=0;ipo<map_po_del.size();ipo++){
			if( map_po_del[ipo] == -2 ) continue;
			unsigned int ipo1 = map_po_del[ipo];
			aPo2D[ipo1] = aPo_tmp[ipo];
		}
	}
	for(unsigned int itri=0;itri<aTri_in.size();itri++){
		for(unsigned int ifatri=0;ifatri<3;ifatri++){
			assert( aTri_in[itri].v[ifatri] != -2 );
			const unsigned int ipo = aTri_in[itri].v[ifatri];
			aTri_in[itri].v[ifatri] = map_po_del[ipo];
			aPo2D[ipo].e = itri;
			aPo2D[ipo].d = ifatri;
		}
	}
	return true;
}

void LaplacianSmoothing( std::vector<CPoint2D>& aPo, const std::vector<STri2D>& aTri,
	const std::vector<unsigned int>& aflg_isnt_move)
{	
	for(unsigned int ipoin=0;ipoin<aPo.size();ipoin++){	// 点周りの点を探索して調べる。
		if( ipoin < aflg_isnt_move.size() ){
			if( aflg_isnt_move[ipoin] == 1 ) continue;
		}
		const unsigned int itri_ini = aPo[ipoin].e;
		const unsigned int inoel_c_ini = aPo[ipoin].d;
		assert( itri_ini < aTri.size() );
		assert( inoel_c_ini < 3 );
		assert( aTri[itri_ini].v[inoel_c_ini] == ipoin );
		unsigned int itri0= itri_ini;
		unsigned int inoel_c0 = inoel_c_ini;
		unsigned int inoel_b0 = noelTriEdge[inoel_c0][0];
		bool is_bound_flg = false;
		CVector2D vec_delta = aPo[ipoin].p;
		unsigned int ntri_around = 1;
		for(;;){
			assert( itri0 < aTri.size() );
			assert( inoel_c0 < 3 );
			assert( aTri[itri0].v[inoel_c0] == ipoin );
			{
				vec_delta.x += aPo[ aTri[itri0].v[inoel_b0] ].p.x;
				vec_delta.y += aPo[ aTri[itri0].v[inoel_b0] ].p.y;
				ntri_around++;
			}
			if( aTri[itri0].g2[inoel_b0] == -2 ){
				unsigned int itri1 = aTri[itri0].s2[inoel_b0];
				const char rel01 = aTri[itri0].r2[inoel_b0];
				unsigned int inoel_c1 = relTriTri[rel01][inoel_c0];
				unsigned int inoel_b1 = relTriTri[rel01][ noelTriEdge[inoel_c0][1] ];
				assert( itri1 < aTri.size() );
				assert( aTri[itri1].s2[ relTriTri[rel01][inoel_b0] ] == itri0 );
				assert( aTri[itri1].v[inoel_c1] == ipoin );
				if( itri1 == itri_ini ) break;
				itri0 = itri1;
				inoel_c0 = inoel_c1;
				inoel_b0 = inoel_b1;
			}
			else{	// この点は境界上の点だから動かしてはならない。
				is_bound_flg = true;
				break;
			}
		}
		if( is_bound_flg ) continue;		
		aPo[ipoin].p.x = vec_delta.x / ntri_around;
		aPo[ipoin].p.y = vec_delta.y / ntri_around;
	}
}


void LaplaceDelaunaySmoothing( std::vector<CPoint2D>& aPo, std::vector<STri2D>& aTri,
	const std::vector<unsigned int>& aflg_isnt_move )
{
	for(unsigned int ipoin=0;ipoin<aPo.size();ipoin++){	// 点周りの点を探索して調べる。
		if( ipoin < aflg_isnt_move.size() ){
			if( aflg_isnt_move[ipoin] == 1 ) continue;
		}
		const unsigned int itri_ini = aPo[ipoin].e;
		const unsigned int inoel_c_ini = aPo[ipoin].d;
		assert( itri_ini < aTri.size() );
		assert( inoel_c_ini < 3 );
		assert( aTri[itri_ini].v[inoel_c_ini] == ipoin );
		unsigned int itri0= itri_ini;
		unsigned int inoel_c0 = inoel_c_ini;
		unsigned int inoel_b0 = noelTriEdge[inoel_c0][0];
		bool is_bound_flg = false;
		CVector2D vec_delta = aPo[ipoin].p;
		unsigned int ntri_around = 1;
		for(;;){	// 点の周りの要素を一回りする
			assert( itri0 < aTri.size() );
			assert( inoel_c0 < 3 );
			assert( aTri[itri0].v[inoel_c0] == ipoin );
			{
				vec_delta.x += aPo[ aTri[itri0].v[inoel_b0] ].p.x;
				vec_delta.y += aPo[ aTri[itri0].v[inoel_b0] ].p.y;
				ntri_around++;
			}
			if( aTri[itri0].g2[inoel_b0] == -2 ){
				unsigned int itri1 = aTri[itri0].s2[inoel_b0];
				const char rel01 = aTri[itri0].r2[inoel_b0];
				unsigned int inoel_c1 = relTriTri[rel01][inoel_c0];
				unsigned int inoel_b1 = relTriTri[rel01][ noelTriEdge[inoel_c0][1] ];
				assert( itri1 < aTri.size() );
				assert( aTri[itri1].s2[ relTriTri[rel01][inoel_b0] ] == itri0 );
				assert( aTri[itri1].v[inoel_c1] == ipoin );
				if( itri1 == itri_ini ) break;
				itri0 = itri1;
				inoel_c0 = inoel_c1;
				inoel_b0 = inoel_b1;
			}
			else{	// この点は境界上の点だから動かしてはならない。
				is_bound_flg = true;
				break;
			}
		}
		if( is_bound_flg ) continue;		
		aPo[ipoin].p.x = vec_delta.x / ntri_around;
		aPo[ipoin].p.y = vec_delta.y / ntri_around;
		DelaunayAroundPoint(ipoin,aPo,aTri);
	}
}

void MeshingInside(std::vector<CPoint2D>& aPo2D, std::vector<STri2D>& aTri, 
				   const std::vector<unsigned int>& aVtxInd, 
				   const double len)
{
	std::vector<unsigned int> aflag_isnt_move;	// フラグが１なら動かさない
	{
		aflag_isnt_move.resize( aPo2D.size(), 0 );
		for(unsigned int iver=0;iver<aVtxInd.size();iver++){
			const unsigned int ivec = aVtxInd[iver];
			aflag_isnt_move[ivec] = 1;
		}
	}

	{	// aTriに節点を追加
		double ratio = 3.0;
		for(;;){
			unsigned int nadd = 0;
			for(unsigned int itri=0;itri<aTri.size();itri++){
				const double area = TriArea(
					aPo2D[aTri[itri].v[0]].p, 
					aPo2D[aTri[itri].v[1]].p, 
					aPo2D[aTri[itri].v[2]].p);
				if( area > len * len * ratio ){
					// itriの重心に新しい節点を追加
					const unsigned int ipo0 = aPo2D.size();	// ipo0は新しい節点番号
					aPo2D.resize( aPo2D.size()+1 );
					aPo2D[ipo0].p.x = (aPo2D[aTri[itri].v[0]].p.x+aPo2D[aTri[itri].v[1]].p.x+aPo2D[aTri[itri].v[2]].p.x)/3.0;
					aPo2D[ipo0].p.y = (aPo2D[aTri[itri].v[0]].p.y+aPo2D[aTri[itri].v[1]].p.y+aPo2D[aTri[itri].v[2]].p.y)/3.0;
					InsertPoint_Elem(ipo0,itri,aPo2D,aTri);
					DelaunayAroundPoint(ipo0,aPo2D,aTri);
					nadd++;
				}
			}
			LaplacianSmoothing(aPo2D,aTri,aflag_isnt_move);
//			LaplaceDelaunaySmoothing(aPo2D,aTri);
			if( nadd != 0 ){ ratio *= 0.8; }
			else{ ratio *= 0.5; }
			if( ratio < 0.65 ) break;
		}
	}

	LaplaceDelaunaySmoothing(aPo2D,aTri,aflag_isnt_move);
}

inline double TriArea2D(const double v1[], const double v2[], const double v3[]){
	return 0.5*( (v2[0]-v1[0])*(v3[1]-v1[1]) - (v3[0]-v1[0])*(v2[1]-v1[1]) );
}

bool IsCrossLines(const double po_s0[], const double po_e0[],
	const double po_s1[], const double po_e1[] )
{
	const double area1 = TriArea2D(po_s0,po_e0,po_s1);
	const double area2 = TriArea2D(po_s0,po_e0,po_e1);
	if( area1 * area2 > 0.0 ) return false;
	const double area3 = TriArea2D(po_s1,po_e1,po_s0);
	const double area4 = TriArea2D(po_s1,po_e1,po_e0);
	if( area3 * area4 > 0.0 ) return false;
	return true;
}

bool IsInclude_Loop(
		const double co[],
		const unsigned int nxys, const double* xys)
{
	unsigned int inum_cross = 0;
	for(unsigned int itr=0;itr<10;itr++){
		const double dir[2] = { cos(1.0*(itr+1)*23), sin(1.0*(itr+1)*23) };
		const double codir[2] = { co[0]+dir[0], co[1]+dir[1] };
		bool is_fail = false;	// 判定不能の場合はフラグが立つ
		inum_cross = 0;
		for(unsigned int ixys=0;ixys<nxys;ixys++){
			const unsigned int ipo0 = ixys;
			unsigned int ipo1 = ixys+1;
			if( ipo1 == nxys ){ ipo1 = 0; }
			const double area0 = TriArea2D(co,codir,xys+ipo0*2);
			const double area1 = TriArea2D(co,xys+ipo1*2,codir);
			double r1 =  area0 / (area0+area1);
			double r0 =  area1 / (area0+area1);
			if( fabs(area0+area1) < 1.0e-20 ){
				is_fail = true;
				break;
			}
			if( fabs(r0) < 1.0e-3 || fabs(r1) < 1.0e-3 ){
				is_fail = true;
				break;
			}
			if( r0*r1 < 0 ){ continue; }
			double po2[2] = { r0*xys[ipo0*2  ]+r1*xys[ipo1*2  ],  r0*xys[ipo0*2+1]+r1*xys[ipo1*2+1] };
			const double area2 = TriArea2D(co,codir,po2);
			double d = (po2[0]-co[0])*dir[0] + (po2[1]-co[1])*dir[1];
			if( d > 0 ) inum_cross++;
		}
		if( is_fail ){ continue; }
		if( inum_cross%2 == 0 ){ return false; }
		else if( inum_cross%2 == 1 ){ return true; }
	}
	return false;
}


////////////////////////////////
int delaunay_triangulation(
	// input
	int nloop,	// number of loops (at least one = outmost loop)  e.g 2 (dounat)
	int* nxys,	// numbers of vertices in each loop (first one is the outmost loop)  e.g. {4, 3} (triangle inside square)
	double* xys,	// xyz values  e.g. {0,0,  0,1, 1,1, 1,0,      0.1,0.1, 0,5,0.9, 0.9,0.1} 
	double max_edge_length, // maximum edge length (if minus then no adding) e.g. {-1}

	// output
	int* ntri,		// number of triangles e.g. 2
	int** atri,		// list of triangles (each triangle is a list of three indices) e.g. {0,1,2, 0,2,3}
	int* nxys_new,	// number of newly added vertices e.g. 1
	double** xys_new)	// coordinates of newly added vertices e.g. {0,0.5,0}
{
	// 辺の頂点へのポインタaIndXYsを作る
	std::vector<unsigned int> aIndXYs;
	aIndXYs.resize(nloop+1);
	aIndXYs[0] = 0;
	for(unsigned int iloop=0;iloop<nloop;iloop++){
		aIndXYs[iloop+1] = aIndXYs[iloop] + nxys[iloop];
	}
	////////////////////////////////
	// enter Input check section

	{ // ループ上の辺の数が３より小さくないかチェック
		for(unsigned int iloop=0;iloop<nloop;iloop++){
			if( nxys[iloop] < 3 ) return 0;
		}
	}
	{ 	// 親ループの中に子ループが入っているか？子ループが他の子ループに含まれていないか
		////////////////////////////////
		// 親ループの中に子ループが入っているか？
		for(unsigned int iloop=1;iloop<nloop;iloop++){
			for(unsigned int ipo=aIndXYs[iloop];ipo<aIndXYs[iloop+1];ipo++){
				if( !IsInclude_Loop(xys+ipo*2, aIndXYs[1]-aIndXYs[0],xys) ) return 0;
			}
		}
		// 子ループが他の子ループに含まれていないか
		for(unsigned int iloop=1;iloop<nloop;iloop++){
			for(unsigned int jloop=0;jloop<nloop;jloop++){
				if( iloop == jloop ) continue;
				for(unsigned int jpo=aIndXYs[jloop];jpo<aIndXYs[jloop+1];jpo++){
					if( IsInclude_Loop(xys+jpo*2, 
						aIndXYs[iloop+1]-aIndXYs[iloop],xys+aIndXYs[iloop]*2) ) return 0;
				}
			}
		}
	}
	{ // 自己交錯をしないかをチェック
		bool is_intersect = false;
		for(unsigned int iloop=0;iloop<nloop;iloop++){
			const unsigned int nbar_i = nxys[iloop];
			for(unsigned int ibar=0;ibar<nbar_i;ibar++){
				const unsigned int ipo0 = aIndXYs[iloop] + ibar;
				unsigned int ipo1 = aIndXYs[iloop] + ibar+1;
				if( ibar == nbar_i-1 ){ ipo1 = aIndXYs[iloop]; }
				const double xmax_i = ( xys[ipo1*2+0] > xys[ipo0*2+0] ) ? xys[ipo1*2+0] : xys[ipo0*2+0];
				const double xmin_i = ( xys[ipo1*2+0] < xys[ipo0*2+0] ) ? xys[ipo1*2+0] : xys[ipo0*2+0];
				const double ymax_i = ( xys[ipo1*2+1] > xys[ipo0*2+1] ) ? xys[ipo1*2+1] : xys[ipo0*2+1];
				const double ymin_i = ( xys[ipo1*2+1] < xys[ipo0*2+1] ) ? xys[ipo1*2+1] : xys[ipo0*2+1];
				for(unsigned int jbar=ibar+2;jbar<nbar_i;jbar++){
					const unsigned int jpo0 = aIndXYs[iloop] + jbar;
					unsigned int jpo1 = aIndXYs[iloop] + jbar+1;
					if( jbar == nbar_i-1 ){
						if( ibar == 0 ) continue;
						jpo1 = aIndXYs[iloop];
					}
					const double xmax_j = ( xys[jpo1*2+0] > xys[jpo0*2+0] ) ? xys[jpo1*2+0] : xys[jpo0*2+0];
					const double xmin_j = ( xys[jpo1*2+0] < xys[jpo0*2+0] ) ? xys[jpo1*2+0] : xys[jpo0*2+0];
					const double ymax_j = ( xys[jpo1*2+1] > xys[jpo0*2+1] ) ? xys[jpo1*2+1] : xys[jpo0*2+1];
					const double ymin_j = ( xys[jpo1*2+1] < xys[jpo0*2+1] ) ? xys[jpo1*2+1] : xys[jpo0*2+1];
					if( xmin_j > xmax_i || xmax_j < xmin_i ) continue;	// 交錯がありえないパターンを除外
					if( ymin_j > ymax_i || ymax_j < ymin_i ) continue;	// 上に同じ
					if( IsCrossLines(xys+ipo0*2,xys+ipo1*2,  xys+jpo0*2,xys+jpo1*2) ){
						is_intersect = true;
						break;
					}
				}
				if( is_intersect ) break;
				for(unsigned int jloop=iloop+1;jloop<nloop;jloop++){
					const unsigned int nbar_j = nxys[jloop];
					for(unsigned int jbar=0;jbar<nbar_j;jbar++){
						const unsigned int jpo0 = aIndXYs[jloop] + jbar;
						unsigned int jpo1 = aIndXYs[jloop] + jbar+1;
						if( jbar == nbar_j-1 ){ jpo1 = aIndXYs[jloop]; }
						const double xmax_j = ( xys[jpo1*2+0] > xys[jpo0*2+0] ) ? xys[jpo1*2+0] : xys[jpo0*2+0];
						const double xmin_j = ( xys[jpo1*2+0] < xys[jpo0*2+0] ) ? xys[jpo1*2+0] : xys[jpo0*2+0];
						const double ymax_j = ( xys[jpo1*2+1] > xys[jpo0*2+1] ) ? xys[jpo1*2+1] : xys[jpo0*2+1];
						const double ymin_j = ( xys[jpo1*2+1] < xys[jpo0*2+1] ) ? xys[jpo1*2+1] : xys[jpo0*2+1];
						if( xmin_j > xmax_i || xmax_j < xmin_i ) continue;	// 交錯がありえないパターンを除外
						if( ymin_j > ymax_i || ymax_j < ymin_i ) continue;	// 上に同じ
						if( IsCrossLines(xys+ipo0*2,xys+ipo1*2,  xys+jpo0*2,xys+jpo1*2) ){
							is_intersect = true;
							break;
						}
					}
					if( is_intersect ) break;
				}
				if( is_intersect ) break;
			}
			if( is_intersect ) break;
		}
		if( is_intersect ) return 0;
	}

	// end of input check section
	////////////////////////////////////////////////

	unsigned int nxys_presum = aIndXYs[nloop];
	std::vector<CPoint2D> aPo2D;
	aPo2D.resize(nxys_presum);
	for(unsigned int ixys=0;ixys<nxys_presum;ixys++){
		aPo2D[ixys].p.x = xys[ixys*2+0];
		aPo2D[ixys].p.y = xys[ixys*2+1];
		aPo2D[ixys].e = -1;
		aPo2D[ixys].d = -1;
	}
	////////////////////////////////
	// 辺に加える点を作る
	std::vector< std::vector<unsigned int> > aPoInEd;
	aPoInEd.resize(nxys_presum);
	if( max_edge_length > 0 ){
		for(unsigned int iloop=0;iloop<nloop;iloop++){
			unsigned int nadd = 0;
			const unsigned int nbar = nxys[iloop];
			for(unsigned int ibar=0;ibar<nbar;ibar++){
				unsigned int ipo0 = aIndXYs[iloop]+ibar;
				unsigned int ipo1 = aIndXYs[iloop]+ibar+1;
				if( ibar == nbar-1 ){ ipo1 = aIndXYs[iloop]; }
				const double len = Length( aPo2D[ipo0].p, aPo2D[ipo1].p );
				nadd = (int)(len / max_edge_length);
				if( nadd == 0 ) continue;
				const unsigned int ndiv = nadd+1;
				const double delx = (aPo2D[ipo1].p.x - aPo2D[ipo0].p.x)/ndiv;
				const double dely = (aPo2D[ipo1].p.y - aPo2D[ipo0].p.y)/ndiv;
				for(unsigned int iadd=0;iadd<nadd;iadd++){
					const unsigned int ipo = aPo2D.size();
					CPoint2D po;
					po.p.x = aPo2D[ipo0].p.x + delx*(iadd+1);
					po.p.y = aPo2D[ipo0].p.y + dely*(iadd+1);
					po.e = -1;
					po.d = -1;
					aPo2D.push_back(po);
					aPoInEd[ aIndXYs[iloop]+ibar ].push_back(ipo);
				}
			}
		}
	}
	////////////////////////////////
	// ループ上の点のIndexを作る
	std::vector<unsigned int> aPtrVtxInd, aVtxInd;
	{
		aPtrVtxInd.resize(nloop+1);
		aPtrVtxInd[0] = 0;
		for(unsigned int iloop=0;iloop<nloop;iloop++){
			const unsigned int nbar0 = nxys[iloop];	// 入力のループの中の点の辺の数
			unsigned int nbar1 = nbar0;	// 辺上に点を追加した後のループの中の辺の数
			for(unsigned int ibar=0;ibar<nbar0;ibar++){
				nbar1 += aPoInEd[ aIndXYs[iloop]+ibar].size();
			}
			aPtrVtxInd[iloop+1] = aPtrVtxInd[iloop] + nbar1;
		}
		////////////////
		aVtxInd.resize(aPtrVtxInd[nloop]);
		{
			unsigned int ivtx0 = 0;
			for(unsigned int iloop=0;iloop<nloop;iloop++){
				double area_loop = 0;
				{
					CVector2D vtmp(0,0);
					const unsigned int nbar = aPtrVtxInd[iloop+1]-aPtrVtxInd[iloop];
					for(unsigned int ibar=0;ibar<nbar;ibar++){
						unsigned int ipo0 = aPtrVtxInd[iloop]+ibar;
						unsigned int ipo1 = aPtrVtxInd[iloop]+ibar+1;
						if( ibar == nbar-1 ){ ipo1 = aPtrVtxInd[iloop]; }
						area_loop += TriArea( vtmp, aPo2D[ipo0].p, aPo2D[ipo1].p );
					}
				}
				const unsigned int nbar0 = nxys[iloop];	// 入力のループの中の点の辺の数
				if( (area_loop > 0) == (iloop == 0) ){	// 左回り
					for(unsigned int ibar=0;ibar<nbar0;ibar++){
						unsigned int ie = aIndXYs[iloop] + ibar;
						const std::vector<unsigned int>& add = aPoInEd[ie];
						aVtxInd[ivtx0] = ie;	ivtx0++;
						for(unsigned int iadd=0;iadd<add.size();iadd++){
							aVtxInd[ivtx0] = add[iadd];	ivtx0++;
						}
					}
				}
				else{	// 右回り
					for(unsigned int ibar=0;ibar<nbar0;ibar++){
						unsigned int ie = aIndXYs[iloop] + nxys[iloop] - 1 - ibar;
						const std::vector<unsigned int>& add = aPoInEd[ie];
						const unsigned int nadd = add.size();
						for(unsigned int iadd=0;iadd<nadd;iadd++){
							aVtxInd[ivtx0] = add[nadd-1-iadd];	ivtx0++;
						}
						aVtxInd[ivtx0] = ie;	ivtx0++;
					}
				}
			}
		}
	}
	////////////////////////////////
	// ここから三角分割ルーティン
	std::vector<STri2D> aTri_in;
	if( !TriangulateOuterLoop(aPo2D,aTri_in,    aPtrVtxInd, aVtxInd) ){
		return 0;
	}
	if( max_edge_length > 0 ){
		MeshingInside(aPo2D,aTri_in, aVtxInd,max_edge_length);
	}

	////////////////
	// ここから出力を作る
	(*ntri) = aTri_in.size();
	(*atri) = new int [(*ntri)*3];
	for(unsigned int itri=0;itri<(*ntri);itri++){
		(*atri)[itri*3+0] = aTri_in[itri].v[0];
		(*atri)[itri*3+1] = aTri_in[itri].v[1];
		(*atri)[itri*3+2] = aTri_in[itri].v[2];
	}
	(*nxys_new) = aPo2D.size()-nxys_presum;
	(*xys_new) = new double [(*nxys_new)*2];
	for(unsigned int ixys=0;ixys<(*nxys_new);ixys++){
		(*xys_new)[ixys*2+0] = aPo2D[nxys_presum+ixys].p.x;
		(*xys_new)[ixys*2+1] = aPo2D[nxys_presum+ixys].p.y;
	}

	return 1;
}

void test_dtri()
{
	int vtx_num = 4;
	int * nxys = new int[ 1 ];
	double xys[ 8 ] = { 0, 0, 1, 0, 1, 1, 0, 1 };

	nxys[ 0 ] = vtx_num;

	// output
	int ntri;
	int* atri; 
	int nxys_new;
	double* xys_new;

	if( delaunay_triangulation( 1, nxys, xys, 0.5, &ntri, &atri, &nxys_new, &xys_new ) ){
		cout << ntri << endl;
	}
	else{
		cout << "fail" << endl;
	}
}