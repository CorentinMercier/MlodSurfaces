// --------------------------------------------------------------------------
// Source code provided FOR REVIEW ONLY, as part of the submission entitled
// "Moving Level-of-Detail Surfaces".
//
// A proper version of this code will be released if the paper is accepted
// with the proper licence, documentation and bug fix.
// Currently, this material has to be considered confidential and shall not
// be used outside the review process.
//
// All right reserved. The Authors
// --------------------------------------------------------------------------

#pragma once

#include <algorithm>
#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>

using Vec3 = glm::vec3;
using Mat3 = glm::mat3;

template<typename Pt>
inline Mat3 tensor(const Pt& p1, const Pt& p2)
{
	return glm::outerProduct(p1, p2);
}

inline Vec3 operator/(const Vec3& v, double x) { return v / float(x); }
inline Vec3 operator*(double x, const Vec3& v) { return v * float(x); }

struct BBOX2
{
	Vec3 bmin,bmax;

	BBOX2() { clear(); }
	BBOX2(Vec3 const& min, Vec3 const& max) : bmin(min), bmax(max) {}

	void clear()
	{
		bmin[0] = std::numeric_limits<float>::infinity();
		bmax[0] = -std::numeric_limits<float>::infinity();
	}
	bool isCleared() const { return bmin[0] <= bmax[0]; }

	Vec3 center() const { return (bmin + bmax) / 2; }

	template<class Pt>
	void set(const Pt& p)
	{
		bmin = Vec3(p[0],p[1],p[2]);
		bmax = Vec3(p[0],p[1],p[2]);
	}

	template<class Pt>
	void set(const Pt& min, const Pt& max )
	{
		bmin = Vec3(min[0],min[1],min[2]);
		bmax = Vec3(max[0],max[1],max[2]);
	}

	template<class Pt>
	void add(const Pt& p)
	{
		bmin = glm::min(bmin,p);
		bmax = glm::max(bmax,p);
	}

	void add(const BBOX2& b)
	{
		bmin = glm::min(bmin,b.bmin);
		bmax = glm::max(bmax,b.bmax);
	}

	float squareDiagonal() const { return length2(bmax - bmin); }
	float diagonal() const { return sqrt(squareDiagonal()); }
	float radius() const { return diagonal() / 2.0; }
	float squareRadius() const { return squareDiagonal() / 4.0; }

	char getLargestExtent() const
	{
		if(  bmax[0] - bmin[0]  >  bmax[1] - bmin[1]  )
		{
			if( bmax[0] - bmin[0] > bmax[2] - bmin[2] )
				return 0;
			return 2;
		}
		else
		{
			if( bmax[1] - bmin[1] > bmax[2] - bmin[2] )
				return 1;
			return 2;
		}
	}
	float getExtentValue(char i) const { return bmax[i] - bmin[i]; }
	float getLargestExtentValue() const { return getExtentValue( getLargestExtent() ); }

	template<class Pt>
	float getPseudoExtentInDirection(Pt const& dir) const
	{
		return getExtentValue(0) * fabs(dir[0])
			+ getExtentValue(1) * fabs(dir[1])
			+ getExtentValue(2) * fabs(dir[2]);
	}

	void splitAlongAxis( char axis , float value , BBOX2 & box1 , BBOX2 & box2 )
	{
		// for safety:
		value = std::max<float>( std::min<float>( value , bmax[axis] ) , bmin[axis] );
		Vec3 bmax1 = bmax;
		bmax1[axis] = value;
		Vec3 bmin2 = bmin;
		bmin2[axis] = value;
		box1.set( bmin , bmax1 );
		box2.set( bmin2 , bmax );
	}

	Vec3 getCorner0() const { return bmin; }
	Vec3 getCorner1() const { return Vec3( bmax[0] , bmin[1] , bmin[2] ); }
	Vec3 getCorner2() const { return Vec3( bmin[0] , bmax[1] , bmin[2] ); }
	Vec3 getCorner3() const { return Vec3( bmax[0] , bmax[1] , bmin[2] ); }
	Vec3 getCorner4() const { return Vec3( bmin[0] , bmin[1] , bmax[2] ); }
	Vec3 getCorner5() const { return Vec3( bmax[0] , bmin[1] , bmax[2] ); }
	Vec3 getCorner6() const { return Vec3( bmin[0] , bmax[1] , bmax[2] ); }
	Vec3 getCorner7() const { return bmax; }
};
