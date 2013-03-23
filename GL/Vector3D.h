/*
 *  Vector3D.h
 *
 *  Created by Yusuf Yilmaz on 4/19/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *	Ripped of from DirectX.
 */
#ifndef _VECTOR3D_H
#define _VECTOR3D_H

#include <math.h>

class Vector3D
{
public:
    Vector3D() {};
    Vector3D( const double * pf)
	{
		x = pf[0];
		y = pf[1];
		z = pf[2];
	}

    Vector3D( const Vector3D& v)
	{
		x = v.x;
		y = v.y;
		z = v.z;
	}

    Vector3D( double fx, double fy, double fz )
	{
		x = fx;
		y = fy;
		z = fz;
	}

    // assignment operators
    Vector3D& operator += ( const Vector3D& v)
	{
		x += v.x;
		y += v.y;
		z += v.z;
		return *this;
	}

    Vector3D& operator -= ( const Vector3D& v)
	{
		x -= v.x;
		y -= v.y;
		z -= v.z;
		return *this;
	}

    Vector3D& operator *= ( double f)
	{
		x *= f;
		y *= f;
		z *= f;
		return *this;
	}

    Vector3D& operator /= ( double f)
	{
		double fInv = 1.0 / f;
		x *= fInv;
		y *= fInv;
		z *= fInv;
		return *this;
	}

    // unary operators
    Vector3D operator + () const
	{
		return *this;
	}

    Vector3D operator - () const
	{
		return Vector3D(-x, -y, -z);
	}

    // binary operators
    Vector3D operator + ( const Vector3D& v) const
	{
		return Vector3D(x + v.x, y + v.y, z + v.z);
	}

    Vector3D operator - ( const Vector3D& v) const
	{
		return Vector3D(x - v.x, y - v.y, z - v.z);
	}

    Vector3D operator * ( double f) const
	{
		return Vector3D(x * f, y * f, z * f);
	}

    Vector3D operator / ( double f) const
	{
		double fInv = 1.0 / f;
		return Vector3D(x * fInv, y * fInv, z * fInv);
	}

	bool operator == ( const Vector3D& v) const
	{
		return x == v.x && y == v.y && z == v.z;
	}

    bool operator != ( const Vector3D& v) const
	{
		return x != v.x || y != v.y || z != v.z;
	}

	double &operator[](int i) {
		switch (i) {
			case 0: return x;
			case 1: return y;
			case 2: return z;
		}
		return z;
	}

	Vector3D Normalize()
	{

		double l = length();
		//if(l != 0)
		//{
			x /= l;
			y /= l;
			z /= l;
		/*}
		else {

			x = 1;
			y = 0;
			z = 0;
		}*/

		return *this;
	}

	Vector3D Cross(const Vector3D &v) const{
		return Vector3D(y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y * v.x);
	}

	double Dot(const Vector3D &v) const{
		return x*v.x + y*v.y + z*v.z;
	}

	double lengthSquared() const
	{
		return x*x + y*y + z*z;
	}

	double length() const
	{
		return sqrt(lengthSquared());
	}

	Vector3D clamp01() {
		if (x>1.f)
			x=1.f;
		else if (x<0.f)
			x=0.f;

		if (y>1.f)
			y=1.f;
		else if (y<0.f)
			y=0.f;

		if (z>1.f)
			z=1.f;
		else if (z<0.f)
			z=0.f;

		return *this;
	}

	// Other DirectX functions
	/*

	 // Hermite interpolation between position V1, tangent T1 (when s == 0)
	 // and position V2, tangent T2 (when s == 1).
	 D3DXVECTOR3* WINAPI D3DXVec3Hermite
	 ( D3DXVECTOR3 *pOut, CONST D3DXVECTOR3 *pV1, CONST D3DXVECTOR3 *pT1,
	 CONST D3DXVECTOR3 *pV2, CONST D3DXVECTOR3 *pT2, FLOAT s );

	 // CatmullRom interpolation between V1 (when s == 0) and V2 (when s == 1)
	 D3DXVECTOR3* WINAPI D3DXVec3CatmullRom
	 ( D3DXVECTOR3 *pOut, CONST D3DXVECTOR3 *pV0, CONST D3DXVECTOR3 *pV1,
	 CONST D3DXVECTOR3 *pV2, CONST D3DXVECTOR3 *pV3, FLOAT s );

	 // Barycentric coordinates.  V1 + f(V2-V1) + g(V3-V1)
	 D3DXVECTOR3* WINAPI D3DXVec3BaryCentric
	 ( D3DXVECTOR3 *pOut, CONST D3DXVECTOR3 *pV1, CONST D3DXVECTOR3 *pV2,
	 CONST D3DXVECTOR3 *pV3, FLOAT f, FLOAT g);

	 // Transform (x, y, z, 1) by matrix.
	 D3DXVECTOR4* WINAPI D3DXVec3Transform
	 ( D3DXVECTOR4 *pOut, CONST D3DXVECTOR3 *pV, CONST D3DXMATRIX *pM );

	 // Transform (x, y, z, 1) by matrix, project result back into w=1.
	 D3DXVECTOR3* WINAPI D3DXVec3TransformCoord
	 ( D3DXVECTOR3 *pOut, CONST D3DXVECTOR3 *pV, CONST D3DXMATRIX *pM );

	 // Transform (x, y, z, 0) by matrix.  If you transforming a normal by a
	 // non-affine matrix, the matrix you pass to this function should be the
	 // transpose of the inverse of the matrix you would use to transform a coord.
	 D3DXVECTOR3* WINAPI D3DXVec3TransformNormal
	 ( D3DXVECTOR3 *pOut, CONST D3DXVECTOR3 *pV, CONST D3DXMATRIX *pM );


	 // Transform Array (x, y, z, 1) by matrix.
	 D3DXVECTOR4* WINAPI D3DXVec3TransformArray
	 ( D3DXVECTOR4 *pOut, UINT OutStride, CONST D3DXVECTOR3 *pV, UINT VStride, CONST D3DXMATRIX *pM, UINT n );

	 // Transform Array (x, y, z, 1) by matrix, project result back into w=1.
	 D3DXVECTOR3* WINAPI D3DXVec3TransformCoordArray
	 ( D3DXVECTOR3 *pOut, UINT OutStride, CONST D3DXVECTOR3 *pV, UINT VStride, CONST D3DXMATRIX *pM, UINT n );

	 // Transform (x, y, z, 0) by matrix.  If you transforming a normal by a
	 // non-affine matrix, the matrix you pass to this function should be the
	 // transpose of the inverse of the matrix you would use to transform a coord.
	 D3DXVECTOR3* WINAPI D3DXVec3TransformNormalArray
	 ( D3DXVECTOR3 *pOut, UINT OutStride, CONST D3DXVECTOR3 *pV, UINT VStride, CONST D3DXMATRIX *pM, UINT n );

	 // Project vector from object space into screen space
	 D3DXVECTOR3* WINAPI D3DXVec3Project
	 ( D3DXVECTOR3 *pOut, CONST D3DXVECTOR3 *pV, CONST D3D10_VIEWPORT *pViewport,
	 CONST D3DXMATRIX *pProjection, CONST D3DXMATRIX *pView, CONST D3DXMATRIX *pWorld);

	 // Project vector from screen space into object space
	 D3DXVECTOR3* WINAPI D3DXVec3Unproject
	 ( D3DXVECTOR3 *pOut, CONST D3DXVECTOR3 *pV, CONST D3D10_VIEWPORT *pViewport,
	 CONST D3DXMATRIX *pProjection, CONST D3DXMATRIX *pView, CONST D3DXMATRIX *pWorld);

	 // Project vector Array from object space into screen space
	 D3DXVECTOR3* WINAPI D3DXVec3ProjectArray
	 ( D3DXVECTOR3 *pOut, UINT OutStride,CONST D3DXVECTOR3 *pV, UINT VStride,CONST D3D10_VIEWPORT *pViewport,
	 CONST D3DXMATRIX *pProjection, CONST D3DXMATRIX *pView, CONST D3DXMATRIX *pWorld, UINT n);

	 // Project vector Array from screen space into object space
	 D3DXVECTOR3* WINAPI D3DXVec3UnprojectArray
	 ( D3DXVECTOR3 *pOut, UINT OutStride, CONST D3DXVECTOR3 *pV, UINT VStride, CONST D3D10_VIEWPORT *pViewport,
	 CONST D3DXMATRIX *pProjection, CONST D3DXMATRIX *pView, CONST D3DXMATRIX *pWorld, UINT n);

	 */

public:
	double x;
	double y;
	double z;
};

Vector3D operator * ( double f, const Vector3D& v )
{
    return Vector3D(f * v.x, f * v.y, f * v.z);
}

double Vector3DLength( const Vector3D *pV )
{
	return sqrtf(pV->x * pV->x + pV->y * pV->y + pV->z * pV->z);
}

Vector3D* Vector3DNormalize( Vector3D *pOut, const Vector3D *pV)
{
	double l = Vector3DLength(pV);
	pOut->x /= l;
	pOut->y /= l;
	pOut->z /= l;
	return pOut;
}

double Vector3DLengthSq( const Vector3D *pV )
{
	return pV->x * pV->x + pV->y * pV->y + pV->z * pV->z;
}

double Vector3DDot( const Vector3D *pV1, const Vector3D *pV2 )
{
	return pV1->x * pV2->x + pV1->y * pV2->y + pV1->z * pV2->z;
}

Vector3D* Vector3DCross( Vector3D *pOut, const Vector3D *pV1, const Vector3D *pV2 )
{
	Vector3D v;

	v.x = pV1->y * pV2->z - pV1->z * pV2->y;
	v.y = pV1->z * pV2->x - pV1->x * pV2->z;
	v.z = pV1->x * pV2->y - pV1->y * pV2->x;

	*pOut = v;
	return pOut;
}

Vector3D* Vector3DAdd( Vector3D *pOut, const Vector3D *pV1, const Vector3D *pV2 )
{
	pOut->x = pV1->x + pV2->x;
	pOut->y = pV1->y + pV2->y;
	pOut->z = pV1->z + pV2->z;
	return pOut;
}

Vector3D* Vector3DSubtract( Vector3D *pOut, const Vector3D *pV1, const Vector3D *pV2 )
{
	pOut->x = pV1->x - pV2->x;
	pOut->y = pV1->y - pV2->y;
	pOut->z = pV1->z - pV2->z;
	return pOut;
}

Vector3D* Vector3DMinimize( Vector3D *pOut, const Vector3D *pV1, const Vector3D *pV2 )
{
	pOut->x = pV1->x < pV2->x ? pV1->x : pV2->x;
	pOut->y = pV1->y < pV2->y ? pV1->y : pV2->y;
	pOut->z = pV1->z < pV2->z ? pV1->z : pV2->z;
	return pOut;
}

Vector3D* Vector3DMaximize( Vector3D *pOut, const Vector3D *pV1, const Vector3D *pV2 )
{
	pOut->x = pV1->x > pV2->x ? pV1->x : pV2->x;
	pOut->y = pV1->y > pV2->y ? pV1->y : pV2->y;
	pOut->z = pV1->z > pV2->z ? pV1->z : pV2->z;
	return pOut;
}

Vector3D* Vector3DScale( Vector3D *pOut, const Vector3D *pV, double s)
{
	pOut->x = pV->x * s;
	pOut->y = pV->y * s;
	pOut->z = pV->z * s;
	return pOut;
}

Vector3D* Vector3DLerp( Vector3D *pOut, const Vector3D *pV1, const Vector3D *pV2, double s )
{
	pOut->x = pV1->x + s * (pV2->x - pV1->x);
	pOut->y = pV1->y + s * (pV2->y - pV1->y);
	pOut->z = pV1->z + s * (pV2->z - pV1->z);
	return pOut;
}

Vector3D* Vector3DClamp01( Vector3D *pOut, const Vector3D *pV) {
	if (pV->x > 1.f)
		pOut->x = 1.f;
	else if(pV->x < 0.f)
		pOut->x = 0.f;

	if (pV->y > 1.f)
		pOut->y = 1.f;
	else if(pV->y < 0.f)
		pOut->y = 0.f;

	if (pV->z > 1.f)
		pOut->z = 1.f;
	else if(pV->z < 0.f)
		pOut->z = 0.f;

	return pOut;
}

#endif
