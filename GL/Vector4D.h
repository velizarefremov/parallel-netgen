/*
 *  Vector4D.h
 *
 *  Created by Yusuf Yilmaz on 4/19/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *	Ripped of from DirectX.
 */

#ifndef _VECTOR4D_H
#define _VECTOR4D_H

#include <math.h>

class Vector4D
{
public:
    Vector4D()
	{
		Vector4D(0.0,0.0,0.0,1.0);
	}
    Vector4D( const double* pf)
	{
		x = pf[0];
		y = pf[1];
		z = pf[2];
		w = pf[3];
	}

    Vector4D( const Vector3D& v, double f )
	{
		x = v.x;
		y = v.y;
		z = v.z;
		w = f;
	}

    Vector4D( double fx, double fy, double fz, double fw )
	{
		x = fx;
		y = fy;
		z = fz;
		w = fw;
	}


    // assignment operators
    Vector4D& operator += ( const Vector4D& v)
	{
		x += v.x;
		y += v.y;
		z += v.z;
		w += v.w;
		return *this;
	}

    Vector4D& operator -= ( const Vector4D& v)
	{
		x -= v.x;
		y -= v.y;
		z -= v.z;
		w -= v.w;
		return *this;
	}

    Vector4D& operator *= ( double f)
	{
		x *= f;
		y *= f;
		z *= f;
		w *= f;
		return *this;
	}

    Vector4D& operator /= ( double f)
	{
		double fInv = 1.0 / f;
		x *= fInv;
		y *= fInv;
		z *= fInv;
		w *= fInv;
		return *this;
	}

    // unary operators
    Vector4D operator + () const
	{
		return *this;
	}

    Vector4D operator - () const
	{
		return Vector4D(-x, -y, -z, -w);
	}


    // binary operators
    Vector4D operator + ( const Vector4D& v) const
	{
		return Vector4D(x + v.x, y + v.y, z + v.z, w + v.w);
	}

    Vector4D operator - ( const Vector4D& v) const
	{
		return Vector4D(x - v.x, y - v.y, z - v.z, w - v.w);
	}

    Vector4D operator * ( double f) const
	{
		return Vector4D(x * f, y * f, z * f, w * f);
	}

    Vector4D operator / ( double f) const
	{
		double fInv = 1.0 / f;
		return Vector4D(x * fInv, y * fInv, z * fInv, w * fInv);
	}

    bool operator == ( const Vector4D& v) const
	{
		return x == v.x && y == v.y && z == v.z && w == v.w;
	}

    bool operator != ( const Vector4D& v) const
	{
		return x != v.x || y != v.y || z != v.z || w != v.w;
	}

	Vector4D clamp01() {
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

		if (w>1.f)
			w=1.f;
		else if (w<0.f)
			w=0.f;

		return *this;
	}

	Vector4D componentMul(const Vector4D &v) const {
		return Vector4D(x*v.x, y*v.y, z*v.z, w*v.w);
	}


	// Other DirectX functions
	/*

	 // Cross-product in 4 dimensions.
	 D3DXVECTOR4* WINAPI D3DXVec4Cross
	 ( D3DXVECTOR4 *pOut, CONST D3DXVECTOR4 *pV1, CONST D3DXVECTOR4 *pV2,
	 CONST D3DXVECTOR4 *pV3);

	 // Hermite interpolation between position V1, tangent T1 (when s == 0)
	 // and position V2, tangent T2 (when s == 1).
	 D3DXVECTOR4* WINAPI D3DXVec4Hermite
	 ( D3DXVECTOR4 *pOut, CONST D3DXVECTOR4 *pV1, CONST D3DXVECTOR4 *pT1,
	 CONST D3DXVECTOR4 *pV2, CONST D3DXVECTOR4 *pT2, FLOAT s );

	 // CatmullRom interpolation between V1 (when s == 0) and V2 (when s == 1)
	 D3DXVECTOR4* WINAPI D3DXVec4CatmullRom
	 ( D3DXVECTOR4 *pOut, CONST D3DXVECTOR4 *pV0, CONST D3DXVECTOR4 *pV1,
	 CONST D3DXVECTOR4 *pV2, CONST D3DXVECTOR4 *pV3, FLOAT s );

	 // Barycentric coordinates.  V1 + f(V2-V1) + g(V3-V1)
	 D3DXVECTOR4* WINAPI D3DXVec4BaryCentric
	 ( D3DXVECTOR4 *pOut, CONST D3DXVECTOR4 *pV1, CONST D3DXVECTOR4 *pV2,
	 CONST D3DXVECTOR4 *pV3, FLOAT f, FLOAT g);

	 // Transform vector by matrix.
	 D3DXVECTOR4* WINAPI D3DXVec4Transform
	 ( D3DXVECTOR4 *pOut, CONST D3DXVECTOR4 *pV, CONST D3DXMATRIX *pM );

	 // Transform vector array by matrix.
	 D3DXVECTOR4* WINAPI D3DXVec4TransformArray
	 ( D3DXVECTOR4 *pOut, UINT OutStride, CONST D3DXVECTOR4 *pV, UINT VStride, CONST D3DXMATRIX *pM, UINT n );

	 */

public:
    double x;
	double y;
	double z;
	double w;
};

Vector4D operator * ( double f, const Vector4D& v )
{
    return Vector4D(f * v.x, f * v.y, f * v.z, f * v.w);
}

double Vector4DLength( const Vector4D *pV )
{
	return sqrtf(pV->x * pV->x + pV->y * pV->y + pV->z * pV->z + pV->w * pV->w);
}

double Vector4DLengthSq( const Vector4D *pV )
{
	return pV->x * pV->x + pV->y * pV->y + pV->z * pV->z + pV->w * pV->w;
}

double Vector4DDot( const Vector4D *pV1, const Vector4D *pV2 )
{
	return pV1->x * pV2->x + pV1->y * pV2->y + pV1->z * pV2->z + pV1->w * pV2->w;
}

Vector4D* Vector4DAdd( Vector4D *pOut, const Vector4D *pV1, const Vector4D *pV2)
{
	pOut->x = pV1->x + pV2->x;
	pOut->y = pV1->y + pV2->y;
	pOut->z = pV1->z + pV2->z;
	pOut->w = pV1->w + pV2->w;
	return pOut;
}

Vector4D* Vector4DSubtract( Vector4D *pOut, const Vector4D *pV1, const Vector4D *pV2)
{
	pOut->x = pV1->x - pV2->x;
	pOut->y = pV1->y - pV2->y;
	pOut->z = pV1->z - pV2->z;
	pOut->w = pV1->w - pV2->w;
	return pOut;
}

Vector4D* Vector4DMinimize( Vector4D *pOut, const Vector4D *pV1, const Vector4D *pV2)
{
	pOut->x = pV1->x < pV2->x ? pV1->x : pV2->x;
	pOut->y = pV1->y < pV2->y ? pV1->y : pV2->y;
	pOut->z = pV1->z < pV2->z ? pV1->z : pV2->z;
	pOut->w = pV1->w < pV2->w ? pV1->w : pV2->w;
	return pOut;
}

Vector4D* Vector4DMaximize( Vector4D *pOut, const Vector4D *pV1, const Vector4D *pV2)
{
	pOut->x = pV1->x > pV2->x ? pV1->x : pV2->x;
	pOut->y = pV1->y > pV2->y ? pV1->y : pV2->y;
	pOut->z = pV1->z > pV2->z ? pV1->z : pV2->z;
	pOut->w = pV1->w > pV2->w ? pV1->w : pV2->w;
	return pOut;
}

Vector4D* Vector4DScale( Vector4D *pOut, const Vector4D *pV, double s)
{
	pOut->x = pV->x * s;
	pOut->y = pV->y * s;
	pOut->z = pV->z * s;
	pOut->w = pV->w * s;
	return pOut;
}

Vector4D* Vector4DLerp( Vector4D *pOut, const Vector4D *pV1, const Vector4D *pV2, double s )
{
	pOut->x = pV1->x + s * (pV2->x - pV1->x);
	pOut->y = pV1->y + s * (pV2->y - pV1->y);
	pOut->z = pV1->z + s * (pV2->z - pV1->z);
	pOut->w = pV1->w + s * (pV2->w - pV1->w);
	return pOut;
}

Vector4D* Vector4DNormalize( Vector4D *pOut, const Vector4D *pV)
{
	double l = Vector4DLength(pV);
	pOut->x /= l;
	pOut->y /= l;
	pOut->z /= l;
	pOut->w /= l;
	return pOut;
}

Vector4D* Vector4DClamp01( Vector4D *pOut, const Vector4D *pV) {
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

	if (pV->w > 1.f)
		pOut->w = 1.f;
	else if(pV->w < 0.f)
		pOut->w = 0.f;

	return pOut;
}

#endif
