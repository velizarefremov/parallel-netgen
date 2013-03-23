/*
 *  Vector2D.h
 *  OpenGLTest
 *
 *  Created by Yusuf Yilmaz on 4/19/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *	Ripped of from DirectX.
 */
#ifndef _VECTOR2D_H
#define _VECTOR2D_H

#include <math.h>

class Vector2D
{
public:
	
	// constructors
	Vector2D() {};
	Vector2D( const double *pf ) 
	{
		x = pf[0];
		y = pf[1];
	}
	
	Vector2D( const Vector2D& v)
	{
		x = v.x;
		y = v.y;
	}
	
	Vector2D( double fx, double fy )
	{
		x = fx;
		y = fy;
	}
	
	// assignment operators
	Vector2D& operator += ( const Vector2D& v)
	{
		x += v.x;
		y += v.y;
		return *this;
	}
	
	Vector2D& operator -= ( const Vector2D& v)
	{
		x -= v.x;
		y -= v.y;
		return *this;
	}
	
	Vector2D& operator *= ( double f)
	{
		x *= f;
		y *= f;
		return *this;
	}
	
	Vector2D& operator /= ( double f)
	{
		double fInv = 1.0 / f;
		x *= fInv;
		y *= fInv;
		return *this;
	}
	
	// unary operators
	Vector2D operator + () const
	{
		return *this;
	}
	Vector2D operator - () const
	{
		return Vector2D(-x, -y);
	}
	
	// binary operators
	Vector2D operator + ( const Vector2D& v) const
	{
		return Vector2D(x + v.x, y + v.y);
	}
	
	Vector2D operator - ( const Vector2D& v) const
	{
		return Vector2D(x - v.x, y - v.y);
	}
	
	Vector2D operator * ( double f) const
	{
		return Vector2D(x * f, y * f);
	}
	
	Vector2D operator / ( double f) const
	{
		double fInv = 1.0 / f;
		return Vector2D(x * fInv, y * fInv);
	}
	
	bool operator == ( const Vector2D& v) const
	{
		return x == v.x && y == v.y;
	}
	
	bool operator != ( const Vector2D& v) const
	{
		return x != v.x || y != v.y;
	}
	
	
		
	// Other DirectX functions
	/* 
	 
	 // Hermite interpolation between position V1, tangent T1 (when s == 0)
	 // and position V2, tangent T2 (when s == 1).
	 D3DXVECTOR2* WINAPI D3DXVec2Hermite
	 ( D3DXVECTOR2 *pOut, CONST D3DXVECTOR2 *pV1, CONST D3DXVECTOR2 *pT1,
	 CONST D3DXVECTOR2 *pV2, CONST D3DXVECTOR2 *pT2, FLOAT s );
	 
	 // CatmullRom interpolation between V1 (when s == 0) and V2 (when s == 1)
	 D3DXVECTOR2* WINAPI D3DXVec2CatmullRom
	 ( D3DXVECTOR2 *pOut, CONST D3DXVECTOR2 *pV0, CONST D3DXVECTOR2 *pV1,
	 CONST D3DXVECTOR2 *pV2, CONST D3DXVECTOR2 *pV3, FLOAT s );
	 
	 // Barycentric coordinates.  V1 + f(V2-V1) + g(V3-V1)
	 D3DXVECTOR2* WINAPI D3DXVec2BaryCentric
	 ( D3DXVECTOR2 *pOut, CONST D3DXVECTOR2 *pV1, CONST D3DXVECTOR2 *pV2,
	 CONST D3DXVECTOR2 *pV3, FLOAT f, FLOAT g);
	 
	 // Transform (x, y, 0, 1) by matrix.
	 D3DXVECTOR4* WINAPI D3DXVec2Transform
	 ( D3DXVECTOR4 *pOut, CONST D3DXVECTOR2 *pV, CONST D3DXMATRIX *pM );
	 
	 // Transform (x, y, 0, 1) by matrix, project result back into w=1.
	 D3DXVECTOR2* WINAPI D3DXVec2TransformCoord
	 ( D3DXVECTOR2 *pOut, CONST D3DXVECTOR2 *pV, CONST D3DXMATRIX *pM );
	 
	 // Transform (x, y, 0, 0) by matrix.
	 D3DXVECTOR2* WINAPI D3DXVec2TransformNormal
	 ( D3DXVECTOR2 *pOut, CONST D3DXVECTOR2 *pV, CONST D3DXMATRIX *pM );
	 
	 // Transform Array (x, y, 0, 1) by matrix.
	 D3DXVECTOR4* WINAPI D3DXVec2TransformArray
	 ( D3DXVECTOR4 *pOut, UINT OutStride, CONST D3DXVECTOR2 *pV, UINT VStride, CONST D3DXMATRIX *pM, UINT n);
	 
	 // Transform Array (x, y, 0, 1) by matrix, project result back into w=1.
	 D3DXVECTOR2* WINAPI D3DXVec2TransformCoordArray
	 ( D3DXVECTOR2 *pOut, UINT OutStride, CONST D3DXVECTOR2 *pV, UINT VStride, CONST D3DXMATRIX *pM, UINT n );
	 
	 // Transform Array (x, y, 0, 0) by matrix.
	 D3DXVECTOR2* WINAPI D3DXVec2TransformNormalArray
	 ( D3DXVECTOR2 *pOut, UINT OutStride, CONST D3DXVECTOR2 *pV, UINT VStride, CONST D3DXMATRIX *pM, UINT n );
	 
	 */
	
	
public:
	double x;
	double y;
};

Vector2D operator * ( double f, const Vector2D& v)
{
	return Vector2D(f * v.x, f * v.y);
}

double Vector2DLength( const Vector2D *pV )
{
	return sqrtf(pV->x * pV->x + pV->y * pV->y);
}

double Vector2DLengthSq( const Vector2D *pV )
{
	return pV->x * pV->x + pV->y * pV->y;
}

double Vector2DDot( const Vector2D *pV1, const Vector2D *pV2 )
{
	return pV1->x * pV2->x + pV1->y * pV2->y;
}

double Vector2DCCW( const Vector2D *pV1, const Vector2D *pV2 )
{
	return pV1->x * pV2->y - pV1->y * pV2->x;
}

Vector2D* Vector2DAdd( Vector2D *pOut, const Vector2D *pV1, const Vector2D *pV2 )
{
	pOut->x = pV1->x + pV2->x;
	pOut->y = pV1->y + pV2->y;
	return pOut;
}

Vector2D* Vector2DSubtract( Vector2D *pOut, const Vector2D *pV1, const Vector2D *pV2 )
{
	pOut->x = pV1->x - pV2->x;
	pOut->y = pV1->y - pV2->y;
	return pOut;
}

Vector2D* Vector2DMinimize( Vector2D *pOut, const Vector2D *pV1, const Vector2D *pV2 )
{
	pOut->x = pV1->x < pV2->x ? pV1->x : pV2->x;
	pOut->y = pV1->y < pV2->y ? pV1->y : pV2->y;
	return pOut;
}

Vector2D* Vector2DMaximize( Vector2D *pOut, const Vector2D *pV1, const Vector2D *pV2 )
{
	pOut->x = pV1->x > pV2->x ? pV1->x : pV2->x;
	pOut->y = pV1->y > pV2->y ? pV1->y : pV2->y;
	return pOut;
}

Vector2D* Vector2DScale( Vector2D *pOut, const Vector2D *pV, double s )
{
	pOut->x = pV->x * s;
	pOut->y = pV->y * s;
	return pOut;
}

Vector2D* Vector2DLerp( Vector2D *pOut, const Vector2D *pV1, const Vector2D *pV2, double s )
{
	pOut->x = pV1->x + s * (pV2->x - pV1->x);
	pOut->y = pV1->y + s * (pV2->y - pV1->y);
	return pOut;
}

Vector2D* Vector2DNormalize( Vector2D *pOut, const Vector2D *pV )
{
	double l = Vector2DLength(pV);
	pOut->x /= l;
	pOut->y /= l;
	return pOut;
}

Vector2D* Vector2DClamp01( Vector2D *pOut, const Vector2D *pV) { 
	if (pV->x > 1.f) 
		pOut->x = 1.f;
	else if(pV->x < 0.f) 
		pOut->x = 0.f;
	
	if (pV->y > 1.f)
		pOut->y = 1.f;
	else if(pV->y < 0.f) 
		pOut->y = 0.f;
	
	return pOut;
}


#endif
