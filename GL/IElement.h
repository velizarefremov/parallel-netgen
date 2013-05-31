/*
 *  IElement.h
 *
 *  Created by Yusuf Yilmaz on 5/23/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _IELEMENT_H
#define _IELEMENT_H

#include "Ray.h"
#include "AABB.h"
#include "IntersectionData.h"

class IElement {

public:
	IElement(void)
	{

	}

	~IElement(void)
	{
	}

	// get bounding box.
	virtual AABB getBB() const
	{
		return AABB();
	}

	// get center of an element
	virtual Vector3D getCentroid() const
	{
		return Vector3D();
	}
};


#endif
