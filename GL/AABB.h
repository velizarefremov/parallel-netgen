/*
 *  AABB.h
 *
 *  Created by Yusuf Yilmaz on 5/22/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _AABB_H
#define _AABB_H

#include "Vector3D.h"

struct AABB {
	AABB()
	{
	}
	AABB(Vector3D l, Vector3D u)
	{
		lowerCorner=l;
		upperCorner=u;
	}

	//info about axis aligned bounding box
	Vector3D lowerCorner;
	Vector3D upperCorner;
};

#endif
