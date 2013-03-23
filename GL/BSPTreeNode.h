/*
 *  BSPTreeNode.h
 *
 *  Created by Yusuf Yilmaz on 5/23/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _BSPTREENODE_H
#define _BSPTREENODE_H

#include <list>
#include <map>
#include <string>
#include "AABB.h"
#include "IElement.h"
#include "Vector3D.h"
#include "ElmerTriangle.h"
#include "ElmerVertex.h"

typedef struct BSPTreeNode_{
	AABB bb;		// Bounding Box of the Node
	int level;		// Depth Level of the tree.
	int axis;		// x=0 y=1 z=2	// x means y-z is the split plane. y means x-z is split plane etc.
	double splitCoordinate;		// value of the splitting coordinate.
	std::list<IElement*> elements; // element list. // NULL or empty if not leaf node.

	std::map<int,ElmerVertex*> verts;
	std::list<ElmerTriangle*> tris;
	std::list<bool> orientation;

	// Plane vertices to draw;
	Vector3D ll;
	Vector3D lr;
	Vector3D ur;
	Vector3D ul;

	Vector3D planeclr; // for visualization purpose

    Vector3D splitPoint;    // plane point.
    Vector3D planeNormal;   // plane normal
    Vector3D parentPlane;   // plane normal of parent.
    double d;               // plane d value

	std::string partname;
	int partId;             // Only available for endnode elements.
	std::string planename;

	// Left and Right nodes.
	BSPTreeNode_ *leftChildNode;
	BSPTreeNode_ *rightChildNode;

} BSPTreeNode;


#endif
