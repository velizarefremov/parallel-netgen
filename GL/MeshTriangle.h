/*
 *  MeshTriangle.h
 *
 *  Created by Yusuf Yilmaz on 5/23/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _MESH_TRIANGLE_H
#define _MESH_TRIANGLE_H

#include "Mesh.h"
#include "MeshVertex.h"
#include "IElement.h"
#include "AABB.h"

class Mesh;

class MeshTriangle : public IElement {

public:
	MeshTriangle()
	{
		parentMesh = 0;
		v[0] = 0;
		v[1] = 0;
		v[2] = 0;
        id = 0;
	}

	MeshTriangle(Mesh *pm, MeshVertex* v0, MeshVertex* v1, MeshVertex* v2)
	{
		parentMesh = pm;

		v[0] = v0;
		v[1] = v1;
		v[2] = v2;
        id = 0;
	}

	~MeshTriangle(void)
	{
	}


	MeshVertex* getVertex(int i)
	{
		return v[i];
	}

	// get center of the triangle.
	Vector3D getCentroid() const
	{
		return ( *(v[0]->getPosition()) + *(v[1]->getPosition()) + *(v[2]->getPosition()) ) / 3;
	}

	unsigned int getId()
	{
	    return id;
	}

	void setId(unsigned int _i)
	{
	    id = _i;
	}

	// get bounding box.
	AABB getBB() const
	{
		Vector3D p1 = *(v[0]->getPosition());
		Vector3D p2 = *(v[1]->getPosition());
		Vector3D p3 = *(v[2]->getPosition());

		Vector3D lower = p1, upper = p1;

		//set to max of p2 & p3
		for (int i = 0; i < 3; i++) {
			if (lower[i] > p2[i])
				lower[i] = p2[i];

			if (upper[i] < p2[i])
				upper[i] = p2[i];
		}

		for (int i = 0; i < 3; i++) {
			if (lower[i] > p3[i])
				lower[i] = p3[i];

			if (upper[i] < p3[i])
				upper[i] = p3[i];
		}

		return AABB(lower, upper);
	}

	void addProc(int id)
    {
        procids.push_back(id);
    }

    bool isShared()
    {
        if(procids.size() > 1)
            return true;
        else
            return false;
    }

private:

	MeshVertex* v[3];		// contains three vertices.
	Mesh* parentMesh;		// belongs to a mesh.
    unsigned int id;
    std::list<int> procids;      // Contains processor ID's the element is in
};

#endif
