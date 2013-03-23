/*
 *  MeshVertex.h
 *
 *  Created by Yusuf Yilmaz on 5/23/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _MESH_VERTEX_H
#define _MESH_VERTEX_H

#include "Vector3D.h"

class MeshVertex
{
public:
	MeshVertex()
	{
		position2 = Vector3D(0.0,0.0,0.0);
		normal = 0;
		id = 0;
	}

	MeshVertex(const Vector3D &pos)
	{
		position2 = pos;
		normal = 0;
		id = 0;
	}

	~MeshVertex()
	{
		if (normal)
			delete normal;
	}

	void setNormal(Vector3D *n1)
	{
		normal = n1;
	}

	Vector3D* getNormal()
	{
		return normal;
	}

	Vector3D* getPosition()
	{
		return &position2;
	}

	void setPosition(Vector3D &point)
	{
		position2 = point;
	}

	unsigned int getIndex()
	{
		return index;
	}

	unsigned int getId()
	{
	    return id;
	}

	void setId(unsigned int _i)
	{
	    id = _i;
	}

	void setIndex(unsigned int i)
	{
		index = i;
	}

	double getValueAt(int index)
	{
		return position2[index];
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

	Vector3D position2;			// position of the point.
	Vector3D* normal;			// normal at that point.
	unsigned int index;			// used to find position inside the mesh.
	unsigned int id;    // ID of the Vertex.
    std::list<int> procids;      // Contains processor ID's the element is in.
};

#endif
