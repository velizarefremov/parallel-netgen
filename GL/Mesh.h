/*
 *  Mesh.h
 *
 *  Created by Yusuf Yilmaz on 5/23/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _MESH_H
#define _MESH_H

#include <vector>
#include "AABB.h"
#include "MeshTriangle.h"
#include "MeshTetrahedron.h"
#include "MeshVertex.h"

class MeshVertex;
class MeshTriangle;
class MeshTetrahedron;

class Mesh {

public:
	Mesh(int n_vertices, int n_triangles, int n_tetrahedrons)		// create mesh with a specific number of vertices and faces
	{
		vertices.clear();
		vertices.resize(n_vertices);

		triangles.clear();
		triangles.resize(n_triangles);

		tetrahedrons.clear();
		tetrahedrons.resize(n_tetrahedrons);
	}

	~Mesh(void)
	{
		for (unsigned int i=0; i < vertices.size(); i++)
			delete vertices[i];
		vertices.clear();

		triangles.clear();		// we do not clear each element here since, each MeshTriangle is also added to the Scene element list.
		tetrahedrons.clear();
	}

    void addTetrahedron(MeshTetrahedron *m)
    {
        tetrahedrons.push_back(m);
    }

	void addTriangle(MeshTriangle *m)
	{
		triangles.push_back(m);
	}

	void addVertex(MeshVertex *v)
	{
		v->setIndex(static_cast<unsigned int>(vertices.size()));
		vertices.push_back(v);
	}

	unsigned int numberOfVertices()
	{
		return (unsigned int)vertices.size();
	}

	unsigned int numberOfFaces()
	{
		return (unsigned int)triangles.size();
	}

	unsigned int numberOfVolumeElems()
	{
	    return (unsigned int)tetrahedrons.size();
	}

	MeshVertex *getVertex(int i)
	{
		return vertices[i];
	}

	MeshTriangle *getFace(int i)
	{
		return triangles[i];
	}

	MeshTetrahedron *getVolumeElement(int i)
	{
	    return tetrahedrons[i];
	}

	void setBB( AABB aabb )
	{
		bb.lowerCorner = aabb.lowerCorner;
		bb.upperCorner = aabb.upperCorner;
	}

	AABB getBB()
	{
		return bb;
	}


private:
	Mesh() {};

	// Contains a vertex and triangle list.
	std::vector<MeshTetrahedron*> tetrahedrons;
	std::vector<MeshTriangle*> triangles;
	std::vector<MeshVertex*> vertices;

	// bounding box of the mesh.
	AABB bb;
};

#endif
