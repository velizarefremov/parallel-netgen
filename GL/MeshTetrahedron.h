#ifndef MESHTETRAHEDRON_H
#define MESHTETRAHEDRON_H

#include "Mesh.h"
#include "MeshVertex.h"
#include "IElement.h"
#include "AABB.h"

class Mesh;

class MeshTetrahedron : public IElement {

public:
	MeshTetrahedron()
	{
		parentMesh = 0;
		v[0] = 0;
		v[1] = 0;
		v[2] = 0;
        v[3] = 0;
        id = 0;
	}

	MeshTetrahedron(Mesh *pm, MeshVertex* v0, MeshVertex* v1, MeshVertex* v2, MeshVertex* v3)
	{
		parentMesh = pm;

		v[0] = v0;
		v[1] = v1;
		v[2] = v2;
		v[3] = v3;
        id = 0;
	}

	~MeshTetrahedron(void)
	{
	}


	MeshVertex* getVertex(int i)
	{
		return v[i];
	}

	// get center of the triangle.
	Vector3D getCentroid() const
	{
		return ( *(v[0]->getPosition()) + *(v[1]->getPosition()) + *(v[2]->getPosition()) + *(v[3]->getPosition()) ) / 4;
	}

    unsigned int getId()
	{
	    return id;
	}

	void setId(unsigned int _i)
	{
	    id = _i;
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

    virtual bool intersect(const Ray &ray, IntersectionData* iData)
	{
		/// Not implemented. Since it is not needed.

		return true;
	}

	// get bounding box.
	AABB getBB() const
	{
        // std::cout << " 1 " << std::endl;

		Vector3D p1 = *(v[0]->getPosition());
		// std::cout << " 2 " << std::endl;

		Vector3D p2 = *(v[2]->getPosition());
		// std::cout << " 3 " << std::endl;

		Vector3D p3 = *(v[3]->getPosition());
		// std::cout << " 4 " << std::endl;

		Vector3D p4 = *(v[1]->getPosition());


		Vector3D lower = p1, upper = p1;

		//set to max of p2 & p3
		for (int i = 0; i < 3; i++) {
			if (lower[i] > p2[i])
				lower[i] = p2[i];

			if (upper[i] < p2[i])
				upper[i] = p2[i];
		}
        // set to max of p1 & p2 & p3
		for (int i = 0; i < 3; i++) {
			if (lower[i] > p3[i])
				lower[i] = p3[i];

			if (upper[i] < p3[i])
				upper[i] = p3[i];
		}
		// set max of all 4.
        for (int i = 0; i < 3; i++) {
			if (lower[i] > p4[i])
				lower[i] = p4[i];

			if (upper[i] < p4[i])
				upper[i] = p4[i];
		}

		return AABB(lower, upper);
	}

private:

	MeshVertex* v[4];		// contains four vertices.
	Mesh* parentMesh;		// belongs to a mesh.
    unsigned int id;
    std::list<int> procids;      // Contains processor ID's the element is in
};


#endif // MESHTETRAHEDRON_H
