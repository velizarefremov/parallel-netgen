//#include<GL/glut.h>
//#include<GL/gl.h>
#include <iostream>

// Include Utility Libraries
#include "AABB.h"
#include "Vector2D.h"
#include "Vector3D.h"
#include "Vector4D.h"

// Mesh and Tree Libraries
#include "IElement.h"
#include "BSPTreeNode.h"
#include "MeshVertex.h"
#include "MeshTriangle.h"
#include "MeshTetrahedron.h"
#include "Mesh.h"
#include "ElmerVertex.h"
#include "ElmerTriangle.h"
#include "TreeBuilder.h"

// Netgen Library Header
namespace nglib{
#include <nglib.h>
}

// Globals
int usedmethod = 2;     // 2 for the metis partitioning and refinement method.

int recursionDepth = 5; //the recursion depth of the tree 2^n processors. Ex: If n=10 processor number is 1024.
int drawDepth = 3; // the recursion depth to draw the split planes.
int minsize = 100000;

bool isglut = false;
bool justonepart = false;
const char* partitname;
const char* inputFileName = "cone.geo";

int createElmer = 0;
int createGeom = 0;

TreeBuilder* treebuilder;

// temporary structures for reading from file.
typedef struct Vertex3D {
	double x;
	double y;
	double z;
} Vertex3D;

typedef struct Triangle2 {
	long a;				// Does not hold Vertex3D but an index to it.
	long b;
	long c;
} Triangle2;

typedef struct Tetrahedron4{
    long a;
    long b;
    long c;
    long d;
}Tetrahedron4;

using namespace nglib;

// Dosyayi okuyan fonksiyon.
int loadNGData(const char * filename, TreeBuilder * treebuilder, double sizing, Vector3D translation, Vector4D col);

// Total vertex ve triangle sayisi.
long nvertices, ntriangles, ntetrahedron;

Vertex3D * verticelist;
Triangle2 * trianglelist;
Tetrahedron4 * tetralist;

double * vertice2list;
int * triangle2list;
int * tetra2list;

long curvertice;
long curtriangle;

unsigned long getTime() //OK
{
	time_t seconds;
	seconds = time(NULL);

	return (unsigned long)seconds;
}

AABB preprocessing(Mesh * _mesh,double _newRadius,Vector3D & _translate) //OK
{
	// compute bounding box
	Vector3D minPosition;
	Vector3D maxPosition;
	Vector3D centerPosition;
	// double radius=0.0;
	maxPosition[0]=-3.4e38;
	maxPosition[1]=-3.4e38;
	maxPosition[2]=-3.4e38;

	minPosition[0]=3.4e38;
	minPosition[1]=3.4e38;
	minPosition[2]=3.4e38;

	for(unsigned int i=0;i<_mesh->numberOfVertices();i++)
	{
		Vector3D point = *_mesh->getVertex(i)->getPosition();
		if(point[0] < minPosition[0])
		{
			minPosition[0] = point[0];
		}
		if(point[1] <minPosition[1])
		{
			minPosition[1] = point[1];
		}
		if(point[2] <minPosition[2])
		{
			minPosition[2] = point[2];
		}

		if(point[0] >maxPosition[0])
		{
			maxPosition[0] = point[0];
		}
		if(point[1] >maxPosition[1])
		{
			maxPosition[1] = point[1];
		}
		if(point[2] >maxPosition[2])
		{
			maxPosition[2] = point[2];
		}

	}

    return AABB(minPosition, maxPosition);
}

void calcVertexNormals(Mesh * _mesh)
{
	// Calculate vertex normals.

	std::vector<Vector3D> vertexNormalList;
	vertexNormalList.clear();
	vertexNormalList.resize(_mesh->numberOfVertices());

	for (unsigned int i = 0; i < _mesh->numberOfVertices(); i++)
	{
		vertexNormalList[i] = Vector3D(0,0,0);
	}

	for (unsigned int i = 0; i < _mesh->numberOfFaces(); i++) {

		MeshTriangle * triangle = _mesh->getFace(i);

		// For each triangle, calculate it's face normal.

		Vector3D point1 = *(triangle->getVertex(0)->getPosition());
		Vector3D point2 = *(triangle->getVertex(1)->getPosition());
		Vector3D point3 = *(triangle->getVertex(2)->getPosition());

		Vector3D vec1 = point2 - point1;
		Vector3D vec2 = point3 - point1;
		Vector3D faceNormal = vec1.Cross(vec2);

		// Add this face normal to the normal vector for each of the triangle's three indices.

		vertexNormalList[triangle->getVertex(0)->getIndex()] += faceNormal;
		vertexNormalList[triangle->getVertex(1)->getIndex()] += faceNormal;
		vertexNormalList[triangle->getVertex(2)->getIndex()] += faceNormal;

	}

	// Finally normalize all of the normals.
	for (unsigned int i = 0; i < _mesh->numberOfVertices(); i++)
	{
		vertexNormalList[i].Normalize();
		Vector3D* normal = new Vector3D(vertexNormalList[i][0],vertexNormalList[i][1],vertexNormalList[i][2]);
		_mesh->getVertex(i)->setNormal(normal);
	}
}

TreeBuilder* createTreeBuilder() // OK.
{
	// create treebuilder
	TreeBuilder* treebuilder = new TreeBuilder();

	loadNGData(inputFileName, treebuilder, 10.0, Vector3D(0.0,0.0,0.0), Vector4D(1.0,0.0,0.0,1.0));

	return treebuilder;
}

int loadNGData(const char * filename, TreeBuilder * treebuilder, double sizing, Vector3D translation, Vector4D col) //OK
{
    using namespace nglib;
    Mesh *m = new Mesh(0,0,0);
	treebuilder->addMesh(m);

    Ng_Init();

    vertice2list = new double[6000000];
    tetra2list = new int[8000000];
    triangle2list = new int[6000000];

    Ng_CSG_GenerateMesh(filename, &nvertices, &ntriangles, &ntetrahedron, vertice2list, triangle2list, tetra2list);
    std::cout << "Successfully loaded .geo File: " << filename << std::endl;

    std::cout << "Hello Numbers Are: " << nvertices << " " << ntriangles << " " << ntetrahedron << std::endl;
    for (int i=0; i<nvertices; i++) {
		MeshVertex *mv = new MeshVertex(Vector3D(vertice2list[i*3], vertice2list[i*3+1], vertice2list[i*3+2]));
		m->addVertex(mv);

		if (i%50000 == 0) {
			std::cout << "Vertex No: " << i <<std::endl;
		}
	}

    for (int i=0; i<ntriangles; i++) {
		MeshTriangle *mt = new MeshTriangle(m, m->getVertex(triangle2list[i*3]-1),
											m->getVertex(triangle2list[i*3+1]-1),
											m->getVertex(triangle2list[i*3+2]-1));

		treebuilder->addSurfaceElem(mt);

		if (i%50000 == 0) {
			std::cout << "Triangle No: " << i <<std::endl;
		}
	}

	for (int i=0; i<ntetrahedron; i++) {
		MeshTetrahedron *mt = new MeshTetrahedron(m, m->getVertex(tetra2list[i*4]-1),
											m->getVertex(tetra2list[i*4+1]-1),
											m->getVertex(tetra2list[i*4+2]-1),
											m->getVertex(tetra2list[i*4+3]-1));
		m->addTetrahedron(mt);
		treebuilder->addElement(mt);

		if (i%50000 == 0) {
			std::cout << "Tetra No: " << i <<std::endl;
		}
	}

	calcVertexNormals(m);
	AABB bb = preprocessing(m,sizing,translation);
	m->setBB(bb);

    delete [] tetra2list;
    delete [] vertice2list;
    delete [] triangle2list;

    return 0;
}

int main(int argc, char **argv) {

    int rank2, size2;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank2);
    MPI_Comm_size(MPI_COMM_WORLD, &size2);

    int time1; // Start Time.
    int time2; // Init Finish Time.
    int time3; // Skin Mesh Time.
    int time4; // Refinement Finish Time.
    int time5; // Migration Finish Time.

    std::string myfile;

        if(argc >= 2)
        {
            inputFileName = argv[1];
            std::cout << "Selected File is: " << inputFileName << std::endl;
        }

        if(argc >= 3)
        {
            usedmethod = atoi(argv[2]);
            std::cout << "Selected Method: " << usedmethod << std::endl;
        }
        if(argc >=4)
        {
            minsize = atoi(argv[3]);
            std::cout << "Minimum Size is set to: " << minsize << std::endl;
        }
        if(argc >=5)
        {
            recursionDepth = atoi(argv[4]);
            std::cout << "Repeat set to" << recursionDepth << std::endl;
        }
        if(argc >=6)
        {
            createElmer = atoi(argv[5]);
            std::cout << "CreateElmer set to: " << createElmer << std::endl;
        }
        if(argc >=7)
        {
            createGeom = atoi(argv[6]);
            std::cout << "CreateGeom set to: " << createGeom << std::endl;
        }

    if(usedmethod == 1)
    {

        if(rank2 == 0)
        {
            srand ( time(NULL) );

            treebuilder = createTreeBuilder();
            std::cout << "TreeBuilder created" << std::endl;

            // Build Tree
            if (treebuilder != NULL)
            {
                unsigned long startTime = getTime();

                treebuilder->buildBSPTree(recursionDepth);

                unsigned long stopTime = getTime();

                std::cout << "Total build time: " << (stopTime - startTime) << " sec \n";
            }

            std::cout << "BSP-Tree built" << std::endl;
            std::cout << "Total Nodes: " << treebuilder->gettotalNodes() << std::endl;
            std::cout << "Total Leaf Nodes: " << treebuilder->gettotalLeafs() << std::endl;

            treebuilder->createGeoFile(inputFileName, "out.geo");
            std::cout << "Partitioned GEO File is created!" << std::endl;

            treebuilder->createPartitions("out.geo", minsize);
            std::cout << "Partitions are created!" << std::endl;

            treebuilder->setStringSend(size2);

            treebuilder->createPartFiles();
            std::cout << "Partitions are saved into files!" << std::endl;
        }

        MPI_Barrier(MPI_COMM_WORLD);

        std::cout << "Calling MPI" << std::endl;
        treebuilder->mpifunc(rank2, size2);

        MPI_Barrier(MPI_COMM_WORLD);

    }
    else        // Method 2 and 3
    {
        MeshMig* localMesh = new MeshMig(rank2, size2, false);


        time1 = getTime();
        localMesh->createMetisMesh(MPI_COMM_WORLD, inputFileName, minsize);
        time2 = getTime();
        if(usedmethod == 3)
            localMesh->skinMeshv2();
        MPI_Barrier(MPI_COMM_WORLD);
        time3 = getTime();
        for(int i=0; i<recursionDepth; i++)
        {
            if(usedmethod == 3)
            {
                localMesh->refineParallelV2(MPI_COMM_WORLD, inputFileName);
                MPI_Barrier(MPI_COMM_WORLD);
            }
            else if(usedmethod == 2)
            {
                localMesh->refineParallel(MPI_COMM_WORLD);
                MPI_Barrier(MPI_COMM_WORLD);
            }
        }

        if(usedmethod == 3)
            localMesh->createVolume();

        time4 = getTime();

        if(createElmer == 1)
            localMesh->createElmerOutput();

        if(createGeom == 1)
            localMesh->createGeomViewOutput();

        if(rank2 == 0)
        {
            std::cout << "INIT TIME TAKEN: " << time3 - time1 << std::endl;
            std::cout << "PARALLEL REFINEMENT TIME TAKEN: " << time4 - time3 << std::endl;
        }
    }

    std::cout << "Finalizing: " << rank2 << std::endl;
    MPI_Finalize();
}

