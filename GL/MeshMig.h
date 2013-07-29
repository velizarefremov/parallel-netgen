/*
 *  Mesh2D.h
 *  MYMIG
 *
 *  Created by Yusuf Yilmaz on 5/2/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

/*
 *  MeshClass.h
 *  MYMIG
 *
 *  Created by Yusuf Yilmaz on 5/2/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <vector>
#include <map>
#include <list>
#include <set>
#include <fstream>
#include <sstream>

namespace nglib {
#include <nglib.h>
}

using namespace std;

#define NUM_VOLM_VERT 4
#define NUM_SURF_VERT 3
#define NUM_VERT_SIZE 3

// Vertex
typedef struct Vertex {
	double         coor[NUM_VERT_SIZE];		// Coordinates of the vertex
	int            idx;						// global ID of the vertex
	// unsigned char  vtype;					// type of the vertex: i-internal; s-surface; p-partition; m-mixed(partition+surface)
	int			   useCount;				// count of parent volume element.
} Vertex;

typedef Vertex * pVertex;

// Boundary Element
typedef struct BoundaryElement {
	int    verts[NUM_SURF_VERT];			// Global ID's of the containing vertices.
	int    geoid;							// ID of the geometric surface the element belongs to.
	int    parentelement;					// Global ID of the Parent Volume Element
	int	   idx;								// Local ID of the surface element.
} BoundaryElement;

typedef BoundaryElement *pBoundary;

// Volume Element
typedef struct VolumeElement{
	int    verts[NUM_VOLM_VERT];			// Global ID's of the containing vertices.
	int    idx;								// Global ID of the volume element.
	// bool   deleted;							// false default. True if the element is deleted.
} VolumeElement;

typedef VolumeElement * pVolume ;

// Some useful typedefs.
typedef map<int,int> MapI2I ;
typedef map<int,list<int> > MapI2Ilist ;
typedef map<int,vector<VolumeElement> > MapI2vVol ;
typedef map<int,vector<Vertex> > MapI2vVertex ;

class Index2
{
    public:
        int x[2];

        Index2()
        {
            x[0] = 0;
            x[1] = 0;
        }

        Index2(int _x, int _y)
        {
            x[0] = _x;
            x[1] = _y;
        }

        void swapel(int a, int b)
        {
            int tmp;
            tmp = x[a];
            x[a] = x[b];
            x[b] = tmp;
        }

        void Sort()
        {
            if(x[1] < x[0]) swapel(1,0);
        }
};

class Index3
{
    public:
        int x[3];

        Index3()
        {
            x[0] = 0;
            x[1] = 0;
            x[2] = 0;
        }

        Index3(int _x, int _y, int _z)
        {
            x[0] = _x;
            x[1] = _y;
            x[2] = _z;
        }

        void swapel(int a, int b)
        {
            int tmp;
            tmp = x[a];
            x[a] = x[b];
            x[b] = tmp;
        }

        void Sort()
        {
            if(x[1] < x[0]) swapel(1,0);
            if(x[2] < x[1]) swapel(2,1);
            if(x[1] < x[0]) swapel(1,0);
        }
};

bool fncomp(Index2 in1, Index2 in2)
{
    if(in1.x[0] < in2.x[0])
        return true;
    else if(in1.x[0] > in2.x[0])
        return false;

    if(in1.x[1] < in2.x[1])
        return true;
    else if(in1.x[1] > in2.x[1])
        return false;

    return false;
}

bool fncomp(Index3 in1, Index3 in2)
{
    if(in1.x[0] < in2.x[0])
        return true;
    else if(in1.x[0] > in2.x[0])
        return false;

    if(in1.x[1] < in2.x[1])
        return true;
    else if(in1.x[1] > in2.x[1])
        return false;

    if(in1.x[2] < in2.x[2])
        return true;
    else if(in1.x[2] > in2.x[2])
        return false;

    return false;
}

// Mesh Class
class MeshMig {
private:

	vector<Vertex>   verts ;
	vector<BoundaryElement>    boundelems ;
	vector<BoundaryElement>    partbound;
	vector<VolumeElement>     volelems ;
	map<int,string> sharedVert;
	map<int,string> sharedTempVert;

    bool boundset;

	int *procid;			// Processor Id's to send the elements.
	int *elids;				// Element ID list.
	int lsize;				// Size of this list.

	vector<Vertex> migVertex;	// Vector containing vertices to be send.
	vector<BoundaryElement> migBoundary;	// Vector containing boundary elements to be send.
	vector<VolumeElement> migVolume;	// Vector containing volume elements to be send.

	std::string absentData;
	std::string excessData;
	std::string excessData2;

	MapI2I vertmap;
	MapI2I volmap;
	MapI2I volbound;

	int		 numvertex;
	int		 numboundary;
	int		 numvolume;
	int		 numshared;
	int      numdirichlet;

	void readMeshHeader(const char* filename)
	{
		ifstream infile(filename);

		infile >> numvertex >> numvolume >> numboundary >> numshared;

		infile.close();
	}

	void readMeshNodes(const char* filename)
	{
		ifstream infile(filename);

		double x, y;
		int glid, tmp;

		verts.resize(numvertex);

		for (int i=0; i<numvertex; i++) {
			infile >> glid >> tmp >> x >> y;

			verts[i].coor[0] = x;
			verts[i].coor[1] = y;
			verts[i].idx = glid;
			verts[i].useCount = 0;
			// verts[i].vtype = 'i';

			vertmap[glid] = i;
		}

		infile.close();
	}

	void readMeshElements(const char* filename)
	{
		ifstream infile(filename);

		int glid, x, y, z;

		volelems.resize(numvolume);

		for (int i=0; i<numvolume; i++) {
			infile >> glid >> x >> y >> z;

			volelems[i].verts[0] = x;
			volelems[i].verts[1] = y;
			volelems[i].verts[2] = z;
			volelems[i].idx = glid;
			// volelems[i].deleted = false;

			// increase usage of vertices.
			verts[vertmap.find(x)->second].useCount++;
			verts[vertmap.find(y)->second].useCount++;
			verts[vertmap.find(z)->second].useCount++;

			volmap[glid] = i;
		}

		infile.close();
	}

	void readMeshBoundary(const char* filename)
	{
		ifstream infile(filename);

		int lid, x, y, par, gfid;

		boundelems.resize(numboundary);

		for (int i=0; i<numboundary; i++) {
			infile >> lid >> x >> y >> par >> gfid;

			boundelems[i].verts[0] = x;
			boundelems[i].verts[1] = y;
			boundelems[i].idx = lid;
			boundelems[i].geoid = gfid;
			boundelems[i].parentelement = par;

			volbound[par] = lid;
		}

		infile.close();
	}

	void readMeshShared(const char* filename)
	{
		ifstream infile(filename);
		stringstream oss;

		int vid;
		int count;
		int prid;

		for (int i=0; i<numshared; i++) {

			infile >> vid >> count;

			// First one is the owner.
			infile >> prid;

			oss << count << " " << prid;

			for (int j=1; j<count; j++) {
				// The other processors.
				infile >> prid;
				oss << " " << prid;
			}

			sharedVert[vid] = oss.str();

			oss.clear();
			oss.str("");
		}

		infile.close();
	}

	void readMeshMigrate(const char* filename)
	{
		ifstream infile(filename);

		int elid, prid;

		lsize = numvolume;
		elids = new int[lsize];
		procid = new int[lsize];

		for (int i=0; i<numvolume; i++) {

			infile >> elid >> prid;
			elids[i] = elid;
			procid[i] = prid - 1;
		}

		infile.close();
	}

public:
	int      mypid ;
	int      numprocs ;

    MeshMig(int rank, int size, bool useFiles)
	{
		mypid = rank;
		numprocs = size;
		numdirichlet = 0;
        boundset = false;

        if(useFiles == true)
        {
            // Here create the filename char *
            char * boundaryfile = new char[132];
            char * elementfile = new char[132];
            char * headerfile = new char[132];
            char * nodefile = new char[132];
            char * sharedfile = new char[132];
            char * migfile = new char[132];

            /// Here create the partitioning.n directory first. Clean if already exists.

            sprintf(boundaryfile, "/root/partitioning.%d/part.%d.boundary", size, rank+1);
            sprintf(elementfile, "/root/partitioning.%d/part.%d.elements", size, rank+1);
            sprintf(headerfile, "/root/partitioning.%d/part.%d.header", size, rank+1);
            sprintf(nodefile, "/root/partitioning.%d/part.%d.nodes", size, rank+1);
            sprintf(sharedfile, "/root/partitioning.%d/part.%d.shared", size, rank+1);
            sprintf(migfile, "/root/partitioning.%d/part.%d.mig", size, rank+1);

            /*
            std::cout << "BOUNDARY FILE: " << boundaryfile << std::endl;
            std::cout << "ELEMENT FILE: " << elementfile << std::endl;
            std::cout << "HEADER FILE: " << headerfile << std::endl;
            std::cout << "NODE FILE: " << nodefile << std::endl;
            std::cout << "SHARED FILE: " << sharedfile << std::endl;
            */

            readMeshHeader(headerfile);
            readMeshNodes(nodefile);
            readMeshElements(elementfile);
            readMeshBoundary(boundaryfile);
            readMeshShared(sharedfile);
            readMeshMigrate(migfile);
        }
	}

	void insertHeaderData(int _numvert, int _numvol, int _numbound, int _numshared)
    {
        numvertex = _numvert;
        numvolume = _numvol;
        numboundary = _numbound;
        numshared = _numshared;
    }

    void insertMeshNodes(int *glidptr, double *xptr, double *yptr, double *zptr)
    {
        verts.resize(numvertex);

		for (int i=0; i<numvertex; i++) {

			verts[i].coor[0] = xptr[i];
			verts[i].coor[1] = yptr[i];
			verts[i].coor[2] = zptr[i];
			verts[i].idx = glidptr[i];
			verts[i].useCount = 0;
			// verts[i].vtype = 'i';

			vertmap[glidptr[i]] = i;
		}
    }

    void insertMeshElements(int *glidptr, int *xptr, int *yptr, int *zptr, int *wptr)
    {
        volelems.resize(numvolume);

		for (int i=0; i<numvolume; i++) {

			volelems[i].verts[0] = xptr[i];
			volelems[i].verts[1] = yptr[i];
			volelems[i].verts[2] = zptr[i];
			volelems[i].verts[3] = wptr[i];
			volelems[i].idx = glidptr[i];

			// increase usage of vertices.
			verts[vertmap.find(xptr[i])->second].useCount++;
			verts[vertmap.find(yptr[i])->second].useCount++;
			verts[vertmap.find(zptr[i])->second].useCount++;
			verts[vertmap.find(wptr[i])->second].useCount++;

			volmap[glidptr[i]] = i;
		}
    }

    void insertMeshBoundary(int *gfidptr, int *xptr, int *yptr, int *zptr, int *parptr)
    {
        boundelems.resize(numboundary);

		for (int i=0; i<numboundary; i++) {

			boundelems[i].verts[0] = xptr[i];
			boundelems[i].verts[1] = yptr[i];
			boundelems[i].verts[2] = zptr[i];
			boundelems[i].idx = i+1;
			boundelems[i].geoid = gfidptr[i];
			boundelems[i].parentelement = parptr[i];

			volbound[parptr[i]] = i;
		}
    }

    void insertMeshShared(int *vidptr, string *sharedptr)
    {
        for (int i=0; i<numshared; i++) {
			sharedVert[vidptr[i]] = sharedptr[i];
		}
    }

    void insertMigrationData(int ne, int *dest)
    {
        lsize = numvolume;
        elids = new int[lsize];
		procid = new int[lsize];

        for (int i=0; i<ne; i++) {
		    elids[i] = volelems[i].idx;     // By default use the same order for ids as volelems
			procid[i] = dest[i];
		}
    }

    void createElmerOutput()
    {
        int rank = mypid;
        int size = numprocs;

        char * boundaryfile = new char[132];
        char * elementfile = new char[132];
        char * headerfile = new char[132];
        char * nodefile = new char[132];
        char * sharedfile = new char[132];
        char * migfile = new char[132];
        char * partfile = new char[132];

        /// Here create the partitioning.n directory first. Clean if already exists.

        sprintf(boundaryfile, "/root/partitioning.%d/part.%d.boundary", size, rank+1);
        sprintf(elementfile, "/root/partitioning.%d/part.%d.elements", size, rank+1);
        sprintf(headerfile, "/root/partitioning.%d/part.%d.header", size, rank+1);
        sprintf(nodefile, "/root/partitioning.%d/part.%d.nodes", size, rank+1);
        sprintf(sharedfile, "/root/partitioning.%d/part.%d.shared", size, rank+1);
        sprintf(migfile, "/root/partitioning.%d/part.%d.mig", size, rank+1);
        sprintf(partfile, "/root/partitioning.%d/part.%d.part", size, rank+1);

        /// Write header file
        ofstream outheader(headerfile);
        outheader << numvertex << " " << numvolume << " " << numboundary + numdirichlet << endl;

        if(numdirichlet == 0)
        {
            outheader << 2 << endl;
            outheader << "504 " << numvolume << endl;
            outheader << "303 " << numboundary << endl;
        }
        else
        {
            outheader << 3 << endl;
            outheader << "504 " << numvolume << endl;
            outheader << "101 " << numdirichlet << endl;
            outheader << "303 " << numboundary << endl;
        }

        if(numshared != 0)
        {
            outheader << numshared << " 0" << endl;
        }

        outheader.close();

        /// Write node file
        ofstream outnode(nodefile);

        for(int i=0; i<numvertex; i++)
        {
            outnode << verts[i].idx << " -1 " << verts[i].coor[0] << " " << verts[i].coor[1] << " " << verts[i].coor[2] << endl;
        }

        outnode.close();

        /// Write elements file.
        ofstream outelements(elementfile);

        for(int i=0; i<numvolume; i++)
        {
            outelements << volelems[i].idx << " 1 504 " << volelems[i].verts[0] << " " << volelems[i].verts[1] << " " << volelems[i].verts[2] << " " << volelems[i].verts[3] << endl;
        }

        outelements.close();

        /// Write boundary file. No dirichlet creation yet.
        ofstream outbound(boundaryfile);

        for(int i=0; i<numboundary; i++)
        {
            outbound << boundelems[i].idx << " " << boundelems[i].geoid << " " << boundelems[i].parentelement << " " << "0" << " 303 " << boundelems[i].verts[0] << " " << boundelems[i].verts[1] << " " << boundelems[i].verts[2] << endl;
        }

        outbound.close();

        /// Write shared file.
        ofstream outshared(sharedfile);

        map<int,string>::iterator sharit;

        for(sharit = sharedVert.begin(); sharit != sharedVert.end(); sharit++)
        {
            outshared << (*sharit).first << " " << (*sharit).second << std::endl;
        }

        outshared.close();

        /*
        /// Write mig file.
        ofstream outmig(migfile);

        for(int i=0; i<numvolume; i++)
        {
            outmig << elids[i] << " " << procid[i] << std::endl;
        }

        outmig.close();
        */

        /// Write part boundary face file.
        ofstream outpart(partfile);
        int partsize = partbound.size();

        outpart << partsize << std::endl;

        for(int i=0; i<partsize; i++)
        {
            outpart << partbound[i].verts[0] << " " << partbound[i].verts[1] << " " << partbound[i].verts[2] << std::endl;
        }
    }

    void createGeomViewOutput()
    {
        int rank = mypid;
        int size = numprocs;

        char * boundaryfile = new char[132];

        sprintf(boundaryfile, "/root/partitioning.%d/geomview.%d.off", size, rank+1);

        ofstream outbound(boundaryfile);

        outbound << "OFF" << std::endl;
        outbound << numvertex << " " << numboundary << " 0 #" << std::endl;

        map<int,int> tmap;

        for(int i=0; i<numvertex; i++)
        {
            outbound << verts[i].coor[0] << " " << verts[i].coor[1] << " " << verts[i].coor[2] << endl;

            tmap[verts[i].idx] = i;
        }

        for(int i=0; i<numboundary; i++)
        {
            outbound << "3 " << tmap.find(boundelems[i].verts[0])->second << " " << tmap.find(boundelems[i].verts[1])->second << " " << tmap.find(boundelems[i].verts[2])->second << " " << (rank+1)*10 << endl;
        }

        outbound.close();

        if(mypid == 0)
        {
            /// Write the main file.
            char * mainfile = new char[132];

            sprintf(mainfile, "/root/partitioning.%d/geomview.off", size);

            ofstream outmain(mainfile);

            outmain << "LIST" << std::endl;

            for(int i=1; i<=numprocs; i++)
            {
                 outmain << " { < geomview." << i << ".off }" << std::endl;
            }
        }
    }

	void dumpData()
	{
		// Vertices Dump.
		std::cout << "RANK: " << mypid << " Vertices Count:\t\t\t" << numvertex << std::endl;

		for (int i=0; i<numvertex; i++) {

			// if (verts[i].useCount != 0) {
				std::cout << verts[i].idx << " " << verts[i].useCount << std::endl;
			// }
		}

		std::cout << "RANK: " << mypid << " Volume Elements Count:\t\t" << numvolume << std::endl;

		for (int i=0; i<numvolume; i++) {
			//if (!volelems[i].deleted) {
				std::cout << volelems[i].idx << std::endl;
			//}
		}

		// std::cout << "RANK: " << mypid << " Boundary Elements Count:\t" << numboundary << std::endl;

		// Shared Vertices Dump.
		std::cout << "RANK: " << mypid << " Shared Vertices Count:\t\t" << numshared << std::endl;

		map<int,string>::iterator it;
		for ( it=sharedVert.begin() ; it != sharedVert.end(); it++ )
		{
			std::cout << (*it).first << "=>  " << (*it).second << std::endl;
		}
	}

	int mymigrate2(MPI_Comm comm, int n, int *eids, int *dest)
	{
		for (int i=0; i<n; i++) {
			elids[i] = eids[i];
			procid[i] = dest[i] - 1;
		}

		mymigrate2(comm);

		return 0;
	}

	int mymigrate2(MPI_Comm comm, vector<int> & elid, vector<int> & dest)
	{
		for (unsigned int i=0; i<dest.size(); i++) {
			elids[i] = elid[i];
			procid[i] = dest[i] - 1;
		}

		mymigrate2(comm);

		return 0;
	}

	int mymigrate2(MPI_Comm comm, int n, int *dest)
	{
		for (int i=0; i<n; i++) {
		    elids[i] = volelems[i].idx;     // By default use the same order for ids as volelems
			procid[i] = dest[i] - 1;
		}

		mymigrate2(comm);

		return 0;
	}

	int mymigrate2(MPI_Comm comm, vector<int> & dest)
	{
		for (unsigned int i=0; i<dest.size(); i++) {
		    elids[i] = volelems[i].idx;     // By default use the same order for ids as volelems
			procid[i] = dest[i] - 1;
		}

		mymigrate2(comm);

		return 0;
	}

	int findEl(int* tosearch, int eltofind, int count)
	{
	    int toReturn = -1;
        for(int i=0; i<count; i++)
        {
            if(eltofind == tosearch[i])
            {
                toReturn = i;
                break;
            }
        }

        return toReturn;
	}

	void setBoundaryElements()
	{
	    bool (*fn_pt)(Index3,Index3) = fncomp;
        std::multimap<Index3,int, bool(*)(Index3, Index3)> face2vol(fn_pt);
        std::multimap<Index3,int, bool(*)(Index3, Index3)>::iterator myit;

        Index3 i3;
        int l;

        for(int i=0; i<numvolume; i++)
        {
            for (int j = 1; j <= 4; j++)   // loop over faces of tet
            {
                l = 0;
                for (int k = 1; k <= 4; k++)
                {
                    if (k != j)
                    {
                        i3.x[l] = volelems[i].verts[k-1];
                        l++;
                    }
                }

                i3.Sort();

                face2vol.insert(pair<Index3,int>(i3,volelems[i].idx));
            }
        }

        for(int i=0; i<numboundary; i++)
        {
            i3.x[0] = boundelems[i].verts[0];
            i3.x[1] = boundelems[i].verts[1];
            i3.x[2] = boundelems[i].verts[2];

            i3.Sort();

            /*
            if(mypid == 0)
            {
                std::cout << i3.x[0] << " " << i3.x[1] << " " << i3.x[2] << endl;
            }
            */

            myit = face2vol.find(i3);

            if(myit == face2vol.end())
            {
                std::cout << "////////////" << endl;
                std::cout << "////////////" << endl;
                std::cout << "It Happens: " << endl;
                std::cout << "////////////" << endl;
                std::cout << "////////////" << endl;
            }

            boundelems[i].parentelement = (*myit).second;
        }

        /*
        if(mypid == 0)
        {
            for(myit = face2vol.begin(); myit != face2vol.end(); myit++)
            {
                std::cout << (*myit).first.x[0] << " " << (*myit).first.x[1] << " " << (*myit).first.x[2] << " " << (*myit).second << endl;
            }
        }
        */

        face2vol.clear();

        boundset = true;
	}

	int mymigrate2(MPI_Comm comm)
	{
	    if(boundset == false)
	    {
            if(mypid == 0)
            {
                std::cout << "Beginning Parent Element Search" << std::endl;
            }
	        setBoundaryElements();

	        if(mypid == 0)
            {
                std::cout << "Ending Parent Element Search" << std::endl;
            }
	    }

		// Volume Datatype
		MPI_Datatype Volumetype;
		MPI_Datatype type[2] = { MPI_INT, MPI_INT };
		int blocklen[2] = { NUM_VOLM_VERT, 1 };
		MPI_Aint disp[2];
		disp[0] = (MPI_Aint)&migVolume[0].verts - (MPI_Aint)&migVolume[0];
		disp[1] = (MPI_Aint)&migVolume[0].idx - (MPI_Aint)&migVolume[0];
		MPI_Type_create_struct(2, blocklen, disp, type, &Volumetype);
		MPI_Type_commit(&Volumetype);

		// Vertex Datatype
		MPI_Datatype Vertextype;
		MPI_Datatype type2[3] = { MPI_DOUBLE, MPI_INT, MPI_INT };
		int blocklen2[3] = { NUM_VERT_SIZE, 1, 1 };
		MPI_Aint disp2[3];
		disp2[0] = (MPI_Aint)&migVertex[0].coor - (MPI_Aint)&migVertex[0];
		disp2[1] = (MPI_Aint)&migVertex[0].idx - (MPI_Aint)&migVertex[0];
		disp2[2] = (MPI_Aint)&migVertex[0].useCount - (MPI_Aint)&migVertex[0];
		MPI_Type_create_struct(3, blocklen2, disp2, type2, &Vertextype);
		MPI_Type_commit(&Vertextype);

		// Boundary Datatype
		MPI_Datatype Boundtype;
		MPI_Datatype type3[4] = { MPI_INT, MPI_INT, MPI_INT, MPI_INT };
		int blocklen3[4] = { NUM_SURF_VERT, 1, 1, 1 };
		MPI_Aint disp3[4];
		disp3[0] = (MPI_Aint)&migBoundary[0].verts - (MPI_Aint)&migBoundary[0];
		disp3[1] = (MPI_Aint)&migBoundary[0].geoid - (MPI_Aint)&migBoundary[0];
		disp3[2] = (MPI_Aint)&migBoundary[0].parentelement - (MPI_Aint)&migBoundary[0];
		disp3[3] = (MPI_Aint)&migBoundary[0].idx - (MPI_Aint)&migBoundary[0];
		MPI_Type_create_struct(4, blocklen3, disp3, type3, &Boundtype);
		MPI_Type_commit(&Boundtype);

        /// Variables needed for migration (new algorithm) START.
        int vsendcount = 0;          // Number of different processors to be send.
        int vrecvcount = 0;          // Number of different processors to be received.
        int *vertsendprids;     // Processor Id's of vertices to be send.
        int *vertsendcounts;    // Count of vertices to be send.
        int *vertrecvprids;     // Processor Id's of vertices to be received.
        int *vertrecvcounts;    // Count of vertices to be received.

        // the same variables for volume migration.
        int volsendcount = 0;
        int volrecvcount = 0;
        int *volsendprids;
        int *volsendcounts;
        int *volrecvprids;
        int *volrecvcounts;

        int * volumesend2 = new int[numprocs];
        for(int i=0; i<numprocs; i++)
        {
            volumesend2[i] = 0;
        }

        // The same variables for boundary migration.

        int facesendcount = 0;
        int facerecvcount = 0;
        int *facesendprids;
        int *facesendcounts;
        int *facerecvprids;
        int *facerecvcounts;

        int * facesend = new int[numprocs];
        for(int i=0; i<numprocs; i++)
        {
            facesend[i] = 0;
        }

        // The same variables for shared1

        int sharsendcount = 0;
        int sharrecvcount = 0;
        int *sharsendprids;
        int *sharsendcounts;
        int *sharrecvprids;
        int *sharrecvcounts;

        int * sharedsend3 = new int[numprocs];
        for(int i=0; i<numprocs; i++)
        {
            sharedsend3[i] = 0;
        }

        // The same variables for shared2

        int shar2sendcount = 0;
        int shar2recvcount = 0;
        int *shar2sendprids;
        int *shar2sendcounts;
        int *shar2recvprids;
        int *shar2recvcounts;

        int * sharedsend4 = new int[numprocs];
        for(int i=0; i<numprocs; i++)
        {
            sharedsend4[i] = 0;
        }

        /// Variables needed for migration (new algorithm) END.

		int* sendpr;
		int* sendpr2;
		int* sendpr3;
		int* sendpr4;
		int* sendpr5;

		int* recvpr;
		int* recvpr2;
		int* recvpr3;
		int* recvpr4;
		int* recvpr5;

		int *sharedtosend;
		int *sharedtorecv;

		int *sharedtosend2;
		int *sharedtorecv2;

		set<int>::iterator setit;
		map<int,string>::iterator myit;
		Vertex tmpVert2;
		MapI2I::iterator vertit2;
		stringstream ss2;
		int curowner = 0;

		// To store vertices to send to each processor.
		set<int>* procsets;

        { // 0) Send Boundary Elements

            int count = 0;
            int tmp34;
			for (int i=0; i<numboundary; i++)
			{
			    tmp34 = boundelems[i].parentelement;
			    tmp34 = volmap.find(tmp34)->second;

			    if(elids[tmp34] != boundelems[i].parentelement)
                {
                    std::cout << "////////////" << endl;
                    std::cout << "////////////" << endl;
                    std::cout << "It Happens: " << endl;
                    std::cout << "////////////" << endl;
                    std::cout << "////////////" << endl;
                }

				if (procid[tmp34] != mypid) {
					count++;
				}
			}

			migBoundary.resize(count);

			for (int i=0; i<numboundary; i++)
			{
			    tmp34 = boundelems[i].parentelement;
			    tmp34 = volmap.find(tmp34)->second;

				if (procid[tmp34] != mypid) {
					facesend[procid[tmp34]]++;
				}
			}


            /// New migration find send counts START
            facesendcount = 0;
            for(int i=0; i<numprocs; i++)
            {
                if(facesend[i] == 0)
                {
                    continue;
                }
                else
                {
                    facesendcount++;
                }
            }

            facesendprids = new int[facesendcount];
            facesendcounts = new int[facesendcount];

            facesendcount = 0;
            for(int i=0; i<numprocs; i++)
            {
                if(facesend[i] == 0)
                {
                    continue;
                }
                else
                {
                    facesendcounts[facesendcount] = facesend[i];
                    facesend[i] = 1;
                    facesendprids[facesendcount] = i;
                    facesendcount++;
                }
            }

            /// Find receive counts
            {
                int *totdiffrecvs = new int[numprocs];  // total number of different processors to receive from.
                MPI_Allreduce(facesend, totdiffrecvs, numprocs,MPI_INT, MPI_SUM, comm);

                facerecvcount = totdiffrecvs[mypid];
                delete [] totdiffrecvs;

                facerecvcounts = new int[facerecvcount];
                facerecvprids = new int[facerecvcount];
            }

            // std::cout << "RANK: " << mypid << " receives " << volrecvcount << " elements" << std::endl;

            MPI_Request * reqs = new MPI_Request[facerecvcount + facesendcount];
            MPI_Status * stats = new MPI_Status[facerecvcount + facesendcount];

            /// Now send to each one how much you will send.
            for(int i=0; i<facesendcount; i++)
            {
                // Send number of volume elements to be send.
                MPI_Isend(&facesendcounts[i], 1, MPI_INT, facesendprids[i], 123, comm, &reqs[i]);
            }
            for(int i=0; i<facerecvcount; i++)
            {
                MPI_Irecv(&facerecvcounts[i], 1, MPI_INT, MPI_ANY_SOURCE, 123, comm, &reqs[facesendcount+i]);
            }
            MPI_Waitall(facerecvcount + facesendcount, reqs, stats);

            for(int i=0; i<facerecvcount; i++)
            {
                facerecvprids[i] = stats[facesendcount + i].MPI_SOURCE;
            }

            delete [] reqs;
            delete [] stats;

            /// New migration find send and recv counts END

            recvpr5 = new int[facerecvcount+1];
			sendpr5 = new int[facesendcount+1];
			//
			sendpr5[0] = 0;
			for (int j=0; j<facesendcount; j++) {
				sendpr5[j+1] = sendpr5[j] + facesendcounts[j];
			}

			recvpr5[0] = 0;
			for (int i=0; i<facerecvcount; i++) {
				recvpr5[i+1] = recvpr5[i] + facerecvcounts[i];
			}

			int *currentel5 = new int[facesendcount];
			for (int i=0; i<facesendcount; i++) {
				currentel5[i] = 0;
			}

			int procindex = 0;

			list<int> toerase;

			for (int i=0; i<numboundary; i++) {

                tmp34 = boundelems[i].parentelement;
			    tmp34 = volmap.find(tmp34)->second;

				if (procid[tmp34] != mypid) {

                    procindex = findEl(facesendprids, procid[tmp34], facesendcount);
					migBoundary[sendpr5[procindex] + currentel5[procindex]] = boundelems[i];

					currentel5[procindex]++;

                    toerase.push_back(i);
					//boundelems[i] = boundelems[numboundary-1];
					//numboundary--;
					//i--;
				}
			}

            list<int>::iterator myite;
            for(myite = toerase.begin(); myite != toerase.end(); myite++)
            {
                boundelems[(*myite)].idx = -1;
                // numboundary--;
            }

            for(int i=0; i<numboundary; i++)
            {
                if(boundelems[i].idx == -1)
                {
                    if(i != (numboundary-1))
                    {
                        do
                        {
                            boundelems[i] = boundelems[numboundary-1];
                            numboundary--;

                        }while(boundelems[i].idx == -1 && i != (numboundary-1));
                    }
                    else
                    {
                        numboundary--;
                    }
                }
            }

			boundelems.resize(numboundary + recvpr5[facerecvcount]);

			MPI_Barrier(comm);

			// Send Elements
			int tmpp2;

			reqs = new MPI_Request[facesendcount + facerecvcount];
			stats = new MPI_Status[facesendcount + facerecvcount];

			for (int i=0; i<facesendcount; i++) {
                tmpp2 = sendpr5[i+1] - sendpr5[i];
				MPI_Isend(&migBoundary[sendpr5[i]], tmpp2, Boundtype, facesendprids[i], 123, comm, &reqs[i]);
			}

			for (int i=0; i<facerecvcount; i++) {
                tmpp2 = recvpr5[i+1] - recvpr5[i];
                MPI_Irecv(&boundelems[numboundary + recvpr5[i]], tmpp2, Boundtype, facerecvprids[i], 123, comm, &reqs[facesendcount + i]);
			}
			MPI_Waitall(facerecvcount + facesendcount, reqs, stats);

            delete [] reqs;
            delete [] stats;

            MPI_Barrier(comm);

            std::cout << "RANK: " << mypid << " boundary receives ";
            for(int i=0; i<facerecvcount; i++)
            {
                std::cout << facerecvcounts[i] << " - " << facerecvprids[i] << "; ";
            }
            std::cout << std::endl;

            numboundary += recvpr5[facerecvcount];
			MPI_Barrier(comm);

		} // 0) Send Boundary Elements


		{ // 1) Send Volume Elements

			int count = 0;
			for (int i=0; i<lsize; i++) {
				if (procid[i] != mypid) {
					count++;
				}
			}

			migVolume.resize(count);

			for (int i=0; i<lsize; i++) {
				if (procid[i] != mypid) {
					volumesend2[procid[i]]++;
				}
			}


            /// New migration find send counts START
            volsendcount = 0;
            for(int i=0; i<numprocs; i++)
            {
                if(volumesend2[i] == 0)
                {
                    continue;
                }
                else
                {
                    volsendcount++;
                }
            }

            volsendprids = new int[volsendcount];
            volsendcounts = new int[volsendcount];

            volsendcount = 0;
            for(int i=0; i<numprocs; i++)
            {
                if(volumesend2[i] == 0)
                {
                    continue;
                }
                else
                {
                    volsendcounts[volsendcount] = volumesend2[i];
                    volumesend2[i] = 1;
                    volsendprids[volsendcount] = i;
                    volsendcount++;
                }
            }

            /// Find receive counts
            {
                int *totdiffrecvs = new int[numprocs];  // total number of different processors to receive from.
                MPI_Allreduce(volumesend2, totdiffrecvs, numprocs,MPI_INT, MPI_SUM, comm);

                volrecvcount = totdiffrecvs[mypid];
                delete [] totdiffrecvs;

                volrecvcounts = new int[volrecvcount];
                volrecvprids = new int[volrecvcount];
            }

            // std::cout << "RANK: " << mypid << " receives " << volrecvcount << " elements" << std::endl;

            MPI_Request * reqs = new MPI_Request[volrecvcount + volsendcount];
            MPI_Status * stats = new MPI_Status[volrecvcount + volsendcount];

            /// Now send to each one how much you will send.
            for(int i=0; i<volsendcount; i++)
            {
                // Send number of volume elements to be send.
                MPI_Isend(&volsendcounts[i], 1, MPI_INT, volsendprids[i], 123, comm, &reqs[i]);
            }
            for(int i=0; i<volrecvcount; i++)
            {
                MPI_Irecv(&volrecvcounts[i], 1, MPI_INT, MPI_ANY_SOURCE, 123, comm, &reqs[volsendcount+i]);
            }
            MPI_Waitall(volrecvcount + volsendcount, reqs, stats);

            for(int i=0; i<volrecvcount; i++)
            {
                volrecvprids[i] = stats[volsendcount + i].MPI_SOURCE;
            }

            delete [] reqs;
            delete [] stats;

            /// New migration find send and recv counts END

			recvpr2 = new int[volrecvcount+1];
			sendpr2 = new int[volsendcount+1];
			//
			sendpr2[0] = 0;
			for (int j=0; j<volsendcount; j++) {
				sendpr2[j+1] = sendpr2[j] + volsendcounts[j];
			}

			recvpr2[0] = 0;
			for (int i=0; i<volrecvcount; i++) {
				recvpr2[i+1] = recvpr2[i] + volrecvcounts[i];
			}

			int *currentel2 = new int[volsendcount];
			for (int i=0; i<volsendcount; i++) {
				currentel2[i] = 0;
			}

            procsets = new set<int>[volsendcount];

			int volelind = 0;
			int procindex = 0;
			for (int i=0; i<lsize; i++) {
				if (procid[i] != mypid) {

                    procindex = findEl(volsendprids, procid[i], volsendcount);

					volelind = volmap.find(elids[i])->second;
					migVolume[sendpr2[procindex] + currentel2[procindex]] = volelems[volelind];

					// Add to the vertices set to send.
					for (int k=0; k<NUM_VOLM_VERT; k++) {
						procsets[procindex].insert(volelems[volelind].verts[k]);
					}

					currentel2[procindex]++;

					// Remove from the end. Update map.
					volmap.erase(volelems[volelind].idx);
					volelems[volelind] = volelems[numvolume-1];
					volmap[volelems[numvolume-1].idx] = volelind;
					numvolume--;
				}
			}

			volelems.resize(numvolume + recvpr2[volrecvcount]);

			MPI_Barrier(comm);

			// Send Elements
			int tmpp2;

			reqs = new MPI_Request[volsendcount + volrecvcount];
			stats = new MPI_Status[volsendcount + volrecvcount];

			for (int i=0; i<volsendcount; i++) {
                tmpp2 = sendpr2[i+1] - sendpr2[i];
				MPI_Isend(&migVolume[sendpr2[i]], tmpp2, Volumetype, volsendprids[i], 123, comm, &reqs[i]);
			}

			for (int i=0; i<volrecvcount; i++) {
                tmpp2 = recvpr2[i+1] - recvpr2[i];
                MPI_Irecv(&volelems[numvolume + recvpr2[i]], tmpp2, Volumetype, volrecvprids[i], 123, comm, &reqs[volsendcount + i]);
			}
			MPI_Waitall(volrecvcount + volsendcount, reqs, stats);

            delete [] reqs;
            delete [] stats;


            std::cout << "RANK: " << mypid << " receives ";
            for(int i=0; i<volrecvcount; i++)
            {
                std::cout << volrecvcounts[i] << " - " << volrecvprids[i] << "; ";
            }
            std::cout << std::endl;

            // numvolume += recvpr2[volrecvcount];

            // createElmerOutput();

			MPI_Barrier(comm);


		} // 1) Send Volume Elements



		{ // 2) Send Vertices

            /// volrecvcount are equal to vertrecvcount. Since the processor we send vertices are the same as we sent elements.
            vrecvcount = volrecvcount;
			vsendcount = volsendcount;
            vertsendcounts = new int[vsendcount];
            vertsendprids = new int[vsendcount];
            vertrecvcounts = new int[vrecvcount];
            vertrecvprids = new int[vrecvcount];

			for (int i=0; i<vsendcount; i++) {
				vertsendcounts[i] = procsets[i].size();
				if(procsets[i].size() == 0)
				{
				    std::cout << "PROBLEM at RANK: " << mypid << std::endl;
				}
				vertsendprids[i] = volsendprids[i];
			}

            /// New migration START

            MPI_Request * reqs = new MPI_Request[vrecvcount + vsendcount];
            MPI_Status * stats = new MPI_Status[vrecvcount + vsendcount];

            /// Now send to each one how much you will send.
            for(int i=0; i<vsendcount; i++)
            {
                // Send number of volume elements to be send.
                MPI_Isend(&vertsendcounts[i], 1, MPI_INT, vertsendprids[i], 123, comm, &reqs[i]);
            }
            for(int i=0; i<vrecvcount; i++)
            {
                MPI_Irecv(&vertrecvcounts[i], 1, MPI_INT, MPI_ANY_SOURCE, 123, comm, &reqs[vsendcount+i]);
            }
            MPI_Waitall(vrecvcount + vsendcount, reqs, stats);

            for(int i=0; i<vrecvcount; i++)
            {
                vertrecvprids[i] = stats[vsendcount + i].MPI_SOURCE;
            }

            delete [] reqs;
            delete [] stats;

            /// New migration END

			recvpr = new int[vrecvcount+1];
			sendpr = new int[vsendcount+1];
			//
			sendpr[0] = 0;
			for (int j=0; j<vsendcount; j++) {
				sendpr[j+1] = sendpr[j] + vertsendcounts[j];
			}

			recvpr[0] = 0;
			for (int i=0; i<vrecvcount; i++) {
				recvpr[i+1] = recvpr[i] + vertrecvcounts[i];
			}
			// Fill the elements to send.

			migVertex.resize(sendpr[vsendcount]);

			set<int>::iterator it;
			int k=0;
			for (int i=0; i<volsendcount; i++) {

				k = 0;

				for (it=procsets[i].begin(); it != procsets[i].end(); it++)
				{
					migVertex[sendpr[i] + k] = verts[vertmap.find(*it)->second];
					k++;
				}
			}
            // Send elements.

			int ind3 = 0;
			// Decrease use count of verts.
			for (unsigned int i=0; i<migVolume.size(); i++) {

				for (int k=0; k<NUM_VOLM_VERT; k++) {
					ind3 = migVolume[i].verts[k];
					ind3 = vertmap.find(ind3)->second;
					verts[ind3].useCount--;
				}
			}

			// Increase use count of old verts after new additions.
			MapI2I::iterator ite;
			for (unsigned int i=numvolume; i<volelems.size(); i++) {

				for (int k=0; k<NUM_VOLM_VERT; k++)
				{
					ind3 = volelems[i].verts[k];
					ite = vertmap.find(ind3);
					if (ite != vertmap.end())
					{
						ind3 = vertmap.find(ind3)->second;
						verts[ind3].useCount++;
					}
				}
			}
			// Calculate shared vertice information count




			for (int i=0; i<volsendcount; i++)
			{
				for (setit=procsets[i].begin(); setit != procsets[i].end(); setit++)
				{
					myit = sharedVert.find(*setit);
					vertit2 = vertmap.find(*setit);

					tmpVert2 = verts[(*vertit2).second];

					if (myit == sharedVert.end())
					{
						if (tmpVert2.useCount == 0)
						{
							sharedsend3[mypid] += 2;
						}
						else
						{
							sharedsend3[mypid] += 2;
						}
					}
					else
					{
						// Find Owner processor here
						ss2 << (*myit).second;

						ss2 >> curowner;
						ss2 >> curowner;

                        // ss2.clear();
						ss2.str("");

						// Warning. Be aware that curowner = prid + 1

						if (tmpVert2.useCount == 0)
						{
							sharedsend3[curowner - 1] += 2;
						}
						else
						{
							sharedsend3[curowner - 1]++;
						}
					}
				}
			}

			/// New migration find send counts START
            sharsendcount = 0;
            for(int i=0; i<numprocs; i++)
            {
                if(sharedsend3[i] == 0)
                {
                    continue;
                }
                else
                {
                    sharsendcount++;
                }
            }

            sharsendprids = new int[sharsendcount];
            sharsendcounts = new int[sharsendcount];

            sharsendcount = 0;
            for(int i=0; i<numprocs; i++)
            {
                if(sharedsend3[i] == 0)
                {
                    continue;
                }
                else
                {
                    sharsendcounts[sharsendcount] = sharedsend3[i];
                    sharedsend3[i] = 1;
                    sharsendprids[sharsendcount] = i;
                    sharsendcount++;
                }
            }

            MPI_Barrier(comm);

            /// Find receive counts
            {
                int *totdiffrecvs = new int[numprocs];  // total number of different processors to receive from.
                MPI_Allreduce(sharedsend3, totdiffrecvs, numprocs,MPI_INT, MPI_SUM, comm);

                sharrecvcount = totdiffrecvs[mypid];
                delete [] totdiffrecvs;

                sharrecvcounts = new int[sharrecvcount];
                sharrecvprids = new int[sharrecvcount];
            }

            reqs = new MPI_Request[sharrecvcount + sharsendcount];
            stats = new MPI_Status[sharrecvcount + sharsendcount];

            /// Now send to each one how much you will send.
            for(int i=0; i<sharsendcount; i++)
            {
                // Send number of volume elements to be send.
                MPI_Isend(&sharsendcounts[i], 1, MPI_INT, sharsendprids[i], 123, comm, &reqs[i]);
            }
            for(int i=0; i<sharrecvcount; i++)
            {
                MPI_Irecv(&sharrecvcounts[i], 1, MPI_INT, MPI_ANY_SOURCE, 123, comm, &reqs[sharsendcount+i]);
            }
            MPI_Waitall(sharrecvcount + sharsendcount, reqs, stats);

            for(int i=0; i<sharrecvcount; i++)
            {
                sharrecvprids[i] = stats[sharsendcount + i].MPI_SOURCE;
            }

            delete [] reqs;
            delete [] stats;

            MPI_Barrier(comm);

            /// New migration find send and recv counts END

			// Fill in shared element information data.

			// Shared Elements
			recvpr3 = new int[sharrecvcount+1];
			sendpr3 = new int[sharsendcount+1];

            sendpr3[0] = 0;
			for (int j=0; j<sharsendcount; j++) {
				sendpr3[j+1] = sendpr3[j] + sharsendcounts[j];
			}

			recvpr3[0] = 0;
			for (int i=0; i<sharrecvcount; i++) {
				recvpr3[i+1] = recvpr3[i] + sharrecvcounts[i];
			}

			sharedtosend = new int[3*sendpr3[sharsendcount]];
			sharedtorecv = new int[3*recvpr3[sharrecvcount]];

			int *currentel3 = new int[sharsendcount];
			for (int i=0; i<sharsendcount; i++) {
				currentel3[i] = 0;
			}

            int maxallowed = sendpr3[sharsendcount];
			int cur3 = 0;
			int procind = 0;

            std::cout << "RANK: " << mypid << " has Max loopsize: " << vsendcount << "_" << volsendcount << std::endl;

			for (int i=0; i<volsendcount; i++)
			{
			    // std::cout << volsendprids[i] << std::endl;
				for (setit=procsets[i].begin(); setit != procsets[i].end(); setit++)
				{
					myit = sharedVert.find(*setit);
					vertit2 = vertmap.find(*setit);

					tmpVert2 = verts[(*vertit2).second];

					if (myit == sharedVert.end())    // If internal node.
					{
					    procind = findEl(sharsendprids, mypid, sharsendcount);

					    if(procind == -1 || procind == 16)
					    {
					        std::cout << "Program Patates!" << std::endl;
					        return 0;
					    }

                        if((sendpr3[procind] + currentel3[procind]) > maxallowed)
                        {
                            std::cout << "Stuff Happens!- " << (sendpr3[procind] + currentel3[procind]) << "--" << maxallowed << std::endl;
                            return 0;
                        }

						cur3 = 3 * (sendpr3[procind] + currentel3[procind]);

						if (tmpVert2.useCount == 0)
						{
							// If only sent to one processor do nothing.
							// If send to multiple processors create a new entry and send it to the new processors.
							// std::cout << "NOTHING" << " (2) ";

							// sharedsend[stride + mypid]+=2;	// 1 tane send gonder. 1 tane de delete gonder. Kendine gonder.

							// Delete one from itself.
							sharedtosend[cur3] = tmpVert2.idx;
							sharedtosend[cur3+1] = mypid;
							sharedtosend[cur3+2] = 0;	// 0 means delete 1 means add.

							//std::cout << tmmm << " " << tmpVert2.idx << " " << mypid << " 0 A" << std::endl;

							// Add one to curselectedproc
							sharedtosend[cur3+3] = tmpVert2.idx;
							sharedtosend[cur3+4] = volsendprids[i];
							sharedtosend[cur3+5] = 1;	// 0 means delete 1 means add.

							//std::cout << tmmm << " " << tmpVert2.idx << " " << curselectedproc << " 1 A" << std::endl;

							currentel3[procind]+=2;

						}
						else
						{
							// Create a new entry and send it to the new processors.
							// sharedsend[stride + mypid] += 2;		// 2 tane send gonder.
							// std::cout << "NOTHING" << " (2) ";

							// add one itself.
							sharedtosend[cur3] = tmpVert2.idx;
							sharedtosend[cur3+1] = mypid;
							sharedtosend[cur3+2] = 1;	// 0 means delete 1 means add.

							//std::cout << tmmm << " " << tmpVert2.idx << " " << mypid << " 1 B" << std::endl;

							// add one send processor
							sharedtosend[cur3+3] = tmpVert2.idx;
							sharedtosend[cur3+4] = volsendprids[i];
							sharedtosend[cur3+5] = 1;	// 0 means delete 1 means add.

							//std::cout << tmmm << " " << tmpVert2.idx << " " << curselectedproc << " 1 B" << std::endl;

							currentel3[procind]+=2;
						}

						// Hold a copy of it.
					}
					else // MISSING find owner prid.
					{
						// Find Owner processor here
						ss2 << (*myit).second;

						ss2 >> curowner;
						ss2 >> curowner;

                    	// ss2.clear();
						ss2.str("");

                        procind = findEl(sharsendprids, curowner-1, sharsendcount);

                        if(procind == -1 || procind == 16)
					    {
					        std::cout << "Program Patates!" << std::endl;
					        return 0;
					    }

                        if((sendpr3[procind] + currentel3[procind]) > maxallowed)
                        {
                            std::cout << "Stuff Happens!- " << (sendpr3[procind] + currentel3[procind]) << "--" << maxallowed << std::endl;
                            return 0;
                        }

						cur3 = 3 * (sendpr3[procind] + currentel3[procind]);

						// Warning. Be aware that curowner = prid + 1

						if (tmpVert2.useCount == 0)
						{
							// Send delete information of self.
							// Send add information of the sent processor.

							// std::cout << curowner << " (2) ";

							// sharedsend[stride + curowner - 1] += 2;			// vertex sahibine 1 adet delete 1 adet ekle gonder.

							// Delete one from itself.
							sharedtosend[cur3] = tmpVert2.idx;
							sharedtosend[cur3+1] = mypid;
							sharedtosend[cur3+2] = 0;	// 0 means delete 1 means add.

							//std::cout << tmmm << " " << tmpVert2.idx << " " << mypid << " 0 C" << std::endl;

							// Add one to curselectedproc
							sharedtosend[cur3+3] = tmpVert2.idx;
							sharedtosend[cur3+4] = volsendprids[i];
							sharedtosend[cur3+5] = 1;	// 0 means delete 1 means add.

							//std::cout << tmmm << " " << tmpVert2.idx << " " << curselectedproc << " 1 C" << std::endl;

							currentel3[procind]+=2;

						}
						else
						{
							// Send add information of the sent processor.
							// std::cout << "NOTHING" << " (1) ";
							// sharedsend[stride + curowner - 1]++;			// vertex sahibine 1 adet ekle gonder.

							// Add one to curselectedproc
							sharedtosend[cur3] = tmpVert2.idx;
							sharedtosend[cur3+1] = volsendprids[i];
							sharedtosend[cur3+2] = 1;	// 0 means delete 1 means add.

							//std::cout << tmmm << " " << tmpVert2.idx << " " << curselectedproc << " 1 D" << std::endl;

							currentel3[procind]++;

						}

					}
				}
				// std::cout << "RANK: " << mypid << " reached loop number: " << i << std::endl;
			}

            MPI_Barrier(comm);

            /// BURAYA KADAR KONTROLLU
			// Delete verts with use count = 0;

			for (int i=0; i<numvertex; i++)
			{
				if (verts[i].useCount == 0)
				{
					// std::cout << "Erasing " << verts[i].idx << " and replacing index of " << verts[numvertex-1].idx << " with " << i << std::endl;

					vertmap.erase(verts[i].idx);
					verts[i] = verts[numvertex-1];
					vertmap[verts[i].idx] = i;

					numvertex--;
					i--;					// Makes sure to make calculation for the new element.
				}
			}

            MPI_Barrier(comm);

			verts.resize(numvertex + recvpr[vrecvcount]);




            // Send volume elements
            reqs = new MPI_Request[vsendcount + vrecvcount];
			stats = new MPI_Status[vsendcount + vrecvcount];

            MPI_Barrier(comm);
			int tmpp2;
			for (int i=0; i<vsendcount; i++) {
				tmpp2 = sendpr[i+1] - sendpr[i];

				MPI_Isend(&migVertex[sendpr[i]], tmpp2, Vertextype, vertsendprids[i], 123, comm, &reqs[i]);
            }

			for (int i=0; i<vrecvcount; i++) {
				tmpp2 = recvpr[i+1] - recvpr[i];

				MPI_Irecv(&verts[numvertex + recvpr[i]], tmpp2, Vertextype, vertrecvprids[i], 123, comm, &reqs[vsendcount+i]);
			}
			MPI_Waitall(vrecvcount + vsendcount, reqs, stats);

			delete [] reqs;
            delete [] stats;

            /*
            std::cout << "RANK: " << mypid << " receives vertices: ";
            for(int i=0; i<vrecvcount; i++)
            {
                std::cout << vertrecvcounts[i] << " - " << vertrecvprids[i] << "; ";
            }
            std::cout << std::endl;
            */

			MPI_Barrier(comm);

			// New vertices are added.

			/// BURDA KALDIM.

			ind3 = 0;
			int curmax = verts.size();
			int max = curmax;

			// Set count to default.
			for (int i=numvertex; i<max; i++) {
				verts[i].useCount = 0;
			}

			// Delete duplicate elements.
			for (int i=numvertex; i<curmax; i++) {
				ind3 = verts[i].idx;
				ite = vertmap.find(ind3);

				if (ite != vertmap.end()) {

					verts[i] = verts[curmax-1];
					i--;
					curmax--;
				}
				else
				{
					vertmap[verts[i].idx] = i;
				}
			}

			// Update newly added vertice counts
			for (unsigned int i=numvolume; i<volelems.size(); i++) {

				for (int k=0; k<NUM_VOLM_VERT; k++)
				{
					ind3 = volelems[i].verts[k];
					ite = vertmap.find(ind3);
					if (ite != vertmap.end())
					{
						ind3 = vertmap.find(ind3)->second;
						if (ind3 >= numvertex) {
							verts[ind3].useCount++;
						}
					}
				}
			}

			numvertex = curmax;
			numvolume += recvpr2[volrecvcount];

			verts.resize(numvertex);

            MPI_Barrier(comm);
		} // 2) Send Vertices

		{ // 3) Shared Vertices



			int tmpp2;

            int tosub = 0;

            /*
            for(int i=0; i<sharsendcount; i++)
            {
                if(sharsendprids[i] == mypid)
                    tosub = 2;
            }

            int cntr = 0;
            int cntr2 = 0;

            for(int i=0; i<sharsendcount; i++)
            {
                if(sharsendprids[i] == mypid)
                    cntr++;
            }

            for(int i=0; i<sharrecvcount; i++)
            {
                if(sharrecvprids[i] == mypid)
                    cntr2++;
            }

            std::cout << "RANK: " << mypid << " sends to: ";
            for(int i=0; i<sharsendcount; i++)
            {
                std::cout << sharsendprids[i] << "-" << sharsendcounts[i] << " ";
            }
            std::cout << std::endl;
            MPI_Barrier(comm);

            std::cout << "RANK: " << mypid << " receives from: ";
            for(int i=0; i<sharrecvcount; i++)
            {
                std::cout << sharrecvprids[i] << "-" << sharrecvcounts[i] << " ";
            }
            std::cout << std::endl;
            MPI_Barrier(comm);
            */

            // std::cout << "RANK "<< mypid << " tosub: " << tosub << " and sharsendcount: " << sharsendcount << " and sharrecvcount: " << sharrecvcount << std::endl;
            // std::cout << "RANK "<< mypid << " and shars: " << cntr << " and recvcount: " << cntr2 << std::endl;

            MPI_Barrier(comm);

            MPI_Request *reqs = new MPI_Request[sharrecvcount + sharsendcount - tosub];
            MPI_Status *stats = new MPI_Status[sharrecvcount + sharsendcount - tosub];

            /*
            int i2 = 0;

			for (int i=0; i<sharsendcount; i++) {
				tmpp2 = sendpr3[i+1] - sendpr3[i];

                if (sharsendprids[i] != mypid) {
                    std::cout << "Sending " << mypid << " to " << sharsendprids[i] << std::endl;
                    MPI_Isend(&sharedtosend[sendpr3[i]*3], tmpp2*3, MPI_INT, sharsendprids[i], 124, comm, &reqs[i2]);
                    i2++;
                }
                else {
                    for (int k=0; k<tmpp2*3; k++) {
                        sharedtorecv[recvpr3[i]*3 + k] = sharedtosend[sendpr3[i]*3 + k];
                    }
                }
			}

			for (int i=0; i<sharrecvcount; i++) {
				tmpp2 = recvpr3[i+1] - recvpr3[i];

                if (sharrecvprids[i] != mypid) {
                    std::cout << "Receiving " << mypid << " from " << sharrecvprids[i] << std::endl;
                    MPI_Irecv(&sharedtorecv[recvpr3[i]*3], tmpp2*3, MPI_INT, sharrecvprids[i], 124, comm, &reqs[i2]);
                    i2++;
                }
			}
			*/

			for (int i=0; i<sharsendcount; i++) {
				tmpp2 = sendpr3[i+1] - sendpr3[i];

                MPI_Isend(&sharedtosend[sendpr3[i]*3], tmpp2*3, MPI_INT, sharsendprids[i], 124, comm, &reqs[i]);

			}

			for (int i=0; i<sharrecvcount; i++) {
				tmpp2 = recvpr3[i+1] - recvpr3[i];

                MPI_Irecv(&sharedtorecv[recvpr3[i]*3], tmpp2*3, MPI_INT, sharrecvprids[i], 124, comm, &reqs[sharsendcount + i]);
			}

			MPI_Waitall(sharsendcount + sharrecvcount - tosub, reqs, stats);

			delete [] reqs;
			delete [] stats;

			MPI_Barrier(comm);


			/*
			std::cout << "RANK: " << mypid << std::endl;

			for (int i=0; i<recvpr3[numprocs]; i++) {
				std::cout << sharedtorecv[3*i+0] << " " << sharedtorecv[3*i+1] << " " << sharedtorecv[3*i+2] << std::endl;
			}
			 */

            MapI2I listmap;

			int vid2;
			int procid2;
			int boo;
			int cur5 = 0;

			//std::cout << "RANK: " << mypid << std::endl;

			for (int i=0; i<recvpr3[sharrecvcount]; i++) {

				vid2 = sharedtorecv[3*i];
				// procid2 = sharedtorecv[3*i+1];		// procid burda MPI ranke esit. ELmer formatindaki gibi degil.
				// boo = sharedtorecv[3*i+2];

				if (listmap.find(vid2) == listmap.end())
				{
					listmap[vid2] = cur5;
					cur5++;
				}

				// std::cout << vid2 << " " << procid2 << " " << boo << std::endl;
			}

			MapI2I::iterator listmapit;

			/*
			std::cout << "RANK: " << mypid << std::endl;
			 for (listmapit = listmap.begin(); listmapit != listmap.end(); listmapit++) {
			 std::cout << (*listmapit).first << " => " << (*listmapit).second << std::endl;
			 }*/


			vector<int> * mylist = new vector<int>[listmap.size()];
			int tmm33;
			bool found2 = false;
			int foundindex = 0;

			// First fill from already existing data.
			int inn2;
			int inn3;
			map<int,string>::iterator sharit;
			stringstream oss;
			string tmpstr;

			int tt = 0;
			int tt2 = 0;	// count;
			int tt3 = 0;	// temp storage.

			//std::cout << "RANK: " << mypid << std::endl;
			for (listmapit = listmap.begin(); listmapit != listmap.end(); listmapit++)
			{
				inn2 = (*listmapit).second;
				inn3 = (*listmapit).first;

				sharit = sharedVert.find(inn3);

				if (sharit != sharedVert.end()) {
					tmpstr = (*sharit).second;
					tt = (*sharit).first;

					oss << tmpstr;
					oss >> tt2;

					//std::cout << tt << " - ";

					for (int i=0; i<tt2; i++) {
						oss >> tt3;
						//std::cout << tt3-1 << " ";
						mylist[inn2].push_back(tt3-1);
					}

					//std::cout << std::endl;

					oss.clear();
					oss.str("");
				}

			}

			// So far so good.


			for (int i=0; i<recvpr3[sharrecvcount]; i++) {
				vid2 = sharedtorecv[3*i];
				procid2 = sharedtorecv[3*i+1];		// procid burda MPI ranke esit. ELmer formatindaki gibi degil.
				boo = sharedtorecv[3*i+2];

				listmapit = listmap.find(vid2);
				tmm33 = (*listmapit).second;

				if (boo == 1)
				{
					// Add element

					found2 = false;
					// Check whether procid2 exists already.
					for (unsigned int j=0; j<mylist[tmm33].size(); j++) {
						if (mylist[tmm33][j] == procid2) {
							found2 = true;
						}
					}

					if (!found2) {
						mylist[tmm33].push_back(procid2);
					}
				}
				else if(boo == 0)
				{
					// Remove element.

					found2 = false;

					// Check whether procid2 exists.
					for (unsigned int j=0; j<mylist[tmm33].size(); j++) {
						if (mylist[tmm33][j] == procid2) {
							found2 = true;
							foundindex = j;
						}
					}

					if (found2) {
						mylist[tmm33].erase(mylist[tmm33].begin()+foundindex);
					}
				}
			}

           /*
			 std::cout << "RANK: " << mypid << std::endl;

			 for (listmapit = listmap.begin(); listmapit != listmap.end(); listmapit++)
			 {
			 inn2 = (*listmapit).second;

			 std::cout << (*listmapit).first << ": ";

			 for (int j=0; j<mylist[inn2].size(); j++) {
			 std::cout << mylist[inn2][j] << " ";
			 }
			 std::cout << std::endl;
			 }*/

			// Delete sharedvertices which are no longer present here.

			map<int, string>::iterator sharit2;
			MapI2I::iterator vertit3;

			list<int> eraselist;

			//std::cout << "RANK: " << mypid << std::endl;
			for (sharit2 = sharedVert.begin(); sharit2 != sharedVert.end(); sharit2++)
			{
				vertit3 = vertmap.find((*sharit2).first);

				// std::cout << (*vertit3).second << " " << (*vertit3).first << " " << verts[(*vertit3).second].useCount << std::endl;

				if (vertit3 == vertmap.end())
				{
					eraselist.push_back((*sharit2).first);

					// sharedVert.erase(sharit2);
					// sharedVert.erase((*vertit3).first);
					numshared--;
				}
			}

			list<int>::iterator eraselistit;

			for (eraselistit = eraselist.begin(); eraselistit != eraselist.end(); eraselistit++) {
				// std::cout << "RANK: " << mypid << " erasing " << (*eraselistit) << std::endl;
				sharedVert.erase((*eraselistit));
			}

			// 11) Yenilenmis shared vertex listesi olustur.
			// std::cout << "RANK11) : " << mypid << std::endl;

			for (listmapit = listmap.begin(); listmapit != listmap.end(); listmapit++)
			{
				inn2 = (*listmapit).second;

				// std::cout << (*listmapit).first << ": ";

				for (unsigned int j=0; j<mylist[inn2].size(); j++) {
					// std::cout << mylist[inn2][j] << " ";

					sharedsend4[mylist[inn2][j]] += (mylist[inn2].size() + 2);

				}
				// std::cout << std::endl;
			}
            /*
			 std::cout << "RANK: " << mypid << " - ";
			 for (int i=0; i<numprocs; i++) {
			 std::cout << sharedsend2[stride+i] << " ";
			 }
			 std::cout << std::endl;
			 */
            /// New migration find send counts START
            shar2sendcount = 0;
            for(int i=0; i<numprocs; i++)
            {
                if(sharedsend4[i] == 0)
                {
                    continue;
                }
                else
                {
                    shar2sendcount++;
                }
            }

            shar2sendprids = new int[shar2sendcount];
            shar2sendcounts = new int[shar2sendcount];

            shar2sendcount = 0;
            for(int i=0; i<numprocs; i++)
            {
                if(sharedsend4[i] == 0)
                {
                    continue;
                }
                else
                {
                    shar2sendcounts[shar2sendcount] = sharedsend4[i];
                    sharedsend4[i] = 1;
                    shar2sendprids[shar2sendcount] = i;
                    shar2sendcount++;
                }
            }

            MPI_Barrier(comm);



            /// Find receive counts
            {
                int *totdiffrecvs = new int[numprocs];  // total number of different processors to receive from.
                MPI_Allreduce(sharedsend4, totdiffrecvs, numprocs,MPI_INT, MPI_SUM, comm);

                shar2recvcount = totdiffrecvs[mypid];
                delete [] totdiffrecvs;

                shar2recvcounts = new int[shar2recvcount];
                shar2recvprids = new int[shar2recvcount];
            }



            reqs = new MPI_Request[shar2recvcount + shar2sendcount];
            stats = new MPI_Status[shar2recvcount + shar2sendcount];

            /// Now send to each one how much you will send.
            for(int i=0; i<shar2sendcount; i++)
            {
                // Send number of volume elements to be send.
                MPI_Isend(&shar2sendcounts[i], 1, MPI_INT, shar2sendprids[i], 123, comm, &reqs[i]);
            }
            for(int i=0; i<shar2recvcount; i++)
            {
                MPI_Irecv(&shar2recvcounts[i], 1, MPI_INT, MPI_ANY_SOURCE, 123, comm, &reqs[shar2sendcount+i]);
            }
            MPI_Waitall(shar2recvcount + shar2sendcount, reqs, stats);

            for(int i=0; i<shar2recvcount; i++)
            {
                shar2recvprids[i] = stats[shar2sendcount + i].MPI_SOURCE;
            }

            delete [] reqs;
            delete [] stats;

            MPI_Barrier(comm);



            /// New migration find send and recv counts END

			// Shared Elements
			recvpr4 = new int[shar2recvcount+1];
			sendpr4 = new int[shar2sendcount+1];
			//
			sendpr4[0] = 0;
			for (int j=0; j<shar2sendcount; j++) {
				sendpr4[j+1] = sendpr4[j] + shar2sendcounts[j];
			}

			recvpr4[0] = 0;
			for (int i=0; i<shar2recvcount; i++) {
				recvpr4[i+1] = recvpr4[i] + shar2recvcounts[i];
			}

			sharedtosend2 = new int[sendpr4[shar2sendcount]];
			sharedtorecv2 = new int[recvpr4[shar2recvcount]+1];
			sharedtorecv2[recvpr4[shar2recvcount]] = -1;


			// Now fill the elements to send.

			int *currentel4 = new int[shar2sendcount];
			for (int i=0; i<shar2sendcount; i++) {
				currentel4[i] = 0;
			}

			int cur6;
			int prind = 0;

			//std::cout << "RANK: " << mypid << std::endl;
			for (listmapit = listmap.begin(); listmapit != listmap.end(); listmapit++)
			{
				inn2 = (*listmapit).second;

				//std::cout << (*listmapit).first << ": ";

				for (unsigned int j=0; j<mylist[inn2].size(); j++) {
					//std::cout << mylist[inn2][j] << " ";

					prind = findEl(shar2sendprids, mylist[inn2][j], shar2sendcount);

					cur6 = sendpr4[prind] + currentel4[prind];

					sharedtosend2[cur6] = -1;
					sharedtosend2[cur6+1] = (*listmapit).first;

					for (unsigned int k=0; k<mylist[inn2].size(); k++) {
						sharedtosend2[cur6+2+k] = mylist[inn2][k];
					}

					currentel4[prind] += (mylist[inn2].size() + 2);
				}
				//std::cout << std::endl;
			}

			// Send sharedsend2.
            MPI_Barrier(comm);
			int tmpp4;

        /*
			for (int i=0; i<numprocs; i++) {
				tmpp4 = sendpr4[i+1] - sendpr4[i];

				if (tmpp4 != 0) {
					if (i != mypid) {
						// std::cout << "Sending from " << mypid << " to " << i << std::endl;
						MPI_Isend(&sharedtosend2[sendpr4[i]], tmpp4, MPI_INT, i, 123, comm, &request);
					}
					else {
						for (int k=0; k<tmpp4; k++) {
							sharedtorecv2[recvpr4[i] + k] = sharedtosend2[sendpr4[i] + k];
						}
					}

				}
			}

			for (int i=0; i<numprocs; i++) {
				tmpp4 = recvpr4[i+1] - recvpr4[i];

				if (tmpp4 != 0) {
					if (i != mypid) {
						// std::cout << "Receiving " << mypid << " from " << i << std::endl;
						MPI_Irecv(&sharedtorecv2[recvpr4[i]], tmpp4, MPI_INT, i, 123, comm, &request2);
					}
				}
			}
			*/

			reqs = new MPI_Request[shar2recvcount + shar2sendcount];
            stats = new MPI_Status[shar2recvcount + shar2sendcount];

            for (int i=0; i<shar2sendcount; i++) {
				tmpp4 = sendpr4[i+1] - sendpr4[i];

                MPI_Isend(&sharedtosend2[sendpr4[i]], tmpp4, MPI_INT, shar2sendprids[i], 124, comm, &reqs[i]);

			}

			for (int i=0; i<shar2recvcount; i++) {
				tmpp4 = recvpr4[i+1] - recvpr4[i];

                MPI_Irecv(&sharedtorecv2[recvpr4[i]], tmpp4, MPI_INT, shar2recvprids[i], 124, comm, &reqs[shar2sendcount + i]);
			}

            MPI_Waitall(shar2recvcount + shar2sendcount, reqs, stats);

			delete [] reqs;
			delete [] stats;

			MPI_Barrier(comm);

			// Now check what is send.
			/*
			 std::cout << "RANK: " << mypid << std::endl;
			 for (int i=0; i<numprocs; i++) {
			 for (int j=0; j<numprocs; j++) {
			 std::cout << sharedsend2[i*numprocs + j] << "\t";
			 }
			 std::cout << std::endl;
			 }
			 */

			int vid3;
			int proid;
			int kounter = 0;
			stringstream oss2;
			map<int,string>::iterator stit;

			// std::cout << "RANK: " << mypid << std::endl;
			// std::cout << "MAX ELEM: " << recvpr4[numprocs] << std::endl;

 			for (int i=0; i<recvpr4[shar2recvcount]; )
			{
				kounter = 0;
				if (sharedtorecv2[i] == -1)
				{
					vid3 = sharedtorecv2[i+1];

					// std::cout << vid3 << " -";

					do
					{
						proid = sharedtorecv2[i+2+kounter];
						oss << " " << proid;

						kounter++;
					}
					while (sharedtorecv2[i+2+kounter] != -1);

					oss2 << kounter << oss.str();

					stit = sharedVert.find(vid3);

					if (stit == sharedVert.end())
					{
						if (kounter > 1) {
							// Add to the shared vertex list.

							// std::cout << "RANK: " << mypid << " inserts " << vid3 << "-" << oss2.str() << " to the list" << std::endl;

							sharedVert[vid3] = oss2.str();
							numshared++;
						}
						else {
							// std::cout << "RANK: " << mypid << " does nothing with " << vid3 << std::endl;

							// If not currently inside the shared list. And kounter == 1
							// do nothing.
						}
					}
					else
					{
						// Shared element already exists. So update it.
						if (kounter > 1) {	// If still shared update
							// std::cout << "RANK: " << mypid << " updates " << vid3 << " with " << oss2.str() << std::endl;

							sharedVert[vid3] = oss2.str();
						}
						else {				// If not shared anymore delete.
							// std::cout << "RANK: " << mypid << " removes " << vid3 << std::endl;

							sharedVert.erase(vid3);
							numshared--;
						}
					}

					// std::cout << oss2.str();
				}

				i = (i + 2 + kounter);

				// std::cout << std::endl;

				oss.clear();
				oss.str("");
				oss2.clear();
				oss2.str("");
			}

		} // 3) Shared Vertices

        std::cout << "RANK: " << mypid << " has reached this position!" << std::endl;
        MPI_Barrier(comm);

		return 0;
	}

    void refineParallel(MPI_Comm comm, const char* inputFile)
    {
        bool (*fn_pt)(Index2,Index2) = fncomp;
        bool (*fn_pt2)(Index2,Index2) = fncomp;
        std::map<Index2,int, bool(*)(Index2, Index2)> edges(fn_pt);
        std::map<Index2,int, bool(*)(Index2, Index2)> edgesonbound(fn_pt2);
        std::map<Index2,int, bool(*)(Index2, Index2)>::iterator myit;

        Index2 i2;
        int isshared = 0;       // 0 for not shared.  1 for shared.
        int totshared = 0;
        int totnonshared = 0;
        int totall = 0;

        int maxelid = 0;
        int maxglobalid = 0;

        /// Fill edgesonbound
        int curgeoid = 0;

        for(int i=0; i<numboundary; i++)
        {
            curgeoid = boundelems[i].geoid - 1; // ElmerID's are (NetgenID + 1)

            // 1
            i2.x[0] = boundelems[i].verts[0];
            i2.x[1] = boundelems[i].verts[1];

            i2.Sort();

            // If edge is not already in the edge list the add it to the list.
            if(edgesonbound.find(i2) == edgesonbound.end())
            {
                edgesonbound.insert(pair<Index2,int>(i2,curgeoid));
            }

            // 2
            i2.x[0] = boundelems[i].verts[0];
            i2.x[1] = boundelems[i].verts[2];

            i2.Sort();

            if(edgesonbound.find(i2) == edgesonbound.end())
            {
                edgesonbound.insert(pair<Index2,int>(i2,curgeoid));
            }

            // 3
            i2.x[0] = boundelems[i].verts[1];
            i2.x[1] = boundelems[i].verts[2];

            i2.Sort();

            if(edgesonbound.find(i2) == edgesonbound.end())
            {
                edgesonbound.insert(pair<Index2,int>(i2,curgeoid));
            }
        }
        /// Fill edgesobbound END

        for(int i=0; i<numvolume; i++)
        {
            if(volelems[i].verts[0] > maxelid)
            {
                maxelid = volelems[i].verts[0];
            }

            if(volelems[i].verts[1] > maxelid)
            {
                maxelid = volelems[i].verts[1];
            }

            if(volelems[i].verts[2] > maxelid)
            {
                maxelid = volelems[i].verts[2];
            }

            if(volelems[i].verts[3] > maxelid)
            {
                maxelid = volelems[i].verts[3];
            }

            // 1
            i2.x[0] = volelems[i].verts[0];
            i2.x[1] = volelems[i].verts[1];

            i2.Sort();

            if(sharedVert.find(i2.x[0]) == sharedVert.end() || sharedVert.find(i2.x[1]) == sharedVert.end())
            {
                isshared = 0;
            }
            else
            {
                isshared = -1;
            }

            // If edge is not already in the edge list the add it to the list.
            if(edges.find(i2) == edges.end())
            {
                edges.insert(pair<Index2,int>(i2,isshared));

                if(isshared == 0)
                {
                    totnonshared++;
                }
                else
                {
                    totshared++;
                }
            }

            // 2
            i2.x[0] = volelems[i].verts[0];
            i2.x[1] = volelems[i].verts[2];

            i2.Sort();

            if(sharedVert.find(i2.x[0]) == sharedVert.end() || sharedVert.find(i2.x[1]) == sharedVert.end())
            {
                isshared = 0;
            }
            else
            {
                isshared = -1;
            }

            // If edge is not already in the edge list the add it to the list.
            if(edges.find(i2) == edges.end())
            {
                edges.insert(pair<Index2,int>(i2,isshared));

                if(isshared == 0)
                {
                    totnonshared++;
                }
                else
                {
                    totshared++;
                }
            }

            // 3
            i2.x[0] = volelems[i].verts[0];
            i2.x[1] = volelems[i].verts[3];

            i2.Sort();

            if(sharedVert.find(i2.x[0]) == sharedVert.end() || sharedVert.find(i2.x[1]) == sharedVert.end())
            {
                isshared = 0;
            }
            else
            {
                isshared = -1;
            }

            // If edge is not already in the edge list the add it to the list.
            if(edges.find(i2) == edges.end())
            {
                edges.insert(pair<Index2,int>(i2,isshared));

                if(isshared == 0)
                {
                    totnonshared++;
                }
                else
                {
                    totshared++;
                }
            }

            // 4
            i2.x[0] = volelems[i].verts[1];
            i2.x[1] = volelems[i].verts[2];

            i2.Sort();

            if(sharedVert.find(i2.x[0]) == sharedVert.end() || sharedVert.find(i2.x[1]) == sharedVert.end())
            {
                isshared = 0;
            }
            else
            {
                isshared = -1;
            }

            // If edge is not already in the edge list the add it to the list.
            if(edges.find(i2) == edges.end())
            {
                edges.insert(pair<Index2,int>(i2,isshared));

                if(isshared == 0)
                {
                    totnonshared++;
                }
                else
                {
                    totshared++;
                }
            }

            // 5
            i2.x[0] = volelems[i].verts[1];
            i2.x[1] = volelems[i].verts[3];

            i2.Sort();

            if(sharedVert.find(i2.x[0]) == sharedVert.end() || sharedVert.find(i2.x[1]) == sharedVert.end())
            {
                isshared = 0;
            }
            else
            {
                isshared = -1;
            }

            // If edge is not already in the edge list the add it to the list.
            if(edges.find(i2) == edges.end())
            {
                edges.insert(pair<Index2,int>(i2,isshared));

                if(isshared == 0)
                {
                    totnonshared++;
                }
                else
                {
                    totshared++;
                }
            }

            // 6
            i2.x[0] = volelems[i].verts[2];
            i2.x[1] = volelems[i].verts[3];

            i2.Sort();

            if(sharedVert.find(i2.x[0]) == sharedVert.end() || sharedVert.find(i2.x[1]) == sharedVert.end())
            {
                isshared = 0;
            }
            else
            {
                isshared = -1;
            }

            // If edge is not already in the edge list the add it to the list.
            if(edges.find(i2) == edges.end())
            {
                edges.insert(pair<Index2,int>(i2,isshared));

                if(isshared == 0)
                {
                    totnonshared++;
                }
                else
                {
                    totshared++;
                }
            }
        }

        totall = totshared + totnonshared;

        /// We have found the edges. And we also know whether they are fully shared by another one.

        std::cout << "RANK: " << mypid << " SIZE OF EDGES : " << totall << std::endl;
        std::cout << "RANK: " << mypid << " SIZE OF SHARED: " << totshared << std::endl;
        // std::cout << "NUMBER OF VOLELEMS: " << numvolume << std::endl;
        // std::cout << "NUMBER OF VERTICES: " << numvertex << std::endl;
        // std::cout << "MAX ELEMENT ID: " << maxelid << std::endl;

        /*
        for(myit = edges.begin(); myit != edges.end(); myit++)
        {
            std::cout << (*myit).first.x[0] << "   " << (*myit).first.x[1] << "  " << (*myit).second * '*' << std::endl;
        }
        */

        MPI_Barrier(comm);
        MPI_Allreduce(&maxelid, &maxglobalid, 1, MPI_INT, MPI_MAX, comm);
        MPI_Barrier(comm);
        std::cout << "MAX GLOBAL ID: " << maxglobalid << std::endl;

        /// Now we know max global id of a vertex. We should calculate how many new vertices each processor will have
        map<int,string>::iterator sharit;
        stringstream oss;
        string tmpstr;

        set<int> procset;
        set<int> combprocset;
        set<int>::iterator procsetit;

        int tt = 0;
        int tt2 = 0;	// count;
        int tt3 = 0;	// temp storage.

        int minid = 10000000;   // Increase this if using more than 10000000 processors. Unlikely.
        int totowner = 0;

        std::cout << "RANK: " << mypid << std::endl;

        for(myit = edges.begin(); myit != edges.end(); myit++)
        {
            // (*myit).first.x[0]
            // (*myit).first.x[1]
            // (*myit).second

            if((*myit).second == -1) // it is shared
            {
                /// here compare the shared processor id's of the two vertices.
                /// And the smallest processor id is the owner and will assign the global id of the new vertex.

                //std::cout << (*myit).first.x[0] << "    " << sharedVert.find((*myit).first.x[0])->second << std::endl;
                //std::cout << (*myit).first.x[1] << "    " << sharedVert.find((*myit).first.x[1])->second << std::endl;

                sharit = sharedVert.find((*myit).first.x[0]);

				if (sharit != sharedVert.end()) // This should always be true.
				{
					tmpstr = (*sharit).second;
					tt = (*sharit).first;

					oss << tmpstr;
					oss >> tt2;

					for (int i=0; i<tt2; i++) {
						oss >> tt3;                     // Processor ID (tt3-1) elmer rank and MPI rank are different.
						procset.insert(tt3-1);  // Map or set usage.
					}

					oss.clear();
					oss.str("");
				}
				else
				{
				    std::cout << "PROBLEM" << std::endl;
				}

				/// Now we have the processor id's from the first vertex.

				sharit = sharedVert.find((*myit).first.x[1]);

				if (sharit != sharedVert.end()) // This should always be true.
				{
					tmpstr = (*sharit).second;
					tt = (*sharit).first;

					oss << tmpstr;
					oss >> tt2;

					for (int i=0; i<tt2; i++) {
						oss >> tt3;                     // Processor ID (tt3-1) elmer rank and MPI rank are different.
						if(procset.find(tt3-1) != procset.end())    // If both contain the same element.
						{
                            combprocset.insert(tt3-1);  // Map or set usage.
						}
					}

					oss.clear();
					oss.str("");
				}
				else
				{
				    std::cout << "PROBLEM" << std::endl;
				}

				/// Now we have the processor id's in combprocset which contain both the vertices.

                minid = 10000000;

				for(procsetit = combprocset.begin(); procsetit != combprocset.end(); procsetit++)
				{
                    if(*procsetit < minid)
                    {
                        minid = *procsetit;
                    }
				}

				///Uncomment this part to assign somewhat uniform owner processors.

                /*
				int detlef = ((*myit).first.x[0] + (*myit).first.x[1]) % combprocset.size();
				for(procsetit = combprocset.begin(); procsetit != combprocset.end(); procsetit++)
				{
				    if(detlef == 0)
				    {
                        minid = *procsetit;
				    }
				    detlef--;
				}
				*/




                if((minid == mypid) && (combprocset.size() > 1) )
                {
                    (*myit).second = -2;         // If equals -2 means it is shared and the owner is this one.
                    totowner++;
                }
                else if( (minid == mypid) && (combprocset.size() == 1))   // If the size of the set is 1. This is the owner
                {
                    (*myit).second = 0;         // Internal vertex
                    totowner++;
                }
                else if(combprocset.size() < 1)
                {
                    std::cout << "PROBLEM - Combprocset Size Error!" << std::endl; // Should never happen.
                }
                else if((minid != mypid) && (combprocset.size() == 1))
                {
                    std::cout << "PROBLEM - Combprocset Size Error Type 2!" << std::endl; // Should never happen.
                }

                /*
                if(minid == mypid)
                {
                    (*myit).second = -2;         // If equals 2 means it is shared and the owner is this one.
                    totowner++;
                }
                else if( combprocset.size() == 1)   // If the size of the set is 1. This is the owner
                {
                    (*myit).second = 0;         // Internal vertex
                    totowner++;
                }*/

                procset.clear();
                combprocset.clear();
            }
        }

        /// Here -2 means shared and owned, -1 means shared but not owned, 0 means internal vertex.


        int myvertices = totowner + totnonshared;
        MPI_Barrier(comm);

        /// Now we have calculated the number of shared vertices which

        std::cout << "RANK: " << mypid << " OWNER VERTEX NUMBER: " << myvertices << std::endl;

        int data1 = myvertices;
        int npindex = 0;

        MPI_Scan(&data1, &npindex, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Barrier(comm);
        int startindex = npindex - myvertices + maxglobalid + 1;
        std::cout << "RANK: " << mypid << " NEW START INDEX: " << startindex << std::endl;

        /// Now we have the startindex for the vertices globalids. We should start creating the id's and then send
        /// The global id's of our owned vertices to the other ones.

        int curindex = startindex;

        int sendcount = 0;
        int recvcount = 0;
        int *sendprids;
        int *sendcounts;
        int *recvprids;
        int *recvcounts;

        int* sendpr5;
        int* recvpr5;

        int * send2 = new int[numprocs];
        for(int i=0; i<numprocs; i++)
        {
            send2[i] = 0;
        }


        for(myit = edges.begin(); myit != edges.end(); myit++)
        {
            // (*myit).first.x[0]
            // (*myit).first.x[1]
            // (*myit).second

            if((*myit).second == -1) // it is shared and i am not the owner
            {
                // Do nothing.
            }
            else if((*myit).second == -2) // it is shared and i am the owner.
            {

                /// Now here we should find out to which processor to send and what.

                sharit = sharedVert.find((*myit).first.x[0]);

				if (sharit != sharedVert.end()) // This should always be true.
				{
					tmpstr = (*sharit).second;
					tt = (*sharit).first;

					oss << tmpstr;
					oss >> tt2;

					for (int i=0; i<tt2; i++) {
						oss >> tt3;                     // Processor ID (tt3-1) elmer rank and MPI rank are different.
						procset.insert(tt3-1);  // Map or set usage.
					}

					oss.clear();
					oss.str("");
				}
				else
				{
				    std::cout << "PROBLEM2" << std::endl;
				}

				/// Now we have the processor id's from the first vertex.

				sharit = sharedVert.find((*myit).first.x[1]);

				if (sharit != sharedVert.end()) // This should always be true.
				{
					tmpstr = (*sharit).second;
					tt = (*sharit).first;

					oss << tmpstr;
					oss >> tt2;

					for (int i=0; i<tt2; i++) {
						oss >> tt3;                     // Processor ID (tt3-1) elmer rank and MPI rank are different.
						if(procset.find(tt3-1) != procset.end())    // If both contain the same element.
						{
                            combprocset.insert(tt3-1);  // Map or set usage.
						}
					}

					oss.clear();
					oss.str("");
				}
				else
				{
				    std::cout << "PROBLEM2" << std::endl;
				}

				/// Now we have the processor id's in combprocset which contain both the vertices.

				for(procsetit = combprocset.begin(); procsetit != combprocset.end(); procsetit++)
				{
                    if(*procsetit != mypid)
                    {
                        // Add to the send list.
                        send2[*procsetit]++;
                    }
				}

                procset.clear();
                combprocset.clear();

            }
            else if((*myit).second == 0)    // It is nonshared.
            {
                (*myit).second =  curindex;     // give a global id.
                curindex++;
            }
        }

        /// Find send counts START
        sendcount = 0;
        for(int i=0; i<numprocs; i++)
        {
            if(send2[i] == 0)
            {
                continue;
            }
            else
            {
                sendcount++;
            }
        }

        sendprids = new int[sendcount];
        sendcounts = new int[sendcount];

        sendcount = 0;
        for(int i=0; i<numprocs; i++)
        {
            if(send2[i] == 0)
            {
                continue;
            }
            else
            {
                sendcounts[sendcount] = send2[i];
                send2[i] = 1;
                sendprids[sendcount] = i;
                sendcount++;
            }
        }

        /// Find receive counts
        {
            int *totdiffrecvs = new int[numprocs];  // total number of different processors to receive from.
            MPI_Allreduce(send2, totdiffrecvs, numprocs,MPI_INT, MPI_SUM, comm);

            recvcount = totdiffrecvs[mypid];
            delete [] totdiffrecvs;

            recvcounts = new int[recvcount];
            recvprids = new int[recvcount];
        }

        MPI_Request * reqs = new MPI_Request[recvcount + sendcount];
        MPI_Status * stats = new MPI_Status[recvcount + sendcount];

        /// Now send to each one how much you will send.
        for(int i=0; i<sendcount; i++)
        {
            // Send number of volume elements to be send.
            MPI_Isend(&sendcounts[i], 1, MPI_INT, sendprids[i], 123, comm, &reqs[i]);
        }
        for(int i=0; i<recvcount; i++)
        {
            MPI_Irecv(&recvcounts[i], 1, MPI_INT, MPI_ANY_SOURCE, 123, comm, &reqs[sendcount+i]);
        }
        MPI_Waitall(recvcount + sendcount, reqs, stats);

        for(int i=0; i<recvcount; i++)
        {
            recvprids[i] = stats[sendcount + i].MPI_SOURCE;
        }

        delete [] reqs;
        delete [] stats;
        /// Now each processor knows how much it will receive from which processor. (recvcounts & recvprids)
        /// And how much it will send to which processor. (sendcounts & sendprids)

        /// Find send and recv counts END

        recvpr5 = new int[recvcount+1];
        sendpr5 = new int[sendcount+1];
        //
        sendpr5[0] = 0;
        for (int j=0; j<sendcount; j++) {
            sendpr5[j+1] = sendpr5[j] + sendcounts[j];
        }

        recvpr5[0] = 0;
        for (int i=0; i<recvcount; i++) {
            recvpr5[i+1] = recvpr5[i] + recvcounts[i];
        }

        int *currentel5 = new int[sendcount];
        for (int i=0; i<sendcount; i++) {
            currentel5[i] = 0;
        }

        /// Now we know which processor will send what, which will get what and how many from whom.
        /// We should fill the contents.

        std::cout << "RANK: " << mypid << " has to receive: " << recvpr5[recvcount] << "   " << totshared - totowner << std::endl;
        std::cout << "RANK: " << mypid << " has to send: " << sendpr5[sendcount] << "   " << totowner << std::endl;

        int* datatosend = new int[3 * sendpr5[sendcount]];
        int curin = 0;
        int curin2 = 0;

        for(myit = edges.begin(); myit != edges.end(); myit++)
        {
            if((*myit).second == -2) // it is shared and i am the owner.
            {

                /// Now here we should find out to which processor to send and what.

                sharit = sharedVert.find((*myit).first.x[0]);

				if (sharit != sharedVert.end()) // This should always be true.
				{
					tmpstr = (*sharit).second;
					tt = (*sharit).first;

					oss << tmpstr;
					oss >> tt2;

					for (int i=0; i<tt2; i++) {
						oss >> tt3;                     // Processor ID (tt3-1) elmer rank and MPI rank are different.
						procset.insert(tt3-1);  // Map or set usage.
					}

					oss.clear();
					oss.str("");
				}
				else
				{
				    std::cout << "PROBLEM2 during actual data fill to send" << std::endl;
				}

				/// Now we have the processor id's from the first vertex.

				sharit = sharedVert.find((*myit).first.x[1]);

				if (sharit != sharedVert.end()) // This should always be true.
				{
					tmpstr = (*sharit).second;
					tt = (*sharit).first;

					oss << tmpstr;
					oss >> tt2;

					for (int i=0; i<tt2; i++) {
						oss >> tt3;                     // Processor ID (tt3-1) elmer rank and MPI rank are different.
						if(procset.find(tt3-1) != procset.end())    // If both contain the same element.
						{
                            combprocset.insert(tt3-1);  // Map or set usage.
						}
					}

					oss.clear();
					oss.str("");
				}
				else
				{
				    std::cout << "PROBLEM2 during actual data fill to send" << std::endl;
				}

				/// Now we have the processor id's in combprocset which contain both the vertices.

				(*myit).second =  curindex;     // give a global id.
                curindex++;

				for(procsetit = combprocset.begin(); procsetit != combprocset.end(); procsetit++)
				{
                    if(*procsetit != mypid)
                    {
                        /// Now fill the content to send.
                        curin = findEl(sendprids, *procsetit, sendcount);
                        curin2 = 3*(sendpr5[curin] + currentel5[curin]);

                        datatosend[curin2] = (*myit).first.x[0];
                        datatosend[curin2+1] = (*myit).first.x[1];
                        datatosend[curin2+2] = (*myit).second;

                        if((*myit).second == -2)
                        {
                            std::cout << "PROBLEM2 - ID to be send should'nt equal -2" << std::endl;
                        }

                        currentel5[curin]++;
                    }
				}

                procset.clear();
                combprocset.clear();
            }
        }

        int maxcurrentlid = curindex;

        int* datattorecv = new int[3 * recvpr5[recvcount]];

        MPI_Barrier(comm);

        // Send Elements
        int tmpp2;

        reqs = new MPI_Request[sendcount + recvcount];
        stats = new MPI_Status[sendcount + recvcount];

        for (int i=0; i<sendcount; i++) {
            tmpp2 = 3*(sendpr5[i+1] - sendpr5[i]);
            MPI_Isend(&datatosend[3*sendpr5[i]], tmpp2, MPI_INT, sendprids[i], 123, comm, &reqs[i]);
        }

        for (int i=0; i<recvcount; i++) {
            tmpp2 = 3* (recvpr5[i+1] - recvpr5[i]);
            MPI_Irecv(&datattorecv[3*recvpr5[i]], tmpp2, MPI_INT, recvprids[i], 123, comm, &reqs[sendcount + i]);
        }
        MPI_Waitall(recvcount + sendcount, reqs, stats);

        delete [] reqs;
        delete [] stats;

        // std::cout << "RANK: " << mypid << " is here." << std::endl;
        MPI_Barrier(comm);

        /// Now fill in the new global id's for nonowned shared vertices.

        std::cout << "RECC COUNT: " << recvpr5[recvcount] << "  and  the expected recv count is: "  << totshared - totowner << std::endl;

        int totproblemm = 0;

        std::list<Index3> excessList;
        Index3 i3;
        std::list<Index3>::iterator excessit;

        for(int i=0; i<recvpr5[recvcount]; i++)
        {
            i2.x[0] = datattorecv[3*i];
            i2.x[1] = datattorecv[3*i+1];

            i2.Sort();

            myit = edges.find(i2);

            if(myit == edges.end())
            {
                /// This creates excess problem. Too many data is sent.
                // std::cout << "EKSTRA: " << mypid << "   " << i2.x[0] << "   "  << i2.x[1] << "   " << datattorecv[3*i+2] << std::endl;

                i3.x[0] = i2.x[0];
                i3.x[1] = i2.x[1];
                i3.x[2] = datattorecv[3*i+2];

                excessList.push_back(i3);

                totproblemm++;

                continue;
            }
            // std::cout << datattorecv[3*i+2] << std::endl;

            if (myit->second != -1)
            {
                std::cout << "PROBLEM - Processor ownership conflict!" << std::endl;
            }
            else
            {
                myit->second = datattorecv[3*i+2];
            }
        }

        int totabsent = 0;

        for(myit = edges.begin(); myit != edges.end(); myit++)
        {
            /// This is the absess problem. Not received data that should be here.
            if((*myit).second == -1)
            {
                totabsent++;
            }
        }

        int *ids1 = new int[totabsent];
        int *ids2 = new int[totabsent];

        int qq = 0;

        for(myit = edges.begin(); myit != edges.end(); myit++)
        {
            if((*myit).second == -1)
            {
                // std::cout << "ABSENT: " << mypid << "   " << myit->first.x[0] << "   "  << myit->first.x[1] << "   " << -1 << std::endl;
                // totproblemm--;
                ids1[qq] = myit->first.x[0];
                ids2[qq] = myit->first.x[1];
                qq++;
            }
        }

        /// For processing of received information from this function.
        /// (Beware; if shared proc list has only one procesor then it is an internal vertex.)
        /// In this case no one gets extra information.
        /// Do not forget (procid = procid + 1) for elmer
        solveAbsentProblem(ids1, ids2, totabsent, maxcurrentlid);
        // std::cout << "RANK: " << mypid << " received the following absent Data: " << std::endl;
        // std::cout << absentData;

        /// Adding the absent global id's.
        std::stringstream oss5(stringstream::in | stringstream::out);
        oss5 << absentData;

        int a1,b1,gl2,c2, tp2;
        while(oss5.good() && oss5.str().length() != 0)
        {
            oss5 >> a1 >> b1 >> gl2 >> c2;

            for(int i=0;i<c2;i++)
            {
                oss5 >> tp2;
            }

            i2.x[0] = a1;
            i2.x[1] = b1;
            i2.Sort();
            edges.find(i2)->second = gl2;

            // std::cout << a1 << " " << b1 << " " << gl2 << std::endl;
        }

        int excesssize = excessList.size();
        int *ids3 = new int[excesssize];
        int *ids4 = new int[excesssize];
        int *idsgl = new int[excesssize];

        qq = 0;

        for(excessit = excessList.begin(); excessit != excessList.end(); excessit++)
        {
            i3 = *excessit;
            ids3[qq] = i3.x[0];
            ids4[qq] = i3.x[1];
            idsgl[qq] = i3.x[2];
            qq++;
        }

        solveExcessProblem(ids3, ids4, idsgl, excesssize);

        // std::cout << "RANK: " << mypid << " received the following excess Data: " << std::endl;
        // std::cout << excessData << excessData2;


        /// After both processes have completed.
        /// First update global ids of those who were -1. (uninitialized)
        /// The following parts can then proceed.
        /// The shared vertex information should be updated with the received data afterwards.

        delete [] datatosend;
        delete [] datattorecv;
        /// Now we have sent everything needed. Lets calculate what's left.

        /// Create the new vertices;
        verts.resize(numvertex + edges.size());
        int cs = 0;
        int a3 = 0;
        int b3 = 0;

        /// Project Edges on Boundary to Geometry

        int boundedgesize = edgesonbound.size();
        double* xdeger = new double[boundedgesize];
        double* ydeger = new double[boundedgesize];
        double* zdeger = new double[boundedgesize];
        int* yuzeydeger = new int[boundedgesize];

        std::cout << "EDGES ON BOUNDARY SIZE OF RANK: " << mypid << " is: " << boundedgesize << " versus " << edges.size() << std::endl;

        /// Now project them on the surface. Get the results on xdeger, ydeger, zdeger. (This will be an NG Call).
        using namespace nglib;

        int pp1 = 0;
        int pa1 = 0;
        int pa2 = 0;

        for(myit = edgesonbound.begin(); myit != edgesonbound.end(); myit++)
        {
            pa1 = myit->first.x[0];
            pa2 = myit->first.x[1];

            pa1 = vertmap.find(pa1)->second;
            pa2 = vertmap.find(pa2)->second;

            xdeger[pp1] = (verts[pa1].coor[0] + verts[pa2].coor[0])/2;
            ydeger[pp1] = (verts[pa1].coor[1] + verts[pa2].coor[1])/2;
            zdeger[pp1] = (verts[pa1].coor[2] + verts[pa2].coor[2])/2;

            yuzeydeger[pp1] = (*myit).second;
            pp1++;
        }

        Ng_CSG_ProjectMesh(inputFile, boundedgesize, xdeger, ydeger, zdeger, yuzeydeger);

        MPI_Barrier(comm);

        pp1 = 0;
        for(myit = edgesonbound.begin(); myit != edgesonbound.end(); myit++)
        {
            (*myit).second = pp1;
            pp1++;
        }

        int totalscheus = 0;
        bool isonface = false;
        int indexxx = 0;

        /// Project Edges on Boundary to Geometry END

        for(myit = edges.begin(); myit != edges.end(); myit++)
        {
            a3 = (*myit).first.x[0];
            b3 = (*myit).first.x[1];

            /// Check if edge is on boundary

            i2.x[0] = a3;
            i2.x[1] = b3;

            i2.Sort();

            if(edgesonbound.find(i2) != edgesonbound.end())
            {
                isonface = true;
                indexxx = edgesonbound.find(i2)->second;
            }
            /// Check if edge is on boundary END

            /// Create shared information
            if(sharedVert.find(a3) != sharedVert.end() || sharedVert.find(b3) != sharedVert.end())
            {
                sharit = sharedVert.find(a3);

				if (sharit != sharedVert.end()) // This should always be true.
				{
					tmpstr = (*sharit).second;
					tt = (*sharit).first;

					oss << tmpstr;
					oss >> tt2;

					for (int i=0; i<tt2; i++) {
						oss >> tt3;                     // Processor ID (tt3-1) elmer rank and MPI rank are different.
						procset.insert(tt3-1);  // Map or set usage.
					}

					oss.clear();
					oss.str("");
				}

				/// Now we have the processor id's from the first vertex.

				sharit = sharedVert.find(b3);

				if (sharit != sharedVert.end()) // This should always be true.
				{
					tmpstr = (*sharit).second;
					tt = (*sharit).first;

					oss << tmpstr;
					oss >> tt2;

					for (int i=0; i<tt2; i++) {
						oss >> tt3;                     // Processor ID (tt3-1) elmer rank and MPI rank are different.
						if(procset.find(tt3-1) != procset.end())    // If both contain the same element.
						{
                            combprocset.insert(tt3-1);  // Map or set usage.
						}
					}

					oss.clear();
					oss.str("");
				}

				/// Now we have the processor id's in combprocset which contain both the vertices.
                if(combprocset.size() > 1)
                {
                    // myit->second;
                    oss << combprocset.size();

                    for(procsetit = combprocset.begin(); procsetit != combprocset.end(); procsetit++)
                    {
                        oss << " " << (*procsetit) + 1;      // Elmer procid = MPIprocid + 1
                    }



                    if((*myit).second == -1)
                    {
                        std::cout << "PROBLEMMPPPP!!!!  Missing data: " << a3 <<  "   "  << b3 << std::endl;
                    }
                    else if((*myit).second <= 0)
                    {
                        std::cout << "ERROR NEVER SHOULD HAVE HAPPENED!" << std::endl;
                    }
                    else
                    {
                        sharedVert[(*myit).second] = oss.str();
                        numshared++;    // This is the private member of the MeshMig class. Not a local variable.
                    }

                    oss.clear();
                    oss.str("");
                }

                procset.clear();
                combprocset.clear();
            }

            /// Create the new vertices

            a3 = vertmap.find(a3)->second;
            b3 = vertmap.find(b3)->second;

            // std::cout << a3 << "  " << b3 << std::endl;

            if(myit->second == 0 || myit->second == -1 || myit->second == -2)
            {
                totalscheus++;
            }

            verts[numvertex + cs].idx = myit->second;
            verts[numvertex + cs].useCount = 0;

            /// Put projected values if edge is on boundary.
            if(isonface != true)
            {
                verts[numvertex + cs].coor[0] = (verts[a3].coor[0] + verts[b3].coor[0])/2;
                verts[numvertex + cs].coor[1] = (verts[a3].coor[1] + verts[b3].coor[1])/2;
                verts[numvertex + cs].coor[2] = (verts[a3].coor[2] + verts[b3].coor[2])/2;
            }
            else
            {
                verts[numvertex + cs].coor[0] = xdeger[indexxx];
                verts[numvertex + cs].coor[1] = ydeger[indexxx];
                verts[numvertex + cs].coor[2] = zdeger[indexxx];
            }

            vertmap[myit->second] = numvertex + cs;

            cs++;
            isonface = false;
        }
        numvertex += edges.size();

        std::cout << "RANK: " << mypid << " has Total Sh.. " << totalscheus << std::endl;

        /// Updating absent Data related shared processor info errors.
        std::stringstream oss6(stringstream::in | stringstream::out);
        std::stringstream oss7(stringstream::in | stringstream::out);

        oss6.clear();
        oss6.str("");
        oss6 << absentData;

        oss7.clear();
        oss7.str("");

        while(!oss6.eof())
        {
            oss6 >> a1 >> b1 >> gl2 >> c2;

            if(oss6.eof())
            {
                break;
            }

            oss7 << c2;

            for(int i=0;i<c2;i++)
            {
                oss6 >> tp2;
                oss7 << " " << tp2 + 1;
            }

            tmpstr = oss7.str();

            if(c2 == 1)
            {
                // delete shared vertex from list.
                sharedVert.erase(gl2);
                // std::cout << "RANK: " << mypid << " has deleted shared info with global id: " << gl2 << std::endl;
            }
            else if(c2 > 1)
            {
                sharedVert[gl2] = tmpstr;
                // std::cout << "RANK: " << mypid << " has updated shared info with global id: " << gl2 << " to: " << tmpstr << std::endl;
            }

            oss7.clear();
            oss7.str("");
        }

        /// Update excessData
        oss6.clear();
        oss6.str("");
        oss6 << excessData;

        oss7.clear();
        oss7.str("");

        while(!oss6.eof())
        {
            oss6 >> a1 >> b1  >> c2;

            if(oss6.eof())
            {
                break;
            }

            i2.x[0] = a1;
            i2.x[1] = b1;
            i2.Sort();
            gl2 = edges.find(i2)->second;

            if(c2 == 1)
            {
                // delete shared vertex from list.
                sharedVert.erase(gl2);
            }

            oss7 << c2;

            for(int i=0;i<c2;i++)
            {
                oss6 >> tp2;
                oss7 << " " << tp2 + 1;
            }

            if(c2 > 1)
            {
                sharedVert[gl2] = oss7.str();
            }

            oss7.clear();
            oss7.str("");
        }

        /// Update excessData 2
        oss6.clear();
        oss6.str("");
        oss6 << excessData2;

        oss7.clear();
        oss7.str("");

        while(!oss6.eof())
        {
            oss6 >> a1 >> b1  >> c2;

            if(oss6.eof())
            {
                break;
            }

            i2.x[0] = a1;
            i2.x[1] = b1;
            i2.Sort();
            gl2 = edges.find(i2)->second;

            if(c2 == 1)
            {
                // delete shared vertex from list.
                sharedVert.erase(gl2);
            }

            oss7 << c2;

            for(int i=0;i<c2;i++)
            {
                oss6 >> tp2;
                oss7 << " " << tp2 + 1;
            }

            if(c2 > 1)
            {
                sharedVert[gl2] = oss7.str();
            }

            oss7.clear();
            oss7.str("");
        }

        /// Create the new boundaryelements.

        boundelems.resize(numboundary*4);

        int cind = numboundary*4-4;

        int a4,b4,c4,d4;
        int ab,ac,ad,bc,bd,cd;

        for(int i=numboundary-1; i>=0; i--)
        {
            a4 = boundelems[i].verts[0];
            b4 = boundelems[i].verts[1];
            c4 = boundelems[i].verts[2];

            i2.x[0] = a4;
            i2.x[1] = b4;
            i2.Sort();
            ab = edges.find(i2)->second;

            i2.x[0] = a4;
            i2.x[1] = c4;
            i2.Sort();
            ac = edges.find(i2)->second;

            i2.x[0] = b4;
            i2.x[1] = c4;
            i2.Sort();
            bc = edges.find(i2)->second;

            boundelems[cind].verts[0] = a4;
            boundelems[cind].verts[1] = ab;
            boundelems[cind].verts[2] = ac;

            boundelems[cind+1].verts[0] = ab;
            boundelems[cind+1].verts[1] = b4;
            boundelems[cind+1].verts[2] = bc;

            boundelems[cind+2].verts[0] = ac;
            boundelems[cind+2].verts[1] = bc;
            boundelems[cind+2].verts[2] = c4;

            boundelems[cind+3].verts[0] = ab;
            boundelems[cind+3].verts[1] = bc;
            boundelems[cind+3].verts[2] = ac;

            boundelems[cind].geoid = boundelems[i].geoid;
            boundelems[cind+1].geoid = boundelems[i].geoid;
            boundelems[cind+2].geoid = boundelems[i].geoid;
            boundelems[cind+3].geoid = boundelems[i].geoid;

            boundelems[cind].parentelement = 0;     // This will be set later.
            boundelems[cind+1].parentelement = 0;
            boundelems[cind+2].parentelement = 0;
            boundelems[cind+3].parentelement = 0;

            boundelems[cind].idx = cind+1;          // ID's are local.
            boundelems[cind+1].idx = cind+2;
            boundelems[cind+2].idx = cind+3;
            boundelems[cind+3].idx = cind+4;

            cind -= 4;
        }

        numboundary *= 4;


        /// Create the new volume elements.
        /// Don't forget updating the useCounts of vertices.
        volelems.resize(numvolume*8);
        cind = numvolume*8-8;

        int newvolsize = numvolume * 8;
        int newvolindex = 0;

        MPI_Scan(&newvolsize, &newvolindex, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Barrier(comm);
        newvolindex = newvolindex - newvolsize;

        std::cout << "RANK: " << mypid << "  has new volsize: " << newvolsize << std::endl;
        std::cout << "RANK: " << mypid << "  has new index: " << newvolindex << std::endl;

        MPI_Barrier(comm);

        volmap.clear();

        for(int i=numvolume-1; i>=0; i--)
        {
            a4 = volelems[i].verts[0];
            b4 = volelems[i].verts[1];
            c4 = volelems[i].verts[2];
            d4 = volelems[i].verts[3];

            i2.x[0] = a4;
            i2.x[1] = b4;
            i2.Sort();
            ab = edges.find(i2)->second;

            i2.x[0] = a4;
            i2.x[1] = c4;
            i2.Sort();
            ac = edges.find(i2)->second;

            i2.x[0] = a4;
            i2.x[1] = d4;
            i2.Sort();
            ad = edges.find(i2)->second;

            i2.x[0] = b4;
            i2.x[1] = c4;
            i2.Sort();
            bc = edges.find(i2)->second;

            i2.x[0] = b4;
            i2.x[1] = d4;
            i2.Sort();
            bd = edges.find(i2)->second;

            i2.x[0] = c4;
            i2.x[1] = d4;
            i2.Sort();
            cd = edges.find(i2)->second;

            /// Create the new elements.

            volelems[cind].verts[0] = a4;
            volelems[cind].verts[1] = ab;
            volelems[cind].verts[2] = ac;
            volelems[cind].verts[3] = ad;

            volelems[cind+1].verts[0] = ab;
            volelems[cind+1].verts[1] = b4;
            volelems[cind+1].verts[2] = bc;
            volelems[cind+1].verts[3] = bd;

            volelems[cind+2].verts[0] = ab;
            volelems[cind+2].verts[1] = bc;
            volelems[cind+2].verts[2] = ac;
            volelems[cind+2].verts[3] = bd;

            volelems[cind+3].verts[0] = ac;
            volelems[cind+3].verts[1] = bc;
            volelems[cind+3].verts[2] = c4;
            volelems[cind+3].verts[3] = cd;

            volelems[cind+4].verts[0] = ac;
            volelems[cind+4].verts[1] = bc;
            volelems[cind+4].verts[2] = cd;
            volelems[cind+4].verts[3] = bd;

            volelems[cind+5].verts[0] = ad;
            volelems[cind+5].verts[1] = cd;
            volelems[cind+5].verts[2] = d4;
            volelems[cind+5].verts[3] = bd;

            volelems[cind+6].verts[0] = ac;
            volelems[cind+6].verts[1] = cd;
            volelems[cind+6].verts[2] = ad;
            volelems[cind+6].verts[3] = bd;

            volelems[cind+7].verts[0] = ab;
            volelems[cind+7].verts[1] = ac;
            volelems[cind+7].verts[2] = ad;
            volelems[cind+7].verts[3] = bd;


            /// Global ID's

            volelems[cind].idx = newvolindex + cind + 1;          // ID's are global.
            volelems[cind + 1].idx = newvolindex + cind + 2;
            volelems[cind + 2].idx = newvolindex + cind + 3;
            volelems[cind + 3].idx = newvolindex + cind + 4;
            volelems[cind + 4].idx = newvolindex + cind + 5;
            volelems[cind + 5].idx = newvolindex + cind + 6;
            volelems[cind + 6].idx = newvolindex + cind + 7;
            volelems[cind + 7].idx = newvolindex + cind + 8;

            volmap[newvolindex + cind + 1] = cind;                    // ID's are global.
            volmap[newvolindex + cind + 2] = cind + 1;
            volmap[newvolindex + cind + 3] = cind + 2;
            volmap[newvolindex + cind + 4] = cind + 3;
            volmap[newvolindex + cind + 5] = cind + 4;
            volmap[newvolindex + cind + 6] = cind + 5;
            volmap[newvolindex + cind + 7] = cind + 6;
            volmap[newvolindex + cind + 8] = cind + 7;

            cind -= 8;
        }
        numvolume *= 8;

        /// Update boundary elements parent ids'
        setBoundaryElements();

        /// Volbound'u doldur silbastan.
        volbound.clear();
        for(int i=0; i<numboundary; i++)
        {
            volbound[boundelems[i].parentelement] = boundelems[i].idx;
        }


        // std::cout << "TOTAL PROBLEM FACTOR OF RANK " << mypid << " IS " << totproblemm << std::endl;

    }

    void createMetisMesh(MPI_Comm comm, const char* inputFile, int minSize) // minSize is the total size of the first generated volume element.
    {
        if(mypid == 0)
        {
            /// Generate the mesh. Partition it using metis. And create the distribution data to be sent.
            /// GENERATION PART
            using namespace nglib;
            Ng_Init();

            long nvertices, ntriangles, nelems;

            double *vertexlist;
            int *trianglelist;
            int *triangleFN;        // Boundary Element Index.
            int *triangleDIN;
            int *triangleDOUT;
            int *triangleParent;
            int *elemlist;

            //vertexlist = new double[12000000];
            //trianglelist = new int[12000000];
            //triangleFN = new int[4000000];

            Ng_CSG_GenerateVolumeMesh(inputFile, minSize, &nvertices, &ntriangles, &nelems);
            // Get first the number to be allocated.

            vertexlist = new double[nvertices*3];
            trianglelist = new int[ntriangles*3];
            triangleFN = new int[ntriangles];
            triangleDIN = new int[ntriangles];
            triangleDOUT = new int[ntriangles];
            triangleParent = new int[ntriangles];
            elemlist = new int[nelems*4];

            Ng_CSG_GenerateVolumeMesh(inputFile, minSize, vertexlist, trianglelist, triangleFN, triangleDIN, triangleDOUT, elemlist);
            std::cout << "Successfully loaded .geo File seond time: " << inputFile << std::endl;

            /*
            for(int i=0; i< nelems; i++)
            {
                std::cout << elemlist[4*i] << " " << elemlist[4*i+1] << " " << elemlist[4*i+2] << " " << elemlist[4*i+3] << " " << std::endl;
            }
            */

            std::cout << "VERTS: " << nvertices << std::endl;
            std::cout << "TRIS: " << ntriangles << std::endl;
            std::cout << "ELEMS: " << nelems << std::endl;



            MPI_Bcast(&nelems, 1, MPI_INT, 0, comm);
            MPI_Barrier(comm);

            /// GENERATION PART END

            /// METIS SETUP

            int *eptr2 = new int[nelems-numprocs+1+1]; // 0 2 5 8 11 13 (ne+1 adet) Degisken olabiliyor.Ama bizim durumumuzda sadece tetrahedral var. Yani 0 4 8 12 16 ...
            int * elemts = new int[4*(nelems-numprocs+1)];        // [12 11 7 3] [11 12 8 9] Sirayla elemenlarin indexleri yaziliyor.
            int *part = new int[nelems];
            float *ubvec = new float[1];     // balance constraint
            float *tpwgts = new float[numprocs]; // normally ncon*nparts.
            int * elmdist = new int[numprocs+1];    // analogue to vtxdist

            // int *part2 = new int[nvertices];

            elmdist[0] = 0;             // Since all elements are at the first processor for now. Send the last numprocs-1 elements one by one to each other.
            elmdist[1] = nelems-numprocs+1;
            for(int i=2; i<=numprocs; i++)
            {
                elmdist[i] = elmdist[i-1] + 1;
            }

            int numflag = 0;    // indexing starts from 0 (set to 1 if index start from 1)
            int wgtflag = 0;    // not weighted
            int ncon = 1;
            int ncommonnodes = 3;   // for tetrahedral.


            int nparts = numprocs;      // size of partition.

            ubvec[0] = 1.05;        // recommended value by parmetis
            int options[10];        // options.
            int edgecut;            // edges cut in total by the new partitioning.

            options[0] = 1;
            options[PMV3_OPTION_DBGLVL] = 7;
            options[PMV3_OPTION_SEED] = 0;

            eptr2[0] = 0;
            for(int i=0; i<(nelems-numprocs+1); i++)
            {
                eptr2[i+1] = eptr2[i] + 4;
            }

            for(int i=0; i<(nelems-numprocs+1); i++)
            {
                elemts[4*i] = elemlist[4*i];
                elemts[4*i+1] = elemlist[4*i+1];
                elemts[4*i+2] = elemlist[4*i+2];
                elemts[4*i+3] = elemlist[4*i+3];
            }

            for(int kk=0; kk<numprocs; kk++)
            {
                tpwgts[kk] = 1.0/(float)(nparts);
            }


            /// METIS SETUP END

            std::cout << "I AM HERE before send" << std::endl;

            /// CALL METIS   (We will be using Parmetis with all elements first at node 0)
            MPI_Barrier(comm);
            // here send one element to each other processor.

            int startindex = nelems - numprocs + 1;

            for(int i=0; i< numprocs-1; i++)
            {
                std::cout << "SENDING TO RANK: " << i+1 << " " << nelems << " " << elemlist[4*(startindex+i)] << " " << elemlist[4*(startindex+i)+1] << " " << elemlist[4*(startindex+i)+2] << " " << elemlist[4*(startindex+i)+3] << " " << std::endl;
                MPI_Send(&elemlist[4*(startindex+i)], 4, MPI_INT, i+1, 123, comm);
            }

            MPI_Barrier(comm);
            ParMETIS_V3_PartMeshKway(elmdist, eptr2, elemts, NULL, &wgtflag, &numflag, &ncon, &ncommonnodes, &nparts, tpwgts, ubvec, options, &edgecut, part, &comm);
            // METIS_PartMeshNodal(nelems, nvertices, eptr2, elemts, NULL, NULL, &nparts, NULL, NULL, &edgecut, part, part2);

            // int nell = nelems;

            // METIS_mCPartGraphKway(&nell, &ncon, eptr2, elemts, NULL, NULL, &nparts, NULL, NULL, NULL, NULL, &edgecut, part);
            /// CALL METIS END

            MPI_Barrier(comm);

            MPI_Status * stats = new MPI_Status[numprocs-1];
            MPI_Request * reqs = new MPI_Request[numprocs-1];

            startindex = nelems - numprocs + 1;

            for(int i=0; i<numprocs-1; i++)
            {
                MPI_Irecv(&part[startindex+i], 1, MPI_INT, i+1, 123, comm, &reqs[i]);
            }

            MPI_Waitall(numprocs-1, reqs, stats);

            MPI_Barrier(comm);

            for(int i=0; i<numprocs-1; i++)
            {
                std::cout << "RECEIVED FROM: " << i+1 << " " << part[startindex+1] << std::endl;
            }

            delete [] reqs;
            delete [] stats;

            /*
            for(int i=0; i<nelems; i++)
            {
                std::cout << part[i] << std::endl;
            }
            */


            std::cout << "RANK: " << mypid << " => I AM HERE3" << std::endl;

            // We can delete the unnecessary data.
            delete [] eptr2;
            delete [] tpwgts;
            delete [] elmdist;
            delete [] elemts;

            MPI_Barrier(comm);          // Barrier 4

            /// DISTRIBUTE DATA
            set<int>* sharedData;
            sharedData = new set<int>[nvertices];       // Set array of size = number of vertices

            // Fill the sharedData;
            int tmpvid = 0; //tmp vertex id.
            int tmppid = 0; //tmp processor id.

            // Important. Filling vertex id's with id-1. So don't forget to set it back.

            for(int i=0; i<nelems; i++)
            {
                tmppid = part[i];

                tmpvid = elemlist[4*i];
                sharedData[tmpvid-1].insert(tmppid);

                tmpvid = elemlist[4*i+1];
                sharedData[tmpvid-1].insert(tmppid);

                tmpvid = elemlist[4*i+2];
                sharedData[tmpvid-1].insert(tmppid);

                tmpvid = elemlist[4*i+3];
                sharedData[tmpvid-1].insert(tmppid);
            }
            // Shared data is filled.
            // Volume Element data is filled.
            // Vertex data is known.
            // So what's left is the surface to volume relation.

            /// CREATE PARENT DATA

            bool (*fn_pt)(Index3,Index3) = fncomp;
            std::multimap<Index3,int, bool(*)(Index3, Index3)> face2vol(fn_pt);
            std::multimap<Index3,int, bool(*)(Index3, Index3)>::iterator myit;

            Index3 i3;
            int l;

            for(int i=0; i<nelems; i++)
            {
                for (int j = 1; j <= 4; j++)   // loop over faces of tet
                {
                    l = 0;
                    for (int k = 1; k <= 4; k++)
                    {
                        if (k != j)
                        {
                            i3.x[l] = elemlist[4*i+k-1];
                            l++;
                        }
                    }

                    i3.Sort();

                    face2vol.insert(pair<Index3,int>(i3,i+1));
                }
            }

            for(int i=0; i<ntriangles; i++)
            {
                i3.x[0] = trianglelist[3*i];
                i3.x[1] = trianglelist[3*i+1];
                i3.x[2] = trianglelist[3*i+2];

                i3.Sort();

                myit = face2vol.find(i3);

                if(myit == face2vol.end())
                {
                    std::cout << "////////////" << endl;
                    std::cout << "////////////" << endl;
                    std::cout << "It Happens: " << endl;
                    std::cout << "////////////" << endl;
                    std::cout << "////////////" << endl;
                }

                triangleParent[i] = (*myit).second;
            }

            face2vol.clear();

            // Now we have everything we need.
            // Let's send each element, vertex, face and shared data info to the appropriate processor.

            /// CREATE PARENT DATA END

            /// SEND DATA

            MPI_Barrier(comm);   // Barrier 5

            std::string* stringtosend;
            stringtosend = new std::string[numprocs];

            std::stringstream* sfile;
            sfile = new std::stringstream[numprocs];

            int * datanumber = new int[numprocs*4];     // i = vert i+1 = face i+2 = elem i+3 = stringsize

            for(int i=0; i<4*numprocs; i++)
            {
                datanumber[i] = 0;
            }

            int sizeofset;
            int curselected;
            std::set<int>::iterator myiter;

            std::stringstream tmpstream;
            int val = 0;                        // O means not shared. 1 means shared.

            for(int i=0; i<nvertices; i++)
            {
                sizeofset = sharedData[i].size();

                if(sizeofset == 1)
                {
                    val = 0;
                }
                else
                {
                    val = 1;
                }

                tmpstream << sizeofset;

                for(myiter = sharedData[i].begin(); myiter != sharedData[i].end(); myiter++)
                {
                    tmpstream << " " << *myiter;
                }

                for(myiter = sharedData[i].begin(); myiter != sharedData[i].end(); myiter++)
                {
                    curselected = *myiter;
                    // Enter vertex data.   Vertex ID begins from 1
                    sfile[curselected] << i+1 << " " << vertexlist[3*i] << " " << vertexlist[3*i+1] << " " << vertexlist[3*i+2] << " " << val << std::endl;

                    // If shared then also enter shared line.
                    if(val == 1)
                    {
                        sfile[curselected] << tmpstream.str() << std::endl;
                    }

                    datanumber[4*curselected]++;
                }

                tmpstream.clear();
                tmpstream.str("");
            }


            for(int i=0; i<ntriangles; i++)
            {
                curselected = triangleParent[i];        // get parent volume element.
                //std::cout << curselected << std::endl;
                curselected = part[curselected-1];      // find out which processor it is send.
                //std::cout << curselected << std::endl;
                // std::cout << "TEST" << std::endl;

                // Ids begin with 1.
                sfile[curselected] << i+1 << " " << triangleFN[i] << " " << triangleParent[i] << " " << trianglelist[3*i] << " " << trianglelist[3*i+1] << " " << trianglelist[3*i+2] << std::endl;

                datanumber[4*curselected+1]++;
            }

            for(int i=0; i<nelems; i++)
            {
                curselected = part[i];

                sfile[curselected] << i+1 << " " << elemlist[4*i] << " " << elemlist[4*i+1] << " " << elemlist[4*i+2] << " " << elemlist[4*i+3] << " " << std::endl;

                datanumber[4*curselected+2]++;
            }

            for(int i=0; i<numprocs; i++)
            {
                stringtosend[i] = sfile[i].str();
                datanumber[4*i+3] = stringtosend[i].length();
                sfile[i].clear();
                sfile[i].str("");
            }


            // Now we have vertex, shared vertex and face information.
            // Only vol element info is left.

            /*
            for(int i=0; i<numprocs; i++)
            {
                std::cout << vertnumber[i] << std::endl;
            }
            */

            //std::cout << sfile[0].str() << std::endl;
            //std::cout << vertnumber[0] << std::endl;
            //std::cout << facenumber[0] << std::endl;
            //std::cout << elemnumber[0] << std::endl;

            // Now we have everything needed. Send the data to the other processors.

            int nve = datanumber[0];    // vertex
            int nel = datanumber[2];    // element
            int nfa = datanumber[1];    // face
            int stringsize = stringtosend[0].length(); // size of receive string

            std::string recvstring = stringtosend[0];

            MPI_Barrier(comm);   // Barrier 6

            for(int i=1; i< numprocs; i++)
            {
                std::cout << "SENDING COUNTS TO RANK: " << i+1 << " " << datanumber[4*i] << " " << datanumber[4*i+1] << " " << datanumber[4*i+2] << " " << datanumber[4*i+3] << std::endl;
                MPI_Send(&datanumber[4*i], 4, MPI_INT, i, 123, comm);
            }

            MPI_Barrier(comm); // Barrier 7

            std::cout << "RANK " << mypid << " has received " << nve << " " << nfa << " " << nel << " " << stringsize << " " << recvstring.length() << std::endl;

            MPI_Barrier(comm); // Barrier 8 Send the strings.


            reqs = new MPI_Request[numprocs-1];
            stats = new MPI_Status[numprocs-1];

            for(int i=1; i<numprocs; i++)
            {
                std::cout << "Trying to send " << stringtosend[i].length() << " to rank " << i << std::endl;
                MPI_Isend(&stringtosend[i][0], stringtosend[i].length(), MPI_CHAR, i, 123, comm, &reqs[i-1]);
            }

            MPI_Waitall(numprocs-1, reqs, stats);

            delete [] reqs;
            delete [] stats;

            MPI_Barrier(comm); // Barrier 9

            /// SEND DATA END

            /// END DISTRIBUTE DATA
            /// Now we have the data. Fill the mesh.


            // Header
            numvertex = nve;
            numvolume = nel;
            numboundary = nfa;
            numshared = 0;      // This will be increased while reading the vertex data.
            numdirichlet = 0;   // Not known. So set to 0 for now.

            std::stringstream strread;

            strread << recvstring;


            double x, y, z;
            int glid, tmp;

            verts.resize(numvertex);
            stringstream oss;

            for (int i=0; i<numvertex; i++)
            {
                strread >> glid >> x >> y >> z >> tmp;

                verts[i].coor[0] = x;
                verts[i].coor[1] = y;
                verts[i].coor[2] = z;
                verts[i].idx = glid;
                verts[i].useCount = 0;

                vertmap[glid] = i;

                if(tmp == 1)
                {
                    // Read shared information.
                    strread >> glid;    // glid used for count of processors.

                    strread >> tmp;

                    oss << glid << " " << tmp+1;

                    for(int k=1; k<glid; k++)
                    {
                        strread >> tmp;

                        oss << " " << tmp+1;
                    }

                    sharedVert[verts[i].idx] = oss.str();

                    oss.clear();
                    oss.str("");
                    numshared++;
                }
            }

            // Header, Vertex and Shared Data is filled. Now boundary and volume data should be entered.

            int lid2, x2, y2, z2, par2, gfid2;

            boundelems.resize(numboundary);

            for (int i=0; i<numboundary; i++) {
                strread >> lid2 >> gfid2 >> par2 >> x2 >> y2 >> z2;

                boundelems[i].verts[0] = x2;
                boundelems[i].verts[1] = y2;
                boundelems[i].verts[2] = z2;
                boundelems[i].idx = lid2;
                boundelems[i].geoid = gfid2 + 1;
                boundelems[i].parentelement = par2;

                volbound[par2] = lid2;
            }

            // Now volume elements.

            int glid3, x3, y3, z3 , w3;

            volelems.resize(numvolume);

            for (int i=0; i<numvolume; i++) {
                strread >> glid3 >> x3 >> y3 >> z3 >> w3;

                volelems[i].verts[0] = x3;
                volelems[i].verts[1] = y3;
                volelems[i].verts[2] = z3;
                volelems[i].verts[3] = w3;
                volelems[i].idx = glid3;

                // increase usage of vertices.
                verts[vertmap.find(x3)->second].useCount++;
                verts[vertmap.find(y3)->second].useCount++;
                verts[vertmap.find(z3)->second].useCount++;
                verts[vertmap.find(w3)->second].useCount++;

                volmap[glid3] = i;
            }

            /// MESH FILLED


        }
        else
        {
            /// En az birer dene eleman lazim. Root gonderecek ilk n-1 elemani.


            // int *eptr2 = new int[nelems+1]; // 0 2 5 8 11 13 (ne+1 adet) Degisken olabiliyor.Ama bizim durumumuzda sadece tetrahedral var. Yani 0 4 8 12 16 ...
            // int *elemts = new int[4*nelems];        // [12 11 7 3] [11 12 8 9] Sirayla elemenlarin indexleri yaziliyor.
            // int *part = new int[nelems];

            int nelems = 0;

            MPI_Bcast(&nelems, 1, MPI_INT, 0, comm);
            MPI_Barrier(comm);

            int *eptr2 = new int[2];
            eptr2[0] = 0;
            eptr2[1] = 4;
            int *elemts = new int[4];
            int *part = new int[1];
            part[0] = mypid;        // Random for testing.

            float *ubvec = new float[1];     // balance constraint
            float *tpwgts = new float[numprocs]; // normally ncon*nparts.
            int * elmdist = new int[numprocs+1];    // analogue to vtxdist

            elmdist[0] = 0;             // Since all elements are at the first processor for now. Send the last numprocs-1 elements one by one to each other.
            elmdist[1] = nelems-numprocs+1;
            for(int i=2; i<=numprocs; i++)
            {
                elmdist[i] = elmdist[i-1] + 1;
            }

            int numflag = 0;    // indexing starts from 0 (set to 1 if index start from 1)
            int wgtflag = 0;    // not weighted
            int ncon = 1;
            int ncommonnodes = 3;   // for tetrahedral.

            int nparts = numprocs;      // size of partition.

            ubvec[0] = 1.05;        // recommended value by parmetis
            int options[10];        // options.
            int edgecut;            // edges cut in total by the new partitioning.

            options[0] = 1;
            options[PMV3_OPTION_DBGLVL] = 7;
            options[PMV3_OPTION_SEED] = 0;

            for(int kk=0; kk<numprocs; kk++)
            {
                tpwgts[kk] = 1.0/(float)(nparts);
            }


            /// CALL METIS   (We will be using Parmetis with all elements first at node 0)

            std::cout << "I AM HERE before send" << std::endl;
            MPI_Barrier(comm);
            // here get one element from root.
            MPI_Request request;
            MPI_Status status;
            MPI_Irecv(&elemts[0], 4, MPI_INT, 0, 123, comm, &request);
            MPI_Wait(&request, &status);

            MPI_Barrier(comm);
            std::cout << "RANK: " << mypid << " " << nelems << " " << elemts[0] << " " << elemts[1] << " " << elemts[2] << " " << elemts[3] << " " << std::endl;

            ParMETIS_V3_PartMeshKway(elmdist, eptr2, elemts, NULL, &wgtflag, &numflag, &ncon, &ncommonnodes, &nparts, tpwgts, ubvec, options, &edgecut, part, &comm);
            /// CALL METIS END
            // Send back the result to root.

            MPI_Barrier(comm);

            std::cout << "RANK: " << mypid << " sending: " << part[0] << std::endl;
            MPI_Send(part, 1, MPI_INT, 0, 123, comm);

            MPI_Barrier(comm);
            std::cout << "RANK: " << mypid << " => I AM HERE3" << std::endl;

            /// Get the data from the root processor.
            MPI_Barrier(comm);          // Barrier 4

            /// RECEIVE DATA
            MPI_Barrier(comm);   // Barrier 5

            int nve = 0;    // vertex
            int nel = 0;    // element
            int nfa = 0;    // face
            int stringsize = 0; // size of receive string
            int * dataa = new int[4];

            std::string recvstring;

            MPI_Barrier(comm);   // Barrier 6
            // Receive count numbers here.

            MPI_Irecv(&dataa[0], 4, MPI_INT, 0, 123, comm, &request);
            MPI_Wait(&request, &status);

            MPI_Barrier(comm); // Barrier 7
            nve = dataa[0];
            nfa = dataa[1];
            nel = dataa[2];
            stringsize = dataa[3];

            recvstring.resize(stringsize);

            std::cout << "RANK " << mypid << " has received " << nve << " " << nfa << " " << nel << " " << stringsize << std::endl;

            MPI_Barrier(comm); // Barrier 8 receive the string.

            MPI_Request reqe;
            MPI_Status state;

            std::cout << "Trying to receive " << stringsize << " My PID = " << mypid << std::endl;
            MPI_Irecv(&recvstring[0], stringsize, MPI_CHAR, 0, 123, comm, &reqe);
            MPI_Wait(&reqe, &state);

            MPI_Barrier(comm); // Barrier 9

            /// RECEIVE DATA END

            /// Now we have the data. Fill the mesh.


            // Header
            numvertex = nve;
            numvolume = nel;
            numboundary = nfa;
            numshared = 0;      // This will be increased while reading the vertex data.
            numdirichlet = 0;   // Not known. So set to 0 for now.

            std::stringstream strread;

            strread << recvstring;


            double x, y, z;
            int glid, tmp;

            verts.resize(numvertex);
            stringstream oss;

            for (int i=0; i<numvertex; i++)
            {
                strread >> glid >> x >> y >> z >> tmp;

                verts[i].coor[0] = x;
                verts[i].coor[1] = y;
                verts[i].coor[2] = z;
                verts[i].idx = glid;
                verts[i].useCount = 0;

                vertmap[glid] = i;

                if(tmp == 1)
                {
                    // Read shared information.
                    strread >> glid;    // glid used for count of processors.

                    strread >> tmp;

                    oss << glid << " " << tmp+1;

                    for(int k=1; k<glid; k++)
                    {
                        strread >> tmp;

                        oss << " " << tmp+1;
                    }

                    sharedVert[verts[i].idx] = oss.str();

                    oss.clear();
                    oss.str("");
                    numshared++;
                }
            }

            // Header, Vertex and Shared Data is filled. Now boundary and volume data should be entered.

            int lid2, x2, y2, z2, par2, gfid2;

            boundelems.resize(numboundary);

            for (int i=0; i<numboundary; i++) {
                strread >> lid2 >> gfid2 >> par2 >> x2 >> y2 >> z2;

                boundelems[i].verts[0] = x2;
                boundelems[i].verts[1] = y2;
                boundelems[i].verts[2] = z2;
                boundelems[i].idx = lid2;
                boundelems[i].geoid = gfid2 + 1;    // BC ID's begin from 1.
                boundelems[i].parentelement = par2;

                volbound[par2] = lid2;
            }

            // Now volume elements.

            int glid3, x3, y3, z3 , w3;

            volelems.resize(numvolume);

            for (int i=0; i<numvolume; i++) {
                strread >> glid3 >> x3 >> y3 >> z3 >> w3;

                volelems[i].verts[0] = x3;
                volelems[i].verts[1] = y3;
                volelems[i].verts[2] = z3;
                volelems[i].verts[3] = w3;
                volelems[i].idx = glid3;

                // increase usage of vertices.
                verts[vertmap.find(x3)->second].useCount++;
                verts[vertmap.find(y3)->second].useCount++;
                verts[vertmap.find(z3)->second].useCount++;
                verts[vertmap.find(w3)->second].useCount++;

                volmap[glid3] = i;
            }

            /// MESH FILLED
        }


        /// Process the data and create the surface mesh for each processor.

        /// Refine the 2D mesh.

        /// Generate volume mesh at each processor.
    }

    void createRandomMigData(int totdiffprocs, float percentage)
    {
        int prid;
        srand((unsigned)time(0));
		lsize = numvolume;
		elids = new int[lsize];
		procid = new int[lsize];

		int* procidss = new int[totdiffprocs];

		for(int i=0; i<totdiffprocs; i++)
		{
            procidss[i] = rand()%numprocs;
		}

        int limit = numvolume*percentage;

		for (int i=0; i<limit; i++) {
            prid = rand()%totdiffprocs;
			elids[i] = volelems[i].idx;
			procid[i] = procidss[prid];
		}

		for(int i=limit; i<numvolume; i++)
		{
			elids[i] = volelems[i].idx;
			procid[i] = mypid;
		}
    }

    void refineParallelV2(MPI_Comm comm, const char* inputFile)
    {
        /// Then refine the surface meshes in parallel.
        bool (*fn_pt)(Index2,Index2) = fncomp;
        bool (*fn_pt2)(Index2,Index2) = fncomp;
        std::map<Index2,int, bool(*)(Index2, Index2)> edges(fn_pt);
        std::map<Index2,int, bool(*)(Index2, Index2)> edgesonbound(fn_pt2);
        std::map<Index2,int, bool(*)(Index2, Index2)>::iterator myit;

        std::cout << "REFINEMENT BASLAMAKTA." << std::endl;

        Index2 i2;
        int isshared = 0;       // 0 for not shared.  1 for shared.
        int totshared = 0;
        int totnonshared = 0;
        int totall = 0;

        int maxelid = 0;
        int maxglobalid = 0;

        int curgeoid = 0;

        for(int i=0; i<numboundary; i++)
        {
            curgeoid = boundelems[i].geoid - 1; // ElmerID's are (NetgenID + 1)

            if(boundelems[i].verts[0] > maxelid)
            {
                maxelid = boundelems[i].verts[0];
            }

            if(boundelems[i].verts[1] > maxelid)
            {
                maxelid = boundelems[i].verts[1];
            }

            if(boundelems[i].verts[2] > maxelid)
            {
                maxelid = boundelems[i].verts[2];
            }

            // 1
            i2.x[0] = boundelems[i].verts[0];
            i2.x[1] = boundelems[i].verts[1];

            i2.Sort();

            /*
            if(i2.x[0] == 18104 && i2.x[1] == 35924)
            {
                std::cout << "ARADIGIMIZI BULDUK: " << mypid << " " << i2.x[0] << " " << i2.x[1] << std::endl;
            }*/

            if(sharedVert.find(i2.x[0]) == sharedVert.end() || sharedVert.find(i2.x[1]) == sharedVert.end())
            {
                isshared = 0;
            }
            else
            {
                isshared = -1;
            }

            // If edge is not already in the edge list the add it to the list.
            if(edges.find(i2) == edges.end())
            {
                edges.insert(pair<Index2,int>(i2,isshared));

                if(isshared == 0)
                {
                    totnonshared++;
                }
                else
                {
                    totshared++;
                }
            }
            if(edgesonbound.find(i2) == edgesonbound.end())
            {
                edgesonbound.insert(pair<Index2,int>(i2,curgeoid));
            }

            // 2
            i2.x[0] = boundelems[i].verts[0];
            i2.x[1] = boundelems[i].verts[2];

            i2.Sort();

            /*
            if(i2.x[0] == 18104 && i2.x[1] == 35924)
            {
                std::cout << "ARADIGIMIZI BULDUK: " << mypid << " " << i2.x[0] << " " << i2.x[1] << std::endl;
            }*/

            if(sharedVert.find(i2.x[0]) == sharedVert.end() || sharedVert.find(i2.x[1]) == sharedVert.end())
            {
                isshared = 0;
            }
            else
            {
                isshared = -1;
            }

            // If edge is not already in the edge list the add it to the list.
            if(edges.find(i2) == edges.end())
            {
                edges.insert(pair<Index2,int>(i2,isshared));

                if(isshared == 0)
                {
                    totnonshared++;
                }
                else
                {
                    totshared++;
                }
            }
            if(edgesonbound.find(i2) == edgesonbound.end())
            {
                edgesonbound.insert(pair<Index2,int>(i2,curgeoid));
            }

            // 3
            i2.x[0] = boundelems[i].verts[1];
            i2.x[1] = boundelems[i].verts[2];

            i2.Sort();

            /*
            if(i2.x[0] == 18104 && i2.x[1] == 35924)
            {
                std::cout << "ARADIGIMIZI BULDUK: " << mypid << " " << i2.x[0] << " " << i2.x[1] << std::endl;
            }*/

            if(sharedVert.find(i2.x[0]) == sharedVert.end() || sharedVert.find(i2.x[1]) == sharedVert.end())
            {
                isshared = 0;
            }
            else
            {
                isshared = -1;
            }

            // If edge is not already in the edge list the add it to the list.
            if(edges.find(i2) == edges.end())
            {
                edges.insert(pair<Index2,int>(i2,isshared));

                if(isshared == 0)
                {
                    totnonshared++;
                }
                else
                {
                    totshared++;
                }
            }
            if(edgesonbound.find(i2) == edgesonbound.end())
            {
                edgesonbound.insert(pair<Index2,int>(i2,curgeoid));
            }
        }
        int numparts = partbound.size();

        for(int i=0; i<numparts; i++)
        {
            if(partbound[i].verts[0] > maxelid)
            {
                maxelid = partbound[i].verts[0];
            }

            if(partbound[i].verts[1] > maxelid)
            {
                maxelid = partbound[i].verts[1];
            }

            if(partbound[i].verts[2] > maxelid)
            {
                maxelid = partbound[i].verts[2];
            }

            // 1
            i2.x[0] = partbound[i].verts[0];
            i2.x[1] = partbound[i].verts[1];

            i2.Sort();

            /*
            if(i2.x[0] == 18104 && i2.x[1] == 35924)
            {
                std::cout << "ARADIGIMIZI BULDUK: " << mypid << " " << i2.x[0] << " " << i2.x[1] << std::endl;
            }*/

            if(sharedVert.find(i2.x[0]) == sharedVert.end() || sharedVert.find(i2.x[1]) == sharedVert.end())
            {
                isshared = 0;
            }
            else
            {
                isshared = -1;
            }

            // If edge is not already in the edge list the add it to the list.
            if(edges.find(i2) == edges.end())
            {
                edges.insert(pair<Index2,int>(i2,isshared));

                if(isshared == 0)
                {
                    totnonshared++;
                }
                else
                {
                    totshared++;
                }
            }

            // 2
            i2.x[0] = partbound[i].verts[0];
            i2.x[1] = partbound[i].verts[2];

            i2.Sort();

            /*
            if(i2.x[0] == 18104 && i2.x[1] == 35924)
            {
                std::cout << "ARADIGIMIZI BULDUK: " << mypid << " " << i2.x[0] << " " << i2.x[1] << std::endl;
            }*/

            if(sharedVert.find(i2.x[0]) == sharedVert.end() || sharedVert.find(i2.x[1]) == sharedVert.end())
            {
                isshared = 0;
            }
            else
            {
                isshared = -1;
            }

            // If edge is not already in the edge list the add it to the list.
            if(edges.find(i2) == edges.end())
            {
                edges.insert(pair<Index2,int>(i2,isshared));

                if(isshared == 0)
                {
                    totnonshared++;
                }
                else
                {
                    totshared++;
                }
            }

            // 3
            i2.x[0] = partbound[i].verts[1];
            i2.x[1] = partbound[i].verts[2];

            i2.Sort();

            /*
            if(i2.x[0] == 18104 && i2.x[1] == 35924)
            {
                std::cout << "ARADIGIMIZI BULDUK: " << mypid << " " << i2.x[0] << " " << i2.x[1] << std::endl;
            }*/

            if(sharedVert.find(i2.x[0]) == sharedVert.end() || sharedVert.find(i2.x[1]) == sharedVert.end())
            {
                isshared = 0;
            }
            else
            {
                isshared = -1;
            }

            // If edge is not already in the edge list the add it to the list.
            if(edges.find(i2) == edges.end())
            {
                edges.insert(pair<Index2,int>(i2,isshared));

                if(isshared == 0)
                {
                    totnonshared++;
                }
                else
                {
                    totshared++;
                }
            }
        }

        totall = totshared + totnonshared;

        /// We have found the edges. And we also know whether they are fully shared by another one.

        // std::cout << "RANK: " << mypid << " SIZE OF EDGES : " << totall << std::endl;
        // std::cout << "RANK: " << mypid << " SIZE OF SHARED: " << totshared << std::endl;
        // std::cout << "NUMBER OF VOLELEMS: " << numvolume << std::endl;
        // std::cout << "NUMBER OF VERTICES: " << numvertex << std::endl;
        // std::cout << "MAX ELEMENT ID: " << maxelid << std::endl;

        /*
        for(myit = edges.begin(); myit != edges.end(); myit++)
        {
            std::cout << (*myit).first.x[0] << "   " << (*myit).first.x[1] << "  " << (*myit).second * '*' << std::endl;
        }
        */

        MPI_Barrier(comm);
        MPI_Allreduce(&maxelid, &maxglobalid, 1, MPI_INT, MPI_MAX, comm);
        MPI_Barrier(comm);
        // std::cout << "MAX GLOBAL ID: " << maxglobalid << std::endl;

        /// Now we know max global id of a vertex. We should calculate how many new vertices each processor will have
        map<int,string>::iterator sharit;
        stringstream oss;
        string tmpstr;

        set<int> procset;
        set<int> combprocset;
        set<int>::iterator procsetit;

        int tt = 0;
        int tt2 = 0;	// count;
        int tt3 = 0;	// temp storage.

        int minid = 10000000;   // Increase this if using more than 10000000 processors. Unlikely.
        int totowner = 0;
        int totnonowner = 0;
        int totpat = 0;
        bool checkit = false;
        int tottosend = 0;

        int myrcv = 0;
        int mysnd = 0;

        int *tosend55 = new int[numprocs];
        int *torecv55 = new int[numprocs];

        for(int i=0; i<numprocs; i++)
        {
            tosend55[i] = 0;
            torecv55[i] = 0;
        }

        int combsize = 0;
        bool iamowner = false;

        // std::cout << "RANK: " << mypid << std::endl;

        for(myit = edges.begin(); myit != edges.end(); myit++)
        {
            // (*myit).first.x[0]
            // (*myit).first.x[1]
            // (*myit).second

            if((*myit).second == -1) // it is shared
            {
                /// here compare the shared processor id's of the two vertices.
                /// And the smallest processor id is the owner and will assign the global id of the new vertex.

                //std::cout << (*myit).first.x[0] << "    " << sharedVert.find((*myit).first.x[0])->second << std::endl;
                //std::cout << (*myit).first.x[1] << "    " << sharedVert.find((*myit).first.x[1])->second << std::endl;

                sharit = sharedVert.find((*myit).first.x[0]);

				if (sharit != sharedVert.end()) // This should always be true.
				{
					tmpstr = (*sharit).second;
					// std::cout << tmpstr << std::endl;

					tt = (*sharit).first;

					oss << tmpstr;
					oss >> tt2;

					for (int i=0; i<tt2; i++) {
						oss >> tt3;                     // Processor ID (tt3-1) elmer rank and MPI rank are different.
						procset.insert(tt3-1);  // Map or set usage.
					}

					oss.clear();
					oss.str("");
				}
				else
				{
				    std::cout << "PROBLEM" << std::endl;
				}

				/// Now we have the processor id's from the first vertex.

				sharit = sharedVert.find((*myit).first.x[1]);

				if (sharit != sharedVert.end()) // This should always be true.
				{
					tmpstr = (*sharit).second;
					// std::cout << tmpstr << std::endl;
					tt = (*sharit).first;

					oss << tmpstr;
					oss >> tt2;

					for (int i=0; i<tt2; i++) {
						oss >> tt3;                     // Processor ID (tt3-1) elmer rank and MPI rank are different.
						if(procset.find(tt3-1) != procset.end())    // If both contain the same element.
						{
                            combprocset.insert(tt3-1);  // Map or set usage.
						}
					}

					oss.clear();
					oss.str("");
				}
				else
				{
				    std::cout << "PROBLEM" << std::endl;
				}

				/// Now we have the processor id's in combprocset which contain both the vertices.

                minid = 10000000;

                checkit = false;

				for(procsetit = combprocset.begin(); procsetit != combprocset.end(); procsetit++)
				{
                    if(*procsetit < minid)
                    {
                        minid = *procsetit;
                    }

                    if(*procsetit == mypid)
                    {
                        checkit = true;
                    }
				}

				///Uncomment this part to make minid back the smallest one again.

				/*
				int detlef = ((*myit).first.x[0] + (*myit).first.x[1]) % combprocset.size();
				for(procsetit = combprocset.begin(); procsetit != combprocset.end(); procsetit++)
				{
				    if(detlef == 0)
				    {
                        minid = *procsetit;
				    }
				    detlef--;
				}
				*/

				if( checkit == false)
				{
                    std::cout << "YOK BOYLE BISEY" << std::endl;
				}

                // std::cout << combprocset.size() << " " << minid << std::endl;

                if(minid < 0 || minid >= numprocs || combprocset.size() == 0)
                {
                    std::cout << "PATTES OLDU" << std::endl;
                }

                combsize = combprocset.size();

                if(minid == mypid)
                {
                    iamowner = true;
                }
                else
                {
                    iamowner = false;
                }

                /*
                if( (*myit).first.x[0] == 18104 && (*myit).first.x[1] == 35924 )
                {
                    std::cout << "COMBSIZE: " << combsize << " RANK: " << mypid << std::endl;
                }
                */

                if( iamowner && combsize > 1 )
                {
                    (*myit).second = -2;         // If equals -2 means it is shared and the owner is this one.
                    // std::cout << "BU OLUYO" << std::endl;
                    tottosend += (combprocset.size() - 1);
                    totowner++;

                    mysnd += (combprocset.size()-1);

                    for(procsetit = combprocset.begin(); procsetit != combprocset.end(); procsetit++)
                    {
                        if(*procsetit != mypid)
                        {
                            tosend55[(*procsetit)]++;
                        }
                    }

                }
                else if( iamowner && combsize == 1 )   // If the size of the set is 1. This is the owner
                {
                    (*myit).second = 0;         // Internal vertex

                    // tottosend++;
                    totowner++;
                    //std::cout << "BU DA OLUYO" << std::endl;
                }
                else if( combprocset.size() == 1)
                {
                    totpat++;           // Should never be.
                }
                else
                {
                    (*myit).second = -1;    // unnecessary
                    /// I am not the owner.
                    totnonowner++;

                    myrcv++;

                    torecv55[minid]++;
                }

                procset.clear();
                combprocset.clear();
            }
        }

        if(totpat != 0)
        {
            std::cout << "BI NOKTADA KOYUVERDI PROGRAM" << std::endl;
        }

        /*
        std::cout << "MY RANK: " << mypid << " sends these values: ";
        for(int i=0; i<numprocs; i++)
        {
            std::cout<< tosend55[i] << "   ";
        }
        std::cout << std::endl;

        std::cout << "MY RANK: " << mypid << " receives these values: ";
        for(int i=0; i<numprocs; i++)
        {
            std::cout<< torecv55[i] << "   ";
        }
        std::cout << std::endl;
        */

        int myvertices = totowner + totnonshared;
        MPI_Barrier(comm);

        /// Now we have calculated the number of shared vertices which

        //std::cout << "RANK: " << mypid << " TOT SENDS: " << mysnd << "   AND TOT RECVS: " << myrcv << std::endl;

        int data1 = myvertices;
        int npindex = 0;

        MPI_Scan(&data1, &npindex, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Barrier(comm);
        int startindex = npindex - myvertices + maxglobalid + 1;
        // std::cout << "RANK: " << mypid << " NEW START INDEX: " << startindex << std::endl;

        /// Now we have the startindex for the vertices globalids. We should start creating the id's and then send
        /// The global id's of our owned vertices to the other ones.

        int curindex = startindex;

        int sendcount = 0;
        int recvcount = 0;
        int *sendprids;
        int *sendcounts;
        int *recvprids;
        int *recvcounts;

        int* sendpr5;
        int* recvpr5;

        int * send2 = new int[numprocs];
        for(int i=0; i<numprocs; i++)
        {
            send2[i] = 0;
        }


        for(myit = edges.begin(); myit != edges.end(); myit++)
        {
            // (*myit).first.x[0]
            // (*myit).first.x[1]
            // (*myit).second

            if((*myit).second == -1) // it is shared and i am not the owner
            {
                // Do nothing.
                myrcv--;
            }
            else if((*myit).second == -2) // it is shared and i am the owner.
            {

                /// Now here we should find out to which processor to send and what.

                sharit = sharedVert.find((*myit).first.x[0]);

				if (sharit != sharedVert.end()) // This should always be true.
				{
					tmpstr = (*sharit).second;
					tt = (*sharit).first;

					oss << tmpstr;
					oss >> tt2;

					for (int i=0; i<tt2; i++) {
						oss >> tt3;                     // Processor ID (tt3-1) elmer rank and MPI rank are different.
						procset.insert(tt3-1);  // Map or set usage.
					}

					oss.clear();
					oss.str("");
				}
				else
				{
				    std::cout << "PROBLEM2" << std::endl;
				}

				/// Now we have the processor id's from the first vertex.

				sharit = sharedVert.find((*myit).first.x[1]);

				if (sharit != sharedVert.end()) // This should always be true.
				{
					tmpstr = (*sharit).second;
					tt = (*sharit).first;

					oss << tmpstr;
					oss >> tt2;

					for (int i=0; i<tt2; i++) {
						oss >> tt3;                     // Processor ID (tt3-1) elmer rank and MPI rank are different.
						if(procset.find(tt3-1) != procset.end())    // If both contain the same element.
						{
                            combprocset.insert(tt3-1);  // Map or set usage.
						}
					}

					oss.clear();
					oss.str("");
				}
				else
				{
				    std::cout << "PROBLEM2" << std::endl;
				}

				/// Now we have the processor id's in combprocset which contain both the vertices.

				for(procsetit = combprocset.begin(); procsetit != combprocset.end(); procsetit++)
				{
                    if(*procsetit != mypid)
                    {
                        // Add to the send list.
                        send2[*procsetit]++;
                        mysnd--;
                    }
				}

                procset.clear();
                combprocset.clear();

            }
            else if((*myit).second == 0)    // It is nonshared.
            {
                (*myit).second =  curindex;     // give a global id.
                curindex++;
            }
        }

        //std::cout << "MY RANK: " << mypid << " has RECEIVES SET BACK TO: " << myrcv << std::endl;
        //std::cout << "MY RANK: " << mypid << " has SENDS SET BACK TO: " << mysnd << std::endl;

        /// New migration find send counts START
        sendcount = 0;
        for(int i=0; i<numprocs; i++)
        {
            if(send2[i] == 0)
            {
                // Do nothing.
            }
            else
            {
                sendcount++;
            }
        }

        sendprids = new int[sendcount];
        sendcounts = new int[sendcount];

        sendcount = 0;
        int sendtotals = 0;

        for(int i=0; i<numprocs; i++)
        {
            if(send2[i] == 0)
            {
                //Do nothing
            }
            else
            {
                sendcounts[sendcount] = send2[i];
                sendtotals += send2[i];
                send2[i] = 1;
                sendprids[sendcount] = i;
                sendcount++;
            }
        }

        /// Find receive counts
        {
            int *totdiffrecvs = new int[numprocs];  // total number of different processors to receive from.
            MPI_Allreduce(send2, totdiffrecvs, numprocs,MPI_INT, MPI_SUM, comm);

            recvcount = totdiffrecvs[mypid];
            delete [] totdiffrecvs;

            recvcounts = new int[recvcount];
            recvprids = new int[recvcount];
        }

        MPI_Request * reqs = new MPI_Request[recvcount + sendcount];
        MPI_Status * stats = new MPI_Status[recvcount + sendcount];

        /// Now send to each one how much you will send.
        for(int i=0; i<sendcount; i++)
        {
            // Send number of volume elements to be send.
            MPI_Isend(&sendcounts[i], 1, MPI_INT, sendprids[i], 123, comm, &reqs[i]);
        }
        for(int i=0; i<recvcount; i++)
        {
            MPI_Irecv(&recvcounts[i], 1, MPI_INT, MPI_ANY_SOURCE, 123, comm, &reqs[sendcount+i]);
        }
        MPI_Waitall(recvcount + sendcount, reqs, stats);

        for(int i=0; i<recvcount; i++)
        {
            recvprids[i] = stats[sendcount + i].MPI_SOURCE;
        }

        delete [] reqs;
        delete [] stats;
        /// Now each processor knows how much it will receive from which processor. (recvcounts & recvprids)
        /// And how much it will send to which processor. (sendcounts & sendprids)

        /// New migration find send and recv counts END

        recvpr5 = new int[recvcount+1];
        sendpr5 = new int[sendcount+1];
        //
        sendpr5[0] = 0;
        for (int j=0; j<sendcount; j++) {
            sendpr5[j+1] = sendpr5[j] + sendcounts[j];
        }

        recvpr5[0] = 0;
        for (int i=0; i<recvcount; i++) {
            recvpr5[i+1] = recvpr5[i] + recvcounts[i];
        }

        int *currentel5 = new int[sendcount];
        for (int i=0; i<sendcount; i++) {
            currentel5[i] = 0;
        }

        /// Now we know which processor will send what, which will get what and how many from whom.
        /// We should fill the contents.

        //std::cout << "RANK: " << mypid << " has to receive: " << recvpr5[recvcount] << "   " << totshared - totowner << "    " << totnonowner << std::endl;
        //std::cout << "RANK: " << mypid << " has to send: " << sendpr5[sendcount] << "   " << tottosend << std::endl;

        int* datatosend = new int[3 * sendpr5[sendcount]];
        int curin = 0;
        int curin2 = 0;

        for(myit = edges.begin(); myit != edges.end(); myit++)
        {
            if((*myit).second == -2) // it is shared and i am the owner.
            {

                /// Now here we should find out to which processor to send and what.

                sharit = sharedVert.find((*myit).first.x[0]);

				if (sharit != sharedVert.end()) // This should always be true.
				{
					tmpstr = (*sharit).second;
					tt = (*sharit).first;

					oss << tmpstr;
					oss >> tt2;

					for (int i=0; i<tt2; i++) {
						oss >> tt3;                     // Processor ID (tt3-1) elmer rank and MPI rank are different.
						procset.insert(tt3-1);  // Map or set usage.
					}

					oss.clear();
					oss.str("");
				}

				/// Now we have the processor id's from the first vertex.

				sharit = sharedVert.find((*myit).first.x[1]);

				if (sharit != sharedVert.end()) // This should always be true.
				{
					tmpstr = (*sharit).second;
					tt = (*sharit).first;

					oss << tmpstr;
					oss >> tt2;

					for (int i=0; i<tt2; i++) {
						oss >> tt3;                     // Processor ID (tt3-1) elmer rank and MPI rank are different.
						if(procset.find(tt3-1) != procset.end())    // If both contain the same element.
						{
                            combprocset.insert(tt3-1);  // Map or set usage.
						}
					}

					oss.clear();
					oss.str("");
				}

				/// Now we have the processor id's in combprocset which contain both the vertices.

				(*myit).second =  curindex;     // give a global id.
                curindex++;

				for(procsetit = combprocset.begin(); procsetit != combprocset.end(); procsetit++)
				{
                    if(*procsetit != mypid)
                    {
                        /// Now fill the content to send.
                        curin = findEl(sendprids, *procsetit, sendcount);
                        curin2 = 3*(sendpr5[curin] + currentel5[curin]);

                        datatosend[curin2] = (*myit).first.x[0];
                        datatosend[curin2+1] = (*myit).first.x[1];
                        datatosend[curin2+2] = (*myit).second;

                        currentel5[curin]++;

                        mysnd++;
                    }
				}

                procset.clear();
                combprocset.clear();
            }
        }

        int maxcurrentlid = curindex;

        //std::cout << "MY RANK: " << mypid << " has SENDS POP UP TO: " << mysnd << std::endl;

        int* datattorecv = new int[3 * recvpr5[recvcount]];

        MPI_Barrier(comm);

        // Send Elements
        int tmpp2;

        reqs = new MPI_Request[sendcount + recvcount];
        stats = new MPI_Status[sendcount + recvcount];

        for (int i=0; i<sendcount; i++) {
            tmpp2 = 3*(sendpr5[i+1] - sendpr5[i]);
            MPI_Isend(&datatosend[3*sendpr5[i]], tmpp2, MPI_INT, sendprids[i], 123, comm, &reqs[i]);
        }

        for (int i=0; i<recvcount; i++) {
            tmpp2 = 3* (recvpr5[i+1] - recvpr5[i]);
            MPI_Irecv(&datattorecv[3*recvpr5[i]], tmpp2, MPI_INT, recvprids[i], 123, comm, &reqs[sendcount + i]);
        }
        MPI_Waitall(recvcount + sendcount, reqs, stats);

        delete [] reqs;
        delete [] stats;

        // std::cout << "RANK: " << mypid << " is here." << std::endl;
        MPI_Barrier(comm);

        /// Now fill in the new global id's for nonowned shared vertices.

        // std::cout << "RECC COUNT: " << recvpr5[recvcount] << "  and  the expected recv count is: "  << totshared - totowner << std::endl;

        int totproblemm = 0;

        std::list<Index3> excessList;
        Index3 i3;
        std::list<Index3>::iterator excessit;

        for(int i=0; i<recvpr5[recvcount]; i++)
        {
            i2.x[0] = datattorecv[3*i];
            i2.x[1] = datattorecv[3*i+1];

            i2.Sort();

            myit = edges.find(i2);

            if(myit == edges.end())
            {
                /// This here finds the number of extra sends to this processor.
                std::cout << "EXTRA: " << mypid << "   " << i2.x[0] << "   "  << i2.x[1] << "   " << datattorecv[3*i+2] << std::endl;

                i3.x[0] = i2.x[0];
                i3.x[1] = i2.x[1];
                i3.x[2] = datattorecv[3*i+2];

                excessList.push_back(i3);

                totproblemm++;
            }
            // std::cout << datattorecv[3*i+2] << std::endl;
            myit->second = datattorecv[3*i+2];
            myrcv++;
        }

        //std::cout << "MY RANK: " << mypid << " has RECVS POP UP TO: " << myrcv << std::endl;

        int totabsent = 0;

        for(myit = edges.begin(); myit != edges.end(); myit++)
        {
            if((*myit).second == -1)
            {
                totabsent++;
            }
        }

        int *ids1 = new int[totabsent];
        int *ids2 = new int[totabsent];

        int qq = 0;



        for(myit = edges.begin(); myit != edges.end(); myit++)
        {
            if((*myit).second == -1)
            {
                std::cout << "ABSENT: " << mypid << "   " << myit->first.x[0] << "   "  << myit->first.x[1] << "   " << -1 << std::endl;
                // totproblemm--;
                ids1[qq] = myit->first.x[0];
                ids2[qq] = myit->first.x[1];
                qq++;
            }
        }

        /// For processing of received information from this function.
        /// (Beware; if shared proc list has only one procesor then it is an internal vertex.)
        /// In this case no one gets extra information.
        /// Do not forget (procid = procid + 1) for elmer
        solveAbsentProblem(ids1, ids2, totabsent, maxcurrentlid);
        std::cout << "RANK: " << mypid << " received the following absent Data: " << std::endl;
        std::cout << absentData;

        /// Adding the absent global id's.
        std::stringstream oss5(stringstream::in | stringstream::out);
        oss5 << absentData;

        int a1,b1,gl2,c2, tp2;
        while(oss5.good() && oss5.str().length() != 0)
        {
            oss5 >> a1 >> b1 >> gl2 >> c2;

            for(int i=0;i<c2;i++)
            {
                oss5 >> tp2;
            }

            i2.x[0] = a1;
            i2.x[1] = b1;
            i2.Sort();
            edges.find(i2)->second = gl2;

            // std::cout << a1 << " " << b1 << " " << gl2 << std::endl;
        }

        int excesssize = excessList.size();
        int *ids3 = new int[excesssize];
        int *ids4 = new int[excesssize];
        int *idsgl = new int[excesssize];

        qq = 0;

        for(excessit = excessList.begin(); excessit != excessList.end(); excessit++)
        {
            i3 = *excessit;
            ids3[qq] = i3.x[0];
            ids4[qq] = i3.x[1];
            idsgl[qq] = i3.x[2];
            qq++;
        }

        solveExcessProblem(ids3, ids4, idsgl, excesssize);

        std::cout << "RANK: " << mypid << " received the following excess Data: " << std::endl;
        std::cout << excessData << excessData2;


        /// After both processes have completed.
        /// First update global ids of those who were -1. (uninitialized)
        /// The following parts can then proceed.
        /// The shared vertex information should be updated with the received data afterwards.

        delete [] datatosend;
        delete [] datattorecv;
        /// Now we have sent everything needed. Lets calculate what's left.

        /// Create the new vertices;
        verts.resize(numvertex + edges.size());
        int cs = 0;
        int a3 = 0;
        int b3 = 0;

        int boundedgesize = edgesonbound.size();
        double* xdeger = new double[boundedgesize];
        double* ydeger = new double[boundedgesize];
        double* zdeger = new double[boundedgesize];
        int* yuzeydeger = new int[boundedgesize];

        std::cout << "EDGES ON BOUNDARY SIZE OF RANK: " << mypid << " is: " << boundedgesize << " versus " << edges.size() << std::endl;

        /// Now project them on the surface. Get the results on xdeger, ydeger, zdeger. (This will be an NG Call).
        using namespace nglib;

        int pp1 = 0;
        int pa1 = 0;
        int pa2 = 0;

        for(myit = edgesonbound.begin(); myit != edgesonbound.end(); myit++)
        {
            pa1 = myit->first.x[0];
            pa2 = myit->first.x[1];

            pa1 = vertmap.find(pa1)->second;
            pa2 = vertmap.find(pa2)->second;

            xdeger[pp1] = (verts[pa1].coor[0] + verts[pa2].coor[0])/2;
            ydeger[pp1] = (verts[pa1].coor[1] + verts[pa2].coor[1])/2;
            zdeger[pp1] = (verts[pa1].coor[2] + verts[pa2].coor[2])/2;

            yuzeydeger[pp1] = (*myit).second;
            pp1++;
        }

        Ng_CSG_ProjectMesh(inputFile, boundedgesize, xdeger, ydeger, zdeger, yuzeydeger);

        MPI_Barrier(comm);

        pp1 = 0;
        for(myit = edgesonbound.begin(); myit != edgesonbound.end(); myit++)
        {
            (*myit).second = pp1;
            pp1++;
        }

        int totalscheus = 0;
        bool isonface = false;
        int indexxx = 0;

        for(myit = edges.begin(); myit != edges.end(); myit++)
        {
            a3 = (*myit).first.x[0];
            b3 = (*myit).first.x[1];

            i2.x[0] = a3;
            i2.x[1] = b3;

            i2.Sort();

            if(edgesonbound.find(i2) != edgesonbound.end())
            {
                isonface = true;
                indexxx = edgesonbound.find(i2)->second;
            }

            /// Create shared information
            if(sharedVert.find(a3) != sharedVert.end() || sharedVert.find(b3) != sharedVert.end())
            {
                sharit = sharedVert.find(a3);

				if (sharit != sharedVert.end()) // This should always be true.
				{
					tmpstr = (*sharit).second;
					tt = (*sharit).first;

					oss << tmpstr;
					oss >> tt2;

					for (int i=0; i<tt2; i++) {
						oss >> tt3;                     // Processor ID (tt3-1) elmer rank and MPI rank are different.
						procset.insert(tt3-1);  // Map or set usage.
					}

					oss.clear();
					oss.str("");
				}

				/// Now we have the processor id's from the first vertex.

				sharit = sharedVert.find(b3);

				if (sharit != sharedVert.end()) // This should always be true.
				{
					tmpstr = (*sharit).second;
					tt = (*sharit).first;

					oss << tmpstr;
					oss >> tt2;

					for (int i=0; i<tt2; i++) {
						oss >> tt3;                     // Processor ID (tt3-1) elmer rank and MPI rank are different.
						if(procset.find(tt3-1) != procset.end())    // If both contain the same element.
						{
                            combprocset.insert(tt3-1);  // Map or set usage.
						}
					}

					oss.clear();
					oss.str("");
				}

				/// Now we have the processor id's in combprocset which contain both the vertices.
                if(combprocset.size() > 1)
                {
                    // myit->second;
                    oss << combprocset.size();

                    for(procsetit = combprocset.begin(); procsetit != combprocset.end(); procsetit++)
                    {
                        oss << " " << *procsetit + 1;      // Elmer procid = MPIprocid + 1
                    }

                    sharedVert[(*myit).second] = oss.str();

                    if((*myit).second == -1)
                    {
                        //std::cout << "PROBLEMMPPPP!!!!  " << a3 <<  "   "  << b3 << std::endl;
                    }

                    numshared++;

                    oss.clear();
                    oss.str("");
                }

                procset.clear();
                combprocset.clear();
            }

            /// Create the new vertices

            a3 = vertmap.find(a3)->second;
            b3 = vertmap.find(b3)->second;

            // std::cout << a3 << "  " << b3 << std::endl;

            if(myit->second == 0 || myit->second == -1 || myit->second == -2)
            {
                totalscheus++;
            }

            verts[numvertex + cs].idx = myit->second;
            verts[numvertex + cs].useCount = 0;

            if(isonface != true)
            {
                verts[numvertex + cs].coor[0] = (verts[a3].coor[0] + verts[b3].coor[0])/2;
                verts[numvertex + cs].coor[1] = (verts[a3].coor[1] + verts[b3].coor[1])/2;
                verts[numvertex + cs].coor[2] = (verts[a3].coor[2] + verts[b3].coor[2])/2;
            }
            else
            {
                verts[numvertex + cs].coor[0] = xdeger[indexxx];
                verts[numvertex + cs].coor[1] = ydeger[indexxx];
                verts[numvertex + cs].coor[2] = zdeger[indexxx];
            }

            vertmap[myit->second] = numvertex + cs;

            cs++;
            isonface = false;
        }
        numvertex += edges.size();

        std::cout << "RANK: " << mypid << " has Total Sh.. " << totalscheus << std::endl;


        /// Updating absent Data related shared processor info errors.
        std::stringstream oss6(stringstream::in | stringstream::out);
        std::stringstream oss7(stringstream::in | stringstream::out);

        oss6.clear();
        oss6.str("");
        oss6 << absentData;

        oss7.clear();
        oss7.str("");

        while(!oss6.eof())
        {
            oss6 >> a1 >> b1 >> gl2 >> c2;

            if(oss6.eof())
            {
                break;
            }

            oss7 << c2;

            for(int i=0;i<c2;i++)
            {
                oss6 >> tp2;
                oss7 << " " << tp2 + 1;
            }

            tmpstr = oss7.str();

            if(c2 == 1)
            {
                // delete shared vertex from list.
                sharedVert.erase(gl2);
                // std::cout << "RANK: " << mypid << " has deleted shared info with global id: " << gl2 << std::endl;
            }
            else if(c2 > 1)
            {
                sharedVert[gl2] = tmpstr;
                // std::cout << "RANK: " << mypid << " has updated shared info with global id: " << gl2 << " to: " << tmpstr << std::endl;
            }

            oss7.clear();
            oss7.str("");
        }

        /// Update excessData
        oss6.clear();
        oss6.str("");
        oss6 << excessData;

        oss7.clear();
        oss7.str("");

        while(!oss6.eof())
        {
            oss6 >> a1 >> b1  >> c2;

            if(oss6.eof())
            {
                break;
            }

            i2.x[0] = a1;
            i2.x[1] = b1;
            i2.Sort();
            gl2 = edges.find(i2)->second;

            if(c2 == 1)
            {
                // delete shared vertex from list.
                sharedVert.erase(gl2);
            }

            oss7 << c2;

            for(int i=0;i<c2;i++)
            {
                oss6 >> tp2;
                oss7 << " " << tp2 + 1;
            }

            if(c2 > 1)
            {
                sharedVert[gl2] = oss7.str();
            }

            oss7.clear();
            oss7.str("");
        }

        /// Update excessData 2
        oss6.clear();
        oss6.str("");
        oss6 << excessData2;

        oss7.clear();
        oss7.str("");

        while(!oss6.eof())
        {
            oss6 >> a1 >> b1  >> c2;

            if(oss6.eof())
            {
                break;
            }

            i2.x[0] = a1;
            i2.x[1] = b1;
            i2.Sort();
            gl2 = edges.find(i2)->second;

            if(c2 == 1)
            {
                // delete shared vertex from list.
                sharedVert.erase(gl2);
            }

            oss7 << c2;

            for(int i=0;i<c2;i++)
            {
                oss6 >> tp2;
                oss7 << " " << tp2 + 1;
            }

            if(c2 > 1)
            {
                sharedVert[gl2] = oss7.str();
            }

            oss7.clear();
            oss7.str("");
        }

        /// Create the new boundaryelements.

        boundelems.resize(numboundary*4);

        int cind = numboundary*4-4;

        int a4,b4,c4,d4;
        int ab,ac,ad,bc,bd,cd;

        for(int i=numboundary-1; i>=0; i--)
        {
            a4 = boundelems[i].verts[0];
            b4 = boundelems[i].verts[1];
            c4 = boundelems[i].verts[2];

            i2.x[0] = a4;
            i2.x[1] = b4;
            i2.Sort();
            ab = edges.find(i2)->second;

            i2.x[0] = a4;
            i2.x[1] = c4;
            i2.Sort();
            ac = edges.find(i2)->second;

            i2.x[0] = b4;
            i2.x[1] = c4;
            i2.Sort();
            bc = edges.find(i2)->second;

            boundelems[cind].verts[0] = a4;
            boundelems[cind].verts[1] = ab;
            boundelems[cind].verts[2] = ac;

            boundelems[cind+1].verts[0] = ab;
            boundelems[cind+1].verts[1] = b4;
            boundelems[cind+1].verts[2] = bc;

            boundelems[cind+2].verts[0] = ac;
            boundelems[cind+2].verts[1] = bc;
            boundelems[cind+2].verts[2] = c4;

            boundelems[cind+3].verts[0] = ab;
            boundelems[cind+3].verts[1] = bc;
            boundelems[cind+3].verts[2] = ac;

            boundelems[cind].geoid = boundelems[i].geoid;
            boundelems[cind+1].geoid = boundelems[i].geoid;
            boundelems[cind+2].geoid = boundelems[i].geoid;
            boundelems[cind+3].geoid = boundelems[i].geoid;

            boundelems[cind].parentelement = 0;     // This will be set later.
            boundelems[cind+1].parentelement = 0;
            boundelems[cind+2].parentelement = 0;
            boundelems[cind+3].parentelement = 0;

            boundelems[cind].idx = cind+1;          // ID's are local.
            boundelems[cind+1].idx = cind+2;
            boundelems[cind+2].idx = cind+3;
            boundelems[cind+3].idx = cind+4;

            cind -= 4;
        }

        numboundary *= 4;

        /// Create the new partitionboundaryelements.

        partbound.resize(numparts*4);

        cind = numparts*4-4;

        for(int i=numparts-1; i>=0; i--)
        {
            a4 = partbound[i].verts[0];
            b4 = partbound[i].verts[1];
            c4 = partbound[i].verts[2];

            i2.x[0] = a4;
            i2.x[1] = b4;
            i2.Sort();
            ab = edges.find(i2)->second;

            i2.x[0] = a4;
            i2.x[1] = c4;
            i2.Sort();
            ac = edges.find(i2)->second;

            i2.x[0] = b4;
            i2.x[1] = c4;
            i2.Sort();
            bc = edges.find(i2)->second;

            partbound[cind].verts[0] = a4;
            partbound[cind].verts[1] = ab;
            partbound[cind].verts[2] = ac;

            partbound[cind+1].verts[0] = ab;
            partbound[cind+1].verts[1] = b4;
            partbound[cind+1].verts[2] = bc;

            partbound[cind+2].verts[0] = ac;
            partbound[cind+2].verts[1] = bc;
            partbound[cind+2].verts[2] = c4;

            partbound[cind+3].verts[0] = ab;
            partbound[cind+3].verts[1] = bc;
            partbound[cind+3].verts[2] = ac;

            cind -= 4;
        }

        numparts *= 4;



        std::cout << "TOTAL PROBLEM FACTOR OF RANK " << mypid << " IS " << totproblemm << std::endl;


        /// At last create the volume mesh from the surface meshes. In this part just do volume element numbering.

        /// createVolume();
        /// Other things should not change.
    }

    void skinMesh()
    {
        /// First skin the already created mesh and create the partboundary face data.
        int x3,y3,z3;
        int count = 0;

        for(int i=0; i<numvolume; i++)
        {
            /// 1

            x3 = volelems[i].verts[0];
            y3 = volelems[i].verts[1];
            z3 = volelems[i].verts[2];

            if(sharedVert.find(x3) == sharedVert.end() || sharedVert.find(y3) == sharedVert.end() || sharedVert.find(z3) == sharedVert.end())
            {
                // do nothing.
            }
            else
            {
                /// We have a partition boundary face.
                partbound.resize(partbound.size()+1);
                partbound[count].verts[0] = x3;
                partbound[count].verts[1] = y3;
                partbound[count].verts[2] = z3;
                count++;
            }

            /// 2

            x3 = volelems[i].verts[0];
            y3 = volelems[i].verts[1];
            z3 = volelems[i].verts[3];

            if(sharedVert.find(x3) == sharedVert.end() || sharedVert.find(y3) == sharedVert.end() || sharedVert.find(z3) == sharedVert.end())
            {
                // do nothing.
            }
            else
            {
                /// We have a partition boundary face.
                partbound.resize(partbound.size()+1);
                partbound[count].verts[0] = x3;
                partbound[count].verts[1] = y3;
                partbound[count].verts[2] = z3;
                count++;
            }

            /// 3

            x3 = volelems[i].verts[0];
            y3 = volelems[i].verts[2];
            z3 = volelems[i].verts[3];

            if(sharedVert.find(x3) == sharedVert.end() || sharedVert.find(y3) == sharedVert.end() || sharedVert.find(z3) == sharedVert.end())
            {
                // do nothing.
            }
            else
            {
                /// We have a partition boundary face.
                partbound.resize(partbound.size()+1);
                partbound[count].verts[0] = x3;
                partbound[count].verts[1] = y3;
                partbound[count].verts[2] = z3;
                count++;
            }

            /// 4

            x3 = volelems[i].verts[1];
            y3 = volelems[i].verts[2];
            z3 = volelems[i].verts[3];

            if(sharedVert.find(x3) == sharedVert.end() || sharedVert.find(y3) == sharedVert.end() || sharedVert.find(z3) == sharedVert.end())
            {
                // do nothing.
            }
            else
            {
                /// We have a partition boundary face.
                partbound.resize(partbound.size()+1);
                partbound[count].verts[0] = x3;
                partbound[count].verts[1] = y3;
                partbound[count].verts[2] = z3;
                count++;
            }
        }
    }

    void skinMeshv2()
    {
        /// Create the Index3 maps and their compare function pointers

        bool (*fn_pt)(Index3,Index3) = fncomp;
        std::map<Index3,int, bool(*)(Index3, Index3)> boundfromvol(fn_pt);
        std::map<Index3,int, bool(*)(Index3, Index3)>::iterator myit;

        bool (*fn_pt2)(Index3,Index3) = fncomp;
        std::map<Index3,int, bool(*)(Index3, Index3)> facefromvol(fn_pt2);
        std::map<Index3,int, bool(*)(Index3, Index3)>::iterator myit2;

        /// Algorithm will be as follows.
        /// First fill the boundfromvol with faces from the boundary

        Index3 i3;
        int count = 0;

        for(int i=0; i<numboundary; i++)
        {
            i3.x[0] = boundelems[i].verts[0];
            i3.x[1] = boundelems[i].verts[1];
            i3.x[2] = boundelems[i].verts[2];

            i3.Sort();

            boundfromvol.insert(pair<Index3,int>(i3,i));
        }

        std::cout << boundfromvol.size() << std::endl;

        /// Next for each face at each volume element do the following.

        for(int i=0; i<numvolume; i++)
        {
            /// 1

            i3.x[0] = volelems[i].verts[0];
            i3.x[1] = volelems[i].verts[1];
            i3.x[2] = volelems[i].verts[2];

            i3.Sort();

            if(boundfromvol.find(i3) != boundfromvol.end())
            {
                /// If it exists in boundfromvol. Don't add the face. (It is boundary face and not partition face)
                // std::cout << "ARADA OLIII" << std::endl;
            }
            else if(facefromvol.find(i3) != facefromvol.end())
            {
                /// If it exists in facefromvol. Then delete that face from there.
                /// (This means the face exists on two different volume elements. So it is an inside face and not a partition face)

                facefromvol.erase(i3);
            }
            else
            {
                /// If it does not exist in facefromvol. Then add that face there.
                facefromvol.insert(pair<Index3,int>(i3,i));
            }

            /// 2

            i3.x[0] = volelems[i].verts[0];
            i3.x[1] = volelems[i].verts[1];
            i3.x[2] = volelems[i].verts[3];

            i3.Sort();

            if(boundfromvol.find(i3) != boundfromvol.end())
            {
                /// If it exists in boundfromvol. Don't add the face. (It is boundary face and not partition face)
            }
            else if(facefromvol.find(i3) != facefromvol.end())
            {
                /// If it exists in facefromvol. Then delete that face from there.
                /// (This means the face exists on two different volume elements. So it is an inside face and not a partition face)

                facefromvol.erase(i3);
            }
            else
            {
                /// If it does not exist in facefromvol. Then add that face there.
                facefromvol.insert(pair<Index3,int>(i3,i));
            }

            /// 3

            i3.x[0] = volelems[i].verts[0];
            i3.x[1] = volelems[i].verts[2];
            i3.x[2] = volelems[i].verts[3];

            i3.Sort();

            if(boundfromvol.find(i3) != boundfromvol.end())
            {
                /// If it exists in boundfromvol. Don't add the face. (It is boundary face and not partition face)
            }
            else if(facefromvol.find(i3) != facefromvol.end())
            {
                /// If it exists in facefromvol. Then delete that face from there.
                /// (This means the face exists on two different volume elements. So it is an inside face and not a partition face)

                facefromvol.erase(i3);
            }
            else
            {
                /// If it does not exist in facefromvol. Then add that face there.
                facefromvol.insert(pair<Index3,int>(i3,i));
            }

            /// 4

            i3.x[0] = volelems[i].verts[1];
            i3.x[1] = volelems[i].verts[2];
            i3.x[2] = volelems[i].verts[3];

            i3.Sort();

            if(boundfromvol.find(i3) != boundfromvol.end())
            {
                /// If it exists in boundfromvol. Don't add the face. (It is boundary face and not partition face)
            }
            else if(facefromvol.find(i3) != facefromvol.end())
            {
                /// If it exists in facefromvol. Then delete that face from there.
                /// (This means the face exists on two different volume elements. So it is an inside face and not a partition face)

                facefromvol.erase(i3);
            }
            else
            {
                /// If it does not exist in facefromvol. Then add that face there.
                facefromvol.insert(pair<Index3,int>(i3,i));
            }
        }

        /// In the end we will have only faces which are not on the boundary and yet belong only two one volume element of
            /// the partition. So this means they are partition faces (they are on the partition boundary).

        for(int i=0; i<numvolume; i++)
        {
            /// 1

            i3.x[0] = volelems[i].verts[0];
            i3.x[1] = volelems[i].verts[1];
            i3.x[2] = volelems[i].verts[2];

            i3.Sort();

            if(facefromvol.find(i3) != facefromvol.end())
            {
                partbound.resize(partbound.size()+1);
                partbound[count].verts[0] = volelems[i].verts[2];
                partbound[count].verts[1] = volelems[i].verts[0];
                partbound[count].verts[2] = volelems[i].verts[1];
                count++;
            }

            /// 2

            i3.x[0] = volelems[i].verts[0];
            i3.x[1] = volelems[i].verts[1];
            i3.x[2] = volelems[i].verts[3];

            i3.Sort();

            if(facefromvol.find(i3) != facefromvol.end())
            {
                partbound.resize(partbound.size()+1);
                partbound[count].verts[0] = volelems[i].verts[0];
                partbound[count].verts[1] = volelems[i].verts[3];
                partbound[count].verts[2] = volelems[i].verts[1];
                count++;
            }

            /// 3

            i3.x[0] = volelems[i].verts[0];
            i3.x[1] = volelems[i].verts[2];
            i3.x[2] = volelems[i].verts[3];

            i3.Sort();

            if(facefromvol.find(i3) != facefromvol.end())
            {
                partbound.resize(partbound.size()+1);
                partbound[count].verts[0] = volelems[i].verts[2];
                partbound[count].verts[1] = volelems[i].verts[3];
                partbound[count].verts[2] = volelems[i].verts[0];
                count++;
            }

            /// 4

            i3.x[0] = volelems[i].verts[1];
            i3.x[1] = volelems[i].verts[2];
            i3.x[2] = volelems[i].verts[3];

            i3.Sort();

            if(facefromvol.find(i3) != facefromvol.end())
            {
                partbound.resize(partbound.size()+1);
                partbound[count].verts[0] = volelems[i].verts[2];
                partbound[count].verts[1] = volelems[i].verts[1];
                partbound[count].verts[2] = volelems[i].verts[3];
                count++;
            }
        }

        /// Now the partition faces together with the boundary faces form a closed surface that can be volume meshed.

    }



    void createVolume()
    {
        /// Create Volume Mesh from Boundary and Partition Face Data.

        using namespace nglib;

        Ng_Mesh * mesh;
        Ng_Init();

        // creates mesh structure
        mesh = Ng_NewMesh ();

        int np, nse, ne;
        double point[3];
        int trig[3], tet[4];

        int totalelems = boundelems.size() + partbound.size();

        std::set<int> vertids;
        std::set<int>::iterator vtit;

        std::map<int,int> tmpvertmap;

        int numpart = partbound.size();

        int a1,a2,a3;

        for(int i=0; i<numpart; i++)
        {
            a1 = partbound[i].verts[0];
            a2 = partbound[i].verts[1];
            a3 = partbound[i].verts[2];

            vertids.insert(a1);
            vertids.insert(a2);
            vertids.insert(a3);
        }

        for(int i=0; i<numboundary; i++)
        {
            a1 = boundelems[i].verts[0];
            a2 = boundelems[i].verts[1];
            a3 = boundelems[i].verts[2];

            vertids.insert(a1);
            vertids.insert(a2);
            vertids.insert(a3);
        }

        np = vertids.size();

        vector<Vertex> tmpVerts;
        tmpVerts.resize(np);

        int xx21 = 0;
        int count = 0;
        for(vtit = vertids.begin(); vtit != vertids.end(); vtit++)
        {
            xx21 = *vtit;

            tmpVerts[count] = verts[vertmap.find(xx21)->second];
            tmpvertmap[xx21] = count+1;
            count++;
        }

        if(mypid != -1)
        {

        std::cout << "START OF VERTEX ADD" << std::endl;
        for (int i = 0; i < np; i++)
        {
            point[0] = tmpVerts[i].coor[0];
            point[1] = tmpVerts[i].coor[1];
            point[2] = tmpVerts[i].coor[2];

            // std::cout << point[0] << " " << point[1] << " " << point[2] << std::endl;

            Ng_AddPoint (mesh, point);
        }
        std::cout << "END OF VERTEX ADD" << std::endl;

        std::cout << "START OF SURFACE ADD" << std::endl;
        nse = totalelems;
        for (int i = 0; i < numpart; i++)
        {
            trig[0] = tmpvertmap.find(partbound[i].verts[0])->second;
            trig[1] = tmpvertmap.find(partbound[i].verts[1])->second;
            trig[2] = tmpvertmap.find(partbound[i].verts[2])->second;

            // std::cout << trig[0] << " " << trig[1] << " " << trig[2] << std::endl;

            Ng_AddSurfaceElement (mesh, NG_TRIG, trig);
        }
        for (int i = numpart; i < nse; i++)
        {
            trig[0] = tmpvertmap.find(boundelems[i-numpart].verts[0])->second;
            trig[1] = tmpvertmap.find(boundelems[i-numpart].verts[1])->second;
            trig[2] = tmpvertmap.find(boundelems[i-numpart].verts[2])->second;

            // std::cout << trig[0] << " " << trig[1] << " " << trig[2] << std::endl;

            Ng_AddSurfaceElement (mesh, NG_TRIG, trig);
        }
        std::cout << "END OF SURFACE ADD" << std::endl;



        // generate volume mesh
        Ng_Meshing_Parameters mp;
        mp.maxh = 1e6;
        mp.fineness = 1;
        mp.secondorder = 0;

        cout << "start meshing" << endl;
        Ng_GenerateVolumeMesh (mesh, &mp);
        cout << "meshing done" << endl;


        char * filenme = new char[132];
        sprintf(filenme, "/root/partitioning.%d/part.%d.vol", numprocs, mypid+1);
        Ng_SaveMesh(mesh, filenme);

        ne = Ng_GetNE(mesh);

        std::cout << "NUMBER OF ELEMENTS BY RANK: " << mypid << " IS: " << ne << std::endl;

        /*
        MPI_Barrier(MPI_COMM_WORLD);

        int globaltotal = 0;
        int globsurftotal = 0;

        MPI_Allreduce(&ne, &globaltotal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&numboundary, &globsurftotal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);

        std::cout << "NUMBER OF ELEMENTS GLOBAL IS: " << globaltotal << std::endl;
        std::cout << "NUMBER OF SURFACE ELEMENTS GLOBAL IS: " << globsurftotal << std::endl;
        */

        }
        /*
        // volume mesh output
        np = Ng_GetNP(mesh);
        cout << "Points: " << np << endl;

        for (int i = 1; i <= np; i++)
        {
            Ng_GetPoint (mesh, i, point);
            cout << i << ": " << point[0] << " " << point[1] << " " << point[2] << endl;
        }

        ne = Ng_GetNE(mesh);
        cout << "Elements: " << ne << endl;

        for (int i = 1; i <= ne; i++)
        {
            Ng_GetVolumeElement (mesh, i, tet);
            cout << i << ": " << tet[0] << " " << tet[1] << " " << tet[2] << " " << tet[3] << endl;
        }
        */
    }

    void solveExcessProblem(int* ids3, int* ids4, int* idsgl, int excessCount)
    {
        int *tosend = new int[numprocs];

        for(int i=0; i<numprocs; i++)
        {
            tosend[i] = 0;
        }

        /*
        for(int i=0; i<excessCount; i++)
        {
            std::cout << ids3[i] << " " << ids4[i] << " " << idsgl[i] << std::endl;
        }*/

        int *destprids = new int[excessCount];
        map<int,string>::iterator sharit;
        stringstream oss;
        string tmpstr;

        set<int> procset;
        set<int> combprocset;
        set<int>::iterator procsetit;

        int tt = 0;
        int tt2 = 0;	// count;
        int tt3 = 0;	// temp storage.

        int minid = 0;

        /// Calculate destination processor ID's.
        for(int i=0; i<excessCount; i++)
        {
            sharit = sharedVert.find(ids3[i]);

            if (sharit != sharedVert.end()) // This should always be true.
			{
				tmpstr = (*sharit).second;
				// std::cout << tmpstr << std::endl;

				tt = (*sharit).first;

				oss << tmpstr;
				oss >> tt2;

				for (int k=0; k<tt2; k++) {
					oss >> tt3;                     // Processor ID (tt3-1) elmer rank and MPI rank are different.
					procset.insert(tt3-1);  // Map or set usage.
				}

				oss.clear();
				oss.str("");
			}
			else
			{
			   std::cout << "RANK: " << mypid << "PROBLEM1" << std::endl;
			}

			/// Now we have the processor id's from the first vertex.

			sharit = sharedVert.find(ids4[i]);

			if (sharit != sharedVert.end()) // This should always be true.
			{
				tmpstr = (*sharit).second;
				// std::cout << tmpstr << std::endl;
				tt = (*sharit).first;

				oss << tmpstr;
				oss >> tt2;

				for (int k=0; k<tt2; k++) {
					oss >> tt3;                     // Processor ID (tt3-1) elmer rank and MPI rank are different.
					if(procset.find(tt3-1) != procset.end())    // If both contain the same element.
					{
                        combprocset.insert(tt3-1);  // Map or set usage.
                    }
                }

				oss.clear();
				oss.str("");
			}
			else
			{
			    std::cout << "RANK: " << mypid << "PROBLEM2" << std::endl;
			}

			/// Now we have the processor id's in combprocset which contain both the vertices.

            minid = 10000000;

			for(procsetit = combprocset.begin(); procsetit != combprocset.end(); procsetit++)
			{
                if(*procsetit < minid)
                {
                    minid = *procsetit;
                }
            }

            ///Uncomment this part to make minid back the smallest one again.
			/*
				int detlef = (ids3[i]+ ids4[i]) % combprocset.size();
				for(procsetit = combprocset.begin(); procsetit != combprocset.end(); procsetit++)
				{
				    if(detlef == 0)
				    {
                        minid = *procsetit;
				    }
				    detlef--;
				}
				*/

            destprids[i] = minid;

            procset.clear();
            combprocset.clear();
        }

        /// Now send to destprids the data. (id1, id2, mypid)
        // First send how much to send.
        for(int i=0; i<excessCount; i++)
        {
            tosend[destprids[i]]++;
        }

        /*
        std::cout << "/////////////////////////////////////" << std::endl;
        std::cout << "RANK: " << mypid << " sends:";
        for(int i=0; i<numprocs; i++)
        {
             std::cout << " " << tosend[i];
        }
        std::cout << std::endl;
        std::cout << "/////////////////////////////////////" << std::endl;
        */

        // We now know who to send how much.
        int sendcount = 0;
        int recvcount = 0;
        int sendtotals = 0;
        int *sendprids;
        int *sendcounts;
        int *recvprids;
        int *recvcounts;

        int* sendpr5;
        int* recvpr5;

        int curin = 0;
        int curin2 = 0;

        sendcount = 0;
        for(int i=0; i<numprocs; i++)
        {
            if(tosend[i] == 0)
            {
                // Do nothing.
            }
            else
            {
                sendcount++;
            }
        }

        sendprids = new int[sendcount];
        sendcounts = new int[sendcount];

        sendcount = 0;


        for(int i=0; i<numprocs; i++)
        {
            if(tosend[i] == 0)
            {
                //Do nothing
            }
            else
            {
                sendcounts[sendcount] = tosend[i];
                sendtotals += tosend[i];
                tosend[i] = 1;
                sendprids[sendcount] = i;
                sendcount++;
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        /// Find receive counts
        {
            int *totdiffrecvs = new int[numprocs];  // total number of different processors to receive from.
            MPI_Allreduce(tosend, totdiffrecvs, numprocs,MPI_INT, MPI_SUM, MPI_COMM_WORLD);

            recvcount = totdiffrecvs[mypid];
            delete [] totdiffrecvs;

            recvcounts = new int[recvcount];
            recvprids = new int[recvcount];
        }

        MPI_Request * reqs = new MPI_Request[recvcount + sendcount];
        MPI_Status * stats = new MPI_Status[recvcount + sendcount];

        /// Now send to each one how much you will send.
        for(int i=0; i<sendcount; i++)
        {
            MPI_Isend(&sendcounts[i], 1, MPI_INT, sendprids[i], 123, MPI_COMM_WORLD, &reqs[i]);
        }
        for(int i=0; i<recvcount; i++)
        {
            MPI_Irecv(&recvcounts[i], 1, MPI_INT, MPI_ANY_SOURCE, 123, MPI_COMM_WORLD, &reqs[sendcount+i]);
        }
        MPI_Waitall(recvcount + sendcount, reqs, stats);

        for(int i=0; i<recvcount; i++)
        {
            recvprids[i] = stats[sendcount + i].MPI_SOURCE;
        }

        delete [] reqs;
        delete [] stats;
        /// Now each processor knows how much it will receive from which processor. (recvcounts & recvprids)

        MPI_Barrier(MPI_COMM_WORLD);

            /*

        std::cout << "/////////////////////////////////////" << std::endl;
        std::cout << "RANK: " << mypid << " receives: " << std::endl;
        for(int i=0; i<recvcount; i++)
        {
             std::cout << recvcounts[i] << " from " << recvprids[i] << std::endl;
        }
        std::cout << "/////////////////////////////////////" << std::endl;
        */

        /// Okay. It works until now.
        /// Now to actually send the data. The receive data should be 3*totreceives. (id1, id2, mypid) * #element
        recvpr5 = new int[recvcount+1];
        sendpr5 = new int[sendcount+1];
        //
        sendpr5[0] = 0;
        for (int j=0; j<sendcount; j++) {
            sendpr5[j+1] = sendpr5[j] + sendcounts[j];
        }

        recvpr5[0] = 0;
        for (int i=0; i<recvcount; i++) {
            recvpr5[i+1] = recvpr5[i] + recvcounts[i];
        }

        int *currentel5 = new int[sendcount];
        for (int i=0; i<sendcount; i++) {
            currentel5[i] = 0;
        }

        int* datatosend = new int[3 * sendpr5[sendcount]];
        int* datattorecv = new int[3 * recvpr5[recvcount]];

        /// Here fill the data to send;

        for(int i=0; i<excessCount; i++)
        {
            /// Now fill the content to send.
            curin = findEl(sendprids, destprids[i], sendcount);
            curin2 = 3*(sendpr5[curin] + currentel5[curin]);

            datatosend[curin2] = ids3[i];
            datatosend[curin2+1] = ids4[i];
            datatosend[curin2+2] = mypid;

            currentel5[curin]++;
        }

        /// Data to send is filled.


        MPI_Barrier(MPI_COMM_WORLD);

        // Send Elements
        int tmpp2;

        reqs = new MPI_Request[sendcount + recvcount];
        stats = new MPI_Status[sendcount + recvcount];

        for (int i=0; i<sendcount; i++) {
            tmpp2 = 3*(sendpr5[i+1] - sendpr5[i]);
            MPI_Isend(&datatosend[3*sendpr5[i]], tmpp2, MPI_INT, sendprids[i], 123, MPI_COMM_WORLD, &reqs[i]);
        }

        for (int i=0; i<recvcount; i++) {
            tmpp2 = 3* (recvpr5[i+1] - recvpr5[i]);
            MPI_Irecv(&datattorecv[3*recvpr5[i]], tmpp2, MPI_INT, recvprids[i], 123, MPI_COMM_WORLD, &reqs[sendcount + i]);
        }
        MPI_Waitall(recvcount + sendcount, reqs, stats);

        delete [] reqs;
        delete [] stats;

        MPI_Barrier(MPI_COMM_WORLD);

        /*

        std::cout << "/////////////////////////////////////" << std::endl;
        std::cout << "RANK: " << mypid << " HAS RECEIVED: " << std::endl;
        for(int i=0; i<recvpr5[recvcount]; i++)
        {
             std::cout << datattorecv[3*i] << " " << datattorecv[3*i+1] << " from " << datattorecv[3*i+2] << std::endl;
        }
        std::cout << "/////////////////////////////////////" << std::endl;
        */

        bool (*fn_pt)(Index2,Index2) = fncomp;
        std::map<Index2,int, bool(*)(Index2, Index2)> edges(fn_pt);
        std::map<Index2,int, bool(*)(Index2, Index2)>::iterator myit;

        Index2 i2;

        /// Calculate number of different new vertices.
        for(int i=0; i<recvpr5[recvcount]; i++)
        {
            i2.x[0] = datattorecv[3*i];
            i2.x[1] = datattorecv[3*i+1];

            i2.Sort();
            if(edges.find(i2) == edges.end())
            {
                edges.insert(pair<Index2,int>(i2,1));
            }
            else
            {
                edges.find(i2)->second++;
            }
        }

        int numnew = edges.size();

        // std::cout << "DIFFERENT VERTICES: " << numnew << " at rank: " << mypid << std::endl;

        MPI_Barrier(MPI_COMM_WORLD);

        /// Now we received all the info about which are send by excess


        std::set<int> proctosend;
        std::set<int>* prsets = new std::set<int>[numnew];
        std::set<int>::iterator setit;
        std::set<int>::iterator setit2;

        int i23 = 0;

        for(myit = edges.begin(); myit != edges.end(); myit++)
        {
            myit->second = i23;
            i23++;
        }

        i23 = 0;

        for(myit = edges.begin(); myit != edges.end(); myit++)
        {
            sharit = sharedVert.find(myit->first.x[0]);

            if (sharit != sharedVert.end()) // This should always be true.
			{
				tmpstr = (*sharit).second;
				// std::cout << tmpstr << std::endl;

				tt = (*sharit).first;

				oss << tmpstr;
				oss >> tt2;

				for (int k=0; k<tt2; k++) {
					oss >> tt3;                     // Processor ID (tt3-1) elmer rank and MPI rank are different.
					procset.insert(tt3-1);  // Map or set usage.
				}

				oss.clear();
				oss.str("");
			}
			else
			{
			   std::cout << "PROBLEM" << std::endl;
			}

			/// Now we have the processor id's from the first vertex.

			sharit = sharedVert.find(myit->first.x[1]);

			if (sharit != sharedVert.end()) // This should always be true.
			{
				tmpstr = (*sharit).second;
				// std::cout << tmpstr << std::endl;
				tt = (*sharit).first;

				oss << tmpstr;
				oss >> tt2;

				for (int k=0; k<tt2; k++) {
					oss >> tt3;                     // Processor ID (tt3-1) elmer rank and MPI rank are different.
					if(procset.find(tt3-1) != procset.end())    // If both contain the same element.
					{
                        combprocset.insert(tt3-1);  // Map or set usage.
                    }
                }

				oss.clear();
				oss.str("");
			}
			else
			{
			    std::cout << "PROBLEM" << std::endl;
			}

			/// Now we have the processor id's in combprocset which contain both the vertices.

            minid = 10000000;

			for(procsetit = combprocset.begin(); procsetit != combprocset.end(); procsetit++)
			{
                prsets[i23].insert(*procsetit);
            }

            procset.clear();
            combprocset.clear();
            i23++;
        }


        /// Surda insert degil delete yap. Yukarda once doldur prsets'i.
        for(int i=0; i<recvpr5[recvcount]; i++)
        {
            i2.x[0] = datattorecv[3*i];
            i2.x[1] = datattorecv[3*i+1];

            i2.Sort();
            curin = edges.find(i2)->second;

            prsets[curin].erase(datattorecv[3*i+2]);
        }

        for(int i=0; i<numnew; i++)
        {
            for(setit = prsets[i].begin(); setit != prsets[i].end(); setit++)
            {
                if(mypid != *setit)
                    proctosend.insert(*setit);
            }
        }


        int numdifprocs = proctosend.size();

        std::string* strlist = new std::string[numdifprocs];
        std::stringstream* osslist = new std::stringstream[numdifprocs];
        std::list<std::string>::iterator listit;
        int* sendprocids = new int[numdifprocs];

        i23 = 0;
        for(setit = proctosend.begin(); setit != proctosend.end(); setit++)
        {
            sendprocids[i23] = *setit;
            i23++;
        }

        int sendid = 0;

        /// He also sets which processor will be the actual owner.
        /// He then sends to each processor the info: (id1, id2, listp1, listp2, listp3 ... , ';' ) listp1 is the new owner.

        int i33 = 0;
        for(myit = edges.begin(); myit != edges.end(); myit++)
        {
            for(setit = prsets[i33].begin(); setit != prsets[i33].end(); setit++)
            {
                sendid = *setit;

                if(sendid != mypid)
                {
                    curin = findEl(sendprocids, sendid, numdifprocs);
                }

                //oss << "Sending to: " << sendid << " the following: " << myit->first.x[0] << " " << myit->first.x[1] << " " << myit->first.x[2]
                //    << " " << prsets[i33].size();

                oss << myit->first.x[0] << " " << myit->first.x[1] << " " << prsets[i33].size();

                for(setit2 = prsets[i33].begin(); setit2 != prsets[i33].end(); setit2++)
                {
                    oss << " " << *setit2;
                }
                oss << std::endl;

                tmpstr = oss.str();
                // std::cout << tmpstr;

                if(sendid != mypid)
                    strlist[curin].append(tmpstr);
                else
                {
                    excessData2.append(tmpstr);
                }

                oss.clear();
                oss.str("");
            }
            i33++;
        }


        /// Buraya geldiysen tamam.



        /*
        for(int i=0; i<numdifprocs; i++)
        {
            std::cout << "Sending to Processor: " << sendprocids[i] << std::endl;
            std::cout << strlist[i] << std::endl;
        }*/

        /*
        for(int i=0; i<numnew; i++)
        {
            std::cout << "RANK: " << mypid << " NEW ELEMENT: " << glids[i] << std::endl;
        }
        */

        /// Actually send the data to
        /// excessData private variable.
        /// CORRECT SO FAR.

        for(int i=0; i<numprocs; i++)
        {
            tosend[i] = 0;
        }

        for(int i=0; i<numdifprocs; i++)
        {
            tosend[sendprocids[i]] = 1;
        }

        MPI_Barrier(MPI_COMM_WORLD);

        /// Find receive counts
        {
            int *totdiffrecvs = new int[numprocs];  // total number of different processors to receive from.
            MPI_Allreduce(tosend, totdiffrecvs, numprocs,MPI_INT, MPI_SUM, MPI_COMM_WORLD);

            recvcount = totdiffrecvs[mypid];
            delete [] totdiffrecvs;

            recvcounts = new int[recvcount];
            recvprids = new int[recvcount];
        }

        reqs = new MPI_Request[recvcount + numdifprocs];
        stats = new MPI_Status[recvcount + numdifprocs];

        int ll = 0;

        /// Now send to each one how much you will send.
        for(int i=0; i<numdifprocs; i++)
        {
            ll = strlist[i].length();
            MPI_Isend(&ll, 1, MPI_INT, sendprocids[i], 123, MPI_COMM_WORLD, &reqs[i]);
        }
        for(int i=0; i<recvcount; i++)
        {
            MPI_Irecv(&recvcounts[i], 1, MPI_INT, MPI_ANY_SOURCE, 123, MPI_COMM_WORLD, &reqs[numdifprocs+i]);
        }
        MPI_Waitall(recvcount + numdifprocs, reqs, stats);

        for(int i=0; i<recvcount; i++)
        {
            recvprids[i] = stats[numdifprocs + i].MPI_SOURCE;
        }

        delete [] reqs;
        delete [] stats;
        /// Now each processor knows how much it will receive from which processor. (recvcounts & recvprids)
        MPI_Barrier(MPI_COMM_WORLD);

        /*
        std::cout << "/////////////////////////////////////" << std::endl;
        std::cout << "RANK: " << mypid << " receives: " << std::endl;
        for(int i=0; i<recvcount; i++)
        {
             std::cout << recvcounts[i] << " from " << recvprids[i] << std::endl;
        }
        std::cout << "/////////////////////////////////////" << std::endl;

        */

        sendcount = numdifprocs;
        /// Okay. It works until now.
        recvpr5 = new int[recvcount+1];
        sendpr5 = new int[sendcount+1];
        //
        sendpr5[0] = 0;
        for (int j=0; j<sendcount; j++) {
            sendpr5[j+1] = sendpr5[j] + strlist[j].length();
        }

        recvpr5[0] = 0;
        for (int i=0; i<recvcount; i++) {
            recvpr5[i+1] = recvpr5[i] + recvcounts[i];
        }

        currentel5 = new int[sendcount];
        for (int i=0; i<sendcount; i++) {
            currentel5[i] = 0;
        }

        char* datatosend2 = new char[sendpr5[sendcount]];
        char* datattorecv2 = new char[recvpr5[recvcount]];

        excessData.resize(recvpr5[recvcount]);

        /// Here fill the data to send;

        for(int i=0; i<numdifprocs; i++)
        {
            /// Now fill the content to send.
            ll = strlist[i].length();
            tmpstr = strlist[i];
            curin2 = sendpr5[i];
            for(int k=0; k<ll; k++)
            {
                datatosend2[curin2+k] = tmpstr[k];
            }

        }

        /// Data to send is filled.

        // std::cout << "TEST1" << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);

        // Send Elements

        reqs = new MPI_Request[sendcount + recvcount];
        stats = new MPI_Status[sendcount + recvcount];

        int counnd = 0;

        for (int i=0; i<sendcount; i++) {
            tmpp2 = sendpr5[i+1] - sendpr5[i];
            // std::cout << mypid << " tries to send: " << tmpp2 << " sized string to: " << sendprocids[i] << std::endl;

            if(sendprocids[i] == mypid)
            {
                for(int k=0; k<tmpp2; k++)
                {
                    excessData[recvpr5[i]+k] = datatosend2[sendpr5[i]+k];
                }
                continue;
            }

            if(tmpp2 == 0)
                continue;

            MPI_Isend(&datatosend2[sendpr5[i]], tmpp2, MPI_CHAR, sendprocids[i], 123, MPI_COMM_WORLD, &reqs[counnd]);
            counnd++;
        }

        for (int i=0; i<recvcount; i++) {
            tmpp2 = recvpr5[i+1] - recvpr5[i];
            // std::cout << mypid << " tries to receive: " << tmpp2 << " sized string from: " << recvprids[i] << std::endl;

            if(recvprids[i] == mypid)
                continue;

            if(tmpp2 == 0)
                continue;

            MPI_Irecv(&excessData[recvpr5[i]], tmpp2, MPI_CHAR, recvprids[i], 123, MPI_COMM_WORLD, &reqs[counnd]);
            counnd++;
        }
        MPI_Waitall(counnd, reqs, stats);

        delete [] reqs;
        delete [] stats;

        MPI_Barrier(MPI_COMM_WORLD);
        // std::cout << "TEST2" << std::endl;

        /// Now test that the data is there.
    }

    void solveAbsentProblem(int *ids1, int *ids2, int absentCount, int maxcurrentlid)
    {
        int *tosend = new int[numprocs];

        for(int i=0; i<numprocs; i++)
        {
            tosend[i] = 0;
        }

        int *destprids = new int[absentCount];
        map<int,string>::iterator sharit;
        stringstream oss;
        string tmpstr;

        set<int> procset;
        set<int> combprocset;
        set<int>::iterator procsetit;

        int tt = 0;
        int tt2 = 0;	// count;
        int tt3 = 0;	// temp storage.

        int minid = 0;

        /// Calculate destination processor ID's.
        for(int i=0; i<absentCount; i++)
        {
            sharit = sharedVert.find(ids1[i]);

            if (sharit != sharedVert.end()) // This should always be true.
			{
				tmpstr = (*sharit).second;
				// std::cout << tmpstr << std::endl;

				tt = (*sharit).first;

				oss << tmpstr;
				oss >> tt2;

				for (int k=0; k<tt2; k++) {
					oss >> tt3;                     // Processor ID (tt3-1) elmer rank and MPI rank are different.
					procset.insert(tt3-1);  // Map or set usage.
				}

				oss.clear();
				oss.str("");
			}
			else
			{
			   std::cout << "PROBLEM" << std::endl;
			}

			/// Now we have the processor id's from the first vertex.

			sharit = sharedVert.find(ids2[i]);

			if (sharit != sharedVert.end()) // This should always be true.
			{
				tmpstr = (*sharit).second;
				// std::cout << tmpstr << std::endl;
				tt = (*sharit).first;

				oss << tmpstr;
				oss >> tt2;

				for (int k=0; k<tt2; k++) {
					oss >> tt3;                     // Processor ID (tt3-1) elmer rank and MPI rank are different.
					if(procset.find(tt3-1) != procset.end())    // If both contain the same element.
					{
                        combprocset.insert(tt3-1);  // Map or set usage.
                    }
                }

				oss.clear();
				oss.str("");
			}
			else
			{
			    std::cout << "PROBLEM" << std::endl;
			}

			/// Now we have the processor id's in combprocset which contain both the vertices.

            minid = 10000000;

			for(procsetit = combprocset.begin(); procsetit != combprocset.end(); procsetit++)
			{
                if(*procsetit < minid)
                {
                    minid = *procsetit;
                }
            }

            ///Uncomment this part to make minid back the smallest one again.

			/*
				int detlef = (ids1[i] + ids2[i]) % combprocset.size();
				for(procsetit = combprocset.begin(); procsetit != combprocset.end(); procsetit++)
				{
				    if(detlef == 0)
				    {
                        minid = *procsetit;
				    }
				    detlef--;
				}
				*/

            destprids[i] = minid;

            procset.clear();
            combprocset.clear();
        }

        /// Now send to destprids the data. (id1, id2, mypid)
        // First send how much to send.
        for(int i=0; i<absentCount; i++)
        {
            tosend[destprids[i]]++;
        }

        /*
        std::cout << "/////////////////////////////////////" << std::endl;
        std::cout << "RANK: " << mypid << " sends:";
        for(int i=0; i<numprocs; i++)
        {
             std::cout << " " << tosend[i];
        }
        std::cout << std::endl;
        std::cout << "/////////////////////////////////////" << std::endl;
        */

        // We now know who to send how much.
        int sendcount = 0;
        int recvcount = 0;
        int sendtotals = 0;
        int *sendprids;
        int *sendcounts;
        int *recvprids;
        int *recvcounts;

        int* sendpr5;
        int* recvpr5;

        int curin = 0;
        int curin2 = 0;

        sendcount = 0;
        for(int i=0; i<numprocs; i++)
        {
            if(tosend[i] == 0)
            {
                // Do nothing.
            }
            else
            {
                sendcount++;
            }
        }

        sendprids = new int[sendcount];
        sendcounts = new int[sendcount];

        sendcount = 0;


        for(int i=0; i<numprocs; i++)
        {
            if(tosend[i] == 0)
            {
                //Do nothing
            }
            else
            {
                sendcounts[sendcount] = tosend[i];
                sendtotals += tosend[i];
                tosend[i] = 1;
                sendprids[sendcount] = i;
                sendcount++;
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        /// Find receive counts
        {
            int *totdiffrecvs = new int[numprocs];  // total number of different processors to receive from.
            MPI_Allreduce(tosend, totdiffrecvs, numprocs,MPI_INT, MPI_SUM, MPI_COMM_WORLD);

            recvcount = totdiffrecvs[mypid];
            delete [] totdiffrecvs;

            recvcounts = new int[recvcount];
            recvprids = new int[recvcount];
        }

        MPI_Request * reqs = new MPI_Request[recvcount + sendcount];
        MPI_Status * stats = new MPI_Status[recvcount + sendcount];

        /// Now send to each one how much you will send.
        for(int i=0; i<sendcount; i++)
        {
            // Send number of volume elements to be send.
            MPI_Isend(&sendcounts[i], 1, MPI_INT, sendprids[i], 123, MPI_COMM_WORLD, &reqs[i]);
        }
        for(int i=0; i<recvcount; i++)
        {
            MPI_Irecv(&recvcounts[i], 1, MPI_INT, MPI_ANY_SOURCE, 123, MPI_COMM_WORLD, &reqs[sendcount+i]);
        }
        MPI_Waitall(recvcount + sendcount, reqs, stats);

        for(int i=0; i<recvcount; i++)
        {
            recvprids[i] = stats[sendcount + i].MPI_SOURCE;
        }

        delete [] reqs;
        delete [] stats;
        /// Now each processor knows how much it will receive from which processor. (recvcounts & recvprids)

        MPI_Barrier(MPI_COMM_WORLD);

        /*
        std::cout << "/////////////////////////////////////" << std::endl;
        std::cout << "RANK: " << mypid << " receives: " << std::endl;
        for(int i=0; i<recvcount; i++)
        {
             std::cout << recvcounts[i] << " from " << recvprids[i] << std::endl;
        }
        std::cout << "/////////////////////////////////////" << std::endl;
        */

        /// Okay. It works until now.
        /// Now to actually send the data. The receive data should be 3*totreceives. (id1, id2, mypid) * #element
        recvpr5 = new int[recvcount+1];
        sendpr5 = new int[sendcount+1];
        //
        sendpr5[0] = 0;
        for (int j=0; j<sendcount; j++) {
            sendpr5[j+1] = sendpr5[j] + sendcounts[j];
        }

        recvpr5[0] = 0;
        for (int i=0; i<recvcount; i++) {
            recvpr5[i+1] = recvpr5[i] + recvcounts[i];
        }

        int *currentel5 = new int[sendcount];
        for (int i=0; i<sendcount; i++) {
            currentel5[i] = 0;
        }

        int* datatosend = new int[3 * sendpr5[sendcount]];
        int* datattorecv = new int[3 * recvpr5[recvcount]];

        /// Here fill the data to send;

        for(int i=0; i<absentCount; i++)
        {
            /// Now fill the content to send.
            curin = findEl(sendprids, destprids[i], sendcount);
            curin2 = 3*(sendpr5[curin] + currentel5[curin]);

            datatosend[curin2] = ids1[i];
            datatosend[curin2+1] = ids2[i];
            datatosend[curin2+2] = mypid;

            currentel5[curin]++;
        }

        /// Data to send is filled.


        MPI_Barrier(MPI_COMM_WORLD);

        // Send Elements
        int tmpp2;

        reqs = new MPI_Request[sendcount + recvcount];
        stats = new MPI_Status[sendcount + recvcount];

        for (int i=0; i<sendcount; i++) {
            tmpp2 = 3*(sendpr5[i+1] - sendpr5[i]);
            MPI_Isend(&datatosend[3*sendpr5[i]], tmpp2, MPI_INT, sendprids[i], 123, MPI_COMM_WORLD, &reqs[i]);
        }

        for (int i=0; i<recvcount; i++) {
            tmpp2 = 3* (recvpr5[i+1] - recvpr5[i]);
            MPI_Irecv(&datattorecv[3*recvpr5[i]], tmpp2, MPI_INT, recvprids[i], 123, MPI_COMM_WORLD, &reqs[sendcount + i]);
        }
        MPI_Waitall(recvcount + sendcount, reqs, stats);

        delete [] reqs;
        delete [] stats;

        MPI_Barrier(MPI_COMM_WORLD);

        /*
        std::cout << "/////////////////////////////////////" << std::endl;
        std::cout << "RANK: " << mypid << " HAS RECEIVED: " << std::endl;
        for(int i=0; i<recvpr5[recvcount]; i++)
        {
             std::cout << datattorecv[3*i] << " " << datattorecv[3*i+1] << " from " << datattorecv[3*i+2] << std::endl;
        }
        std::cout << "/////////////////////////////////////" << std::endl;
        */

        bool (*fn_pt)(Index2,Index2) = fncomp;
        std::map<Index2,int, bool(*)(Index2, Index2)> edges(fn_pt);
        std::map<Index2,int, bool(*)(Index2, Index2)>::iterator myit;

        Index2 i2;

        /// Calculate number of different new vertices.
        for(int i=0; i<recvpr5[recvcount]; i++)
        {
            i2.x[0] = datattorecv[3*i];
            i2.x[1] = datattorecv[3*i+1];

            i2.Sort();
            if(edges.find(i2) == edges.end())
            {
                edges.insert(pair<Index2,int>(i2,1));
            }
            else
            {
                edges.find(i2)->second++;
            }
        }

        int numnew = edges.size();
        int maxglobalid = 0;

        // std::cout << "DIFFERENT VERTICES: " << numnew << " at rank: " << mypid << std::endl;

        /// Now calculate the new start index for globalids.

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(&maxcurrentlid, &maxglobalid, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        // std::cout << "MAX GLOBAL ID: " << maxglobalid << std::endl;

        MPI_Barrier(MPI_COMM_WORLD);

        int data1 = numnew;
        int npindex = 0;

        MPI_Scan(&data1, &npindex, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        int startindex = npindex - numnew + maxglobalid + 1;
        // std::cout << "RANK: " << mypid << " NEW START INDEX: " << startindex << std::endl;

        int* glids = new int[numnew];

        int i22=0;
        /// Global ID calculation.
        for( myit = edges.begin(); myit != edges.end(); myit++)
        {
            glids[i22] = startindex;

            myit->second = startindex;

            startindex++;
            i22++;
        }

        std::set<int> proctosend;
        std::set<int>* prsets = new std::set<int>[numnew];
        std::set<int>::iterator setit;
        std::set<int>::iterator setit2;

        for(int i=0; i<recvpr5[recvcount]; i++)
        {
            i2.x[0] = datattorecv[3*i];
            i2.x[1] = datattorecv[3*i+1];

            i2.Sort();
            curin = findEl(glids, edges.find(i2)->second, numnew);

            prsets[curin].insert(datattorecv[3*i+2]);
            proctosend.insert(datattorecv[3*i+2]);
        }

        int numdifprocs = proctosend.size();

        std::string* strlist = new std::string[numdifprocs];
        std::stringstream* osslist = new std::stringstream[numdifprocs];
        std::list<std::string>::iterator listit;
        int* sendprocids = new int[numdifprocs];

        i22 = 0;
        for(setit = proctosend.begin(); setit != proctosend.end(); setit++)
        {
            sendprocids[i22] = *setit;
            i22++;
        }

        int sendid = 0;

        /// He also sets which processor will be the actual owner.
        /// He then sends to each processor the info: (id1, id2, glid, listp1, listp2, listp3 ... , ';' ) listp1 is the new owner.

        int i33 = 0;
        for(myit = edges.begin(); myit != edges.end(); myit++)
        {
            for(setit = prsets[i33].begin(); setit != prsets[i33].end(); setit++)
            {
                sendid = *setit;

                curin = findEl(sendprocids, sendid, numdifprocs);

                //oss << "Sending to: " << sendid << " the following: " << myit->first.x[0] << " " << myit->first.x[1] << " " << myit->first.x[2]
                //    << " " << prsets[i33].size();

                oss << myit->first.x[0] << " " << myit->first.x[1] << " " << myit->second << " " << prsets[i33].size();

                for(setit2 = prsets[i33].begin(); setit2 != prsets[i33].end(); setit2++)
                {
                    oss << " " << *setit2;
                }
                oss << std::endl;

                tmpstr = oss.str();
                // std::cout << tmpstr << " " << curin;

                strlist[curin].append(tmpstr);

                oss.clear();
                oss.str("");
            }
            i33++;
        }

        /*
        for(int i=0; i<numdifprocs; i++)
        {
            std::cout << "Sending to Processor: " << sendprocids[i] << std::endl;
            std::cout << strlist[i] << std::endl;
        }
        */


        /*
        for(int i=0; i<numnew; i++)
        {
            std::cout << "RANK: " << mypid << " NEW ELEMENT: " << glids[i] << std::endl;
        }
        */

        /// Actually send the data to
        /// absentData private variable.

        for(int i=0; i<numprocs; i++)
        {
            tosend[i] = 0;
        }

        for(int i=0; i<numdifprocs; i++)
        {
            tosend[sendprocids[i]] = 1;
        }

        MPI_Barrier(MPI_COMM_WORLD);

        /// Find receive counts
        {
            int *totdiffrecvs = new int[numprocs];  // total number of different processors to receive from.
            MPI_Allreduce(tosend, totdiffrecvs, numprocs,MPI_INT, MPI_SUM, MPI_COMM_WORLD);

            recvcount = totdiffrecvs[mypid];
            delete [] totdiffrecvs;

            recvcounts = new int[recvcount];
            recvprids = new int[recvcount];
        }

        reqs = new MPI_Request[recvcount + numdifprocs];
        stats = new MPI_Status[recvcount + numdifprocs];

        int ll = 0;
        /// Now send to each one how much you will send.
        for(int i=0; i<numdifprocs; i++)
        {
            ll = strlist[i].length();
            // Send number of volume elements to be send.
            MPI_Isend(&ll, 1, MPI_INT, sendprocids[i], 123, MPI_COMM_WORLD, &reqs[i]);
        }
        for(int i=0; i<recvcount; i++)
        {
            MPI_Irecv(&recvcounts[i], 1, MPI_INT, MPI_ANY_SOURCE, 123, MPI_COMM_WORLD, &reqs[numdifprocs+i]);
        }
        MPI_Waitall(recvcount + numdifprocs, reqs, stats);

        for(int i=0; i<recvcount; i++)
        {
            recvprids[i] = stats[numdifprocs + i].MPI_SOURCE;
        }

        delete [] reqs;
        delete [] stats;
        /// Now each processor knows how much it will receive from which processor. (recvcounts & recvprids)

        MPI_Barrier(MPI_COMM_WORLD);

        /*
        std::cout << "/////////////////////////////////////" << std::endl;
        std::cout << "RANK: " << mypid << " receives: " << std::endl;
        for(int i=0; i<recvcount; i++)
        {
             std::cout << recvcounts[i] << " from " << recvprids[i] << std::endl;
        }
        std::cout << "/////////////////////////////////////" << std::endl;
        */


        sendcount = numdifprocs;
        /// Okay. It works until now.
        /// Now to actually send the data. The receive data should be 3*totreceives. (id1, id2, mypid) * #element
        recvpr5 = new int[recvcount+1];
        sendpr5 = new int[sendcount+1];
        //
        sendpr5[0] = 0;
        for (int j=0; j<sendcount; j++) {
            sendpr5[j+1] = sendpr5[j] + strlist[j].length();
        }

        recvpr5[0] = 0;
        for (int i=0; i<recvcount; i++) {
            recvpr5[i+1] = recvpr5[i] + recvcounts[i];
        }

        currentel5 = new int[sendcount];
        for (int i=0; i<sendcount; i++) {
            currentel5[i] = 0;
        }

        char* datatosend2 = new char[sendpr5[sendcount]];
        char* datattorecv2 = new char[recvpr5[recvcount]];

        absentData.resize(recvpr5[recvcount]);

        /// Here fill the data to send;

        for(int i=0; i<numdifprocs; i++)
        {
            /// Now fill the content to send.
            ll = strlist[i].length();
            tmpstr = strlist[i];
            curin2 = sendpr5[i];
            for(int k=0; k<ll; k++)
            {
                datatosend2[curin2+k] = tmpstr[k];
            }

        }

        /// Data to send is filled.


        MPI_Barrier(MPI_COMM_WORLD);

        // Send Elements

        reqs = new MPI_Request[sendcount + recvcount];
        stats = new MPI_Status[sendcount + recvcount];

        for (int i=0; i<sendcount; i++) {
            tmpp2 = sendpr5[i+1] - sendpr5[i];
            MPI_Isend(&datatosend2[sendpr5[i]], tmpp2, MPI_CHAR, sendprocids[i], 123, MPI_COMM_WORLD, &reqs[i]);
        }

        for (int i=0; i<recvcount; i++) {
            tmpp2 = recvpr5[i+1] - recvpr5[i];
            MPI_Irecv(&absentData[recvpr5[i]], tmpp2, MPI_CHAR, recvprids[i], 123, MPI_COMM_WORLD, &reqs[sendcount + i]);
        }
        MPI_Waitall(recvcount + sendcount, reqs, stats);

        delete [] reqs;
        delete [] stats;

        MPI_Barrier(MPI_COMM_WORLD);

        /// Now test that the data is there.
    }

} ;
