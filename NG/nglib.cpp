/**************************************************************************/
/* File:   nglib.cpp                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   7. May. 2000                                                   */
/**************************************************************************/

/*
  
  Interface to the netgen meshing kernel
  
*/


#include <mystdlib.h>
#include <myadt.hpp>

#include <linalg.hpp>
#include <csg.hpp>
#include <stlgeom.hpp>
#include <occgeom.hpp>
#include <geometry2d.hpp>
#include <meshing.hpp>

#ifdef OCCGEOMETRY
#include <occgeom.hpp>
#endif


namespace netgen {
   extern void MeshFromSpline2D (SplineGeometry2d & geometry,
                                 Mesh *& mesh, 
                                 MeshingParameters & mp);
   extern CSGeometry * ParseCSG (istream & istr);
}


#ifdef PARALLEL
#include <mpi.h>

namespace netgen
{
   int id, ntasks;

   MPI_Group MPI_HIGHORDER_WORLD;
   MPI_Comm MPI_HIGHORDER_COMM;
}
#endif


/*
// should not be needed (occ currently requires it)
namespace netgen {
#include "../libsrc/visualization/vispar.hpp"
  VisualizationParameters vispar;
  VisualizationParameters :: VisualizationParameters() { ; }
}
*/


namespace nglib {
#include "nglib.h"
}

using namespace netgen;

// constants and types:

namespace nglib
{
   // initialize, deconstruct Netgen library:
   DLL_HEADER void Ng_Init ()
   {
      mycout = &cout;
      myerr = &cerr;
      // netgen::testout->SetOutStream (new ofstream ("test.out"));
      testout = new ofstream ("test.out");
   }



   // Clean-up functions before ending usage of nglib
   DLL_HEADER void Ng_Exit ()
   {
      ;
   }



   // Create a new netgen mesh object
   DLL_HEADER Ng_Mesh * Ng_NewMesh ()
   {
      Mesh * mesh = new Mesh;  
      mesh->AddFaceDescriptor (FaceDescriptor (1, 1, 0, 1));
      return (Ng_Mesh*)(void*)mesh;
   }



   // Delete an existing netgen mesh object
   DLL_HEADER void Ng_DeleteMesh (Ng_Mesh * mesh)
   {
      if(mesh != NULL)
      {
         // Delete the Mesh structures
         ((Mesh*)mesh)->DeleteMesh();

         // Now delete the Mesh class itself
         delete (Mesh*)mesh;

         // Set the Ng_Mesh pointer to NULL
         mesh = NULL;
      }
   }



   // Save a netgen mesh in the native VOL format 
   DLL_HEADER void Ng_SaveMesh(Ng_Mesh * mesh, const char* filename)
   {
      ((Mesh*)mesh)->Save(filename);
   }



   // Load a netgen native VOL mesh from a given file
   DLL_HEADER Ng_Mesh * Ng_LoadMesh(const char* filename)
   {
      Mesh * mesh = new Mesh;
      mesh->Load(filename);
      return ( (Ng_Mesh*)mesh );
   }



   // Manually add a point to an existing mesh object
   DLL_HEADER void Ng_AddPoint (Ng_Mesh * mesh, double * x)
   {
      Mesh * m = (Mesh*)mesh;
      m->AddPoint (Point3d (x[0], x[1], x[2]));
   }

   void Ng_SetPointIndex(Ng_Mesh * mesh, const int ei, const int index)
   {
      Mesh * m = (Mesh*)mesh;
      m->SurfaceElement(ei).SetIndex(index);
   }

   // Manually add a surface element of a given type to an existing mesh object
   DLL_HEADER void Ng_AddSurfaceElement (Ng_Mesh * mesh, Ng_Surface_Element_Type et,
                                         int * pi)
   {
      Mesh * m = (Mesh*)mesh;
      Element2d el (3);
      el.SetIndex (1);
      el.PNum(1) = pi[0];
      el.PNum(2) = pi[1];
      el.PNum(3) = pi[2];
      m->AddSurfaceElement (el);
   }

   void Ng_SetSurfaceElementIndex(Ng_Mesh * mesh, const int ei, const int index)
   {
      Mesh * m = (Mesh*)mesh;
      m->SurfaceElement(ei).SetIndex(index);
   }

   // Manually add a volume element of a given type to an existing mesh object
   DLL_HEADER void Ng_AddVolumeElement (Ng_Mesh * mesh, Ng_Volume_Element_Type et,
                                        int * pi)
   {
      Mesh * m = (Mesh*)mesh;
      Element el (4);
      el.SetIndex (1);
      el.PNum(1) = pi[0];
      el.PNum(2) = pi[1];
      el.PNum(3) = pi[2];
      el.PNum(4) = pi[3];
      m->AddVolumeElement (el);
   }

   void Ng_SetVolumeElementIndex(Ng_Mesh * mesh, const int ei, const int index)
   {
      Mesh * m = (Mesh*)mesh;
      m->VolumeElement(ei).SetIndex(index);
   }


   // Obtain the number of points in the mesh
   DLL_HEADER int Ng_GetNP (Ng_Mesh * mesh)
   {
      return ((Mesh*)mesh) -> GetNP();
   }



   // Obtain the number of surface elements in the mesh
   DLL_HEADER int Ng_GetNSE (Ng_Mesh * mesh)
   {
      return ((Mesh*)mesh) -> GetNSE();
   }



   // Obtain the number of volume elements in the mesh
   DLL_HEADER int Ng_GetNE (Ng_Mesh * mesh)
   {
      return ((Mesh*)mesh) -> GetNE();
   }



   //  Return point coordinates of a given point index in the mesh
   DLL_HEADER void Ng_GetPoint (Ng_Mesh * mesh, int num, double * x)
   {
      const Point3d & p = ((Mesh*)mesh)->Point(num);
      x[0] = p.X();
      x[1] = p.Y();
      x[2] = p.Z();
   }

   DLL_HEADER void GetMyPoint (Mesh * mesh, int num, double * x)
   {
      const Point3d & p = mesh->Point(num);
      x[0] = p.X();
      x[1] = p.Y();
      x[2] = p.Z();
   }

   // Return the surface element at a given index "pi"
   DLL_HEADER Ng_Surface_Element_Type 
      Ng_GetSurfaceElement (Ng_Mesh * mesh, int num, int * pi)
   {
      const Element2d & el = ((Mesh*)mesh)->SurfaceElement(num);
      for (int i = 1; i <= el.GetNP(); i++)
         pi[i-1] = el.PNum(i);
      Ng_Surface_Element_Type et;
      switch (el.GetNP())
      {
      case 3: et = NG_TRIG; break;
      case 4: et = NG_QUAD; break;
      case 6: et = NG_TRIG6; break;
      default:
         et = NG_TRIG; break; // for the compiler
      }
      return et;
   }



   // Return the volume element at a given index "pi"
   DLL_HEADER Ng_Volume_Element_Type
      Ng_GetVolumeElement (Ng_Mesh * mesh, int num, int * pi)
   {
      const Element & el = ((Mesh*)mesh)->VolumeElement(num);
      for (int i = 1; i <= el.GetNP(); i++)
         pi[i-1] = el.PNum(i);
      Ng_Volume_Element_Type et;
      switch (el.GetNP())
      {
      case 4: et = NG_TET; break;
      case 5: et = NG_PYRAMID; break;
      case 6: et = NG_PRISM; break;
      case 10: et = NG_TET10; break;
      default:
         et = NG_TET; break; // for the compiler
      }
      return et;
   }

   // Return the volume element at a given index "pi"
   DLL_HEADER void
      GetMyVolumeElement (Mesh * mesh, int num, int * pi)
   {
      const Element & el = mesh->VolumeElement(num);
      for (int i = 1; i <= el.GetNP(); i++)
         pi[i-1] = el.PNum(i);
   }

   DLL_HEADER void
      GetMySurfaceElement (Mesh * mesh, int num, int * pi)
   {
      const Element2d & el = mesh->SurfaceElement(num);
      for (int i = 1; i <= el.GetNP(); i++)
         pi[i-1] = el.PNum(i);
   }

   DLL_HEADER int GetBoundaryID(Mesh *mesh, int id)
   {
      const Element2d & el = mesh->SurfaceElement(id);
      return mesh->GetFaceDescriptor(el.GetIndex()).SurfNr();
   }

   DLL_HEADER int GetDomIn(Mesh *mesh, int id)
   {
      const Element2d & el = mesh->SurfaceElement(id);
      return mesh->GetFaceDescriptor(el.GetIndex()).DomainIn();
   }

   DLL_HEADER int GetDomOut(Mesh *mesh, int id)
   {
      const Element2d & el = mesh->SurfaceElement(id);
      return mesh->GetFaceDescriptor(el.GetIndex()).DomainOut();
   }

   DLL_HEADER void SetBoundaryID(Mesh *mesh, int triid, int newid)
   {
      const Element2d & el = mesh->SurfaceElement(triid);
      mesh->GetFaceDescriptor(el.GetIndex()).SetSurfNr(newid);
   }

   // Set a global limit on the maximum mesh size allowed
   DLL_HEADER void Ng_RestrictMeshSizeGlobal (Ng_Mesh * mesh, double h)
   {
      ((Mesh*)mesh) -> SetGlobalH (h);
   }



   // Set a local limit on the maximum mesh size allowed around the given point
   DLL_HEADER void Ng_RestrictMeshSizePoint (Ng_Mesh * mesh, double * p, double h)
   {
      ((Mesh*)mesh) -> RestrictLocalH (Point3d (p[0], p[1], p[2]), h);
   }



   // Set a local limit on the maximum mesh size allowed within a given box region
   DLL_HEADER void Ng_RestrictMeshSizeBox (Ng_Mesh * mesh, double * pmin, double * pmax, double h)
   {
      for (double x = pmin[0]; x < pmax[0]; x += h)
         for (double y = pmin[1]; y < pmax[1]; y += h)
            for (double z = pmin[2]; z < pmax[2]; z += h)
               ((Mesh*)mesh) -> RestrictLocalH (Point3d (x, y, z), h);
   }

   DLL_HEADER void Ng_GetParentElement2 (Ng_Mesh* mesh, int ei, int *dd)
   {
      *dd = 1;
      if (ei <= ((Mesh*)mesh)->mlparentsurfaceelement.Size())
         *dd = ((Mesh*)mesh)->mlparentsurfaceelement.Get(ei);
      else
         *dd = 2;
   }

   DLL_HEADER void Ng_GetVertex_SurfaceElements(Ng_Mesh* mesh, int selnr, int *elems1, int* elems2 )
   {
	const MeshTopology& topology = ((Mesh*)mesh)->GetTopology();
	int x1, x2;	

	topology.GetSurface2VolumeElement( selnr, x1, x2 );

	*elems1 = x1;
	*elems2 = x2;
   } 
   DLL_HEADER int Ng_GetVertex_SurfaceElements2(Ng_Mesh* mesh, int selnr)
   {
	const MeshTopology& topology = ((Mesh*)mesh)->GetTopology();
	int x1, x2;	
	const Element2d & el = ((Mesh*)mesh)->SurfaceElement(selnr);	

	cout << "HERE I AM " << el.GetIndex() << " " << selnr << endl;
	topology.GetSurface2VolumeElement( el.GetIndex(), x1, x2 );
	
	return x1;
   }


   // Generates volume mesh from an existing surface mesh
   DLL_HEADER Ng_Result Ng_GenerateVolumeMesh (Ng_Mesh * mesh, Ng_Meshing_Parameters * mp)
   {
      Mesh * m = (Mesh*)mesh;

      // Philippose - 30/08/2009
      // Do not locally re-define "mparam" here... "mparam" is a global 
      // object 
      //MeshingParameters mparam;
      mparam.maxh = mp->maxh;
      mparam.meshsizefilename = mp->meshsize_filename;

      double fineness = min(1., max(0., mp->fineness));
      mparam.curvaturesafety = 0.3 + 5 * fineness;
      mparam.segmentsperedge = 0.3 + 5 * fineness;

      m->CalcLocalH();

      MeshVolume (mparam, *m);
      RemoveIllegalElements (*m);
      OptimizeVolume (mparam, *m);

      return NG_OK;
   }

Array<SpecialPoint> specpoints;
static Array<MeshPoint> spoints;

DLL_HEADER Ng_Result Ng_CSG_FindPoints(CSGeometry* geom, Mesh * mesh)
{
	// Find Points
	PrintMessage (1, "Start Findpoints");
	for (int i = 0; i < geom->GetNUserPoints(); i++)
	{
		mesh->AddPoint(geom->GetUserPoint (i));
		mesh->Points().Last().Singularity (geom->GetUserPointRefFactor(i));
		mesh->AddLockedPoint (PointIndex (i+1));
	}
	
	SpecialPointCalculation spc;
	
	spc.SetIdEps(geom->GetIdEps());
	
	if (spoints.Size() == 0)
	spc.CalcSpecialPoints (*geom, spoints);
	
	PrintMessage (2, "Analyze spec points");
	spc.AnalyzeSpecialPoints (*geom, spoints, specpoints);
	
	PrintMessage (5, "done");
	// Find Points Done.

	
	// FindEdges (geom, *mesh, true);
	// MeshSurface (geom, *mesh);
}
DLL_HEADER Ng_Result Ng_CSG_GenerateSurfaceMesh (const char * filename, int minSize, long *nvertices, long *ntriangles, 
double *& vertexlist, int *& trianglelist, int *& triangleBC, int *&triangleDIN, int *&triangleDOUT)
{
	AutoPtr<CSGeometry> geometry (new CSGeometry(""));
	ifstream infile (filename);
	geometry.Reset( netgen::ParseCSG(infile) );
	geometry -> FindIdenticSurfaces(1e-6);
	Box<3> box (geometry->BoundingBox());
	double detail = 0.001;
      	double facets = 20;
     	geometry->CalcTriangleApproximation(box, detail, facets);
	
	// Geometry Created.
	AutoPtr<Mesh> mesh;
	
	// Now call the generate mesh function here.
	
	int res = geometry-> GenerateMesh (mesh.Ptr(), 1, 3, NULL);
	
	// Mesh generated.

	int np, nse;
   	double point[3];
   	int trig[3];

	// Refine until element number > minSize

	cout << "Hello World!" << endl;
	cout << "surface elements before refinement: " << mesh->GetNSE() << endl;
	cout << "elements before refinement: " << mesh->GetNE() << endl;
	cout << "points   before refinement: " << mesh->GetNP() << endl;
	while( mesh->GetNSE() < minSize)
	{
		geometry->GetRefinement().Refine ( *mesh.Ptr());
	}
	cout << "surface elements after refinement: " << mesh->GetNSE() << endl;
	cout << "points   after refinement: " << mesh->GetNP() << endl;

	// Now get that data as output.	

	// Output Mesh Data
	np = mesh->GetNP();
	*nvertices = np;
	// vertexlist = new double[np*3];
	cout << "Points Data Entered: " << np << endl;

	vertexlist = new double[np*3];
	
	for (int i = 1; i <= np; i++)
	{
		GetMyPoint (mesh.Ptr(), i, point);
		// cout << i << ": " << point[0] << " " << point[1] << " " << point[2] << endl;
		vertexlist[(i-1)*3] = point[0];
		vertexlist[(i-1)*3+1] = point[1];
		vertexlist[(i-1)*3+2] = point[2];
	}

	nse = mesh->GetNSE();
	*ntriangles = nse;
	// trianglelist = new int[np*3];
	trianglelist = new int[nse*3];
	triangleBC = new int[nse];
	triangleDIN = new int[nse];
	triangleDOUT = new int[nse];

        bool inverse = false;

	cout << "Surface Element Data Entered: " << nse << endl;
	for (int i = 1; i <= nse; i++)
	{
		GetMySurfaceElement (mesh.Ptr(), i, trig);
		// cout << i << ": " << trig[0] << " " << trig[1] << " " << trig[2] << endl;
                
                trianglelist[(i-1)*3] = trig[0];
                trianglelist[(i-1)*3+1] = trig[1];
                trianglelist[(i-1)*3+2] = trig[2];
                
                triangleBC[i-1] = GetBoundaryID(mesh.Ptr(), i); 		// Here put Boundary Element Index
		triangleDIN[i-1] = GetDomIn(mesh.Ptr(),i);
		triangleDOUT[i-1] = GetDomOut(mesh.Ptr(),i);

	}
	// Mesh Data Output Complete.

	return NG_OK;
}

DLL_HEADER Ng_Result Ng_CSG_GenerateVolumeMesh (const char * filename, int minSize, long *nvertices, long *ntriangles, long *nelems,
double *& vertexlist, int *& trianglelist, int *& triangleBC, int *&triangleDIN, int *&triangleDOUT, int *&elemlist)
{
	
	AutoPtr<CSGeometry> geometry (new CSGeometry(""));
	ifstream infile (filename);
	geometry.Reset( netgen::ParseCSG(infile) );
	geometry -> FindIdenticSurfaces(1e-6);
	Box<3> box (geometry->BoundingBox());
	double detail = 0.001;
      	double facets = 20;
     	geometry->CalcTriangleApproximation(box, detail, facets);

	
	/// New Part
	
	//AutoPtr<OCCGeometry> occgeometry;
	//OCCGeometry * occgeometry;
	//occgeometry = LoadOCC_IGES (filename);
	// occgeometry = LoadOCC_STEP (filename);
	
	
	
	/// New Part
	
	// Geometry Created.
	AutoPtr<Mesh> mesh;
	
	// Now call the generate mesh function here.
	
	int res = geometry-> GenerateMesh (mesh.Ptr(), 1, 6, NULL);
	
	// Mesh generated.

	int np, nse, ne;
   	double point[3];
   	int trig[3];
	int tet[4];

	// Refine until element number > minSize

	cout << "Hello World!" << endl;
	cout << "surface elements before refinement: " << mesh->GetNSE() << endl;
	cout << "elements before refinement: " << mesh->GetNE() << endl;
	cout << "points   before refinement: " << mesh->GetNP() << endl;
	while( mesh->GetNE() < minSize)
	{
		geometry->GetRefinement().Refine ( *mesh.Ptr());
	}
	cout << "surface elements after refinement: " << mesh->GetNSE() << endl;
        cout << "elements after refinement: " << mesh->GetNE() << endl;
	cout << "points   after refinement: " << mesh->GetNP() << endl;

	// Now get that data as output.	

	// Output Mesh Data
	np = mesh->GetNP();
	*nvertices = np;
	// vertexlist = new double[np*3];
	cout << "Points Data Entered: " << np << endl;

	vertexlist = new double[np*3];
	
	for (int i = 1; i <= np; i++)
	{
		GetMyPoint (mesh.Ptr(), i, point);
		// cout << i << ": " << point[0] << " " << point[1] << " " << point[2] << endl;
		vertexlist[(i-1)*3] = point[0];
		vertexlist[(i-1)*3+1] = point[1];
		vertexlist[(i-1)*3+2] = point[2];
	}

	nse = mesh->GetNSE();
	*ntriangles = nse;
	// trianglelist = new int[np*3];
	trianglelist = new int[nse*3];
	triangleBC = new int[nse];
	triangleDIN = new int[nse];
	triangleDOUT = new int[nse];

        bool inverse = false;

	cout << "Surface Element Data Entered: " << nse << endl;
	for (int i = 1; i <= nse; i++)
	{
		GetMySurfaceElement (mesh.Ptr(), i, trig);
		// cout << i << ": " << trig[0] << " " << trig[1] << " " << trig[2] << endl;
                
                trianglelist[(i-1)*3] = trig[0];
                trianglelist[(i-1)*3+1] = trig[1];
                trianglelist[(i-1)*3+2] = trig[2];
                
                triangleBC[i-1] = GetBoundaryID(mesh.Ptr(), i); 		// Here put Boundary Element Index
		triangleDIN[i-1] = GetDomIn(mesh.Ptr(),i);
		triangleDOUT[i-1] = GetDomOut(mesh.Ptr(),i);

	}
	
	ne = mesh->GetNE();
	*nelems = ne;
        elemlist = new int[ne*4];

        cout << "Volume Element Data Entered: " << ne << endl;

	for(int i=1; i<=ne; i++)
	{
		GetMyVolumeElement(mesh.Ptr(), i, tet);		

                elemlist[(i-1)*4] = tet[0];
                elemlist[(i-1)*4+1] = tet[1];
                elemlist[(i-1)*4+2] = tet[2];
                elemlist[(i-1)*4+3] = tet[3];
	}

	// Mesh Data Output Complete.
	
	return NG_OK;
}

DLL_HEADER Ng_Result Ng_CSG_ProjectMesh (const char * filename, int nelems, double *xdeger, double *ydeger, double *zdeger, int *yuzeydeger)
{
	AutoPtr<CSGeometry> geometry (new CSGeometry(""));
	ifstream infile (filename);
	geometry.Reset( netgen::ParseCSG(infile) );
	geometry -> FindIdenticSurfaces(1e-6);
	Box<3> box (geometry->BoundingBox());
	double detail = 0.001;
      	double facets = 20;
     	geometry->CalcTriangleApproximation(box, detail, facets);
	
	// Geometry Created.
	AutoPtr<Mesh> mesh;
	
	// Now call the generate mesh function here.
	
	int res = geometry-> GenerateMesh (mesh.Ptr(), 1, 3, NULL);

	// Surface mesh is generated. Now project the elements.
	Point<3> hnewp;	

	double tmp;

	for (int i=0; i<nelems; i++)
	{
  		hnewp(0) = xdeger[i];
		hnewp(1) = ydeger[i];
		hnewp(2) = zdeger[i];

		// tmp = xdeger[i];		

  		if (yuzeydeger[i] != -1)
    		{
      			geometry->GetSurface (yuzeydeger[i]) -> Project (hnewp);
    		}
		xdeger[i] = hnewp(0);
		ydeger[i] = hnewp(1);
		zdeger[i] = hnewp(2);
		
		// if(tmp != xdeger[i])
		//	std::cout << tmp << " becomes " << xdeger[i] << std::endl;
	}
}


DLL_HEADER Ng_Result Ng_CSG_GenerateVolumeMesh (const char * filename, int minSize, long *nvertices, long *ntriangles, long *nelems)
{
	AutoPtr<CSGeometry> geometry (new CSGeometry(""));
	ifstream infile (filename);
	geometry.Reset( netgen::ParseCSG(infile) );
	geometry -> FindIdenticSurfaces(1e-6);
	Box<3> box (geometry->BoundingBox());
	double detail = 0.001;
      	double facets = 20;
     	geometry->CalcTriangleApproximation(box, detail, facets);
	
	// Geometry Created.
	AutoPtr<Mesh> mesh;
	
	// Now call the generate mesh function here.
	
	int res = geometry-> GenerateMesh (mesh.Ptr(), 1, 6, NULL);
	
	// Mesh generated.

	int np, nse, ne;
   	double point[3];
   	int trig[3];
	int tet[4];

	// Refine until element number > minSize

	cout << "Hello World!" << endl;
	cout << "surface elements before refinement: " << mesh->GetNSE() << endl;
	cout << "elements before refinement: " << mesh->GetNE() << endl;
	cout << "points   before refinement: " << mesh->GetNP() << endl;
	while( mesh->GetNE() < minSize)
	{
		geometry->GetRefinement().Refine ( *mesh.Ptr());
	}
	cout << "surface elements after refinement: " << mesh->GetNSE() << endl;
        cout << "elements after refinement: " << mesh->GetNE() << endl;
	cout << "points   after refinement: " << mesh->GetNP() << endl;

	// Now get that data as output.	

	// Output Mesh Data
	np = mesh->GetNP();
	*nvertices = np;

	nse = mesh->GetNSE();
	*ntriangles = nse;
	
	ne = mesh->GetNE();
	*nelems = ne;
	
	return NG_OK;
}

DLL_HEADER Ng_Result Ng_CSG_GenerateVolumeMesh (const char * filename, int minSize, double *& vertexlist, int *& trianglelist, int *& triangleBC, int *&triangleDIN, int *&triangleDOUT, int *&elemlist)
{
	AutoPtr<CSGeometry> geometry (new CSGeometry(""));
	ifstream infile (filename);
	geometry.Reset( netgen::ParseCSG(infile) );
	geometry -> FindIdenticSurfaces(1e-6);
	Box<3> box (geometry->BoundingBox());
	double detail = 0.001;
      	double facets = 20;
     	geometry->CalcTriangleApproximation(box, detail, facets);
	
	// Geometry Created.
	AutoPtr<Mesh> mesh;
	
	// Now call the generate mesh function here.
	
	int res = geometry-> GenerateMesh (mesh.Ptr(), 1, 6, NULL);
	
	// Mesh generated.

	int np, nse, ne;
   	double point[3];
   	int trig[3];
	int tet[4];

	// Refine until element number > minSize

	cout << "Hello World!" << endl;
	cout << "surface elements before refinement: " << mesh->GetNSE() << endl;
	cout << "elements before refinement: " << mesh->GetNE() << endl;
	cout << "points   before refinement: " << mesh->GetNP() << endl;
	while( mesh->GetNE() < minSize)
	{
		geometry->GetRefinement().Refine ( *mesh.Ptr());
	}
	cout << "surface elements after refinement: " << mesh->GetNSE() << endl;
        cout << "elements after refinement: " << mesh->GetNE() << endl;
	cout << "points   after refinement: " << mesh->GetNP() << endl;

	// Now get that data as output.	

	// Output Mesh Data
	np = mesh->GetNP();
	cout << "Points Data Entered: " << np << endl;

	vertexlist = new double[np*3];
	
	for (int i = 1; i <= np; i++)
	{
		GetMyPoint (mesh.Ptr(), i, point);
		// cout << i << ": " << point[0] << " " << point[1] << " " << point[2] << endl;
		vertexlist[(i-1)*3] = point[0];
		vertexlist[(i-1)*3+1] = point[1];
		vertexlist[(i-1)*3+2] = point[2];
	}

	nse = mesh->GetNSE();

        bool inverse = false;

	cout << "Surface Element Data Entered: " << nse << endl;
	for (int i = 1; i <= nse; i++)
	{
		GetMySurfaceElement (mesh.Ptr(), i, trig);
		// cout << i << ": " << trig[0] << " " << trig[1] << " " << trig[2] << endl;
                
                trianglelist[(i-1)*3] = trig[0];
                trianglelist[(i-1)*3+1] = trig[1];
                trianglelist[(i-1)*3+2] = trig[2];
                
                triangleBC[i-1] = GetBoundaryID(mesh.Ptr(), i); 		// Here put Boundary Element Index
		triangleDIN[i-1] = GetDomIn(mesh.Ptr(),i);
		triangleDOUT[i-1] = GetDomOut(mesh.Ptr(),i);

	}
	
	ne = mesh->GetNE();

        cout << "Volume Element Data Entered: " << ne << endl;

	for(int i=1; i<=ne; i++)
	{
		GetMyVolumeElement(mesh.Ptr(), i, tet);		

                elemlist[(i-1)*4] = tet[0];
                elemlist[(i-1)*4+1] = tet[1];
                elemlist[(i-1)*4+2] = tet[2];
                elemlist[(i-1)*4+3] = tet[3];
	}

	// Mesh Data Output Complete.
	
	return NG_OK;
}


DLL_HEADER Ng_Result Ng_CSG_GenerateMesh (const char * filename, long *nvertices, long *ntriangles, long *ntetrahedrons,
double * vertexlist, int * trianglelist, int * tetralist)
{
	AutoPtr<CSGeometry> geometry (new CSGeometry(""));
	ifstream infile (filename);
	geometry.Reset( netgen::ParseCSG(infile) );
	geometry -> FindIdenticSurfaces(1e-6);
	Box<3> box (geometry->BoundingBox());
	double detail = 0.001;
      	double facets = 20;
     	geometry->CalcTriangleApproximation(box, detail, facets);
	
	// Geometry Created.
	AutoPtr<Mesh> mesh;

	// Now call the generate mesh function.	
	int res = geometry-> GenerateMesh (mesh.Ptr(), 1, 6, NULL);

	// Mesh generated.

	int np, ne, nse;
   	double point[3];
   	int trig[3], tet[4];

	// Refine until element number > 100000


	cout << "elements before refinement: " << mesh->GetNE() << endl;
	cout << "points   before refinement: " << mesh->GetNP() << endl;
	while( mesh->GetNE() < 80000)
	{
		geometry->GetRefinement().Refine ( *mesh.Ptr());
	}
	cout << "elements after refinement: " << mesh->GetNE() << endl;
	cout << "points   after refinement: " << mesh->GetNP() << endl;

	// Now get that data as output.	

	// Output Mesh Data
	np = mesh->GetNP();
	*nvertices = np;
	// vertexlist = new double[np*3];
	cout << "Points Data Entered: " << np << endl;
	
	for (int i = 1; i <= np; i++)
	{
		GetMyPoint (mesh.Ptr(), i, point);
		// cout << i << ": " << point[0] << " " << point[1] << " " << point[2] << endl;
		vertexlist[(i-1)*3] = point[0];
		vertexlist[(i-1)*3+1] = point[1];
		vertexlist[(i-1)*3+2] = point[2];
	}

	nse = mesh->GetNSE();
	*ntriangles = nse;
	// trianglelist = new int[np*3];
	cout << "Surface Element Data Entered: " << nse << endl;
	for (int i = 1; i <= nse; i++)
	{
		GetMySurfaceElement (mesh.Ptr(), i, trig);
		// cout << i << ": " << trig[0] << " " << trig[1] << " " << trig[2] << endl;
		trianglelist[(i-1)*3] = trig[0];
		trianglelist[(i-1)*3+1] = trig[1];
		trianglelist[(i-1)*3+2] = trig[2];
	}
	
	ne = mesh->GetNE();
	*ntetrahedrons = ne;
	// tetralist = new int[ne*4];
	cout << "Volume Element Data Entered: " << ne << endl;
	for (int i = 1; i <= ne; i++)
	{
		GetMyVolumeElement (mesh.Ptr(), i, tet);
		// cout << i << ": " << tet[0] << " " << tet[1] << " " << tet[2] << " " << tet[3] << endl;
		tetralist[(i-1)*4] = tet[0];
		tetralist[(i-1)*4+1] = tet[1];
		tetralist[(i-1)*4+2] = tet[2];
		tetralist[(i-1)*4+3] = tet[3];	
	}
	// Mesh Data Output Complete.
	
	
	return NG_OK;
}

   /* ------------------ 2D Meshing Functions ------------------------- */



   DLL_HEADER void Ng_AddPoint_2D (Ng_Mesh * mesh, double * x)
   {
      Mesh * m = (Mesh*)mesh;

      m->AddPoint (Point3d (x[0], x[1], 0));
   }

   DLL_HEADER void Ng_AddBoundarySeg_2D (Ng_Mesh * mesh, int pi1, int pi2)
   {
      Mesh * m = (Mesh*)mesh;

      Segment seg;
      seg[0] = pi1;
      seg[1] = pi2;
      m->AddSegment (seg);
   }


   DLL_HEADER int Ng_GetNP_2D (Ng_Mesh * mesh)
   {
      Mesh * m = (Mesh*)mesh;
      return m->GetNP();
   }

   DLL_HEADER int Ng_GetNE_2D (Ng_Mesh * mesh)
   {
      Mesh * m = (Mesh*)mesh;
      return m->GetNSE();
   }

   DLL_HEADER int Ng_GetNSeg_2D (Ng_Mesh * mesh)
   {
      Mesh * m = (Mesh*)mesh;
      return m->GetNSeg();
   }

   DLL_HEADER void Ng_GetPoint_2D (Ng_Mesh * mesh, int num, double * x)
   {
      Mesh * m = (Mesh*)mesh;

      Point<3> & p = m->Point(num);
      x[0] = p(0);
      x[1] = p(1);
   }

   DLL_HEADER void Ng_GetElement_2D (Ng_Mesh * mesh, int num, int * pi, int * matnum)
   {
      const Element2d & el = ((Mesh*)mesh)->SurfaceElement(num);
      for (int i = 1; i <= 3; i++)
         pi[i-1] = el.PNum(i);
      if (matnum)
         *matnum = el.GetIndex();
   }


   DLL_HEADER void Ng_GetSegment_2D (Ng_Mesh * mesh, int num, int * pi, int * matnum)
   {
      const Segment & seg = ((Mesh*)mesh)->LineSegment(num);
      pi[0] = seg[0];
      pi[1] = seg[1];

      if (matnum)
         *matnum = seg.edgenr;
   }




   DLL_HEADER Ng_Geometry_2D * Ng_LoadGeometry_2D (const char * filename)
   {
      SplineGeometry2d * geom = new SplineGeometry2d();
      geom -> Load (filename);
      return (Ng_Geometry_2D *)geom;
   }



   
   DLL_HEADER Ng_Result Ng_GenerateMesh_2D (Ng_Geometry_2D * geom,
                                            Ng_Mesh ** mesh,
                                            Ng_Meshing_Parameters * mp)
   {
      // use global variable mparam
      //  MeshingParameters mparam;  
      mparam.maxh = mp->maxh;
      mparam.meshsizefilename = mp->meshsize_filename;
      mparam.quad = mp->quad_dominated;

      Mesh * m;
      MeshFromSpline2D (*(SplineGeometry2d*)geom, m, mparam);

      cout << m->GetNSE() << " elements, " << m->GetNP() << " points" << endl;

      *mesh = (Ng_Mesh*)m;
      return NG_OK;
   }



   
   DLL_HEADER void Ng_HP_Refinement (Ng_Geometry_2D * geom,
      Ng_Mesh * mesh,
      int levels)
   {
      Refinement2d ref(*(SplineGeometry2d*)geom);
      HPRefinement (*(Mesh*)mesh, &ref, levels);
   }




   DLL_HEADER void Ng_HP_Refinement (Ng_Geometry_2D * geom,
      Ng_Mesh * mesh,
      int levels, double parameter)
   {
      Refinement2d ref(*(SplineGeometry2d*)geom);
      HPRefinement (*(Mesh*)mesh, &ref, levels, parameter);
   }








   Array<STLReadTriangle> readtrias; //only before initstlgeometry
   Array<Point<3> > readedges; //only before init stlgeometry

   // loads geometry from STL file
   DLL_HEADER Ng_STL_Geometry * Ng_STL_LoadGeometry (const char * filename, int binary)
   {
      int i;
      STLGeometry geom;
      STLGeometry* geo;
      ifstream ist(filename);

      if (binary)
      {
         geo = geom.LoadBinary(ist);
      }
      else
      {
         geo = geom.Load(ist);
      }

      readtrias.SetSize(0);
      readedges.SetSize(0);

      Point3d p;
      Vec3d normal;
      double p1[3];
      double p2[3];
      double p3[3];
      double n[3];

      Ng_STL_Geometry * geo2 = Ng_STL_NewGeometry();

      for (i = 1; i <= geo->GetNT(); i++)
      {
         const STLTriangle& t = geo->GetTriangle(i);
         p = geo->GetPoint(t.PNum(1));
         p1[0] = p.X(); p1[1] = p.Y(); p1[2] = p.Z(); 
         p = geo->GetPoint(t.PNum(2));
         p2[0] = p.X(); p2[1] = p.Y(); p2[2] = p.Z(); 
         p = geo->GetPoint(t.PNum(3));
         p3[0] = p.X(); p3[1] = p.Y(); p3[2] = p.Z();
         normal = t.Normal();
         n[0] = normal.X(); n[1] = normal.Y(); n[2] = normal.Z();

         Ng_STL_AddTriangle(geo2, p1, p2, p3, n);
      }

      return geo2;
   }



   // generate new STL Geometry
   DLL_HEADER Ng_STL_Geometry * Ng_STL_NewGeometry ()
   {
      return (Ng_STL_Geometry*)(void*)new STLGeometry;
   } 

   // generate new CSG Geometry
   DLL_HEADER Ng_CSG_Geometry * Ng_CSG_LoadGeometry(const char* filename)
   {
	return (Ng_CSG_Geometry*)(void*)new CSGeometry(filename);
   }


   // after adding triangles (and edges) initialize
   DLL_HEADER Ng_Result Ng_STL_InitSTLGeometry (Ng_STL_Geometry * geom)
   {
      STLGeometry* geo = (STLGeometry*)geom;
      geo->InitSTLGeometry(readtrias);
      readtrias.SetSize(0);

      if (readedges.Size() != 0)
      {
         /*
         for (int i = 1; i <= readedges.Size(); i+=2)
         {
         cout << "e(" << readedges.Get(i) << "," << readedges.Get(i+1) << ")" << endl;
         }
         */
         geo->AddEdges(readedges);
      }

      if (geo->GetStatus() == STLTopology::STL_GOOD || geo->GetStatus() == STLTopology::STL_WARNING) return NG_OK;
      return NG_SURFACE_INPUT_ERROR;
   }



   // automatically generates edges:
   DLL_HEADER Ng_Result Ng_STL_MakeEdges (Ng_STL_Geometry * geom,
                                          Ng_Mesh* mesh,
                                          Ng_Meshing_Parameters * mp)
   {
      STLGeometry* stlgeometry = (STLGeometry*)geom;
      Mesh* me = (Mesh*)mesh;

      // Philippose - 27/07/2009
      // Do not locally re-define "mparam" here... "mparam" is a global 
      // object 
      //MeshingParameters mparam;

      mparam.maxh = mp->maxh;
      mparam.meshsizefilename = mp->meshsize_filename;

      me -> SetGlobalH (mparam.maxh);
      me -> SetLocalH (stlgeometry->GetBoundingBox().PMin() - Vec3d(10, 10, 10),
                       stlgeometry->GetBoundingBox().PMax() + Vec3d(10, 10, 10),
                       0.3);

      me -> LoadLocalMeshSize (mp->meshsize_filename);
      /*
      if (mp->meshsize_filename)
      {
      ifstream infile (mp->meshsize_filename);
      if (!infile.good()) return NG_FILE_NOT_FOUND;
      me -> LoadLocalMeshSize (infile);
      }
      */

      STLMeshing (*stlgeometry, *me);

      stlgeometry->edgesfound = 1;
      stlgeometry->surfacemeshed = 0;
      stlgeometry->surfaceoptimized = 0;
      stlgeometry->volumemeshed = 0;

      return NG_OK;
   }


   // generates mesh, empty mesh be already created.
   DLL_HEADER Ng_Result Ng_STL_GenerateSurfaceMesh (Ng_STL_Geometry * geom,
                                                    Ng_Mesh* mesh,
                                                    Ng_Meshing_Parameters * mp)
   {
      STLGeometry* stlgeometry = (STLGeometry*)geom;
      Mesh* me = (Mesh*)mesh;

      // Philippose - 27/07/2009
      // Do not locally re-define "mparam" here... "mparam" is a global 
      // object
      //MeshingParameters mparam;

      mparam.maxh = mp->maxh;
      mparam.meshsizefilename = mp->meshsize_filename;

      /*
      me -> SetGlobalH (mparam.maxh);
      me -> SetLocalH (stlgeometry->GetBoundingBox().PMin() - Vec3d(10, 10, 10),
      stlgeometry->GetBoundingBox().PMax() + Vec3d(10, 10, 10),
      0.3);
      */
      /*
      STLMeshing (*stlgeometry, *me);

      stlgeometry->edgesfound = 1;
      stlgeometry->surfacemeshed = 0;
      stlgeometry->surfaceoptimized = 0;
      stlgeometry->volumemeshed = 0;
      */  
      int retval = STLSurfaceMeshing (*stlgeometry, *me);
      if (retval == MESHING3_OK)
      {
         (*mycout) << "Success !!!!" << endl;
         stlgeometry->surfacemeshed = 1;
         stlgeometry->surfaceoptimized = 0;
         stlgeometry->volumemeshed = 0;
      } 
      else if (retval == MESHING3_OUTERSTEPSEXCEEDED)
      {
         (*mycout) << "ERROR: Give up because of too many trials. Meshing aborted!" << endl;
      }
      else if (retval == MESHING3_TERMINATE)
      {
         (*mycout) << "Meshing Stopped!" << endl;
      }
      else
      {
         (*mycout) << "ERROR: Surface meshing not successful. Meshing aborted!" << endl;
      }


      STLSurfaceOptimization (*stlgeometry, *me, mparam);

      return NG_OK;
   }


   // fills STL Geometry
   // positive orientation
   // normal vector may be null-pointer
   DLL_HEADER void Ng_STL_AddTriangle (Ng_STL_Geometry * geom, 
                                       double * p1, double * p2, double * p3, 
                                       double * nv)
   {
      Point<3> apts[3];
      apts[0] = Point<3>(p1[0],p1[1],p1[2]);
      apts[1] = Point<3>(p2[0],p2[1],p2[2]);
      apts[2] = Point<3>(p3[0],p3[1],p3[2]);

      Vec<3> n;
      if (!nv)
         n = Cross (apts[0]-apts[1], apts[0]-apts[2]);
      else
         n = Vec<3>(nv[0],nv[1],nv[2]);

      readtrias.Append(STLReadTriangle(apts,n));
   }

   // add (optional) edges:
   DLL_HEADER void Ng_STL_AddEdge (Ng_STL_Geometry * geom, 
      double * p1, double * p2)
   {
      readedges.Append(Point3d(p1[0],p1[1],p1[2]));
      readedges.Append(Point3d(p2[0],p2[1],p2[2]));
   }





#ifdef OCCGEOMETRY
   // --------------------- OCC Geometry / Meshing Utility Functions -------------------

   // Create new OCC Geometry Object
   DLL_HEADER Ng_OCC_Geometry * Ng_OCC_NewGeometry ()
   {
      return (Ng_OCC_Geometry*)(void*)new OCCGeometry;
   } 


   
   // Delete the OCC Geometry Object
   DLL_HEADER Ng_Result Ng_OCC_DeleteGeometry(Ng_OCC_Geometry * geom)
   {
      if (geom != NULL)
      {
         delete (OCCGeometry*)geom;
         geom = NULL;
         return NG_OK;
      }
      
      return NG_ERROR;
   }


   
   // Loads geometry from STEP File
   DLL_HEADER Ng_OCC_Geometry * Ng_OCC_Load_STEP (const char * filename)
   {
      // Call the STEP File Load function. Note.. the geometry class 
      // is created and instantiated within the load function
      OCCGeometry * occgeo = LoadOCC_STEP(filename);

      return ((Ng_OCC_Geometry *)occgeo);
   }


   
   // Loads geometry from IGES File
   DLL_HEADER Ng_OCC_Geometry * Ng_OCC_Load_IGES (const char * filename)
   {
      // Call the IGES File Load function. Note.. the geometry class 
      // is created and instantiated within the load function
      OCCGeometry * occgeo = LoadOCC_IGES(filename);

      return ((Ng_OCC_Geometry *)occgeo);
   }


   
   // Loads geometry from BREP File
   DLL_HEADER Ng_OCC_Geometry * Ng_OCC_Load_BREP (const char * filename)
   {
      // Call the BREP File Load function. Note.. the geometry class 
      // is created and instantiated within the load function
      OCCGeometry * occgeo = LoadOCC_BREP(filename);

      return ((Ng_OCC_Geometry *)occgeo);
   }



   // Locally limit the size of the mesh to be generated at various points 
   // based on the topology of the geometry
   DLL_HEADER Ng_Result Ng_OCC_SetLocalMeshSize (Ng_OCC_Geometry * geom,
                                                 Ng_Mesh * mesh,
                                                 Ng_Meshing_Parameters * mp)
   {
      OCCGeometry * occgeom = (OCCGeometry*)geom;
      Mesh * me = (Mesh*)mesh;

      me->geomtype = Mesh::GEOM_OCC;

      mparam.uselocalh = mp->uselocalh;

      mparam.maxh = mp->maxh;
      mparam.minh = mp->minh;

      mparam.segmentsperedge = mp->elementsperedge;
      mparam.curvaturesafety = mp->elementspercurve;

      mparam.grading = mp->grading;
      mparam.meshsizefilename = mp->meshsize_filename;
      
      occparam.resthcloseedgeenable = mp->closeedgeenable;
      occparam.resthcloseedgefac = mp->closeedgefact;

      // Delete the mesh structures in order to start with a clean 
      // slate
      me->DeleteMesh();

      OCCSetLocalMeshSize(*occgeom, *me);

      return(NG_OK);
   }


   
   // Mesh the edges and add Face descriptors to prepare for surface meshing
   DLL_HEADER Ng_Result Ng_OCC_GenerateEdgeMesh (Ng_OCC_Geometry * geom,
                                                 Ng_Mesh * mesh,
                                                 Ng_Meshing_Parameters * mp)
   {
      OCCGeometry * occgeom = (OCCGeometry*)geom;
      Mesh * me = (Mesh*)mesh;

      mparam.uselocalh = mp->uselocalh;

      OCCFindEdges(*occgeom, *me);

      if((me->GetNP()) && (me->GetNFD()))
      {
         return NG_OK;
      }
      else
      {
         return NG_ERROR;
      }
   }



   
   // Mesh the edges and add Face descriptors to prepare for surface meshing
   DLL_HEADER Ng_Result Ng_OCC_GenerateSurfaceMesh (Ng_OCC_Geometry * geom,
                                                    Ng_Mesh * mesh,
                                                    Ng_Meshing_Parameters * mp)
   {
      int numpoints = 0;

      OCCGeometry * occgeom = (OCCGeometry*)geom;
      Mesh * me = (Mesh*)mesh;

      mparam.uselocalh = mp->uselocalh;

      // Only go into surface meshing if the face descriptors have already been added
      if(!me->GetNFD())
         return NG_ERROR;

      numpoints = me->GetNP();

      // Initially set up only for surface meshing without any optimisation
      int perfstepsend = MESHCONST_MESHSURFACE;

      // Check and if required, enable surface mesh optimisation step
      if(mp->optsurfmeshenable)
      {
         perfstepsend = MESHCONST_OPTSURFACE;
      }

      OCCMeshSurface(*occgeom, *me, perfstepsend);

      me->CalcSurfacesOfNode();
      
      if(me->GetNP() <= numpoints)
         return NG_ERROR;

      if(me->GetNSE() <= 0)
         return NG_ERROR;

      return NG_OK;
   }


   // Extract the face map from the OCC geometry
   // The face map basically gives an index to each face in the geometry, 
   // which can be used to access a specific face
   DLL_HEADER Ng_Result Ng_OCC_GetFMap(Ng_OCC_Geometry * geom, 
                                       Ng_OCC_TopTools_IndexedMapOfShape * FMap)
   {
      OCCGeometry* occgeom = (OCCGeometry*)geom;
      TopTools_IndexedMapOfShape *occfmap = (TopTools_IndexedMapOfShape *)FMap;

      // Copy the face map from the geometry to the given variable
      occfmap->Assign(occgeom->fmap);

      if(occfmap->Extent())
      {
         return NG_OK;
      }
      else
      {
         return NG_ERROR;
      }
   }

   // ------------------ End - OCC Geometry / Meshing Utility Functions ----------------
#endif






   DLL_HEADER Ng_Meshing_Parameters :: Ng_Meshing_Parameters()
   {
      uselocalh = 1;

      maxh = 1000;
      minh = 0.0;

      fineness = 0.5;
      grading = 0.3;

      elementsperedge = 2.0;
      elementspercurve = 2.0;

      closeedgeenable = 0;
      closeedgefact = 2.0;

      secondorder = 0;
      quad_dominated = 0;

      meshsize_filename = 0;

      optsurfmeshenable = 1;
      optvolmeshenable = 1;

      optsteps_2d = 3;
      optsteps_3d = 3;
   }



   DLL_HEADER void Ng_Meshing_Parameters :: Reset_Parameters()
   {
      uselocalh = 1;

      maxh = 1000;
      minh = 0;

      fineness = 0.5;
      grading = 0.3;

      elementsperedge = 2.0;
      elementspercurve = 2.0;

      closeedgeenable = 0;
      closeedgefact = 2.0;

      secondorder = 0;
      quad_dominated = 0;

      meshsize_filename = 0;

      optsurfmeshenable = 1;
      optvolmeshenable = 1;

      optsteps_2d = 3;
      optsteps_3d = 3;
   }












  DLL_HEADER void Ng_Uniform_Refinement (Ng_Mesh * mesh)
  {
    Refinement ref;
    ref.Refine ( * (Mesh*) mesh );
  }

  DLL_HEADER void Ng_2D_Uniform_Refinement (Ng_Geometry_2D * geom,
					 Ng_Mesh * mesh)
  {
    ( (SplineGeometry2d*)geom ) -> GetRefinement().Refine ( * (Mesh*) mesh );
  }

  DLL_HEADER void Ng_STL_Uniform_Refinement (Ng_STL_Geometry * geom,
					     Ng_Mesh * mesh)
  {
    ( (STLGeometry*)geom ) -> GetRefinement().Refine ( * (Mesh*) mesh );
  }

  DLL_HEADER void Ng_CSG_Uniform_Refinement (Ng_CSG_Geometry * geom,
					     Ng_Mesh * mesh)
  {
    ( (CSGeometry*)geom ) -> GetRefinement().Refine ( * (Mesh*) mesh );
  }


#ifdef OCCGEOMETRY
  DLL_HEADER void Ng_OCC_Uniform_Refinement (Ng_OCC_Geometry * geom,
					     Ng_Mesh * mesh)
  {
    ( (OCCGeometry*)geom ) -> GetRefinement().Refine ( * (Mesh*) mesh );
  }
#endif



} // End of namespace nglib








// compatibility functions:

namespace netgen 
{

   char geomfilename[255];

   DLL_HEADER void MyError (const char * ch)
   {
      cerr << ch;
   }

   //Destination for messages, errors, ...
   DLL_HEADER void Ng_PrintDest(const char * s)
   {
      (*mycout) << s << flush;
   }

   DLL_HEADER double GetTime ()
   {
      return 0;
   }

   void ResetTime ()
   {
      ;
   }

   void MyBeep (int i)
   {
      ;
   }

   void Render()
   {
      ; 
   }

} // End of namespace netgen

