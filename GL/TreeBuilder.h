/*
 *  TreeBuilder.h
 *
 *  Created by Yusuf Yilmaz on 5/23/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#define SPLITTYPE 3 //1: OctreeType splitting, 2: mean splitting(not implemented yet) 3:median splitting
#define BSPTYPE 1 // 1: Game BSP Type, 2: iterative PCA based
#define REPSIZE 1000 // Repetition for iterative PCA

#define smoothOut 1 // Smooth Out the bisection areas
#define smoothFact 8 // Factor to smooth
#define ROUNDFACT 3

#ifndef _TREEBUILDER_H
#define _TREEBUILDER_H

#include <mpi.h>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <cstring>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <algorithm>
#include <time.h>
#include <fstream>

#include "IElement.h"
#include "Mesh.h"
#include "BSPTreeNode.h"
#include "parmetisbin.h"
#include "MeshMig.h"

namespace nglib {
#include <nglib.h>
}

enum CLASSIFY{
    CLASSIFY_LEFT,
    CLASSIFY_RIGHT,
    CLASSIFY_ONPLANE
};

class TreeBuilder {

public:

    int tot;

	TreeBuilder(void)
	{
	    tot = 0;
		totalnode = 0;
		totalleafnode = 0;
	}

	~TreeBuilder(void)	// destructor
	{
		unsigned int i;

		// delete elements. We delete elements here and not in the mesh destructor.
		std::list<IElement*>::iterator elementListIt;
		for (elementListIt = elementList.begin(); elementListIt != elementList.end(); elementListIt++) {
			delete(*elementListIt);
		}

        for (elementListIt = surfaceList.begin(); elementListIt != elementList.end(); elementListIt++) {
            delete(*elementListIt);
        }

		// delete meshes
		for (i=0; i<meshList.size(); i++)
			delete(meshList[i]);
		meshList.clear();

	}

	void setStringSend(int parsize)
	{
		sindex = 0;
		stringtosend = new std::string[parsize];
	}



	// Build BSP-Tree Start
	int buildBSPTree(unsigned int recdepth) //OK
	{
		rootNode = new BSPTreeNode;
		totalnode++;
		rootNode->level = 0;

		Vector3D l(9999999, 9999999, 9999999);	// set max numbers.
		Vector3D u(-9999999,-9999999,-9999999);
		rootNode->bb = AABB(l, u);

		std::list<IElement*>::const_iterator element;
		// Fill the root node.

		for(element = elementList.begin(); element != elementList.end(); element++)
		{
			rootNode->elements.push_back((*element));
		}
		// find the bounding box of the whole Tree.

		if(SPLITTYPE != 5)
		{
            AABB temp;
            std::vector<Mesh*>::const_iterator meshs;
            for(meshs = meshList.begin(); meshs != meshList.end(); meshs++)
            {
                temp = (*meshs)->getBB();
                if(temp.lowerCorner.x < rootNode->bb.lowerCorner.x)
                    rootNode->bb.lowerCorner.x = temp.lowerCorner.x;
                if(temp.lowerCorner.y < rootNode->bb.lowerCorner.y)
                    rootNode->bb.lowerCorner.y = temp.lowerCorner.y;
                if(temp.lowerCorner.z < rootNode->bb.lowerCorner.z)
                    rootNode->bb.lowerCorner.z = temp.lowerCorner.z;

                if(temp.upperCorner.x > rootNode->bb.upperCorner.x)
                    rootNode->bb.upperCorner.x = temp.upperCorner.x;
                if(temp.upperCorner.y > rootNode->bb.upperCorner.y)
                    rootNode->bb.upperCorner.y = temp.upperCorner.y;
                if(temp.upperCorner.z > rootNode->bb.upperCorner.z)
                    rootNode->bb.upperCorner.z = temp.upperCorner.z;
            }

            int biggest = 0;

            if((rootNode->bb.upperCorner.x - rootNode->bb.lowerCorner.x) > (rootNode->bb.upperCorner.y - rootNode->bb.lowerCorner.y))
            {
                biggest = 0;

            }
            else
            {
                biggest = 1;
            }

            if((rootNode->bb.upperCorner[biggest] - rootNode->bb.lowerCorner[biggest]) > (rootNode->bb.upperCorner.z - rootNode->bb.lowerCorner.z))
            {
                biggest = biggest;
            }
            else
            {
                biggest = 2;
            }



            rootNode->axis = biggest;
		}


		if(SPLITTYPE == 1)  // Octree.
			rootNode->splitCoordinate = (rootNode->bb.upperCorner[rootNode->axis] + rootNode->bb.lowerCorner[rootNode->axis])/2;
		else if(SPLITTYPE == 2)
		{
			// here make mean splitting.
		}
		else if(SPLITTYPE == 3) //KD-Tree
		{
			std::vector<double> centers;
			std::list<IElement*>::const_iterator element2;


			for(element2 = rootNode->elements.begin(); element2 != rootNode->elements.end(); element2++)
			{
				centers.push_back((*element2)->getCentroid()[rootNode->axis]);
			}

			std::sort(centers.begin(),centers.end());
			if(rootNode->elements.size() % 2 != 0)
				rootNode->splitCoordinate = centers.at(int(rootNode->elements.size()/2));
			else
				rootNode->splitCoordinate = (centers.at(int(rootNode->elements.size()/2)) + centers.at(int(rootNode->elements.size()/2 - 1/2))) / 2;
		}
		else if(SPLITTYPE == 5) //  BSP-Tree
		{
		    if(BSPTYPE == 1) // Random Game style
		    {
                // BSP Splitting Here.
                std::vector<IElement*> myVect;
                std::list<IElement*>::const_iterator element3;

                for(element3 = rootNode->elements.begin(); element3 != rootNode->elements.end(); element3++)
                {
                    myVect.push_back(*element3);
                }

                int mySize = myVect.size();
                int score = 0;
                int minScore = mySize / 10;

                Vector3D pt1;
                Vector3D norm1;
                double d1;
                MeshVertex* vt1;
                MeshVertex* vt2;
                MeshVertex* vt3;

                Vector3D vt1t2;
                Vector3D vt1t3;

                do
                {
                    score = 0;

                    int selected = rand() % mySize;

                    IElement* myElem = myVect[selected];

                    int ret2 = 0;

                    vt1 = ((MeshTetrahedron*)myElem)->getVertex(0);
                    vt2 = ((MeshTetrahedron*)myElem)->getVertex(1);
                    vt3 = ((MeshTetrahedron*)myElem)->getVertex(2);

                    vt1t2 = *(vt2->getPosition()) - *(vt1->getPosition());
                    vt1t3 = *(vt3->getPosition()) - *(vt1->getPosition());

                    pt1 = (*(vt1->getPosition()) + *(vt2->getPosition()) + *(vt3->getPosition()))/3;
                    Vector3DCross(&norm1, &vt1t2, &vt1t3);

                    d1 = (-1) * (pt1[0]*norm1[0] + pt1[1]*norm1[1] + pt1[2]*norm1[2]);

                    std::list<IElement*>::const_iterator element4;

                    for(element4 = rootNode->elements.begin(); element4 != rootNode->elements.end(); element4++)
                    {
                        ret2 = classifyPointPlane((*element4)->getCentroid(), norm1, d1);

                        if(ret2 == CLASSIFY_LEFT)
                        {
                            score++;
                        }
                        else
                        {
                            score--;
                        }
                    }

                    std::cout << "SCORE: " << score << " MINSCORE: " << minScore << " SIZE: " << mySize << std::endl;
                }
                while(abs(score) > minScore);

                rootNode->splitPoint = pt1;
                rootNode->planeNormal = norm1;
                rootNode->d = d1;
		    }
		    else if(BSPTYPE == 2)   // Iterative PCA
		    {
                int sz = rootNode->elements.size();

                double* points = new double[sz * 3];
                Vector3D elem3;
                double totx = 0, toty = 0, totz = 0;
                double* orig = new double[3];
                double* dir = new double[3];

                int k = 0;
                std::list<IElement*>::const_iterator elementt;
                for(elementt = rootNode->elements.begin(); elementt != rootNode->elements.end(); elementt++)
                {
                    elem3 = (*elementt)->getCentroid();

                    points[3*k] = elem3[0];
                    totx += points[3*k];
                    points[3*k+1] = elem3[1];
                    toty += points[3*k+1];
                    points[3*k+2] = elem3[2];
                    totz += points[3*k+2];
                    k++;
                }

                totx /= sz;
                toty /= sz;
                totz /= sz;

                orig[0] = totx;
                orig[1] = toty;
                orig[2] = totz;

                // calculate inertia axis here and put it into orig[3] and dir[3]

                double px,py,pz = 0;
                double tx,ty,tz = 0;

                // std::cout << px << " " << py << " " << pz << std::endl;

                for(int i=0; i<sz; i++)
                {
                    points[3*i] -= orig[0];
                    points[3*i+1] -= orig[1];
                    points[3*i+2] -= orig[2];
                }


                double xp = 0;

                px = rand()%10000;
                py = rand()%10000;
                pz = rand()%10000;

                double tot = sqrt(px*px + py*py + pz*pz);

                px /= tot;
                py /= tot;
                pz /= tot;

                for(int i=0; i<REPSIZE; i++)
                {
                    tx = 0;
                    ty = 0;
                    tz = 0;

                    for(int j=0; j<sz; j++)
                    {
                        xp = px*points[3*i] + py*points[3*i+1] + pz*points[3*i+2];

                        tx += xp * points[3*i];
                        ty += xp * points[3*i+1];
                        tz += xp * points[3*i+2];
                    }

                    tot = sqrt(tx*tx + ty*ty + tz*tz);
                    px = tx/tot;
                    py = ty/tot;
                    pz = tz/tot;
                }

                dir[0] = px;
                dir[1] = py;
                dir[2] = pz;

                orig[0] = round(orig[0], ROUNDFACT);
                orig[1] = round(orig[1], ROUNDFACT);
                orig[2] = round(orig[2], ROUNDFACT);
                dir[0] = round(dir[0], ROUNDFACT);
                dir[1] = round(dir[1], ROUNDFACT);
                dir[2] = round(dir[2], ROUNDFACT);

                rootNode->splitPoint = Vector3D(orig[0], orig[1], orig[2]);
                rootNode->planeNormal = Vector3D(dir[0], dir[1], dir[2]);

                rootNode->d = (-1) * (dir[0]*orig[0] + dir[1]*orig[1] + dir[2]*orig[2]);

                delete [] points;
		    }
		}

        if(SPLITTYPE != 5)
        {
            if(smoothOut)
            {
                /// update split coordinate. 10'a bol en fazla 10da 1'i kadar cozunurluk olsun.
                double width2 = rootNode->bb.upperCorner[rootNode->axis] - rootNode->bb.lowerCorner[rootNode->axis];
                double temp3 = rootNode->splitCoordinate - rootNode->bb.lowerCorner[rootNode->axis];
                double pos2 = temp3/width2;
                int tmp2 = (int)(pos2 * smoothFact + 0.5);
                rootNode->splitCoordinate = rootNode->bb.lowerCorner[rootNode->axis] + (tmp2 * width2 / smoothFact);
            }
            rootNode->splitCoordinate = round(rootNode->splitCoordinate, ROUNDFACT);

            // Calculate the vectors for the split plane.
            rootNode->ll = Vector3D(rootNode->splitCoordinate, rootNode->bb.lowerCorner[1], rootNode->bb.lowerCorner[2]);
            rootNode->lr = Vector3D(rootNode->splitCoordinate, rootNode->bb.lowerCorner[1], rootNode->bb.upperCorner[2]);
            rootNode->ur = Vector3D(rootNode->splitCoordinate, rootNode->bb.upperCorner[1], rootNode->bb.upperCorner[2]);
            rootNode->ul = Vector3D(rootNode->splitCoordinate, rootNode->bb.upperCorner[1], rootNode->bb.lowerCorner[2]);

            switch(rootNode->axis)
            {
                case 0:
                    rootNode->planeNormal = Vector3D(1,0,0);
                    break;
                case 1:
                    rootNode->planeNormal = Vector3D(0,1,0);
                    break;
                case 2:
                    rootNode->planeNormal = Vector3D(0,0,1);
                    break;
            }
        }
		rootNode->planeclr = Vector3D((rand() % 100/100.0), (rand() % 100/100.0), (rand() % 100/100.0));

        // Put the node name here.
        rootNode->partname = "main_";
        rootNode->planename = "plane_";

		rootNode->leftChildNode = NULL;
		rootNode->rightChildNode = NULL;
		// call the recursive function.
		buildBSPTreeRec(rootNode, recdepth);

		return 0;
	}

	// Build Tree kisminin recursive fonksiyonu.
	void buildBSPTreeRec(BSPTreeNode *parentNode,int recdepth) //OK
	{
		// depth levela eristiysek yada bir node da minimum istedigimiz eleman sayisina eristiysek duruyoruz.
        int biggest = 0;
        if( parentNode->level == recdepth)
		{
            int rank = 0;

            for(int i=5; i<(recdepth+5) ; i++)
            {
                if(parentNode->partname.at(i) == 'l')
                {
                    rank = rank*2;
                }
                else if(parentNode->partname.at(i) == 'r')
                {
                    rank = rank*2;
                    rank += 1;
                }
            }
            parentNode->partId = rank;
            std::cout << parentNode->partId << std::endl;

			totalleafnode++;
			return;
		}

        if(SPLITTYPE !=5)
        {
            switch(parentNode->axis)
            {
                case 0:
                    parentNode->planeNormal = Vector3D(1,0,0);
                    break;
                case 1:
                    parentNode->planeNormal = Vector3D(0,1,0);
                    break;
                case 2:
                    parentNode->planeNormal = Vector3D(0,0,1);
                    break;
            }
        }


		if(parentNode->level  == 1 || parentNode->level == 2)	// to see progress in the console.
			std::cout << "Level " << parentNode->level;

		// sag node yaratiliyor.
		parentNode->leftChildNode = new BSPTreeNode;
		totalnode++;
		parentNode->leftChildNode->leftChildNode = NULL;
		parentNode->leftChildNode->rightChildNode = NULL;

		if(SPLITTYPE != 5)
		{
            // parentNode->leftChildNode->axis = (parentNode->axis+1)%3;	// alternate axis
            parentNode->leftChildNode->bb = parentNode->bb;				// set bb of parent node
            parentNode->leftChildNode->bb.upperCorner[parentNode->axis] = parentNode->splitCoordinate;	// only change uppercorner.
            biggest = 0;

            if((parentNode->leftChildNode->bb.upperCorner.x - parentNode->leftChildNode->bb.lowerCorner.x) > (parentNode->leftChildNode->bb.upperCorner.y - parentNode->leftChildNode->bb.lowerCorner.y))
            {
                biggest = 0;

            }
            else
            {
                biggest = 1;
            }

            if((parentNode->leftChildNode->bb.upperCorner[biggest] - parentNode->leftChildNode->bb.lowerCorner[biggest]) > (parentNode->leftChildNode->bb.upperCorner.z - parentNode->leftChildNode->bb.lowerCorner.z))
            {
                biggest = biggest;
            }
            else
            {
                biggest = 2;
            }
            parentNode->leftChildNode->axis = biggest;
            biggest = 0;
		}

		parentNode->leftChildNode->level = parentNode->level + 1;	// set level
        parentNode->leftChildNode->partname = parentNode->partname + "l";
        parentNode->leftChildNode->planename = parentNode->planename + "l";

		// sol node yaratiliyor.
		parentNode->rightChildNode = new BSPTreeNode;
		totalnode++;
		parentNode->rightChildNode->leftChildNode = NULL;
		parentNode->rightChildNode->rightChildNode = NULL;

		if(SPLITTYPE != 5)
		{
            // parentNode->rightChildNode->axis = (parentNode->axis+1)%3;	// alternate axis.
            parentNode->rightChildNode->bb = parentNode->bb;			// set bb of parent node
            parentNode->rightChildNode->bb.lowerCorner[parentNode->axis] = parentNode->splitCoordinate;	// only change lowercorner.
            if((parentNode->rightChildNode->bb.upperCorner.x - parentNode->rightChildNode->bb.lowerCorner.x) > (parentNode->rightChildNode->bb.upperCorner.y - parentNode->rightChildNode->bb.lowerCorner.y))
            {
                biggest = 0;
            }
            else
            {
                biggest = 1;
            }

            if((parentNode->rightChildNode->bb.upperCorner[biggest] - parentNode->rightChildNode->bb.lowerCorner[biggest]) > (parentNode->rightChildNode->bb.upperCorner.z - parentNode->rightChildNode->bb.lowerCorner.z))
            {
                biggest = biggest;
            }
            else
            {
                biggest = 2;
            }
            parentNode->rightChildNode->axis = biggest;
            biggest = 0;
		}
		parentNode->rightChildNode->level = parentNode->level + 1;	// set level
        parentNode->rightChildNode->partname = parentNode->partname + "r";
        parentNode->rightChildNode->planename = parentNode->planename + "r";

		fillChildNodes( parentNode );

		// sag ve sol nodelar icin split coordinate hesabi yap.
		parentNode->rightChildNode->splitCoordinate = computeSplitPlanePos(parentNode,false);

		if(SPLITTYPE != 5)
		{
            if(smoothOut)
            {
                /// update split coordinate. 10'a bol en fazla 10da 1'i kadar cozunurluk olsun.
                double width2 = parentNode->rightChildNode->bb.upperCorner[parentNode->rightChildNode->axis] - parentNode->rightChildNode->bb.lowerCorner[parentNode->rightChildNode->axis];
                double temp3 = parentNode->rightChildNode->splitCoordinate - parentNode->rightChildNode->bb.lowerCorner[parentNode->rightChildNode->axis];
                double pos2 = temp3/width2;
                int tmp2 = (int)(pos2 * smoothFact + 0.5);
                parentNode->rightChildNode->splitCoordinate = parentNode->rightChildNode->bb.lowerCorner[parentNode->rightChildNode->axis] + (tmp2 * width2 / smoothFact);
            }
            parentNode->rightChildNode->splitCoordinate = round(parentNode->rightChildNode->splitCoordinate, ROUNDFACT);

            switch (parentNode->rightChildNode->axis) {
                case 0:
                    parentNode->rightChildNode->ll = Vector3D(parentNode->rightChildNode->splitCoordinate, parentNode->rightChildNode->bb.lowerCorner[1], parentNode->rightChildNode->bb.lowerCorner[2]);
                    parentNode->rightChildNode->lr = Vector3D(parentNode->rightChildNode->splitCoordinate, parentNode->rightChildNode->bb.lowerCorner[1], parentNode->rightChildNode->bb.upperCorner[2]);
                    parentNode->rightChildNode->ur = Vector3D(parentNode->rightChildNode->splitCoordinate, parentNode->rightChildNode->bb.upperCorner[1], parentNode->rightChildNode->bb.upperCorner[2]);
                    parentNode->rightChildNode->ul = Vector3D(parentNode->rightChildNode->splitCoordinate, parentNode->rightChildNode->bb.upperCorner[1], parentNode->rightChildNode->bb.lowerCorner[2]);
                    break;

                case 1:
                    parentNode->rightChildNode->ll = Vector3D(parentNode->rightChildNode->bb.lowerCorner[0], parentNode->rightChildNode->splitCoordinate, parentNode->rightChildNode->bb.lowerCorner[2]);
                    parentNode->rightChildNode->lr = Vector3D(parentNode->rightChildNode->bb.lowerCorner[0], parentNode->rightChildNode->splitCoordinate, parentNode->rightChildNode->bb.upperCorner[2]);
                    parentNode->rightChildNode->ur = Vector3D(parentNode->rightChildNode->bb.upperCorner[0], parentNode->rightChildNode->splitCoordinate, parentNode->rightChildNode->bb.upperCorner[2]);
                    parentNode->rightChildNode->ul = Vector3D(parentNode->rightChildNode->bb.upperCorner[0], parentNode->rightChildNode->splitCoordinate, parentNode->rightChildNode->bb.lowerCorner[2]);
                    break;

                case 2:
                    parentNode->rightChildNode->ll = Vector3D(parentNode->rightChildNode->bb.lowerCorner[0], parentNode->rightChildNode->bb.lowerCorner[1], parentNode->rightChildNode->splitCoordinate);
                    parentNode->rightChildNode->lr = Vector3D(parentNode->rightChildNode->bb.lowerCorner[0], parentNode->rightChildNode->bb.upperCorner[1], parentNode->rightChildNode->splitCoordinate);
                    parentNode->rightChildNode->ur = Vector3D(parentNode->rightChildNode->bb.upperCorner[0], parentNode->rightChildNode->bb.upperCorner[1], parentNode->rightChildNode->splitCoordinate);
                    parentNode->rightChildNode->ul = Vector3D(parentNode->rightChildNode->bb.upperCorner[0], parentNode->rightChildNode->bb.lowerCorner[1], parentNode->rightChildNode->splitCoordinate);
                    break;

                default:
                    break;
            }
		}
		parentNode->rightChildNode->planeclr = Vector3D((rand() % 100/100.0), (rand() % 100/100.0), (rand() % 100/100.0));

		parentNode->leftChildNode->splitCoordinate = computeSplitPlanePos(parentNode,true);
		if(SPLITTYPE != 5)
		{
            if(smoothOut)
            {
                /// update split coordinate. 10'a bol en fazla 10da 1'i kadar cozunurluk olsun.
                double width2 = parentNode->leftChildNode->bb.upperCorner[parentNode->leftChildNode->axis] - parentNode->leftChildNode->bb.lowerCorner[parentNode->leftChildNode->axis];
                double temp3 = parentNode->leftChildNode->splitCoordinate - parentNode->leftChildNode->bb.lowerCorner[parentNode->leftChildNode->axis];
                double pos2 = temp3/width2;
                int tmp2 = (int)(pos2 * smoothFact + 0.5);
                parentNode->leftChildNode->splitCoordinate = parentNode->leftChildNode->bb.lowerCorner[parentNode->leftChildNode->axis] + (tmp2 * width2 / smoothFact);
            }
            parentNode->leftChildNode->splitCoordinate = round(parentNode->leftChildNode->splitCoordinate, ROUNDFACT);


            switch (parentNode->leftChildNode->axis) {
                case 0:
                    parentNode->leftChildNode->ll = Vector3D(parentNode->leftChildNode->splitCoordinate, parentNode->leftChildNode->bb.lowerCorner[1], parentNode->leftChildNode->bb.lowerCorner[2]);
                    parentNode->leftChildNode->lr = Vector3D(parentNode->leftChildNode->splitCoordinate, parentNode->leftChildNode->bb.lowerCorner[1], parentNode->leftChildNode->bb.upperCorner[2]);
                    parentNode->leftChildNode->ur = Vector3D(parentNode->leftChildNode->splitCoordinate, parentNode->leftChildNode->bb.upperCorner[1], parentNode->leftChildNode->bb.upperCorner[2]);
                    parentNode->leftChildNode->ul = Vector3D(parentNode->leftChildNode->splitCoordinate, parentNode->leftChildNode->bb.upperCorner[1], parentNode->leftChildNode->bb.lowerCorner[2]);
                    break;

                case 1:
                    parentNode->leftChildNode->ll = Vector3D(parentNode->leftChildNode->bb.lowerCorner[0], parentNode->leftChildNode->splitCoordinate, parentNode->leftChildNode->bb.lowerCorner[2]);
                    parentNode->leftChildNode->lr = Vector3D(parentNode->leftChildNode->bb.lowerCorner[0], parentNode->leftChildNode->splitCoordinate, parentNode->leftChildNode->bb.upperCorner[2]);
                    parentNode->leftChildNode->ur = Vector3D(parentNode->leftChildNode->bb.upperCorner[0], parentNode->leftChildNode->splitCoordinate, parentNode->leftChildNode->bb.upperCorner[2]);
                    parentNode->leftChildNode->ul = Vector3D(parentNode->leftChildNode->bb.upperCorner[0], parentNode->leftChildNode->splitCoordinate, parentNode->leftChildNode->bb.lowerCorner[2]);
                    break;

                case 2:
                    parentNode->leftChildNode->ll = Vector3D(parentNode->leftChildNode->bb.lowerCorner[0], parentNode->leftChildNode->bb.lowerCorner[1], parentNode->leftChildNode->splitCoordinate);
                    parentNode->leftChildNode->lr = Vector3D(parentNode->leftChildNode->bb.lowerCorner[0], parentNode->leftChildNode->bb.upperCorner[1], parentNode->leftChildNode->splitCoordinate);
                    parentNode->leftChildNode->ur = Vector3D(parentNode->leftChildNode->bb.upperCorner[0], parentNode->leftChildNode->bb.upperCorner[1], parentNode->leftChildNode->splitCoordinate);
                    parentNode->leftChildNode->ul = Vector3D(parentNode->leftChildNode->bb.upperCorner[0], parentNode->leftChildNode->bb.lowerCorner[1], parentNode->leftChildNode->splitCoordinate);
                    break;

                default:
                    break;
            }
		}
		parentNode->leftChildNode->planeclr = Vector3D((rand() % 100/100.0), (rand() % 100/100.0), (rand() % 100/100.0));

		buildBSPTreeRec(parentNode->leftChildNode, recdepth);	// sol tarafi yap.
		buildBSPTreeRec(parentNode->rightChildNode, recdepth); // sag tarafi yap.
	}

	// Split coordinatini belirleyen fonksiyon. Suanda sadece nodeu tam ortadan bolen bir yapi var. Istenirse bu fonksiyonun ici
	// degistirilerek median splitting veya degisik heuristic algoritmalar kullanilip denenebilinir.
	double computeSplitPlanePos(BSPTreeNode * parentNode, bool isLeft) //OK
	{
		if(SPLITTYPE == 1)
		{
			return (parentNode->bb.upperCorner[parentNode->rightChildNode->axis] + parentNode->bb.lowerCorner[parentNode->rightChildNode->axis])/2;
		}
		else if(SPLITTYPE == 2)
		{
			// here make mean splitting.
			return 0;
		}
		else if(SPLITTYPE == 3)
		{
			std::vector<double> centers;
			std::list<IElement*>::const_iterator element;
			BSPTreeNode * curNode;
			if(isLeft)
			{
				curNode = parentNode->leftChildNode;
				for(element = parentNode->leftChildNode->elements.begin(); element != parentNode->leftChildNode->elements.end(); element++)
				{
					centers.push_back((*element)->getCentroid()[parentNode->leftChildNode->axis]);
				}
			}
			else
			{
				curNode = parentNode->rightChildNode;
				for(element = parentNode->rightChildNode->elements.begin(); element != parentNode->rightChildNode->elements.end(); element++)
				{
					centers.push_back((*element)->getCentroid()[parentNode->rightChildNode->axis]);
				}
			}
			std::sort(centers.begin(),centers.end());
			if(curNode->elements.size() % 2 != 0)
				return centers.at(int(curNode->elements.size()/2));
			else
				return (centers.at(int(curNode->elements.size()/2)) + centers.at(int(curNode->elements.size()/2 - 1/2))) / 2;
		}
		else if(SPLITTYPE == 5)
		{
		    if(BSPTYPE == 1)
		    {
                // BSP Split
                BSPTreeNode * curNode;
                if(isLeft)
                {
                    curNode = parentNode->leftChildNode;
                }
                else
                {
                    curNode = parentNode->rightChildNode;
                }

                // BSP Splitting Here.
                std::vector<IElement*> myVect;
                std::list<IElement*>::const_iterator element3;

                for(element3 = curNode->elements.begin(); element3 != curNode->elements.end(); element3++)
                {
                    myVect.push_back(*element3);
                }

                int mySize = myVect.size();
                int score = 0;
                int minScore = mySize / 10;

                Vector3D pt1;
                Vector3D norm1;
                double d1;
                MeshVertex* vt1;
                MeshVertex* vt2;
                MeshVertex* vt3;

                Vector3D vt1t2;
                Vector3D vt1t3;

                do
                {
                    score = 0;

                    int selected = rand() % mySize;

                    IElement* myElem = myVect[selected];

                    int ret2 = 0;

                    vt1 = ((MeshTetrahedron*)myElem)->getVertex(0);
                    vt2 = ((MeshTetrahedron*)myElem)->getVertex(1);
                    vt3 = ((MeshTetrahedron*)myElem)->getVertex(2);

                    vt1t2 = *(vt2->getPosition()) - *(vt1->getPosition());
                    vt1t3 = *(vt3->getPosition()) - *(vt1->getPosition());

                    pt1 = (*(vt1->getPosition()) + *(vt2->getPosition()) + *(vt3->getPosition()))/3;
                    Vector3DCross(&norm1, &vt1t2, &vt1t3);

                    d1 = (-1) * (pt1[0]*norm1[0] + pt1[1]*norm1[1] + pt1[2]*norm1[2]);

                    std::list<IElement*>::const_iterator element4;

                    for(element4 = curNode->elements.begin(); element4 != curNode->elements.end(); element4++)
                    {
                        ret2 = classifyPointPlane((*element4)->getCentroid(), norm1, d1);

                        if(ret2 == CLASSIFY_LEFT)
                        {
                            score++;
                        }
                        else
                        {
                            score--;
                        }
                    }

                    // std::cout << "SCORE: " << score << " MINSCORE: " << minScore << " SIZE: " << mySize << std::endl;
                }
                while(abs(score) > minScore);

                curNode->splitPoint = pt1;
                curNode->planeNormal = norm1;
                curNode->d = d1;
		    }
		    else if(BSPTYPE == 2)
		    {
		        BSPTreeNode * curNode;
                if(isLeft)
                {
                    curNode = parentNode->leftChildNode;
                }
                else
                {
                    curNode = parentNode->rightChildNode;
                }

                int sz = curNode->elements.size();

                double* points = new double[sz * 3];
                Vector3D elem3;
                double totx = 0, toty = 0, totz = 0;
                double* orig = new double[3];
                double* dir = new double[3];

                int k = 0;
                std::list<IElement*>::const_iterator elementt;
                for(elementt = curNode->elements.begin(); elementt != curNode->elements.end(); elementt++)
                {
                    elem3 = (*elementt)->getCentroid();

                    points[3*k] = elem3[0];
                    totx += points[3*k];
                    points[3*k+1] = elem3[1];
                    toty += points[3*k+1];
                    points[3*k+2] = elem3[2];
                    totz += points[3*k+2];
                    k++;
                }

                totx /= sz;
                toty /= sz;
                totz /= sz;

                orig[0] = totx;
                orig[1] = toty;
                orig[2] = totz;

                // calculate inertia axis here and put it into orig[3] and dir[3]

                double px,py,pz = 0;
                double tx,ty,tz = 0;

                // std::cout << px << " " << py << " " << pz << std::endl;

                for(int i=0; i<sz; i++)
                {
                    points[3*i] -= orig[0];
                    points[3*i+1] -= orig[1];
                    points[3*i+2] -= orig[2];
                }


                double xp = 0;

                px = rand()%10000;
                py = rand()%10000;
                pz = rand()%10000;

                double tot = sqrt(px*px + py*py + pz*pz);

                px /= tot;
                py /= tot;
                pz /= tot;

                for(int i=0; i<REPSIZE; i++)
                {
                    tx = 0;
                    ty = 0;
                    tz = 0;

                    for(int j=0; j<sz; j++)
                    {
                        xp = px*points[3*i] + py*points[3*i+1] + pz*points[3*i+2];

                        tx += xp * points[3*i];
                        ty += xp * points[3*i+1];
                        tz += xp * points[3*i+2];
                    }

                    tot = sqrt(tx*tx + ty*ty + tz*tz);
                    px = tx/tot;
                    py = ty/tot;
                    pz = tz/tot;
                }

                dir[0] = px;
                dir[1] = py;
                dir[2] = pz;

                orig[0] = round(orig[0], ROUNDFACT);
                orig[1] = round(orig[1], ROUNDFACT);
                orig[2] = round(orig[2], ROUNDFACT);
                dir[0] = round(dir[0], ROUNDFACT);
                dir[1] = round(dir[1], ROUNDFACT);
                dir[2] = round(dir[2], ROUNDFACT);

                curNode->splitPoint = Vector3D(orig[0], orig[1], orig[2]);
                curNode->planeNormal = Vector3D(dir[0], dir[1], dir[2]);

                curNode->d = (-1) * (dir[0]*orig[0] + dir[1]*orig[1] + dir[2]*orig[2]);

                delete [] points;
		    }

            return 0;
		}

	}

	// Parent nodeun splitting coordinate ine ve axisine gore elemanlari right and left node'a koyuyor. Aslinda pointerlari sadece.
	void fillChildNodes(BSPTreeNode * parentNode) //OK
	{
	    if(SPLITTYPE == 5)
	    {
	        Vector3D res ;
	        int res2;
	        if (parentNode->elements.size() > 0)
	        {
                std::list<IElement*>::const_iterator element;
                for (element = parentNode->elements.begin(); element != parentNode->elements.end(); element++)
                {

                    res = (*element)->getCentroid();
                    res2  = classifyPointPlane(res, parentNode->planeNormal, parentNode->d );

                    if( res2 == CLASSIFY_LEFT)
                    {
                        parentNode->leftChildNode->elements.push_back( (*element) );
                    }
                    else if( res2 == CLASSIFY_RIGHT )
                    {
                        parentNode->rightChildNode->elements.push_back( (*element) );
                    }
                    else	// if bounding box of the element is on both sides (if it intersects the split plane) add to both
                    {
                        parentNode->leftChildNode->elements.push_back( (*element) );
                        parentNode->rightChildNode->elements.push_back( (*element) );
                    }
                }

                parentNode->elements.clear();								// parent node da hic eleman birakmiyoruz.
            }
	    }
	    else
	    {
            AABB bb;
            int axis = parentNode->axis;
            int i=0;

            if (parentNode->elements.size() > 0) {
                std::list<IElement*>::const_iterator element;
                for (element = parentNode->elements.begin(); element != parentNode->elements.end(); element++) {

                    // std::cout << i+1 << std::endl;
                    bb = (*element)->getBB();
                    // std::cout << "=====-===" << std::endl;

                    if( bb.upperCorner[axis] < parentNode->splitCoordinate )
                    {
                        parentNode->leftChildNode->elements.push_back( (*element) );
                    }
                    else if( bb.lowerCorner[axis] > parentNode->splitCoordinate )
                    {
                        parentNode->rightChildNode->elements.push_back( (*element) );
                    }
                    else	// if bounding box of the element is on both sides (if it intersects the split plane) add to both
                    {


                        parentNode->leftChildNode->elements.push_back( (*element) );
                        parentNode->rightChildNode->elements.push_back( (*element) );
                    }
                    i++;
                }

                parentNode->elements.clear();								// parent node da hic eleman birakmiyoruz.
            }
	    }
	}

    void createGeoFile(const char *inFileName, const char* outFileName) // OK
    {
        /// Edit here.

        fileName = inFileName;

        std::string line;
        std::string foundname;
        size_t found;
        std::ifstream myfile (inFileName);
        std::ofstream myfile2 (outFileName);
        if (myfile.is_open())
        {
            while ( myfile.good() )
            {
                /// Here find the main object. Name after tlo command.
                /// Comment it out.

                getline (myfile,line);
                found = line.find("tlo");
                if(found != std::string::npos)
                {
                    std::cout << line << std::endl;
                    foundname = &line[found+3];
                    line.insert(0,"#");
                }
                myfile2 << line << std::endl;
            }

            /// here write solid main_ = buldugumuz; Gerekirse sonuna ; koy.
            myfile2 << "solid " << rootNode->partname << " = " << foundname <<  " " << std::endl << std::endl;

            /// Fill the new parts here. Create cutplanes and partitions.

            writeCutPlanes(rootNode, myfile2);

            myfile2 << std::endl;
            /// Now write tlo command for each of the partitions.
            /// Traverse the Tree. It's as easy as that.

            writeTlo(rootNode, myfile2);


            /// Close the files.
            myfile.close();
            myfile2.close();
        }
        else
        {
            std::cerr << "ERROR ENCOUNTERED" << std::endl;
            std::cout << "Unable to open file" << std::endl;
        }

        return;
    }

    void writeCutPlanes(BSPTreeNode* node2, std::ofstream & myfile2) //OK
    {

        /// Fill Here the cutplanes and create the partitions from them.
        /// Make breadth first search.
        if(node2->elements.empty())
        {
            if(SPLITTYPE == 5)
            {
                myfile2 << "solid " << node2->planename << " = plane("<< node2->splitPoint[0] << ","<< node2->splitPoint[1] <<","<< node2->splitPoint[2] <<";"
                << node2->planeNormal[0] <<","<< node2->planeNormal[1] <<","<< node2->planeNormal[2] <<");" << std::endl;
            }
            else
            {
                switch (node2->axis) {
                    case 0:
                        myfile2 << "solid " << node2->planename << " = plane("<< node2->splitCoordinate << ",0,0;1,0,0);" << std::endl;
                        break;

                    case 1:
                        myfile2 << "solid " << node2->planename << " = plane(0,"<< node2->splitCoordinate << ",0;0,1,0);" << std::endl;
                        break;

                    case 2:
                        myfile2 << "solid " << node2->planename << " = plane(0,0,"<< node2->splitCoordinate << ";0,0,1);" << std::endl;
                        break;

                    default:
                        break;
                }
            }
            myfile2 << "solid " << node2->rightChildNode->partname << " = " << node2->partname << " and not " << node2->planename << ";" << std::endl;
            myfile2 << "solid " << node2->leftChildNode->partname << " = " << node2->partname << " and " << node2->planename << ";" << std::endl;
            myfile2 << std::endl;
        }

        /*
        1  procedure BFS(Graph,source):
        2      create a queue Q
        3      enqueue source onto Q
        4      mark source
        5      while Q is not empty:
        6          dequeue an item from Q into v
        7          for each edge e incident on v in Graph:
        8              let w be the other end of e
        9              if w is not marked:
        10                 mark w
        11                 enqueue w onto Q

        */

        //if(!node2->elements.empty())
          //  myfile2 << "tlo " << node2->partname << std::endl;

        if(node2->leftChildNode != NULL)
            writeCutPlanes(node2->leftChildNode, myfile2);

        if(node2->rightChildNode != NULL)
            writeCutPlanes(node2->rightChildNode, myfile2);
    }

    void writeTlo(BSPTreeNode* node2, std::ofstream & myfile2) //OK
    {
        if(!node2->elements.empty())
            myfile2 << "tlo " << node2->partname << ";" << std::endl;

        if(node2->leftChildNode != NULL)
            writeTlo(node2->leftChildNode, myfile2);

        if(node2->rightChildNode != NULL)
            writeTlo(node2->rightChildNode, myfile2);
    }

    void createPartitions(const char* inputFile, int minSize) //OK
    {
        using namespace nglib;
        Ng_Init();

        long nvertices, ntriangles;

        double *vertexlist;
        int *trianglelist;
        int *triangleFN;        // Boundary Element Index.
        int *triangleDIN;
        int *triangleDOUT;

        //vertexlist = new double[12000000];
        //trianglelist = new int[12000000];
        //triangleFN = new int[4000000];

        Ng_CSG_GenerateSurfaceMesh(inputFile, minSize, &nvertices, &ntriangles, vertexlist, trianglelist, triangleFN, triangleDIN, triangleDOUT);
        std::cout << "Successfully loaded .geo File: " << inputFile << std::endl;

        /// Now partition these elements to the tree nodes;

        std::cout << "Reading " << nvertices  << " points..." << std::endl;
        numvertexused = nvertices;
        elmVertexList = new ElmerVertex[nvertices];

        for (long int i = 0; i < nvertices; i++)
        {
		//std::cout << "Element Number: " << i+1 << " " << vertexlist[i*3] << std::endl;
            Vector3D position = Vector3D(vertexlist[i*3], vertexlist[i*3+1], vertexlist[i*3+2]);
            elmVertexList[i].setid(i+1);
            elmVertexList[i].setPosition(position);
            intersectBSPTreeVert(i);
            // i+1 stands for id.

        }
        std::cout << "done" << std::endl;

        std::cout << "Reading " << ntriangles  << " faces..." << std::endl;
        numtriangleused = ntriangles;
        elmTriangleList = new ElmerTriangle[ntriangles * 2];
        tricount = ntriangles;
        addTriNode = 0;

        for (long int i = 0; i < ntriangles; i++)
        {
            int npos = 3;
            Vector3D* positions = new Vector3D[npos]();
            int*pos = new int[npos];

            pos[0] = trianglelist[i*3];
            pos[1] = trianglelist[i*3+1];
            pos[2] = trianglelist[i*3+2];

            // std::cout << "ID: " << triangleFN[i] << std::endl;

            positions[0] = Vector3D(vertexlist[pos[0]*3-3], vertexlist[pos[0]*3-2], vertexlist[pos[0]*3-1]);
            positions[1] = Vector3D(vertexlist[pos[1]*3-3], vertexlist[pos[1]*3-2], vertexlist[pos[1]*3-1]);
            positions[2] = Vector3D(vertexlist[pos[2]*3-3], vertexlist[pos[2]*3-2], vertexlist[pos[2]*3-1]);

            elmTriangleList[i].setid(i+1);
            elmTriangleList[i].setPositions(Vector3D(pos[0],pos[1],pos[2]));
            elmTriangleList[i].setbid(triangleFN[i]);
            elmTriangleList[i].setDIN(triangleDIN[i]);
            elmTriangleList[i].setDOUT(triangleDOUT[i]);

            intersectBSPTreeFace(positions, pos, npos, i);
        }
        std::cout << "done" << std::endl;

        delete [] vertexlist;
        delete [] trianglelist;
        delete [] triangleFN;

        return;
    }

    /// This function should find which node the element belongs and puts it in there.
    /// It is also possible that the element will be put into more than one node.
    /// Elements can be vertices, faces or volume elements.
    /// in this case only vertices are present.
    void intersectBSPTreeVert(int id) // OK
    {
        intersectBSPTreeVertRec(rootNode, id);
        return;
    }

    void intersectBSPTreeVertRec(BSPTreeNode* node2, int id) //OK
    {
        if(node2->leftChildNode == NULL && node2->rightChildNode == NULL)
        {
            /// If we reached a leafNode then add this element to the list of verts.
            node2->verts[id+1] = &elmVertexList[id];
            /// Insert procs info.
            elmVertexList[id].addProc(node2->partId);
            tot++;
            // std::cout << node2->partname << "  " << id << ": " << pos[0] << "  " << pos[1] << "  " << pos[2] << "  " << node2->axis << "  " << node2->splitCoordinate << std::endl;
        }
        else    /// If we are not at a terminal node.
        {
            /// Check whether the element is to the right or the left of the cutplane.
            /// If it is on the plane call both childnodes.
            if(SPLITTYPE == 5)
            {
                int result = 0;
                result = classifyPointPlane(elmVertexList[id].getPosition(), node2->planeNormal, node2->d);

                if(result == CLASSIFY_ONPLANE)
                {
                    intersectBSPTreeVertRec(node2->leftChildNode, id);
                    intersectBSPTreeVertRec(node2->rightChildNode, id);
                }
                else if(result == CLASSIFY_LEFT)
                {
                    intersectBSPTreeVertRec(node2->leftChildNode, id);
                }
                else
                {
                    intersectBSPTreeVertRec(node2->rightChildNode, id);
                }
            }
            else
            {
                switch (node2->axis) {
                    case 0:     /// Check X value against the plane.(splitCoordinate)
                        if(abs(node2->splitCoordinate - elmVertexList[id].getPosition()[0]) < 0.00001)
                        {
                            // Point is on the plane. Put them into both the nodes.
                            //std::cout << " LR " << node2->partname << "  " << id << ": " << pos[0] << "  " << pos[1] << "  " << pos[2] << "  " << node2->splitCoordinate - pos[0] << std::endl;
                            intersectBSPTreeVertRec(node2->leftChildNode, id);
                            intersectBSPTreeVertRec(node2->rightChildNode, id);
                        }
                        else if((node2->splitCoordinate - elmVertexList[id].getPosition()[0]) >= 0.00001)
                        {
                            /// Put to the left.
                            //std::cout << node2->partname << " L " << id << ": " << pos[0] << "  " << pos[1] << "  " << pos[2] << "  " << node2->axis << "  " << node2->splitCoordinate - pos[0] << std::endl;
                            intersectBSPTreeVertRec(node2->leftChildNode, id);
                        }
                        else
                        {
                            /// Put to the right.
                            //std::cout << node2->partname << " R " << id << ": " << pos[0] << "  " << pos[1] << "  " << pos[2] << "  " << node2->axis << "  " << node2->splitCoordinate - pos[0] << std::endl;
                            intersectBSPTreeVertRec(node2->rightChildNode, id);
                        }
                        break;

                    case 1:     /// Check Y value against the plane.(splitCoordinate)
                        if(abs(node2->splitCoordinate - elmVertexList[id].getPosition()[1]) < 0.00001)
                        {
                            // Point is on the plane. Put them into both the nodes.
                            //std::cout << " LR " << node2->partname << "  " << id << ": " << pos[0] << "  " << pos[1] << "  " << pos[2] << "  " << node2->splitCoordinate - pos[0] << std::endl;
                            intersectBSPTreeVertRec(node2->leftChildNode, id);
                            intersectBSPTreeVertRec(node2->rightChildNode, id);
                        }
                        else if((node2->splitCoordinate - elmVertexList[id].getPosition()[1]) >= 0.00001)
                        {
                            /// Put to the left.
                            //std::cout << node2->partname << " L " << id << ": " << pos[0] << "  " << pos[1] << "  " << pos[2] << "  " << node2->axis << "  " << node2->splitCoordinate - pos[0] << std::endl;
                            intersectBSPTreeVertRec(node2->leftChildNode, id);
                        }
                        else
                        {
                            /// Put to the right.
                            //std::cout << node2->partname << " R " << id << ": " << pos[0] << "  " << pos[1] << "  " << pos[2] << "  " << node2->axis << "  " << node2->splitCoordinate - pos[0] << std::endl;
                            intersectBSPTreeVertRec(node2->rightChildNode, id);
                        }
                        break;

                    case 2:     /// Check Z value against the plane.(splitCoordinate)
                        if(abs(node2->splitCoordinate - elmVertexList[id].getPosition()[2]) < 0.00001)
                        {
                            // Point is on the plane. Put them into both the nodes.
                            //std::cout << " LR " << node2->partname << "  " << id << ": " << pos[0] << "  " << pos[1] << "  " << pos[2] << "  " << node2->splitCoordinate - pos[0] << std::endl;
                            intersectBSPTreeVertRec(node2->leftChildNode, id);
                            intersectBSPTreeVertRec(node2->rightChildNode, id);
                        }
                        else if((node2->splitCoordinate - elmVertexList[id].getPosition()[2]) >= 0.00001)
                        {
                            /// Put to the left.
                            //std::cout << node2->partname << " L " << id << ": " << pos[0] << "  " << pos[1] << "  " << pos[2] << "  " << node2->axis << "  " << node2->splitCoordinate - pos[0] << std::endl;
                            intersectBSPTreeVertRec(node2->leftChildNode, id);
                        }
                        else
                        {
                            /// Put to the right.
                            //std::cout << node2->partname << " R " << id << ": " << pos[0] << "  " << pos[1] << "  " << pos[2] << "  " << node2->axis << "  " << node2->splitCoordinate - pos[0] << std::endl;
                            intersectBSPTreeVertRec(node2->rightChildNode, id);
                        }
                        break;

                    default:
                        break;
                }
            }
        }
    }

    int classifyPointPlane(Vector3D pt, Vector3D norm, double d1) // OK
    {
        double result = pt[0]*norm[0] + pt[1]*norm[1] + pt[2]*norm[2] + d1;

        if(result <= -1e-6f)
        {
            return CLASSIFY_LEFT;
        }
        else if(result >= 1e-6f)
        {
            return CLASSIFY_RIGHT;
        }
        else
        {
            return CLASSIFY_ONPLANE;
        }
    }

    void intersectBSPTreeFace(Vector3D* positions, int*ids, int npos, int procid ) // OK
    {
        intersectBSPTreeFaceRec(rootNode, positions, ids, npos, procid, true);
        return;
    }

    void intersectBSPTreeFaceRec(BSPTreeNode* node2, Vector3D* pos, int*ids, int npos ,int procid, bool orient ) // OK
    {
        // int tmp;
        if(node2->leftChildNode == NULL && node2->rightChildNode == NULL)
        {
            /// If we reached a leafNode then add this element to the list of verts.
            node2->tris.push_back(&elmTriangleList[procid]);
            node2->orientation.push_back(orient);
            elmTriangleList[procid].addProc(node2->partId);
            // std::cout << node2->partname << "  " << ids[0] << "  " << ids[1] << "  " << ids[2] << std::endl;
        }
        else    /// If we are not at a terminal node.
        {
            /// Check whether the element is to the right or the left of the cutplane.
            /// If it is on the plane call both childnodes.
            /// In this case we check each of the vertices which constitute the face.
            /// Currently coded just to work with triangles.
            if(SPLITTYPE == 5)
            {
                int result1,result2,result3;
                result1 = classifyPointPlane(pos[0], node2->planeNormal, node2->d);
                result2 = classifyPointPlane(pos[1], node2->planeNormal, node2->d);
                result3 = classifyPointPlane(pos[2], node2->planeNormal, node2->d);

                if(result1 == CLASSIFY_ONPLANE && result2 == CLASSIFY_ONPLANE && result3 == CLASSIFY_ONPLANE)
                {
                    intersectBSPTreeFaceRec(node2->leftChildNode, pos, ids, npos, procid, orient);

                    //tmp = ids[0];
                    //ids[0] = ids[2];
                    //ids[2] = tmp;

                    // elmTriangleList[tricount + addTriNode] = new ElmerTriangle();
                    // elmTriangleList[tricount + addTriNode].setPositions(Vector3D(ids[0],ids[1],ids[2]));
                    // elmTriangleList[tricount + addTriNode].setid(tricount + addTriNode);
                    // addTriNode++;
                    // procid = tricount + addTriNode;

                    intersectBSPTreeFaceRec(node2->rightChildNode, pos, ids, npos, procid, !orient);
                }
                else if(result1 == CLASSIFY_LEFT || result2 == CLASSIFY_LEFT || result3 == CLASSIFY_LEFT)
                {
                    intersectBSPTreeFaceRec(node2->leftChildNode, pos, ids, npos, procid, orient);
                }
                else
                {
                    intersectBSPTreeFaceRec(node2->rightChildNode, pos, ids, npos, procid, orient);
                }
            }
            else
            {
                switch (node2->axis) {
                    case 0:     /// Check X value against the plane.(splitCoordinate)
                        if(abs(node2->splitCoordinate - pos[0][0]) < 0.00001 && abs(node2->splitCoordinate - pos[1][0]) < 0.00001 && abs(node2->splitCoordinate - pos[2][0]) < 0.00001)
                        {
                            // Point is on the plane. Put them into both the nodes.
                            //std::cout << node2->partname << " LR " << ids[0] << "  " << ids[1] << "  " << ids[2] << std::endl;
                            intersectBSPTreeFaceRec(node2->leftChildNode, pos, ids, npos, procid, orient);

                            //Vector3D vec1 = pos[1] - pos[0];
                            //Vector3D vec2 = pos[2] - pos[0];
                            //Vector3D faceNormal = vec1.Cross(vec2).Normalize();

                            //std::cout << "ID: " << ids[0] << "  " << ids[1] << "  " << ids[2] << std::endl;
                            //std::cout << "Face Normal: " << faceNormal[0] << " " << faceNormal[1] << " " << faceNormal[2] << std::endl;
                            //std::cout << "Plane Normal: " << node2->planeNormal[0] << " " << node2->planeNormal[1] << " " << node2->planeNormal[2] << std::endl;

                            //tmp = ids[0];
                            //ids[0] = ids[2];
                            //ids[2] = tmp;

                            //elmTriangleList[tricount + addTriNode] = new ElmerTriangle();
                            //elmTriangleList[tricount + addTriNode].setPositions(Vector3D(ids[0],ids[1],ids[2]));
                            //elmTriangleList[tricount + addTriNode].setid(tricount + addTriNode);
                            //addTriNode++;
                            //procid = tricount + addTriNode;

                            intersectBSPTreeFaceRec(node2->rightChildNode, pos, ids, npos, procid, !orient);
                        }
                        else if((node2->splitCoordinate - pos[0][0]) >= 0.00001 || (node2->splitCoordinate - pos[1][0]) >= 0.00001 || (node2->splitCoordinate - pos[2][0]) >= 0.00001)
                        {
                            /// Put to the left.
                            //std::cout << node2->partname << " L " << ids[0] << "  " << ids[1] << "  " << ids[2] << std::endl;
                            intersectBSPTreeFaceRec(node2->leftChildNode, pos, ids, npos, procid, orient);
                        }
                        else
                        {
                            /// Put to the right.
                            //std::cout << node2->partname << " R " << ids[0] << "  " << ids[1] << "  " << ids[2] << std::endl;
                            intersectBSPTreeFaceRec(node2->rightChildNode, pos, ids, npos, procid, orient);
                        }
                        break;

                    case 1:     /// Check Y value against the plane.(splitCoordinate)
                        if(abs(node2->splitCoordinate - pos[0][1]) < 0.00001 && abs(node2->splitCoordinate - pos[1][1]) < 0.00001 && abs(node2->splitCoordinate - pos[2][1]) < 0.00001)
                        {
                            // Point is on the plane. Put them into both the nodes.
                            //std::cout << node2->partname << " LR " << ids[0] << "  " << ids[1] << "  " << ids[2] << std::endl;
                            intersectBSPTreeFaceRec(node2->leftChildNode, pos, ids, npos, procid, orient);

                            //tmp = ids[0];
                            //ids[0] = ids[2];
                            //ids[2] = tmp;

                            //Vector3D vec1 = pos[1] - pos[0];
                            //Vector3D vec2 = pos[2] - pos[0];
                            //Vector3D faceNormal = vec1.Cross(vec2).Normalize();

                            //std::cout << "ID: " << ids[0] << "  " << ids[1] << "  " << ids[2] << std::endl;
                            //std::cout << "Face Normal: " << faceNormal[0] << " " << faceNormal[1] << " " << faceNormal[2] << std::endl;
                            //std::cout << "Plane Normal: " << node2->planeNormal[0] << " " << node2->planeNormal[1] << " " << node2->planeNormal[2] << std::endl;

                            //elmTriangleList[tricount + addTriNode] = new ElmerTriangle();
                            //elmTriangleList[tricount + addTriNode].setPositions(Vector3D(ids[0],ids[1],ids[2]));
                            //elmTriangleList[tricount + addTriNode].setid(tricount + addTriNode);
                            //addTriNode++;
                            //procid = tricount + addTriNode;

                            intersectBSPTreeFaceRec(node2->rightChildNode, pos, ids, npos, procid, !orient);
                        }
                        else if((node2->splitCoordinate - pos[0][1]) >= 0.00001 || (node2->splitCoordinate - pos[1][1]) >= 0.00001 || (node2->splitCoordinate - pos[2][1]) >= 0.00001)
                        {
                            /// Put to the left.
                            //std::cout << node2->partname << " L " << ids[0] << "  " << ids[1] << "  " << ids[2] << std::endl;
                            intersectBSPTreeFaceRec(node2->leftChildNode, pos, ids, npos, procid, orient);
                        }
                        else
                        {
                            /// Put to the right.
                            //std::cout << node2->partname << " R " << ids[0] << "  " << ids[1] << "  " << ids[2] << std::endl;
                            intersectBSPTreeFaceRec(node2->rightChildNode, pos, ids, npos, procid, orient);
                        }
                        break;

                    case 2:     /// Check Z value against the plane.(splitCoordinate)
                        if(abs(node2->splitCoordinate - pos[0][2]) < 0.00001 && abs(node2->splitCoordinate - pos[1][2]) < 0.00001 && abs(node2->splitCoordinate - pos[2][2]) < 0.00001)
                        {
                            // Point is on the plane. Put them into both the nodes.
                            //std::cout << node2->partname << " LR " << ids[0] << "  " << ids[1] << "  " << ids[2] << std::endl;
                            intersectBSPTreeFaceRec(node2->leftChildNode, pos, ids, npos, procid, orient);

                            // tmp = ids[0];
                            // ids[0] = ids[2];
                            // ids[2] = tmp;

                            //Vector3D vec1 = pos[1] - pos[0];
                            //Vector3D vec2 = pos[2] - pos[0];
                            //Vector3D faceNormal = vec1.Cross(vec2).Normalize();


                            //std::cout << "ID: " << ids[0] << "  " << ids[1] << "  " << ids[2] << std::endl;
                            //std::cout << "Face Normal: " << faceNormal[0] << " " << faceNormal[1] << " " << faceNormal[2] << std::endl;
                            //std::cout << "Plane Normal: " << node2->planeNormal[0] << " " << node2->planeNormal[1] << " " << node2->planeNormal[2] << std::endl;

                            //elmTriangleList[tricount + addTriNode] = new ElmerTriangle();
                            //elmTriangleList[tricount + addTriNode].setPositions(Vector3D(ids[0],ids[1],ids[2]));
                            //elmTriangleList[tricount + addTriNode].setid(tricount + addTriNode);
                            //addTriNode++;
                            //procid = tricount + addTriNode;

                            intersectBSPTreeFaceRec(node2->rightChildNode, pos, ids, npos, procid, !orient);
                        }
                        else if((node2->splitCoordinate - pos[0][2]) >= 0.00001 || (node2->splitCoordinate - pos[1][2]) >= 0.00001 || (node2->splitCoordinate - pos[2][2]) >= 0.00001)
                        {
                            /// Put to the left.
                            //std::cout << node2->partname << " L " << ids[0] << "  " << ids[1] << "  " << ids[2] << std::endl;
                            intersectBSPTreeFaceRec(node2->leftChildNode, pos, ids, npos, procid, orient);
                        }
                        else
                        {
                            /// Put to the right.
                            //std::cout << node2->partname << " R " << ids[0] << "  " << ids[1] << "  " << ids[2] << std::endl;
                            intersectBSPTreeFaceRec(node2->rightChildNode, pos, ids, npos, procid, orient);
                        }
                        break;

                    default:
                        break;
                }
            }
        }
    }


    void createPartFiles() // OK
    {
        createRecursive(rootNode);
    }

    void createRecursive(BSPTreeNode* node2) // OK
    {

        if(node2->leftChildNode == NULL && node2->rightChildNode == NULL)
        {
            /// We reached a terminal node. So Write its file.
            /// Now we will also write to a string.
            std::stringstream sfile;

            // const char* outFile = node2->partname.c_str();
            // std::ofstream myfile2 (outFile);

            // myfile2 << node2->verts.size() << std::endl;
            sfile << node2->verts.size() << std::endl;

            std::map<int,ElmerVertex*>::iterator vectlistit;

            for (vectlistit = node2->verts.begin(); vectlistit != node2->verts.end(); vectlistit++)
            {
                // myfile2 << (*vectlistit).second->getid() << " " << (*vectlistit).second->getPosition()[0] << " " << (*vectlistit).second->getPosition()[1] << " "
                // << (*vectlistit).second->getPosition()[2] << " " << (*vectlistit).second->getProcs() << std::endl;

                sfile << (*vectlistit).second->getid() << " " << (*vectlistit).second->getPosition()[0] << " " << (*vectlistit).second->getPosition()[1] << " "
                << (*vectlistit).second->getPosition()[2] << " " << (*vectlistit).second->getProcs() << std::endl;
            }

            // myfile2 << node2->tris.size() << std::endl;
            sfile << node2->tris.size() << std::endl;

            std::list<ElmerTriangle*>::iterator triit;
            std::list<bool>::iterator orit;
            std::string proclist;
            int x1,y1,z1, id1, fid;
            x1 = y1 = z1 = id1 = fid = 0;
            int tmp1;
            // bool isfound;

            std::map<std::string, int> edgemap;
            std::map<std::string, int>::iterator edgeit;
            std::string tmpstr;
            std::stringstream ss;

            /// Here find the region which is inside of the object.
            std::map<int, int> regmap;
            int tmmm33;

            // Count the number.
            for(triit = node2->tris.begin(); triit != node2->tris.end(); triit++)
            {
                tmmm33 = (*triit)->getDIN();
                if(regmap.find(tmmm33) != regmap.end())
                {
                    regmap[tmmm33] += 1;
                }
                else
                {
                    regmap[tmmm33] = 1;
                }

                tmmm33 = (*triit)->getDOUT();
                if(regmap.find(tmmm33) != regmap.end())
                {
                    regmap[tmmm33] += 1;
                }
                else
                {
                    regmap[tmmm33] = 1;
                }
            }

            std::map<int, int>::iterator regmapit;
            int maxelem = 0;
            int maxcount = 0;

            for(regmapit = regmap.begin(); regmapit != regmap.end(); regmapit++)
            {
                if((*regmapit).second > maxcount)
                {
                    maxelem = (*regmapit).first;
                    maxcount = (*regmapit).second;
                }
            }

            std::cout << "///////////////////////////////////////////" << std::endl;
            std::cout << "HERE THE MAX ELEM: " << maxelem << " " << maxcount << " " << regmap.size() << std::endl;
            std::cout << "///////////////////////////////////////////" << std::endl;

            /// Now we have the domain index. Lets decide whether to change face orientation.

            for (triit = node2->tris.begin(), orit = node2->orientation.begin(); triit != node2->tris.end(); triit++, orit++)
            {
                tmmm33 = (*triit)->getDIN();
                if(tmmm33 != maxelem)
                {
                    (*orit) = false;
                }
                else
                {
                    (*orit) =  true;
                }
            }

            for (triit = node2->tris.begin(), orit = node2->orientation.begin(); triit != node2->tris.end(); triit++, orit++)
            {
                id1 = (*triit)->getid();
                x1 = (int)(*triit)->getPositions()[0];
                y1 = (int)(*triit)->getPositions()[1];
                z1 = (int)(*triit)->getPositions()[2];
                fid = (*triit)->getbid();
                proclist = (*triit)->getProcs();

                /*
                /// Check Edge Map. If exists change orientation.
                isfound = false;

                ss.clear();
                ss.str("");
                ss << x1 << ";" << y1;
                // std::cout << ss.str() << std::endl;
                if(edgemap.find(ss.str()) == edgemap.end())
                {
                    /// Add to the edgelist.
                    edgemap[ss.str()] = id1;
                }
                else
                {
                    std::cout << ss.str() << std::endl;
                    isfound = true;
                    /// Change orientation and add to list.
                }
                ss.clear();
                ss.str("");
                ss << y1 << ";" << z1;
                // std::cout << ss.str() << std::endl;
                if(edgemap.find(ss.str()) == edgemap.end())
                {
                    /// Add to the edgelist.
                    edgemap[ss.str()] = id1;
                }
                else
                {
                    std::cout << ss.str() << std::endl;
                    isfound = true;
                    /// Change orientation and add to list.
                }
                ss.clear();
                ss.str("");
                ss << z1 << ";" << x1;
                // std::cout << ss.str() << std::endl;
                if(edgemap.find(ss.str()) == edgemap.end())
                {
                    /// Add to the edgelist.
                    edgemap[ss.str()] = id1;
                }
                else
                {
                    std::cout << ss.str() << std::endl;
                    isfound = true;
                    /// Change orientation and add to list.
                }
                */

                if(!(*orit))
                {
                    tmp1 = x1;
                    x1 = z1;
                    z1 = tmp1;

                    /*
                    ss.clear();
                    ss.str("");
                    ss << x1 << ";" << y1;
                    edgemap[ss.str()] = id1;

                    ss.clear();
                    ss.str("");
                    ss << y1 << ";" << z1;
                    edgemap[ss.str()] = id1;

                    ss.clear();
                    ss.str("");
                    ss << z1 << ";" << x1;
                    edgemap[ss.str()] = id1;
                    */
                }

                // myfile2 << id1 << " " << x1 << " " << y1 << " " << z1 << " " << fid << " " << proclist << std::endl;
                sfile << id1 << " " << x1 << " " << y1 << " " << z1 << " " << fid << " " << proclist << std::endl;
                // std::cout << id1 << " " << x1 << " " << y1 << " " << z1 << std::endl;
            }

            // myfile2.close();

            /// Here write sfile to its string in memory.

            stringtosend[sindex] = sfile.str();
            sindex++;
            sfile.str(std::string());	// Clear content

            return;
        }

        if(node2->leftChildNode != NULL)
            createRecursive(node2->leftChildNode);

        if(node2->rightChildNode != NULL)
            createRecursive(node2->rightChildNode);
    }

    void mpifunc(int rank, int size)
    {
        using namespace nglib;
        using namespace std;

        std::stringstream sfile;
        std::string mydata4;
        int datasize = 0;

        if(rank == 0)
        {
            MPI_Request* reqs = new MPI_Request[size-1];
            MPI_Status* stats = new MPI_Status[size-1];

            datasize = stringtosend[0].length();

            int* sizetmp = new int[size-1] ;
            for(int i=1; i<size; i++)
            {
                sizetmp[i-1] = stringtosend[i].length();
                MPI_Isend(&sizetmp[i-1], 1, MPI_INT, i, 123, MPI_COMM_WORLD, &reqs[i-1]);
            }

            MPI_Waitall(size-1, reqs, stats);

            delete [] sizetmp;
            delete [] reqs;
            delete [] stats;
        }
        else
        {
            MPI_Request reqe;
            MPI_Status state;
            MPI_Irecv(&datasize, 1, MPI_INT, 0, 123, MPI_COMM_WORLD, &reqe);
            MPI_Wait(&reqe, &state);
        }

        /// Here we have the datasizes. Now we need to get the data.
        mydata4.resize(datasize);

        if(rank == 0)
        {
            MPI_Request* reqs = new MPI_Request[size-1];
            MPI_Status* stats = new MPI_Status[size-1];

            for(int i=1; i<size; i++)
            {
                MPI_Isend(&stringtosend[i][0], stringtosend[i].length(), MPI_CHAR, i, 123, MPI_COMM_WORLD, &reqs[i-1]);
            }

            MPI_Waitall(size-1, reqs, stats);

            mydata4 = stringtosend[0];

            delete [] reqs;
            delete [] stats;
        }
        else
        {
            MPI_Request reqe;
            MPI_Status state;
            MPI_Irecv(&mydata4[0], datasize, MPI_CHAR, 0, 123, MPI_COMM_WORLD, &reqe);
            MPI_Wait(&reqe, &state);
        }

        /// Now we also have the data in mydata4

        int totallevel = (int)ceil(log10f(size)/log10f(2)) + 6;

        // partitioning.2
        // part.1.boundary
        /// part.1.elements     ..............DONE
        /// part.1.header       ..............DONE
        /// part.1.nodes        ..............DONE
        /// part.1.shared       ..............DONE

        // Here create the filename char *
        char * boundaryfile = new char[32];
        char * elementfile = new char[32];
        char * headerfile = new char[32];
        char * nodefile = new char[32];
        char * sharedfile = new char[32];
        char * repartfile = new char[32];


        /// Here create the partitioning.n directory first. Clean if already exists.

        sprintf(boundaryfile, "./partitioning.%d/part.%d.boundary", size, rank+1);
        sprintf(elementfile, "./partitioning.%d/part.%d.elements", size, rank+1);
        sprintf(headerfile, "./partitioning.%d/part.%d.header", size, rank+1);
        sprintf(nodefile, "./partitioning.%d/part.%d.nodes", size, rank+1);
        sprintf(sharedfile, "./partitioning.%d/part.%d.shared", size, rank+1);
        sprintf(repartfile, "./partitioning.%d/part.%d.repart", size, rank+1);

//        std::cout << "BOUNDARY FILE: " << boundaryfile << std::endl;
//        std::cout << "ELEMENT FILE: " << elementfile << std::endl;
//        std::cout << "HEADER FILE: " << headerfile << std::endl;
//        std::cout << "NODE FILE: " << nodefile << std::endl;
//        std::cout << "SHARED FILE: " << sharedfile << std::endl;
//        std::cout << "REPARTITION FILE: " << repartfile << std::endl;


        // std::cin >> totallevel;

        // cout << size << " " << log10f(size)/log10f(2) << " " << totallevel << endl;

        char * filename = new char[totallevel];
        // char * outFile = new char[totallevel + 4];
        filename[0] = 'm';
        filename[1] = 'a';
        filename[2] = 'i';
        filename[3] = 'n';
        filename[4] = '_';

        int tempp = size / 2;
        int tempp2 = rank;

        for(int i=5; i<totallevel-1; i++)
        {
            if(tempp2/tempp == 0)
                filename[i] = 'l';
            else
                filename[i] = 'r';

            tempp2 = tempp2%tempp;
            tempp = tempp/2;
        }
        filename[totallevel-1] = '\0';

        Ng_Mesh * mesh;

        Ng_Init();

        // creates mesh structure
        mesh = Ng_NewMesh ();

        int numpoint;
        int numface;
        int numboundary;

        int i, np, npshared, npd, nse, ne, proccnt, prccnttmp;
        npshared = 0;
        npd = 0;
        double point[3];
        int trig[3], tet[4] , ptid, triid;
        int fid;
        // int tetid;

        std::map<int,int> ind1;
        std::map<int,int> ind2;

        if(rank != -1)
        {
            // reads surface mesh from file
            ifstream in(filename);
            std::stringstream ifile;
            ifile.str(std::string());
            ifile << mydata4;

            ofstream outshared(sharedfile);
            ofstream outbound(boundaryfile);
            std::stringstream ss;

            //in >> np;
            ifile >> np;

            numpoint = np;
            cout << "PROCESSOR " << rank << " is reading " << np  << " points..." << std::endl;

            for (i = 1; i <= np; i++)
            {
                //in >> ptid >> point[0] >> point[1] >> point[2] >> proccnt;
                ifile >> ptid >> point[0] >> point[1] >> point[2] >> proccnt;

                ss.clear();
                ss.str("");

                ss << ptid << " " << proccnt;

                for(int j=0; j<proccnt; j++)
                {
                    //in >> prccnttmp;
                    ifile >> prccnttmp;
                    ss << " " << (prccnttmp+1);
                }

                if((proccnt > 1))
                {
                    /// Shared File is written here.

                    outshared << ss.str() << endl;
                    npshared++;
                }

                // cout << i << " " << ptid << " " << point[0] << "  " << point[1] << "  " << point[2] << endl;
                Ng_AddPoint (mesh, point);
                ind1[ptid] = i;
                ind2[i] = ptid;
                // Ng_SetPointIndex(mesh, i, ptid);
            }
            outshared.close();
            cout << "done" << endl;

            //in >> nse;
            ifile >> nse;

            numface = nse;
            numboundary = 0;
            int x11,y11,z11;
            int parid1, parid2;

            cout << "Reading " << nse  << " faces..."; cout.flush();
            for (i = 1; i <= nse; i++)
            {
                // cout << "ANOTHER ONE:" << endl;
                //in >> triid >> trig[0] >> trig[1] >> trig[2] >> fid >> proccnt;
                ifile >> triid >> trig[0] >> trig[1] >> trig[2] >> fid >> proccnt;

                for(int j=0; j<proccnt; j++)
                {
                    ifile >> prccnttmp;
                    //in >> prccnttmp;
                }

                // cout << trig[0] << "  " << trig[1] << "  " << trig[2] << std::endl;
                x11 = trig[0];
                y11 = trig[1];
                z11 = trig[2];
                trig[0] = ind1[x11];
                trig[1] = ind1[y11];
                trig[2] = ind1[z11];

                parid1 = 0;
                parid2 = 0;

                /// Format is:  ID  FACEID PARENTID1 PARENTID2 TYPE N1 ... N2
                /// Change parids
                /// Make sure the trig[n]'s are correct.
                /// Correct faceids index, add only nonshared faces to the file.
                if(proccnt == 1)
                {
                    outbound << numboundary+1 << " " << fid << " " << parid1 << " " << parid2 << " 303 " << trig[0] << " " << trig[1] << " " << trig[2] << endl;
                    numboundary++;
                }

                // cout << trig[0] << "  " << trig[1] << "  " << trig[2] << std::endl;
                Ng_AddSurfaceElement (mesh, NG_TRIG, trig);
                // Ng_SetSurfaceElementIndex(mesh, i, triid);
            }
            outbound.close();
            cout << "done" << endl;
        }
        /// Burdan Gerekirse cikar.
        // MPI_Barrier(MPI_COMM_WORLD);

        if(rank != -1)
        {
            // generate volume mesh
            Ng_Meshing_Parameters mp;
            mp.maxh = 1e8;
            mp.fineness = 1;
            mp.secondorder = 0;

            /// Volume Meshing START
            cout << "start meshing" << endl;
            Ng_GenerateVolumeMesh (mesh, &mp);
            cout << "meshing done" << endl;

            /// Volume Meshing END


            /// here the numbering of the nodes, elements etc is done.

            np = Ng_GetNP(mesh);
            cout << "Nodes: " << np << endl;

            ne = Ng_GetNE(mesh);
            cout << "Elements: " << ne << endl;

            int data1 = np - numpoint;
            int npindex = 0;
            int neindex = 0;
            int totvertt = 0;
            int tottrian = 0;
            int* vbuffer = new int[2];

            // If root node send numvertices and numtriangles used to the other processors.

            if(rank == 0)
            {
                vbuffer[0] = numvertexused;
                vbuffer[1] = numtriangleused;
            }

            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(vbuffer, 2 , MPI_INT, 0, MPI_COMM_WORLD);
            totvertt = vbuffer[0];
            tottrian = vbuffer[1];
            // std::cout << "RANK: " << rank << " : " << vbuffer[0] << " - " << vbuffer[1] << std::endl;

            // Then MPI_Scan the newly added np (np - numpoint) into npindex.
            MPI_Scan(&data1, &npindex, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            // Then MPI_Scan the new elements. (ne) into neindex.
            MPI_Scan(&ne, &neindex, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

            neindex -= ne;
            npindex -= data1;
            npindex += totvertt;
            npindex++;
            int npindex2 = npindex;

            std::cout << "RANK: " << rank << " : " << npindex << " - " << neindex << std::endl;

            /// PARMETIS SETUP START

            int * elmdist = new int[size+1];     // 0 5 10 15 (np+1 adet) element indexleri tutuyor.
            int * elmdist2 = new int[size+1];
            elmdist2[0] = 0;
            elmdist[0] = 0;
            elmdist2[rank+1] = ne;

            MPI_Allgather(&elmdist2[rank+1], 1, MPI_INT, &elmdist[1], 1, MPI_INT, MPI_COMM_WORLD);

            elmdist[0] = 1;
            for(int z=1; z<=size; z++)
            {
                elmdist[z] = elmdist[z] + elmdist[z-1];
            }


            if(rank == 0)
            {
                for(int z=0; z<=size; z++)
                {
                    cout << "ELMDIST: " << elmdist[z] << endl;
                }
            }


            int * eptr = new int[ne+1];          // 0 2 5 8 11 13 (ne+1 adet) Degisken olabiliyor.Ama bizim durumumuzda sadece tetrahedral var. Yani 0 4 8 12 16 ...
            int * elemts = new int[4*ne];        // [12 11 7 3] [11 12 8 9] Sirayla elemenlarin indexleri yaziliyor.
            eptr[0] = 0;

            // elmweight = NULL;
            int numflag = 0;    // indexing starts from 0 (set to 1 if index start from 1)
            int wgtflag = 0;    // not weighted
            int ncon = 1;
            int ncommonnodes = 3;   // for tetrahedral.
            int nparts = size;      // size of partition.
            float *tpwgts = new float[size]; // normally ncon*nparts.
            float *ubvec = new float[1];     // balance constraint
            int options[10];        // options.
            int edgecut;            // edges cut in total by the new partitioning.
            int *part = new int[ne];
            MPI_Comm comm;
            MPI_Comm_dup(MPI_COMM_WORLD, &comm);

            ubvec[0] = 1.05;

            for(int kk=0; kk<size; kk++)
            {
                tpwgts[kk] = 1.0/(float)(nparts);
            }

            options[0] = 1;
            options[PMV3_OPTION_DBGLVL] = 7;
            options[PMV3_OPTION_SEED] = 0;

            /// Fill in elmdist, eptr, elems
            /// PARMETIS SETUP END

            /// MIGRATION PART
            MeshMig* localMesh = new MeshMig(rank, size, false); // false for not reading data from files.
            localMesh->insertHeaderData(np, ne, numboundary, npshared);

            /// Here the nodes list START
            // ofstream outnode(nodefile);
            int iddd = 0;

            int *nodeidptr = new int[np];
            double *xptr = new double[np];
            double *yptr = new double[np];
            double *zptr = new double[np];

            for (i = 1; i <= np; i++)
            {
                Ng_GetPoint (mesh, i, point);

                if(i <= numpoint)
                {
                    iddd = ind2[i];
                }
                else
                {
                    iddd = npindex;
                    npindex++;
                }

                nodeidptr[i-1] = iddd;
                xptr[i-1] = point[0];
                yptr[i-1] = point[1];
                zptr[i-1] = point[2];

                // outnode << iddd << " -1 " << point[0] << " " << point[1] << " " << point[2] << endl;
            }
            localMesh->insertMeshNodes(nodeidptr, xptr, yptr, zptr);

            delete [] nodeidptr;
            delete [] xptr;
            delete [] yptr;
            delete [] zptr;
            // outnode.close();
            /// Here the nodes list END

            /// Here the elements list START
            // ofstream outelement(elementfile);
            int x11,y11,z11,w11;

            int *eglidptr = new int[ne];
            int *exptr = new int[ne];
            int *eyptr = new int[ne];
            int *ezptr = new int[ne];
            int *ewptr = new int[ne];

            for (i = 1; i <= ne; i++)
            {
                Ng_GetVolumeElement (mesh, i, tet);

                x11 = tet[0];
                y11 = tet[1];
                z11 = tet[2];
                w11 = tet[3];

                /// Global Numbering for the volume elements.
                /// Global Numbering for the tet[n]'s.

                if(tet[0] <= numpoint)
                    tet[0] = ind2[x11];
                else
                    tet[0] = npindex2 + tet[0] - numpoint - 1;

                if(tet[1] <= numpoint)
                    tet[1] = ind2[y11];
                else
                    tet[1] = npindex2 + tet[1] - numpoint - 1;

                if(tet[2] <= numpoint)
                    tet[2] = ind2[z11];
                else
                    tet[2] = npindex2 + tet[2] - numpoint - 1;

                if(tet[3] <= numpoint)
                    tet[3] = ind2[w11];
                else
                    tet[3] = npindex2 + tet[3] - numpoint - 1;


                eglidptr[i-1] = neindex + i;
                exptr[i-1] = tet[0];
                eyptr[i-1] = tet[1];
                ezptr[i-1] = tet[2];
                ewptr[i-1] = tet[3];

                // outelement << neindex+i << " 1 504 " << tet[0] << " " << tet[1] << " " << tet[2] << " " << tet[3] << endl;

                // Parmetis stuff
                eptr[i] = eptr[i-1] + 4;
                int w = i-1;
                elemts[w*4] = tet[0];
                elemts[w*4+1] = tet[1];
                elemts[w*4+2] = tet[2];
                elemts[w*4+3] = tet[3];
            }
            localMesh->insertMeshElements(eglidptr, exptr, eyptr, ezptr, ewptr);
            delete [] eglidptr;
            delete [] exptr;
            delete [] eyptr;
            delete [] ezptr;
            delete [] ewptr;
            // outelement.close();
            /// Here the elements list END

            /*
            /// Here the header file START
            ofstream outheader(headerfile);

            // npd contains dirichlet addition.
            // numboundary, np, npshared and ne contain number of nodes, shared nodes, elements, boundary elements.
            outheader << np << " " << ne << " " << numboundary + npd << endl;

            if(npd == 0)
            {
                outheader << 2 << endl;
                outheader << "504 " << ne << endl;
                outheader << "303 " << numboundary << endl;
            }
            else
            {
                outheader << 3 << endl;
                outheader << "504 " << ne << endl;
                outheader << "101 " << npd << endl;
                outheader << "303 " << numboundary << endl;
            }

            if(npshared != 0)
            {
                outheader << npshared << " 0" << endl;
            }

            outheader.close();
            /// Here the header file END
            */
            // std::cout << "Elmer Files Created by RANK: " << rank << std::endl;
            /// Boundary and Shared Files START

            ifstream in(filename);
            std::stringstream ifile;
            ifile.str(std::string());
            ifile << mydata4;

            std::stringstream ss;
            vector<string> sharedstring;
            vector<int> sharedids;

            int npx ,nsex;

            // in >> npx;
            ifile >> npx;
            for (i = 1; i <= npx; i++)
            {
                // in >> ptid >> point[0] >> point[1] >> point[2] >> proccnt;
                ifile >> ptid >> point[0] >> point[1] >> point[2] >> proccnt;

                ss.clear();
                ss.str("");

                ss << proccnt;

                for(int j=0; j<proccnt; j++)
                {
                    //in >> prccnttmp;
                    ifile >> prccnttmp;
                    ss << " " << (prccnttmp+1);
                }

                if((proccnt > 1))
                {
                    // Add to the shared list.
                    sharedstring.push_back(ss.str());
                    sharedids.push_back(ptid);
                    // npshared++;
                }
            }
            localMesh->insertMeshShared(&sharedids[0], &sharedstring[0]);
            sharedids.clear();
            sharedstring.clear();

            // in >> nsex;
            ifile >> nsex;
            numboundary = 0;     // Beginning local ID for boundary elements
            int parid1, parid2;

            int *gfidptr = new int[nsex];
            int *sxptr = new int[nsex];
            int *syptr = new int[nsex];
            int *szptr = new int[nsex];
            int *parptr = new int[nsex];

            for (i = 1; i <= nsex; i++)
            {
                // cout << "ANOTHER ONE:" << endl;
                // in >> triid >> trig[0] >> trig[1] >> trig[2] >> fid >> proccnt;
                ifile >> triid >> trig[0] >> trig[1] >> trig[2] >> fid >> proccnt;

                for(int j=0; j<proccnt; j++)
                {
                    //in >> prccnttmp;
                    ifile >> prccnttmp;
                }

                x11 = trig[0];
                y11 = trig[1];
                z11 = trig[2];

                if(trig[0] <= numpoint)
                    trig[0] = ind2[x11];
                else
                    trig[0] = npindex2 + trig[0] - numpoint - 1;

                if(trig[1] <= numpoint)
                    trig[1] = ind2[y11];
                else
                    trig[1] = npindex2 + trig[1] - numpoint - 1;

                if(trig[2] <= numpoint)
                    trig[2] = ind2[z11];
                else
                    trig[2] = npindex2 + trig[2] - numpoint - 1;

                parid1 = 0;
                parid2 = 0;

                /// Format is:  ID  FACEID PARENTID1 PARENTID2 TYPE N1 ... N2
                /// Change parids
                /// Make sure the trig[n]'s are correct.
                /// Correct faceids index, add only nonshared faces to the file.
                if(proccnt == 1)
                {
                    gfidptr[numboundary] = fid;
                    parptr[numboundary] = parid1;
                    sxptr[numboundary] = x11;
                    syptr[numboundary] = y11;
                    szptr[numboundary] = z11;

                    numboundary++;
                }
            }
            localMesh->insertMeshBoundary(gfidptr, sxptr, syptr, szptr, parptr);
            delete [] gfidptr;
            delete [] sxptr;
            delete [] syptr;
            delete [] szptr;
            delete [] parptr;

            /// Boundary and Shared Files END

            /// PARMETIS START
            /*
            MPI_Barrier(MPI_COMM_WORLD);
            ParMETIS_V3_PartMeshKway(elmdist, eptr, elemts, NULL, &wgtflag, &numflag, &ncon, &ncommonnodes, &nparts, tpwgts, ubvec, options, &edgecut, part, &comm);
            */
            /// PARMETIS END


            /// Migrate Mesh according to Parmetis result.

            // localMesh->insertMigrationData(ne, part);
            localMesh->setBoundaryElements();
            // localMesh->createElmerOutput();

            std::cout << "RANK: " << rank << " has reached the position" << std::endl;
            MPI_Barrier(MPI_COMM_WORLD);

            /*
            std::cout << "Migration begins. RANK: " << rank << std::endl;
            localMesh->mymigrate2(MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
            std::cout << "Migration Ends. RANK: " << rank << std::endl;
            MPI_Barrier(MPI_COMM_WORLD);
            */

            localMesh->refineParallel(MPI_COMM_WORLD, fileName);
            std::cout << "RANK: " << rank << " ended refinement." << std::endl;
            MPI_Barrier(MPI_COMM_WORLD);

            localMesh->refineParallel(MPI_COMM_WORLD, fileName);
            std::cout << "RANK: " << rank << " ended refinement." << std::endl;
            MPI_Barrier(MPI_COMM_WORLD);

            /// Write output like that.
            localMesh->createElmerOutput();

            /// MIGRATION PART END
            // volume mesh output

            /*
            np = Ng_GetNP(mesh);
            cout << "Points: " << np << endl;

            for (i = 1; i <= np; i++)
            {
                Ng_GetPoint (mesh, i, point);
                cout << i << ": " << point[0] << " " << point[1] << " " << point[2] << endl;
            }

            */

            /*
            MPI_Barrier(MPI_COMM_WORLD);

            int id21, id22;
            for(i = 1; i <= nse; i++)
            {
                std::cout << "Girmeden Once" << std::endl;
                id21 = id22 = Ng_GetVertex_SurfaceElements2(mesh,i);
                std::cout << "Element: " << i << " ParentID: " << id21 << "    " << id22 << std::endl;
            }

            // Create output file name

            for(int k=0; k<totallevel-1; k++)
            {
                outFile[k] = filename[k];
            }
            outFile[totallevel-1] = '.';
            outFile[totallevel] = 'v';
            outFile[totallevel+1] = 'o';
            outFile[totallevel+2] = 'l';
            outFile[totallevel+3] = '\0';

            Ng_SaveMesh(mesh, outFile);
            */
        }
        /* Quit */
        return;
    }


	void addElement(IElement* element) // OK
	{
		elementList.push_back(element);
	}

	void addSurfaceElem(IElement* element) // OK
	{
	    surfaceList.push_back(element);
	}

	BSPTreeNode * getRootNode() // OK
	{
		return rootNode;
	}

	unsigned long getTime()
	{
		time_t seconds;
		seconds = time(NULL);

		return (unsigned long)seconds;
	}

	void addMesh(Mesh* mesh)
	{
		meshList.push_back(mesh);
	}

	int gettotalNodes()
	{
		return totalnode;
	}

	int gettotalLeafs()
	{
		return totalleafnode;
	}

	void getPlane(double *liste, BSPTreeNode* node2)
	{
		liste[0] = node2->ll[0];
		liste[1] = node2->ll[1];
		liste[2] = node2->ll[2];
		liste[3] = node2->lr[0];
		liste[4] = node2->lr[1];
		liste[5] = node2->lr[2];
		liste[6] = node2->ur[0];
		liste[7] = node2->ur[1];
		liste[8] = node2->ur[2];
		liste[9] = node2->ul[0];
		liste[10] = node2->ul[1];
		liste[11] = node2->ul[2];

	}

	int getNumMesh()
	{
		return (unsigned int)meshList.size();
	}

	Mesh* getMeshAt(int index)
	{
		return meshList.at(index);
	}

	Vector3D getPlaneColor(BSPTreeNode* node2)
	{
		return node2->planeclr;
	}

	BSPTreeNode* getRoot()
	{
		return rootNode;
	}

	double abs(double input)
	{
	    if (input < 0)
            input *= -1;

        return input;
	}

	double round(double num, int places)
    {
        double temp, mult;
        mult = pow(10.0, places);
        temp = floor(num * mult + 0.5);
        temp = temp / mult;
        return temp;
    }

private:

	std::list<IElement*> elementList;			// list of all primitives in the scene
	std::list<IElement*> surfaceList;
	std::vector<Mesh*> meshList;				// list of all meshes in the scene.

	std::string *stringtosend;
	int sindex;

    ElmerVertex* elmVertexList;

    ElmerTriangle* elmTriangleList;
    int tricount;
    int addTriNode;

	int totalnode;
	int totalleafnode;

	unsigned int numvertexused;
	unsigned int numtriangleused;

	BSPTreeNode * rootNode;						// root node of the tree.
    const char *fileName;
};

#endif
