#ifndef ELMERTRIANGLE_H_INCLUDED
#define ELMERTRIANGLE_H_INCLUDED

#include <sstream>

class ElmerTriangle
{
    public:
        ElmerTriangle()
        {
            positions = Vector3D(0,0,0);
            id = 0;
        }

        ElmerTriangle(Vector3D _pos, int _id)
        {
            positions = _pos;
            id = _id;
        }

        ~ElmerTriangle()
        {
            procs.clear();
        }

        Vector3D getPositions()
        {
            return positions;
        }

        void setPositions(Vector3D pos)
        {
            positions = pos;
        }

        int getid()
        {
            return id;
        }

        void setid(int _id)
        {
            id = _id;
        }

        int getbid()        // boundary Id set and get.
        {
            return boundId;
        }

        void setbid(int _id)
        {
            boundId = _id;
        }

        int getDIN()
        {
            return domin;
        }

        void setDIN(int _din)
        {
            domin = _din;
        }

        int getDOUT()
        {
            return domout;
        }

        void setDOUT(int _dout)
        {
            domout = _dout;
        }

        void addProc(int procId)
        {
            procs.push_back(procId);
        }

        bool isShared()
        {
            if(procs.size() > 1)
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        std::string getProcs()
        {
            std::string s;
            std::stringstream myout;

            myout << procs.size();
            std::list<int>::iterator it;

            for ( it=procs.begin() ; it != procs.end(); it++ )
            {
                myout << " " << *it;
            }

            s = myout.str();
            return s;
        }

    private:

        Vector3D positions;
        int id;
        std::list<int> procs;
        int boundId;
        int domin;
        int domout;

};


#endif // ELMERTRIANGLE_H_INCLUDED
