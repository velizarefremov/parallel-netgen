#ifndef ELMERVERTEX_H_INCLUDED
#define ELMERVERTEX_H_INCLUDED

#include <sstream>

class ElmerVertex
{
    public:
        ElmerVertex()
        {
            position = Vector3D(0,0,0);
            id = 0;
        }
        ElmerVertex(Vector3D _pos, int _id)
        {
            position = _pos;
            id = _id;
        }
        ~ElmerVertex()
        {
            procs.clear();
        }

        Vector3D getPosition()
        {
            return position;
        }

        void setPosition(Vector3D pos)
        {
            position = pos;
        }

        int getid()
        {
            return id;
        }

        void setid(int _id)
        {
            id = _id;
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

        Vector3D position;
        int id;
        std::list<int> procs;
};


#endif // ELMERVERTEX_H_INCLUDED
