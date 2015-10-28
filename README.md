Download the code from
http://code.google.com/p/parallel-netgen/source/

Follow the instruction on the following page. http://sourceforge.net/apps/mediawiki/netgen-mesher/index.php?title=Installing_on_Ubuntu But replace the nglib.cpp and nglib.h files in the netgen directory with the ones in "NG" folder. Compile with the replaced files.

Download Parmetis from http://glaros.dtc.umn.edu/gkhome/metis/parmetis/download and install it using make - make install

Locate the TreeBuilder.h main.cpp and MeshMig.h files in "GL" folder. Change the netgen path's inside these files according to your installation.

Open parmetisbin.h and change the parmetis paths according to your installation.

compile the program.

call the program as follows.

mpirun -n 16 a.exe cube.geo 2 100000 1 1 1

where after a.exe

filename method minsize recursiondepth elmerout geomout



////

filename.geo // Name of the input geo file.
Method 3 // Method 1,2 or 3

// The following are for method 2 and 3

minstartsize 100000 // Min size of the starting coarse mesh.
refinerepeat 3 // Number of refinement steps. Increase to get higher # of elements.
geomout 1 // Get geomview output. 0 for no. 1 for yes.
elmerout 1 // Get Elmer output. 0 for no. 1 for yes.

**Code is tested with Netgen 4.9.13 and Parmetis 3.2.0**