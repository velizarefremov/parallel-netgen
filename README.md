# parallel-netgen
Parallel Mesh Generator using NETGEN and C++/MPI.

To Run the Program.

* First install the following preliminaries.

  a) sudo apt-get install build-essential tcl8.5 tcl8.5-dev tk8.5 tk8.5-dev tix tix-dev libtogl1 libtogl-dev glutg3 glutg3-dev libxmu-dev liblapack-dev

* Download the source files for NETGEN from http://sourceforge.net/projects/netgenmesher/files/ and extract it as follows
  
  a) tar xfs netgen-4.9.11.tar.gz
  
  b) cd netgen-4.9.11

* Replace nglib.cpp and nglib.h files with the ones contained in ng.zip

* Install it as follows.
  
  a) ./configure --with-tcl=/usr/lib/tcl8.5/ --with-tk=/usr/lib/tk8.5/ --withtogl=/usr/lib/Togl1.7/
  
  b) Make
  
  c) Sudo make install

* export LD_LIBRARY_PATH=/usr/lib/Togl1.7:/opt/netgen/lib

* download Parmetis version 3.2 from http://glaros.dtc.umn.edu/gkhome/metis/parmetis/download and install it using make and make install

* unzip gl.zip and locate the TreeBuilder.h main.cpp and MeshMig.h. Change the netgen path's inside these files.

* Open parmetisbin.h and change the parmetis paths according to your installation.

* compile the program using g++.

* call the program as follows: “mpirun -n 16 a.exe cube.geo 3 100000 1 1 1” where the parameters are
 
  a) a.exe – the compiled program
  
  b) filename.geo – Name of the input geo file.
  
  c) 3 – Method name 1,2 or 3
  
  d) The following parameters are for method 2 and 3
  
  e) 100000 – Min size of the starting coarse mesh.
  
  f) 3 – Number of refinement steps. Increase to get higher # of elements.
  
  g) 1 – Get geomview output. 0 for no. 1 for yes.
  
  h) 1 – Get Elmer output. 0 for no. 1 for yes.
