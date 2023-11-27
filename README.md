Nurex
=======
Nuclear Reaction Cross Section calculator.
Nurex is a library for Glauber Model nuclear reaction caclulation. The library is written in C++ and Python binding are available.


Installation
============
Get the copy and type:
~~~~
> mkdir build
> cd build
> cmake ../
> make
~~~~

cmake options
-------------
compile options, enable or disable with cmake:
> cmake ../ -D[OPTION]

available options:
  * THREADS - enable/disable multi-threading support, default ON
  * PYTHON_MODULE - enable/disable building of the python bindigs, pybind11 library required to build the module, default OFF
  * TESTS - build tests, default OFF
  * EXAMPLES - build examples, default OFF
  * APPS - build nurex applications, ie command line nurex app  

ie:
> cmake -DPYTHON_MODULE=ON ../


after the compilation the libraries and headers must be either installed system-wide by make install. 
(The default install path can be change, ie: cmake -DCMAKE_INSTALL_PREFIX=/opt/nurex)

Or alternatively the library can be installed in any directory and variables PATH, LD_LIBRARY_PATH and later include directories must be set accordingly.
This can be done sourcing the init.sh file, which is generated in the build directory:
> source init.sh


Manual
========
[Tutorial](docs/manual.md)

How To
======
[How to](docs/howto.md)



PyNurex
=======
Nuclear Reaction Cross Section calculator.
Nurex is a library for Glauber Model nuclear reaction caclulation. 
The library is written in C++ with Python bindings.
Note: This is work in progress. At the moment only basic documentation is available here. The rest will come.

#### Web Interface
WebUI to this libary is running here: https://www.isotopea.com/nurex/

#### Physics Manual
https://www.isotopea.com/nurex/nurex.pdf


	
Build and Installation
======================

Requirements python version __3.7-3.11__, OS: __Linux__ and __Windows__ via pip:
```
pip install pynurex
```

For basic usage see [Python](docs/python.md) and [JSON format](docs/json_format.md) documentation.
