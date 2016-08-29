# MagBottlePy

Developer: Tom Gade

MagBottlePy is a library for defining magnetic bottles and modeling motion of charged particles within those bottles.  VPython is used to visualize motion, if visualization is desired.  Visual code has been developed using OpenGL as well, however is still a work in progress, and should not be relied upon for most visualization uses.  Visualizing in MinVR2 (UMN Virtual Reality Project) is an exception (and the reason OpenGL is used).

Additionally, some code is written in C to optimize a few computationally intensive functions.  See below for info on building for your platform, as well as developing your own C functions and integrating them into your implementation of MagBottlePy.

## Compatibility:
MagBottlePy is compatible with all flavors of Python, OS platforms (Darwin is untested, but should work), and 32/64 bit architectures as long as you don't need to visualize the results.  If visualization is desired, at the moment, there must be a version of VPython 5 or 6 available to Python.  This means Python 3 is out on all platforms, as well as