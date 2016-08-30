# MagBottlePy

Developer: Tom Gade

MagBottlePy is a library for defining magnetic bottles and modeling motion of charged particles within those bottles.  VPython is used to visualize motion, if visualization is desired.  Visual code has been developed using OpenGL as well, however it is still a work in progress, and should not be relied upon for most visualization uses.  Visualizing in MinVR2 (UMN Virtual Reality Project) is an exception (and the reason OpenGL is used).

Additionally, some code is written in C to optimize a few computationally intensive functions.  See below for info on building for your platform, as well as developing your own C functions to extend the functionality of MagBottlePy.

## Compatibility
MagBottlePy is compatible with all flavors of Python, OS platforms (Darwin is untested, but should work), and 32/64 bit architectures as long as you don't need to visualize the results.  If visualization is desired, at the moment, there must be a version of VPython 5 or 6 available to Python.  This means Python 3 is out on all platforms, as well as VPython visualization on Linux (to the knowledge of the author at the time of writing), unfortunately.  Perhaps other visualization engines will be examined in the future and will be added.  We'll see.

## Getting Started

#### Download Repository

  ```
  git clone https://github.com/tomgade09/MagBottlePy
  cd MagBottlePy
  ```

Use as is or compile the C optimizations (instructions below).

#### Build C Optimization Library (optional)

*Note: This is only necessary if you desire to utilize the code developed in C (used to speed up the calculations by a little over 2x).  This requires a compatible compiler.  For Windows, CMake (3.6.0-rc4 as of this writing) and Visual Studio (Community 2015) are used.  VS can be obtained for free from Microsoft's website.  Google search both to easily find.  On Linux, gcc is used.*

* Windows (CMake and Visual Studio)

1. Run CMake gui.
2. Click "Browse Source" and point to c directory in MagBottlePy root.
3. Click "Browse Build" and point to c/build directory in MagBottlePy root.  If the folder doesn't exist, create it.
4. Click "Configure".
5. Select the appropriate version of VS ("Visual Studio 14 2015" at the time of writing) and click "Finish".  Ignore warnings.
*Note: Ensure to select the applicable processor architecture.  A library compiled for 64 bit applications will not work in 32 bit Python and vice versa!  If you intend to use a 64 bit version of Python AND have a 64 bit version of Windows installed, click the version followed by "Win64" ("Visual Studio ## YEAR Win64").  Otherwise, click the appropriate version of VS without anything following the year ("Visual Studio ## YEAR").  You can ignore "ARM" unless you are going to run on a mobile device or Raspberry Pi (unlikely).*
6. Click "Generate".
7. Open the build folder and double click the *.sln file.  This file should open in VS.  In the menu bar go to "Build" > "Build Solution".
8. Once complete, right click on "Install" in the Solution Explorer (usually on the left side of the screen.  Click "Build".  DLL file should be installed to the lib folder in the project root.  Python code will point to the library automatically.

* Linux (much easier)

  ```
  # Run from c directory in the project root
  gcc -shared -o WireCoilB.so -fPIC WireCoilB.c
  mv WireCoilB.so ./../lib/
  ```
  
* Mac OS X
Untested

## Classes and functions: A brief documentation

Of course, this documentation will assume a working knowledge of Python.  Any arguments defined as ```argument=something``` are optional.

#### Classes

```
BField(windObj, BObjList=[], PartList=[])
```
*Define a Magnetic Field Object with the below arguments, variables, and methods.  Multiple can be defined in one program, allowing for simultaneous execution of non-interacting fields (could be useful for comparing two fields side by side).*

* Use

  ```
  *someVariable* = BField(windObj, BObjList=[*someList*], PartList=[*someList*])
  ```

* Arguments (__init__)

  1. ```windObj``` is a reference to the object that represents the window you want to draw to.  If using MagBottlePy in a purely computational way (no visualization is desired), set this to None.
  2. ```BObjList``` is a list of Magnetic Field producing objects.  This list is iterated over in the TotalBatP method.
  3. ```PartList``` is a list of particles placed in the magnetic field.  These objects may also produce a magnetic field (such as the included classes: electrons, protons, and positrons).
  
* Returns (__init__)
  *No Returns*

* Internal variables

  1. windObj - see above
  2. BObjList - see above
  3. particleList - see above ```PartList```
  
* Methods

  ```totalBatP(p)```
  
  *Calculate the magnetic field at a point, ```p```, due to all objects in BObjList and particleList.  ```p``` is defined as a tuple or list of three values: [x position, y position, z position].  Returns a list of x, y, and z values of a vector representing the B field in the form [Bx, By, Bz].*
  
  ```drawBlines(self, windObj, p, pupbound=[None,None,None], plobound=[None,None,None], numiter=None, linelength=None, multlng=None)```
  *Calculate and draw B field lines starting at a point ```p``` represented as a list or tuple of three values: [x position, y position, z position].  ```[pup|plo]bound``` represents the x,y,z values of the upper|lower bounds.  ```numiter``` allows for the specification of a number of iterations instead of a bound.  ```linelength``` allows specification of a length of B Field line instead of a bound or number of iterations.  ```multlng``` allows for scaling the B vector by a constant.*

## Extending the functionality of MagBottlePy

#### Adding additional B Field producing elements

#### Speeding up B Field iterative calculations with C
