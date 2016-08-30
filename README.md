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
*Note: Ensure to select the applicable processor architecture.  A library compiled for 64 bit applications will not work in 32 bit Python!  If you intend to use a 64 bit version of Python AND have a 64 bit version of Windows installed, click the version followed by "Win64" ("Visual Studio ## YEAR Win64").  Otherwise, click the appropriate version of VS without anything following the year ("Visual Studio ## YEAR").  You can ignore "ARM" unless you are going to run on a mobile device or Raspberry Pi (unlikely).*
6. Click "Generate".
7. Open the build folder and double click the *.sln file.  This file should open in VS.  In the menu bar go to "Build" > "Build Solution".
8. Once complete, right click on "Install" in the Solution Explorer (usually on the left side of the screen.  Click "Build".  DLL file should be installed to the lib folder in the project root.  Python code will point to the library automatically.

* Linux (much easier)

1. Change to c directory in the project root from the command line.  Type the following.

  ```
  gcc -shared -o WireCoilB.so -fPIC WireCoilB.c
  mv WireCoilB.so ./../lib/
  ```

.so Library file should be moved to the appropriate directory with these commands.  Python code will point to the library automatically.

## Classes and functions: A brief documentation

## Extending the functionality of MagBottlePy

#### Adding additional B Field producing elements

#### Speeding up B Field iterative calculations with C
