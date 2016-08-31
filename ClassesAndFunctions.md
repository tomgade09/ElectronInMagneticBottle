## Classes and functions: A brief documentation

Of course, this documentation will assume a working knowledge of Python.  Any arguments defined as ```argument=something``` are optional.

### Classes

#### BField
```
BField(windObj, BObjList=[], PartList=[])
```
*Define a Magnetic Field Object with the below arguments, variables, and methods.  Multiple can be defined in one program, allowing for simultaneous execution of non-interacting fields (could be useful for comparing two fields side by side).*

* Use

  ```
  *someVariable* = BField(windObj, BObjList=[*someList*], PartList=[*someList*])
  *someVariable*.totalBatP(p)
  *someVariable*.drawBlines(windObj, p, pupbound=[None,None,None], plobound=[None,None,None], numiter=None, linelength=None, multlng=None)
  ```

* Arguments (\_\_init\_\_)

  1. ```windObj``` is a reference to the object that represents the window you want to draw to.  If using MagBottlePy in a purely computational way (no visualization is desired), set this to None.
  2. ```BObjList``` is a list of Magnetic Field producing objects.  This list is iterated over in the TotalBatP method.
  3. ```PartList``` is a list of particles placed in the magnetic field.  These objects may also produce a magnetic field (such as the included classes: electrons, protons, and positrons).

* Internal variables

  1. windObj - see above
  2. BObjList - see above
  3. particleList - see above ```PartList```
  
* Methods

  ```
  totalBatP(p)
  ```
  
  *Calculate the magnetic field at a point, ```p```, due to all objects in BObjList and particleList.  ```p``` is defined as a tuple or list of three values: [x position, y position, z position].
  Returns a list of x, y, and z values of a vector representing the B field in the form [Bx, By, Bz].*
  
  ```
  drawBlines(windObj, p, pupbound=[None,None,None], plobound=[None,None,None], numiter=None, linelength=None, multlng=None)
  ```
  
  *Calculate and draw B field lines starting at a point ```p``` represented as a list or tuple of three values: [x position, y position, z position].  ```[pup|plo]bound``` represents the x,y,z values of the upper|lower bounds.  ```numiter``` allows for the specification of a number of iterations instead of a bound.  ```linelength``` allows specification of a length of B Field line instead of a bound or number of iterations.  ```multlng``` allows for scaling the B vector by a constant.
  Returns nothing.*
  
#### Particle (and subclasses Electron, Proton, and Positron)
```
Particle(windObj, charge, mass, po, vo)
Electron(windObj, po, vo) #Charge, mass values are filled in with appropriate data
Positron(windObj, po, vo)
Proton(windObj, po, vo)
```
*Define a Particle object with the below arguments, variables, and methods.  Define as many of each as you like.*

* Use

  ```
  *someVariable* = Particle(windObj, charge, mass, [pox,poy,poz], [vox,voy,voz]) OR Electron|Positron|Proton(windObj, [pox,poy,poz], [vox,voy,voz])
  *someVariable*.initDraw(intrvl, traillng, Dcolor=(#R,#G,#B)
  *someVariable*.updDraw()
  *someVariable*.calcBatP(pB)
  *someVariable*.updP(b, dt)
  *someVariable*.foRKvCrossB(BFieldObj, h)
  ```

* Arguments (\_\_init\_\_)

  1. ```windObj``` see above.
  2. ```charge``` is the particle's charge in [C].
  3. ```mass``` is the particle's mass in [kg].
  4. ```po``` is a list of initial position values of the particle in the form [px, py, pz].
  5. ```vo``` is a list of initial position values of the particle in the form [vx, vy, vz].

* Internal variables

  1. wind - see above
  2. q - charge
  3. mass - mass
  4. p - particle position, set initially by po and updated as appropriate methods are called
  5. v - particle velocity, set initially by po and updated as appropriate methods are called
  6. eom - charge over mass (constant)
  7. pic - picture object of the particle (what is drawn on the screen)
  
  
* Methods

  ```
  initDraw(intrvl, traillng, Dcolor=(0,1,0))
  ```
  
  *Initiates an object representing the sphere drawn in ```wind``` and stores it in the variable ```pic```.  Checks to see whether this particle has a ```pic``` object already, and does nothing if it finds one (except print a nagging message).  ```intrvl``` sets the interval of points to draw (if set to 10, will draw every 10th call to ```updDraw()```).  traillng sets how many points to leave on the screen.  Dcolor sets the color with a tuple(only!) of the form (Red Value, Green Value, Blue Value) with values between 0 and 1 (1 being fully "on" and 0 being "off").  Default is green.
  Returns the ```pic``` object.*
  
  ```
  updDraw()
  ```
  
  *Updates the position of ```pic``` in ```windObj``` according to the value of ```p```.
  Returns nothing.*
  
  ```
  calcBatP(pB)
  ```
  
  *Calculates the Magnetic Field at a position ```pB``` due to the particle.
  Returns 3 values: Bx, By, Bz*
  
  ```
  updP(b, dt)
  ```
  
  *Updates ```p``` and ```v``` due to the force exerted on the particle by a Magnetic Field ```b``` of the form [bx,by,bz].  This is iteratively calculated by the change in time ```dt```.
  Returns nothing.*

  ```
  foRKvCrossB(BFieldObj, h)
  ```
  
  *Note: This method is a bit slower than above (calcBatP and updP) and is still experimental.  This code needs (much) more testing before it is relied upon at all.*
  *Uses Fourth Order Runge Kutta to provide (much) higher precision updates to ```p``` and ```v``` due to the magnetic field producing elements of the object of the ```BField class``` represented by the argument ```BFieldObj```.  ```h``` here represents the dt value.
  Returns nothing.*
  
#### WireCoilPair
```
WireCoilPair(windObj, Cpair, axis, N, I, R, d)
```

*To be added*

## Functions

*To be added*