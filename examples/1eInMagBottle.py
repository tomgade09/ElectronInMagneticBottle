from __future__ import division

import os, sys, inspect, time
a = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
sys.path.append(a + '/../')
sys.path.append(a + '/../vis/')
os.chdir(a)

from VPyDraw import *

def main():
    ind = 0                            #Optional - Index (Calculation Counter)
    dt = 1*10**-7                      #Loses sufficient resolution at any faster (larger dt) than 1e-7, Set smaller for greater accuracy/larger for quicker simulation
    
    e1center = [-4,3,0]
    e1vel = [500,4000,4000]            #Need high vparallel to vperpendicular (to B field) ratio to confine
    #e2center = [-4,0,0]
    wccenter = [0,0,0]
    loopaxis = [1,0,0]

    windObj1 = drawWindow(1920, 1080, e1center) #Mandatory for visualization - create a window Object, arguments (width pixels, height pixels, where display is initially centered)
    relclockObj1 = drawTimeClock(windObj1, [-6.5,0,0], 0) #Optional - if you want to see how much time has elapsed, arguments are (windObj from above, location, initial time)

    wireCoils = WireCoilPair(windObj1, wccenter, loopaxis, 1, 1, 5, 5, useC=True) #this creates an object (wire coils) that will be used to calculate the B that shapes the particle's trajectory
    #Arguments are: (window to draw to - use None if not visualizing, center of wire coils, axis the loop lies along, number of turns, current, radius of coil, distance from wccenter)
    wireCoils.initDraw() #Mandatory if you want to visualize
    
    electron1 = Electron(windObj1, e1center, e1vel) #Mandatory - create a particle to interact with the B Bottle - doesn't have to be an electron
    electron1.initDraw(1, 5000) #Mandatory if visualizing - arguments: (interval - draw every ___ point, trail length - keep a history of ___ points on the screen)
    points(pos=e1center, size=5, color=color.red) #Optional - draw a point at electron1 start point for reference

    #electron2 = Electron(windObj1, e2center, e1vel) #Create another particle if you want to model two - note, particles interacting via E field not implemented yet
    #electron2.initDraw(1, 5000)
    #points(pos=e2center, size=5, color=color.red)

    Field = Fields(windObj1, dt, [wireCoils], [electron1])#, electron2]) #Mandatory - create a "wrapper" for everything
    #Arguments: (window to draw to, time step, B field producing objects, particles)
    
    for i in [[-5,3.5,0],[-5,0,3.5],[-5,-3.5,0],[-5,0,-3.5]]: #Optional - draws B Field Lines starting at four points listed to the left
        j = rotateVector(i,wireCoils.axiscf_theta,wireCoils.axiscf_phi) + wireCoils.Cpair
        Field.drawBlines(j, numiter=5750, multlng=500)
    
    wireCoils.drawBfromGC(electron1, Field) #Optional - draw the B line extending from the particle's guiding center to see how well the particle follows it

    start = time.time() #Optional - for tracking real time elapsed - useful for figuring out how fast your computer does all this stuff
    while ind < 2500000: #Mandatory - Everything needs to be nested in a loop - that is unless you only want to do one calculation
    #while ((-10 + wccenter[0]) <= electron1.p[0] <= (10 + wccenter[0])) and ((-10 +  #Optional - Specify a max number of iterations or loop would never end
            #wccenter[1]) <= electron1.p[1] <= (10 + wccenter[1])) and ((-10 + 
            #wccenter[2]) <= electron1.p[2] <= (10 + wccenter[2])):
        FPSrate(10000) #Mandatory for visualization - set frames per second rate for the visualization engine.  Set some ridiculous number like the example to remove the "limit" on how fast the sim runs.
        
        Field.updateParticleP_V() #Mandatory - updates the particles' Position and Velocity due to the sum of forces on it
        electron1.updDraw() #Mandatory for visualization - update the visualization
        #electron2.updDraw()
        
        ind += 1 #Increment Index #Optional - if you specified an index, better increment it
        
        updateTimeClock(windObj1, relclockObj1, Field.t) #Optional - if you created a time clock, better update it
        pauseOnKey('\n', windObj1) #Pause on enter key #Optional - allows you to pause the simulation
        if ind % 1000 == 0: #Optional - updates the center of the visualization so you can follow the particle
            windObj1.center = (electron1.p[0], electron1.p[1], electron1.p[2])
            #print(ind, time.time() - start)
    
    print("Done!") #Optional - if you like your sims talking to you
    while True: #Optional but recommended.  Done with calculations - don't close window, leave it open to continue looking at it
        FPSrate(30)

if __name__ == "__main__": #Python formalism - if you don't like it, remove these two lines and remove the main(): line from above, and remove one tabbed indent per line
    main() #But honestly it's just easier to leave it