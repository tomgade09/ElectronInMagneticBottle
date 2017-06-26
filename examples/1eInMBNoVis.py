from __future__ import division

import os, sys, inspect, time
a = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
sys.path.append(a + '/../')
sys.path.append(a + '/../vis/')
os.chdir(a)

from classes import *

def main():
    ind = 0                            #Optional - Index (Calculation Counter)
    dt = 1*10**-10                      #Loses sufficient resolution at any faster (larger dt) than 1e-7, Set smaller for greater accuracy/larger for quicker simulation
    
    e1center = [-4,3,0]
    e1vel = [500,4000,4000]            #Need high vparallel to vperpendicular (to B field) ratio to confine
    #e2center = [-4,0,0]
    wccenter = [0,0,0]
    loopaxis = [1,0,0]
    positionHistory = []               #Used in example of what to do with the data

    wireCoils = WireCoilPair(None, wccenter, loopaxis, 1, 1, 5, 5, useC=True) #this creates an object (wire coils) that will be used to calculate the B that shapes the particle's trajectory
    #Arguments are: (window to draw to - use None if not visualizing, center of wire coils, axis the loop lies along, number of turns, current, radius of coil, distance from wccenter)
    
    electron1 = Electron(None, e1center, e1vel) #Mandatory - create a particle to interact with the B Bottle - doesn't have to be an electron

    #electron2 = Electron(None, e2center, e1vel) #Create another particle if you want to model two - note, particles interacting via E field not implemented yet

    Field = Fields(None, dt, [wireCoils], [electron1])#, electron2]) #Mandatory - create a "wrapper" for everything
    #Arguments: (window to draw to, time step, B field producing objects, particles)

    start = time.time() #Optional - for tracking real time elapsed - useful for figuring out how fast your computer does all this stuff
    while (Field.t <= 5.0001e-6): #Mandatory - Everything needs to be nested in a loop - that is unless you only want to do one calculation
    #while ((-10 + wccenter[0]) <= electron1.p[0] <= (10 + wccenter[0])) and ((-10 +  #Optional - Specify a max number of iterations or loop would never end
            #wccenter[1]) <= electron1.p[1] <= (10 + wccenter[1])) and ((-10 + 
            #wccenter[2]) <= electron1.p[2] <= (10 + wccenter[2])):
        
        Field.updateParticleP_V() #Mandatory - updates the particles' Position and Velocity due to the sum of forces on it
        
        ind += 1 #Increment Index #Optional - if you specified an index, better increment it
        
        if ind % 1000 == 0: #Optional but you probably want to know the results of the calculations - otherwise, what's the point?
            #print(ind, time.time() - start)
            #print(ind, electron1.p, electron1.v)
            print('Location: { '+str(electron1.p[0])+', '+str(electron1.p[1])+', '+str(electron1.p[2])+' } | Index: '+str(ind))
            print('Velocity: { '+str(electron1.v[0])+', '+str(electron1.v[1])+', '+str(electron1.v[2])+' } | Time: '+str(Field.t))
            print('')
        
        #Save results here, or do something with the data - whatever you like - an example is below
        positionHistory.append(electron1.p)
        #Now you can save to a csv or other file, print the results on the screen (good luck!), or massage the data some other way
    
    print("Done!") #Optional - if you like your sims talking to you

if __name__ == "__main__": #Python formalism - if you don't like it, remove these two lines and remove the main(): line from above, and remove one tabbed indent per line
    main() #But honestly it's just easier to leave it