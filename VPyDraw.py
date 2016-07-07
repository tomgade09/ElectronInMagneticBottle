from visual import *

def drawWindow(wd, ht, centx, centy, centz):
    # Call this first to create a window that other draw functions will draw to.
    windObj = display(title='Electron in Magnetic Bottle', autocenter=0, width=wd, 
        height=ht, center=(centx,centy,centz), exit=0, range=(15,15,15))
    xaxpt=[0,1,2,3,4,5,6,7,8,9,10]
    yaxpt=[0,1,2,3,4,5,6,7,8,9,10]
    zaxpt=[0,1,2,3,4,5,6,7,8,9,10]
    xlbl = label(pos=(10,1,0), text='x')
    ylbl =  label(pos=(1,10,0), text='y')
    zlbl =  label(pos=(0,1,10), text='z')
    for i in range(0,10,1):
        points(pos=(xaxpt[i],0,0), size=5, color=color.cyan)
        points(pos=(0,yaxpt[i],0), size=5, color=color.cyan)
        points(pos=(0,0,zaxpt[i]), size=5, color=color.cyan)
    
    return windObj
    
def drawParticlePic(windObj, po_x, po_y, po_z, intrvl, traillng):
    windObj.select()
    particle = sphere(pos=(po_x,po_y,po_z), radius=0.0000001, color=color.green,
        make_trail=True, trail_type="points", interval=intrvl, retain=traillng)
    #particle.trail_object.color = color.yellow
        
    return particle
    
def updateParticlePic(windObj, partObj, px, py, pz):
    windObj.select()
    partObj.pos = (px,py,pz)

def drawTimeClock(windObj, to):
    windObj.select()
    relclockObj = label(pos=(-6.5, 0, 0), text='t = ' + str(to) + ' s')
    
    return relclockObj

def updateTimeClock(windObj, relclockObj, relt):
    windObj.select()
    relclockObj.text = 't = ' + str(relt) + ' s'
    
def drawLine(windObj, px, py, pz, axis_x, axis_y, axis_z):
    windObj.select()
    arrow(pos=(px, py, pz), axis=(axis_x, axis_y, axis_z), color=color.red,
        headwidth=0.005, headlength=0.001, shaftwidth=0.01)
        
def drawWireCoilPair(windObj, d, R):
    windObj.select()
    ring(pos=(-d,0,0), axis=(1,0,0), radius=R, thickness=0.01)
    ring(pos=(d,0,0), axis=(1,0,0), radius=R, thickness=0.01)
    
def FPSrate(fps):
    rate(fps)