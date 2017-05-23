from classes import *
from visual import *

#def setupDisplay(wd=1920, ht=1080, center=[0,0,0], axistype="points"):
# Could do this to make setup easier...put everything in a wrapper function/class
#def updateDisplay():

def pauseOnKey(keystr, windObj):
    if windObj.kb.keys:
        key = windObj.kb.getkey()
        if key == keystr:
            print("Press Enter to resume.")
            windObj.waitfor('\n keydown')
            print("Resumed.")
            windObj.kb.queue = []

def drawWindow(wd, ht, cent, axistype="points"):
    centx = cent[0]
    centy = cent[1]
    centz = cent[2]
    # Call this first to create a window that other draw functions will draw to.
    windObj = display(title='Electron in Magnetic Bottle', autocenter=0, width=wd, 
        height=ht, center=(centx,centy,centz), exit=0, range=(15,15,15))
    if axistype == "points":
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
    elif axistype == "lines":
        for i in [[10,0,0], [0,10,0], [0,0,10]]:
            drawLine(windObj, [0,0,0], i, shaftwd=0.05, headwd=0.1, headln=0.15,
                col=color.cyan)
    
    return windObj
    
def drawParticlePic(windObj, po, intrvl, traillng, col):
    windObj.select()
    particle = sphere(pos=(po[0],po[1],po[2]), radius=0.0000001, color=col,
        make_trail=True, trail_type="points", interval=intrvl, retain=traillng)
    
    return particle
    
def updateParticlePic(windObj, partObj, p):
    windObj.select()
    partObj.pos = (p[0],p[1],p[2])

def drawTimeClock(windObj, po, to):
    windObj.select()
    relclockObj = label(pos=po, text='t = ' + str(to) + ' s')
    
    return relclockObj

def updateTimeClock(windObj, relclockObj, relt):
    windObj.select()
    relclockObj.text = 't = ' + str(relt) + ' s'
    
def drawLine(windObj, p, laxis, headwd=0.005, headln=0.001, shaftwd=0.01,
    col=color.red):
    px = p[0]
    py = p[1]
    pz = p[2]
    axis_x = laxis[0]
    axis_y = laxis[1]
    axis_z = laxis[2]
    windObj.select()
    arrow(pos=(px, py, pz), axis=(axis_x, axis_y, axis_z), color=col, 
        headwidth=headwd, headlength=headln, shaftwidth=shaftwd)
        
def drawWireCoilPair(windObj, C, norm, cntlf, cntrt, R):
    windObj.select()
    objlf = ring(pos=cntlf, axis=norm, radius=R, thickness=0.01)
    objrt = ring(pos=cntrt, axis=norm, radius=R, thickness=0.01)
    
    return [objlf, objrt]
    
def FPSrate(fps):
    rate(fps)
    
def displayFlip(): #Placeholder
    return
    
#Needs work
def keyInput(evt):
    #global key
    #key = evt.key
    return

def inputHandler():
    scene.bind('keydown', keyInput)
    
#Methods for "Particle" Class
def PinitDraw(self, intrvl, traillng, Dcolor=(0,1,0)):
    """Draw a point object representing the particle in the object specified by self.wind.  intrvl represents how often to 'draw' a point.  traillng represents how long of a 'trail' to leave behind the current position of the particle."""
    if self.pic is not None:
        print("Pic has already been initialized. Use updDraw to change the position.")
        return
    self.pic = drawParticlePic(self.wind, self.p, intrvl, traillng, Dcolor)
    return self.pic
        
def PupdDraw(self):
    """Update the location of the point object drawn with initDraw.  Obviously, it can't be updated if it hasn't been initialized.  Use the initDraw function first."""
    if self.pic is None:
        print("Pic has not been initialized.  Use initDraw to create a picture of the particle first.")
        return
    updateParticlePic(self.wind, self.pic, self.p)

Particle.initDraw = PinitDraw
Particle.updDraw = PupdDraw

#Methods for "WireCoilPair" Class
def WCinitDraw(self):
    """Draw the pair of Wire Coils."""
    self.pic = drawWireCoilPair(self.wind, self.Cpair, self.axis, self.Cleft,
        self.Cright, self.R)

WireCoilPair.initDraw = WCinitDraw
        
#Methods for "BField" Class
def drawBlines(self, p, pupbound=[None,None,None], plobound=[None,None,None], 
    numiter=None, linelength=None, multlng=None):
    """Draw B field lines starting at po and ending at ####."""
    loopind = 0
    BoundBool = True
    pts=[]
    while BoundBool:
        for i in [pupbound, plobound]:
            #ToDo Write code to check lower bound - maybe negative both sides
            if i[0] is not None:
                BoundBool = BoundBool and p[0] <= i[0]
            if i[1] is not None:
                BoundBool = BoundBool and p[1] <= i[1]
            if i[2] is not None:
                BoundBool = BoundBool and p[2] <= i[2]
        if numiter is not None:
            loopind += 1
            BoundBool = BoundBool and loopind < numiter
            
        bx, by, bz = self.totalBatP(p)
            
        if linelength is not None:
            Blen, Bth, Bphi = cartesianToSpherical([bx,by,bz])
            bx, by, bz = sphericalToCartesian(linelength, Bth, Bphi)
        elif multlng is not None:
            bx *= multlng; by *= multlng; bz *= multlng
            
        #drawLine(self.windObj, p, [bx,by,bz])
        p[0] += bx; p[1] += by; p[2] += bz
        pts.append((p[0],p[1],p[2]))
    #print(pts)
    self.windObj.select()
    curve(pos=pts, color=color.red)
        
Fields.drawBlines = drawBlines