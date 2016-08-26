from __future__ import division,print_function

import os, sys, inspect
a = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
sys.path.append(a + '/../')
sys.path.append(a + '/../vis/')

from classes import *
from VPyDraw import *

windObj1 = drawWindow(1920, 1080, [0,0,0], axistype="lines")
electron1 = Electron(windObj1, [0.1,0.1,0.1], [1000,1000,1000])
electron1.initDraw(1,0)
B = BField(windObj1)
B.BObjList.append(electron1)

drawLine(windObj1, [0.5,0.5,0.5], [0.5,0.5,0.5], shaftwd=0.05, headwd=0.1, headln=0.15,
    col=color.white)
label(pos=[0.5,1,1],text='v')

FPSrate(30)
B.drawBlines(windObj1, [0.11,0.1,0.1], numiter=100, multlng=1e16)
B.drawBlines(windObj1, [0.09,0.1,0.1], numiter=100, multlng=1e16)
B.drawBlines(windObj1, [0.1,0.11,0.1], numiter=100, multlng=1e16)
B.drawBlines(windObj1, [0.1,0.09,0.1], numiter=100, multlng=1e16)
B.drawBlines(windObj1, [0.1,0.1,0.11], numiter=100, multlng=1e16)
B.drawBlines(windObj1, [0.1,0.1,0.09], numiter=100, multlng=1e16)

print("done")

while True:
    FPSrate(30)