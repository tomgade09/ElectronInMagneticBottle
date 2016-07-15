from __future__ import division

import os, sys, inspect
a = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
sys.path.append(a + '/../')
sys.path.append(a + '/../vis/')
sys.path.append(a + '/../../test/')

import transformations as t
from classes import *

x = [[0,0,0],[1,0,0]]
y = [[0,0,0],[0,1,0]]

w = t.superimposition_matrix(x, y)

print w