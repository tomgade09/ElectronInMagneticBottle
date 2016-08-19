from __future__ import division

import os, sys, inspect
a = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
os.chdir(a)

import csv
import numpy as np
from math import sqrt

data = []
base = []
subtract = []
errdist = []

for j in range(-9,-14,-1):
    for i in [5,2,1]:
        filedata = []
        resultfile = open(str(i*(10**j)) + '.csv','a+')
        rd = csv.reader(resultfile)
        for row in rd:
            tmp = []
            for elem in row:
                tmp.append(float(elem))
            filedata.append(tmp)
        data.append(filedata)

filedata = []
resultfile = open('5e-14.csv','a+')
rd = csv.reader(resultfile)
for row in rd:
    tmp = []
    for elem in row:
        tmp.append(float(elem))
    filedata.append(tmp)
base.append(filedata)

for k in range(15):
    subtract.append(np.array(data[k])-np.array(base[0]))
    subtract.append([])

ind = 0

for l in subtract:
    tmp = []
    for rows in l:
        tmp.append([rows[0],sqrt(rows[1]**2+rows[2]**2+rows[3]**2)])
    errdist.append(tmp)
    
resultfile = open('subtract.csv','w+')
wr = csv.writer(resultfile)
for rows in subtract:
    wr.writerows(rows)
    wr.writerow([])

resultfile = open('errdist.csv','w+')
wr = csv.writer(resultfile)
for rows in errdist:
    wr.writerows(rows)
    wr.writerow([])