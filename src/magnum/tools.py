# Copyright 2012, 2013 by the Micromagnum authors.
#
# This file is part of MicroMagnum.
#
# MicroMagnum is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MicroMagnum is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MicroMagnum.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import absolute_import, print_function

import os, sys
import numpy as np
import random

def flush():
    import gc
    gc.collect()
    from .magneto import flush
    flush()


def cpu_count():
    try:
        import multiprocessing
        return multiprocessing.cpu_count()
    except:
        from magnum.logger import logger
        logger.warn("Could not find out number of processors")
        return 1

## Fancy colors #########################################

# Enable colors if not windows and console is interactive (and thus hopefully supports ansi escape codes)
# Fixme: Maybe check $TERM variable.
if os.name != "nt" and hasattr(sys.stdout, "isatty") and sys.stdout.isatty():
    def color(c):
        return "\033[" + str(30+c) + "m"
    def nocolor():
        return "\033[0m"
else:
    def color(c):
        return ""
    def nocolor():
        return ""

## Interactive menus ####################################

if sys.version_info < (3, 0):
    getline = raw_input # Python 2.x
else:
    getline = input # Python 3.x

def print_header(header, width):
    pad = (width - len(header)) - 2
    hdr = "="*pad + "[" + header + "]" + "="*pad
    print(hdr)

def interactive_menu(header, text, options):
    print_header(header, 60)
    print(text)
    for idx, opt in enumerate(options):
        print("  %i. %s" % (idx+1, opt))
    for tempNumberToAvoidInfLoops in range(0,10):
        print("Choice: ", end="")
        try:
            ans = int(getline())
            if ans < 1 or ans > len(options)+1:
                raise ValueError()
        except:
            print("Type a number between 1 and %i." % len(options))
            continue
        break
    return ans

## Generate a list of floats from a range ###############

def frange(*args):
    """
    A float range generator. Usage:

    frange(start [,stop [,step]]) -> generator

    Examples:
      list(frange(4.2))            -> [0.0, 1.0, 2.0, 3.0, 4.0]
      list(frange(2.2, 5.6))       -> [2.2, 2.3, 4.3, 5.3]
      list(frange(2.2, 5.6, 0.25)) -> [2.2, 2.45, 2.7, 2.95, 3.2, 3.45, 3.7, 3.95, 4.2, 4.45, 4.7, 4.95, 5.2, 5.45]
    """
    start = 0.0
    step = 1.0

    l = len(args)
    if l == 1:
        end, = args
    elif l == 2:
        start, end = args
    elif l == 3:
        start, end, step = args
        if step == 0.0:
            raise ValueError("step must not be zero")
    else:
        raise TypeError("frange expects 1-3 arguments, got %d" % l)

    v = start
    while True:
        if (step > 0 and v >= end) or (step < 0 and v <= end):
            raise StopIteration
        yield v
        v += step

def plot2Dsimple(M, path = ".", text="2d"):
	import matplotlib.pyplot as plt
	import matplotlib.cm as cmx
	# plot main data
	(nx,ny,nz) = M.mesh.getNumNodes()
	(dx,dy,dz) = M.mesh.getDelta()
	DATA = []
	X = []
	Y = []
	for j in range(ny):
		for i in range(nx):
			X.append((i+0.5)*dx*1e9)
			Y.append((j+0.5)*dy*1e9)
			DATA.append(M.get(i,j,0))
			
	DATA = np.array(DATA)
	X = np.array(X)
	Y = np.array(Y)
	grid = DATA.reshape((ny, nx))		
	
	plt.imshow(grid, extent=(X.min(), X.max(), Y.min(), Y.max()),interpolation='nearest', cmap=cmx.gist_rainbow)
	
	plt.xlabel("x [nm]")
	plt.ylabel("y [nm]")
	filedest = path + "/" + "plot_" + text + ".png"
	plt.savefig(filedest)
	print("saved:", filedest)
	plt.close()

def update_progress(progress, barLength,text="progress"):
    #barLength = 30 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "FRAMEWORK: error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress > 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\r" + text + ": [{0}] {1}% ".format( "#"*block + "-"*(barLength-block), int((progress)*100)) 
    sys.stdout.write(text)
    sys.stdout.flush()
	
def makeGrains2D(K, numcenters, variationPercent, distMin, seed, showProgress=False):
	from scipy.spatial import Voronoi
	from magnum.logger import logger
	
	# determine if a point is inside a given polygon or not
	# Polygon is a list of (x,y) pairs.
	#http://www.ariel.com.au/a/python-point-int-poly.html
	def point_inside_polygon(x,y,poly):
		n = len(poly)
		inside =False
		p1x,p1y = poly[0]
		for i in range(n+1):
			p2x,p2y = poly[i % n]
			if y > min(p1y,p2y):
				if y <= max(p1y,p2y):
					if x <= max(p1x,p2x):
						if p1y != p2y:
							xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
						if p1x == p2x or x <= xinters:
							inside = not inside
			p1x,p1y = p2x,p2y
		return inside

	def getDistance(M,a,b,c,x,y,z):

		(nx,ny,nz) = M.mesh.getNumNodes()
		(dx,dy,dz) = M.mesh.getDelta()
		
		px = M.mesh.periodic_bc[0].find("x") != -1
		py = M.mesh.periodic_bc[0].find("y") != -1
		pz = M.mesh.periodic_bc[0].find("z") != -1
		
		if not (px or py or pz):
			res = (a-x,b-y,c-z)
		elif px:
			res = ((a-x)%(nx*dx),b-y,c-z)
		elif py:
			res = (a-x,(b-y)%(ny*dy),c-z)
		elif pz:
			res = (a-x,b-y,(c-z)%(nz*dz))	
		elif py and py:
			res = ((a-x)%(nx*dx),(b-y)%(ny*dy),c-z)	
		elif px and pz:
			res = ((a-x)%(nx*dx),b-y,(c-z)%(nz*dz))	
		elif py and pz:
			res = (a-x,(b-y)%(ny*dy),(c-z)%(nz*dz))		
		else: 
			print("encountered unknown periodic boundary settings. exiting.")
			exit(2)
		distance = (res[0]**2 + res[1]**2 + res[2]**2)**0.5
				
		return distance
		
		
	logger.info("Creating grains...")
	if seed != "": 	random.seed(seed)
	else:			random.seed(1234568756489)
	
	(nx,ny,nz) = K.mesh.getNumNodes()
	(dx,dy,dz) = K.mesh.getDelta()
	
	centers = []
	idc = 0
	logger.info("  -- Building Centers...")
	while idc < numcenters:
		xx = (random.random()*nx+0.5)*dx
		yy = (random.random()*ny+0.5)*dy
		zz = (random.random()*nz+0.5)*dz
		continueO = False
		for oldcen in centers:
			if getDistance(K, (xx,yy,zz), (oldcen[0],oldcen[1],oldcen[2])) < distMin: continueO = True
		if continueO: continue
		centers.append((xx, yy, zz))
		idc += 1
		if showProgress: update_progress(float(idc)/float(numcenters),20, "building centers")
	if showProgress: print()
	
	points = []
	for cen in centers:
		points.append((cen[0], cen[1]))
	
	vor = Voronoi(points)
	
	regions = []
	for reg in vor.regions:
		entry = []
		for elm in reg:
			if elm < 0: continue
			entry.append(vor.vertices[elm])
		if len(entry) < 3: continue
		if entry != []: regions.append((entry, 1.+variationPercent*(2.*random.random()-1.)))
	
	logger.info("  -- Filling Patterns...")
	for i in range(nx):
		for j in range(ny):
			for k in range(nz):
				x,y,z = (i+0.5)*dx, (j+0.5)*dy, (k+0.5)*dz
				for reg in regions:
					if point_inside_polygon(x,y,reg[0]): 
						K.set(i,j,k, K.get(i,j,k)*reg[1])
						break
		if showProgress: update_progress(float(i)/float(nx),20, "filling patterns")
		#else: 
		#	if int(float(i)/float(nx) * 100.) % 10 == 0: logger.info("     " + str(int(float(i)/float(nx) * 100.)) + "...")
	if showProgress: print()
	return K		
