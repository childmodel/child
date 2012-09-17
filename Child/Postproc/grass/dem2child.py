#! /usr/bin/python

# Written August--September 2012 by Andrew D. Wickert

"""
Copyright 2012 Andrew D. Wickert

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""


import argparse
import sys
import numpy as np
import os

# Parse some command line input to direct the rest of the script, which should
# automate GRASS export

# See http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit() # (2)

parser = MyParser(
  prog = "dem2child.py",
  description = """
  This program uses GRASS GIS 7.0 and Python 2.7 to generate a triangular grid 
  for input to the landscape evolution model CHILD. This grid is written to the  
  output file grassdem.pts.
  
  Sample Usage: ./dem2child.py --dem=topo --resolution=100000 --n=1 --s=1 --w=1 --e=1
  """,
  epilog = """
  Current usage is limited to setting boundary conditions on 
  sides of the whole region, and large cliffs are created at all boundaries. 
  Future work should fix this.
  
  Also is currently limited to only projected coordinate systems; don't know if 
  this will ever change.
  """
)

parser.add_argument('--dem', type=str, help='name of elevation raster')
parser.add_argument('--resolution', type=float, help='output grid resolution (in map units) PROJECTED ONLY')
parser.add_argument('--n', type=int, help='northern boundary: closed (Neumann, 1) or open (Dirichlet, 2)')
parser.add_argument('--s', type=int, help='southern boundary: closed (1) or open (2)')
parser.add_argument('--w', type=int, help='western boundary: closed (1) or open (2)')
parser.add_argument('--e', type=int, help='eastern boundary: closed (1) or open (2)')
#parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")

args = parser.parse_args()

if args.dem and args.resolution and args.n and args.s and args.w and args.e:
  pass
else:
  sys.exit(
  """
  In current version of program, must define:
  dem ---------> elevation data set name in current GRASS location
  resolution --> of CHILD grid
  n, s, w, e --> boundary conditions
  """
  )

try:
  from grass import script as g
  print "---------------------------"
  print "GRASS GIS RUNNING. VERSION:"
  version = g.version()
  for row in version.items():
    for item in row:
      print item,
    print ""
  print "ENVIRONMENT:"
  gisenv = g.gisenv()
  for row in gisenv.items():
    for item in row:
      print item,
    print ""
  print "---------------------------"
except:
  sys.exit("Must be run inside GRASS GIS; tested only with GRASS 7.0")
  
# Create vector map around boundaries and classify them

# Compuational region
print "Obtaining geographical region from GRASS."
reg = g.region()
n = reg['n']
s = reg['s']
w = reg['w']
e = reg['e']
dns = n-s
dew = e-w

dx = args.resolution
dy = args.resolution * np.cos( np.radians(30) )

nsres = dx
ewres = dy

print "Classifying map boundaries..."

if args.n == 1:
  top = n + nsres/2.
else:
  top = n
if args.n == 1:
  bottom = s - nsres/2.
else:
  bottom = s
if args.w == 1:
  left = w - ewres
else:
  left = w - ewres/2.
if args.e == 1:
  right = e + ewres
else:
  right = e + ewres/2.

"""
if args.n == 1:
  top = n
else:
  top = n - nsres/2.
if args.s == 1:
  bottom = s
else:
  bottom = s + nsres/2.
if args.w == 1:
  left = w - ewres/2.
else:
  left = w
if args.e == 1:
  right = e + ewres/2.
else:
  right = e
"""
  
hexns = top-bottom
hexew = right-left

print "Building hex grid nodes..."

ny = int(np.ceil(hexns/args.resolution))
nx_more = int(np.ceil(hexew/args.resolution))
nx_fewer = nx_more - 1
resolution_y = hexns / ny
resolution_x = hexew / nx_more

# Values: loop over rows of constant y with correct x progression
y_values = np.linspace(bottom, top, ny)
x_more_values = np.linspace(left, right, nx_more)
x_fewer_values = np.linspace(left + args.resolution/2., right - args.resolution/2., nx_fewer)

x = []
y = []

# Create x and y
for i in range(ny):
  # Start with a short row, then do a longer one

  if i % 2 == 0:
    x.extend(x_fewer_values)
    y.extend(y_values[i] + 0 * x_fewer_values)
  else:
    x.extend(x_more_values)
    y.extend(y_values[i] + 0 * x_more_values)

# Export x,y list and bring into GRASS
xyarray = np.array([x,y]).transpose()
np.savetxt('.xyarray.txt', xyarray, delimiter='|')
g.run_command('v.in.ascii', input='.xyarray.txt', output='CHILDgrid', overwrite=True, quiet=True)
os.remove('.xyarray.txt') # Clean up after oneself - perhaps this should go in a tmp directory

# For fun, if we want to
# v.voronoi in=CHILDgrid out=CHILDvoronoi
# v.delaunay in=CHILDgrid out=CHILDdelaunay

# Deprecated
## Extrapolate from DEM edges to get elevation values at sides
#g.run_command('g.region', vect='CHILDgrid', flags='p')
#g.run_command('r.fillnulls', input=args.dem, output=args.dem + '_extended_CHILDgrid')

# Create and update attribute table to get z values
g.run_command('v.db.addtable', map='CHILDgrid', columns='x double precision, y double precision, z double precision, b integer', quiet=True)
g.run_command('v.to.db', map='CHILDgrid', option='coor', columns='x, y', quiet=True) # x, y
g.run_command('g.region', vect='CHILDgrid', quiet=True)
g.run_command('v.what.rast', map='CHILDgrid', raster=args.dem, column='z', quiet=True) # z
g.run_command('g.region', rast=args.dem, quiet=True)

# Including boundary types: do open (Dirichlet) boundaries first, and then closed (Neumann) ones
# so corners are all closed
# LATER!
"""
open_boundaries = []
closed_boundaries = []
boundary_names = ['n', 's', 'e', 'w']

i = 0
for boundary in args.n,args.s,args.w,args.e:
  if boundary == 1:
    #closed_boundary_names.append[boundary_names]
    #closed_boundary_names.append[boundary_names]
  elif boundary == 2:
    closed_boundaries.append[boundary]
"""

# Start with 0 everywhere - interior - and then overwrite
# THIS CAUSES BIG CLIFFS OFF SIDES OF MODEL - SO PROBLEMS WITH OPEN BOUNDARIES
g.run_command('v.db.update', map='CHILDgrid', column='b', value=0, quiet=True)
#for boundary in open_boundaries:
g.run_command('v.db.update', map='CHILDgrid', column='b', value=args.n, where="y >= " + str(n), quiet=True)
g.run_command('v.db.update', map='CHILDgrid', column='b', value=args.s, where="y <= " + str(s), quiet=True)
g.run_command('v.db.update', map='CHILDgrid', column='b', value=args.e, where="x >= " + str(e), quiet=True)
g.run_command('v.db.update', map='CHILDgrid', column='b', value=args.w, where="x <= " + str(w), quiet=True)
# d.vect.thematic map=CHILDgrid icon=basic/circle column=b type=point nint=3
# (above) working on some display method

print "***"
print "WRITING GRID..."
# Next, export to array for CHILD
# (could use g.vector_db_select to import into Python if we want to do this in memory)
# Not sure what Greg likes for null values, so I will just use "NULL"
# Could also make equal to closest point with v.distance and v.extract steps
g.run_command('v.db.select', map='CHILDgrid', fs=' ', nv='0', columns='x,y,z,b', file='grassdem.pts', flags='c', overwrite=True)

npoints = xyarray.shape[0] # Will need to change if restricting to a single basin
  

f = open('grassdem.pts','r')
temp = f.read()
f.close()

f = open('grassdem.pts', 'w')
f.write(str(npoints)+'\n') # might not work with Windows
f.write(temp)
f.close()

print "***"
print "DONE. File 'grassdem.pts' should be availble in this directory for CHILD input."
print "Thank you for using grass2child.py :)"
print "***"

