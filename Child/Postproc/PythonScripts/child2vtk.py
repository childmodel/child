#!/usr/env/python

"""
child2vtk.py: Python version of Vincent Godard's child2vtk.f95 utility for
converting CHILD model output to vtk format.

Created GT Oct 2013
"""

import sys    # for command-line arguments
import numpy
from readchild import *

_CHILD_REQUIRED_FILES = {
    'node' : 'nodes',
    'elevation' : 'z',
    'triangle' : 'tri',
    'slope' : 'slp',
    'area' : 'area',
    'discharge' : 'q'
}


_CHILD_OPTIONAL_FILES = {
    'shear_stress' : 'tau',
    'network' : 'net',
    'grain_size' : 'tx',
    'voronoi_area' : 'varea'
}



def get_name_of_child_run():
    
    # If there are at least two command-line arguments, use the second one as
    # the file name ...
    if len(sys.argv)>1:
        run_name = sys.argv[1]
    # ... otherwise, ask for it
    else:
        run_name = u''
        while run_name==u'':
            run_name = raw_input('Enter base name of CHILD run: ')
            
    return run_name
    
    
def close_all_files(child_files):
    
    for file_type in child_files:
        child_files[file_type].close()
    
    
def open_required_files(run_name, child_files):
    
    for file_type in _CHILD_REQUIRED_FILES:
        try:
            file_name = run_name + '.' + _CHILD_REQUIRED_FILES[file_type]
            f = open(file_name, 'r')
        except IOError:
            print 'Unable to open <'+file_name+'>'
            close_all_files(child_files)
            raise
        else:
            child_files[file_type] = f
            
    
def open_optional_files(run_name, child_files):
    
    for file_type in _CHILD_OPTIONAL_FILES:
        try:
            file_name = run_name + '.' + _CHILD_OPTIONAL_FILES[file_type]
            f = open(file_name, 'r')
        except IOError:
            print 'Warning: unable to open <'+file_name+'>'
        else:
            child_files[file_type] = f
            
    
def open_child_files(run_name):
    
    child_files = {}
    
    open_required_files(run_name, child_files)
    
    open_optional_files(run_name, child_files)
    
    return child_files
    
    
#def read_node_data(child_files, child_data):
    
def read_time_and_number(the_file):
    
    try:
        t = float(the_file.readline())
        n = int(the_file.readline())
    except ValueError:
        print 'Problem with file format'
        raise
    else:        
        return t, n
        
        
#def read_gridded_data_for_one_time_slice(the_file, nrows, ncols, datatype=float):
    
 #   if ncols>1:
  #      data = numpy.zeros((nrows, ncols), dtype=datatype)
   # else:
    #    data = numpy.zeros(nrows, dtype=datatype)
            
    #for row in range(0, nrows):
        
     #   if ncols==1:
      #      data[row] = datatype(the_file.readline())
       # else:
        #    row_text_items = the_file.readline().split()
         #   for col in range(0, ncols):
          #      data[row,col] = datatype(row_text_items[col])
                
    #return data
    
    
def read_file_to_string_list(myfile):
    
    data = myfile.read() # Read into one long string
    data = data.split()  # Split into a list of strings
    myfile.close()  # Close the file
    
    return data
    
    
def read_and_write_files(child_files, run_name):
    
    # Create data object (a dictionary) and initialize
    child_data = {}
    child_data['time'] = []  # Time values @ each time slice
    child_data['num_nodes'] = []  # Number of nodes @ "
    child_data['num_active_nodes'] = [] # Number of active nodes @ "
    
    # Read contents of files into lists of strings
    nrd = read_file_to_string_list(child_files['node'])
    zrd = read_file_to_string_list(child_files['elevation'])
    trd = read_file_to_string_list(child_files['triangle'])
    srd = read_file_to_string_list(child_files['slope'])
    ard = read_file_to_string_list(child_files['area'])
    qrd = read_file_to_string_list(child_files['discharge'])
    taurd = read_file_to_string_list(child_files['shear_stress'])
    
    time_slice = 1
 
    while nrd:
    
        # Read node-based data for this time slice
        t = float(nrd.pop(0))
        print 'Time=',t
        n = int(nrd.pop(0))
        t = float(zrd.pop(0))
        nz = int(zrd.pop(0))
        assert (n==nz), 'Unequal number of nodes in node and elevation file'
        srd.pop(0)
        srd.pop(0)
        qrd.pop(0)
        qrd.pop(0)
        taurd.pop(0)
        taurd.pop(0)
        x = numpy.zeros(n)
        y = numpy.zeros(n)
        z = numpy.zeros(n)
        eid = numpy.zeros(n, dtype=int)
        b = numpy.zeros(n, dtype=int)
        s = numpy.zeros(n)
        q = numpy.zeros(n)
        tau = numpy.zeros(n)
        for i in range(n):
            x[i] = float(nrd.pop(0))
            y[i] = float(nrd.pop(0))
            eid[i] = int(nrd.pop(0))
            b[i] = int(nrd.pop(0))
            z[i] = float(zrd.pop(0))
            s[i] = float(srd.pop(0))
            q[i] = float(qrd.pop(0))
            tau[i] = float(taurd.pop(0))
            #print x[i], y[i], z[i], eid[i], b[i], s[i], q[i], tau[i]
            
        # Optional data
        if layers_exist(run_name):
            layer_data = creadlayers(run_name, time_slice)
            reg = calc_reg_thickness(layer_data, n)
            tot_lay_thick = calc_total_layer_thickness(layer_data, n)
            if time_slice > 1: 
                dep_ero_rate = (prev_tot_lay_thick-tot_lay_thick)/(t-prev_t)
                tot_ero = starting_lay_thick-tot_lay_thick
            else:
                dep_ero_rate = numpy.zeros(n)
                tot_ero = numpy.zeros(n)
                starting_lay_thick = numpy.copy(tot_lay_thick)
            prev_tot_lay_thick = numpy.copy(tot_lay_thick)
            prev_t = t
        else:
            reg = None
            dep_ero_rate = None
    
        # Extract triangle data
        t = float(trd.pop(0))
        nt = int(trd.pop(0))
        print nt,'triangles'
        tri = numpy.zeros((nt, 3), dtype=int)
        for i in range(nt):  # for each triangle
            for j in range(3):  # get the 3 node IDs
                tri[i,j] = trd.pop(0)
            for j in range(6):  # get and discard the other 6 data elements
                trd.pop(0)
                
        # Extract drainage area data (and any other active-node-based data)
        t = float(ard.pop(0))
        na = int(ard.pop(0))
        print na,'active nodes'
        area = numpy.zeros(n)  # still need to have n entries, but last will be zeros
        for i in range(na):
            area[i] = ard.pop(0)
        
        # Create the VTK output file for this time slice
        time_string = str(time_slice).zfill(4) #time_slice_to_string(time_slice)
        vtkfile = open(run_name+time_string+'.vtk', 'w')
        
        # Write VTK header
        vtkfile.write('# vtk DataFile Version 3.0\n')
        vtkfile.write('CHILD\n')
        vtkfile.write('ASCII\n')
        vtkfile.write('DATASET UNSTRUCTURED_GRID\n')
        
        # Write node x,y,z data
        vtkfile.write('POINTS {0:10d} float\n'.format(n))
        for i in range(n):
            line_out = '{0:15.5f}{1:15.5f}{2:15.5f}\n'.format(x[i], y[i], z[i])
            vtkfile.write(line_out)
            
        # Write triangle data
        vtkfile.write('CELLS{0:10d}{1:10d}\n'.format(nt, 4*nt))
        for i in range(nt):
            line_out = '         3{0:10d}{1:10d}{2:10d}\n'.format(tri[i,0], tri[i,1], tri[i,2])
            vtkfile.write(line_out)
            
        # Write cell type info
        vtkfile.write('CELL_TYPES {0:10d}\n'.format(nt))
        for i in range(nt):
            vtkfile.write(' 5\n')
        
        vtkfile.write('POINT_DATA{0:10d}\n'.format(n))
            
        # Elevation (altitude)
        vtkfile.write('SCALARS Altitude float 1\n')
        vtkfile.write('LOOKUP_TABLE default\n')
        for i in range(n):
            vtkfile.write(str(z[i])+'\n')
        
        # Shear stress
        vtkfile.write('SCALARS Shear_stress float 1\n')
        vtkfile.write('LOOKUP_TABLE default\n')
        for i in range(n):
            vtkfile.write(str(tau[i])+'\n')
        
        # Discharge
        vtkfile.write('SCALARS Discharge float 1\n')
        vtkfile.write('LOOKUP_TABLE default\n')
        for i in range(n):
            vtkfile.write(str(q[i])+'\n')
        
        # Slope
        vtkfile.write('SCALARS Slope float 1\n')
        vtkfile.write('LOOKUP_TABLE default\n')
        for i in range(n):
            vtkfile.write(str(s[i])+'\n')
        
        # Regolith thickness
        if reg is not None:
            vtkfile.write('SCALARS Regolith float 1\n')
            vtkfile.write('LOOKUP_TABLE default\n')
            for i in range(n):
                vtkfile.write(str(reg[i])+'\n')
        
        # Deposition/erosion rate
        if dep_ero_rate is not None:
            vtkfile.write('SCALARS Depo-Ero float 1\n')
            vtkfile.write('LOOKUP_TABLE default\n')
            for i in range(n):
                vtkfile.write(str(dep_ero_rate[i])+'\n')
        
         # Total erosion depth since start of run
        if tot_ero is not None:
            vtkfile.write('SCALARS Total_erosion_depth float 1\n')
            vtkfile.write('LOOKUP_TABLE default\n')
            for i in range(n):
                vtkfile.write(str(tot_ero[i])+'\n')
        
       # Drainage area
        vtkfile.write('SCALARS Drainage_area float 1\n')
        vtkfile.write('LOOKUP_TABLE default\n')
        for i in range(n):
            vtkfile.write(str(area[i])+'\n')
        
        vtkfile.close()
        
        time_slice += 1
        
        
def main():
    
    run_name = get_name_of_child_run()
    
    child_files = open_child_files(run_name)
    
    read_and_write_files(child_files, run_name)  
      
    
if __name__=='__main__':
    main()
