#! /usr/env/python
"""
readchild.py: contains routines to read data from a CHILD run.

created oct 2013 gt
"""

import numpy


class Layer():
    
    def __init__(self, thick, ctime, rtime, etime, erody, is_regolith, dgrade=[1.0] ):
        
        self.thickness = thick
        self.creation_time = ctime
        self.activation_time = rtime
        self.exposure_time = etime
        self.erodibility = erody
        self.is_regolith = is_regolith 
        self.grain_size_fractions = dgrade
        
        
    def report(self):
        
        if self.is_regolith:
            print 'LAYER INFO: layer is made of REGOLITH'
        else:
            print 'LAYER INFO: layer is made of ROCK'
        print 'Thickness:',self.thickness
        print 'Creation time:',self.creation_time
        print 'Reworking time:',self.activation_time
        print 'Exposure time:',self.exposure_time
        print 'Erodibility:',self.erodibility
        print 'Grain size fractions:',self.grain_size_fractions
        
        
def layers_exist(run_name):
    
    #print 'layers_exist here ...'
    try:
        open(run_name+'.lay0', 'r')
    except IOError:
        #print 'returning False'
        return False
    else:
        #print 'returning True'
        return True
        

def creadlayers(basenm, ts, numg=1, format_version=0):
    
    layerfile = open(basenm+'.lay'+str(ts-1),'r')
    layinfo = layerfile.read().split()
    layerfile.close()
    
    maxlayers = 600
    
    model_time = float(layinfo.pop(0))
    print 'creadlayers: reading data for time',model_time
    n = int(layinfo.pop(0))
    
    num_layers = numpy.zeros(n, dtype=int)
    
    layer_data = []
    
    # Read info on layer(s) for each node
    for i in range(n):
        
        # Create an empty layer list
        layer_list = []
        
        # Read the number of layers at this node
        num_layers[i] = int(layinfo.pop(0))
        #print 'node',i,'has',num_layers[i],'layers'
        assert (num_layers[i]>0), 'Node '+i+' has only '+num_layers[i]+' layers'
        assert (num_layers[i]<maxlayers), 'Node '+i+' claims to have '+num_layers[i]+' layers'
    
        # For each layer, read its properties
        for k in range(num_layers[i]):
            
            # Read properties for the current layer
            ctime = float(layinfo.pop(0))
            rtime = float(layinfo.pop(0))
            etime = float(layinfo.pop(0))
            thick = float(layinfo.pop(0))
            erody = float(layinfo.pop(0))
            sedrock = int(layinfo.pop(0))
            is_regolith = (sedrock==1)
            dgrade = []
            for d in range(numg):
                dgrade.append(float(layinfo.pop(0))/thick)
            
            # Put them in a layer object
            lay = Layer(thick, ctime, rtime, etime, erody, is_regolith, dgrade)
            
            # Put the layer object on the list for this node
            layer_list.append(lay)
            
        # Add the layer list to the layer_data
        layer_data.append(layer_list)
        
    return layer_data
    
    
def calc_total_layer_thickness(layer_data, num_nodes):
    
    thickness = numpy.zeros(num_nodes)
    
    i=0
    for layerlist in layer_data:
        for layer in layerlist:
            thickness[i] += layer.thickness
        i+=1
        
    return thickness
    
    
def calc_total_layer_thickness_from_file(run_name, time_slice, numgsize=1,
                                         num_nodes=None):
    
    if num_nodes is None:
        num_nodes = find_number_of_nodes(run_name, time_slice)
    layer_data = creadlayers(run_name, time_slice, num_rain_sizes=numgsize)
    lay_thick = calc_total_layer_thickness(layer_data, num_nodes)
    
    return lay_thick
    
    
def find_number_of_nodes(run_name, time_slice):
    
    f = open(run_name+'.z', 'r')
    for i in range(time_slice):
        t = f.readline()
        n = int(f.readline())
        if i<(time_slice-1):
            for j in range(n):
                f.readline()
    f.close()

    return n
        
    
def calc_reg_thickness(layer_data, num_nodes):
    
    regolith_thickness = numpy.zeros(num_nodes)
    
    node_num = 0
    for layer_list in layer_data:
        for layer in layer_list:
            if layer.is_regolith:
                regolith_thickness[node_num] += layer.thickness
            else:
                break
        node_num += 1
            
    return regolith_thickness
    
        
def calc_reg_thickness_from_file(run_name, time_slice, num_grain_size=1):
    
    layer_data = creadlayers(run_name, time_slice, numg=num_grain_size)
    num_nodes = find_number_of_nodes(run_name, time_slice)
    reg_thick = calc_reg_thickness(layer_data, num_nodes)
    return reg_thick
    
    
# Here's some temporary testing stuff

#ld = creadlayers('/Users/gtucker/Runs/Test/readlaytest', 21)

#for node in ld:
#    print '\n**** LAYERS AT NODE ****'
#    for lay in node:
#        lay.report()