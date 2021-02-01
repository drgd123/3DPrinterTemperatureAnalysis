#!/usr/bin/env python3
import time
from Vertex import Vertex

__author__ = 'Greg Dreifus'

in_file_path = './toolPath_128.txt'
out_file_path = './toolPath_128-40percentEval.gcode'
layers = []

# Define Print Settings.
num_layers = 1
epew = 0.50 #external perimeters extrusion width
pew = 0.72 #perimeters extrusion width
iew = 0.72 #infill extrusion width
tiew = 0.72 #top infill extrusion width

FlayerZ = 7800.000 #Feedrate to next layer?
Flayer1 = 2400.00000 #Feedrate to purge?

E_purge = -2.00000 #purge extrusion?
distance = 0 #the distance traveled at each time step
E_ratio = 1.06 / 40 #the ratio of material extruded to distance traveled in mm? 0.0265
E_val = distance * E_ratio * 0.4 #eparameter


layer_thickness = 2.54  # incremental increase in height per Layer
tableHeight = 0
# z_lift = 2
feed = 4000
start_height = 0  # a temp z_lift to be updated for absolute coordinates
# start_x = 1071.428571 for complex geometry
scale_down_factor = 1#100
start_x = 622.0441539999999 / scale_down_factor
start_y = 700 / scale_down_factor

# Read in vertex path file formatted as:
# xx.xxx yy.yyy zz.zzz print
# where print is 1 if printing when leaving node and 0 otherwise.
#with open(in_file_path) as f, open('check_file.txt', 'w') as check_file:
with open(in_file_path) as f:
    previous_z_value = -1.0
    is_airtime = True # Initial move to the first vertex is an airtime move.
    reading_interior = False # Initially reading boundary vertices.

    for line in f:
        if line[0] == '#': # Skip comment lines.
            continue
        
        line_elements = line.split()

        if len(line_elements[0]) > 5:
            reading_interior = True

        # Handle start of a new layer indicated by change in z-value.
        if float(line_elements[2]) != previous_z_value:
            new_layer = []
            layers.append(new_layer)
            reading_interior = False
            previous_z_value = float(line_elements[2])

        #                    xx.xxx           yy.yyy              AirTime
        tmp_vertex = Vertex(line_elements[0], line_elements[1], is_airtime)

        # Change units from the input decimeters to millimeters.
        tmp_vertex = Vertex(tmp_vertex.x / scale_down_factor, tmp_vertex.y  / scale_down_factor, tmp_vertex.isAirtime)

        if not reading_interior: # Scale up boundary vertices by 1.4 to fit interior.
            tmp_vertex = Vertex(tmp_vertex.x , tmp_vertex.y , tmp_vertex.isAirtime)
            #tmp_vertex = Vertex(tmp_vertex.x * 1.4, tmp_vertex.y * 1.4, tmp_vertex.isAirtime)


        #if len(layers) == 5:
        #    check_file.write(str(tmp_vertex.x) + ' ' + str(tmp_vertex.y) + '\r\n')

        # Print parameter of this vertex determines the airtime parameter
        # of the next vertex.
        if line_elements[3] == '0':
            is_airtime = True
        else:
            is_airtime = False

        # Delete the previous vertex if it and the current vertex are both
        # airtime. This allows the printer head to travel directly to the next
        # extruded vertex, instead of backtracking along already extruded edges.
        if tmp_vertex.isAirtime and layers[-1] and layers[-1][-1].isAirtime:
            del layers[-1][-1]

        layers[-1].append(tmp_vertex)

with open(out_file_path, 'w') as outfile:
    # Write Initial GCode Commands to file.
    outfile.write(';START_OF_HEADER\n')
    outfile.write(';HEADER_VERSION:0.1\n')
    outfile.write(';FLAVOR:Griffin\n')
    outfile.write(';GENERATOR.NAME:Cura_SteamEngine\n')
    outfile.write(';GENERATOR.VERSION:2.4.0\n')
    outfile.write(';GENERATOR.BUILD_DATE:2017-02-13\n')
    outfile.write(';TARGET_MACHINE.NAME:Ultimaker 3 Extended\n')
    outfile.write(';EXTRUDER_TRAIN.0.INITIAL_TEMPERATURE:205\n')
    outfile.write(';EXTRUDER_TRAIN.0.MATERIAL.VOLUME_USED:705\n')
    outfile.write(';EXTRUDER_TRAIN.0.MATERIAL.GUID:0e01be8c-e425-4fb1-b4a3-b79f255f1db9\n')
    outfile.write(';GENERATOR.VERSION:2.4.0\n')
    outfile.write(';EXTRUDER_TRAIN.0.NOZZLE.DIAMETER:0.4\n')
    outfile.write(';BUILD_PLATE.INITIAL_TEMPERATURE:60\n')
    outfile.write(';PRINT.TIME:749\n')
    outfile.write(';PRINT.SIZE.MIN.X:104.675\n')
    outfile.write(';PRINT.SIZE.MIN.Y:6\n')
    outfile.write(';PRINT.SIZE.MIN.Z:0.27\n')
    outfile.write(';PRINT.SIZE.MAX.X:170\n')
    outfile.write(';PRINT.SIZE.MAX.Y:119.325\n')
    outfile.write(';PRINT.SIZE.MAX.Z:9.97\n')
    outfile.write(';END_OF_HEADER\n')
    outfile.write(';Generated with Cura_SteamEngine 2.4.0\n\n')

    outfile.write('T0\n')
    outfile.write('G92 E0\n\n')
    
    outfile.write('M109 S205\n')
    outfile.write('G0 F15000 X170 Y6 Z2\n')
    outfile.write('G280\n')
    outfile.write('G1 F1500 E-6.5\n')
    outfile.write(';LAYER_COUNT:{0}\n'.format(len(layers)))
    
    outfile.write('M109 S200 ; wait for temperature to be reached\n')
    outfile.write('G21 ; set units to millimeters\n')
    outfile.write('G90 ; use absolute coordinates\n')
    outfile.write('M82 ; use absolute distance for extrusion\n')
    outfile.write('G92 E0\n')

    for layer in range(1, len(layers) + 1):
        outfile.write(';LAYER:{0}\n\n'.format(layer))

        if layer == 1:
            outfile.write('G1 E{0} F{1}\n'.format(E_purge,Flayer1))
            outfile.write('G92 E0\n')
            outfile.write('G1 X3 Y0\n')
            outfile.write('G1 Z-1.33\n')
            outfile.write('G1 F{1} E{0}\n'.format(FlayerZ,E_val* 0.4))
        E_val = 0
        for (index, vertex) in enumerate(layers[layer - 1]):
            # Continue to second loop iteration if on the first
            # so we have Vertices[index] and Vertices[index - 1] defined.
            #if index == 0:
                #outfile.write('G1 X{0} Y{1} ; TRAVEL\n'.format(vertex.x - start_x, vertex.y - start_y))
                #continue
            

            airTime = layers[layer-1][index].isAirtime
            #if layer == 1:
                #airTime=1
            distance = ((vertex.x)**2+(vertex.y)**2)**(1/2)
            E_val = E_val + distance * E_ratio #eparameter
            

            #outfile.write(airTime)
            #if layer != 1:
                #if index + 1 >= len(layers[layer - 1]) or layers[layer - 1][index + 1].isAirtime:
                    #outfile.write('G1 Z{0} F{1}\n'.format(layer_thickness,FlayerZ))
                    #outfile.write('G1 X{0} Y{1} E{2} F{3}\n'.format(vertex.x, vertex.y,E_val,FlayerZ))
                    

            # Turn off extruder and raise if airtime.
            if airTime:
                outfile.write('G1 Z{0} F{1}\n'.format(layer_thickness,FlayerZ))
                outfile.write('G1 X{0} Y{1} E0\n'.format(vertex.x, vertex.y))
                outfile.write('G1 Z0\n')
            else:
                if index == 0:
                    outfile.write('G1 F{0} X{1} Y{2} E{3}\n'.format(FlayerZ,vertex.x,vertex.y,E_val* 0.4))
                else:
                    outfile.write('G1 X{0} Y{1} E{2}\n'.format(vertex.x,vertex.y,E_val* 0.4))
            

            # If we are on the last edge or the next edge is airtime.
            
        

        # Turn extruder off and raise to next layer height.
        #outfile.write('G1 F{0} Z{1} \n'.format(FlayerZ,layer_thickness))

    #outfile.write('G1 Z3\n')
    #outfile.write('M104 S0\n')
    #outfile.write('T1 S0\n')
    #outfile.write(';End of Gcode\n')
    #outfile.write(';SETTING_3 {"extruder_quality": ["[general]\\nversion = 2\\nname = empty\\ndefin\n')
    #outfile.write(';SETTING_3 ition = ultimaker3_extended\\n\\n[metadata]\\ntype = quality_changes\\n')
    #outfile.write(';SETTING_3 \nquality_type = normal\\nextruder = ultimaker3_extended_extruder_lef\n')
    #outfile.write(';SETTING_3 t\\n\\n[values]\\n\\n", "[general]\\nversion = 2\\nname = empty\\ndef\n')
    #outfile.write(';SETTING_3 inition = ultimaker3_extended\\n\\n[metadata]\\ntype = quality_change\n')
    #outfile.write(';SETTING_3 s\\nquality_type = normal\\nextruder = ultimaker3_extended_extruder_r\n')
    #outfile.write(';SETTING_3 ight\\n\\n[values]\\n\\n"], "global_quality": "[general]\\nversion = \n')
    #outfile.write(';SETTING_3 2\\nname = empty\\ndefinition = ultimaker3_extended\\n\\n[metadata]\\\n')
    #outfile.write(';SETTING_3 ntype = quality_changes\\nquality_type = normal\\n\\n[values]\\n\\n"}\n')
