# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 17:17:10 2019

@author: Home
"""

"""
Writing the solution and the mesh into a file readable by paraview
"""

import numpy as np
import os
import shutil

def write_file( mesh_,U,it):
    
	os.makedirs("output..vtu", exist_ok=True)
	out_file="output..vtu/output_" + str(it) + ".vtu"
	file = open(out_file,"w")
	file.write('<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">\n')
	file.write("<UnstructuredGrid>\n")
	file.write('<Piece NumberOfPoints="' + str(mesh_.Ns)+'" NumberOfCells="'+str(mesh_.Nt )+'">\n')
	file.write("<Points>\n")
	file.write('<DataArray NumberOfComponents="3" type="Float64">\n')
	for i in mesh_.Nodes:
		file.write(str(i.x) + " " + str(i.y) + " " + str(i.z))
		file.write('\n')

	file.write('</DataArray>\n</Points>\n<Cells>\n<DataArray type="Int32" Name="connectivity">\n')

	for i in range(0, mesh_.Nt):
		file.write(str(mesh_.Triangles[i].sommets[0]-1)+ " " + str(mesh_.Triangles[i].sommets[1]-1)+ " " + str(mesh_.Triangles[i].sommets[2]-1))
		file.write('\n')

	file.write('</DataArray>\n<DataArray type="Int32" Name="offsets">\n')

	for i in range(mesh_.Nt):
		file.write(str((i+1)*3) + '\n')

	file.write('</DataArray>\n<DataArray type="UInt8" Name="types">\n')

	for i in range(mesh_.Nt):
		file.write('5\n')

	file.write('</DataArray>\n</Cells>\n<PointData Scalars="solution">\n<DataArray type="Float64" Name="Real part" format="ascii">\n')
	for u in U:
		file.write(str(np.real(u)) + "\n")
	file.write('</DataArray>\n')
	file.write('<DataArray type="Float64" Name="Imag part" format="ascii">\n')
	for u in U:
		file.write(str(np.imag(u)) + "\n")
	file.write('</DataArray>\n')

	file.write('</PointData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n')
	file.close()
    
    
def erase_files():
    repertoire='output..vtu'
    shutil.rmtree(repertoire,ignore_errors=True)
#    for filename in os.listdir(repertoire) :
#        os.remove(repertoire + "/" + filename)
    return