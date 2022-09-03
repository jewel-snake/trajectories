from evtk.hl import pointsToVTK
import numpy
import sys
data = numpy.transpose(numpy.genfromtxt(sys.argv[1],delimiter=','))
#print(data[1])
x = numpy.array([round(xi,9) for xi in data[0]])
y = numpy.array([round(xi,9) for xi in data[1]])
z = numpy.array([round(xi,9) for xi in data[2]])
pointsToVTK(sys.argv[1],x,y,z)
