from pathlib import Path
import numpy as np


class Gmsh:
    pcount = 0
    lcount = 0
    loopcount = 0
    scount = 0
    phycount = 0
    gfile = ''

    def addPoint(self, x, y, z):
        self.pcount += 1
        self.gfile += 'Point('+str(self.pcount)+') = {'+str(x)+','+str(y)+','+str(z)+', cl__1};\n'
        return self.pcount

    def extrudeMesh(self, direction, line, layers):
        d = direction
        self.lcount += 3
        self.gfile += 'Extrude { ' + str(d[0]) + ', ' + str(d[1]) + ', ' + str(d[2]) + '} { Curve{'+ str(line) +'}; Layers {' + str(layers) + '}; Recombine;}\n'

        return self.lcount-2, self.lcount-1, self.lcount

    def addLine(self, points):
        self.lcount += 1
        self.gfile += 'Line('+str(self.lcount)+') = {'
        for i,p in enumerate(points):
            self.gfile += str(p)
            if i < np.size(points)-1:
                self.gfile += ', '
        self.gfile += '};\n'
        
        return self.lcount

    def addSpline(self, points):
        self.lcount += 1
        self.gfile += 'Spline('+str(self.lcount)+') = {'
        for i,p in enumerate(points):
            self.gfile += str(p)
            if i < np.size(points)-1:
                self.gfile += ', '
        self.gfile += '};\n'
        
        return self.lcount

    def addLoop(self, lines):
        self.loopcount += 1
        self.gfile += 'Line Loop('+str(self.loopcount)+') = {'
        for i,l in enumerate(lines):
            self.gfile += str(l)
            if i < np.size(lines)-1:
                self.gfile += ', '
        self.gfile += '};\n'
        
        return self.loopcount

    def addPlaneSurface(self, loops):
        self.scount += 1
        self.gfile += 'Plane Surface('+str(self.scount)+') = {'
        for i,l in enumerate(loops):
            self.gfile += str(l)
            if i < np.size(loops)-1:
                self.gfile += ', '
        self.gfile += '};\n'
        
        return self.scount
            
    def transfiniteLine(self, line, points, progression=1, using='Progression'):
        self.gfile += 'Transfinite Line {' + str(line) + '} = ' + str(points) + f' Using {using} ' + str(progression) + ';\n'

    def transfiniteSurface(self, surface, points):
        self.gfile += 'Transfinite Surface {' + str(surface) + '} = {'
        for i,p in enumerate(points):
            self.gfile += str(p)
            if i < np.size(points)-1:
                self.gfile += ', '
        self.gfile += '};\n'
        
    def recombineSurface(self, surface):
        self.gfile += 'Recombine Surface {' + str(surface) + '};\n'

    def addPhysicalLine(self, boundName, line):
        self.phycount +=1
        self.gfile += 'Physical Line(\"' + boundName +'\") = {' + str(line) + '};\n'

    def addPhysicalSurface(self, surface):
        self.phycount +=1
        self.gfile += 'Physical Surface(' + str(self.phycount) +') = {' + str(surface) + '};\n'
        
    def writeLine(self, s):
        self.gfile += s

    def print(self):
        print(self.gfile)

    def writeFile(self, filename):
        f = open(filename,'w')
        f.write(self.gfile)
        f.close()

    def addSymmetryCurve(self, direction, line):
        self.lcount += 1
        d = direction
        self.gfile += f"Symmetry {{{d[0]}, {d[1]}, {d[2]}, 0}} {{ Duplicata {{ Curve{{{line}}}; }} }}\n"

        return self.lcount