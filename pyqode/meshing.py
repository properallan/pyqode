from pathlib import Path
import numpy as np


class Gmsh:
    pcount = 0
    lcount = 0
    loopcount = 0
    scount = 0
    phycount = 0
    gfile = ''
    blcount = 0

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
    
    def addBoundaryLayer(self, edges, nodes, thickness=0.02, ratio=1.2, hwall_n=0.00001, beta=1.001, nblayers=10):
        self.blcount += 1

        edges_list = '{'
        for i, edge in enumerate(edges):
            edges_list += str(edge)
            if i < np.size(edges) -1 :
                edges_list +=  ', '
        edges_list += '}'

        nodes_list = '{'
        for i, node in enumerate(nodes):
            nodes_list += str(node)
            if i < np.size(nodes) -1 :
                nodes_list +=  ', '
        nodes_list += '}'

        self.gfile += f"""//Boundary Layer
Field[{self.blcount}] = BoundaryLayer;
//Field[{self.blcount}].AnisoMax = 1000;
Field[{self.blcount}].Quads = 1;
Field[{self.blcount}].BetaLaw = 1;
Field[{self.blcount}].Beta = {beta};
//Field[{self.blcount}].Thickness = {thickness};
Field[{self.blcount}].EdgesList = {edges_list};
Field[{self.blcount}].NbLayers = {nblayers};
Field[{self.blcount}].NodesList = {nodes_list};
//Field[{self.blcount}].Ratio = {ratio};
Field[{self.blcount}].hwall_n = {hwall_n};
//Field[{self.blcount}].hfar = 1;
Field[{self.blcount}].IntersectMetrics = 0;

BoundaryLayer Field = {self.blcount};
//Mesh.Algorithm = 6;
//Recombine Surface {1};"""