import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import math
from scipy.spatial import ConvexHull
import cmath



class ThreeDShape:
    # type = ""
    # color = ""


    def __init__(self, type, color, x, y, z):
        self.type = type
        self.color = color
        self.x = x
        self.y = y
        self.z = z
        self.object = None

    def volume(self):
        return None
    
    def surfaceArea(self):
        return None
    
    def getType(self):
        return f"This shape is a {self.type}"

    def setColor(self, color):
        self.color = color
    
    def getColor(self):
        return f"The color of this shape is {self.color}"
    
    def display(self, ax):
        pass

    def __del__(self):
        print(f"Deleting {self.type}...")
        # object.remove()
        

class Sphere(ThreeDShape):
    def __init__(self, color="red", radius=50, x=0, y=0, z=0, type="Sphere"):
        super().__init__(type, color, x, y, z)
        self.radius = radius
    
    def volume(self):
        return ((4/3)*math.pi*float(self.radius)**3)
    
    def surfaceArea(self):
        return 4*math.pi*self.radius**2
    
    def setRadius(self, radius):
        self.radius = radius
    
    def display(self, ax):
        r = self.radius
        theta = np.linspace(0, 2 * np.pi, 100)
        phi = np.linspace(0, np.pi, 100)

        x = self.x + r * np.outer(np.cos(theta), np.sin(phi))
        y = self.y + r * np.outer(np.sin(theta), np.sin(phi))
        z = self.z + r * np.outer(np.ones(100), np.cos(phi))

        ax.plot_surface(x, y, z, color=self.color, alpha=0.6)

    

class Rhombicosidodecahedron(ThreeDShape):
    def __init__(self, color="red", eLength=50, x=0, y=0, z=0, type="Rhombicosidodecahedron"):
        super().__init__(type, color, x, y, z)
        self.edgeLength = eLength
    
    def volume(self):
        return (float(self.edgeLength)**3/3*(60+29*math.sqrt(5)))
    
    def surfaceArea(self):
        return (float(self.edgeLength)**2 * (30 + 5*math.sqrt(3) + 3*math.sqrt(25+10*math.sqrt(5))))
    
    def setEdgeLength(self, eLength):
        self.edgeLength = eLength

    def display(self, ax):
        def create_rhombicosidodecahedron():
            phi = (1 + np.sqrt(5)) / 2
            vertices = np.array([
                [0, 1, 3 * phi], [0, -1, 3 * phi], [0, 1, -3 * phi], [0, -1, -3 * phi],
                [1, 3 * phi, 0], [-1, 3 * phi, 0], [1, -3 * phi, 0], [-1, -3 * phi, 0],
                [3 * phi, 0, 1], [3 * phi, 0, -1], [-3 * phi, 0, 1], [-3 * phi, 0, -1],
            ])
            hull = ConvexHull(vertices)
            return vertices, hull.simplices

        vertices, faces = create_rhombicosidodecahedron()
        ax.add_collection3d(Poly3DCollection([vertices[face] for face in faces], facecolors=self.color, linewidths=1, edgecolors=self.color, alpha=0.7))

class Heart(ThreeDShape):
    def __init__(self, color="red", eLength=50, height=100, x=0, y=0, z=0, type="Heart"):
        super().__init__(type, color, x, y, z)
        self.edgeLength = eLength
        self.height = height
    
    def volume(self):
        return(((1+math.pi/4)*float(self.edgeLength)**2)*float(self.height))
    
    def surfaceArea(self):
        return(((1+math.pi/4)*float(self.edgeLength)**2)*2) + ((self.edgeLength*self.height)*2) + ((1+math.pi/4)*self.height*2) #this looks weird because i made it myself, (it should technically work)
    
    def display(self, ax):
        def heart_2d(t, edgeLength):
            x = self.x + 16 * np.sin(t)**3 * edgeLength / 16
            y = self.y + (13 * np.cos(t) - 5 * np.cos(2*t) - 2 * np.cos(3*t) - np.cos(4*t)) * edgeLength / 16
            return x, y

        def create_heart_box(edgeLength=1, height=1, resolution=100):
            t = np.linspace(0, 2 * np.pi, resolution)
            x_top, y_top = heart_2d(t, edgeLength)
            z_top = np.full_like(x_top, height / 2)
            x_bottom, y_bottom = heart_2d(t, edgeLength)
            z_bottom = np.full_like(x_bottom, -height / 2)
            return x_top, y_top, z_top, x_bottom, y_bottom, z_bottom

        x_top, y_top, z_top, x_bottom, y_bottom, z_bottom = create_heart_box(self.edgeLength, self.height)
        verts_top = [list(zip(x_top, y_top, z_top))]
        verts_bottom = [list(zip(x_bottom, y_bottom, z_bottom))]
        
        ax.add_collection3d(Poly3DCollection(verts_top, facecolors=self.color, linewidths=1, edgecolors='r', alpha=0.7))
        ax.add_collection3d(Poly3DCollection(verts_bottom, facecolors=self.color, linewidths=1, edgecolors='r', alpha=0.7))

        for i in range(len(x_top) - 1):
            verts_side = [
                [x_top[i], y_top[i], z_top[i]],
                [x_top[i+1], y_top[i+1], z_top[i+1]],
                [x_bottom[i+1], y_bottom[i+1], z_bottom[i+1]],
                [x_bottom[i], y_bottom[i], z_bottom[i]]
            ]
            ax.add_collection3d(Poly3DCollection([verts_side], facecolors=self.color, linewidths=1, edgecolors='r', alpha=0.7))


class RoundedCone(ThreeDShape):
    def __init__(self, color="blue", bottomRadius=100, topRadius=50, height=100, x=0, y=0, z=0, type="RoundedCone"):
        super().__init__(type, color, x, y, z)
        self.bottomRadius = bottomRadius
        self.topRadius = topRadius
        self.height = height
    
    def volume(self):
        R = self.bottomRadius
        r = self.topRadius
        j = self.height
        h = R - cmath.sqrt(R**2 - r**2)
        sqrt_term = h**2 + (R - r)**2
        a = cmath.acos((R - r) / cmath.sqrt(sqrt_term))
        # a = np.arccos(( ( R - r ) / math.sqrt(h**2 + ( R - r )**2)))
        i =  math.sqrt(r**2 / ( 2 / ( 1 - np.cos(a) ) - 1 ))
        return(math.pi/3 * ( j * ( R**2 + R*r + r**2 ) + i**2 * ( 3 * ( r**2 + i**2 ) / (2*i) - i )))

    def surfaceArea(self):
        R = self.bottomRadius
        r = self.topRadius
        j = self.height
        h = R - cmath.sqrt(R**2 - r**2)
        sqrt_term = h**2 + (R - r)**2
        a = cmath.acos((R - r) / cmath.sqrt(sqrt_term))
        i =  math.sqrt(r**2 / ( 2 / ( 1 - np.cos(a) ) - 1 ))
        return ( R + r ) * math.pi * math.sprt(( R - r )**2 + j**2) + math.pi*R**2 + math.pi * ( r**2 + i**2 )
    
    def display(self, ax):
        def create_rounded_cone(radius=5, height=10, resolution=100):
            theta = np.linspace(0, 2 * np.pi, resolution)
            z_cone = np.linspace(-height, 0, resolution)
            z_cone = z_cone[:, None]
            x_cone = radius * (1 - z_cone / -height) * np.cos(theta)
            y_cone = radius * (1 - z_cone / -height) * np.sin(theta)

            phi = np.linspace(0, np.pi / 2, resolution)
            theta_sphere = np.linspace(0, 2 * np.pi, resolution)
            phi, theta_sphere = np.meshgrid(phi, theta_sphere)
            x_sphere = radius * np.sin(phi) * np.cos(theta_sphere)
            y_sphere = radius * np.sin(phi) * np.sin(theta_sphere)
            z_sphere = radius * np.cos(phi)

            return x_cone, y_cone, z_cone, x_sphere, y_sphere, z_sphere

        x_cone, y_cone, z_cone, x_sphere, y_sphere, z_sphere = create_rounded_cone(self.bottomRadius, self.height)

        ax.plot_surface(x_cone, y_cone, z_cone, color=self.color, alpha=0.7)
        ax.plot_surface(x_sphere, y_sphere, z_sphere, color=self.color, alpha=0.7)


class Box(ThreeDShape):
    def __init__(self, type="Box", color="blue", x=0, y=0, z=0, width=50, height=50, length=50):
        super().__init__(type, color, x, y, z)
        self.height = height
        self.width = width
        self.length = length

    def setPlacement(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def volume(self):
        return self.width*self.length*self.height
    
    def surfaceArea(self):
        return 2*self.width*self.length + 2*self.length+self.height + 2*self.height*self.width
    
    def display(self, ax):
        x = np.arange(self.length + 1) + self.x   
        y = np.arange(self.height + 1) + self.y  
        z = np.arange(self.width + 1) + self.z  
        x, y, z = np.meshgrid(x, y, z, indexing='ij')
        data = np.ones((self.length, self.height, self.width), dtype=bool)
        self.object = ax.voxels(x, y, z, data, facecolors=self.color, linewidth=0.1)