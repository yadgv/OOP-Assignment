from ThreeDObjClass import *

x_min, x_max = -500,500
y_min, y_max = -500,500
z_min, z_max = -500,500

def plot_all_shapes():
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim([x_min, x_max])
    ax.set_ylim([y_min, y_max])
    ax.set_zlim([z_min, z_max])

    wall = Box(width=1,height=20,length=20, color="purple")
    sphere = Sphere("blue", 100, 200, 0, 0)
    rhombic = Rhombicosidodecahedron("cyan", 100, 200, 500, 0)
    heart = Heart("red", 50, 100, -200, 0, 0)
    cone = RoundedCone("orange", 50, 20, 100, 0, -200, 0)

    sphere.display(ax)
    rhombic.display(ax)
    heart.display(ax)
    cone.display(ax)
    wall.display(ax)

    plt.show()


plot_all_shapes()

