import matplotlib.pyplot as plt
import matplotlib.patches as ptc
import numpy as np
from shapely.geometry import Point
from shapely.ops import unary_union

n = 100
size = 0.02
alpha = 0.5

def points():
    x = np.random.uniform(size=n)
    y = np.random.uniform(size=n)
    return x, y

x1, y1 = points()
x2, y2 = points()
polygons1 = [Point(x1[i], y1[i]).buffer(size) for i in range(n)]
polygons2 = [Point(x2[i], y2[i]).buffer(size) for i in range(n)]
print(len(polygons1))
print(type(polygons1[0]))
print(polygons1[0].is_valid)
polygons1 = unary_union(polygons1)
print(type(polygons1))
polygons2 = unary_union(polygons2)

fig = plt.figure(figsize=(4,4))
ax = fig.add_subplot(111, title="Test scatter")
for polygon1 in polygons1:
    polygon1 = ptc.Polygon(np.array(polygon1.exterior), facecolor="red", lw=0, alpha=alpha)
    ax.add_patch(polygon1)
for polygon2 in polygons2:
    polygon2 = ptc.Polygon(np.array(polygon2.exterior), facecolor="blue", lw=0, alpha=alpha)
    ax.add_patch(polygon2)
ax.axis([-0.2, 1.2, -0.2, 1.2])

fig.savefig("test_scatter.png")