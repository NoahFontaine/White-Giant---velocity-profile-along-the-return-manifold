import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm

V_max = 2e-2

# Calcuates the velocity at any point along the cross-section of the manifold
def v(x,y,theta):
    # Calculate R (m) and Y_bar (m)
    R = (2.17e-3-27e-3)/(2*np.pi)*theta + 27e-3
    Y_bar = 4*R/(3*np.pi)

    if y <= Y_bar:
        vel = V_max * y/Y_bar * (1 - (x**2+y**2)/R**2) * (1 - (y-Y_bar)**2/Y_bar**2)
    else:
        vel = V_max * y/Y_bar * (1 - (x**2+y**2)/R**2) * (1 - (y-Y_bar)**2/(R-Y_bar)**2)

    if vel < 0:
        return float(0)
    else:
        return vel

#Velocity profile at entrance
angle = 0
radius = (2.17e-3-27e-3)/(2*np.pi)*angle + 27e-3

x_range = np.arange(-radius, radius, 1e-3)
y_range = np.arange(0, radius, 5e-4)
list = []

for x0 in x_range:
    for y0 in y_range:
        if x0**2 + y0**2 <= radius**2:
            list.append(x0)
            list.append(y0)
            list.append(v(x0,y0,angle))

X = np.array(list[::3])
Y = np.array(list[1::3])
Z = np.array(list[2::3])

#Plot the return manifold. (x,y,z)_m are points along the inside wall of the manifold, touching the nozzle

u_m = np.linspace(0, 2.0 * np.pi, endpoint=True, num=200) #Angle along the manifold (0 -> 2pi)
z_m = np.array([])
x_m = np.array([])
y_m = np.array([])
for u in u_m:
    j = np.arange(-27e-3, 2*(2.17e-3-27e-3)/(2*np.pi)*u + 27e-3, 1e-3) #Finds the radius of the manifold at every angle and plots equally spaced points along it
    z_m = np.append(z_m, j) #High of the manifold
    #Plot the points
    for i in range(len(j)):
        x_m = np.append(x_m, 150e-3*np.cos(u))
        y_m = np.append(y_m, 150e-3*np.sin(u))

# Points on curved side of the manifold. z coordinates do not change
x_t = np.array([])
y_t = np.array([])

for u in u_m:
    rad = (2.17e-3 - 27e-3)/(2*np.pi)*u + 27e-3
    z_t = np.arange(-27e-3, 2*rad - 27e-3, 1e-3)
    a_t = np.arange(-rad, rad, 1e-3) #relative position of point on coordinate centered about semicircle
    x_t = np.append(x_t, np.cos(u)*np.sqrt(rad**2-a_t**2))
    y_t = np.append(y_t, np.sin(u)*np.sqrt(rad**2-a_t**2))

#Plot detailed velocity profile
fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
t = ax.plot_trisurf(X,Y,Z, cmap=plt.cm.CMRmap, antialiased=True)
ax.set_xlim(-30e-3, 30e-3)
ax.set_ylim(-15e-3, 30e-3)
ax.set_zlim(0, 0.03)
fig.colorbar(t, shrink=0.5, aspect=5)
plt.show()

#Plot velocity profile on the manifold
fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
ax.scatter(Y+x_m[0],Z,X, color='purple', s=0.1)
ax.scatter(x_m,y_m,z_m, color='grey', alpha=0.1, s=1)
ax.scatter(x_m+x_t, y_m+y_t, z_m, color='grey', alpha=0.1, s=1)
ax.set_xlim(-0.20, 0.20)
ax.set_ylim(-0.20, 0.20)
ax.set_zlim(-0.20, 0.20)
plt.show()