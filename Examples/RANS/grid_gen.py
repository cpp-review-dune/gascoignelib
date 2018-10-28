import numpy as np
import pathlib

nelem = 30
np_x = 10
np_y = 3
nodes = 17
n_bd = 24
x = np.array([0, 1, 1.75, 2, 2.25, 2.75, 3, 3.25, 4, 5, 6])
y = np.array([0, 0.5, 1])
grid = np.meshgrid(x,y)
coords = np.hstack((grid[0].reshape(-1,1), grid[1].reshape(-1,1)))
grid_string =''
for i in range(nelem):
    grid_string += str(i) + ' ' + ' '.join(str(cell) for cell in coords[i]) +' 0\n'


print(grid_string)
nidx = np.arange(nodes).reshape((-1,1))
seccol = np.zeros((nodes,1))
quad_string = ''
g = np.arange(nelem)
g = g.reshape((3,10))
g = np.flip(g,1)
idx = 0
for i in range((np_x-1)*(np_y-1) + 1):
    if i == 3 or i == (np_x-1):
        pass
    else:
        quad_string += str(idx) + ' 0 quad '\
                        + str(i) + ' ' + str(i+1) + ' ' + str(i+1+np_x) + ' ' + str(i+np_x) + '\n'
        idx = idx + 1
line_string = '\n'.join((str(i)+' line ') for i in np.arange(n_bd))
f_name = pathlib.Path('wall_mounted.inp')
if f_name.exists():
    print('File exists, will be overwritten!')
with f_name.open(mode='wt') as file:
    file.write(str(nelem) + ' ' + str(nodes + n_bd) + ' 0 0 0\n')
    file.write(grid_string)
    file.write(quad_string)
    file.write(line_string)
