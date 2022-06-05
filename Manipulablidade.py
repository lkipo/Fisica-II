import numpy as np
import sympy as sp
import FisicaII as f2
import optparse

# Pedimos torques e momentos
parser = optparse.OptionParser(description = 'Velocidades do elemento terminal a partir das velocidades nas articulacións')
parser.add_option('-j', '--joints', help = 'Coordenadas ou distancias (arts. 2, 5) das articulacións', action = 'store')
parser.add_option('-v', '--velocity', help = 'Velocidades das articulacións', action = 'store')
parser.set_defaults(joints = '0, 0, 0, 0, 0, 0, 0, 0, 0, 0', velocity = '0, 0, 0, 0, 0, 0, 0, 0, 0, 0')
options, arguments = parser.parse_args()

j = (str(options.joints).split(', '))
v = (str(options.velocity).split(', '))

joints = []
velocity = []
for i in range(10):
    joints.append(np.float(j[i]))
    velocity.append(np.float(v[i]))

vel = f2.velET(joints, velocity)
vlin = vel[3: 6]
vang = vel[0: 3]
print('Velocidade lineal do elemento terminal (x, y, z): \n', vlin)
print('Velocidade angular do elemento terminal (x, y, z): \n', vang)