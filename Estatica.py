import numpy as np
import sympy as sp
import FisicaII as f2
import optparse

# Pedimos torques e momentos
parser = optparse.OptionParser(description = 'Torques necesarios para contrarrestar momento e forza para LD1.')
parser.add_option('-j', '--joints', help = 'Coordenadas ou distancias (arts. 2, 5) das articulacións', action = 'store')
parser.add_option('-f', '--force', help = 'Forzas do elemento terminal (x, y, z)')
parser.add_option('-m', '--moment', help = 'Momentos do elemento terminal (x, y, z)')
parser.set_defaults(joints = '0, 0, 0, 0, 0, 0, 0, 0, 0, 0', force = '0, 0, 0', moment = '0, 0, 0')
options, arguments = parser.parse_args()

j = (str(options.joints).split(', '))
f = (str(options.force).split(', '))
m = (str(options.moment).split(', '))

joints = []
force = []
moment = []

for i in range(10):
    joints.append(np.float(j[i]))

for i in range(3):
    force.append(np.float(f[i]))
    moment.append(np.float(f[i]))

tq_vec = np.r_[moment, force] # Creamos o vector de momentos e forzas
torques = f2.Torq(tq_vec, joints) # Chamamos a funcion de torques
print('Torques necesarios para contrarrestar as forzas e o momento aplicadas ao elemento terminal:\n', torques)