import FisicaII as f2
import numpy as np
import sympy as sp
import optparse

parser = optparse.OptionParser(description="Resolucion de Cinematica Directa para LD1.")
parser.add_option('-s', '--seed', help = 'Coordenadas ou distancias (arts. 2, 5) iniciais das articulacións (seed)', action = 'store')
parser.add_option('-r', '--pos_xyz', help = 'Posición do elemento terminal (x, y, z)', action = 'store')
parser.add_option('-a', '--ang', help = 'Ángulos de Euler finais do elemento terminal', action = 'store')
parser.add_option('-o', '--er_or', help='Error na orientación do elemento terminal (0.01 por defecto)', action='store')
parser.add_option('-p', '--er_pos', help='Error na posición do elemento terminal (0.001 por defecto)', action='store')
parser.add_option('-i', '--max_iter', help='Número máximo de iteracións (20 por defecto)', action='store')

parser.set_defaults(seed = '0, 0, 0, 0, 0, 0, 0, 0, 0, 0', pos_xyz = '0.1, 0.1, 0.1', ang = '0, 0, 0', er_or = '0.01', er_pos = '0.001', max_iter = '20')
options, arguments = parser.parse_args()

# Pasamos os str a float e gardamolos nos arrays corresoindentes
sj = str(options.seed).split(', ')
seed_joints = []
for i in range (0,10,1):
    seed_joints.append(np.float(sj[i]))

xyz = str(options.pos_xyz).split(', ')
ang = str(options.ang).split(', ')
position = []
orientation = []
for i in range (0,3,1): 
    position.append(np.float(xyz[i]))
    orientation.append(np.float(ang[i]))

eomg = float(options.er_or)
ev = float(options.er_pos)
max_iter = float(options.max_iter)

# Matriz de transformacion homoxenea no elemento terminal na posicion 0
M = np.array([[-1, 0, 0, 0],
              [0, -1, 0, 0],
              [0, 0, -1, 0.05],
              [0, 0, 0, 1]])
T = f2.getT(orientation, position)

i = 0
Tsb = f2.CinematicaDirectaLD1(seed_joints)
a = np.dot(np.linalg.inv(Tsb), T)
print(a)
Vb = f2.MatrixLog6(np.dot(np.linalg.inv(Tsb), T))


Vs = np.dot(f2.Adj(Tsb), f2.se3ToVec(Vb))
err = (np.linalg.norm([Vs[0], Vs[1], Vs[2]]) > eomg) or (np.linalg.norm([Vs[3], Vs[4], Vs[5]]) > ev)
thetalist = np.array(seed_joints).copy()

while err and i < max_iter:
    print('while')
    J = f2.Jacobian(thetalist)
    
    thetalist = thetalist + np.dot(np.linalg.pinv(J), Vs)
    print(thetalist)

    Tsb = f2.CinematicaDirectaLD1(thetalist)
    print(Tsb)
    a = np.dot(np.linalg.inv(Tsb), T)
    print(a)
    Vb = f2.MatrixLog6(np.dot(np.linalg.inv(Tsb), T))
    Vs = np.dot(f2.Adj(Tsb), f2.se3ToVec(Vb))
    err = (np.linalg.norm([Vs[0], Vs[1], Vs[2]]) > eomg) or (np.linalg.norm([Vs[3], Vs[4], Vs[5]]) > ev)
    i += 1
    
print ("\n\nCoordenadas de las articulaciones:\n", thetalist)
print ("Error en w:", np.round(np.linalg.norm([Vs[0], Vs[1], Vs[2]]),8))
print ("Error en v:", np.round(np.linalg.norm([Vs[3], Vs[4], Vs[5]]),8))
print ("Número de iteraciones:", i)
