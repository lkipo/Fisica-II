import FisicaII as f2
import numpy as np
import optparse

parser = optparse.OptionParser(description="Resolucion de Cinematica Directa para LD1.")
parser.add_option('-j', '--joints', help = 'Coordenadas ou distancias (arts. 2, 5) das articulacións', action = 'store')
parser.set_defaults(joints = '0, 0, 0, 0, 0, 0, 0, 0, 0, 0')
options, arguments = parser.parse_args()

# Pasamos os str a float e gardamolos nos arrays corresoindentes
joints = (str(options.joints).split(', '))

coordenadas = []
for i in range (10):
    coordenadas.append(np.float(joints[i]))

# Chamamos a funcion de cinematica directa
T = f2.CinematicaDirectaLD1(coordenadas)
print('Matriz de transformación homoxénea do elemento terminal:\n', np.round(T, 3))
print("\nCoordenadas (x,y,z) do TCP: ", np.round(T[0: 3, 3],3))