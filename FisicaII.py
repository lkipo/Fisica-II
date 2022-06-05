import numpy as np
import sympy as sp

def rotg(omg, theta): 
    # Matriz de rotacion para vector + angulo
    

    c=np.cos(theta)
    s=np.sin(theta)
    w=omg/np.linalg.norm(omg)

    a = c + ((w[0])**2) * (1-c)
    b = w[0] * w[1] * (1-c) - (w[2]*s)
    c = w[0] * w[2] * (1-c) + (w[1]*s)
    d = w[0] * w[1] * (1-c) + (w[2]*s)
    e = c + ((w[1])**2) * (1-c)
    f = w[1] * w[2] * (1-c) - (w[0]*s)
    g = w[0] * w[2] * (1-c) - (w[1]*s)
    h = w[1] * w[2] * (1-c) + (w[0]*s)
    i = c + ((w[2])**2) * (1-c)
    
    return np.array([[a, b, c],
                     [d, e, f],
                     [g, h, i]])

def normalize(V):

    return V / np.linalg.norm(V)

def VecToso3(v):
    # Matriz antisimetrica para vector dado
    return np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])

def so3ToVec(so3mat):
    # Vector para matrix antisimetrica dada
    return np.array([so3mat[2][1], so3mat[0][2], so3mat[1][0]])

def VecTose3(V): # convierte un vector giro (1x6) o eje helicoidal en matriz 4x4 se3
    return np.r_[np.c_[VecToso3([V[0], V[1], V[2]]), [V[3], V[4], V[5]]], np.zeros((1, 4))]

def se3ToVec(se3mat): # Convierte una matriz se3 en un vector giro 1x6
    return np.r_[[se3mat[2][1], se3mat[0][2], se3mat[1][0]],
                 [se3mat[0][3], se3mat[1][3], se3mat[2][3]]]   

def MatrixExp3(w, theta):
    # Formula de Rodrigues
    w=w/np.linalg.norm(w)
    matriz = VecToso3(w)
    
    return np.eye(3) + np.sin(theta) * matriz + (1 - np.cos(theta)) * np.dot(matriz, matriz)

def MatrixLog3(R):
    traza = np.trace(R)
    if traza >= 3:
        return np.zeros(3)
    elif traza <= -1:
        if not abs((1 + R[2][2])) < 1e-6:
            w = (1.0 / np.sqrt(2 * (1 + R[2][2]))) * np.array([R[0][2], R[1][2], 1 + R[2][2]])
        elif not abs((1 + R[1][1])) < 1e-6:
            w = (1.0 / np.sqrt(2 * (1 + R[1][1]))) * np.array([R[0][1], 1 + R[1][1], R[2][1]])
        else:
            w = (1.0 / np.sqrt(2 * (1 + R[0][0]))) * np.array([1 + R[0][0], R[1][0], R[2][0]])
        return w, np.pi
    else:
        theta = np.arccos(0.5*(traza-1))
        so3 = (1/(2*np.sin(theta)))*(R-np.transpose(R))
        return so3ToVec(so3), theta
'''def MatrixLog3(R):
    traza = np.trace(R)
    if traza >= 3:
        return np.zeros(3)
    elif traza <= -1:
        if not abs((1 + R[2][2])) < 1e-6:
            w = (1.0 / np.sqrt(2 * (1 + R[2][2]))) * np.array([R[0][2], R[1][2], 1 + R[2][2]])
        elif not abs((1 + R[1][1])) < 1e-6:
            w = (1.0 / np.sqrt(2 * (1 + R[1][1]))) * np.array([R[0][1], 1 + R[1][1], R[2][1]])
        else:
            w = (1.0 / np.sqrt(2 * (1 + R[0][0]))) * np.array([1 + R[0][0], R[1][0], R[2][0]])
        return w, np.pi
    else:
        theta = np.arccos(0.5*(traza-1))
        so3 = (1/(2*np.sin(theta)))*(R-np.transpose(R))
        return so3ToVec(so3), theta'''

def RpToTrans(R, p):
    return np.array([[R[0, 0], R[0, 1], R[0, 2], p[0]],
                     [R[1, 0], R[1, 1], R[1, 2], p[1]],
                     [R[2, 0], R[2, 1], R[2, 2], p[2]],
                     [0, 0, 0, 1]])

def TransToRp(R):
    M = np.array([[R[0, 0], R[0, 1], R[0, 2]],
                  [R[1, 0], R[1, 1], R[1, 2]],
                  [R[2, 0], R[2, 1], R[2, 2]]])
    p = np.array([R[0, 3], R[1, 3], R[2, 3]])
    return (M, p)

def TransInv(T):
    (r, p) = TransToRp(T)

    inv = np.hstack((r.T, np.dot(-r.T, p).reshape(3, 1)))
    inv = np.vstack((inv, np.array([0, 0, 0, 1])))
    return inv

def Adj(T): # Calcula la matriz adjunta de una MTH
    R=T[0: 3, 0: 3]; p = T[0: 3, 3]
    return np.r_[np.c_[R, np.zeros((3, 3))], np.c_[np.dot(VecToso3(p), R), R]]

def VecTose3(V): # vector a antisimetrica 4x4
    return np.r_[np.c_[VecToso3([V[0], V[1], V[2]]), [V[3], V[4], V[5]]], np.zeros((1, 4))]

def ScrewToAxis(q, s, h):
    return([np.concatenate([s, np.cross(q, s) + np.dot(h, s)])])

def MatrixExp6(se3mat):
    se3mat = np.array(se3mat) 
    v = se3mat[0: 3, 3]
    omgmattheta = se3mat[0: 3, 0: 3]
    omgtheta = so3ToVec(omgmattheta) 

    if (np.linalg.norm(omgtheta))<1.e-6: 
        return np.r_[np.c_[np.eye(3), v], [[0, 0, 0, 1]]] 

    else:
        theta = np.linalg.norm(omgtheta)
        omgmat = omgmattheta / theta

        G_theta = np.eye(3)*theta + (1-np.cos(theta))*omgmat + (theta-np.sin(theta))*np.dot(omgmat,omgmat)
        R = np.eye(3) + np.sin(theta)*omgmat + (1.-np.cos(theta))*np.dot(omgmat,omgmat)
        return np.r_[np.c_[R,np.dot(G_theta,v)/theta],[[0, 0, 0, 1]]]

def MatrixLog6(T): # Calcula la matriz logaritmo de una MTH
    R=T[0: 3, 0: 3]; p = T[0: 3, 3]  # separa la MTH en matriz de rotación y vector traslación
    omgmat = MatrixLog3(R) # coordenadas exponenciales de la matriz de rotación
                           # o sea, un vector de rotación como matriz antisimétrica so3 (3x3)
    if np.array_equal(omgmat, np.zeros((3, 3))): # Si no hay rotación, es una matriz de ceros 
        return np.r_[np.c_[np.zeros((3, 3)),p],[[0, 0, 0, 0]]]
    else:
        omgvec = so3ToVec(omgmat) # expresa la rotación como un vector en la dirección del eje por el ángulo
        omgmat=omgmat/np.linalg.norm(omgvec) # el vector en el eje de rotación normalizado y en forma matricial
        theta = np.linalg.norm(omgvec) # también se puede calcular como np.arccos((np.trace(R)-1)/2.0)
        # a continuación aplicamos la definición que vimos en clase (ver diapositivas)
        invG_theta=np.eye(3)/theta-omgmat*0.5+(1.0/theta-0.5/np.tan(theta*0.5))*np.dot(omgmat,omgmat)
        v=np.dot(invG_theta,p)
        return np.r_[np.c_[omgmat,v],[[0, 0, 0, 0]]]*theta # primero concatena columnas y luego filas

def R2Euler(R):
    sy=np.sqrt(R[0,0]*R[0,0]+R[1,0]*R[1,0])
    singular=sy<1.e-6
    if not singular:
        x=np.arctan2(R[2,1],R[2,2])
        y=np.arctan2(-R[2,0],sy)
        z=np.arctan2(R[1,0],R[0,0])
        x=np.arctan2(-R[1,2],R[1,1])
        y=np.arctan2(-R[2,0],sy)
        z=0.
    return np.array([x,y,z])

def Hel6(w, q):
    if (w[0] == 0 and w[1] == 0 and w[2] == 0):
        hel = [0, 0, 0]
        for i in range(3):
            if q[i] == 0:
                hel.append(0)
            else:
                hel.append(q[i]/abs(q[i]))
        return(np.array(hel))
    else:
        return((np.r_[w, np.cross(q, w)]))

def EixosLD1():
    w = [] # Array de eixos
    w.append(np.array((0, 0, 1))) # Articulacion 1
    w.append(np.array((0, 0, 0))) # Articulacion 2
    w.append(np.array((0, 0, 1))) # Articulacion 3
    w.append(np.array((0, -1, 0))) # Articulacion 4
    w.append(np.array((0, 0, 0))) # Articulacion 5
    w.append(np.array((-1, 0, 0))) # Articulacion 6
    w.append(np.array((-1, 0, 0))) # Articulacion 7
    w.append(np.array((-1, 0, 0))) # Articulacion 8
    w.append(np.array((0, -1, 0))) # Articulacion 9
    w.append(np.array((0, 0, -1))) # Articulacion 10
    return w

def VectoresLD1():
    L = np.array([1.1, 1, 1.1, 0.4, 0.4, 0.05]) # Eslabons
    q = [] # array de vectores de cada eixo ao seguinte (q)
    q.append(np.array([0, L[0], 0])) #1
    q.append(np.array([0, 0, L[1]])) #2
    q.append(np.array([0, 0, 0])) #3
    q.append(np.array([0, 0, 0])) #4
    q.append(np.array([0, -L[2], 0])) #5
    q.append(np.array([0, 0, -L[3]])) #6
    q.append(np.array([0, 0, -L[4]])) #7
    q.append(np.array([0, 0, 0])) #8
    q.append(np.array([0, 0, 0])) #9
    q.append(np.array([0, 0, 0])) #10
    return q

def Jacobian(theta):
    w = EixosLD1()
    q = VectoresLD1()

    t=sp.symbols('t0, t1, t2, t3, t4, t5, t6, t7, t8, t9')  #Coordenadas de las articulaciones

    # Calculamos las matrices de rotación a partir de los ejes w, utilizando la fórmula de Rodrigues
    R=[]
    for i in range(10):
        wmat=sp.Matrix(VecToso3(w[i]))
        R.append(sp.eye(3)+sp.sin(t[i])*wmat+(1-sp.cos(t[i]))*(wmat*wmat))

    # Aplicamos rotaciones a los vectores q y w para llevarlos a la configuración del robot que queremos
    qs=[]; ws=[]; Ri=R[0]
    qs.append(sp.Matrix(q[0]))
    ws.append(sp.Matrix(w[0]))
    for i in range(1,10,1):
        ws.append(Ri*sp.Matrix(w[i]))
        qs.append(Ri*sp.Matrix(q[i])+qs[i-1])
        Ri=Ri*R[i]

    # Calculamos las velocidades lineales, los vectores giro correspondientes y la matriz Jacobiana
    vs=[]; Ji=[]; i=0
    vs.append(qs[i].cross(ws[i]))
    Ji.append(ws[i].row_insert(3,vs[i]))
    J=Ji[0]
    for i in range(1,10,1):
        vs.append(qs[i].cross(ws[i]))
        Ji.append(ws[i].row_insert(3,vs[i]))
        J=J.col_insert(i,Ji[i])

    return np.array(J.subs({t[0]:theta[0], t[1]:theta[1], t[2]:theta[2], t[3]:theta[3], t[4]:theta[4], t[5]:theta[5], t[6]:theta[6], t[7]:theta[7], t[8]:theta[8]}))

def velET(thetalist, vel):
    Jaco = Jacobian(thetalist)
    return np.dot(Jaco, vel)

def Torq(MomentTorq, thetalist):
    J = Jacobian(thetalist)
    Js = []
    print('Matriz Jacobiana:\n', J)
    torques = np.dot(MomentTorq, J)
    return(torques)

def CinematicaDirectaLD1(coordenadas):
    w = EixosLD1() # Eixos de rotacion (w)
    q = VectoresLD1() # Vectores de cada eixo ao seguinte (q)

    # Calcular velocidades lineais   
    qs=[]; ws=[]
    qs.append(np.array(q[0]))
    ws.append(np.array(w[0]))
    for i in range(1,10,1):
        ws.append(np.array(w[i]))
        qs.append(np.array(q[i])+qs[i-1])

    # Calcular eixos helicoidais
    v = []
    #print('\nEixos helicoidais: ')
    for i in range(10):
        v.append(Hel6(ws[i], qs[i]))
        #print(v[i], 'Eixo #%i' %(i+1))

    mth = []
    T = np.eye(4)
    # Matriz m con respecto ao sistema de referencia da base
    M = np.array([[-1, 0, 0, 0],
                  [0, -1, 0, 0],
                  [0, 0, -1, 0.05],
                  [0, 0, 0, 1]])

    #print('\nMatrices de transformación homoxéneas: ')

    for n in range(10):
        mth.append(MatrixExp6(VecTose3(v[n] * coordenadas[n])))
        # print('Matriz no elemento %i: \n' %(n+1), np.round(mth[n], 3)) # DESCOMENTAR PARA VER MATRICES INDIVIDUALMENTE
        T = np.dot(T, mth[n])
        #print('Matriz despois do elemento %i: \n' %(n+1), np.round(T, 3))
    T = np.dot(T, M)

    return T

def getT(orientation, r):
    Ri=MatrixExp3(np.array([1, 0, 0]), orientation[0])
    Rj=MatrixExp3(np.array([0, 1, 0]), orientation[1])
    Rk=MatrixExp3(np.array([0, 0, 1]), orientation[2])
    R = np.matmul(Rk, np.matmul(Rj, Ri))
    aux = np.array([[0, 0, 0, 1]])
    return np.r_[np.c_[R, r], aux]