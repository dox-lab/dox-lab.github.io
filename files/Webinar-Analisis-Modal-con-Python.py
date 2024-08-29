import openseespy.opensees as ops
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import opsvis as opsv

m = 1
kg = 1
s = 1
# Otras Unidades
cm = 0.01*m
N = kg*m/s**2
kgf = 9.81*N
Pa = N/m**2
MPa = 10**6*Pa
inch = 2.54*cm
ft = 12*inch
ksi = 6894757.2932*Pa
kip = ksi*inch**2
psi = 6894.76*Pa
# Constantes Físicas
g = 9.81*m/s**2


fce = 210*kg/cm**2
E = 15100*fce**0.5
G = 0.5*E/(1+0.2)
# Viga
b,h = 30*cm, 60*cm
Av = b*h
Izv = b*h**3/12
Iyv = b**3*h/12
aa, bb = max(b,h),min(b,h)
β= 1/3-0.21*bb/aa*(1-(bb/aa)**4/12)
Jxxv = β*bb**3*aa

# Columna
a = 40*cm
Ac = a**2
Izc = a**4/12
Iyc = a**4/12
β= 1/3-0.21*1.*(1-(1.)**4/12)
Jxxc = β*a**4

# Iniciar modelo
ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 6)

# Generar malla
dx, dz = 7 * m, 3.1 * m
nx, nz = 3, 3

# Crear nodos y elementos
Nodes, Elems = [], []
for i in range(nx+1):
    for j in range(nz+1):
        nodeId = (i+1) * 10 + j
        Nodes.append([nodeId, i*dx, 0., j*dz])
        
        # No crear vigas en z = 0
        if j > 0:  # Elementos verticales
            Elems.append([f'{nodeId-1}{nodeId}', nodeId-1, nodeId, 1])
        
        if i > 0 and j != 0:  # Elementos horizontales, excluyendo z=0
            Elems.append([f'{nodeId-10}{nodeId}', nodeId-10, nodeId, 2])

# Crear nodos en OpenSees
for Ni in Nodes:
    ops.node(int(Ni[0]), *Ni[1:4])
    
    # Restringimos solo los nodos que tienen z=0
    if Ni[3] == 0.:
        ops.fix(int(Ni[0]), *[1,1,1,1,1,1])
    else:
        ops.fix(int(Ni[0]), *[0,0,0,0,0,0])  # Los demás nodos tienen libertad completa

# Transformaciones geométricas
ops.geomTransf('PDelta', 1, *[0, 1, 0])
ops.geomTransf('Linear', 2, *[0, 1, 0])

# Crear elementos en OpenSees
for Ele in Elems:
    tag = int(Ele[0])  # Usamos el tag como identificador, aunque ahora es el nombre de la conexión de nodos
    if int(Ele[3]) == 1:  # Columna
        ops.element('elasticBeamColumn', tag, int(Ele[1]), int(Ele[2]), Ac, E, G*1000., Jxxc*1000., Iyc*1000., Izc, int(Ele[3]))
    else:  # Viga
        ops.element('elasticBeamColumn', tag, int(Ele[1]), int(Ele[2]), Av, E, G*1000., Jxxv*1000., Iyv*1000., Izv, int(Ele[3]))

# Visualizar el modelo
opsv.plot_model(fig_wi_he=(40., 32.), az_el=(-90, 0))
plt.show()

# Aplicando Cargas vivas y muertas
wLive = 250*kg/m**2
wLosa = 300*kg/m**2
wAcab = 100*kg/m**2
wTabi = 250*kg/m**2
wVigaColumna = 287.5*kg/m**2
wTotal = 1.0*(wLosa+wAcab+wTabi+wVigaColumna)+0.25*wLive
# print(wTotal)

# Asignación de masa
dy = dz  # Asunción, ajustar si es necesario
Carga = wTotal * dx * dy * m**2
for Ni in Nodes:
    ops.mass(int(Ni[0]), Carga, Carga, 0.0, 0.0, 0.0, 0.0)

# Obtenemos los modos
Nmodes = 9

vals = ops.eigen(Nmodes)
Tmodes = np.zeros(len(vals))
for i in range(Nmodes):
    Tmodes[i] = 2*np.pi/vals[i]**0.5
    print("T[%i]: %.5f" % (i+1, Tmodes[i]))

# Ahora que ya se han calculado los modos, intentamos visualizar las formas modales
opsv.plot_mode_shape(1)
opsv.plot_mode_shape(2)
opsv.plot_mode_shape(3)
plt.show()