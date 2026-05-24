from math import sqrt
import numpy as np

def defactor(linha: str, f, tipos=None):
    a = []
    try:
        linha = next(f)
        while linha.strip():
            row = []
            for i, elemento in enumerate(linha.strip().split(',')):
                s = elemento.strip()
                t = tipos[i] if tipos and i < len(tipos) else None
                try:
                    if t is str:
                        row.append(s)
                    elif t in (int, float):
                        row.append(t(s))
                    else:
                        val = float(s)
                        row.append(int(val) if val == int(val) else val)
                except ValueError:
                    row.append(s)
            a.append(row)
            linha = next(f)
    except StopIteration:
        pass
    return a

Entrada = {}

with open("entrada_Ex_Viga.txt", "r", encoding="utf-8") as f:

    for linha in f:

        if linha == "*POINTS (id,x_point,y_point)\n":
            Entrada["POINTS"] = defactor(linha, f, [int, float, float])

        if linha == "*CURVES (id,point_ini,point_fim)\n":
            Entrada["CURVES"] = defactor(linha, f, [int, int, int])

        if linha == "*MATERIALS (id,E,nu)\n":
            Entrada["MATERIALS"] = defactor(linha, f, [int, float, float])

        if linha == "*PROPERTIES (id,material_id,A,I)\n":
            Entrada["PROPERTIES"] = defactor(linha, f, [int, int, float, float])

        if linha == "*MESH (curve_id,property_id,num_elem)\n":
            Entrada["MESH"] = defactor(linha, f, [int, int, int])

        if linha == "*POINT_LOADS (id,point_id,related_gdl,value)\n":
            Entrada["POINT_LOADS"] = defactor(linha, f, [int, int, int, float])

        if linha == "*DIST_LOADS (id,curve_id,value_ini,value_fim,direction:(x,y,l,t))\n":
            Entrada["DIST_LOADS"] = defactor(linha, f, [int, int, float, float, str])

        if linha == "*SPRING (point_id,related_gdl,value)\n":
            Entrada["SPRING"] = defactor(linha, f, [int, int, float])

        if linha == "*BC (point_id,related_gdl,value)\n":
            Entrada["BC"] = defactor(linha, f, [int, int, float])

# Preciso montar a matriz de rigidez
    # Preciso do angulo da curva
        # Preciso do tamanho da curva
""" 
A curva vai me dá os pontos então tenho n e m
"""

def geometria_curva(curve_id:int,Entrada:list):
    """
    curve_id: É o número da curva 
    Entrada: É o Dicionário de todos os dados da entrada
    return:: lista de propriedades geométricas da curva
    return: 
    """

    curve = Entrada["CURVES"][curve_id]
    a_y = Entrada["POINTS"][curve[1]-1][2]
    a_x = Entrada["POINTS"][curve[1]-1][1]
    b_y = Entrada["POINTS"][curve[2]-1][2]
    b_x = Entrada["POINTS"][curve[2]-1][1]

    d = sqrt((a_y-b_y)**2+(a_x-b_x)**2)
    c = b_x - a_x
    c = c/d

    s = b_y - a_y
    s = s/d

    return  d, c, s

#To precisando da Rigidez. M

def material_property(property_id:int, Entrada:list):
    A = Entrada["PROPERTIES"][property_id-1][2]
    I = Entrada["PROPERTIES"][property_id-1][3]
    material_id = Entrada["PROPERTIES"][property_id-1][1]
    E =  Entrada["MATERIALS"][material_id-1][1]
    mu =  Entrada["MATERIALS"][material_id-1][2]

    return A, I, E, mu

tamanho_matriz =0
for i in range(len((Entrada["MESH"]))): 
    tamanho_matriz += Entrada["MESH"][i][2]

tamanho_matriz = int(tamanho_matriz) + 1
Entrada["TAMANHO_MATRIZ"] = 3*tamanho_matriz

def matriz_rigidez(curve_id:int, Entrada:list):

    point_1 = Entrada["CURVES"][curve_id][1]-1
    point_2 = Entrada["CURVES"][curve_id][2]-1
    L, c,s = geometria_curva(curve_id, Entrada)
    
    tamanho_matriz = Entrada["TAMANHO_MATRIZ"]
    K_temp = np.zeros((tamanho_matriz,tamanho_matriz))
    property_id = Entrada["MESH"][curve_id][1]
    A, I, E, mu = material_property(property_id, Entrada)

    cte = A*L**2/(2*I)

    K_temp[3*point_1][3*point_1] = cte*c**2 + 6*s**2
    K_temp[3*point_1 + 1][3*point_1] = c*s*(cte -6)
    K_temp[3*point_1][3*point_1 + 1] = c*s*(cte -6)
    K_temp[3*point_1 + 1][3*point_1+1] = cte*s**2 + 6*c**2

    K_temp[3*point_1][3*point_2] = -(cte*c**2 + 6*s**2)
    K_temp[3*point_1 + 1][3*point_2] = -(c*s*(cte -6))
    K_temp[3*point_1][3*point_2 + 1] = -(c*s*(cte -6))
    K_temp[3*point_1 + 1][3*point_2+1] = -(cte*s**2 + 6*c**2)

    K_temp[3*point_2][3*point_1] = -(cte*c**2 + 6*s**2)
    K_temp[3*point_2 + 1][3*point_1] = -(c*s*(cte -6))
    K_temp[3*point_2][3*point_1 + 1] = -(c*s*(cte -6))
    K_temp[3*point_2 + 1][3*point_1+1] = -(cte*s**2 + 6*c**2)

    K_temp[3*point_2][3*point_2] = cte*c**2 + 6*s**2
    K_temp[3*point_2 + 1][3*point_2] = c*s*(cte -6)
    K_temp[3*point_2][3*point_2 + 1] = c*s*(cte -6)
    K_temp[3*point_2 + 1][3*point_2+1] = cte*s**2 + 6*c**2

    #---
    K_temp[3*point_1][3*point_1+2] = -3*L*s
    K_temp[3*point_1+2][3*point_1] = -3*L*s

    K_temp[3*point_1+1][3*point_1+2] = 3*L*c
    K_temp[3*point_1+2][3*point_1+1] = 3*L*c

    K_temp[3*point_1+2][3*point_1+2] = 2*L**2

    #---
    K_temp[3*point_2][3*point_1+2] = 3*L*s
    K_temp[3*point_2+2][3*point_1] = -3*L*s

    K_temp[3*point_2+1][3*point_1+2] = -3*L*c
    K_temp[3*point_2+2][3*point_1+1] = 3*L*c

    K_temp[3*point_2+2][3*point_1+2] = L**2

    #---
    K_temp[3*point_1][3*point_2+2] = -3*L*s
    K_temp[3*point_1+2][3*point_2] = 3*L*s

    K_temp[3*point_1+1][3*point_2+2] = 3*L*c
    K_temp[3*point_1+2][3*point_2+1] = -3*L*c

    K_temp[3*point_1+2][3*point_2+2] = L**2

    #---
    K_temp[3*point_2][3*point_2+2] = +3*L*s
    K_temp[3*point_2+2][3*point_2] = +3*L*s

    K_temp[3*point_2+1][3*point_2+2] = -3*L*c
    K_temp[3*point_2+2][3*point_2+1] = -3*L*c

    K_temp[3*point_2+2][3*point_2+2] = 2*L**2

    K_temp = K_temp*2*E*I/L**3
    return K_temp

#FUNCIONA PARA APENAS UMA MESH

tamanho_matriz = Entrada["TAMANHO_MATRIZ"]
K = np.zeros((tamanho_matriz,tamanho_matriz))

for i in range(len(Entrada["CURVES"])):
    K += matriz_rigidez(i,Entrada)

#Força distribuida

def forca_distribuida(force_dist_id:int, Entrada:list):
    direcao = Entrada["DIST_LOADS"][force_dist_id][4]
    curve_id = Entrada["DIST_LOADS"][force_dist_id][1] - 1
    L, c,s = geometria_curva(curve_id, Entrada)

    point_1 = Entrada["CURVES"][curve_id][1]-1
    point_2 = Entrada["CURVES"][curve_id][2]-1

    tamanho_matriz = Entrada["TAMANHO_MATRIZ"]
    f_temp = np.zeros(tamanho_matriz)

    if direcao == 'x':
        q_1 = Entrada["DIST_LOADS"][force_dist_id][2]
        q_2 = Entrada["DIST_LOADS"][force_dist_id][3]
        f_temp[3*point_1] = L*(2*q_1+q_2)/6
        f_temp[3*point_2] = L*(q_1+2*q_2)/6
    
    if direcao == 'y':
        q_1 = Entrada["DIST_LOADS"][force_dist_id][2]
        q_2 = Entrada["DIST_LOADS"][force_dist_id][3]
        f_temp[3*point_1+1] = L*(7*q_1+3*q_2)/20
        f_temp[3*point_1+2] = L**2*(3*q_1+2*q_2)/60
        
        f_temp[3*point_2+1] = L*(3*q_1+7*q_2)/20
        f_temp[3*point_2+2] = -L**2*(2*q_1+3*q_2)/60

    #   FALTA MUDAR A DIREÇÃO 
    return f_temp

f = np.zeros(tamanho_matriz)
for i in range(len(Entrada["DIST_LOADS"])):
    f += forca_distribuida(i,Entrada)

#Força Concetrada
def forca_concentrada(force_conc_id:int, Entrada:list):
    point_id = Entrada["POINT_LOADS"][force_conc_id][1] -1 
    value = Entrada["POINT_LOADS"][force_conc_id][3] 
    related_gdl = Entrada["POINT_LOADS"][force_conc_id][2] 
    
    tamanho_matriz = Entrada["TAMANHO_MATRIZ"]
    f_conc_temp = np.zeros(tamanho_matriz)

    if related_gdl == 1:
        f_conc_temp[3*point_id] = value
    if related_gdl == 2:
        f_conc_temp[3*point_id+1] = value
    if related_gdl == 3:
        f_conc_temp[3*point_id+2] = value
 
    return f_conc_temp

b = np.zeros(tamanho_matriz)
for i in range(len(Entrada["POINT_LOADS"])):
    b += forca_concentrada(i,Entrada)

#Condições de Contorno

def boudary_condition(BC_id: int, Entrada: list, BC:list, forca: list):
    
    point_id = Entrada["BC"][BC_id][0] -1 

    value = Entrada["BC"][BC_id][2] 
    related_gdl = Entrada["BC"][BC_id][1] 

    if related_gdl == 1:
        for i in range(tamanho_matriz):
            BC[3*point_id][i] = 0
        for i in range(tamanho_matriz):
            BC[i][3*point_id] = 0
        BC[3*point_id][3*point_id] = 1
        forca[3*point_id] = 0
        
    if related_gdl == 2:
        for i in range(tamanho_matriz):
            BC[3*point_id+1][i] = 0
        for i in range(tamanho_matriz):
            BC[i][3*point_id+1] = 0
        BC[3*point_id+1][3*point_id+1] = 1
        forca[3*point_id+1] = 0
    
    if related_gdl == 3:
        for i in range(tamanho_matriz):
            BC[3*point_id+2][i] = 0
        for i in range(tamanho_matriz):
            BC[i][3*point_id+2] = 0
        BC[3*point_id+2][3*point_id+2] = 1
        forca[3*point_id+2] = 0

    return BC, forca

forca_final = b + f
K_star = K.copy()
force_star = forca_final.copy()
for i in range(len(Entrada["BC"])):
    K_star, force_star = boudary_condition(i,Entrada,K_star, force_star)

u = np.linalg.solve(K_star, force_star)
reactions_force = K@u - forca_final

# PEGAR O ELEMENTO DA CURVA
# CALCULAR OS ESFORÇOS DAQUELE ELEMENTO
# COLOCAR EM UM ARRAY

N = [] #PONTO INICIAL E FINAL
V = []
M = []

for curve_id in range(len(Entrada["CURVES"])):
    point_1 = Entrada["CURVES"][curve_id][1]-1
    point_2 = Entrada["CURVES"][curve_id][2]-1
    L, c,s = geometria_curva(curve_id, Entrada)

    property_id = Entrada["MESH"][curve_id][1]
    A, I, E, mu = material_property(property_id, Entrada)

    R_T = np.array([
        [ c, s,0],
        [-s, c,0],
        [ 0, 0,1]
    ])
    u_local_1 = R_T@u[[3*point_1,3*point_1+1,3*point_1+2]]
    u_local_2 = R_T@u[[3*point_2,3*point_2+1,3*point_2+2]]
    
    N.append([E*A/L*(u_local_2[0] - u_local_1[0]),
              E*A/L*(u_local_2[0] - u_local_1[0])])
    V.append([E*I*( 12/L**3*u_local_1[1] + 6/L**2*u_local_1[2] - 12/L**3*u_local_2[1] + 6/L**2*u_local_2[2]),
              E*I*( 12/L**3*u_local_1[1] + 6/L**2*u_local_1[2] - 12/L**3*u_local_2[1] + 6/L**2*u_local_2[2])])
    M.append([E*I*(-6/L**2*u_local_1[1] - 4/L*u_local_1[2] + 6/L**2*u_local_2[1] - 2/L*u_local_2[2]),
              E*I*( 6/L**2*u_local_1[1] + 2/L*u_local_1[2] - 6/L**2*u_local_2[1] + 4/L*u_local_2[2])])

#Escrita da saída

with open("saida_viga.txt", "w", encoding="utf-8") as arquivo:
    arquivo.write("Resultado da Simulação\n\n")
    arquivo.write("Arquivo de Entrada: entrada_Ex_Portico.txt\n\n\n")

    arquivo.write("-"*18 + " Resultados Nodais " + "-"*18 + "\n")
    arquivo.write("|%5s|%15s|%15s|%15s|\n" % ("Nó", "u", "v", "theta"))
    arquivo.write("-"*55 + "\n")
    for no in range(len(u)//3):
        arquivo.write("|%5d|%15s|%15s|%15s|\n" % (
            no+1,
            str(round(float(u[3*no]), 6)),
            str(round(float(u[3*no+1]), 6)),
            str(round(float(u[3*no+2]), 6))))
    arquivo.write("-"*55 + "\n\n\n")

    arquivo.write("-"*54 + " Resultados por Elemento " + "-"*54 + "\n")
    arquivo.write("|%9s|%12s|%15s|%15s|%15s|%12s|%15s|%15s|%15s|\n" % (
        "Elemento", "Nó Inicial", "Força Axial", "Cisalhamento", "Momento Fletor",
        "Nó Final", "Força Axial", "Cisalhamento", "Momento Fletor"))
    arquivo.write("-"*133 + "\n")
    for i in range(len(Entrada["CURVES"])):
        no_ini = Entrada["CURVES"][i][1]
        no_fim = Entrada["CURVES"][i][2]
        arquivo.write("|%9d|%12d|%15s|%15s|%15s|%12d|%15s|%15s|%15s|\n" % (
            i+1, no_ini,
            str(round(float(N[i][0]), 3)), str(round(float(V[i][0]), 3)), str(round(float(M[i][0]), 3)),
            no_fim,
            str(round(float(N[i][1]), 3)), str(round(float(V[i][1]), 3)), str(round(float(M[i][1]), 3))))
    arquivo.write("-"*133 + "\n\n\n")

    arquivo.write("-"*6 + " Forças de Reação " + "-"*6 + "\n")
    arquivo.write("|%5s|%6s|%15s|\n" % ("Nó", "GDL", "Value"))
    arquivo.write("-"*30 + "\n")
    for no in range(len(Entrada["BC"])):
        point_id = int(Entrada["BC"][no][0])
        gdl = int(Entrada["BC"][no][1])
        dof = 3*(point_id-1) + (gdl-1)
        arquivo.write("|%5d|%6d|%15s|\n" % (point_id, gdl, str(round(float(reactions_force[dof]), 4))))
    arquivo.write("-"*30 + "\n")
    

















# f_local = np.zeros(6)
    
#     if direcao == 'l': # Carga Axial Local
#         q_1 = Entrada["DIST_LOADS"][force_dist_id][2]
#         q_2 = Entrada["DIST_LOADS"][force_dist_id][3]
#         f_local[0] = L*(2*q_1+q_2)/6
#         f_local[3] = L*(q_1+2*q_2)/6
    
#     elif direcao == 't': # Carga Transversal Local
#         p_1 = Entrada["DIST_LOADS"][force_dist_id][2]
#         p_2 = Entrada["DIST_LOADS"][force_dist_id][3]
#         f_local[1] = L*(7*p_1+3*p_2)/20
#         f_local[2] = L**2*(3*p_1+2*p_2)/60
#         f_local[4] = L*(3*p_1+7*p_2)/20
#         f_local[5] = -L**2*(2*p_1+3*p_2)/60

#     # Matriz de transformação R transposta
#     R_T = np.array([
#         [c, -s, 0, 0, 0, 0],
#         [s,  c, 0, 0, 0, 0],
#         [0,  0, 1, 0, 0, 0],
#         [0,  0, 0, c, -s, 0],
#         [0,  0, 0, s,  c, 0],
#         [0,  0, 0, 0,  0, 1]
#     ])

#     f_global_elem = R_T @ f_local

#     # Alocação no vetor global
#     f_temp[3*point_1 : 3*point_1+3] += f_global_elem[0:3]
#     f_temp[3*point_2 : 3*point_2+3] += f_global_elem[3:6]