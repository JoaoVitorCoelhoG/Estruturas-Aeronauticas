from math import sqrt
import numpy as np

def defactor(linha: str, f):
    """
    n_mesh:
    Entrada:
    """
    linha = next(f)
    a = []
    while linha.strip():
        row = []
        for elemento in linha.strip().split(','):
            try:
                row.append(float(elemento))
            except ValueError:
                row.append(elemento.strip())
        a.append(row)
        linha = next(f)
    return a

Entrada = {}

with open("entrada_Ex_Trelica.txt", "r", encoding="utf-8") as f:

    for linha in f:

        if linha == "*POINTS (id,x_point,y_point)\n":
            Entrada["POINTS"] = defactor(linha,f)

        if linha == "*CURVES (id,point_ini,point_fim)\n":
            Entrada["CURVES"] = defactor(linha,f)

        if linha == "*MATERIALS (id,E,nu)\n":
            Entrada["MATERIALS"] = defactor(linha,f)

        if linha == "*PROPERTIES (id,material_id,A)\n":
            Entrada["PROPERTIES"] = defactor(linha,f)

        if linha == "*MESH (curve_id,property_id,num_elem)\n":
            Entrada["MESH"] = defactor(linha,f)

        if linha == "*POINT_LOADS (id,point_id,related_gdl,value)\n":
            Entrada["POINT_LOADS"] = defactor(linha,f)

        if linha == "*DIST_LOADS (id,curve_id,value_ini,value_fim,direction:(x,y,l,t))\n":
            Entrada["DIST_LOADS"] = defactor(linha,f)

        if linha == "*BC (point_id,related_gdl,value)\n":
            Entrada["BC"] = defactor(linha,f)

tamanho_matriz =0
for i in range(len((Entrada["MESH"]))): 
    tamanho_matriz += Entrada["MESH"][i][2]

tamanho_matriz = int(tamanho_matriz) + 1

Entrada["TAMANHO_MATRIZ"] = 2*tamanho_matriz

def tamanho_mesh(n_mesh:int,Entrada:list)->float:
    """
    n_mesh: É o número do mesh definido pela ordem de meshs da entrada (Ordem da barra)
    Entrada: É o Dicionário de todos os dados da entrada
    return:: float
    return: O comprimento do n_mesh-esimo elemento.
    """

    curve = Entrada["CURVES"][int(Entrada["MESH"][n_mesh-1][0])-1]
    a_y = Entrada["POINTS"][int(curve[1])-1][2]
    a_x = Entrada["POINTS"][int(curve[1])-1][1]
    b_y = Entrada["POINTS"][int(curve[2])-1][2]
    b_x = Entrada["POINTS"][int(curve[2])-1][1]

    d = sqrt((a_y-b_y)**2+(a_x-b_x)**2)
    return  d

def angulo_mesh(n_mesh:int,Entrada:list)->float:
    """
    n_mesh: É o número do mesh definido pela ordem de meshs da entrada (Ordem da barra)
    Entrada: É o Dicionário de todos os dados da entrada
    return:: float
    return: O comprimento do n_mesh-esimo elemento.
    """

    curve = Entrada["CURVES"][int(Entrada["MESH"][n_mesh-1][0])-1]
    a_y = Entrada["POINTS"][int(curve[1])-1][2]
    a_x = Entrada["POINTS"][int(curve[1])-1][1]
    b_y = Entrada["POINTS"][int(curve[2])-1][2]
    b_x = Entrada["POINTS"][int(curve[2])-1][1]

    d = sqrt((a_y-b_y)**2+(a_x-b_x)**2)
    s= (b_y - a_y)/d
    c= (b_x - a_x)/d
    return  c,s

def phi_integral(curve_id:int, Entrada:list)->float: # Tá linear ainda
    """
    curve_id: é o id que representa a curva
    Entrada: É o Dicionário de todos os dados da entrada
    return:: float
    return: Vetor de força distribuída do n_mesh-esimo elemento.
    """

    mesh = 0
    count = 0
    while curve_id != int(Entrada["MESH"][mesh][0]):
        count+=int(Entrada["MESH"][mesh][2])
        mesh+=1

    element_in_mesh = int(Entrada["MESH"][mesh][2]) 
    cte = (tamanho_mesh(mesh,Entrada))/element_in_mesh

    tamanho_matriz = Entrada["TAMANHO_MATRIZ"] 
    b = np.zeros(tamanho_matriz) # Inicializando 
     
    point_force = 0
    while curve_id != int(Entrada["DIST_LOADS"][point_force-1][1]): #Caso a força não esteja alinhada com o vetor
        point_force+=1

    q_initial = Entrada["DIST_LOADS"][point_force-1][2] #Dá para remover
    q_final = Entrada["DIST_LOADS"][point_force-1][3]

    #vetor = np.linspace(curve_init[1], curve_end[1], element_in_mesh)
    q = np.linspace(q_initial, q_final, element_in_mesh+1)
    
    for k in range(element_in_mesh):
        q_initial = q[k]
        q_final = q[k+1]
        if(Entrada["DIST_LOADS"][point_force-1][4]=='x'):
            b[2*count] += (2*q_initial + q_final)/6*cte
            b[2*count+2] += (1*q_initial + 2*q_final)/6*cte
        if(Entrada["DIST_LOADS"][point_force-1][4]=='y'):
            b[2*count+1] += (2*q_initial + q_final)/6*cte
            b[2*count+3] += (1*q_initial + 2*q_final)/6*cte
        count+=1

    return b

def phi_derivate(n_mesh:int,node:int,Entrada:list) -> list:
    """
    n_mesh: É o número do mesh definido pela ordem de mesh's da entrada.
    Entrada: É o Dicionário de todos os dados da entrada
    return: float
    return: A Matriz de Rigidez do n_mesh-esimo elemento, para cada malha.
    """
    
    cte = property_area(n_mesh,Entrada)*property_elasticity(n_mesh,Entrada)*Entrada["MESH"][n_mesh-1][2]/(tamanho_mesh(n_mesh,Entrada))

    tamanho_matriz = Entrada["TAMANHO_MATRIZ"] 

    K = np.zeros((tamanho_matriz,tamanho_matriz))
    A = np.zeros((2,2))
    c, s = angulo_mesh(n_mesh,Entrada)
    A = cte*np.array([
        [ c**2,  c*s, -c**2, -c*s],
        [ c*s,   s**2, -c*s, -s**2],
        [-c**2, -c*s,  c**2,  c*s],
        [-c*s,  -s**2,  c*s,  s**2]
    ])
    for i in range(len(A)):
        for j in range(len(A)):
            K[node + i][node + j] = A[i][j]
    return K
            
def property_area(n_mesh: int, Entrada: list) ->float:
    """
    n_mesh: É o número do mesh definido pela ordem de meshs da entrada, começa do 0
    Entrada: É o Dicionário de todos os dados da entrada
    return:: float
    return: A área do n_mesh-esimo elemento.
    """
    A = Entrada["PROPERTIES"][int(Entrada["MESH"][n_mesh-1][1])-1][2] # Área do i-esimo elemento
    return A

def property_elasticity(i: int, Entrada:list) -> float:
    """
    n_mesh: É o número do mesh definido pela ordem de meshs da entrada
    Entrada: É o Dicionário de todos os dados da entrada
    return:: float
    return: E o módulo de elasticidade do n_mesh-esimo elemento.
    """
    E = Entrada["MATERIALS"][int(Entrada["PROPERTIES"][int(Entrada["MESH"][i-1][1])-1][1]) - 1][1]#Modulo de Young do i-esimo elemento
    return E

## CÁLCULO DA MATRIZ DE RIGIDEZ

tamanho_matriz = Entrada["TAMANHO_MATRIZ"] 
    
K = np.zeros((tamanho_matriz,tamanho_matriz))

node = 0
for curve in range(1,int(len(Entrada["MESH"]))+1):
    mesh = 1
    while curve != int(Entrada["MESH"][mesh-1][0]):
        mesh+=1
    count = 0
    while (count < Entrada["MESH"][mesh-1][2]):
        K += phi_derivate(mesh,2*node,Entrada) 
        node +=1
        count +=1
        
#print(K)

#Cargas Distribuídas

b = np.zeros(tamanho_matriz)

dist_curve_ids = {int(dl[1]) for dl in Entrada.get("DIST_LOADS", [])}
for curve in range(1, int(len(Entrada["MESH"]))+1):
    if curve in dist_curve_ids:
        b += phi_integral(curve, Entrada)

#print(b)

#Cargas Pontuais

F = np.zeros(tamanho_matriz)

for force_id in range(0,int(len(Entrada["POINT_LOADS"]))):
    curve = 0
    a=0
    while Entrada["POINT_LOADS"][force_id][1] != Entrada["CURVES"][curve][1]:
        curve+=1
    
    mesh = 0
    while curve+1 != Entrada["MESH"][mesh][0]:
        mesh+=1        
        a += int(Entrada["MESH"][mesh][2])
    gdl = int(Entrada["POINT_LOADS"][force_id][2])
    F[2*a + (gdl - 1)] += Entrada["POINT_LOADS"][force_id][3]

#print(F)

Forca_final = F+b

def restricao(K,ponto:int,Forca_Final):
    
    A = K.copy()
    for i in range(len(K)):
        for j in range(len(K)):
            if i == ponto or j == ponto:
                A[i,j] = 0
    A[ponto][ponto] = 1
    b = Forca_Final.copy()
    b[ponto] = 0
    return A, b

Kcc = K.copy()
Fcc = Forca_final.copy()
for restriction_point in range(0, int(len(Entrada["BC"]))):
    point_id = int(Entrada["BC"][restriction_point][0])
    gdl      = int(Entrada["BC"][restriction_point][1])
    dof = 2 * (point_id - 1) + (gdl - 1)
    Kcc, Fcc = restricao(Kcc, dof, Fcc)

for i in range(tamanho_matriz):
    if Kcc[i, i] == 0:
        Kcc, Fcc = restricao(Kcc, i, Fcc)

u = np.linalg.solve(Kcc, Fcc)

f_cc = Kcc@u
reactions_force = K@u - Forca_final

d = []
a = 0
for curve in range(len(Entrada["CURVES"])):
    mesh = 0
    while curve + 1 != Entrada["MESH"][mesh][0]:
        mesh += 1
    d.append(u[a])      
    d.append(u[a+1])    
    a += 2 * int(Entrada["MESH"][mesh][2])

d.append(u[a])      
d.append(u[a+1])    

N = []
for curve in range(len(Entrada["CURVES"])):
    mesh = 0
    while curve + 1 != Entrada["MESH"][mesh][0]:
        mesh += 1
    c, s = angulo_mesh(mesh+1, Entrada)
    cte = property_elasticity(mesh+1, Entrada) * property_area(mesh+1, Entrada) / tamanho_mesh(mesh+1, Entrada)
    dux = d[2*curve+2]   - d[2*curve]
    duy = d[2*curve+3] - d[2*curve+1]
    N.append(cte * (c*dux + s*duy))

with open("saida_barra.txt", "w") as arquivo:
    arquivo.write("Resultado da Simulação\n\n")
    
    arquivo.write("Arquivo de Entrada: entrada_Ex_Trelica.txt\n\n\n")
    
    
    arquivo.write("---------- Resultados Nodais ----------\n")
    arquivo.write("|   Nó|              u|              v|\n")
    arquivo.write("---------------------------------------\n")

    for no in range(len(d)//2):
        arquivo.write("|%5d|%15.4f|%15.4f|\n" % (no+1, d[2*no], d[2*no+1]))

    arquivo.write("---------------------------------------\n\n\n")


    arquivo.write("---- Resultados por Elemento ----\n")
    arquivo.write("|     Elemento|          Força x|\n")
    arquivo.write("---------------------------------\n")

    for no in range(len(N)):
        arquivo.write("|%13d|%17.2f|\n" % (no+1, N[no]))
    arquivo.write("---------------------------------\n\n\n")

    arquivo.write("------ Forças de Reação -----\n")
    arquivo.write("|   Nó|   GDL|          Value|\n")
    arquivo.write("-----------------------------\n")

    for no in range(len(Entrada["BC"])):
        point_id = int(Entrada["BC"][no][0])
        gdl      = int(Entrada["BC"][no][1])
        dof = 2 * (point_id - 1) + (gdl - 1)
        arquivo.write("|%5d|%6d|%15.2f|\n" % (point_id, gdl, reactions_force[dof]))
