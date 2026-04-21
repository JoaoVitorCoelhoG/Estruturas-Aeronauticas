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

with open("entrada_Ex_Barra.txt", "r", encoding="utf-8") as f:

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

Entrada["TAMANHO_MATRIZ"] = tamanho_matriz


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

def tamanho_curve(curve_id:int,Entrada:list)->float:
    """
    n_mesh: É o número do mesh definido pela ordem de meshs da entrada (Ordem da barra)
    Entrada: É o Dicionário de todos os dados da entrada
    return:: float
    return: O comprimento do n_mesh-esimo elemento.
    """

    curve = Entrada["CURVES"][curve_id-1]
    a_y = Entrada["POINTS"][int(curve[1])-1][2]
    a_x = Entrada["POINTS"][int(curve[1])-1][1]
    b_y = Entrada["POINTS"][int(curve[2])-1][2]
    b_x = Entrada["POINTS"][int(curve[2])-1][1]

    d = sqrt((a_y-b_y)**2+(a_x-b_x)**2)
    return  d

def phi(x:float,n_mesh:int,node:int, Entrada:list)->float: # Tá linear ainda
    """
    x: é o ponto que queremos analisar a função
    n_mesh: É o número do mesh definido pela ordem de meshs da entrada
    node: É o nó definido pela barra, podendo ser n_mesh -1 ou n_mesh.
    point: 
    Entrada: É o Dicionário de todos os dados da entrada
    return:: float
    return: A área do n_mesh-esimo elemento.
    """
    
    #Vou desconsiderar os pontos em y. Barra linear
    #Curva Referente a posição i
    point_init = int(Entrada["CURVES"][int(Entrada["MESH"][n_mesh-1][0])-1][1])
    x_init = Entrada["POINTS"][point_init-1][1]
    if(node==n_mesh-1):
        a = 1 - (x - x_init)/(tamanho_mesh(n_mesh,Entrada))
        return a
    if node == n_mesh:
        a = (x - x_init)/(tamanho_mesh(n_mesh, Entrada))
        return a
    return 0

def phi_integral(curve_id:int, Entrada:list)->float: # Tá linear ainda
    """
    curve_id: é o id que representa a curva
    Entrada: É o Dicionário de todos os dados da entrada
    return:: float
    return: Vetor de força distribuída do n_mesh-esimo elemento.
    """

    mesh = 0
    while curve_id != int(Entrada["MESH"][mesh-1][0]):
        mesh+=1

    element_in_mesh = int(Entrada["MESH"][mesh-1][2]) 
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
    

    mesh = 1
    count = 0
    while curve_id != int(Entrada["MESH"][mesh-1][0]):
        mesh+=1
        a = 0
        while a < int(Entrada["MESH"][mesh-1][2]):
            a +=1
        count+=a
    for k in range(element_in_mesh):
        q_initial = q[k]
        q_final = q[k+1]
        b[count] += (2*q_initial + q_final)/6*cte
        b[count+1] += (1*q_initial + 2*q_final)/6*cte
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
    K[node-1][node-1]=cte
    K[node][node]=cte
    K[node-1][node]=-cte
    K[node][node-1]=-cte
    # for i in range(tamanho_matriz):
    #     for j in range(tamanho_matriz): # Dá para colocar condição de para de acordo com o número de meshs
    #         if((i==node-1 and j == node-1) or (i==node and j == node)):
    #             K[i][j] = cte
    #         if((i==node and j == node-1) or (i==node-1 and j == node)):
    #             K[i][j] = -cte #COLOCAR ISSO EM UMA ENTRADA SÓ
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
        node +=1
        count +=1
        K += phi_derivate(mesh,node,Entrada) 
        
#print(K)

#Cargas Distribuídas

b = np.zeros(tamanho_matriz)

for curve in range(1,int(len(Entrada["MESH"]))+1):
    b+=phi_integral(curve,Entrada)

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
    F[a] += Entrada["POINT_LOADS"][force_id][3]

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

    
for restriction_point in range(0,int(len(Entrada["BC"]))):
    curve = 0
    while Entrada["BC"][restriction_point][0] != Entrada["CURVES"][curve][1]:
        if curve + 1 >= len(Entrada["CURVES"]) and Entrada["BC"][restriction_point][0] == Entrada["CURVES"][curve][2]:
            curve+=1
            break
        curve+=1
    Kcc, Fcc = restricao(K,curve, Forca_final)

u = np.linalg.solve(Kcc, Fcc)
#print(u)

f_cc = Kcc@u
reactions_force = K@u - Forca_final
d = u[::int(Entrada["MESH"][0][2])] #Isso não está correto, só funciona para quando o número de meshs é igual
#print(u[::int(Entrada["MESH"][0][2])])

#Falta só calcular as forças normais

N = np.zeros(tamanho_matriz - 1)

#Fazer generalizações para quando tiver muitas meshs
N = np.ones(tamanho_matriz - 1)
a =0 
for curve in range(len(Entrada["CURVES"])):
    mesh = 0
    while curve + 1 != Entrada["MESH"][mesh][0]:
        mesh +=1
    for i in range(int(Entrada["MESH"][mesh][2]) -1):
        a+=1
        N[a] = property_elasticity(mesh+1,Entrada)*property_area(mesh+1,Entrada)/tamanho_mesh(mesh,Entrada)*(u[a+1] - u[a])#problema nas meshs também

# Forma recomendada (fecha automaticamente)
with open("saida_barra.txt", "w") as arquivo:
    arquivo.write("Resultado da Simulação\n\n")
    
    arquivo.write("Arquivo de Entrada: entrada_Ex_Trelica.txt\n\n\n")
    
    
    arquivo.write("---------- Resultados Nodais ----------\n")
    arquivo.write("|   Nó|              u|              v|\n")
    arquivo.write("---------------------------------------\n")

    for no in range(tamanho_matriz):
        arquivo.write("|%5d|%15.4f|%15s|\n" % (no+1, u[no], "0.0"))

    arquivo.write("---------------------------------------\n\n\n")


    arquivo.write("---- Resultados por Elemento ----\n")
    arquivo.write("|     Elemento|          Força x|\n")
    arquivo.write("---------------------------------\n")

    for no in range(tamanho_matriz-1):
        arquivo.write("|%13d|%17.1f|\n" % (no+1, N[no]))
    arquivo.write("---------------------------------\n\n\n")

    arquivo.write("------ Forças de Reação -----\n")
    arquivo.write("|   Nó|   GDL|          Value|\n")
    arquivo.write("-----------------------------\n")

    for no in range(len(Entrada["BC"])):
        arquivo.write("|%13d|%14.1f|\n" % (Entrada["BC"][no][0], reactions_force[no]))
