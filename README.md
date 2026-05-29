# Estruturas Aeronáuticas — EST-40

Implementação em Python do **Método dos Elementos Finitos (MEF)** para análise estática de estruturas reticuladas planas (barras, vigas e pórticos 2D), desenvolvida para a disciplina **EST-40 — Estruturas Aeronáuticas**.

O programa monta a matriz de rigidez global a partir de uma descrição textual da geometria, propriedades e carregamentos, aplica as condições de contorno e resolve o sistema linear $K \cdot u = f$, fornecendo deslocamentos nodais, esforços internos por elemento (normal, cortante, momento fletor) e reações de apoio.

---

## Estrutura do repositório

```
ESTRUTURAS AERONÁUTICAS/
└── MÉTODOS DOS ELEMENTOS FINITOS/
    ├── LAB_1_EST40.py           # Código-fonte do solver
    ├── entrada_Ex3_Lista2.txt   # Exemplo: viga apoiada + mola
    ├── entrada_Ex_Portico.txt   # Exemplo: pórtico plano biengastado
    ├── entrada_Ex_Viga.txt      # Exemplo: viga contínua
    ├── saida_Ex3_Lista2.txt     # Saída correspondente
    ├── saida_Ex_Portico.txt
    └── saida_Ex_Viga.txt
```

---

## Pré-requisitos

- Python ≥ 3.8
- [NumPy](https://numpy.org/)

```bash
pip install numpy
```

---

## Como executar

1. Entre na pasta do solver:

   ```bash
   cd "MÉTODOS DOS ELEMENTOS FINITOS"
   ```

2. Selecione o arquivo de entrada desejado editando a linha em [LAB_1_EST40.py:31](MÉTODOS%20DOS%20ELEMENTOS%20FINITOS/LAB_1_EST40.py#L31):

   ```python
   arquivo_entrada = "entrada_Ex3_Lista2.txt"
   ```

3. Rode o script:

   ```bash
   python LAB_1_EST40.py
   ```

4. O arquivo de saída é gerado automaticamente com o nome `saida_<...>.txt` (mesma raiz do arquivo de entrada).

---

## Formato do arquivo de entrada

O arquivo de entrada é dividido em blocos, cada um iniciado por um cabeçalho `*NOME (campos)`. Os blocos são lidos em qualquer ordem, mas o conteúdo de cada bloco deve seguir exatamente os campos indicados no cabeçalho. Linhas em branco terminam o bloco.

| Bloco | Descrição | Campos |
|---|---|---|
| `*POINTS` | Pontos geométricos | `id, x, y` |
| `*CURVES` | Curvas (barras) conectando dois pontos | `id, point_ini, point_fim` |
| `*MATERIALS` | Materiais | `id, E, nu` |
| `*PROPERTIES` | Seções transversais | `id, material_id, A, I` |
| `*MESH` | Discretização de cada curva | `curve_id, property_id, num_elem` |
| `*POINT_LOADS` | Cargas concentradas | `id, point_id, gdl, value` |
| `*DIST_LOADS` | Cargas distribuídas (lineares) | `id, curve_id, value_ini, value_fim, direção` |
| `*SPRING` | Molas aterradas | `point_id, gdl, k` |
| `*BC` | Condições de contorno (apoios) | `point_id, gdl, value` |

**Convenções:**
- `gdl = 1` → deslocamento em $x$
- `gdl = 2` → deslocamento em $y$
- `gdl = 3` → rotação $\theta$
- Direção da carga distribuída:
  - `x`, `y` → eixos globais
  - `l` (longitudinal), `t` (transversal) → eixos locais da barra

### Exemplo mínimo

```text
*POINTS (id,x_point,y_point)
1,0.0,0.0
2,5000.,0.0
3,10000.,0.0

*CURVES (id,point_ini,point_fim)
1,1,2
2,2,3

*MATERIALS (id,E,nu)
1,200000,0.3

*PROPERTIES (id,material_id,A,I)
1,1,109544.5,1000000000

*MESH (curve_id,property_id,num_elem)
1,1,1
2,1,1

*POINT_LOADS (id,point_id,related_gdl,value)
1,2,3,-1250000.

*DIST_LOADS (id,curve_id,value_ini,value_fim,direction:(x,y,l,t))
1,1,-1.,-1.,y

*SPRING (point_id,related_gdl,value)
3,2,20.

*BC (point_id,related_gdl,value)
1,1,0.
1,2,0.
1,3,0.
2,2,0.
```

---

## Saída

O arquivo de saída contém três tabelas:

1. **Resultados Nodais** — deslocamentos $u, v$ e rotação $\theta$ de cada nó.
2. **Resultados por Elemento** — esforço normal $N$, cortante $V$ e momento fletor $M$ nas extremidades inicial e final de cada elemento.
3. **Forças de Reação** — reações nos graus de liberdade restringidos.

Exemplo (recorte de [saida_Ex3_Lista2.txt](MÉTODOS%20DOS%20ELEMENTOS%20FINITOS/saida_Ex3_Lista2.txt)):

```
------------------ Resultados Nodais ------------------
|   Nó|              u|              v|          theta|
-------------------------------------------------------
|    1|           -0.0|           -0.0|            0.0|
|    2|           -0.0|            0.0|       -7.2e-05|
|    3|            0.0|      -0.879007|      -0.000228|
-------------------------------------------------------
```

---

## Visão geral do solver

O fluxo principal de [LAB_1_EST40.py](MÉTODOS%20DOS%20ELEMENTOS%20FINITOS/LAB_1_EST40.py) é:

1. **Leitura** — `defactor()` faz o parsing tipado de cada bloco do arquivo de entrada.
2. **Geração da malha** — cada curva é subdividida em `num_elem` elementos; nós compartilhados entre curvas são reutilizados via `POINT_TO_NODE`.
3. **Montagem da rigidez global** — `matriz_rigidez()` monta a contribuição de cada elemento de pórtico plano (3 GDL por nó: $u$, $v$, $\theta$) já no referencial global usando os cossenos diretores $c, s$. Molas aterradas são somadas na diagonal.
4. **Vetor de forças** — `forca_distribuida()` aplica cargas trapezoidais nos eixos $x$, $y$, $l$ ou $t$; `forca_concentrada()` aplica cargas pontuais.
5. **Condições de contorno** — `boudary_condition()` aplica restrições zerando linha/coluna e fixando a diagonal.
6. **Solução** — `numpy.linalg.solve(K*, f*)` fornece os deslocamentos; as reações são recuperadas por $R = K u - f$.
7. **Pós-processamento** — esforços internos são calculados em coordenadas locais elemento a elemento.

---

## Limitações conhecidas

- Apenas estruturas planas (2D).
- Elemento de pórtico Euler-Bernoulli (sem cisalhamento).
- Cargas distribuídas trapezoidais lineares.
- O coeficiente de Poisson $\nu$ é lido mas não usado (reservado para futura extensão).

---

## Autor

Desenvolvido para o **LAB 1 de EST-40 — Estruturas Aeronáuticas**.
