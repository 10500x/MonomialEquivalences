import numpy as np
import itertools
import galois
from sympy.polys.fields import xfield
from sympy import symbols
from sympy.polys.domains import GF as SympyGF

################### 1 
################### 1 
################### 1 
def Crear_cuerpo(q):
    """
    Crea el cuerpo finito GF(q) usando la biblioteca galois.
        GF: la clase de campo finito GF(q).
    """
    GF = galois.GF(q)
    return GF
#################### FIN 1
#################### FIN 1
#################### FIN 1




################### 2 
################### 2 
################### 2 
def generar_PMI2(q):
    """ Genera polinomios irreducibles monicos de grado 2 sobre GF(q).
    Devuelve una lista de polinomios irreducibles (objetos Poly de galois) de grado 2."""
    
    PMI2 = list(galois.irreducible_polys(q, 2))
    return PMI2

def generar_GL2(GF):
    """ Genera todas las matrices 2x2 invertibles sobre GF.
    Devuelve una lista de arrays numpy de dtype GF."""
    
    elementos = GF.elements.tolist()
    GL2 = []
    for a, b, c, d in itertools.product(elementos, repeat=4):
        M = np.array([[a, b], [c, d]], dtype=object)
        det = M[0,0]*M[1,1] - M[0,1]*M[1,0]
        if det != GF(0):
            # ignoramos escalares multiplicativos
            GL2.append(M)
    return GL2

def generar_PGL2(GL2, GF):
    """
    Genera representantes únicos de las clases de equivalencia de PGL(2,q),
    identificando matrices que difieren por un escalar no nulo.
    """
    PGL = []
    vistos = []

    for M in GL2:
        # Buscamos el primer elemento no nulo
        for val in M.flatten():
            if val != GF(0):
                escalar = val
                break
        else:
            continue  # (no debería pasar, GL2 solo tiene matrices invertibles)

        # Normalizamos
        M_norm = M / escalar

        # Evitamos duplicados (matrices equivalentes)
        if not any(np.all(M_norm == N) for N in vistos):
            vistos.append(M_norm)
            PGL.append(M)

    return PGL
################# Esto se usa solo para IMS
def aplica_transformacion(M, z, GF):
    """
    Aplica la transformación lineal fraccional definida por M a un punto z.
    z = None representa el punto en el infinito.
    """
    a, b = M[0,0], M[0,1]
    c, d = M[1,0], M[1,1]
    if z is None:
        # M(inf) = a/c si c != 0, sino inf
        return None if c == GF(0) else a / c
    # z es finito en GF
    denom = ((c * z) + GF(d))
    if denom == GF(0):
        return None  # mapea a infinito
    return ( ( (a * z) + GF(b) ) / denom)

def fija_lugar_racional(GF,IMS):
    """
    Verifica si M fija algún lugar racional (alfa o infinito).
    Retorna True si existe alfa tal que M(alfa)=alfa o M(inf)=inf.
    """
    for M in IMS:
        # Verificar infinito
        inf_image = aplica_transformacion(M, None, GF)
        if inf_image is None:
            return True
        # Probar cada elemento finito
        for z in GF.elements:
            if aplica_transformacion(M, GF(z), GF) == GF(z):
                return True
        return False
################# Esto se usa solo para IMS 

def generar_IMS(GF,GL2):
    """Para cada matriz en GL2, nos fijamos si fija lugar racional, si no, lo metemos en GL2_Fq"""
    IMS=[]
    for matriz in GL2:
       if not(fija_lugar_racional(GF,IMS)):
           IMS.append(matriz)
    return IMS
########################### FIN 2
########################### FIN 2
########################### FIN 2       


################### 3 

                
                
##################### Fin 3


################### 5

def orden_matriz(PGL, q):
    """
    Calcula el orden de cada matriz en PGL(2,q).
    """
    GF = galois.GF(q)
    matord = []
    for matriz in PGL:
        n = 1
        while True:
            M_n = np.linalg.matrix_power(matriz.astype(object), n)# Elevamos matriz^n
            I = np.array([[GF(1), GF(0)], [GF(0), GF(1)]], dtype=object) #identidad del cuerpo
            escalar = M_n[0,0] # Escalar candidato
            if np.array_equal(M_n, escalar * I):# Comparamos en el cuerpo
                matord.append(n)
                break
            else:
                n += 1
                if n > q+1:
                    break

    return matord


        

############# FIN 5

def crear_cuerpo_funciones_sympy(p, m):
    """
    Crea simbólicamente el campo de funciones GF(p^m)(x) usando Sympy.
    Esto crea el dominio base y el generador a.
    (Nota: Sympy.GF solo acepta primos, así que usamos GF(primo).)
    """
    x = symbols('x')
    # Simpy GF para primos solo; si q = p^m con m>1, podríamos considerar GF(p) como base.
    domain = SympyGF(p)
    field, x = xfield(x, domain)  # xfield retorna (campo, generador)
    print(f"Campo de funciones (simbólico): {field}, generador {x}")
    return field, x






def calcular_orbitas(M, GF):
    """
    Calcula las órbitas de los puntos de la recta proyectiva GF(q) u {∞} bajo M.
    Retorna una lista de órbitas.
    """
    lugares = [GF(a) for a in GF.elements] + [None]  # None representa ∞
    orbitas = []
    vistos = set()
    for z in lugares:
        if z in vistos:
            continue
        orb = []
        current = z
        while current not in orb:
            orb.append(current)
            vistos.add(current)
            current = aplica_transformacion(M, current, GF)
        orbitas.append(orb)
    return orbitas

def construir_codigo_sigma(orbitas, divisor=None):
    """
    Construye un código sigma-cíclico dado un conjunto de órbitas y un divisor.
    """
    # Ejemplo de salida: 
    return [len(orb) for orb in orbitas]


def son_equivalentes_perm(M1, M2):
    """Verifica si dos matrices generadoras M1 y M2 son equivalentes por permutación de columnas."""
    return False  # Placeholder

def agrupar_codigos_equivalentes(codigos):
    """
    Agrupa códigos bajo equivalencia de permutación/monomio.
    Parámetro:
        codigos: lista de códigos (cada uno puede representarse por su matriz generadora).
    Retorna:
        lista de grupos, donde cada grupo contiene códigos equivalentes entre sí.
    """
    grupos = []
    for C in codigos:
        colocado = False
        for g in grupos:
            if son_equivalentes_perm(g[0], C):
                g.append(C)
                colocado = True
                break
        if not colocado:
            grupos.append([C])
    return grupos

def long_sigma(q,matord):
    N = []
    for n in set(matord):
        if (n < q + 1) and (n > 3) and ((q + 1) % n == 0):
            N.append(n)
    print(31*'-')
    return N
