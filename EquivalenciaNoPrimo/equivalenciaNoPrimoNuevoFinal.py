import sage.all
from sage.parallel.decorate import parallel
import time
from itertools import product

# ==============================================================================
# CONFIGURACIÓN GLOBAL
# ==============================================================================
# Lista de potencias de primos (q) a analizar.
# El código iterará sobre cada uno de estos valores para encontrar códigos sigma-cíclicos.
muestra = [9]

# ==============================================================================
# 1. EL WORKER (PROCESO PARALELO)
#    Esta función realiza el trabajo pesado de hallar los codigos.
#    Se ejecuta en paralelo en múltiples núcleos del CPU.
# ==============================================================================
@parallel
def worker_analizar_matriz(B_list, n, q, monicirred2_coeffs, elem_coeffs, menoselem_coeffs):
    """
    Busca códigos sigma-cíclicos asociados a una matriz B específica.
    """
    
    # --- RECONSTRUCCIÓN DEL ENTORNO MATEMÁTICO ---
    # Al ser un proceso paralelo independiente, no tiene acceso a las variables globales.
    # Debemos reconstruir los objetos algebraicos (cuerpos, anillos) aquí mismo.
    
    k = GF(q, 'a'); a = k.gen()            # Cuerpo finito GF(q) con generador 'a'
    R = PolynomialRing(k, 'x')             # Anillo de polinomios k[x]
    K = FunctionField(k, 'x'); x = K.gen() # Cuerpo de funciones racionales k(x)
    O = K.maximal_order()                  # Anillo de enteros del cuerpo de funciones
    
    # Obtenemos los Lugares Finitos del cuerpo de funciones.
    # Son los puntos sobre los cuales evaluaremos las funciones para generar el código.
    # Usamos la función nativa para asegurar el orden canónico de Sage.
    places_finite = K.places_finite()
    
    # Reconstruimos las listas de elementos y polinomios a partir de los coeficientes recibidos.
    elem = [k(c) for c in elem_coeffs]     
    menoselem = [k(c) for c in menoselem_coeffs]
    monicirred2 = [R(coeffs) for coeffs in monicirred2_coeffs]
    
    # Reconstruimos la matriz B (que define el automorfismo sigma).
    B = matrix(k, 2, B_list)
    # Calculamos la matriz adjunta para la transformación inversa.
    # (Necesaria para calcular las órbitas hacia atrás o verificar consistencia).
    Binv = matrix(k, 2, [B[1,1], -B[0,1], -B[1,0], B[0,0]])

    # --- DEFINICIÓN DEL AUTOMORFISMO SIGMA ---
    try:
        # Sigma es una Transformación de Möbius sobre el cuerpo de funciones.
        # sigma(x) = (ax + b) / (cx + d)
        sigma = K.hom([(B[0,0]*x + B[0,1]) / (B[1,0]*x + B[1,1])])
    except: return []

    # --- BÚSQUEDA DE LUGARES FIJOS ---
    # Buscamos polinomios irreducibles de grado 2 que sean "fijos" bajo sigma.
    # Un lugar P es fijo si sigma(P) = P (como conjuntos de ceros).
    # Algebraicamente: P debe dividir al numerador de sigma(P).
    lg2fijos = []
    for pol in monicirred2:
        im = sigma(pol)
        num = R(im.numerator())
        den = R(im.denominator())
        if num % pol == 0 and den % pol != 0:
            lg2fijos.append(pol)
    
    # Si no hay lugares fijos de grado 2, no podemos construir el código deseado.
    if not lg2fijos: return []

    # --- GENERACIÓN DE ÓRBITAS ---
    # Buscamos conjuntos de puntos {P1, ..., Pn} que formen un ciclo bajo sigma.
    orbitas = []
    for i in range(q):
        # Lógica de búsqueda de semilla:
        # Buscamos un punto inicial basado en la lista 'menoselem' para replicar
        # la lógica de permutación del algoritmo original.
        j = -1 #Asignamos un valor que no está en val o en enumerate(elem).
        for idx, val in enumerate(elem): #Si el valor que nos da el recorrer es igual a un menoselem[i],
            if val == menoselem[i]:       #asignamos al valor ese elemento.
                j = idx; break
        
        P_start = places_finite[j]
        
        
        ya_esta = False # Verificamos si este punto ya está en alguna órbita encontrada.
        for orb in orbitas:
            if P_start in orb: ya_esta = True; break
            
        if not ya_esta:         #Si no es cierto que ya está en alguna orbita,
            orbita = [P_start]  #es decir no está en ninguna orbita, construimos la orbita del punto.
            
            # Semilla de iteración:
            # Usamos potencias a^i para generar diversidad de puntos de inicio.
            if i == 0: val_aux = k(1) 
            else: val_aux = a**i
            aux = [val_aux]
            
            brk2 = 0
            # Generamos la órbita iterando la transformación
            while brk2 == 0 and len(orbita) < n:
                last_val = aux[-1]
                num = Binv[0,0]*last_val + Binv[0,1]
                deno = Binv[1,0]*last_val + Binv[1,1]
                
                if deno != 0:
                    new_val = num/deno
                    aux.append(new_val)
                    
                    # Mapeamos el valor numérico de vuelta a un índice de lugar
                    indice = -1
                    for k_idx, m_val in enumerate(menoselem):
                        if m_val == new_val: indice = k_idx; break
                    
                    if indice != -1: orbita.append(places_finite[indice])
                    else: orbita = []; brk2 = 1 # Fallo en búsqueda
                else:
                    orbita = []; brk2 = 1 # Polo encontrado (infinito)
            
            if len(orbita) == n:
                orbitas.append(orbita)

    # --- CONSTRUCCIÓN DE CÓDIGOS ---
    matrices_halladas = []
    for Q in lg2fijos:
        # Construimos el espacio de funciones L(G) usando Riemann-Roch.
        # El divisor G está asociado al lugar fijo Q.
        I = O.ideal(Q)
        Base = I.divisor_of_zeros().basis_function_space()
        
        for orbita in orbitas:
            # Evaluamos la base en los puntos de la órbita para obtener la matriz generadora.
            lista = []
            for base in Base:
                for p in orbita: lista.append(base.evaluate(p))
            try:
                # Matriz 'cruda': tal cual sale de la evaluación.
                M_cruda = matrix(k, len(Base), len(orbita), lista)
                C_temp = LinearCode(M_cruda)
                
                # Filtramos solo códigos de dimensión 3.
                if C_temp.dimension() == 3: 
                    # Convertimos a Forma Sistemática: [ I | A ]
                    # Esto es crucial para estandarizar las matrices antes de devolverlas.
                    M_sys = C_temp.systematic_generator_matrix()
                    # Devolvemos los datos necesarios para reconstruir la matriz.
                    matrices_halladas.append((M_sys.nrows(), M_sys.ncols(), M_sys.list()))
            except: pass
            
    return matrices_halladas

# ==============================================================================
# 2. Optimización para el calculo monomial
# ==============================================================================
def obtener_reduccion(G):
    """
    reduce una matriz generadora sistemática para facilitar la comparación monomial.
    
    Método:
        Dada M = [I | A], recorre las columnas de A (parte no identidad).
        Para cada columna, busca el primer elemento no nulo (pivote).
        Multiplica toda la columna por el inverso de ese pivote para que se vuelva 1.
    
    Esto reduce el problema de Equivalencia Monomial a un problema de Permutación.
    """
    M = copy(G)                        # Trabajamos sobre una copia
    n = M.ncols()                      # hallamos la cantidad de columnas
    k = M.nrows()                      # hallamos la cantidad de filas
    
    for j in range(k, n):              # Iteramos desde la columna k (saltando la identidad)
        for i in range(0, k):          # Recorremos filas
            if M[i, j] != 0:           # Encontramos pivote
                inv = M[i, j].inverse()# Calculamos la inversa de dicho elemento.
                M.rescale_col(j, inv)  # Multiplicamos la columna por el inverso del pivote.
                break                  # Pasamos a la siguiente columna
    return M                           # Devolvemos la matriz ya reducida.

# ==============================================================================
# 3. MANAGER PRINCIPAL (CONTROLADOR)
# ==============================================================================
def run_dynamic_program(q):       
    t_global_start = time.time()
    print(90*'=')
    print('INICIO DEL PROGRAMA GENÉRICO (PARALELO + DINÁMICO)')
    print(f'>>> INICIO DEL CASO Fq con q = {q}')
    k = GF(q, 'a'); a = k.gen()         # Inicialización del cuerpo de tamaño q y generador a
    print(f'Polinomio característico: {a.charpoly("y")}')
    
    # --- PREPARACIÓN DE DATOS ---
    
    elem = []                            # Generamos listas de elementos para búsquedas,
    if is_prime(q): elem = [i for i in k]# si q primo, es tan solo los elementos.
    else:
        for i in range(q):               # Si q no es primo, guardamos los elementos
            if i == 0: elem.append(a*0)  # como potencias del generador.
            else: elem.append(a**i)      
    menoselem = [-e for e in elem]       # Generamos la lista de menos los elementos, usada más adelante para el isomorfismo sigma.
    
    ### Esto lo hizo la ia, dios sabrá si está bien, pero los resultados son suficientemente concluyentes para asegurar que lo está.
    # Serialización: Convertimos a listas de coeficientes simples
    # Esto evita errores de 'Pickling' al pasar objetos complejos a procesos paralelos.
    elem_coeffs = [e.polynomial().list() for e in elem]          
    menoselem_coeffs = [e.polynomial().list() for e in menoselem]
    
    # Generación de Polinomios Irreducibles Mónicos de Grado 2
    R = PolynomialRing(k, 'x')    #Construimos el anillo de polinomios asociado al cuerpo, con variable 'x'.  
    monicirred2 = [p for p in R.polynomials(2) if p.is_irreducible() and p.is_monic()] #conjunto de los polinomios mónicos irreducibles de grado 2.
    monicirred2_coeffs = [p.list() for p in monicirred2]                               #guardamos los coeficientes, por lo mismo que arriba.
    print(f'Hay {len(monicirred2)} polinomios mónicos irreducibles de grado 2')  
    ###
    
    # --- CLASIFICACIÓN DE MATRICES ---
    # Buscamos las longitudes 'n' válidas para códigos sigma-cíclicos.
    print("Clasificando matrices PGL por orden...")
    matrices_por_orden = {} 
    
    for pmi in monicirred2:                               # Para cada polinomios mónicos irreducibles (pmi), 
        coeffs = pmi.list(); a0, a1 = coeffs[0], coeffs[1]# llamamos a0 y a1 al termino independiente y el lineal respetivamente.
        for t in k: #para cada elemento del cuerpo,
            B = matrix(k, 2, [t, -a0, 1, a1+t])  # Construimos matriz candidata asociada al polinomio
            if B.det() != 0:                    
                # Calculamos orden en PGL (B^n = escalar * I)
                n_gl = B.multiplicative_order() 
                for d in divisors(n_gl):
                    B_d = B**d
                    if B_d == B_d[0,0] * identity_matrix(k, 2):
                        n_gl = d; break
                
                # Filtramos según criterios: n>3, n<q+1, n|q+1. El abordado en la tesis.
                if n_gl > 3 and n_gl < q+1 and (q+1) % n_gl == 0:
                    if n_gl not in matrices_por_orden: matrices_por_orden[n_gl] = []
                    matrices_por_orden[n_gl].append(B)
    
    longitudes_validas = sorted(matrices_por_orden.keys())
    print(f'Las longitudes posibles (n) son: {longitudes_validas}')

    # --- BUCLE PRINCIPAL POR LONGITUD ---
    for n in longitudes_validas:
        t_n_start = time.time()
        matrices_candidatas = matrices_por_orden[n]
        print(30*'=')
        print(f'ESTUDIAREMOS LONGITUD n = {n}')
        print(f'Hay {len(matrices_candidatas)} matrices candidatas de orden {n}')
        
        ##################################
        #Aquí, hacemos la llamda a la funcion que busca los códgios sigma ciclicos.
        #Lo que se obtiene son todas las matrices que generan códigos AG-SC, en forma estandar.
        # Preparamos tareas y lanzamos workers
        inputs = [(B.list(), n, q, monicirred2_coeffs, elem_coeffs, menoselem_coeffs) 
                  for B in matrices_candidatas]

        MATGEN_RAW = []
        # Recolectamos resultados de todos los núcleos
        for _, res in worker_analizar_matriz(inputs):
            if res: MATGEN_RAW.extend(res)
        
        ###################################
        
        # 2. REDUCCIÓN DE DUPLICADOS
        #Una ves que hemos encontrado todos las matrices, 
        #eliminamos matrices idénticas generadas por órbitas equivalentes.
        Theta_set = set()
        Theta = []
        for (rows, cols, entries) in MATGEN_RAW:
            M_clean = matrix(k, rows, cols, entries)
            M_clean.set_immutable()     # Necesario para usar en set.
            if M_clean not in Theta_set:# Si la matriz no está en el conjunto, la metemos y tambien la metemos
                Theta_set.add(M_clean)  # a la lista theta. si está en el conjunto, ya está en la lista,
                Theta.append(M_clean)   # por lo tanto, no guardamos esa matriz.
                                      
        print(30*'-')
        print(f'HAY {len(Theta)} MATRICES GENERADORAS DISTINTAS (n={n})')
        if len(Theta) == 0: continue

        # 3. FILTRADO POR PERMUTACIÓN
        # Buscamos clases de equivalencia bajo permutación de columnas.
        print("Filtrando por permutaciones...")
        t_perm_start = time.time()
        codigos_perm = []
        reps_perm = [] # Representantes de cada clase
        
        for mu in Theta:
            C1 = LinearCode(mu)
            es_nuevo = True
            for C2 in codigos_perm:
                if C1.is_permutation_equivalent(C2):
                    es_nuevo = False; break
            if es_nuevo:
                codigos_perm.append(C1)
                reps_perm.append(mu)
        
        print(f'Existen {len(reps_perm)} códigos NO equivalentes por permutación.')
        for r in reps_perm: print(r); print('')
        print(f"Tiempo Permutación: {time.time() - t_perm_start:.2f} s")
        
        # 4. FILTRADO MONOMIAL (VIA NORMALIZACIÓN)
        # Verificamos si las clases colapsan al permitir escalado de columnas.
        print('Filtrando por equivalencia monomial...')
        t_mono_start = time.time()
        
        lista_reducida = [] # Paso 1: Tomar cada matriz que pasó el filtro de permutaciones
        for G in reps_perm:    # y reducirlas a traves de la funcion obtener_reducción()
            G_norm = obtener_reduccion(G)   #Llamamos a la matriz reducida G_norm y la metemos en la lista 
            lista_reducida.append(G_norm)
            
        # Paso B: Comparar las versiones reducidas
        codigos_finales = []
        matrices_finales = []
        
        for G_norm_cand in lista_reducida:
            C_norm_cand = LinearCode(G_norm_cand)
            
            es_nuevo = True
            for C_aceptado in codigos_finales:
                # Al estar reducidas, Monomial se reduce a Permutación
                if C_norm_cand.is_permutation_equivalent(C_aceptado):
                    es_nuevo = False
                    break
            
            if es_nuevo:
                codigos_finales.append(C_norm_cand)
                matrices_finales.append(G_norm_cand) 

        # --- REPORTE FINAL ---
        print(30*'*')
        print(f'RESULTADO DEFINITIVO (n={n}): Existen {len(matrices_finales)} códigos NO equivalentes (Monomial).')
        print('Las matrices generadoras finales son:')
    
        for M in matrices_finales:
            print(M)
            print("-" * 20)

        print(f"Tiempo Monomial: {time.time() - t_mono_start:.2f} s")
        print(f"Tiempo TOTAL n={n}: {time.time() - t_n_start:.2f} s")

    print(90*'=')
    print(f"Tiempo Global: {time.time() - t_global_start:.2f} s")
    print('FIN DEL PROGRAMA')

# Ejecución principal sobre la lista de muestras
for x in muestra: 
    run_dynamic_program(x)
