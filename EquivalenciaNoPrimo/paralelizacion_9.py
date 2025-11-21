import sage.all
from sage.parallel.decorate import parallel
import time
from itertools import product

muestra = [9, 16, 25]

# ==============================================================================
# 1. EL WORKER DE MINERÍA (Tu código original)
# ==============================================================================
@parallel
def worker_analizar_matriz(B_list, n, q, monicirred2_coeffs, elem_coeffs, menoselem_coeffs):
    # Reconstrucción
    k = GF(q, 'a'); a = k.gen()
    R = PolynomialRing(k, 'x')
    K = FunctionField(k, 'x'); x = K.gen()
    O = K.maximal_order()
    
    places_finite = K.places_finite()
    
    elem = [k(c) for c in elem_coeffs]
    menoselem = [k(c) for c in menoselem_coeffs]
    monicirred2 = [R(coeffs) for coeffs in monicirred2_coeffs]
    
    B = matrix(k, 2, B_list)
    Binv = matrix(k, 2, [B[1,1], -B[0,1], -B[1,0], B[0,0]])

    try:
        sigma = K.hom([(B[0,0]*x + B[0,1]) / (B[1,0]*x + B[1,1])])
    except: return []

    lg2fijos = []
    for pol in monicirred2:
        im = sigma(pol)
        num = R(im.numerator())
        den = R(im.denominator())
        if num % pol == 0 and den % pol != 0:
            lg2fijos.append(pol)
    
    if not lg2fijos: return []

    orbitas = []
    for i in range(q):
        val_buscado = menoselem[i]
        j = -1
        for idx, val in enumerate(elem):
            if val == val_buscado:
                j = idx; break
        
        P_start = places_finite[j]
        ya_esta = False
        for orb in orbitas:
            if P_start in orb: ya_esta = True; break
            
        if not ya_esta:
            orbita = [P_start]
            if i == 0: val_aux = k(1) 
            else: val_aux = a**i
            aux = [val_aux]
            
            brk2 = 0
            while brk2 == 0 and len(orbita) < n:
                last_val = aux[-1]
                num = Binv[0,0]*last_val + Binv[0,1]
                deno = Binv[1,0]*last_val + Binv[1,1]
                if deno != 0:
                    new_val = num/deno
                    aux.append(new_val)
                    indice = -1
                    for k_idx, m_val in enumerate(menoselem):
                        if m_val == new_val: indice = k_idx; break
                    if indice != -1: orbita.append(places_finite[indice])
                    else: orbita = []; brk2 = 1
                else: orbita = []; brk2 = 1
            
            if len(orbita) == n:
                orbitas.append(orbita)

    matrices_found = []
    for Q in lg2fijos:
        I = O.ideal(Q)
        Base = I.divisor_of_zeros().basis_function_space()
        for orbita in orbitas:
            lista = []
            for base in Base:
                for p in orbita: lista.append(base.evaluate(p))
            try:
                M_raw = matrix(k, len(Base), len(orbita), lista)
                C_temp = LinearCode(M_raw)
                if C_temp.dimension() == 3: 
                    M_sys = C_temp.systematic_generator_matrix()
                    matrices_found.append((M_sys.nrows(), M_sys.ncols(), M_sys.list()))
            except: pass
            
    return matrices_found

# ==============================================================================
# 2. NUEVO: WORKER PARA EQUIVALENCIA MONOMIAL
# ==============================================================================
@parallel
def worker_check_monomial(mat_cand_list, mat_piv_list, rows, cols, q):
    """
    Verifica si mat_cand es equivalente a mat_pivote (Monomial) en un hilo separado.
    """
    k = GF(q, 'a')
    # Reconstruimos matrices
    M_cand = matrix(k, rows, cols, mat_cand_list)
    M_piv = matrix(k, rows, cols, mat_piv_list)
    
    C1 = LinearCode(M_cand)
    C2 = LinearCode(M_piv) # Pivote

    # 1. Permutación (Rápido)
    if C1.is_permutation_equivalent(C2): return True

    # 2. Monomial (Lento - Fuerza Bruta)
    n = cols
    G2 = C2.systematic_generator_matrix()
    unidades = [u for u in k if u != 0]
    
    # Iteramos diagonales
    for escalares in product(unidades, repeat=n-1):
        diag = [k(1)] + list(escalares)
        cols_scaled = [G2.column(j) * diag[j] for j in range(n)]
        G_scaled = matrix(k, cols_scaled).transpose()
        C2_scaled = LinearCode(G_scaled)
        if C1.is_permutation_equivalent(C2_scaled): return True
            
    return False

# ==============================================================================
# 3. NUEVO: FUNCIÓN GESTORA DE FILTRADO PARALELO
# ==============================================================================
def filtrar_monomial_paralelo(lista_matrices, q):
    """
    Recibe una lista de matrices (objetos Sage) y filtra las monomialmente equivalentes
    usando worker_check_monomial.
    """
    # Convertimos a formato serializable para pasar a workers: (lista, rows, cols)
    pendientes_data = [(M.list(), M.nrows(), M.ncols()) for M in lista_matrices]
    
    unicos_data = []
    
    while pendientes_data:
        # Tomamos el primero como pivote y lo guardamos
        pivote_data = pendientes_data.pop(0)
        unicos_data.append(pivote_data)
        
        if not pendientes_data: break
        
        # Preparamos inputs: comparar todos los pendientes contra este pivote
        # Args del worker: (cand_list, piv_list, rows, cols, q)
        r, c = pivote_data[1], pivote_data[2]
        inputs = [(cand[0], pivote_data[0], r, c, q) for cand in pendientes_data]
        
        # Ejecutamos en paralelo
        resultados_bool = [] # Guardará True/False en orden
        
        # IMPORTANTE: @parallel devuelve resultados desordenados.
        # Debemos asociar cada resultado a su candidato original.
        # Sage parallel devuelve: (((args), kwargs), result)
        
        # Mapeamos 'tuple(lista_matriz)' -> booleano
        mapa_res = {}
        for (args, _), es_equiv in worker_check_monomial(inputs):
            # args[0] es la lista de la matriz candidata
            mapa_res[tuple(args[0])] = es_equiv
            
        # Reconstruimos la lista de supervivientes (los que NO son equivalentes)
        nuevos_pendientes = []
        for p_data in pendientes_data:
            mat_list_tuple = tuple(p_data[0])
            # Si es equivalente (True), lo descartamos. Si es False, sobrevive.
            if not mapa_res.get(mat_list_tuple, False):
                nuevos_pendientes.append(p_data)
        
        pendientes_data = nuevos_pendientes

    # Reconstruimos objetos Matrix de Sage para devolver
    k = GF(q, 'a')
    matrices_finales = []
    for (l, r, c) in unicos_data:
        matrices_finales.append(matrix(k, r, c, l))
        
    return matrices_finales

# ==============================================================================
# 4. MANAGER DINÁMICO
# ==============================================================================
def run_dynamic_program(q):
    t_global_start = time.time()
    print(90*'=')
    print('INICIO DEL PROGRAMA GENÉRICO (PARALELO TOTAL)')
    
    print(f'>>> INICIO DEL CASO Fq con q = {q}')
    k = GF(q, 'a'); a = k.gen()
    print(f'Polinomio característico: {a.charpoly("y")}')
    
    # Setup listas
    elem = []
    if is_prime(q): elem = [i for i in k]
    else:
        for i in range(q):
            if i == 0: elem.append(a*0)
            else: elem.append(a**i)
    menoselem = [-e for e in elem]
    
    elem_coeffs = [e.polynomial().list() for e in elem]
    menoselem_coeffs = [e.polynomial().list() for e in menoselem]
    
    R = PolynomialRing(k, 'x')
    monicirred2 = [p for p in R.polynomials(2) if p.is_irreducible() and p.is_monic()]
    monicirred2_coeffs = [p.list() for p in monicirred2]
    print(f'Hay {len(monicirred2)} polinomios mónicos irreducibles de grado 2')

    # 1. BÚSQUEDA DE LONGITUDES
    print("Clasificando matrices PGL por orden...")
    matrices_por_orden = {} 
    
    for pmi in monicirred2:
        coeffs = pmi.list()
        a0, a1 = coeffs[0], coeffs[1]
        for t in k:
            B = matrix(k, 2, [t, -a0, 1, a1+t])
            if B.det() != 0:
                n_gl = B.multiplicative_order()
                n = n_gl
                for d in divisors(n_gl):
                    B_d = B**d
                    if B_d == B_d[0,0] * identity_matrix(k, 2):
                        n = d; break
                
                if n > 3 and n < q+1 and (q+1) % n == 0:
                    if n not in matrices_por_orden: matrices_por_orden[n] = []
                    matrices_por_orden[n].append(B)
    
    longitudes_validas = sorted(matrices_por_orden.keys())
    print(f'Las longitudes posibles (n) son: {longitudes_validas}')

    # 2. BUCLE POR LONGITUD
    for n in longitudes_validas:
        t_n_start = time.time()
        matrices_candidatas = matrices_por_orden[n]
        print(30*'=')
        print(f'ESTUDIAREMOS LONGITUD n = {n}')
        print(f'Hay {len(matrices_candidatas)} matrices candidatas de orden {n}')
        
        # MINERÍA PARALELA
        inputs = [(B.list(), n, q, monicirred2_coeffs, elem_coeffs, menoselem_coeffs) 
                  for B in matrices_candidatas]
        
        MATGEN_RAW = []
        for _, res in worker_analizar_matriz(inputs):
            if res: MATGEN_RAW.extend(res)
            
        Theta_set = set()
        Theta = []
        for (rows, cols, entries) in MATGEN_RAW:
            M_clean = matrix(k, rows, cols, entries)
            M_clean.set_immutable()
            if M_clean not in Theta_set:
                Theta_set.add(M_clean)
                Theta.append(M_clean)
        
        print(30*'-')
        print(f'HAY {len(Theta)} MATRICES GENERADORAS DISTINTAS (n={n})')
        
        if len(Theta) == 0: continue

        # FILTRADO PERMUTACIÓN (Secuencial - suele ser rápido, pero se podría paralelizar también)
        print("Filtrando por permutaciones...")
        t_perm_start = time.time()
        codigos_perm = []
        reps_perm = [] # Matrices representantes de permutación
        
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
        
        # FILTRADO MONOMIAL (PARALELO)
        print('Filtrando por equivalencia monomial (PARALELO)...')
        t_mono_start = time.time()
        
        # Llamamos a la nueva función paralela pasando las matrices representantes de la etapa anterior
        reps_finales_matrices = filtrar_monomial_paralelo(reps_perm, q)
                
        print(30*'*')
        print(f'RESULTADO FINAL (n={n}): Existen {len(reps_finales_matrices)} códigos NO equivalentes (Monomial).')
        for r in reps_finales_matrices: print(r); print('')
        
        print(f"Tiempo Monomial: {time.time() - t_mono_start:.2f} s")
        print(f"Tiempo TOTAL n={n}: {time.time() - t_n_start:.2f} s")

    print(90*'=')
    print(f"Tiempo TOTAL Global: {time.time() - t_global_start:.2f} s")
    print('FIN DEL PROGRAMA')

for x in muestra:
    run_dynamic_program(x)