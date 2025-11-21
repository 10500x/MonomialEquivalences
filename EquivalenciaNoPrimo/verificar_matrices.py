from sage.all import *
from itertools import product

# ==============================================================================
# FUNCIÓN DE AYUDA PARA EQUIVALENCIA MONOMIAL
# ==============================================================================
def check_monomial_equivalence(C1, C2):
    """
    Verifica si C1 es equivalente a C2 bajo permutación y escalado de columnas.
    """
    # 1. Chequeo rápido de permutación
    if C1.is_permutation_equivalent(C2):
        return True
    
    # 2. Chequeo exhaustivo de escalado diagonal
    k = C1.base_field()
    n = C1.length()
    G2 = C2.systematic_generator_matrix()
    
    # Unidades no nulas para escalar
    unidades = [u for u in k if u != 0]
    
    # Fijamos la primera columna a 1 para optimizar (factor q-1 menos)
    for escalares in product(unidades, repeat=n-1):
        # Vector diagonal [1, d2, ..., dn]
        diag = [k(1)] + list(escalares)
        
        # Escalar columnas de G2
        cols_scaled = [G2.column(j) * diag[j] for j in range(n)]
        G_scaled = matrix(k, cols_scaled).transpose()
        
        # Crear código temporal y comparar permutación
        C2_scaled = LinearCode(G_scaled)
        if C1.is_permutation_equivalent(C2_scaled):
            return True
            
    return False

# ==============================================================================
# PROGRAMA PRINCIPAL
# ==============================================================================
def verificar_equivalencias():
    print(80 * "=")
    print("VERIFICACIÓN MATEMÁTICA DE EQUIVALENCIA DE CÓDIGOS (CORREGIDO)")
    print(80 * "=")

    # 1. Definir el cuerpo GF(9)
    k = GF(9, 'a')
    a = k.gen()
    print(f"Cuerpo base: {k}")

    # --- MATRICES ORIGINALES (SET 1 - DEL CÓDIGO VIEJO) ---
    M_old_1 = matrix(k, 3, 5, [
        [1, 0, 0, 2,     2*a],
        [0, 1, 0, a+1,   2*a+2],
        [0, 0, 1, 2*a+1, 2*a+2]
    ])

    M_old_2 = matrix(k, 3, 5, [
        [1, 0, 0, 2*a, 1],
        [0, 1, 0, a,   0],
        [0, 0, 1, 1,   0]
    ])

    # --- MATRICES NUEVAS (SET 2 - DEL CÓDIGO PARALELO) ---
    M_new_1 = matrix(k, 3, 5, [
        [1, 0, 0, 1,   a],
        [0, 1, 0, 2*a, 2*a],
        [0, 0, 1, a,   1]
    ])

    M_new_2 = matrix(k, 3, 5, [
        [1, 1, 0, 0, 1],
        [0, 0, 1, 0, 2*a],
        [0, 0, 0, 1, a]
    ])

    # Crear objetos LinearCode
    C_old_1 = LinearCode(M_old_1)
    C_old_2 = LinearCode(M_old_2)
    C_new_1 = LinearCode(M_new_1)
    C_new_2 = LinearCode(M_new_2)

    # Listas para iterar
    viejos = [("Viejo_1", C_old_1), ("Viejo_2", C_old_2)]
    nuevos = [("Nuevo_1", C_new_1), ("Nuevo_2", C_new_2)]

    print("\nIniciando comparación cruzada (Monomial Equivalence)...\n")

    match_found_count = 0

    for nombre_viejo, codigo_viejo in viejos:
        encontrado = False
        for nombre_nuevo, codigo_nuevo in nuevos:
            # USAMOS LA FUNCIÓN MANUAL CORREGIDA
            if check_monomial_equivalence(codigo_viejo, codigo_nuevo):
                print(f"✅ ÉXITO: {nombre_viejo} es equivalente a {nombre_nuevo}")
                encontrado = True
                match_found_count += 1
        
        if not encontrado:
            print(f"❌ ALERTA: {nombre_viejo} no tiene pareja en el conjunto nuevo.")

    print("\n" + 80 * "-")
    if match_found_count >= 2:
        print("CONCLUSIÓN: Los resultados son MATEMÁTICAMENTE IDÉNTICOS.")
        print("Las diferencias visuales se deben solo a cambio de base, permutación y escalado.")
    else:
        print("CONCLUSIÓN: Hay discrepancias matemáticas.")
    print(80 * "=")

if __name__ == "__main__":
    verificar_equivalencias()