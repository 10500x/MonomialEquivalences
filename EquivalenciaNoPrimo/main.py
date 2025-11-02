import numpy as np
import sympy as sp
import galois
from sec import *

def main():
    q=[8,9,16,25,27]
    for x in q:
        GF = Crear_cuerpo(x)           #Construimos el cuerpo finito (1,2,3,4,...q)
        PMI2=generar_PMI2(x)   #Buscamos los polinomios irreducibles y decimos cuantos hay.
        IMS = generar_IMS(GF)      #Buscamos las IMS dado el cuerpo
        matord=orden_matriz(q)         #Buscamos los ordenes de cada matriz de IMS
        sigma=long_sigma(q)
        
        print(f"Hay {len(PMI2)} pol irred grado 2 en GF({x})")
        print(f"El polinomio caracteristico es: ({GF.irreducible_poly})")
        print(f"Matrices en PGL(2,{x}) sin fijar lugares: {len(IMS)}")
        print('Las longitudes posibles para codigos sigma ciclicos son',sigma)
    
if __name__ == "__main__":
    main()
