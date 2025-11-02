from sec import *

def main():
    q=[8,9,16,25,27]
    for x in q:
        GF = Crear_cuerpo(x)           #Construimos el cuerpo finito (1,2,3,4,...q)
        PMI2=generar_PMI2(x)   #Buscamos los polinomios irreducibles y decimos cuantos hay.
        GL2=generar_GL2(GF)    # Construimos las matrices invertibles.
        PGL=generar_PGL2(GL2,GF)# Subconjunto de GL2 que no fija lugares racionales
        IMS = generar_IMS(GF,GL2)      #Buscamos las IMS dado el cuerpo
        matord = orden_matriz(IMS, x)#Buscamos los ordenes de cada matriz de IMS    
        sigma=long_sigma(x,matord)
        print(f"Hay {len(PMI2)} pol irred grado 2 en GF({x})")
        print(f"El polinomio caracteristico es: ({GF.irreducible_poly})")
        print(f"Matrices en PGL(2,{x}) sin fijar lugares: {len(IMS)} racionales")
        print('Las longitudes posibles para codigos sigma ciclicos son',sigma)
        print(f"fin caso {x}")
        print("")

if __name__ == "__main__":
    main()
