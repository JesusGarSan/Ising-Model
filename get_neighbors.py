# Función para calcular los vecinos de cada uno de los nodos de la red. SÓLO SE HACE UNA VEZ
def get_neighbors(L):
    n1 = [0] * (L * L)
    n2 = [0] * (L * L)
    n3 = [0] * (L * L)
    n4 = [0] * (L * L)
    
    for iy in range(L):
        for ix in range(L):
            i = iy * L + ix
            ix1 = (ix + 1) % L
            n1[i] = iy * L + ix1
            iy2 = (iy + 1) % L
            n2[i] = iy2 * L + ix
            ix3 = (ix - 1) % L
            n3[i] = iy * L + ix3
            iy4 = (iy - 1) % L
            n4[i] = iy4 * L + ix
    
    return n1, n2, n3, n4