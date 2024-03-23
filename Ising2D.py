# # MODELO DE ISING 2D - Algoritmo de Metrópolis
# importación de librerías
import numpy as np
import csv
from datetime import datetime
from get_neighbors import get_neighbors

global L, N
L = 100 # Longitud de un lado del retículo cuadrado
N = L*L # Número de nodos
    
n1, n2, n3, n4 = get_neighbors(L)

# Condición inicial aleatoria
s = np.random.choice([1,-1], N)

T_step = 0.1 # Paso entre cada temperatura para la que se va a ejecutar
T_fine_step = 0.01 # Paso más fino para temperaturas cercanas a la crítica

T_range = np.concatenate((
    np.arange(4.0, 2.4, -T_step),
    np.arange(2.4, 2.0, -T_fine_step),
    np.arange(2.0, 0.0, -T_step)
))

for T in  T_range:
    # T = 2 # Temperatura
    T = round(T, 3)
    M = 500 # Número de medidas a realizar
    M0 = 100 # Número de pasos Monte Carlo hasta termalizar
    mc = 5 # Número de pasos Monte Carlo entre cada medida
    
    if T > 2.0 and T < 2.4: M0 = 500 # Más pasos de termalización cerca de T_c
    
    
    print(f"### Simulando para: ###")
    print(f"L={L}, T={T}, M={M}, M0={M0}, mc={mc}")
    print(f"Número de pasos Monte Carlo a realizar: {M*mc}")

    # Cálculo del parámetro de aceptación
    h = np.zeros(5)
    for j in np.linspace(-4, 4, 5, dtype=int):
        h[j] = np.min([1.0, np.exp(-2*j/T)])

    # Termalización
    print(f"## Termalizamos con {M0} pasos MC ({M0*N} iteraciones) ##")
    for _ in range(M0*N):
        i = np.random.randint(N) # Escogemos un índice al azar
        ib = s[i] *( s[n1[i]] + s[n2[i]] + s[n3[i]] + s[n4[i]] )
        if np.random.rand() < h[ib]: s[i] = -s[i]

    #### Inicialización de los promedios
    RM0 = np.zeros(M) # magnetización a lo largo de la evolución del sistema
    c = 0. # Correlación media
    rm = 0. # Magnetización media
    rm2 = 0. # magnetización media ^ 2
    rm4 = 0. # magnetización media ^ 4
    rm1 = np.abs(np.sum(s))/N

    # ### Evolución del sistema
    print(f"## Comienzo de la toma de medidas ##\n\n")
    for j in range(M): # Número de medidas a realizar
        # print(f"Medida {j}")
        for _ in range(mc*N): # Número de pasos Monte Carlo entre medidas
            i = np.random.randint(N) # Escogemos un índice al azar
            ib = s[i] *( s[n1[i]] + s[n2[i]] + s[n3[i]] + s[n4[i]] )
            if np.random.rand() < h[ib]: s[i] = -s[i]
            
        rm0 = np.abs(np.sum(s))/N # Magnetización actual
        RM0[j] = rm0
        rm = rm + rm0
        rm2 = rm2 + rm0*rm0
        rm4 = rm4 + rm0*rm0*rm0*rm0
        c = c + rm0*rm1
        rm1 = rm0

    # Cálculo y guardado de medidas
    rm = rm/M
    rm2 = rm2/M
    rm4 = rm4/M
    c = (c/M - rm*rm)/rm2
    tau = 0
    if c != 1.0: tau = c/(1.-c)
    error = np.sqrt((rm2 - rm*rm)*(2*tau+1)/M)

    # Guardado de resultados
    with open('./Resultados/medidas.csv', 'a', newline='') as file:
        writer = csv.writer(file, delimiter=',')
        writer.writerows([[T, L, M, M0, mc, rm, rm2, rm4, error, mc*tau, c, datetime.now()]])

    # with open(f'./Resultados/rm0_{L}_{T}_{M0}_{M}_{mc}.csv', 'w') as file:
    #     writer = csv.writer(file, delimiter=',')
    #     writer.writerows([RM0])
        
    # with open(f'./Resultados/sfinal_{L}_{T}_{M0}_{M}_{mc}.csv', 'w') as file:
    #     writer = csv.writer(file, delimiter=',')
    #     writer.writerows([s])


