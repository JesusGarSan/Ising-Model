# # MODELO DE ISING 2D - Algoritmo de Metrópolis
# importación de librerías
import numpy as np
import csv
from datetime import datetime
from get_neighbors import get_neighbors
from tqdm import tqdm

global L, N
L = 10 # Longitud de un lado del retículo cuadrado
# for L in [10, 20, 30, 40, 60, 80]:
for L in [30, 40]:
    print(f"""
             #############################################
             ########### EJECUTANDO PARA L = {L} #########
             #############################################  
          """)
    N = L*L # Número de nodos
        
    n1, n2, n3, n4 = get_neighbors(L)

    # Condición inicial aleatoria
    np.random.seed(20082022)
    s = np.random.choice([1,-1], N)

    T_cutoffs = [4.0, 2.4, 2.32, 2.25, 2.2, 0.0]
    T_steps =   [0.1, 0.01, 0.001, 0.01, 0.1]

    T_range = np.concatenate((
        np.arange(T_cutoffs[1], T_cutoffs[2], -T_steps[1]),
        np.arange(T_cutoffs[2], T_cutoffs[3], -T_steps[2])[1:],
        np.arange(T_cutoffs[3], T_cutoffs[4], -T_steps[3])[1:],
        ))



    # print(f"Número de simulaciones a realizar para L = {L}: {len(T_range)}")
    pbar = tqdm(T_range)
    for T in pbar:
        T = round(T, 3)
        mc = 1


        if   T > T_cutoffs[1]:
            M0 = 1000 # Número de pasos Monte Carlo hasta termalizar
            M  = 5000 # Número de medidas a realizar

        elif T > T_cutoffs[2]:
            M0 = 2000 # Número de pasos Monte Carlo hasta termalizar
            M  = 8000 # Número de medidas a realizar

        elif T > T_cutoffs[3]:
            M0 = 3000 # Número de pasos Monte Carlo hasta termalizar
            M  = 10000 # Número de medidas a realizar

        elif T > T_cutoffs[4]:
            M0 = 2000 # Número de pasos Monte Carlo hasta termalizar
            M  = 8000 # Número de medidas a realizar

        elif T > T_cutoffs[5]:
            M0 = 100 # Número de pasos Monte Carlo hasta termalizar
            M  = 100 # Número de medidas a realizar


        # Damos más tiempo para termalizar si partimos de estado inicial aleatorio
        if T == max(T_range): M0 *= 2

        pbar.set_description(f'Simulando para: L={L}, T={T}, M={M}, M0={M0}, mc={mc}')
        
        
        # print(f"###### Simulando para: ######")
        # print(f"L={L}, T={T}, M={M}, M0={M0}, mc={mc}")
        # print(f"Número de pasos MC a realizar: {M*mc}")

        # Cálculo del parámetro de aceptación
        h = np.zeros(5)
        for j in np.linspace(-4, 4, 5, dtype=int):
            h[j] = np.min([1.0, np.exp(-2*j/T)])

        # Termalización
        current_time = datetime.now()
        # print(f"## Termalizamos con {M0} pasos MC ({M0*N} iteraciones) ##")
        pbar.set_description(f'Termalizando para: L={L}, T={T}, M={M}, M0={M0}, mc={mc}')
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
        
        previous_time = current_time
        current_time = datetime.now()
        elapsed_time = current_time-previous_time
        # print(f"Tiempo transcurrido: {elapsed_time}\n\n")

        # ### Evolución del sistema
        # print(f"## Comienzo de la toma de medidas ##")
        pbar.set_description(f'Midiendo para: L={L}, T={T}, M={M}, M0={M0}, mc={mc}')
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
            previous_time = current_time
            current_time = datetime.now()
            elapsed_time = current_time-previous_time
            # print(f"Tiempo transcurrido: {elapsed_time}\n\n")
            writer.writerows([[T, L, M, M0, mc, rm, rm2, rm4, error, mc*tau, c, current_time]])
            

        # with open(f'./Resultados/rm0_{L}_{T}_{M0}_{M}_{mc}.csv', 'w') as file:
        #     writer = csv.writer(file, delimiter=',')
        #     writer.writerows([RM0])
            
        # with open(f'./Resultados/sfinal_{L}_{T}_{M0}_{M}_{mc}.csv', 'w') as file:
        #     writer = csv.writer(file, delimiter=',')
        #     writer.writerows([s])


