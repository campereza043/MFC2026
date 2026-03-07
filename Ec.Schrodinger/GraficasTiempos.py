import numpy as np
import matplotlib.pyplot as plt

# 1. Cargar datos del Problema 2
try:
    data = np.loadtxt('resultados_p1_m20.dat')
except FileNotFoundError:
    print("Error: No se encontró 'resultados_p1_m20.dat'.")
    exit()

# 2. Definir los tiempos que quieres graficar
tiempos_deseados = [1.0, 0.75, 0.5, 0.25]
tiempos_reales_en_archivo = np.unique(data[:, 0])

# Encontrar los tiempos más cercanos disponibles en el archivo .dat
tiempos_a_graficar = []
for td in tiempos_deseados:
    idx = (np.abs(tiempos_reales_en_archivo - td)).argmin()
    tiempos_a_graficar.append(tiempos_reales_en_archivo[idx])

# 3. Configuración de la figura
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
colores = ['darkred', 'deeppink', 'cyan', 'purple']

for i, t_val in enumerate(tiempos_a_graficar):
    # Filtrar datos por tiempo
    mask = data[:, 0] == t_val
    x = data[mask, 1]
    u_real = data[mask, 2]
    u_imag = data[mask, 3]
    
    # Solución analítica para este tiempo exacto
    u_exacta = np.exp(1j * t_val) * np.sin(x)
    
    # --- Subtrama 1: Parte Real ---
    ax1.plot(x, u_exacta.real, '-', color='black', alpha=0.2)
    ax1.plot(x, u_real, '*', color=colores[i], markersize=7, label=f't ≈ {t_val:.2f}')
    
    # --- Subtrama 2: Parte Imaginaria ---
    ax2.plot(x, u_exacta.imag, '-', color='black', alpha=0.2)
    ax2.plot(x, u_imag, '*', color=colores[i], markersize=7, label=f't ≈ {t_val:.2f}')

# Estética Subtrama 1 (Real)
ax1.set_title('Parte Real (Re(ψ))')
ax1.set_xlabel('Posición x')
ax1.set_ylabel('Amplitud')
ax1.legend()
ax1.grid(True, linestyle=':', alpha=0.6)

# Estética Subtrama 2 (Imaginaria)
ax2.set_title('Parte Imaginaria (Im(ψ))')
ax2.set_xlabel('Posición x')
ax2.set_ylabel('Amplitud')
ax2.legend()
ax2.grid(True, linestyle=':', alpha=0.6)



plt.tight_layout()

# 4. Guardado automático
nombre_img = "Tiempos_P1.png"
plt.savefig(nombre_img, dpi=300)
print(f"Gráfica guardada como {nombre_img}")
plt.show()