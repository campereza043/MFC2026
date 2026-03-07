import numpy as np
import matplotlib.pyplot as plt

# 1. Cargar datos del Problema 2
try:
    data = np.loadtxt('resultados_p2_m1000.dat')
except FileNotFoundError:
    print("Error: No se encontró 'resultados_p2_m1000.dat'.")
    exit()

# 2. Identificar el último instante de tiempo guardado
t_final = np.max(data[:, 0])
print(f"Graficando el último tiempo disponible: t = {t_final:.4f}")

# 3. Filtrar datos para ese tiempo específico
mask = data[:, 0] == t_final
x = data[mask, 1]
u_real = data[mask, 2]
u_imag = data[mask, 3]

# Solución analítica para t_final: u = exp(it) * cos(x)
u_exacta = np.exp(1j * t_final) * np.cos(x)

# 4. Crear la figura
plt.figure(figsize=(12, 6))

# Subtrama 1: Parte Real
plt.subplot(1, 2, 1)
plt.plot(x, u_exacta.real, '-', color='mediumblue', label='Analítica (Real)', alpha=0.3)
plt.plot(x, u_real, '*', color='deeppink', markersize=8, label=f'C++ (t={t_final:.2f})')
plt.title('Parte Real de la Función de Onda')
plt.xlabel('x')
plt.ylabel('Re(ψ)')
plt.legend()
plt.grid(True, linestyle=':')

# Subtrama 2: Parte Imaginaria
plt.subplot(1, 2, 2)
plt.plot(x, u_exacta.imag, '-', color='mediumpurple', label='Analítica (Imag)', alpha=0.3)
plt.plot(x, u_imag, '*', color='purple', markersize=8, label=f'C++ (t={t_final:.2f})')
plt.title('Parte Imaginaria de la Función de Onda')
plt.xlabel('x')
plt.ylabel('Im(ψ)')
plt.legend()
plt.grid(True, linestyle=':')



plt.tight_layout()

# 5. Guardado automático
nombre_salida = "Graficasm1000_P2.png"
plt.savefig(nombre_salida, dpi=300)
print(f"Gráfica guardada como: {nombre_salida}")

plt.show()