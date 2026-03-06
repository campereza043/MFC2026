import numpy as np
import matplotlib.pyplot as plt

# Cargar los datos del archivo .dat
# Skiprows=1 para saltar la cabecera con '#'
try:
    data = np.loadtxt('resultados.dat')
    x = data[:, 0]
    u_real = data[:, 1]
    u_imag = data[:, 2]
    u_abs = data[:, 3]
except FileNotFoundError:
    print("Error: No se encontró el archivo 'resultados.dat'. Ejecuta primero el código C++.")
    exit()

# Solución analítica para validar (t = 1.0)
u_exacta = np.exp(1j * 1.0) * np.sin(x)

# Crear la figura
plt.figure(figsize=(12, 6))

# Subtrama 1: Parte Real
plt.subplot(1, 2, 1)
plt.plot(x, u_exacta.real, 'k-', label='Analítica (Real)', alpha=0.5)
plt.plot(x, u_real, 'r*', label='C++ (Crank-Nicolson)')
plt.title('Parte Real de la Función de Onda')
plt.xlabel('Posición x')
plt.ylabel('Re(ψ)')
plt.legend()
plt.grid(True)

# Subtrama 2: Parte Imaginaria
plt.subplot(1, 2, 2)
plt.plot(x, u_exacta.imag, 'k-', label='Analítica (Imag)', alpha=0.5)
plt.plot(x, u_imag, 'b*', label='C++ (Crank-Nicolson)')
plt.title('Parte Imaginaria de la Función de Onda')
plt.xlabel('Posición x')
plt.ylabel('Im(ψ)')
plt.legend()
plt.grid(True)



plt.tight_layout()
plt.show()