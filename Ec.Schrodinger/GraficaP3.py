import numpy as np
import matplotlib.pyplot as plt

try:
    data = np.loadtxt('resultados_p3_m50.dat')
except FileNotFoundError:
    print("Error: No se encontró el archivo .dat")
    exit()

t_final = np.max(data[:, 0])
mask = data[:, 0] == t_final
x = data[mask, 1]
u_real = data[mask, 2]
u_imag = data[mask, 3]

# Solución analítica corregida: exp(i*t + x)
u_exacta = np.exp(1j * t_final + x)

plt.figure(figsize=(12, 5))

# Parte Real
plt.subplot(1, 2, 1)
plt.plot(x, u_exacta.real, 'k-', alpha=0.3, label='Analítica')
plt.plot(x, u_real, 'r*', label=f'C++ (t={t_final:.2f})')
plt.title('Parte Real: $e^x \cos(t)$')
plt.grid(True, ls=':')
plt.legend()

# Parte Imaginaria
plt.subplot(1, 2, 2)
plt.plot(x, u_exacta.imag, 'k-', alpha=0.3, label='Analítica')
plt.plot(x, u_imag, 'b*', label=f'C++ (t={t_final:.2f})')
plt.title('Parte Imaginaria: $e^x \sin(t)$')
plt.grid(True, ls=':')
plt.legend()



plt.tight_layout()
plt.savefig('grafica_final_problema3.png')
plt.show()