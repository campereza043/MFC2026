import numpy as np
import matplotlib.pyplot as plt

# 1. Cargar los datos generados por C++
try:
    data = np.loadtxt('resultados_cpp.dat')
except FileNotFoundError:
    print("Error: No se encontró 'resultados_cpp.dat'. Ejecuta el código C++ primero.")
    exit()

# 2. Extraer tiempos únicos
tiempos_totales = np.unique(data[:, 0])

# 3. Configuración de la gráfica
plt.figure(figsize=(12, 7))

# Seleccionamos una cantidad moderada de curvas para que los puntos no se encimen
num_curvas = 6 
indices = np.linspace(0, len(tiempos_totales) - 1, num_curvas, dtype=int)
tiempos_a_graficar = tiempos_totales[indices]

# Lista de colores estándar para evitar el mapa "plasma"
colores_estandar = ['blue', 'red', 'green', 'magenta', 'orange', 'black']

print(f"Graficando puntos para {len(tiempos_a_graficar)} instantes de tiempo...")

for i, t in enumerate(tiempos_a_graficar):
    # Filtrar datos para el tiempo t
    mask = data[:, 0] == t
    x = data[mask, 1]
    u_real = data[mask, 2]
    
    # Graficamos con puntos ('o') y una línea muy tenue de fondo
    plt.plot(x, u_real, '*', color=colores_estandar[i % len(colores_estandar)], 
             markersize=4, alpha=0.7, label=f't = {t:.4f}')
    
    # Opcional: una línea punteada muy suave para guiar la vista
    plt.plot(x, u_real, color=colores_estandar[i % len(colores_estandar)], 
             linestyle='--', linewidth=0.5, alpha=0.3)

# 4. Estética
plt.title('Evolución de la Parte Real (Puntos de la Malla C++)', fontsize=14)
plt.xlabel('Posición (x)', fontsize=12)
plt.ylabel('Re(ψ)', fontsize=12)
plt.axhline(0, color='black', linewidth=0.8, linestyle='-')
plt.grid(True, linestyle=':', alpha=0.5)



plt.legend(title="Instantes de Tiempo", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# Guardar y mostrar
plt.savefig('TiemposR.png', dpi=300)
plt.show()