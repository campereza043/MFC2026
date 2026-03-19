# Método de Crank-Nicolson

Este directorio contiene implementaciones del método de Crank-Nicolson para resolver la **Ecuación de Difusión (o de Calor) en 1D**:

$$\frac{\partial u}{\partial t} = \alpha^2 \frac{\partial^2 u}{\partial x^2}$$

## Descripción del Método

El método de Crank-Nicolson es un esquema de diferencias finitas de segundo orden tanto en el espacio como en el tiempo. Es un método implícito que resulta ser **incondicionalmente estable** para la ecuación de calor.

La discretización promedia el operador de diferencia espacial en el paso de tiempo actual $n$ y el siguiente $n+1$:

$$\frac{u_i^{n+1} - u_i^n}{\Delta t} = \frac{\alpha^2}{2} \left[ \frac{u_{i+1}^{n+1} - 2u_i^{n+1} + u_{i-1}^{n+1}}{\Delta x^2} + \frac{u_{i+1}^{n} - 2u_i^{n} + u_{i-1}^{n}}{\Delta x^2} \right]$$

Esto conduce a un sistema tridiagonal que se resuelve eficientemente en cada paso de tiempo.

## Contenidos

1.  **`crank_nicolson.py`**: Implementación en Python utilizando `numpy` y `scipy.linalg.solve_banded`.
2.  **`crank_nicolson.cpp`**: Implementación en C++ utilizando el **Algoritmo de Thomas** para resolver el sistema tridiagonal de manera eficiente ($O(N)$).

## Cómo ejecutar

### Python
Asegúrate de tener instaladas las librerías `numpy`, `scipy` y `pandas`.
```bash
python3 crank_nicolson.py
```

### C++
Compila con un compilador estándar de C++ (como `g++`):
```bash
g++ -O2 -o crank_nicolson crank_nicolson.cpp
./crank_nicolson
```

## Resultados
Ambos códigos generan archivos de datos con la evolución temporal de la solución:
- `results_python.csv`
- `results_cpp.dat`

Además, comparan el resultado final con la solución analítica para validar la precisión del método.
