# Benchmark de la reimplementación (fase 7)

Fuera del build del paquete (`^bench$` en `.Rbuildignore`). Los scripts se
ejecutan desde la raíz del paquete, p. ej. `Rscript bench/bench-engines.R`, y
guardan sus resultados en `bench/results/`. Equipo: Windows 11, 16 núcleos,
R 4.5.0, CmdStan vía cmdstanr, stanli 0.25.0. 2026-11-01.

## 1. Coste del gradiente por motor (`bench-engines.R`, `engines.csv`)

Modelo `cor` con prior unit-vector, N = 1000, J ∈ {10, 15, 20, 30, 60, 100},
M ∈ {2, 3, 5, 8} (con (J − M)² ≥ J + M) y 0/10/30 % de valores perdidos MCAR.
Programas compilados con cmdstanr; NUTS sin adaptación, con paso diminuto y
profundidad máxima 5, de modo que cada iteración da 31 pasos leapfrog. El
tiempo de muestreo de CmdStan entre los pasos es el coste de un gradiente.

**Datos completos** (µs por gradiente):

| J | M = 2 Chol / Wood | M = 3 | M = 5 | M = 8 |
|---|---|---|---|---|
| 10 | 10.8 / 9.1 | 11.8 / 10.2 | 12.6 / 12.4 | — |
| 30 | 51.8 / 18.9 | 53.2 / 20.8 | 56.2 / 27.5 | 60.2 / 33.2 |
| 100 | 770 / 57 | 779 / 67 | 786 / 90 | 831 / 119 |

**10 % perdidos** (≈ un patrón por fila desde J = 30; µs por gradiente):

| J | M = 2 Chol / Wood | M = 3 | M = 5 | M = 8 |
|---|---|---|---|---|
| 10 | 401 / 163 | 381 / 194 | 373 / 281 | — |
| 30 | 18 634 / 1 001 | 18 495 / 1 254 | 17 376 / 1 942 | 18 559 / 3 424 |
| 100 | 261 710 / 1 496 | 274 355 / 1 903 | 256 081 / 3 120 | 255 887 / 6 147 |

- Woodbury nunca es más lento: con J = 10 y datos completos gana entre 1,02 y
  1,18 veces, y con J = 100 entre 7 y 14 veces. Con valores perdidos gana entre 1,3 y
  175 veces. **Regla de `auto`: siempre Woodbury.** Cholesky queda como motor
  interno (`options(BayesianEFA.engine = "cholesky")`) y como referencia de
  los tests.
- **FIML vectorizado.** La primera versión de Woodbury recorría los patrones
  con unas doce operaciones pequeñas por patrón (copias de `B[idx]` y
  `W[idx]`), y con un patrón por fila ese coste fijo dominaba: un ajuste de
  30 ítems, 3 factores y 10 % perdidos tardó 47 min (10 ms por paso). Ahora
  todo lo que crece con J se calcula para todos los patrones a la vez con
  productos de matrices (indicadores de patrón × productos de columnas de B;
  raíces apiladas × W), y por patrón solo queda la Cholesky de M × M. Ese
  gradiente baja de 4,6 ms a 1,25 ms; frente a Cholesky, de 14,8 veces.
- Con 30 % perdidos, la versión vectorizada solo se volvió a medir en
  J = 10, M = 2 (1 210 frente a 563 µs); el resto de esa rejilla es de la
  versión anterior y se descartó.

## 2. Precisión de Woodbury (`bench-precision.R`, `precision.csv`)

Prior normal (ψ es un parámetro con cota inferior, cuya transformación
logarítmica es exacta hasta 1e−16 relativo), J = 10, M = 2, N = 500, con un
patrón y con 82 patrones. El ψ más pequeño va de 1e−1 a 1e−12. Referencia: la
log-verosimilitud en aritmética MPFR de 120 bits.

| ψ mín | error Cholesky | error Woodbury | dif. relativa del gradiente |
|---|---|---|---|
| 1e−3 | 1e−11 | 1e−10 | 1e−13 |
| 1e−6 | 7e−10 | 1.4e−7 | 3e−13 |
| 1e−9 | 1.4e−6 | 1.6e−4 | 3e−13 |
| 1e−12 | 5e−4 | 5e−2 | 4e−13 |

El error de Woodbury crece como 1e−16/ψ, pero con ψ = 1e−12 la densidad vale
−1,2e12, así que el error relativo es 4e−14; los gradientes coinciden hasta
1e−12 en todo el rango. Los priors de ψ (y la cota 1 − h² < 1) dejan esos
valores en zonas de probabilidad despreciable. Woodbury es seguro por defecto.

## 3. NUTS (`bench-nuts.R`, `nuts.csv`)

Tres conjuntos de datos (HS_data con M = 3; simulados J = 30, M = 3, N = 500;
J = 12, M = 4, N = 150), cmdstanr, cuatro cadenas de 1000 + 1000, dos semillas.

| datos | ESS/s con δ = 0,80 | 0,90 | 0,95 | profundidad máx. |
|---|---|---|---|---|
| HS | 201 | 171 | 168 | 9 |
| J30 | 163 | 257 | 101 | 10 |
| J12 | 192 | 154 | 134 | 9 |

- Ninguna transición divergente con ningún δ.
- δ = 0,95 da el menor ESS por segundo en los tres casos (1,2–2,6 veces menos
  que el mejor). **`adapt_delta` pasa a 0,9 por defecto**, que mantiene margen
  para posteriores más difíciles.
- La profundidad del árbol nunca supera 10. **`max_treedepth` sigue en 15**:
  el tope no cuesta nada mientras no se alcanza.

## 4. Backends y métodos (`bench-backends.R`, `backends.csv`, `methods.csv`, `bench-backends-large.R`, `backends-large.csv`)

Valores por defecto de `befa()`, caché vacía en cada backend: el primer ajuste
incluye la compilación (stanli no compila).

| backend | HS: primero / en caché (s) | J30: primero / en caché (s) | ESS por s de cadena (HS / J30) |
|---|---|---|---|
| stanli | 3.5 / 3.4 | 7.7 / 7.7 | 346 / 244 |
| cmdstanr | 37.6 / 6.5 | 19.1 / 13.9 | 146 / 114 |
| stanr | 35.9 / 5.1 | 13.4 / 13.3 | 190 / 112 |
| rstan | 66.8 / 5.9 | 12.1 / 11.6 | 144 / 134 |

Con un modelo mayor (J = 60, M = 4, N = 1000), cmdstanr tarda 98 µs por paso
y stanli 113 µs (31 s frente a 36 s en total).

**Orden de `auto`: stanli → cmdstanr → stanr → rstan.** stanli no necesita
compilador, su primer ajuste tarda segundos, y hasta 30 ítems es igual o más
rápido que los backends compilados ya en caché; con 60 ítems cmdstanr es un
13 % más rápido por paso.

**ADVI y Pathfinder frente a MCMC** (error máximo en la media de h²):

- ADVI: 0,04–0,05, con k de Pareto 1,05–1,22 (> 0,7: la aproximación no es
  fiable para estos modelos, y `befa()` lo avisa).
- Pathfinder: 0,01–0,08.
- El k = `Inf` de Pathfinder con cmdstanr en `methods.csv` era un fallo
  nuestro: CmdStan devuelve draws remuestreados con reemplazo, y los empates en
  la cola dejan el ajuste de Pareto sin definir. Ahora cmdstanr y stanr
  devuelven todos los candidatos, el k se calcula sobre ellos y BayesianEFA
  hace el remuestreo PSIS (`.pathfinder_resample()`). Con HS_data, el k pasa a
  1,50.

## Ideas para más adelante

- Los árboles de NUTS son profundos (100–230 pasos por iteración). Una métrica
  densa podría reducirlos; habría que medirla.
- FIML con muchos patrones sigue costando unas 60 veces más por gradiente que
  los datos completos (1,25 ms frente a 21 µs con J = 30). El bucle por patrón
  que queda podría agruparse por número de columnas de la raíz.
