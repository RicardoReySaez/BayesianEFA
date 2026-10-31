# BayesianEFA v0.2.0 — Plan de reimplementación

Borrador para discutir · 2026-10-31 · base: `24c15ac` (v0.1.0) · rama de trabajo: `dev`

**Estado:** fase 0 ✔ (`3ad095b`) · fase 1 ✔ (`8e3c442`) · fase 2 ✔ (`9b6eaf1`) · fase 3 ✔ (`1923239`) · fase 4 ✔ (`a0d721c`) · fase 5 ✔ (`64b830c`) · fase 6 ✔ (este commit) · siguiente: fase 7.

**Fase 5 (2026-10-31).** `befa(HS_data, 3)` con stanli reproduce v0.1.0. Las medias de h² coinciden hasta 0,02 y las cargas alineadas hasta 0,05 con la tabla de la documentación de v0.1.0. El ajuste no da ningún aviso de convergencia.

**Fase 6 (2026-10-31).** Puntuaciones factoriales, fiabilidad, índices de ajuste y `summary()` leen los draws finales.
- El C++ recibe Λ, ψ y ν por draw (cadenas apiladas) y reproduce las fixtures de v0.1.0 hasta 1e−10, números aleatorios incluidos.
- El modelo nulo se ajusta con el backend y el método del ajuste. Con `HS_data` y stanli, el ajuste completo (puntuaciones, omega e índices) tarda 7 s.
- Los tests unitarios usan un objeto `befa` falso (`helper-befa.R`), sin Stan, para que covr los mida.

### Lo aprendido en la fase 4 (2026-10-31)

- **Procesos aparte.** Los backends no se cargan nunca en la sesión del usuario. Compilan y muestrean en procesos `Rscript` aparte, lanzados con `callr::rscript_process`, que escriben su salida en ficheros que la sesión principal lee para la barra. En Windows, rstan no carga si stanr ya está cargado en la misma sesión (choque de DLLs de TBB).
- **rstan no compilaba en este equipo.** CRAN distribuye rstan 2.32.7 con StanHeaders 2.39.1, una pareja que falla en el propio C++ de rstan. El registro lo detecta (`.rstan_mismatch()`) y marca rstan como no disponible, con el motivo. Se resolvió instalando rstan 2.39.0.9000 desde el r-universe de stan-dev y manteniendo StanHeaders 2.39.1 de CRAN. rstanarm, blavaan y metaBMA siguen funcionando.
- **El ADVI de rstan y el de stanr no devuelven `lp_approx__`** (stanr tampoco `lp__`). En la fase 5, el k de Pareto solo se podrá calcular cuando estén las dos.
- **Progreso por backend.** rstan, cmdstanr y stanli imprimen `Iteration: N`. stanr no imprime nada en un proceso de fondo, así que su barra avanza cadena a cadena.
- **Valores iniciales.** stanli quiere inits sin restringir; el trabajo los convierte con `stanli::unconstrain()`. El Pathfinder de cmdstanr pide un init por camino. stanli ejecuta un solo camino.
- **Nombres en ADVI y Pathfinder.** cmdstanr y stanli devuelven `lp__` y `lp_approx__`. El ADVI de stanr solo trae `lp_approx__`. `.standardize_draws()` unifica los nombres (`log_p__` y `log_g__` de rstan) y quita `path__`.
- **Tiempos con un modelo de 6 ítems y 2 factores** (`cor-uv-single`, N = 200):

  | Backend | Compilación (una vez) | 2 cadenas × 600 iteraciones | ADVI | Pathfinder |
  |---|---|---|---|---|
  | stanli | — | 1,2 s | 0,8 s | 0,8 s |
  | cmdstanr | 33 s | 2,7 s | 2,3 s | 2,3 s |
  | stanr | 30 s | 2,1 s | 1,8 s | 1,6 s |
  | rstan | 58 s | 2,9 s | 2,1 s | — |

---

## 0. Decisiones tomadas

| # | Decisión |
|---|---|
| D1 | Sin modelos Stan precompilados: fuera `rstantools`, `inst/stan/`, `src/stanExports_*`. El código Stan se genera en R según la configuración del modelo. |
| D2 | Cada variante se compila una sola vez y se reutiliza entre sesiones (caché persistente). |
| D3 | Backends: `rstan`, `cmdstanr`, `stanr`, `stanli`. |
| D4 | Métodos: MCMC (NUTS), ADVI (meanfield / fullrank) y Pathfinder. Sin ML ni Laplace. |
| D5 | Todas las densidades llevan sus constantes aditivas, incluida la del término de repulsión. |
| D6 | Los parámetros se devuelven en la métrica del modelo y la otra se calcula a demanda. Sufijo `_std` = métrica estandarizada. |
| D7 | La formulación de Stan se elige por benchmark: una sola Cholesky con raíz de la matriz de dispersión, Woodbury con control interno, FIML agrupado por patrones, `rows_dot_self`. |
| D8 | En C++: RSP (Hungarian), varimax por draw, puntuaciones factoriales, métricas de ajuste, log-verosimilitud del modelo nulo y EM del saturado. |
| D9 | `rsp_align()` vuelve a la versión v0.1.0 (la del paper, solo varimax). Único cambio: el varimax por draw pasa a C++ con resultados idénticos. |
| D10 | CRAN como objetivo futuro: dependencias y caché compatibles desde ya. |
| D11 | Estilo de código: el de uSDT, tu último paquete (decidido el 2026-10-31, sustituye al de v0.1.0). Ver §15. |
| D12 | RSP siempre sobre `Lambda_std`, también con el prior normal. |
| D13 | ω siempre en métrica estandarizada. |
| D14 | `backend = "auto"` por defecto; el orden de preferencia se fija tras el benchmark. |
| D15 | Se trabaja en la rama `dev` (hoy igual a `main` y a `origin/dev`, en `24c15ac`): un commit por fase, siempre con visto bueno, y push solo cuando se pida. |
| D16 | Desaparece `befa$stanfit`. Versión 0.1.0.9000 durante el desarrollo y 0.2.0 al publicar. |
| D17 | El prior normal de las cargas se escala por la SD observada de cada ítem. |
| D18 | `ordered = TRUE` se mantiene como argumento y da error de "no implementado". |
| D19 | Valores iniciales con la misma estrategia que lavaan y `stats::factanal()`: EFA por ML a partir de ψ⁰ = (1 − M/2J)/diag(R⁻¹) (§4.5). |
| D20 | Barra de progreso y mensajes de convergencia como en uSDT (rama `dev`), iguales en los cuatro backends. Los diagnósticos se calculan solo sobre los parámetros finales, ya alineados (§4.6). |
| D21 | Cada funcionalidad entra con sus tests en el mismo commit. El paquete pasa por Codecov, con la misma configuración que uSDT (§10.1). |
| D22 | Con el prior normal, la escala de la InvGamma de cada unicidad se multiplica por s_j² (decidido el 2026-10-31). Así todo el prior describe la solución estandarizada. |

---

## 1. Arquitectura

### 1.1 Flujo de `befa()`

```
befa(data, n_factors, model, lambda_prior, missing, method, backend, ...)
 │
 ├─ 1. validate_*()            entradas, priors, backend, método, argumentos
 ├─ 2. befa_spec()             estructura del modelo (lo único que cambia el código Stan)
 ├─ 3. build_stan_code(spec)   genera el código Stan → texto determinista
 ├─ 4. build_stan_data(...)    patrones, raíces, priors → lista de datos
 ├─ 5. backend_compile()       caché: solo compila si el hash no existe
 ├─ 6. backend_run(method)     MCMC / ADVI / Pathfinder → befa_draws (formato común)
 ├─ 7. reconstruct_params()    Z, h2 → Lambda_std, Psi_std (C++)
 ├─ 8. align_draws()           rsp_align() sobre Lambda_std
 ├─ 9. to_model_metric()       métrica del modelo y nombres finales
 └─ 10. opcionales             puntuaciones factoriales, fiabilidad, índices de ajuste
```

### 1.2 Ficheros en `R/`

| Fichero | Contenido |
|---|---|
| `befa.R` | Función principal (orquesta el flujo). |
| `utils_input.R`, `utils_priors.R` | Validaciones (se conservan y se adaptan). |
| `model-spec.R` | `.befa_spec()`: lista con `model_type`, `loadings` (`"uv"`, `"uni"`, `"normal"`), `patterns` (`"single"`, `"multiple"`), `engine` (`"cholesky"`, `"woodbury"`) y `null_model`. |
| `stan-codegen.R` | Componentes (cada uno aporta trozos a los bloques de Stan) y ensamblador. |
| `stan-data.R` | Patrones de datos perdidos, raíces de dispersión, priors. |
| `backend-registry.R` | Capacidades, disponibilidad, normalización de argumentos. |
| `backend-rstan.R`, `backend-cmdstanr.R`, `backend-stanr.R`, `backend-stanli.R` | Adaptadores. |
| `stan-cache.R` | Directorio, claves, `befa_compile_models()`, `befa_clear_cache()`. |
| `fit-chains.R` | Cadenas en sesiones de fondo, barra de progreso e informe de convergencia (§4.6). |
| `postprocess.R` | Reconstrucción, alineamiento y conversión de métrica. |
| `befa_factor_scores.R`, `befa_reliability.R`, `befa_fit_measures.R`, `summary_befa.R`, `utils_posterior.R` | Funciones de post-procesado, que leen `object$draws`. |
| `rsp_align.R` | v0.1.0 con el varimax en C++. |

Los ficheros que ya existen se renombran con guiones al reescribirlos (por ejemplo, `befa_factor_scores.R` → `befa-factor-scores.R`), según el estilo de uSDT (§15).

Se eliminan: `R/stanmodels.R`, `inst/stan/`, `inst/include/stan_meta_header.hpp`, `src/stanExports_*`, `configure`, `configure.win` y el truco de `Makevars.win`.

---

## 2. El modelo Stan generado

### 2.1 Notación

J ítems, M factores, P patrones de datos observados. Para el patrón p:

- n_p filas e índices observados o_p, con J_p = |o_p|;
- una raíz R_p (J_p × K_p) tal que R_p R_pᵀ = Σ_{i∈p} (y_i − c_p)(y_i − c_p)ᵀ, con c_p = ȳ_p en `raw` y c_p = 0 en `cov`/`cor`.

### 2.2 Componentes de cargas y escala

| Componente | Parámetros Stan | Construcción (local, no se guarda) | Prior (con constantes) |
|---|---|---|---|
| `uv` (M ≥ 2) | `h2` ∈ (0,1)^J, `Z` (J×M) | r²_j = ‖z_j‖²; Λ_std = diag(√(h2 ./ r²))·Z; Ψ_std = 1 − h2 | Beta(h2 \| a, b); Σ_j [½ξ·log r²_j − ½ r²_j] − J·log C(ξ, M) |
| `uni` (M = 1 con UV) | `lambda_std` ∈ (−1, 1)^J | Ψ_std = 1 − λ² | Beta escalada (ya lleva constantes) |
| `normal` (solo `cov`/`raw`) | `Lambda` (J×M), `psi` > 0 | Σ = ΛΛᵀ + diag(ψ), **sin σ** | Normal(λ_jm \| m·s_j, sd·s_j); InvGamma(ψ_j \| a, b·s_j²) |
| `scale` (`uv`/`uni` en `cov`/`raw`) | `sigma` > 0 | Σ = D_σ(Λ_std Λ_stdᵀ + Ψ_std)D_σ | half-t (con `lccdf`) |
| `mean` (`raw`) | `nu` | — | Normal(ν \| m, s) |

**Constante de repulsión.** La densidad de cada fila es p(z) ∝ ‖z‖^ξ·exp(−‖z‖²/2), es decir, ‖z‖² ~ χ²_{ξ+M} con dirección uniforme. Su constante es:

log C(ξ, M) = (M/2)·log π + ((ξ+M)/2)·log 2 + lgamma((ξ+M)/2) − lgamma(M/2)

- Verificada con integración numérica para M = 2…5 y ξ ∈ {0, 3, 100}: diferencias < 1e−13.
- Con ξ = 0 se recupera (M/2)·log(2π), la de la normal estándar.
- Con ξ = 100 vale ≈ 185–193 por fila, así que `lp__` se desplaza ≈ 185·J respecto a v0.1.0. Es lo esperable.

**Prior normal.** Estima Σ = ΛΛᵀ + Ψ directamente. Esto corrige la identificación de v0.1.0, que además reescalaba con σ.

- s_j es la SD observada del ítem j (con datos perdidos, la del modelo saturado por EM). Los hiperparámetros por ítem se pasan como vectores de datos, así que no hace falta recompilar.
- **El prior de ψ también se escala** con b·s_j² (D22): si ψ ~ IG(a, b), entonces c·ψ ~ IG(a, c·b), así que el prior completo no depende de las unidades del ítem.

### 2.3 Verosimilitud por patrones

ℓ_p = −½·[ n_p·J_p·log 2π + n_p·log|Σ_p| + tr(Σ_p⁻¹ R_p R_pᵀ) + n_p·(ȳ_p − ν_p)ᵀ Σ_p⁻¹ (ȳ_p − ν_p) ]

- El último término solo existe en `raw`.
- Listwise y estadísticos resumen son el caso de un único patrón.
- La raíz se calcula en R: si n_p ≥ J_p, es la Cholesky de la dispersión (con respaldo por autovalores si es singular); si n_p < J_p, son las filas centradas (R_p = Y_cᵀ).

**Motor `cholesky`.** L = chol(Σ_p); log|Σ_p| = 2·Σ log diag(L); tr = ‖L⁻¹R_p‖²_F. Coste O(J_p³/3 + J_p²·K_p).

**Motor `woodbury`.** Se escribe Σ = AAᵀ + diag(ψ), con A = D_σΛ_std y ψ = σ²Ψ_std (o A = Λ y ψ directamente con el prior normal). Entonces:

- B = Ψ^{−1/2}A y C = I_M + BᵀB = L_C L_Cᵀ;
- log|Σ_p| = Σ log ψ + 2·Σ log diag(L_C);
- tr = Σ_j ‖R_{p,j}‖²/ψ_j − ‖L_C⁻¹ AᵀΨ⁻¹ R_p‖²_F. Las normas por fila de R_p son datos y se calculan una vez en `transformed data`.
- Coste O(J_p·M·K_p + J_p·M²). Nunca forma una matriz J×J.

Los dos motores están verificados en R contra `mvtnorm::dmvnorm`: coinciden en 10 decimales, también con raíz de filas cuando n < J (diferencia 1.4e−14).

**Control del motor.** Se decide al generar el código. Por defecto lo elige una regla en J y M calibrada en el benchmark (§9); se puede forzar con `options(BayesianEFA.engine = "auto" | "cholesky" | "woodbury")`, que es interno, para tests y benchmark.

### 2.4 Salida mínima

- Stan solo escribe parámetros (`h2`, `Z`, `sigma`, `nu`, o `Lambda` y `psi`) y `lp__`.
- Λ_std, Ψ, Σ y R son variables locales del bloque `model`.
- No hay `generated quantities`: la log-verosimilitud puntual para LOO se calcula después en C++.

### 2.5 Ejemplo: `cor` + `uv` + un patrón + motor `cholesky`

```stan
data {
  int<lower=1> N;
  int<lower=1> J;
  int<lower=2> M;
  matrix[J, J] S_root;          // S_root * S_root' = sum_i y_i y_i' (computed in R)
  vector<lower=0>[2] pr_h2;     // Beta(a, b) on communalities
  real<lower=0> pr_xi;          // repulsion
}
transformed data {
  real ll_const = -0.5 * N * J * log(2 * pi());
  real log_C = 0.5 * M * log(pi()) + 0.5 * (pr_xi + M) * log2()
               + lgamma(0.5 * (pr_xi + M)) - lgamma(0.5 * M);
}
parameters {
  vector<lower=0, upper=1>[J] h2;
  matrix[J, M] Z;
}
model {
  vector[J] r2 = rows_dot_self(Z);
  matrix[J, M] Lambda_std = diag_pre_multiply(sqrt(h2 ./ r2), Z);
  matrix[J, J] L = cholesky_decompose(add_diag(tcrossprod(Lambda_std), 1 - h2));

  // Priors (normalised)
  target += beta_lpdf(h2 | pr_h2[1], pr_h2[2]);
  target += 0.5 * pr_xi * sum(log(r2)) - 0.5 * sum(r2) - J * log_C;

  // Likelihood (normalised)
  target += ll_const - N * sum(log(diagonal(L)))
            - 0.5 * sum(square(mdivide_left_tri_low(L, S_root)));
}
```

Comprobado el 2026-10-31:

- compila con el stanc de rstan (Stan 2.32.2);
- evaluado con stanli (`log_prob_grad`) en un punto aleatorio, su `lp` coincide exactamente con una referencia en R (`dbeta` + repulsión + `dmvnorm` + Jacobiano logit). Diferencia: 0.

### 2.6 Modelo nulo

- **`cov`/`raw`:** se genera con el mismo sistema, Σ = diag(σ²) (más ν en `raw`). Con las normas por fila de las raíces, cada patrón cuesta O(J_p).
- **`cor`:** es analítico, sin Stan, como ahora.

---

## 3. Datos para Stan (`stan-data.R`)

- **Preprocesado por tipo de modelo**, como en v0.1.0:
  - `cor`: estandarizar las columnas, ignorando los NA;
  - `cov`: centrar las columnas;
  - `raw`: datos sin tocar.
- **Agrupación por patrón de datos observados.** Las filas sin ningún dato se descartan, con aviso.
- **Por patrón:** n_p, J_p, índices, ȳ_p (en `raw`), la raíz y sus normas por fila.
- **Varios patrones:** vectores apilados (`obs_idx`, `root_flat`, `root_r2`, `ybar_flat`) cuyas longitudes se precalculan en R. El código recorre los patrones con `segment()` y `to_matrix()`.
- **Estadísticos resumen:** un único patrón con raíz chol(N·S), o chol(N·R) en `cor`.
- **Hiperparámetros:** van siempre como datos, así que cambiar un prior nunca obliga a recompilar.

---

## 4. Backends

### 4.1 Interfaz común

Cada adaptador implementa:

- `available()`, con la versión mínima;
- `compile(code, cache_key)` → objeto modelo (nada en stanli);
- `run(model, data, method, args, inits)` → `befa_draws`, una lista con:
  - `draws`: `draws_array` de iteraciones × cadenas × variables;
  - `sampler`: diagnósticos, o `NULL`;
  - `approx`: `lp_approx` o log-ratios, si los hay;
  - `method`, `backend`, `time` y `native` (el objeto nativo del backend).

### 4.2 Capacidades

| Backend | MCMC | ADVI | Pathfinder | Compila | Toolchain | Instalado aquí |
|---|---|---|---|---|---|---|
| rstan | `sampling()` | `vb()` | ✗ (error claro) | sí (lento) | sí | 2.32.7 (Stan 2.32.2) |
| cmdstanr | `$sample()` | `$variational()` | `$pathfinder()`, varios caminos + PSIS | sí | sí | 0.99.0 + CmdStan 2.40 |
| stanr | `$sample()` | `$variational()` | `$pathfinder()` | sí (C++20) | sí | no instalado |
| stanli | `sample_cstan()` | `variational_model()` | `pathfinder_model()`, 1 camino | no | no | 0.19.1 (sin ADVI ni Pathfinder: hay que pasar a ≥ 0.25) |

En cmdstanr 0.99 la API ha cambiado: por ejemplo, `cmdstan_model()` ya no acepta `compile = FALSE`. El adaptador tiene que fijar una versión mínima y probarse contra ella.

rstan compila en tiempo de ejecución contra las cabeceras de `StanHeaders`, `BH` y `RcppEigen`, que una instalación binaria de rstan no siempre trae. uSDT comprueba esto antes de compilar y nosotros haremos lo mismo.

### 4.3 Argumentos

- **Nombres canónicos**, al estilo de cmdstanr: `chains`, `parallel_chains`, `iter_warmup`, `iter_sampling`, `adapt_delta`, `max_treedepth`, `seed`, `refresh`, `init`.
- **Por compatibilidad**, también se aceptan los de rstan: `iter`, `warmup`, `cores`, `control`.
- Cada adaptador traduce estos nombres a los suyos.
- **Pathfinder:** fijamos `num_paths` y el remuestreo PSIS cuando el backend lo admite.
- **Inits:** se generan en R para cada cadena, en la escala restringida (§4.5). stanr solo admite un init común a todas las cadenas (ver riesgos).
- **`backend = "auto"`:** se usa el primer backend disponible según un orden de preferencia, que se fija tras el benchmark. Se informa con un mensaje de cuál se ha usado.

### 4.4 Compilar una sola vez

- **Clave:** md5 del código, calculado con `tools::md5sum()` para no añadir dependencias, más el backend y su versión (y la de CmdStan o Stan).
- **Directorio:** `tools::R_user_dir("BayesianEFA", "cache")/<backend>/<versión>/<md5>.stan`, junto al binario.
  - cmdstanr: el ejecutable va junto al `.stan`. Si el `.stan` ya existe no se reescribe, porque cmdstanr compara las fechas.
  - rstan: `saveRDS()` del `stanmodel`.
  - stanr: `.stan` en la caché; así su `.so` persiste entre sesiones.
  - stanli: no hace falta nada.
- **Compilación atómica**, en un directorio temporal y luego renombrando, por si dos sesiones compilan la misma variante a la vez.
- **`befa_compile_models(backend)`:** precompila todas las variantes cuando el usuario quiera. Es el equivalente a `BUILD_ALL` de hBayesDM, pero opcional y en tiempo de ejecución.
- **`befa_clear_cache()`** y poda de las versiones antiguas. CRAN lo exige para usar `R_user_dir()`.

Número de variantes: tipo de modelo (3) × cargas (`uv`/`uni`/`normal`, con restricciones) × patrones (1 o varios) × motor (2), más los nulos ≈ 30.

### 4.5 Valores iniciales

Es la misma estrategia que usan lavaan, `stats::factanal()` y `psych::fa(fm = "ml")`. lavaan, para los bloques EFA, parte de `lav_efa_extraction(method = "ML")`. Lo he comprobado en el código instalado: lavaan 0.7.3, psych 2.6.9 y stats de R 4.5.0.

1. **Matriz de correlaciones R:**
   - con datos completos, la muestral;
   - con datos perdidos, la del modelo saturado por EM, como hace lavaan con su `h1`;
   - con estadísticos resumen, la proporcionada.
2. **Unicidades iniciales:** ψ⁰_j = (1 − M/(2J)) / (R⁻¹)_jj.
   - Es la fórmula exacta de `factanal()` y de `lav_efa_ml_init_cache()`.
   - Equivale a (1 − M/2J)·(1 − SMC_j).
3. **EFA por ML:**
   - se optimiza ψ con Λ en forma cerrada, a partir de la descomposición espectral de Ψ^{−1/2}RΨ^{−1/2};
   - cotas [0.005, 1], como `factanal()`;
   - se reimplementa en R (unas 30 líneas), porque `stats:::factanal.fit.mle` no está exportada y CRAN no admite `:::`.
4. **Traducción a cada parametrización**, con un ruido pequeño distinto por cadena:
   - `uv`: h2 = clamp(1 − ψ, 0.05, 0.95), y Z_j = r*·λ_j/‖λ_j‖ con r* = √(ξ + M − 1).
     - r* es la moda de ‖z‖, cuya densidad es ∝ r^{ξ+M−1}·e^{−r²/2}; con ξ = 100 y M = 3 vale ≈ 10.1.
     - v0.1.0 inicializaba Z ~ N(0, 1), con norma ≈ √M, lejos de esa zona, y el warmup tenía que recorrer esa distancia.
   - `uni`: lambda_std = clamp(λ, −0.95, 0.95).
   - `normal`: Lambda = D_s·Λ y psi = s²·ψ, en métrica original.
   - `sigma` = SD observadas; `nu` = medias observadas.
5. **Respaldo:** si la EFA por ML falla o no converge, se recurre al eje principal de un paso con las SMC como comunalidades, como en psych. Los casos Heywood no hacen falta como respaldo, porque h2 ya se recorta a [0,05; 0,95].

### 4.6 Barra de progreso y diagnósticos (como uSDT)

Se copia el diseño de uSDT (rama `dev`, `R/fit-bayes.R`: `.bayes_header()`, `.run_chains()`, `.bar_style()`, `.bayes_report()`, `.bayes_diagnostics()`).

**Cabecera.**
- `cli::cli_rule()` con el título y el backend.
- `cli_alert_info()` con cadenas, draws, warmup y cuántas cadenas corren a la vez.
- Sustituye al banner de v0.1.0 (líneas ‗‗‗).

**Una cadena por sesión de R en segundo plano**, con `callr::r_bg()`.
- Como mucho `parallel_chains` a la vez. `.default_cores()` usa hasta 4 núcleos, deja 2 libres y respeta `_R_CHECK_LIMIT_CORES_`.
- Los mensajes de Stan se quedan en esas sesiones. El proceso principal lee sus líneas, extrae `Iteration: N` y mueve **una sola barra** (`cli_progress_bar`):
  - spinner y fase (Warmup / Sampling);
  - cuadrados ■ azul claro (`#6CB4EE`) para lo hecho y línea ─ gris (`#9AA0A6`) para lo pendiente;
  - porcentaje y ETA.
- Al terminar: "✔ Sampled 4 chains in 12.3s". En consolas no dinámicas (knitr, logs) solo aparece esa línea.
- Si se interrumpe, se matan las sesiones.

**Ventajas añadidas de lanzar cada cadena aparte:**
- la barra es la misma en los cuatro backends;
- cada cadena tiene su propio init, lo que resuelve la limitación de stanr;
- cada cadena tiene su semilla (`seed` + número de cadena);
- las cadenas se combinan como `draws_array` con `posterior::bind_draws(along = "chain")`.

**Coste:** cada sesión tarda ~1–2 s en arrancar en Windows y tiene que cargar el modelo compilado. Se notará en modelos pequeños; lo mediremos en el benchmark.

**ADVI y Pathfinder** no tienen un total de iteraciones fiable, porque paran al converger. En vez de barra con ETA, muestran un spinner con estado y tiempo transcurrido.

**Formato de progreso por backend:** rstan y cmdstanr imprimen `Iteration: N / T`. El de stanr y stanli se comprueba en la fase 4; si no se puede leer, se usa el spinner.

**Informe tras el ajuste.** Se calcula **solo sobre los parámetros finales**, los de `object$draws` ya alineados: `Lambda` o `Lambda_std`, `Psi` o `Psi_std`, `h2`, `sigma` y `nu`. Nunca sobre `Z` ni sobre las cargas sin alinear, que no están identificadas.
- ✔/! divergencias;
- ✔/! transiciones que llegan a la profundidad máxima;
- ✔/! "Largest R-hat 1.003, below 1.01";
- ✔/! "Smallest effective sample size 812, above 400" (el mínimo de ESS bulk y tail);
- en ADVI y Pathfinder, en lugar de R-hat y ESS, el k de Pareto de la aproximación (✔ si ≤ 0.7);
- si algo falla, un `warning()` que repite los problemas y su remedio (subir `adapt_delta`, `max_treedepth` o `iter_sampling`). Así lo ve también quien ejecuta con `verbose = FALSE`.

Los diagnósticos se guardan en `object$diagnostics` y `summary()` los reutiliza. El umbral de R-hat pasa a 1.01; v0.1.0 usaba 1.1.

**Dependencias nuevas en Imports:** `cli`, `callr` y `rlang` (`check_installed()`, `hash()`), todas en CRAN.

**La caché también sigue a uSDT:**
- una carpeta por modelo, que se refresca con `Sys.setFileTime()` cada vez que se usa;
- `.prune_cache()` conserva los N más recientes. uSDT usa 20; aquí habrá que subirlo, porque hay ~30 variantes por backend;
- `options(BayesianEFA.cache_dir = ...)` para moverla.

---

## 5. Objeto `befa` y nombres de parámetros

```
befa
 ├─ draws        draws_array con los parámetros finales (alineados, métrica del modelo)
 ├─ fit          befa_draws: diagnósticos, método, backend, tiempos, objeto nativo
 ├─ spec, stan_code, stan_data, priors_used, options, call
 ├─ rsp          objetivo y vectores de signo y permutación
 └─ factor_scores / reliability / fit_indices (opcionales)
```

**Cambio rompedor (decidido):** desaparece `befa$stanfit`; el objeto nativo pasa a `befa$fit$native`. Versión 0.1.0.9000 durante el desarrollo y 0.2.0 al publicar.

| Modelo | Se guardan (por defecto) | A demanda |
|---|---|---|
| `cor` + `uv`/`uni` | `Lambda_std`, `Psi_std`, `h2` | — |
| `cov` + `uv`/`uni` | `Lambda`, `Psi`, `sigma`, `h2` | `Lambda_std`, `Psi_std` |
| `raw` + `uv`/`uni` | `Lambda`, `Psi`, `sigma`, `nu`, `h2` | `Lambda_std`, `Psi_std` |
| `cov`/`raw` + `normal` | `Lambda`, `Psi` (+ `nu`) | `Lambda_std`, `Psi_std`, `h2` |

- `Sigma` y `Rho` solo se calculan a demanda.
- `h2` no lleva sufijo, porque no depende de la escala.
- `eta` solo se guarda si `factor_scores = TRUE`.
- La métrica se elige al extraer: `extract_posterior_draws(object, pars, metric = c("model", "std", "raw"))`.

---

## 6. Post-procesado

- **Reconstrucción (C++):** Λ_std se reconstruye desde Z y h2 exactamente como en Stan.
- **Alineamiento:** `rsp_align()` siempre sobre Λ_std. Después, Λ = D_s·Λ_std alineada, draw a draw, con D_s = σ_s en `uv`/`uni` y D_s = √diag(Σ_s) con el prior normal. No hace falta guardar las matrices de rotación.
- **Puntuaciones factoriales:** el C++ recibe Λ, ψ y ν por draw, en la misma métrica que los datos. Esto corrige la mezcla de escalas de v0.1.0. FIML por patrón.
- **Fiabilidad:** ω por draw, siempre con Λ_std y Ψ_std.
- **Índices de ajuste:**
  - log-verosimilitud puntual en C++ a partir de Λ, ψ y ν;
  - LOO con r_eff en MCMC y r_eff = 1 en ADVI y Pathfinder;
  - modelo nulo generado con el mismo backend y el mismo método.
- **`summary()`:** h² y u² según la métrica.
- **ADVI y Pathfinder:** no hay R-hat. Se muestra el k de Pareto de la aproximación cuando haya `lp_approx`.

---

## 7. C++

| Función | Estado |
|---|---|
| `rsp_exact_eff()` | Sin cambios. |
| `varimax_draws_cpp()` | Nueva. Mismo algoritmo que `stats::varimax(normalize = FALSE, eps = 1e-5)`. |
| `reconstruct_uv_cpp()` | Nueva. |
| `compute_factor_scores_cpp()` | Adaptada: recibe Λ, ψ y ν. |
| `compute_posterior_metrics_cpp()` | Adaptada: recibe Λ, ψ y ν, y construye Σ en cada draw. |
| `null_loglik_cpp()` | Nueva; sustituye el bucle en R. |
| `em_saturated_cpp()` | Nueva; sustituye `em_E_step()` en R. |

`Makevars` y `Makevars.win` pasan a ser los estándar de RcppArmadillo (`$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)`).

---

## 8. `rsp_align()`

- Se restaura la v0.1.0; ya lo está tras limpiar el repo.
- Único cambio: el `apply()` con `stats::varimax` pasa a `varimax_draws_cpp()`.
- Tests:
  - igualdad con `stats::varimax` hasta 1e−10;
  - una fixture con la salida de v0.1.0, generada antes de tocar nada.

---

## 9. Benchmark (`bench/`, fuera del build)

**Datos:** `HS_data` y simulados con J ∈ {10, 30, 60, 100}, M ∈ {2, 3, 5, 8}, N ∈ {300, 1000, 5000} y datos perdidos MCAR al 0/10/30 %.

**Niveles:**

1. **Coste del gradiente por formulación:** el modelo unificado de v0.1.0 frente al generado; Cholesky frente a Woodbury; un patrón frente a varios. Se mide con `log_prob_grad` en puntos aleatorios. Sale de aquí la regla de `engine = "auto"`.
2. **NUTS completo:**
   - ESS bulk y tail por segundo en h², Ψ y `lp__`, que son invariantes a la rotación;
   - divergencias y tiempo total, con la compilación aparte;
   - de aquí sale también la revisión de los valores por defecto (`adapt_delta = 0.95`, `max_treedepth = 15`).
3. **Backends y métodos:** el mismo modelo en los cuatro backends; ADVI y Pathfinder frente a MCMC (error en h², k de Pareto).
4. **Precisión:** `lp` idéntico entre motores y frente a la referencia en R; error relativo de Woodbury en función de min ψ.

Los resultados se guardan en CSV en `bench/results/`.

---

## 10. Tests

- **Código generado:** snapshots por variante (`expect_snapshot`) y comprobación de sintaxis con `rstan::stanc()` (que además vigila el suelo de Stan 2.32), `stanli_check_syntax()` o `$check_syntax()` de cmdstanr. Se saltan si no hay backend.
- **Datos para Stan:** R_p R_pᵀ coincide con la dispersión y los patrones son correctos.
- **Constantes:** log C frente a integración numérica.
- **Densidad completa:** el `lp` de Stan, obtenido con el `log_prob` del backend, frente a una referencia en R (`dbeta`, `dnorm`, `dmvnorm`, Jacobianos). Tolerancia 1e−8.
- **Equivalencia entre backends:** mismo `lp` en el mismo punto.
- **C++:** varimax frente a `stats::varimax`; fixtures de v0.1.0 para `rsp_align`, las puntuaciones y las métricas.
- **Integración:** modelos minúsculos por backend y método (`skip_if_not_installed`, `skip_on_cran`). stanli es el backend rápido de los tests.
- **Caché:** la segunda llamada no compila y una sesión nueva reutiliza el binario.

### 10.1 Cobertura (Codecov)

**Regla de trabajo:**
- Ninguna función nueva entra sin tests en el mismo commit.
- Al cerrar cada fase se ejecuta `covr::package_coverage()` en local. El resumen por fichero se reporta junto con el commit.
- Objetivo: cubrir prácticamente todas las líneas nuevas, tanto en `R/` como en `src/`. covr mide también el C++ (gcov), así que las funciones C++ necesitan tests directos.
- `# nocov` solo para ramas que de verdad no se pueden probar (por ejemplo, las que dependen de una consola interactiva), con un comentario que explique por qué.

**Cómo cubrir el código que depende de un backend:**
- **Sin backend:** traducción de argumentos, conversión de la salida a `draws_array`, capacidades, errores y caché. Se prueban con `testthat::local_mocked_bindings()` y objetos falsos que imitan la salida de cada backend.
- **Con backend real:** modelos minúsculos, con `skip_if_not_installed()` y `skip_on_cran()`. stanli no compila, así que estos tests son rápidos y se pueden ejecutar siempre en CI.
- **Barra de progreso e informe:** `cli` con `cli.dynamic = FALSE` y `expect_snapshot()` para la salida; la ruta de consola dinámica se prueba forzando `cli.dynamic = TRUE` con `withr::local_options()`.

**Infraestructura**, copiada de uSDT (fase 1):
- `codecov.yml`, con los checks de project y patch informativos y sin comentarios en los PR. Va en `.Rbuildignore`.
- `.github/workflows/test-coverage.yaml`, derivado de r-lib/actions, con `covr` → `to_cobertura()` → `codecov/codecov-action`.
  - Añadimos `dev` a las ramas que lo disparan, para ver la cobertura mientras trabajamos.
  - El workflow instala stanli desde r-universe y ejecuta `stanli::stanli_install()`, para que los tests con backend real cuenten en la cobertura.
- Badge de Codecov en el README.
- **Tienes que hacerlo tú en GitHub:** añadir el secreto `CODECOV_TOKEN` al repositorio, como en uSDT.

---

## 11. Dependencias y CRAN

- **Imports:** `posterior`, `loo`, `Rcpp`, `methods`, `stats`, `utils`, `tools`. Revisar `mvnfast` y `withr` cuando el FIML pase a C++.
- **LinkingTo:** `Rcpp`, `RcppArmadillo`.
- **Suggests:** `rstan`, `cmdstanr`, `stanr`, `stanli`, `testthat`, `knitr`, `rmarkdown`.
- **`Additional_repositories`:** los r-universe de stan-dev, andrjohns y seantalts.
- **Fuera:** `rstantools`, `StanHeaders`, `BH`, `RcppEigen`, `RcppParallel`, `brms`.
- **Reglas de CRAN:**
  - todo uso de un backend va detrás de `requireNamespace()`;
  - los ejemplos van en `\donttest{}`;
  - los tests llevan `skip_on_cran()`;
  - la caché solo se escribe en `R_user_dir()`, y en los tests se redirige a `tempdir()`.
- **Repo:** trackear `.gitignore`, y añadir `^notes$`, `^bench$`, `^\.vscode$` y `^\.Rhistory$` a `.Rbuildignore`. La carpeta de notas no se llama `dev/` para no chocar con la rama `dev` (git da "ambiguous argument").

---

## 12. Fases

| Fase | Contenido | Terminada cuando |
|---|---|---|
| 0. Preparación | Pasar a la rama `dev`; `.gitignore` y `.Rbuildignore`; instalar stanr y actualizar stanli a ≥ 0.25; fixtures de v0.1.0. | Fixtures guardadas y backends disponibles. |
| 1. Quitar el Stan precompilado | Borrar `stanmodels`, `inst/stan`, `stanExports` y `configure`; actualizar DESCRIPTION y Makevars; `codecov.yml` y `test-coverage.yaml` (§10.1). | El paquete instala en segundos y carga; los tests que no usan Stan pasan; primera medida de cobertura. |
| 2. C++ base | Varimax en C++ dentro de `rsp_align`; modelo nulo y EM en C++. | Tests de igualdad en verde. |
| 3. Generador y datos | Spec, componentes, ensamblador, patrones y raíces, constantes; motor Cholesky. | Snapshots, sintaxis y `lp` frente a la referencia en verde. |
| 4. Backends y caché | Los 4 adaptadores, capacidades, argumentos, inits y caché; cadenas en sesiones de fondo con la barra de progreso de uSDT. | Un modelo pequeño corre en los 4 backends con cada método, con la misma barra, y solo compila una vez. |
| 5. Nuevo `befa()` | Objeto nuevo, reconstrucción, alineamiento y métrica; `extract_posterior_draws()` y `posterior_summaries()`; informe de convergencia sobre los parámetros finales. | `befa(HS_data, 3)` reproduce v0.1.0: cargas alineadas iguales dentro del error de Monte Carlo. |
| 6. Post-procesado | Puntuaciones, fiabilidad, ajuste y `summary()`. | Los tests de esas funciones en verde. |
| 7. Woodbury y benchmark | Motor Woodbury; benchmark; regla `auto`; valores por defecto de NUTS; backend por defecto. | Informe del benchmark y valores por defecto fijados. |
| 8. Documentación y check | Roxygen, viñetas, README (con el badge de Codecov), NEWS; `R CMD check --as-cran`. | Check limpio. |

En todas las fases, "terminada" incluye además que las líneas nuevas estén cubiertas por tests y que el resumen de `covr` acompañe al commit (§10.1).

---

## 13. Decisiones abiertas

Las siete del primer borrador quedaron resueltas el 2026-10-31 (D12–D18). Quedan:

1. **Orden de preferencia de `backend = "auto"`:** se decide con los resultados del benchmark (fase 7).
(El estilo de código quedó resuelto: el de uSDT, D11 y §15.)

---

## 14. Riesgos

- **Woodbury** pierde precisión cuando ψ → 0 (casi Heywood). Hay que medirlo antes de activarlo por defecto.
- **rstan** fija el suelo del lenguaje en Stan 2.32 y tarda ~1 min en compilar cada variante (solo la primera vez).
- **stanli y stanr** son pre-1.0, y la API de **cmdstanr 0.99** también ha cambiado. Hay que fijar versiones mínimas y tener tests por adaptador.
- **stanr** no admite inits distintos por cadena: o se lanza cada cadena en una llamada aparte o se usa un init común.
- **Compilaciones concurrentes** de la misma variante: escritura atómica en la caché.
- **FIML en `cov`/`cor`** centra o estandariza con las medias y SD observadas de cada columna (heredado de v0.1.0). Bajo MAR es una aproximación.
- **ADVI y Pathfinder** pueden sesgar p_loo y los índices DevM: documentarlo y mostrar el k de Pareto.
- **La rama `dev` está en GitHub** (`origin/dev`): todo lo que se suba ahí es público. Nada de material privado en el árbol al hacer commit.

---

## 15. Estilo de código (el de uSDT)

Tomado de uSDT (rama `dev`: `R/utils.R`, `R/fit-bayes.R`, `R/fit-hsdt.R`, `tests/testthat/`).

**Ficheros:**
- Nombres con guiones: `fit-stan.R`, `stan-codegen.R`, `backend-cmdstanr.R`.
- Cabecera de cuatro líneas:
  ```r
  # stan-codegen.R
  # Build the Stan code of a BEFA model from its specification
  # Author: Ricardo Rey-Sáez
  # Last modified: DD-MM-YYYY
  ```
- Secciones con un comentario simple (`# Public functions`, `# Internal functions`, `# Messages`), sin banners de caja.

**Funciones públicas:** roxygen en markdown, con `@param` en prosa y `@return` como lista de `$elementos`.

**Funciones internas:**
- Nombre con punto (`.stan_block_data()`).
- Una línea encima que empieza por "This function…" y explica qué hace y por qué, en frases completas.
- Sin roxygen ni `@noRd`.

**Mensajes:**
- `.befa_stop()` y `.befa_warn()`, que equivalen a `stop()` y `warning()` con `paste0(...)` y `call. = FALSE`.
- `.befa_msg()`, que antepone "BayesianEFA: ".
- Los argumentos se citan entre comillas invertidas: "`n_factors` must be …".
- Lo que se muestra durante el ajuste (cabecera, barra, informe) va con `cli`.

**Comprobaciones:** helpers `.check_*()` que devuelven `invisible(TRUE)`.

**Idioma del código:**
- Enteros con `L` (`1L`) y `if (length(x))` en lugar de `length(x) > 0`.
- `vapply()` y `seq_len()`, `%||%`, y siempre `pkg::fun`.
- Las funciones devuelven su última expresión, sin `return()` final.
- Líneas de ≤ 80 caracteres, sangría de 2 espacios y continuaciones alineadas.

**Tests:**
- Descripciones en frase completa: `test_that("the generated code compiles under Stan 2.32", …)`.
- Se llama directamente a las funciones internas.
- Comentarios que justifican las tolerancias.
- `helper-*.R` para los datos compartidos y `_snaps/` para los snapshots.

**`rsp_align()`:** su código se mantiene como en el paper (D9). Solo recibe la cabecera de fichero y el cambio del varimax en C++.
