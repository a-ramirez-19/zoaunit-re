# Simulación 3 en GitHub Actions: diseño de simulación del director

Implementa las correcciones del director a la propuesta (13/10/2026):

- **Generación:** una distribución generadora con parámetros conocidos por escenario, común a los órdenes cuantílicos que se comparan. Es la distribución de la componente interior: Unit-Weibull, Unit-Gamma o Kumaraswamy. Las masas en 0 y 1 se generan aparte con el logit multinomial. Sobre cada dataset se ajustan τ = 0.25, 0.50 y 0.75.
- **Matriz de escenarios:**
  - sujetos 50, 100, 200 y 300, con 4 observaciones cada uno (extendido a n = 30 y n = 500 el 18/10/2026);
  - familia generadora (componente interior) Unit-Weibull y Unit-Gamma, con Kumaraswamy como segunda distribución;
  - enlaces logit, probit y cloglog;
  - inflación baja y moderada, y heterogeneidad baja y alta;
  - valores generadores controlados y estimados con los datos reales;
  - piloto de 50 réplicas; estudio principal con R = 500 (decisión del 18/10/2026).
- **Previas:** β y γ ~ N(0, 2.5²); σ_b ~ half-normal(0, 1); φ con una previa positiva específica de cada familia.
- **Convergencia:** al menos 4 cadenas, R̂ ≤ 1.01, ESS bulk y tail ≥ 400.

Los valores numéricos que el director no fijó están en el bloque `POR_DEFINIR` de `funciones_sim3.R`. `run_bloque_sim3.R` se detiene si falta alguno.

## Archivos

- `funciones_sim3.R`: generador, modelo JAGS, ajuste con caché y resúmenes. Es el mismo archivo que usa el Notebook `ZOAUNIT_RE_simulacion_v5_comun.qmd`.
- `run_bloque_sim3.R`: corre un bloque de ajustes.
- `.github/workflows/simulacion_sim3.yml`: al terminar publica los `.rds` en un release con la etiqueta indicada, en partes de hasta 4000 ajustes. Si se relanza con la misma etiqueta, se saltan los ajustes ya publicados.

## Piloto (50 réplicas de la matriz)

Se lanza una corrida por cada distribución interior que tenga todos sus valores definidos. Cada corrida usa los enlaces logit, probit, cloglog y loglog, τ = 0.25, 0.50 y 0.75, n = 30, 50, 100, 200, 300 y 500 y los 8 escenarios de la matriz, es decir, 28 800 ajustes (7 200 por enlace).

| Etiqueta | Distribución interior |
|---|---|
| `sim3-piloto-uw-logit-probit-cloglog-loglog` | Unit-Weibull (primera etapa del comentario 30) |
| `sim3-piloto-ug-logit-probit-cloglog` | Unit-Gamma |
| `sim3-piloto-ku-logit-probit-cloglog` | Kumaraswamy |

## Niveles elegidos (18/10/2026, opción A)

- Inflación baja (p₀, p₁) = (0.05, 0.04) y moderada = (0.10, 0.08), que es lo observado en los datos.
- Heterogeneidad baja σ_b = 0.30 y alta σ_b = 1.00; la estimada con los datos (0.63) queda en medio.
- Valores estimados con los datos reales (β, φ por familia): ajuste log(−log) de `ajuste_real_loglog.R`.

## Análisis de sensibilidad

`subir_y_lanzar_sim3.ps1 -Sensibilidad` (o `-SoloSensibilidad`) lanza, con n = 100, réplicas 1–20, los tres τ
y los ocho escenarios, una corrida por variante (entrada `variante` del workflow):

| Variante | Qué cambia | Release |
|---|---|---|
| `previa_alt<k>` | previas alternativas `sens_sd_coef`, `sens_escala_sigma`, `sens_previa_phi_*` | `sim3-sens-<fam>-previa_alt<k>` |
| `masas_alt` | inflación (0.10, 0.08) y (0.20, 0.15) (opción B) | `sim3-sens-<fam>-masas_alt` |
| `sigma_alt` | heterogeneidad 0.63 y 1.26 (opción C) | `sim3-sens-<fam>-sigma_alt` |

La variante `director` (la base) reutiliza los ajustes del piloto. El Notebook compara las variantes en sesgo,
RMSE y cobertura, el cambio de la media posterior en DE posteriores y la convergencia, con trazas.

## Estudio principal (R = 500)

`subir_y_lanzar_sim3.ps1 -Etapa principal` lanza las réplicas 51 a 500 (las 1 a 50 son las del piloto, con las
mismas claves) en 256 trabajos, publicadas en `sim3-principal-<familia>-logit`: 64 800 ajustes por familia y
enlace. En GitHub se guarda cada ajuste en modo liviano: el resumen (estimaciones, intervalos, R̂ y ESS
calculados con todas las extracciones) y 250 extracciones por cadena para las trazas.

## Enlaces (18/10/2026)

El piloto ajusta logit, probit y cloglog sobre los mismos datasets (release `sim3-piloto-<familia>-logit-probit-cloglog`,
64 trabajos). El estudio principal queda con logit hasta decidir con el piloto. Ninguno de los tres enlaces coincide
con el del generador, log(−log θ): la comparación entre enlaces se hace en θ_τ y Q_Y(τ).

Desde el 18/10/2026 el piloto incluye además el enlace `loglog`, g(θ) = log(−log θ), que es el del generador común:
con él el modelo está bien especificado y β_τ = β + c_d(τ, φ)·e₁ y σ_b también tienen valor verdadero. Logit,
probit y cloglog miden la robustez frente a un enlace mal elegido.
