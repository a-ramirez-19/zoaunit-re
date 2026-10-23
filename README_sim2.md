# Simulación 2: escenarios del director sobre el diseño de la simulación 1

Reemplaza la simulación con generador común (antigua simulación 3, `_archivo/generador_comun/funciones_sim3.R`, 18/10/2026). Es una extensión de la simulación 1 (en las diapositivas del 21/10, «simulación 2»): para cada familia, enlace y τ se genera un conjunto de datos propio con ese mismo modelo y se ajusta el modelo correcto, así que todos los parámetros tienen valor verdadero. A ese diseño se le agregan los escenarios que pidió el director el 21/10/2026, en un diseño factorial completo.

**Numeración y nombres de archivos (22/10/2026).** Los archivos siguen la numeración del trabajo de grado. Esta es la simulación 2 (antes `sim3e`). La simulación 1 es la que en las diapositivas del 21/10 se llamó simulación 2 (`ZOAUNIT_RE_simulacion1*.qmd`, antes `_v4`; resultados en `resultados_sim1*`). El diseño descartado de generador común (antigua simulación 3) está en `_archivo/generador_comun/`.

Por qué se eligió este diseño y no el de generador común: documento "Simulación 2: decisiones de diseño" (22/10/2026).

## Matriz

| Factor | Niveles |
|---|---|
| Familia | Unit-Gamma (`ug`), Unit-Weibull (`uw`), Kumaraswamy (`ku`) |
| Enlace | logit, probit, cloglog |
| τ | 0.25, 0.50, 0.75 (cada uno con sus propios datos) |
| Sujetos (4 observaciones cada uno) | 30, 50, 100, 300, 500 (n = 200 se quitó el 22/10/2026 para ahorrar cómputo) |
| Inflación (p₀, p₁) | `i1008` = (0.10, 0.08), la de la simulación 1 (datos: 0.098 y 0.081) · `i2030` = (0.20, 0.30) · `i3020` = (0.30, 0.20) |
| σ_b | `s05` = 0.5 · `s10` = 1.0 (simulación 1) · `s15` = 1.5 |
| Forma φ | `fB` (poca dispersión) · `fM` (la de la simulación 1) · `fA` (mucha dispersión) |
| Réplicas | R = 500, sin piloto (decisión del 22/10/2026), lanzadas en tramos de 100 |

Son 27 escenarios por familia, enlace y τ. El escenario `i1008_s10_fM` es el de la simulación 1.

## Generador

Para el sujeto i y la observación j:

1. ηᵢⱼ = β₀ + β₁·sexo + β₂·edad + bᵢ, con bᵢ ~ N(0, σ_b²) y β = (0.50, 0.60, 0.40), los de la simulación 1.
2. θᵢⱼ = g⁻¹(ηᵢⱼ), con el enlace g del ajuste. Es el cuantil τ verdadero de la componente interior.
3. Componente interior: log(−log Y) = log(−log θᵢⱼ) − c_d(τ, φ) + c_d(U, φ), con U ~ U(0, 1). Así el cuantil τ de Y es exactamente θᵢⱼ.
4. Las masas en 0 y 1 salen del logit multinomial, con γ₁ = 0.35 y δ₁ = −0.35 para z.

**Valores verdaderos:** β, γ, δ, φ y σ_b, y los cuantiles θ_τ y Q_Y(τ) en 4 perfiles de covariables.

**Semillas:** las de la simulación 1 (`20260927 + 100000 · réplica + 10000 · τ + n`), en el mismo orden de sorteos. Son iguales para todas las familias, enlaces y escenarios (números aleatorios comunes), de modo que las diferencias entre escenarios no se confunden con diferencias entre muestras.

## Por qué esos valores de φ

φ no significa lo mismo en las tres familias. Por eso los niveles se fijan con una medida común: la desviación estándar de log(−log Y) dentro del sujeto (dado el predictor), que solo depende de φ. Los φ de la simulación 1 vienen de los ajustes a los datos periodontales y coinciden con los del ajuste con tipo de diente (3.13 / 1.97 / 3.66). Dan una DE de 0.61 (Unit-Gamma), 0.64 (Unit-Weibull) y 0.60 (Kumaraswamy). Los otros dos niveles quedan a la misma distancia en escala logarítmica (× 2/3 y × 3/2): DE de 0.40 y 0.90.

| Familia | fB (DE 0.40) | fM (simulación 1) | fA (DE 0.90) | Previa de φ |
|---|---|---|---|---|
| Unit-Gamma | 6.7 | 3.2 | 1.7 | Gamma(2, 0.50) |
| Unit-Weibull | 3.2 | 2.0 | 1.4 | Gamma(2, 0.80) |
| Kumaraswamy | 11.1 | 3.8 | 1.7 | Gamma(2, 0.39) |

`niveles_phi()` en `funciones_sim2.R` reproduce esta tabla con la DE exacta y el percentil de cada φ en su previa.

## Previas

- **Fijadas por el director (comentario 32):**
  - β, γ y δ ~ N(0, 2.5²);
  - σ_b ~ half-normal(0, 1).
- **φ:** Gamma(2, tasa), con la mediana en la media geométrica de los niveles fB y fA. Así los tres niveles caen entre los percentiles 14 y 93 de la previa.
- **Punto a vigilar en los resultados:** σ_b = 1.5 queda en el percentil 87 de la half-normal(0, 1). Si con n pequeño aparece sesgo hacia abajo en σ_b, la variante de sensibilidad `previa_alt1` (escalas × 2) permite ver si viene de la previa.

## Protocolo de convergencia

Es el del director (comentario 23):
- 4 cadenas;
- R̂ ≤ 1.01;
- ESS bulk y tail ≥ 400.

Si no se cumple, las cadenas se extienden hasta 4 veces. La simulación 1 usaba 2 cadenas.

## Archivos

- **`funciones_sim2.R`:** escenarios, generador, valores verdaderos, modelo JAGS, ajuste con caché y métricas.
- **`run_bloque_sim2.R`:** corre un bloque de ajustes (GitHub Actions) y deja en `registro/` el estado de cada ajuste y el entorno de cómputo.
- **`.github/workflows/simulacion_sim2.yml`:** corre en `ubuntu-24.04` y publica los `.rds` en un release con la etiqueta indicada. Si se relanza con la misma etiqueta, se saltan los ajustes ya publicados.
- **`subir_y_lanzar_sim2.ps1`** (portátil): sube los archivos y lanza una corrida por familia y enlace.
- **`traer_sim2_github.ps1`** (portátil): descarga los releases `sim2-*` a `_cache_sim2/`.
- **`resumir_sim2.R`** (portátil): resume la caché en `resultados_sim2/`. Las tablas son convergencia, estabilidad, métricas, criterio de paso y efecto de cada factor.

## Qué guarda cada ajuste

El resumen se calcula con todas las extracciones:
- media, mediana, DE e intervalo de credibilidad al 95 %;
- R̂ y ESS bulk y tail;
- los 9 parámetros, y θ_τ y Q_Y(τ) en 4 perfiles.

Además guarda las extensiones de las cadenas, el error de JAGS si lo hubo, el chequeo del truco de los ceros y las semillas.

Las trazas (250 extracciones por cadena) solo se guardan en las réplicas 1–5. Con trazas, cada ajuste ocupa ~70 KB; sin ellas, ~2.5 KB. El estudio completo ocupa unos 7 GB.

## Cómputo estimado

Los tiempos medios por ajuste (4 cadenas) salen del piloto de la simulación con generador común en GitHub Actions: 0.34, 0.51, 0.97, 3.12 y 5.64 minutos para n = 30, 50, 100, 300 y 500. Cada trabajo de GitHub Actions usa 4 núcleos; la tabla da el reloj con 20 y con 80 trabajos simultáneos (80 y 320 núcleos).

| Corrida | Ajustes | Horas-núcleo | Reloj, 20 trabajos | Reloj, 80 trabajos |
|---|---|---|---|---|
| Cada tramo de 100 réplicas (9 corridas) | 364 500 | ~12 900 | ~6.7 días | ~1.7 días |
| Total, R = 500 (5 tramos, 45 corridas) | 1 822 500 | ~64 300 | ~33 días | ~8.4 días |
| Sensibilidad (n = 30 y 100, réplicas 1–20) | 29 160 | ~320 | ~4 horas | ~1 hora |

Sin n = 200 el cómputo baja un 16 %. n = 500 pesa el 53 % del tiempo. Al reloj hay que sumarle un margen de uno o dos días, porque al final de cada tramo quedan núcleos libres mientras terminan los últimos trabajos.

Las corridas se lanzan en orden de tramo, así que cada ~6.7 días (20 trabajos) o ~1.7 días (80 trabajos) se completan 100 réplicas más en todo el diseño.

Si un trabajo llega al límite de 5 horas, los ajustes que no alcanzó quedan como `sin_tiempo`. Basta relanzar el mismo tramo (`-Replicas 101-200`, por ejemplo): se saltan los ajustes ya publicados.

## Criterio de calidad (propuesta)

Lo calcula `resumir_sim2.R` con las réplicas que haya en ese momento (`criterio_paso.csv`). Se evalúa para cada familia, enlace y escenario:
1. converge al menos el 90 % de las réplicas en cada τ y n;
2. no hay réplicas inestables;
3. se cumple el truco de los ceros;
4. desde n = 100, la cobertura de los 9 parámetros queda dentro de la banda Monte Carlo de 0.95 o por encima.

## Pasos

Las corridas se lanzan en el repositorio de la organización, `Sinu-Labs/zoaunit-re` (80 trabajos simultáneos). Los dos scripts de PowerShell aceptan `-Repo a-ramirez-19/zoaunit-re` para usar el repositorio personal. Requisitos en el portátil: git con permiso de push y GitHub CLI con sesión (`gh auth login`).

1. **Prueba corta** (unos 15 minutos): `powershell -NoProfile -ExecutionPolicy Bypass -File subir_y_lanzar_sim2.ps1 -Prueba`. Sube los archivos a un clon en `%USERPROFILE%\zoaunit-re-sinu` y lanza Unit-Gamma + logit con dos escenarios, n = 30 y réplicas 1–2 (release `sim2-prueba`). Revisar con `gh release view sim2-prueba -R Sinu-Labs/zoaunit-re` que tenga `lista_*.txt` con 12 ajustes y un `cache_*.tar`.
2. **Todo el diseño**: `powershell -NoProfile -ExecutionPolicy Bypass -File subir_y_lanzar_sim2.ps1`. Muestra el resumen (45 corridas, 1 822 500 ajustes, ~8.4 días con 80 trabajos) y pide confirmación. Las etiquetas son `sim2-<familia>-<enlace>-r<ini>-<fin>`. Para agregar la sensibilidad a las previas, añadir `-Sensibilidad`.
3. A medida que terminan los tramos: `powershell -NoProfile -ExecutionPolicy Bypass -File traer_sim2_github.ps1` (descarga todos los releases `sim2-*` menos la prueba) y luego `Rscript resumir_sim2.R`. Es incremental: solo lee los ajustes nuevos.
4. Revisar `resultados_sim2/criterio_paso.csv` y `estabilidad.csv` con el primer tramo (R = 100). Si algo falla, se pueden cancelar las corridas pendientes con `gh run cancel` antes de gastar el resto del cómputo.

GitHub cancela una corrida que pasa 35 días desde que se lanza, contando la espera en cola. Con 80 trabajos simultáneos todo termina en ~8.4 días; con 20 habría que lanzar los tramos en dos veces (`-Replicas 1-300` y luego `-Replicas 301-500`).
