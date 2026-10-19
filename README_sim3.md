# Simulación 3 en GitHub Actions: diseño reformulado

Implementa las correcciones del director a la propuesta (13/10/2026) con el diseño reformulado del 18/10/2026 (documento "Simulación 3: diseño reformulado"). El objetivo es la recuperación de parámetros.

- **Generación:** una distribución generadora con parámetros conocidos por escenario, común a los órdenes cuantílicos que se comparan. Es la distribución de la componente interior. La primera etapa es Unit-Gamma; después vienen Unit-Weibull y Kumaraswamy, cada una generada y ajustada consigo misma. Las masas en 0 y 1 se generan aparte con el logit multinomial. Sobre cada dataset se ajustan τ = 0.25, 0.50 y 0.75 con los enlaces logit, probit y cloglog.
- **Matriz de escenarios:**
  - sujetos 30, 50, 100, 200, 300 y 500, con 4 observaciones cada uno;
  - inflación baja (p₀, p₁) = (0.05, 0.04) y moderada = (0.10, 0.08);
  - heterogeneidad baja σ_b = 0.30 y alta σ_b = 1.00;
  - valores generadores controlados y estimados con los datos reales;
  - piloto de 50 réplicas; estudio principal con R = 500 (réplicas 51–500 sumadas a las del piloto).
- **Previas:** β, γ y δ ~ N(0, 2.5²); σ_b ~ half-normal(0, 1); φ ~ Gamma(2, tasa) propia de cada familia (Unit-Gamma: Gamma(2, 0.65)).
- **Convergencia:** 4 cadenas, R̂ ≤ 1.01, ESS bulk y tail ≥ 400; si no se cumple, las cadenas se extienden hasta 4 veces.

## Archivos

- `funciones_sim3.R`: generador, modelo JAGS, ajuste con caché, resúmenes y métricas de estabilidad. Es el mismo archivo que usa el Notebook `ZOAUNIT_RE_simulacion_v5_comun.qmd`.
- `run_bloque_sim3.R`: corre un bloque de ajustes y deja en `registro/` el estado de cada ajuste y el entorno de cómputo del bloque.
- `.github/workflows/simulacion_sim3.yml`: corre en `ubuntu-24.04` (versión fija) y al terminar publica los `.rds` en un release con la etiqueta indicada, en partes de hasta 4000 ajustes. Si se relanza con la misma etiqueta, se saltan los ajustes ya publicados.

## Qué guarda cada ajuste

En GitHub cada ajuste se guarda en modo liviano: el resumen calculado con todas las extracciones (media, mediana, DE, intervalo, R̂ y ESS) y 250 extracciones por cadena para las trazas. Además:

- **Estabilidad:** número de extensiones, error de JAGS (la réplica queda como inestable, no se descarta) y la comprobación del truco de los ceros (cota de la log-densidad interior frente a Cz, comentario 22).
- **Reproducibilidad (comentario 34):** semilla de los datos (20260927 + 100000 · réplica + n), semilla de las cadenas, configuración MCMC, previas y entorno (versiones de R, JAGS, rjags, coda y posterior; sistema; imagen del ejecutor; repositorio, commit y corrida).

## Piloto (50 réplicas)

`subir_y_lanzar_sim3.ps1 -Familias ug -Sensibilidad` lanza el piloto de Unit-Gamma: logit, probit y cloglog, τ = 0.25, 0.50 y 0.75, n = 30 a 500 y los 8 escenarios, es decir, 21 600 ajustes (unas 6 horas de reloj), en el release `sim3-piloto-ug-logit-probit-cloglog`.

## Análisis de sensibilidad

Con `-Sensibilidad` (o `-SoloSensibilidad`) se lanza la variante `previa_alt1`: todas las escalas de las previas × 2 (N(0, 5²) para β, γ y δ; half-normal(0, 2) para σ_b; Gamma(2, 0.325) para φ en Unit-Gamma). Usa n = 30 y 100, las réplicas 1–20 del piloto (los mismos datasets), los 8 escenarios, los 3 τ y los 3 enlaces: 2 880 ajustes, en `sim3-sens-<familia>-previa_alt1`. La variante `director` reutiliza los ajustes del piloto.

## Estudio principal (R = 500)

Si el piloto cumple el criterio de paso del Notebook, `subir_y_lanzar_sim3.ps1 -Familias ug -Etapa principal` lanza las réplicas 51 a 500 con una corrida por enlace (`sim3-principal-ug-logit`, `-probit`, `-cloglog`; 64 800 ajustes cada una, 256 trabajos). En total son 194 400 ajustes, unas 54 horas de reloj.

## Criterio de paso (propuesta para el director)

El piloto pasa al estudio principal, y la familia a la etapa siguiente, si en todos los escenarios y enlaces: (a) converge al menos el 90 % de las réplicas; (b) no hay réplicas inestables; (c) se cumple el truco de los ceros; y (d) desde n = 100, la cobertura de γ, δ y φ queda dentro de la banda Monte Carlo. Ninguno de los tres enlaces coincide con la escala log(−log) del generador, así que β_τ y σ_b no tienen valor verdadero y el sesgo de θ_τ por el enlace no se cuenta como falla.
