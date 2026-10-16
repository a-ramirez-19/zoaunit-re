# Simulación 2 (v3 Unit-Weibull y v4 Unit-Gamma) en GitHub Actions

Corre en paralelo los ajustes Monte Carlo que faltan de `ZOAUNIT_RE_simulacion_v3.qmd` y `ZOAUNIT_RE_simulacion_v4.qmd`.

- `funciones_v3.R`, `funciones_v4.R`: las funciones de cada qmd, copiadas sin cambios (generador, modelo JAGS, ajuste con caché).
- `calibracion.rds`: media y DE de la edad y proporción de mujeres, calculadas desde los datos reales
  (los datos no se publican).
- `run_bloque.R`: corre un bloque de réplicas; cada ajuste se guarda como `.rds` con la misma clave
  que usa el qmd, de modo que al copiarlos a `_cache_simulacion_asesor/` el render los reutiliza.
- `hechos.txt`: claves que ya existen en el portátil y que se saltan.

Se lanza desde **Actions → Simulacion → Run workflow**. Cada bloque sube un artefacto
`cache-<versión>-<enlace>-<run>-bloque<k>` con sus `.rds`.

Verificado: con la misma semilla, los ajustes reproducen los del portátil (Windows, R 4.2.2) con
diferencias del orden de 1e-14.
