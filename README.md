# Simulación 2 v4 (Unit-Gamma), n = 500, en GitHub Actions

Corre en paralelo los ajustes Monte Carlo que faltan de `ZOAUNIT_RE_simulacion_v4.qmd`.

- `funciones_v4.R`: las funciones del qmd, copiadas sin cambios (generador, modelo JAGS, ajuste con caché).
- `calibracion.rds`: media y DE de la edad y proporción de mujeres, calculadas desde los datos reales
  (los datos no se publican).
- `run_bloque.R`: corre un bloque de réplicas; cada ajuste se guarda como `.rds` con la misma clave
  que usa el qmd, de modo que al copiarlos a `_cache_simulacion_asesor/` el render los reutiliza.
- `hechos.txt`: claves que ya existen en el portátil y que se saltan.

Se lanza desde la pestaña **Actions → Simulacion n500 → Run workflow**. Cada bloque sube un artefacto
`cache-<enlace>-n<n>-bloque<k>` con sus `.rds`.

Verificado: con la misma semilla, los ajustes reproducen los del portátil (Windows, R 4.2.2) con
diferencias del orden de 1e-14.
