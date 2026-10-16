# Simulación 2, diseño v4, en GitHub Actions

Ajustes Monte Carlo de ZOAUNIT-RE para Unit-Weibull y Kumaraswamy con el mismo diseño de la
simulación v4 con Unit-Gamma: 3 enlaces (logit, probit, cloglog) × 3 cuartiles (0.25, 0.50, 0.75)
× 5 tamaños (n = 30, 50, 100, 200, 500) × 500 réplicas = 22.500 ajustes por familia.

- `funciones_v4uw.R`, `funciones_v4kum.R`, `funciones_v4.R`: funciones de `ZOAUNIT_RE_simulacion_v4_uw.qmd`,
  `ZOAUNIT_RE_simulacion_v4_kum.qmd` y `ZOAUNIT_RE_simulacion_v4.qmd`, copiadas sin cambios de sus chunks.
- `calibracion.rds`: media y DE de la edad y proporción de mujeres de los datos reales (los datos no se publican).
- `run_bloque.R`: corre un bloque de tareas; cada ajuste se guarda como `.rds` con la clave de caché del qmd
  (prefijos `uwv4`, `kuv4`, `ughc`), así que al copiarlos a `_cache_simulacion_asesor/` el render los reutiliza.
- `hechos.txt`: claves ya calculadas que se saltan (para relanzar lo que falte).

Se lanza desde **Actions → Simulacion → Run workflow**. Cada trabajo sube un artefacto
`cache-<versión>-<enlace>-<run>-bloque<k>`.

Verificado: con la misma semilla, los ajustes reproducen los del portátil (Windows, R 4.2.2) con diferencias
del orden de 1e-14; las log-densidades de JAGS de Unit-Weibull y Kumaraswamy coinciden con las fórmulas en
escala θ con los tres enlaces.
