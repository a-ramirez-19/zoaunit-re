# Simulación 2 en máquinas virtuales de Google Cloud

Desde el 30/10/2026 la simulación 2 ya no corre en GitHub Actions. Estos scripts la corren en VMs spot de Google Cloud, con el diseño reducido de 9 escenarios (inflación × σ_b, φ en el nivel `fM`).

- `crear_vms_sim2.sh`: crea el bucket y las VMs. Se corre en Cloud Shell. `MODO=prueba` crea 1 VM pequeña (~40 min) para probar el flujo y medir tiempos.
- `arranque_sim2.sh`: startup-script de cada VM. Instala R y JAGS de apt (Ubuntu 24.04: R 4.3.3, JAGS 4.3.2, como en Actions), clona este repo en el commit fijado, salta los ajustes ya publicados en los releases `sim2-*` y en `gs://<bucket>/hechos/`, corre `run_bloque_sim2.R` en todos los núcleos, copia cada 10 minutos los `.rds` nuevos a `gs://<bucket>/cache/` y apaga la VM al terminar. Si Google detiene la VM, al encenderla retoma donde quedó.
- `estado_sim2.sh`: VMs, ajustes en el bucket y registros; `REINICIAR=1` enciende las VMs detenidas.

El procedimiento completo está en el runbook «Runbook: corridas de la simulación 2 en máquinas virtuales».
