#!/bin/bash
# Estado de la simulación 2 en Google Cloud. Se corre en Cloud Shell, desde la carpeta vm/ del repo.
#   bash estado_sim2.sh              # VMs, ajustes en el bucket y últimas líneas de cada registro
#   REINICIAR=1 bash estado_sim2.sh  # además enciende las VMs spot que Google detuvo y no han terminado
# Variables: PROYECTO, ZONA (us-central1-a), BUCKET (sim2-zoaunit-<proyecto>), PREFIJO (sim2-)
set -uo pipefail
PROYECTO=${PROYECTO:-$(gcloud config get-value project 2>/dev/null)}
ZONA=${ZONA:-us-central1-a}
BUCKET=${BUCKET:-sim2-zoaunit-$PROYECTO}
PREFIJO=${PREFIJO:-sim2-}

echo "== VMs ($(date -u +%FT%TZ))"
gcloud compute instances list --project "$PROYECTO" --filter "name~^$PREFIJO" \
  --format "table(name,zone.basename(),machineType.basename(),status,scheduling.provisioningModel)"

echo; echo "== Ajustes en gs://$BUCKET/cache por familia"
gcloud storage ls "gs://$BUCKET/cache/" 2>/dev/null | sed 's#.*/##' | cut -d_ -f1 | sort | uniq -c
echo "Total del diseño reducido: 607 500 (67 500 por familia y enlace). Hay que sumar lo ya publicado en los releases."

echo; echo "== VMs que terminaron su parte"
gcloud storage ls "gs://$BUCKET/estado/" 2>/dev/null || echo "  ninguna todavía"

echo; echo "== Últimas líneas del registro de cada VM"
for f in $(gcloud storage ls "gs://$BUCKET/registro/*/sim2.log" 2>/dev/null); do
  echo "-- $f"; gcloud storage cat "$f" | tail -4
done

if [ -n "${REINICIAR:-}" ]; then
  echo; echo "== Reinicio de VMs detenidas"
  for v in $(gcloud compute instances list --project "$PROYECTO" --filter "name~^$PREFIJO AND status=TERMINATED" --format "value(name)"); do
    if gcloud storage ls "gs://$BUCKET/estado/$v.TERMINADO" >/dev/null 2>&1; then
      echo "  $v ya terminó: no se reinicia"
    else
      echo "  encendiendo $v"; gcloud compute instances start "$v" --zone "$ZONA" --project "$PROYECTO" -q
    fi
  done
fi
