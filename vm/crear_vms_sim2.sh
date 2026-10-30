#!/bin/bash
# Crea las VMs de la simulación 2 en Google Cloud. Se corre en Cloud Shell, desde la carpeta vm/ del repo:
#   git clone https://github.com/a-ramirez-19/zoaunit-re && cd zoaunit-re/vm
#   MODO=prueba bash crear_vms_sim2.sh      # 1 VM pequeña, ~40 min: prueba el flujo y mide tiempos
#   bash crear_vms_sim2.sh                  # diseño reducido completo (9 escenarios, R = 500)
# Variables (todas opcionales):
#   MODO       prueba | completo (por defecto completo)
#   PROYECTO   proyecto de Google Cloud (por defecto el activo en gcloud)
#   ZONA       por defecto us-central1-a
#   TIPO       tipo de máquina (por defecto c3d-highcpu-180; prueba: c3d-highcpu-8)
#   N_VMS      número de VMs que se reparten el trabajo (por defecto 2; prueba: 1)
#   PREFIJO    prefijo del nombre de las VMs (por defecto sim2-vm; prueba: sim2-prueba)
#   BUCKET     bucket de resultados (por defecto sim2-zoaunit-<proyecto>)
#   REPLICAS, ESCENARIOS, FAMILIAS   para cambiar el alcance
#   SI=1       no pide confirmación
# Las VMs son spot (mucho más baratas; Google puede detenerlas). Si se detienen, se reinician con
# estado_sim2.sh y retoman donde quedaron. Cada VM se apaga sola al terminar su parte.
set -euo pipefail
cd "$(dirname "$0")"

MODO=${MODO:-completo}
PROYECTO=${PROYECTO:-$(gcloud config get-value project 2>/dev/null)}
ZONA=${ZONA:-us-central1-a}
REGION=${ZONA%-*}
BUCKET=${BUCKET:-sim2-zoaunit-$PROYECTO}
COMMIT=$(git rev-parse HEAD)
if [ "$MODO" = "prueba" ]; then
  TIPO=${TIPO:-c3d-highcpu-8}; N_VMS=${N_VMS:-1}; PREFIJO=${PREFIJO:-sim2-prueba}
  ESCENARIOS=${ESCENARIOS:-i1008_s10_fM}; FAMILIAS=${FAMILIAS:-ku ug}; REPLICAS=${REPLICAS:-1-1}
else
  TIPO=${TIPO:-c3d-highcpu-180}; N_VMS=${N_VMS:-2}; PREFIJO=${PREFIJO:-sim2-vm}
  ESCENARIOS=${ESCENARIOS:-i1008_s05_fM,i1008_s10_fM,i1008_s15_fM,i2030_s05_fM,i2030_s10_fM,i2030_s15_fM,i3020_s05_fM,i3020_s10_fM,i3020_s15_fM}
  FAMILIAS=${FAMILIAS:-ku uw ug}; REPLICAS=${REPLICAS:-1-500}
fi
[ -z "$PROYECTO" ] && { echo "No hay proyecto activo: gcloud config set project <ID>"; exit 1; }

echo "Proyecto $PROYECTO | zona $ZONA | $N_VMS x $TIPO (spot) | bucket gs://$BUCKET | commit ${COMMIT:0:7}"
echo "Modo $MODO | familias $FAMILIAS | réplicas $REPLICAS"
echo "Escenarios: $ESCENARIOS"

gcloud services enable compute.googleapis.com storage.googleapis.com --project "$PROYECTO" -q

echo; echo "Cuotas de la región $REGION (límite / en uso):"
gcloud compute regions describe "$REGION" --project "$PROYECTO" --format=json \
  | python3 -c "import json,sys; q=json.load(sys.stdin)['quotas']; [print(f\"  {x['metric']}: {x['limit']:.0f} / {x['usage']:.0f}\") for x in q if x['metric'] in ('CPUS','C3D_CPUS','PREEMPTIBLE_CPUS','SSD_TOTAL_GB','DISKS_TOTAL_GB')]"
NUC=$(echo "$TIPO" | grep -o '[0-9]*$')
echo "Se necesitan $(( NUC * N_VMS )) vCPU de C3D (y de CPUS o PREEMPTIBLE_CPUS)."

if [ -z "${SI:-}" ]; then read -r -p "¿Crear las VMs? (s/n) " r; [[ "$r" =~ ^[sS] ]] || { echo "No se creó nada."; exit 0; }; fi

gcloud storage buckets describe "gs://$BUCKET" >/dev/null 2>&1 || \
  gcloud storage buckets create "gs://$BUCKET" --project "$PROYECTO" --location "$REGION" --uniform-bucket-level-access
# La cuenta de servicio por defecto de las VMs no siempre tiene permisos sobre Cloud Storage: se le dan
# sobre este bucket (si ya los tiene, no cambia nada).
SA="$(gcloud projects describe "$PROYECTO" --format 'value(projectNumber)')-compute@developer.gserviceaccount.com"
gcloud storage buckets add-iam-policy-binding "gs://$BUCKET" --member "serviceAccount:$SA" \
  --role roles/storage.objectAdmin --format 'value(etag)' >/dev/null

T=$(mktemp -d)
printf '%s' "$ESCENARIOS" > "$T/escenarios"; printf '%s' "$FAMILIAS" > "$T/familias"
printf '%s' "${ENLACES:-logit,probit,cloglog}" > "$T/enlaces"; printf '%s' "${TAUS:-0.25,0.50,0.75}" > "$T/taus"
printf '%s' "${N_SUJETOS:-30,50,100,300,500}" > "$T/n-sujetos"

for i in $(seq 1 "$N_VMS"); do
  NOMBRE=$PREFIJO$i
  echo "Creando $NOMBRE ..."
  gcloud compute instances create "$NOMBRE" --project "$PROYECTO" --zone "$ZONA" \
    --machine-type "$TIPO" --provisioning-model SPOT --instance-termination-action STOP \
    --image-family ubuntu-2404-lts-amd64 --image-project ubuntu-os-cloud \
    --boot-disk-size 50GB --boot-disk-type pd-balanced \
    --scopes storage-rw --labels proyecto=zoaunit-re,simulacion=sim2 \
    --metadata "bucket=$BUCKET,vm-index=$i,vm-total=$N_VMS,commit=$COMMIT,replicas=$REPLICAS,apagar-al-terminar=si" \
    --metadata-from-file "startup-script=arranque_sim2.sh,escenarios=$T/escenarios,familias=$T/familias,enlaces=$T/enlaces,taus=$T/taus,n-sujetos=$T/n-sujetos"
done
rm -rf "$T"
echo
echo "Listo. Cada VM instala R y JAGS (~3 min) y empieza. Para seguirlas:  bash estado_sim2.sh"
echo "Registro de cada VM: gs://$BUCKET/registro/<nombre de la VM>/sim2.log (se copia cada 10 minutos)."
