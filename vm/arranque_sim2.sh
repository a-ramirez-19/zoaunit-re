#!/bin/bash
# Arranque de una VM de la simulación 2 de ZOAUnit-RE en Google Cloud.
# Es el startup-script de cada VM: corre en cada encendido (también después de una interrupción de una
# VM spot) y retoma donde quedó, porque cada ajuste se guarda en su propio .rds y run_bloque_sim2.R
# reutiliza los que ya existen.
#
# Lee la configuración de los metadatos de la VM (los pone crear_vms_sim2.sh):
#   bucket        bucket de Cloud Storage donde se copian los ajustes (gs://<bucket>/cache/)
#   vm-index      número de esta VM (1..vm-total)
#   vm-total      número de VMs que se reparten el trabajo
#   commit        commit del repositorio con el que se corre (el mismo en todas las VMs)
#   escenarios    escenarios separados por coma (por defecto, los 9 del diseño reducido)
#   familias      orden de las familias (por defecto "ku uw ug": la más lenta primero)
#   enlaces, taus, n-sujetos, replicas   factores (por defecto, el diseño completo)
#   apagar-al-terminar   "si" apaga la VM cuando termina su parte (deja de cobrar CPU)
# Fuera de Google Cloud (prueba local) toma los mismos valores de variables de entorno en mayúsculas
# (BUCKET, VM_INDEX, ...), y BUCKET vacío desactiva la copia al bucket.

set -u
export PATH="$PATH:/snap/bin:/usr/local/bin"   # gcloud viene como snap en las imágenes de Ubuntu
DIR=${SIM2_DIR:-/opt/sim2}
REPO=${SIM2_REPO:-https://github.com/a-ramirez-19/zoaunit-re}
mkdir -p "$DIR"

# El startup-script de Google Cloud debe terminar pronto: el trabajo largo se lanza como servicio aparte.
if [ -z "${SIM2_DESACOPLADO:-}" ] && command -v systemd-run >/dev/null 2>&1 && [ -d /run/systemd/system ]; then
  cp "$0" "$DIR/arranque_sim2.sh"
  systemctl stop sim2 2>/dev/null; systemctl reset-failed sim2 2>/dev/null
  systemd-run --unit=sim2 --setenv=SIM2_DESACOPLADO=1 --setenv=SIM2_DIR="$DIR" --setenv=SIM2_REPO="$REPO" \
    /bin/bash "$DIR/arranque_sim2.sh"
  exit 0
fi
LOG=$DIR/sim2.log
exec >> "$LOG" 2>&1
echo "==== arranque $(date -u +%FT%TZ)"

meta() {  # meta <clave> <variable de entorno> <valor por defecto>
  local v
  v=$(curl -s -f -m 2 -H "Metadata-Flavor: Google" \
      "http://metadata.google.internal/computeMetadata/v1/instance/attributes/$1" 2>/dev/null) || v=""
  [ -z "$v" ] && v=${!2:-$3}
  echo "$v"
}
BUCKET=$(meta bucket BUCKET "")
VM_INDEX=$(meta vm-index VM_INDEX 1)
VM_TOTAL=$(meta vm-total VM_TOTAL 1)
COMMIT=$(meta commit COMMIT "")
ESCENARIOS=$(meta escenarios ESCENARIOS "i1008_s05_fM,i1008_s10_fM,i1008_s15_fM,i2030_s05_fM,i2030_s10_fM,i2030_s15_fM,i3020_s05_fM,i3020_s10_fM,i3020_s15_fM")
FAMILIAS=$(meta familias FAMILIAS "ku uw ug")
ENLACES=$(meta enlaces ENLACES "logit,probit,cloglog")
TAUS=$(meta taus TAUS "0.25,0.50,0.75")
N_SUJETOS=$(meta n-sujetos N_SUJETOS "30,50,100,300,500")
REPLICAS=$(meta replicas REPLICAS "1-500")
APAGAR=$(meta apagar-al-terminar APAGAR "no")
NOMBRE=$(hostname -s)
echo "vm $VM_INDEX/$VM_TOTAL | bucket ${BUCKET:-(ninguno)} | commit ${COMMIT:-(último)} | familias $FAMILIAS | réplicas $REPLICAS"
echo "escenarios $ESCENARIOS"

# 1. R y JAGS, los mismos paquetes de apt que en GitHub Actions (Ubuntu 24.04: R 4.3.3, JAGS 4.3.2).
if ! Rscript -e 'suppressMessages({library(rjags); library(posterior)})' >/dev/null 2>&1; then
  apt-get update -qq
  DEBIAN_FRONTEND=noninteractive apt-get install -y -qq --no-install-recommends \
    r-base-core jags r-cran-rjags r-cran-coda r-cran-posterior git curl
fi

# 2. Código, fijado en el commit indicado.
if [ ! -d "$DIR/repo/.git" ]; then git clone -q "$REPO" "$DIR/repo"; fi
cd "$DIR/repo"
git fetch -q origin || true
if [ -n "$COMMIT" ]; then git checkout -q "$COMMIT"; else git checkout -q origin/HEAD 2>/dev/null || true; fi
echo "código en $(git rev-parse --short HEAD)"
TRAB=$DIR/trabajo
mkdir -p "$TRAB/_cache_sim2" "$TRAB/registro"
cp funciones_sim2.R run_bloque_sim2.R calibracion.rds "$TRAB/"
cd "$TRAB"

# 3. Si el disco es nuevo (VM recreada), se recupera del bucket lo ya calculado.
if [ -n "$BUCKET" ] && [ -z "$(ls -A _cache_sim2 | head -1)" ]; then
  echo "caché local vacía: copiando la del bucket"
  gcloud storage rsync -q "gs://$BUCKET/cache" _cache_sim2 || true
fi

# 4. Ajustes ya hechos en otra parte: listas de los releases sim2-* del repositorio público y, si hay,
#    listas extra que se suban a gs://<bucket>/hechos/ (por ejemplo, las de Sinu-Labs).
if [ ! -s hechos_sim2.txt ]; then
  : > hechos_sim2.txt
  curl -s -f "https://api.github.com/repos/${REPO#https://github.com/}/releases?per_page=100" \
    | grep -o '"browser_download_url": *"[^"]*/sim2-[^"]*/lista_[^"]*\.txt"' | sed 's/.*"\(https[^"]*\)"/\1/' \
    | while read -r u; do curl -s -f -L "$u" >> hechos_sim2.txt || true; done
  if [ -n "$BUCKET" ]; then gcloud storage cat "gs://$BUCKET/hechos/*.txt" >> hechos_sim2.txt 2>/dev/null || true; fi
  echo "claves ya hechas en otra parte: $(wc -l < hechos_sim2.txt)"
fi

# 5. Copia periódica al bucket de los ajustes nuevos (cada 10 minutos), en segundo plano.
copiar() {
  [ -z "$BUCKET" ] && return 0
  touch .marca_nueva
  if [ -f .marca ]; then find _cache_sim2 -name '*.rds' -newer .marca > .nuevos; else find _cache_sim2 -name '*.rds' > .nuevos; fi
  if [ -s .nuevos ]; then gcloud storage cp -q -I "gs://$BUCKET/cache/" < .nuevos && mv .marca_nueva .marca; else mv .marca_nueva .marca; fi
  gcloud storage rsync -q -r registro "gs://$BUCKET/registro/$NOMBRE" >/dev/null 2>&1 || true
  gcloud storage cp -q "$LOG" "gs://$BUCKET/registro/$NOMBRE/sim2.log" >/dev/null 2>&1 || true
}

# 6. Trabajo: K procesos de R por VM (R 4.3 admite ~125 conexiones, así que cada proceso usa como mucho
#    100 núcleos). Las tareas se reparten por turnos entre los VM_TOTAL * K bloques; dentro de cada
#    bloque, n grande primero. Las familias van en el orden de FAMILIAS.
NUCLEOS=$(nproc)
K=$(( (NUCLEOS + 99) / 100 ))
W=$(( NUCLEOS / K ))
NB=$(( VM_TOTAL * K ))
echo "núcleos $NUCLEOS -> $K procesos de $W trabajadores; $NB bloques en total"

if [ -f TERMINADO ]; then
  echo "esta VM ya había terminado su parte"
else
  for k in $(seq 1 "$K"); do
    B=$(( (VM_INDEX - 1) * K + k ))
    (
      for f in $FAMILIAS; do
        echo "bloque $B/$NB: familia $f, inicio $(date -u +%FT%TZ)"
        FAMILIA=$f ENLACES=$ENLACES ESCENARIOS=$ESCENARIOS TAUS=$TAUS N_SUJETOS=$N_SUJETOS REPLICAS=$REPLICAS \
          VARIANTE=director BLOQUE=$B N_BLOQUES=$NB N_WORKERS=$W HORAS_MAX=10000 \
          CACHE_DIR=_cache_sim2 GUARDAR=liviano HECHOS=hechos_sim2.txt \
          Rscript run_bloque_sim2.R > "registro/bloque_${B}_$f.out" 2>&1
        echo "bloque $B/$NB: familia $f, fin $(date -u +%FT%TZ) (estado $?)"
      done
    ) &
  done
  ( while true; do sleep 600; copiar; done ) &
  COPIADOR=$!
  wait $(jobs -p | grep -v "^$COPIADOR$")
  kill "$COPIADOR" 2>/dev/null
  date -u +%FT%TZ > TERMINADO
fi

copiar
[ -n "$BUCKET" ] && gcloud storage cp -q TERMINADO "gs://$BUCKET/estado/$NOMBRE.TERMINADO" || true
echo "==== fin $(date -u +%FT%TZ): $(find _cache_sim2 -name '*.rds' | wc -l) ajustes en la caché local"
if [ "$APAGAR" = "si" ]; then shutdown -h now; fi
