"""Despachador de la simulación 2 (lo corre .github/workflows/despachador_sim2.yml cada 30 minutos).

GitHub cancela un trabajo que espera más de 24 h en cola sin conseguir runner. Por eso las corridas no se
lanzan todas de una vez: este script mira cuántos trabajos de simulacion_sim2.yml esperan runner y, si son
menos que el umbral del repo, lanza la siguiente unidad pendiente de la cola.

Cola: cola_sim2_<dueño del repo>.txt (por ejemplo cola_sim2_Sinu-Labs.txt). Una unidad por línea:
    familia enlace replicas etiqueta
    ug logit 101-150 sim2-ug-logit-r101-200
Líneas de configuración (opcionales):
    # umbral=120      lanzar mientras haya menos de 120 trabajos esperando runner
    # max_lanzar=2    como mucho 2 unidades por ejecución
    # bloques=128     trabajos por unidad
    # intentos=3      veces que se relanza una unidad cancelada o fallida

Cada corrida lleva el título "sim2 <etiqueta> <replicas>" (run-name del workflow), así se sabe qué unidades
ya se lanzaron. Una unidad se relanza (con la misma etiqueta) si su última corrida terminó cancelada o con
falla; el workflow se salta los ajustes ya publicados, así que solo completa lo que falta.
"""
import json
import math
import os
import sys
import urllib.request

REPO = os.environ["GITHUB_REPOSITORY"]
TOKEN = os.environ["GH_TOKEN"]
WORKFLOW = "simulacion_sim2.yml"
API = "https://api.github.com"
ACTIVOS = {"queued", "in_progress", "waiting", "requested", "pending"}
FALLIDOS = {"cancelled", "failure", "timed_out", "startup_failure"}


def api(ruta, metodo="GET", datos=None):
    req = urllib.request.Request(API + ruta, method=metodo,
                                 data=None if datos is None else json.dumps(datos).encode())
    req.add_header("Authorization", "Bearer " + TOKEN)
    req.add_header("Accept", "application/vnd.github+json")
    req.add_header("X-GitHub-Api-Version", "2022-11-28")
    with urllib.request.urlopen(req) as r:
        cuerpo = r.read()
        return json.loads(cuerpo) if cuerpo else None


def leer_cola():
    dueno = REPO.split("/")[0]
    archivo = os.environ.get("COLA") or f"cola_sim2_{dueno}.txt"
    conf = {"umbral": 60, "max_lanzar": 1, "bloques": 128, "intentos": 3}
    unidades = []
    with open(archivo, encoding="utf-8") as f:
        for linea in f:
            linea = linea.strip()
            if not linea:
                continue
            if linea.startswith("#"):
                for parte in linea.lstrip("#").split():
                    if "=" in parte:
                        k, v = parte.split("=", 1)
                        if k in conf:
                            conf[k] = int(v)
                continue
            familia, enlace, replicas, etiqueta = linea.split()
            unidades.append({"familia": familia, "enlace": enlace, "replicas": replicas, "etiqueta": etiqueta,
                             "titulo": f"sim2 {etiqueta} {replicas}"})
    return archivo, conf, unidades


def corridas():
    todas, pagina = [], 1
    while True:
        r = api(f"/repos/{REPO}/actions/workflows/{WORKFLOW}/runs?per_page=100&page={pagina}")
        todas += r["workflow_runs"]
        if len(r["workflow_runs"]) < 100:
            return todas
        pagina += 1


def trabajos_en_espera(activas, bloques):
    espera = 0
    for run in activas:
        sim, pagina = [], 1
        while True:
            r = api(f"/repos/{REPO}/actions/runs/{run['id']}/jobs?filter=latest&per_page=100&page={pagina}")
            sim += [j for j in r["jobs"] if j["name"].startswith("sim2")]
            if len(r["jobs"]) < 100:
                break
            pagina += 1
        if not sim:                      # la matriz todavía no se expandió
            espera += bloques
        else:
            espera += sum(1 for j in sim if j["status"] in ("queued", "waiting", "pending", "requested"))
    return espera


def main():
    archivo, conf, unidades = leer_cola()
    runs = corridas()
    activas = [r for r in runs if r["status"] in ACTIVOS]
    por_titulo = {}
    for r in sorted(runs, key=lambda r: r["created_at"]):
        por_titulo.setdefault(r["display_title"], []).append(r)

    nuevas, reintentos, hechas, en_curso = [], [], 0, 0
    for u in unidades:
        previas = por_titulo.get(u["titulo"], [])
        if not previas:
            nuevas.append(u)
            continue
        ultima = previas[-1]
        if ultima["status"] in ACTIVOS:
            en_curso += 1
        elif ultima["conclusion"] in FALLIDOS and len(previas) < conf["intentos"]:
            reintentos.append(u)
        else:
            hechas += 1
    pendientes = nuevas + reintentos

    espera = trabajos_en_espera(activas, conf["bloques"])
    print(f"Cola {archivo}: {len(unidades)} unidades | terminadas {hechas} | en curso {en_curso} | "
          f"nuevas {len(nuevas)} | por reintentar {len(reintentos)}")
    print(f"Corridas activas: {len(activas)} | trabajos esperando runner: {espera} | umbral {conf['umbral']}")

    if espera >= conf["umbral"] or not pendientes:
        print("No se lanza nada.")
        return
    n = min(conf["max_lanzar"], len(pendientes), math.ceil((conf["umbral"] - espera) / conf["bloques"]))
    for u in pendientes[:n]:
        entradas = {"familia": u["familia"], "enlaces": u["enlace"], "escenarios": "todos",
                    "taus": "0.25,0.50,0.75", "n_sujetos": "30,50,100,300,500", "replicas": u["replicas"],
                    "variante": "director", "n_bloques": str(conf["bloques"]), "horas_max": "5",
                    "etiqueta": u["etiqueta"]}
        api(f"/repos/{REPO}/actions/workflows/{WORKFLOW}/dispatches", "POST", {"ref": "main", "inputs": entradas})
        print(f"Lanzada: {u['titulo']}")


if __name__ == "__main__":
    try:
        main()
    except urllib.error.HTTPError as e:
        print(f"Error de la API: {e.code} {e.read().decode()[:500]}")
        sys.exit(1)
