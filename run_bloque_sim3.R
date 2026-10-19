## Corre un bloque de ajustes de la simulación 3 (un mismo proceso generador para todos los tau).
## Usa funciones_sim3.R, el mismo archivo que el Notebook ZOAUNIT_RE_simulacion_v5_comun.qmd, así que las
## claves de caché coinciden y los .rds se pueden copiar a _cache_sim3/ del portátil.
## Variables de entorno:
##   FAMILIA      ug (primera etapa) | uw | ku
##   ENLACES      p. ej. "logit" o "logit,probit,cloglog"
##   ESCENARIOS   "todos" o lista, p. ej. "iM_sA_vE"
##   TAUS         "0.25,0.50,0.75"
##   N_SUJETOS    "30,50,100,200,300,500"
##   REPLICAS     rango de réplicas "ini-fin" (o REP_INI y REP_FIN)
##   GUARDAR      liviano (por defecto: resumen + 250 extracciones por cadena) | completo
##   VARIANTE     director | previa_alt1 (análisis de sensibilidad: previas con escalas x 2)
##   BLOQUE, N_BLOQUES    este trabajo toma las tareas BLOQUE, BLOQUE + N_BLOQUES, ...
##   HORAS_MAX    no empieza ajustes nuevos pasado este tiempo (los que están en curso terminan)
##   CACHE_DIR    carpeta de los .rds (por defecto _cache_sim3)
##   HECHOS       archivo con claves ya publicadas (una por línea, con o sin .rds) que se saltan
t0 <- Sys.time()
env <- function(x, def) Sys.getenv(x, def)
familia <- env("FAMILIA", "ug")
enlaces <- trimws(strsplit(env("ENLACES", "logit"), ",")[[1]])
cache_dir <- env("CACHE_DIR", "_cache_sim3")
variante <- env("VARIANTE", "director")
params_sim3 <- list(cache_dir = cache_dir)
source("funciones_sim3.R")
guardar <- env("GUARDAR", "liviano")
P3 <- aplicar_variante(variante, familia); P3$guardar <- guardar
validar_configuracion(familia, incluir_perfiles = TRUE)   # se detiene si falta algún valor por definir

esc_arg <- env("ESCENARIOS", "todos")
escenarios <- if (esc_arg == "todos") escenarios_sim3()$id else trimws(strsplit(esc_arg, ",")[[1]])
taus <- leer_numeros(env("TAUS", "0.25,0.50,0.75"))
n_run <- sort(as.integer(leer_numeros(env("N_SUJETOS", "30,50,100,200,300,500"))), decreasing = TRUE)
reps <- env("REPLICAS", paste0(env("REP_INI", "1"), "-", env("REP_FIN", "50")))
rep_ini <- as.integer(sub("-.*", "", reps)); rep_fin <- as.integer(sub(".*-", "", reps))
bloque <- as.integer(env("BLOQUE", "1")); n_bloques <- as.integer(env("N_BLOQUES", "1"))
horas_max <- as.numeric(env("HORAS_MAX", "5"))

tareas <- expand.grid(rep = rep_ini:rep_fin, tau = taus, n = n_run, esc = escenarios, enlace = enlaces,
                      stringsAsFactors = FALSE)
tareas <- tareas[order(-tareas$n, tareas$esc, tareas$enlace, tareas$rep, tareas$tau), ]   # n grande primero
tareas$clave <- mapply(clave_sim3, familia, tareas$enlace, tareas$esc, tareas$tau, tareas$n, tareas$rep)
hechos_arch <- env("HECHOS", "hechos_sim3.txt")
hechos <- if (file.exists(hechos_arch)) sub("\\.rds$", "", trimws(readLines(hechos_arch))) else character(0)
tareas <- tareas[!tareas$clave %in% hechos, ]
tareas <- tareas[seq_len(nrow(tareas)) %% n_bloques == (bloque %% n_bloques), ]
cat(sprintf("Variante %s\n", variante))
cat(sprintf("Familia %s | enlaces %s | escenarios %s | n %s | réplicas %d-%d | bloque %d/%d | tareas: %d\n",
            familia, paste(enlaces, collapse = ","), paste(escenarios, collapse = ","), paste(n_run, collapse = ","),
            rep_ini, rep_fin, bloque, n_bloques, nrow(tareas)))
dir.create("registro", showWarnings = FALSE)
## Reproducibilidad (comentario 34): versiones, entorno y configuración de este bloque.
writeLines(c(capture.output(str(info_entorno())), "", capture.output(print(config_mcmc())), "",
             capture.output(sessionInfo())),
           sprintf("registro/entorno_sim3_%s_%s_%03d.txt", familia, variante, bloque))
if (!nrow(tareas)) quit(save = "no")

n_workers <- as.integer(env("N_WORKERS", parallel::detectCores()))
cl <- parallel::makeCluster(n_workers)
invisible(parallel::clusterCall(cl, setwd, getwd()))
parallel::clusterExport(cl, c("tareas", "t0", "horas_max", "familia", "cache_dir", "variante", "guardar"))
invisible(parallel::clusterEvalQ(cl, {
  params_sim3 <- list(cache_dir = cache_dir); source("funciones_sim3.R"); P3 <- aplicar_variante(variante, familia); P3$guardar <- guardar
}))
res <- tryCatch(
  parallel::parLapplyLB(cl, seq_len(nrow(tareas)), function(k) {
    g <- tareas[k, ]
    if (as.numeric(difftime(Sys.time(), t0, units = "hours")) > horas_max) {
      return(data.frame(clave = g$clave, estado = "sin_tiempo", minutos = NA_real_))
    }
    tryCatch({
      sim <- simular_dataset_comun(definir_escenario(g$esc, familia), g$n, g$rep)
      fit <- ajustar_sim3(sim, g$tau, g$enlace)
      estado <- if (!is.null(fit$error)) paste("error_jags:", fit$error) else fit$cache_status
      data.frame(clave = g$clave, estado = estado, minutos = fit$elapsed / 60)
    }, error = function(e) data.frame(clave = g$clave, estado = paste("error:", conditionMessage(e)), minutos = NA_real_))
  }, chunk.size = 1L),
  finally = parallel::stopCluster(cl)
)
res <- do.call(rbind, res)
write.csv(res, sprintf("registro/bloque_sim3_%s_%03d.csv", familia, bloque), row.names = FALSE)
print(table(sub(":.*", "", res$estado)))
cat(sprintf("Minutos por ajuste (mediana): %.2f | horas de reloj: %.2f\n",
            median(res$minutos, na.rm = TRUE), as.numeric(difftime(Sys.time(), t0, units = "hours"))))
