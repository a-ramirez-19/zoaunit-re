## Corre un bloque de ajustes Monte Carlo de la simulación 2.
## Variables de entorno:
##   VERSION     v4 (Unit-Gamma) | v4uw (Unit-Weibull) | v4kum (Kumaraswamy) | v3 (Unit-Weibull, v3)
##   ENLACE      logit | probit | cloglog
##   TAUS        p. ej. "0.50,0.75"
##   N_SUJETOS   p. ej. "500" o "30,50,100,200,500"
##   R_MC        réplicas por tau (500)
##   BLOQUE, N_BLOQUES   este trabajo toma las tareas BLOQUE, BLOQUE + N_BLOQUES, ...
##   HORAS_MAX   no empieza ajustes nuevos pasado este tiempo (los que están en curso terminan)
##   CACHE_DIR   carpeta donde se guardan los .rds (mismo nombre que en el portátil)
##   HECHOS      archivo con claves ya calculadas (una por línea, sin .rds) que se saltan
t0 <- Sys.time()
version_sim <- Sys.getenv("VERSION", "v4")
archivo_funciones <- sprintf("funciones_%s.R", version_sim)
source(archivo_funciones)

env_num <- function(x, def) as.numeric(Sys.getenv(x, def))
taus_run <- leer_numeros(Sys.getenv("TAUS", "0.25,0.50,0.75"))
n_run <- sort(as.integer(leer_numeros(Sys.getenv("N_SUJETOS", "500"))), decreasing = TRUE)
R_mc <- as.integer(env_num("R_MC", 500))
bloque <- as.integer(env_num("BLOQUE", 1))
n_bloques <- as.integer(env_num("N_BLOQUES", 1))
horas_max <- env_num("HORAS_MAX", 5)
rep_ini <- as.integer(env_num("REP_INI", 1))
rep_fin <- as.integer(env_num("REP_FIN", R_mc))

source("claves.R")

tareas <- expand.grid(rr = rep_ini:rep_fin, tt = taus_run, nn = n_run)  # n grande primero
tareas$clave <- mapply(clave, tareas$tt, tareas$nn, tareas$rr)
hechos_archivo <- Sys.getenv("HECHOS", "hechos.txt")
hechos <- if (file.exists(hechos_archivo)) trimws(readLines(hechos_archivo)) else character(0)
tareas <- tareas[!tareas$clave %in% hechos, ]
tareas <- tareas[seq_len(nrow(tareas)) %% n_bloques == (bloque %% n_bloques), ]
cat(sprintf("Versión %s | enlace %s | taus %s | n %s | bloque %d/%d | tareas: %d\n",
            version_sim, enlace_nombre, paste(taus_run, collapse = ","), paste(n_run, collapse = ","), bloque, n_bloques, nrow(tareas)))
if (!nrow(tareas)) quit(save = "no")

n_workers <- as.integer(env_num("N_WORKERS", parallel::detectCores()))
cl <- parallel::makeCluster(n_workers)
invisible(parallel::clusterCall(cl, setwd, getwd()))
parallel::clusterExport(cl, c("tareas", "t0", "horas_max"))
parallel::clusterExport(cl, "archivo_funciones")
invisible(parallel::clusterEvalQ(cl, source(archivo_funciones)))
res <- tryCatch(
  parallel::parLapplyLB(cl, seq_len(nrow(tareas)), function(k) {
    g <- tareas[k, ]
    if (as.numeric(difftime(Sys.time(), t0, units = "hours")) > horas_max) {
      return(data.frame(clave = g$clave, estado = "sin_tiempo", minutos = NA_real_))
    }
    r <- tryCatch({
      sim <- simular_dataset_tau(g$tt, n_sujetos = g$nn, replica = g$rr)
      fit <- ajustar_simulacion(sim)
      data.frame(clave = g$clave, estado = fit$cache_status, minutos = fit$elapsed / 60)
    }, error = function(e) data.frame(clave = g$clave, estado = paste("error:", conditionMessage(e)), minutos = NA_real_))
    r
  }, chunk.size = 1L),
  finally = parallel::stopCluster(cl)
)
res <- do.call(rbind, res)
dir.create("registro", showWarnings = FALSE)
write.csv(res, sprintf("registro/bloque_%s_%s_n%d_%03d.csv", version_sim, enlace_nombre, max(n_run), bloque), row.names = FALSE)
print(table(sub(":.*", "", res$estado)))
cat(sprintf("Minutos por ajuste (mediana): %.1f | horas de reloj: %.2f\n",
            median(res$minutos, na.rm = TRUE), as.numeric(difftime(Sys.time(), t0, units = "hours"))))
