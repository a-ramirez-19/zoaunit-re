## Corre un bloque de ajustes Monte Carlo de la simulación 2 v4 (Unit-Gamma).
## Variables de entorno:
##   ENLACE      logit | probit | cloglog
##   TAUS        p. ej. "0.50,0.75"
##   N_SUJETOS   p. ej. 500
##   R_MC        réplicas por tau (500)
##   BLOQUE, N_BLOQUES   este trabajo toma las tareas BLOQUE, BLOQUE + N_BLOQUES, ...
##   HORAS_MAX   no empieza ajustes nuevos pasado este tiempo (los que están en curso terminan)
##   CACHE_DIR   carpeta donde se guardan los .rds (mismo nombre que en el portátil)
##   HECHOS      archivo con claves ya calculadas (una por línea, sin .rds) que se saltan
t0 <- Sys.time()
source("funciones_v4.R")

env_num <- function(x, def) as.numeric(Sys.getenv(x, def))
taus_run <- leer_numeros(Sys.getenv("TAUS", "0.25,0.50,0.75"))
n_run <- as.integer(env_num("N_SUJETOS", 500))
R_mc <- as.integer(env_num("R_MC", 500))
bloque <- as.integer(env_num("BLOQUE", 1))
n_bloques <- as.integer(env_num("N_BLOQUES", 1))
horas_max <- env_num("HORAS_MAX", 5)
rep_ini <- as.integer(env_num("REP_INI", 1))
rep_fin <- as.integer(env_num("REP_FIN", R_mc))

## Misma clave que ajustar_simulacion(); se usa para saltar lo ya hecho sin generar datos.
clave <- function(tau, n, replica) {
  huella <- paste(c(beta_true(tau), gamma_true, delta_true, as.numeric(params$phi_true),
                    as.numeric(params$sigma_b_true)), collapse = "_")
  cod <- utf8ToInt(huella)
  huella <- sprintf("%07d", sum(cod * seq_along(cod)) %% 9999991)
  sprintf("%s_%s_tau%03d_n%04d_rep%04d_%d-%d-%d-%d-%d_g%s", "ughc", enlace_nombre, round(100 * tau), n, replica,
          config_mcmc$n_adapt, config_mcmc$n_burn, config_mcmc$n_iter, config_mcmc$n_thin, config_mcmc$n_chains, huella)
}

tareas <- expand.grid(rr = rep_ini:rep_fin, tt = taus_run, nn = n_run)
tareas$clave <- mapply(clave, tareas$tt, tareas$nn, tareas$rr)
hechos_archivo <- Sys.getenv("HECHOS", "hechos.txt")
hechos <- if (file.exists(hechos_archivo)) trimws(readLines(hechos_archivo)) else character(0)
tareas <- tareas[!tareas$clave %in% hechos, ]
tareas <- tareas[seq_len(nrow(tareas)) %% n_bloques == (bloque %% n_bloques), ]
cat(sprintf("Enlace %s | taus %s | n %d | bloque %d/%d | tareas: %d\n",
            enlace_nombre, paste(taus_run, collapse = ","), n_run, bloque, n_bloques, nrow(tareas)))
if (!nrow(tareas)) quit(save = "no")

n_workers <- as.integer(env_num("N_WORKERS", parallel::detectCores()))
cl <- parallel::makeCluster(n_workers)
invisible(parallel::clusterCall(cl, setwd, getwd()))
parallel::clusterExport(cl, c("tareas", "t0", "horas_max"))
invisible(parallel::clusterEvalQ(cl, source("funciones_v4.R")))
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
write.csv(res, sprintf("registro/bloque_%s_%03d.csv", enlace_nombre, bloque), row.names = FALSE)
print(table(sub(":.*", "", res$estado)))
cat(sprintf("Minutos por ajuste (mediana): %.1f | horas de reloj: %.2f\n",
            median(res$minutos, na.rm = TRUE), as.numeric(difftime(Sys.time(), t0, units = "hours"))))
