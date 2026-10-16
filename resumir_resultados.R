## Resume los ajustes de una versión y un enlace igual que el chunk `ejecutar-montecarlo` del qmd:
## resultados por réplica y métricas (RB, MSE, RMSE, CP, ...), con los mismos nombres de archivo.
## Solo lee la caché; los ajustes que falten se listan en faltantes_<enlace>.txt.
## Variables de entorno: VERSION, ENLACE, CACHE_DIR, RESULTADOS_DIR
version_sim <- Sys.getenv("VERSION", "v4uw")
archivo_funciones <- sprintf("funciones_%s.R", version_sim)
source(archivo_funciones)
source("claves.R")
tamanos_mc <- c(30L, 50L, 100L, 200L, 500L)
R_mc <- 500L
## Mismas carpetas que el YAML de cada qmd.
resultados_dir <- Sys.getenv("RESULTADOS_DIR", switch(version_sim, v4uw = "resultados_sim2_v4_uw", v4kum = "resultados_sim2_v4_kum",
                                                      sprintf("resultados_sim2_%s", version_sim)))

grid_mc <- expand.grid(rr = seq_len(R_mc), tt = taus, nn = sort(tamanos_mc, decreasing = TRUE))
grid_mc$ruta <- file.path(params$cache_dir, paste0(mapply(clave, grid_mc$tt, grid_mc$nn, grid_mc$rr), ".rds"))
falta <- !file.exists(grid_mc$ruta)
dir.create(resultados_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(basename(grid_mc$ruta[falta]), file.path(resultados_dir, sprintf("faltantes_%s.txt", enlace_nombre)))
cat(sprintf("%s %s: %d ajustes en caché, %d faltan\n", version_sim, enlace_nombre, sum(!falta), sum(falta)))
grid_ok <- grid_mc[!falta, ]

cl <- parallel::makeCluster(parallel::detectCores())
invisible(parallel::clusterCall(cl, setwd, getwd()))
parallel::clusterExport(cl, c("grid_ok", "archivo_funciones"))
invisible(parallel::clusterEvalQ(cl, source(archivo_funciones)))
resultados_mc <- dplyr::bind_rows(parallel::parLapply(cl, seq_len(nrow(grid_ok)), function(k) {
  fit <- readRDS(grid_ok$ruta[k]); fit$cache_status <- "reutilizado"; resumir_fit(fit)
}))
parallel::stopCluster(cl)

## Mismo resumen que el qmd.
metricas_mc <- resultados_mc |>
  group_by(n, Tau, Parametro) |>
  summarise(
    R = n(),
    RB = mean(Error_rel, na.rm = TRUE),
    MSE = mean(Error2),
    RMSE = sqrt(mean(Error2)),
    CP = mean(Cubre),
    Sesgo_abs = mean(Error),
    MAE = mean(Error_abs),
    Longitud_ICr = mean(Longitud_ICr),
    Rhat_max = max(Rhat, na.rm = TRUE),
    ESS_mediana = median(ESS, na.rm = TRUE),
    Tiempo_medio_min = mean(Tiempo_min),
    .groups = "drop"
  )
base_res <- file.path(resultados_dir, sprintf("sim2_%s_nmax%d", enlace_nombre, max(tamanos_mc)))
saveRDS(resultados_mc, paste0(base_res, "_resultados_por_replica.rds"))
saveRDS(metricas_mc, paste0(base_res, "_metricas.rds"))
write.csv(resultados_mc, paste0(base_res, "_resultados_por_replica.csv"), row.names = FALSE)
write.csv(metricas_mc, paste0(base_res, "_metricas.csv"), row.names = FALSE)
print(as.data.frame(metricas_mc |> filter(n == max(n)) |> select(Tau, Parametro, RB, RMSE, CP)), digits = 3)
