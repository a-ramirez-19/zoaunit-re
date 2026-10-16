## Clave de caché idéntica a la que arma ajustar_simulacion() en cada qmd (requiere haber cargado funciones_<versión>.R).
prefijo <- switch(version_sim, v4 = "ughc", v4uw = "uwv4", v4kum = "kuv4", v3 = if (escala_log) "uwlg" else "uwhc")
clave <- function(tau, n, replica) {
  huella <- paste(c(beta_true(tau), gamma_true, delta_true, as.numeric(params$phi_true),
                    as.numeric(params$sigma_b_true)), collapse = "_")
  cod <- utf8ToInt(huella)
  huella <- sprintf("%07d", sum(cod * seq_along(cod)) %% 9999991)
  sprintf("%s_%s_tau%03d_n%04d_rep%04d_%d-%d-%d-%d-%d_g%s", prefijo, enlace_nombre, round(100 * tau), n, replica,
          config_mcmc$n_adapt, config_mcmc$n_burn, config_mcmc$n_iter, config_mcmc$n_thin, config_mcmc$n_chains, huella)
}
