## Simulación 2 v4: Kumaraswamy (diseño de la v4, escala log(-log) con los tres enlaces).
## Funciones de ZOAUNIT_RE_simulacion_v4_kum.qmd, copiadas sin cambios de los chunks
## setup, valores-beta, valores-masas, funcion-generadora, modelo-jags y funciones-ajuste.
## Los parámetros son los del YAML del qmd; solo el enlace y la carpeta de caché vienen del entorno.
## La calibración (media y DE de la edad, proporción de mujeres) se lee de calibracion.rds,
## calculada en el portátil desde PropData_workDiana.csv (el CSV no se publica).

suppressPackageStartupMessages({
  library(dplyr)
  library(coda)
  library(rjags)
})
load.module("glm", quiet = TRUE)

params <- list(
  taus = "0.25,0.50,0.75",
  m_repeticiones = 4,
  semilla = 20260927,
  beta_q10 = "0.50,0.60,0.40",
  beta_q25 = "0.50,0.60,0.40",
  beta_q50 = "0.50,0.60,0.40",
  beta_q75 = "0.50,0.60,0.40",
  beta_q90 = "0.50,0.60,0.40",
  p0_ref = 0.10,
  p1_ref = 0.08,
  gamma_z = 0.35,
  delta_z = -0.35,
  phi_true = 3.8,
  sigma_b_true = 1.0,
  enlace = Sys.getenv("ENLACE", "logit"),
  n_adapt = 750,
  n_burn = 1200,
  n_iter = 3000,
  n_thin = 3,
  n_chains = 2,
  icr = 0.95,
  cache_dir = Sys.getenv("CACHE_DIR", "_cache_simulacion_asesor"),
  force_refit = FALSE
)

calibracion <- readRDS(Sys.getenv("CALIBRACION", "calibracion.rds"))

## ---- setup (helpers) ----

dir.create(params$cache_dir, recursive = TRUE, showWarnings = FALSE)

leer_numeros <- function(x, n_esperado = NULL) {
  z <- as.numeric(trimws(strsplit(as.character(x), ",", fixed = TRUE)[[1]]))
  if (any(!is.finite(z))) stop("Se encontraron valores no numéricos en: ", x)
  if (!is.null(n_esperado) && length(z) != n_esperado) {
    stop("Se esperaban ", n_esperado, " valores y se recibieron ", length(z), ".")
  }
  z
}

clip01 <- function(x, eps = 1e-10) pmin(pmax(x, eps), 1 - eps)
inv_logit <- function(x) plogis(x)

## Enlace del predictor cuantilico. El asesor pidio barrer las tres funciones
## despues de fijar la familia interior; se elige con params$enlace.
enlaces_disponibles <- list(
  logit   = list(inv = function(x) plogis(x),               jags = "ilogit(eta[m])"),
  probit  = list(inv = function(x) pnorm(x),                jags = "phi(eta[m])"),
  cloglog = list(inv = function(x) 1 - exp(-exp(x)),        jags = "icloglog(eta[m])")
)
enlace_nombre <- tolower(trimws(as.character(params$enlace)))
if (!enlace_nombre %in% names(enlaces_disponibles)) {
  stop("params$enlace debe ser 'logit', 'probit' o 'cloglog'; se recibio: ", enlace_nombre)
}
inv_enlace <- enlaces_disponibles[[enlace_nombre]]$inv
## v4: los tres enlaces se trabajan en escala log(-log theta).
escala_log <- TRUE
jags_inv_enlace <- enlaces_disponibles[[enlace_nombre]]$jags

softmax_masas <- function(eta0, eta1) {
  mm <- pmax(eta0, eta1, 0)
  den <- exp(-mm) + exp(eta0 - mm) + exp(eta1 - mm)
  list(
    p0 = exp(eta0 - mm) / den,
    p1 = exp(eta1 - mm) / den,
    pc = exp(-mm) / den
  )
}

## Kumaraswamy: F(y) = 1 - (1 - y^a)^phi, con a tal que Q(tau) = theta:
## a = log(qb) / log theta, qb = 1 - (1 - tau)^(1/phi). En escala llt = log(-log theta):
## log a = log(-log qb) - llt y log(-log y) = llt - log(-log qb) + log(-log(1 - (1 - u)^(1/phi))).
lqb_kum <- function(phi, tau) log(-log(-expm1(log1p(-tau) / phi)))
q_kum_llt <- function(u, llt, phi, tau) {
  exp(-exp(llt - lqb_kum(phi, tau) + log(-log1p(-exp(log1p(-u) / phi)))))
}

## v4: Unit-Gamma. Y = exp(-T), T ~ Gamma(phi, tasa b), b = q_Gamma(1 - tau; phi, 1) / (-log theta).
## Todo se escribe con llt = log(-log theta): -log theta = exp(llt).
q_ug_llt <- function(u, llt, phi, tau) {
  b <- qgamma(1 - tau, phi, 1) * exp(-llt)
  exp(-qgamma(1 - u, phi, rate = b))
}
ld_ug_llt <- function(ly, llt, phi, tau) {
  lb <- log(qgamma(1 - tau, phi, 1)) - llt
  phi * lb - lgamma(phi) + (phi - 1) * log(-ly) + (exp(lb) - 1) * ly
}

## Corrección de cloglog (07/10/2026). Con theta = 1 - exp(-exp(eta)), 1 - theta decrece como doble
## exponencial y theta llega a valer 1 en doble precisión. Se trabaja con l_theta = log(-log theta),
## que es exacto para cualquier eta: log(-log1p(-exp(-x))) con x = exp(eta), y -x cuando exp(-x) se anula.
lth_cloglog <- function(eta) {
  x <- exp(eta)
  ifelse(x < 700, log(-log1p(-exp(-pmin(x, 700)))), -x)
}

## v4: log(-log theta) estable con cada enlace.
lth_enlace <- function(eta, enlace = enlace_nombre) {
  switch(
    enlace,
    logit = log(log1p(exp(-pmax(eta, -700)))),
    probit = log(-pnorm(eta, log.p = TRUE)),
    cloglog = {
      x <- exp(pmin(eta, 700))
      ## x pequeño: -log(theta) = -log(-expm1(-x)) es exacto; x grande: la fórmula de lth_cloglog().
      suppressWarnings(ifelse(x > 1, lth_cloglog(eta), log(-log(-expm1(-x)))))
    },
    stop("Enlace no reconocido: ", enlace)
  )
}
## ---- valores-beta ----
taus <- leer_numeros(params$taus)
if (!all(taus > 0 & taus < 1)) stop("Todos los valores de tau deben estar en (0,1).")

beta_por_tau <- list(
  `0.10` = leer_numeros(params$beta_q10, 3),
  `0.25` = leer_numeros(params$beta_q25, 3),
  `0.50` = leer_numeros(params$beta_q50, 3),
  `0.75` = leer_numeros(params$beta_q75, 3),
  `0.90` = leer_numeros(params$beta_q90, 3)
)

beta_true <- function(tau) {
  key <- sprintf("%.2f", tau)
  if (!key %in% names(beta_por_tau)) {
    stop("No hay beta verdadero definido para tau = ", key,
         ". Añada el vector correspondiente antes de ejecutar.")
  }
  setNames(beta_por_tau[[key]], c("Intercepto", "Sexo", "Edad_s"))
}


## ---- valores-masas ----
p0_ref <- as.numeric(params$p0_ref)
p1_ref <- as.numeric(params$p1_ref)
pc_ref <- 1 - p0_ref - p1_ref
if (pc_ref <= 0) stop("Se requiere p0_ref + p1_ref < 1.")

gamma_true <- c(
  Intercepto = log(p0_ref / pc_ref),
  z = as.numeric(params$gamma_z)
)
delta_true <- c(
  Intercepto = log(p1_ref / pc_ref),
  z = as.numeric(params$delta_z)
)


## ---- funcion-generadora ----
simular_dataset_tau <- function(
    tau,
    n_sujetos,
    replica = 1L,
    m = as.integer(params$m_repeticiones),
    seed_base = as.integer(params$semilla)
) {
  tau <- as.numeric(tau)
  n_sujetos <- as.integer(n_sujetos)
  replica <- as.integer(replica)
  m <- as.integer(m)
  if (n_sujetos < 5) stop("n_sujetos debe ser al menos 5.")
  if (m < 2) stop("Se requieren al menos dos observaciones por sujeto para el efecto aleatorio.")

  seed_tau <- seed_base + 100000L * replica + as.integer(round(10000 * tau)) + n_sujetos
  set.seed(seed_tau)

  ## Covariables de sujeto sugeridas por el asesor.
  sexo <- rbinom(n_sujetos, size = 1, prob = calibracion$p_mujer)
  edad <- rnorm(n_sujetos, mean = calibracion$edad_media, sd = calibracion$edad_sd)
  edad_s <- (edad - calibracion$edad_media) / calibracion$edad_sd

  ## Una covariable auxiliar para los submodelos de masas.
  z <- rnorm(n_sujetos, mean = 0, sd = 1)

  ## Efecto aleatorio por sujeto.
  sigma_b <- as.numeric(params$sigma_b_true)
  b <- rnorm(n_sujetos, mean = 0, sd = sigma_b)

  ## Repeticiones dentro de sujeto.
  id <- rep(seq_len(n_sujetos), each = m)
  X_suj <- cbind(Intercepto = 1, Sexo = sexo, Edad_s = edad_s)
  Z_suj <- cbind(Intercepto = 1, z = z)
  X <- X_suj[id, , drop = FALSE]
  Z <- Z_suj[id, , drop = FALSE]

  beta <- beta_true(tau)
  phi <- as.numeric(params$phi_true)

  eta_q <- as.vector(X %*% beta) + b[id]
  theta <- inv_enlace(eta_q)

  eta0 <- as.vector(Z %*% gamma_true)
  eta1 <- as.vector(Z %*% delta_true)
  pr <- softmax_masas(eta0, eta1)

  ## Nueva distribución interior para ESTE tau.
  ## Kumaraswamy en escala log(-log y) = llt - log(-log qb) + log(-log(1 - (1 - U)^(1/phi))), U ~ U(0,1).
  ## Exacto con los tres enlaces aunque theta o y valgan 1 en doble precisión.
  lth <- lth_enlace(eta_q)
  u_interior <- runif(length(id))
  lly_interior <- lth - lqb_kum(phi, tau) + log(-log1p(-exp(log1p(-u_interior) / phi)))
  y_interior <- exp(-exp(lly_interior))

  ## Asignación a 0, 1 o componente interior.
  u_masa <- runif(length(id))
  if (escala_log) {
    ## El grupo sale del sorteo de la masa, no de comparar y con 1: un interior que vale 1 en doble
    ## precisión sigue siendo interior.
    grp <- ifelse(u_masa < pr$p0, 1L, ifelse(u_masa < pr$p0 + pr$p1, 2L, 3L))
    y <- ifelse(grp == 1L, 0, ifelse(grp == 2L, 1, y_interior))
  } else {
    y <- ifelse(
      u_masa < pr$p0,
      0,
      ifelse(u_masa < pr$p0 + pr$p1, 1, y_interior)
    )
    grp <- ifelse(y == 0, 1L, ifelse(y == 1, 2L, 3L))
  }
  interior <- which(grp == 3L)

  datos <- tibble(
    id = id,
    repeticion = rep(seq_len(m), times = n_sujetos),
    sexo = sexo[id],
    edad = edad[id],
    edad_s = edad_s[id],
    z = z[id],
    b_true = b[id],
    theta_true = theta,
    p0_true = pr$p0,
    p1_true = pr$p1,
    pc_true = pr$pc,
    y_interior = y_interior,
    y = y,
    grp = grp,
    tau_generador = tau
  )

  list(
    tau = tau,
    replica = replica,
    seed = seed_tau,
    datos = datos,
    X = X,
    Z = Z,
    id = id,
    y = y,
    grp = grp,
    interior = interior,
    Xc = X[interior, , drop = FALSE],
    idc = id[interior],
    yc = y[interior],
    llyc = lly_interior[interior],          # log(-log y), interiores
    lyc = -exp(lly_interior[interior]),     # log y, interiores (exacto aunque y valga 1)
    n = n_sujetos,
    N = length(y),
    Nc = length(interior),
    beta_true = beta,
    gamma0_true = gamma_true,
    gamma1_true = delta_true,
    phi_true = phi,
    sigma_b_true = sigma_b,
    b_true_subject = b,
    X_subject = X_suj,
    Z_subject = Z_suj
  )
}
## ---- modelo-jags ----
## v4: Unit-Gamma en escala llt = log(-log theta) con los tres enlaces, sin recortes sobre theta.
## Cada enlace calcula llt de forma estable (JAGS evalúa las dos ramas de `ifelse`, así que cada una
## recibe argumentos acotados). eta se acota en [-35, 35] (cloglog: 20) solo para que una propuesta
## extrema no produzca Inf; esos valores están lejos de cualquiera que los datos puedan sostener.
llt_jags <- list(
  logit = "
    ec[m] <- max(min(eta[m], 35), -35)
    u[m] <- exp(-ec[m])
    llt[m] <- ifelse(step(ec[m] - 5),
                     -ec[m] + log(1 - u[m] / 2 + pow(u[m], 2) / 3 - pow(u[m], 3) / 4 + pow(u[m], 4) / 5),
                     log(log(1 + u[m])))
",
  probit = "
    ec[m] <- max(min(eta[m], 35), -35)
    pe[m] <- phi(-ec[m])
    llt[m] <- ifelse(step(ec[m] - 3.1),
                     log(pe[m] + pow(pe[m], 2) / 2 + pow(pe[m], 3) / 3),
                     log(-log(phi(min(ec[m], 3.1)))))
",
  cloglog = "
    x[m] <- exp(min(eta[m], 20))
    u[m] <- exp(-x[m])
    xs[m] <- max(min(x[m], 4), 1.0E-5)
    llt[m] <- ifelse(step(x[m] - 4),
                     -x[m] + log(1 + u[m] / 2 + pow(u[m], 2) / 3 + pow(u[m], 3) / 4 + pow(u[m], 4) / 5),
                     ifelse(step(1.0E-5 - x[m]),
                            log(-eta[m] + x[m] / 2),
                            log(-log(1 - exp(-xs[m])))))
"
)

modelo_jags_ug <- "
model {
  for (k in 1:N) {
    grp[k] ~ dcat(pr[k, 1:3])
    e0[k] <- inprod(Z[k, ], gamma0[])
    e1[k] <- inprod(Z[k, ], gamma1[])
    mm[k] <- max(0, max(e0[k], e1[k]))
    den[k] <- exp(-mm[k]) + exp(e0[k] - mm[k]) + exp(e1[k] - mm[k])
    pr[k, 1] <- exp(e0[k] - mm[k]) / den[k]
    pr[k, 2] <- exp(e1[k] - mm[k]) / den[k]
    pr[k, 3] <- exp(-mm[k]) / den[k]
  }

  ## qb = 1 - (1 - tau)^(1/phi), acotado en (1e-12, 1 - 1e-12) solo para que una propuesta extrema
  ## de varphi no produzca log(0); en valores plausibles no actúa.
  lqb <- log(-log(min(max(1.0E-12, 1 - pow(1 - tau, 1 / varphi)), 1 - 1.0E-12)))
  for (m in 1:Nc) {
    eta[m] <- mu[idc[m]]
__LLT__
    ## Kumaraswamy: log a = log(-log qb) - llt; t = a (-log y) = exp(log a + log(-log y)).
    ## log(1 - y^a) = log(1 - exp(-t)): serie en t pequeño, fórmula directa si t >= 1e-4.
    la[m] <- lqb - llt[m]
    lt[m] <- la[m] + llyc[m]
    tk[m] <- exp(min(lt[m], 700))
    ts[m] <- min(tk[m], 1.0E-4)
    tg[m] <- max(tk[m], 1.0E-4)
    l1my[m] <- ifelse(step(tk[m] - 1.0E-4),
                      log(1 - exp(-tg[m])),
                      lt[m] + log(1 - ts[m] / 2 + pow(ts[m], 2) / 6))
    logdc[m] <- la[m] + log(varphi) + (exp(min(la[m], 700)) - 1) * lyc[m] + (varphi - 1) * l1my[m]
    zeros.c[m] ~ dpois(Cz - logdc[m])
  }

  ## Centrado jerárquico: mu[i] = X_i beta + b_i con b_i ~ N(0, sigma.b^2).
  for (i in 1:n) {
    mu[i] ~ dnorm(inprod(Xs[i, ], beta[]), tau.b)
  }

  for (j in 1:px) {
    beta[j] ~ dnorm(0, 0.01)
  }
  for (j in 1:pz) {
    gamma0[j] ~ dnorm(0, 0.01)
    gamma1[j] ~ dnorm(0, 0.01)
  }

  varphi ~ dgamma(0.1, 0.01)
  sigma.b ~ dunif(0, 100)
  tau.b <- pow(sigma.b, -2)
}
"
modelo_jags_ug <- sub("__LLT__", llt_jags[[enlace_nombre]], modelo_jags_ug, fixed = TRUE)
## Nombre que usa el resto del Notebook.
modelo_jags_uw <- modelo_jags_ug

## ---- funciones-ajuste ----
config_mcmc <- list(
  n_adapt = as.integer(params$n_adapt),
  n_burn = as.integer(params$n_burn),
  n_iter = as.integer(params$n_iter),
  n_thin = as.integer(params$n_thin),
  n_chains = as.integer(params$n_chains)
)

truth_vector <- function(sim) {
  c(
    `beta[1]` = unname(sim$beta_true[1]),
    `beta[2]` = unname(sim$beta_true[2]),
    `beta[3]` = unname(sim$beta_true[3]),
    `gamma0[1]` = unname(sim$gamma0_true[1]),
    `gamma0[2]` = unname(sim$gamma0_true[2]),
    `gamma1[1]` = unname(sim$gamma1_true[1]),
    `gamma1[2]` = unname(sim$gamma1_true[2]),
    varphi = sim$phi_true,
    `sigma.b` = sim$sigma_b_true
  )
}

ajustar_simulacion <- function(sim, usar_cache = TRUE) {
  ## La clave incluye el enlace y una huella de los valores generadores:
  ## sin ellos, cambiar la configuracion reutilizaria ajustes que no corresponden.
  huella <- paste(
    c(sim$beta_true, sim$gamma0_true, sim$gamma1_true, sim$phi_true, sim$sigma_b_true),
    collapse = "_"
  )
  cod <- utf8ToInt(huella)
  huella <- sprintf("%07d", sum(cod * seq_along(cod)) %% 9999991)
  key <- sprintf(
    ## "hc": modelo con centrado jerárquico; no reutiliza ajustes de la versión no centrada.
    ## "lg": cloglog en escala log(-log); no reutiliza los ajustes de cloglog con recortes (inválidos).
    "%s_%s_tau%03d_n%04d_rep%04d_%d-%d-%d-%d-%d_g%s",
    "kuv4", enlace_nombre, round(100 * sim$tau), sim$n, sim$replica,   # v4: Kumaraswamy
    config_mcmc$n_adapt, config_mcmc$n_burn, config_mcmc$n_iter,
    config_mcmc$n_thin, config_mcmc$n_chains, huella
  )
  ruta_cache <- file.path(params$cache_dir, paste0(key, ".rds"))

  if (usar_cache && file.exists(ruta_cache) && !isTRUE(params$force_refit)) {
    out <- readRDS(ruta_cache)
    out$cache_status <- "reutilizado"
    return(out)
  }

  if (sim$Nc < 10) stop("Muy pocas observaciones interiores para ajustar el modelo.")

  dat <- list(
    N = sim$N,
    Nc = sim$Nc,
    n = sim$n,
    px = ncol(sim$X),
    pz = ncol(sim$Z),
    Xs = sim$X_subject,
    Z = sim$Z,
    idc = sim$idc,
    yc = sim$yc,
    grp = sim$grp,
    tau = sim$tau,
    Cz = 1e4,
    zeros.c = rep(0, sim$Nc)
  )
  if (escala_log) {
    ## v4 (los tres enlaces): log y y log(-log y) en lugar de y; Cz mayor porque la log-densidad puede
    ## ser grande cuando theta está muy cerca de 1.
    dat$yc <- NULL
    dat$lyc <- sim$lyc
    dat$llyc <- sim$llyc
    dat$Cz <- 1e6
  }

  inits <- lapply(seq_len(config_mcmc$n_chains), function(ch) {
    set.seed(sim$seed + 5000L * ch)
    beta0 <- rnorm(ncol(sim$X), 0, 0.10)
    list(
      beta = beta0,
      gamma0 = c(rnorm(1, -2, 0.20), rnorm(ncol(sim$Z) - 1, 0, 0.10)),
      gamma1 = c(rnorm(1, -2, 0.20), rnorm(ncol(sim$Z) - 1, 0, 0.10)),
      varphi = runif(1, 2, 6),
      sigma.b = runif(1, 0.5, 1.5),
      mu = as.vector(sim$X_subject %*% beta0) + rnorm(sim$n, 0, 0.10),
      .RNG.name = "base::Wichmann-Hill",
      .RNG.seed = sim$seed + 10000L * ch
    )
  })

  inicio <- proc.time()[[3]]
  jm <- rjags::jags.model(
    textConnection(modelo_jags_uw),
    data = dat,
    inits = inits,
    n.chains = config_mcmc$n_chains,
    n.adapt = config_mcmc$n_adapt,
    quiet = TRUE
  )
  update(jm, n.iter = config_mcmc$n_burn, progress.bar = "none")
  mcmc <- rjags::coda.samples(
    jm,
    variable.names = c("beta", "gamma0", "gamma1", "varphi", "sigma.b"),
    n.iter = config_mcmc$n_iter,
    thin = config_mcmc$n_thin,
    progress.bar = "none"
  )
  elapsed <- proc.time()[[3]] - inicio

  out <- list(
    tau = sim$tau,
    n = sim$n,
    replica = sim$replica,
    mcmc = mcmc,
    truth = truth_vector(sim),
    elapsed = elapsed,
    cache_status = "calculado"
  )
  saveRDS(out, ruta_cache)
  out
}

resumir_fit <- function(fit, nivel = as.numeric(params$icr)) {
  mat <- as.matrix(fit$mcmc)
  alpha <- (1 - nivel) / 2
  pars <- names(fit$truth)
  mat <- mat[, pars, drop = FALSE]

  est <- colMeans(mat)
  sd_post <- apply(mat, 2, sd)
  li <- apply(mat, 2, quantile, probs = alpha)
  ls <- apply(mat, 2, quantile, probs = 1 - alpha)

  rhat <- tryCatch({
    ps <- coda::gelman.diag(fit$mcmc[, pars], multivariate = FALSE)$psrf[, 1]
    ps[pars]
  }, error = function(e) setNames(rep(NA_real_, length(pars)), pars))

  ess <- coda::effectiveSize(fit$mcmc[, pars])
  ess <- ess[pars]

  tibble(
    Tau = fit$tau,
    n = fit$n,
    Replica = fit$replica,
    Parametro = pars,
    Verdadero = as.numeric(fit$truth[pars]),
    Estimacion = as.numeric(est[pars]),
    DE_post = as.numeric(sd_post[pars]),
    LI = as.numeric(li[pars]),
    LS = as.numeric(ls[pars]),
    Error = Estimacion - Verdadero,
    Error_abs = abs(Error),
    Error2 = Error^2,
    Error_rel = ifelse(abs(Verdadero) < 1e-8, NA_real_, Error / Verdadero),
    Cubre = LI <= Verdadero & Verdadero <= LS,
    Longitud_ICr = LS - LI,
    Rhat = as.numeric(rhat[pars]),
    ESS = as.numeric(ess[pars]),
    Tiempo_min = fit$elapsed / 60,
    Cache = fit$cache_status
  )
}
