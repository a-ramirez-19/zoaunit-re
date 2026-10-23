## =====================================================================================================
## Simulación 2 de ZOAUNIT-RE (22/10/2026): escenarios pedidos por el director en la
## reunión del 21/10/2026. Reemplaza la simulación con generador común (antigua simulación 3, _archivo/generador_comun/funciones_sim3.R, 18/10/2026).
## Fuente única de funciones: la usan run_bloque_sim2.R (GitHub Actions) y resumir_sim2.R (portátil),
## con las mismas claves de caché.
##
## Diseño de generación (decisión del 22/10/2026; documento "Simulación 2: decisiones de diseño"): el de la
## simulación 1 (en las diapositivas del 21/10, «simulación 2»). Para cada familia, enlace y orden cuantílico se genera un conjunto de datos propio
## con ese mismo modelo, de modo que todos los parámetros tienen valor verdadero, y se ajusta el modelo
## correcto (misma familia, enlace y tau). Se descartó el generador común de la antigua simulación 3: con
## logit, probit y cloglog deja sin valor verdadero a beta, sigma_b y phi.
##
## Lo que pidió el director el 21/10/2026:
##   * Quedarse con tres distribuciones (Unit-Gamma, Unit-Weibull y Kumaraswamy) y ampliar los escenarios
##     en lugar de agregar distribuciones:
##       - porcentajes de inflación: el de la simulación 1 (0.10 en cero, 0.08 en uno) y dos extremos, 20 % / 30 %
##         y 30 % / 20 % ("situaciones extremas que uno debería considerar");
##       - parámetro de forma/precisión phi de cada familia: tres configuraciones que cambian la forma de
##         la distribución;
##       - dispersión del efecto aleatorio sigma_b: baja, la de la simulación 1 (1.0) y alta (1.5).
##     "Había que hacer una combinación entre todos los escenarios": diseño factorial completo
##     3 x 3 x 3 = 27 escenarios por familia, enlace y tau.
##   * Justificar los valores escogidos (ver el bloque de parámetros).
##
## Matriz (decisión del 22/10/2026): familias ug, uw, ku; enlaces logit, probit, cloglog; tau = 0.25, 0.50
## y 0.75; n = 30, 50, 100, 300 y 500 sujetos (n = 200 se quitó el 22/10/2026) con 4 observaciones cada uno; R = 500 sin piloto
## (decisión del 22/10/2026), lanzadas en tramos de 100 réplicas.
##
## Generador (el de la simulación 1, ZOAUNIT_RE_simulacion1*.qmd): para el sujeto i y la observación j,
##   eta_ij = x_ij' beta + b_i,  b_i ~ N(0, sigma_b^2),  theta_ij = g^{-1}(eta_ij)  (enlace g del ajuste),
##   componente interior: log(-log Y) = log(-log theta_ij) - c_d(tau, phi) + c_d(U, phi),  U ~ U(0, 1),
## de modo que Q_d(tau) = theta_ij exactamente; las masas en 0 y 1 salen del logit multinomial.
## Semillas de los datos: las de la simulación 1 (semilla + 100000 * réplica + 10000 * tau + n), iguales para todos
## los escenarios, enlaces y familias (números aleatorios comunes): las diferencias entre escenarios no
## se mezclan con diferencias de las muestras, y el escenario base usa los mismos sorteos que la simulación 1.
## =====================================================================================================

suppressPackageStartupMessages({
  library(coda)
  library(rjags)
})
load.module("glm", quiet = TRUE)

params_sim2_defecto <- list(
  ## ---------------- Factores del director (21/10/2026) ----------------
  ## Inflación (p0, p1) en el perfil de referencia z = 0: la de la simulación 1 (cercana a los datos, 9.8 % y 8.1 %)
  ## y los dos extremos que pidió el director. En los extremos la mitad de las observaciones está en las
  ## fronteras.
  masas_i1008 = "0.10,0.08",
  masas_i2030 = "0.20,0.30",
  masas_i3020 = "0.30,0.20",
  ## Desviación estándar del intercepto aleatorio b_i (escala del predictor eta).
  sigma_s05 = 0.5,
  sigma_s10 = 1.0,               # la de la simulación 1
  sigma_s15 = 1.5,               # la que propuso el director
  ## Forma phi_d. Los tres niveles se definen por la dispersión de la componente interior dentro del
  ## sujeto, medida como la DE de log(-log Y) dado el predictor (escala común a las tres familias, porque
  ## phi no significa lo mismo en todas). Esa DE solo depende de phi:
  ##   Unit-Gamma:   sqrt(trigamma(phi));  Unit-Weibull: pi / (sqrt(6) phi);  Kumaraswamy: numérica.
  ## fM = el phi de la simulación 1 (calibrado con los datos; DE 0.605 / 0.641 / 0.600 en ug / uw / ku);
  ## fB = DE 0.40 (componente interior concentrada, alta precisión);
  ## fA = DE 0.90 (componente interior dispersa, baja precisión).
  ## Los niveles bajo y alto están a la misma distancia de la DE de los datos en escala logarítmica
  ## (x 2/3 y x 3/2). Valores de phi redondeados a un decimal (ver niveles_phi() más abajo).
  phi_fB_ug = 6.7, phi_fM_ug = 3.2, phi_fA_ug = 1.7,
  phi_fB_uw = 3.2, phi_fM_uw = 2.0, phi_fA_uw = 1.4,
  phi_fB_ku = 11.1, phi_fM_ku = 3.8, phi_fA_ku = 1.7,
  ## Coeficientes del predictor cuantílico (los de la simulación 1, iguales para los tres tau): intercepto, sexo y
  ## edad estandarizada. Son valores ilustrativos: con los datos periodontales el efecto de sexo es negativo
  ## (logit: -0.76 a -0.51 entre tau = 0.10 y 0.90). Pendiente de confirmar con el director.
  beta = "0.50,0.60,0.40",
  ## Pendientes de z en las masas (las de la simulación 1).
  gamma_z = 0.35, delta_z = -0.35,
  ## ---------------- Previas ----------------
  ## Fijadas por el director (comentario 32): beta, gamma y delta ~ N(0, 2.5^2); sigma_b ~ half-normal(0, 1).
  previa_sd_coef = 2.5,
  escala_sigma = 1,
  ## Previa de phi propia de cada familia: Gamma(2, tasa), con la mediana en la media geométrica de los
  ## niveles fB y fA, de modo que los tres niveles queden dentro del intervalo central del 90 % de la previa
  ## (con la previa de la simulación con generador común, el nivel fB de Kumaraswamy quedaba en el percentil 99).
  previa_phi_ug = "2,0.5",
  previa_phi_uw = "2,0.8",
  previa_phi_ku = "2,0.39",
  ## Alternativa para el análisis de sensibilidad: escalas de las previas x 2 (variante previa_alt1).
  sens_sd_coef = "5",
  sens_escala_sigma = "2",
  sens_previa_phi_ug = "2,0.25",
  sens_previa_phi_uw = "2,0.4",
  sens_previa_phi_ku = "2,0.195",
  ## Perfiles de covariables donde se mide el error de los cuantiles (sexo, edad_s): hombre y mujer de
  ## edad media, y mujeres con edad -1 y +1 DE (los de la simulación con generador común).
  perfiles = "0,0;1,0;1,-1;1,1",
  ## ---------------- Protocolo de convergencia del director (comentario 23) ----------------
  n_chains = 4,
  rhat_max = 1.01,
  ess_min = 400,
  m_repeticiones = 4,            # observaciones por sujeto
  n_adapt = 750, n_burn = 1200, n_iter = 3000, n_thin = 3,
  ## Si no se cumple el protocolo, las mismas cadenas continúan n_iter iteraciones más, hasta ext_max veces.
  ext_max = 4,
  ## "completo" guarda todas las extracciones; "liviano" (GitHub Actions) guarda el resumen calculado con
  ## todas las extracciones y, en las réplicas 1 a `replicas_traza`, `extracciones_traza` por cadena para
  ## las trazas (con ~2 millones de ajustes en el estudio principal, guardar trazas de todos ocuparía
  ## más de 150 GB).
  guardar = "completo",
  extracciones_traza = 250,
  replicas_traza = 5,            # en modo liviano, réplicas 1-5 con trazas; las demás solo con el resumen
  semilla = 20260927,            # la de la simulación 1
  icr = 0.95,
  ## ---------------- Ejecución ----------------
  familia = "ug",
  enlace = "logit",
  cache_dir = "_cache_sim2",
  force_refit = FALSE
)
## Configuración efectiva P3: los valores por defecto, modificados por `params_sim2` si existe.
P3 <- if (exists("params_sim2") && is.list(params_sim2)) {
  modifyList(params_sim2_defecto, params_sim2[intersect(names(params_sim2), names(params_sim2_defecto))])
} else params_sim2_defecto

if (!exists("calibracion", inherits = FALSE)) {
  arch_cal <- Sys.getenv("CALIBRACION", "calibracion.rds")
  calibracion <- if (file.exists(arch_cal)) readRDS(arch_cal) else
    list(edad_media = 55.31, edad_sd = 10.64, p_mujer = 0.7621, fuente = "valores de respaldo")
}

## Se mantiene por compatibilidad con los scripts de la simulación con generador común: aquí no hay valores por definir.
validar_configuracion <- function(...) invisible(TRUE)

## ---------------------------------------------------------------------------------------------------
## Variantes del análisis de sensibilidad
## ---------------------------------------------------------------------------------------------------
## "director": la configuración base. "previa_alt<k>": la alternativa k de las previas (sens_*). Las previas
## entran en la huella de la clave de caché, así que cada variante queda con claves distintas.
partes_sens <- function(x) trimws(strsplit(as.character(x), ";", fixed = TRUE)[[1]])
variantes_sim2 <- function(familia, base = P3) {
  claves <- c("sens_sd_coef", "sens_escala_sigma", paste0("sens_previa_phi_", familia))
  k_prev <- max(vapply(claves, function(k) length(partes_sens(base[[k]])), 1L))
  c("director", paste0("previa_alt", seq_len(k_prev)))
}
aplicar_variante <- function(variante, familia, base = P3) {
  if (identical(variante, "director")) return(base)
  out <- base
  if (grepl("^previa_alt[0-9]+$", variante)) {
    k <- as.integer(sub("previa_alt", "", variante))
    tomar <- function(clave) { v <- partes_sens(base[[clave]]); v[min(k, length(v))] }
    out$previa_sd_coef <- as.numeric(tomar("sens_sd_coef"))
    out$escala_sigma <- as.numeric(tomar("sens_escala_sigma"))
    out[[paste0("previa_phi_", familia)]] <- tomar(paste0("sens_previa_phi_", familia))
  } else stop("Variante no reconocida: ", variante, call. = FALSE)
  out
}
con_variante <- function(variante, familia, codigo) {
  viejo <- P3
  on.exit(P3 <<- viejo, add = TRUE)
  P3 <<- aplicar_variante(variante, familia, viejo)
  force(codigo)
}

## ---------------------------------------------------------------------------------------------------
## Utilidades
## ---------------------------------------------------------------------------------------------------
leer_numeros <- function(x, n_esperado = NULL) {
  z <- as.numeric(trimws(strsplit(as.character(x), ",", fixed = TRUE)[[1]]))
  if (any(!is.finite(z))) stop("Se encontraron valores no numéricos en: ", x)
  if (!is.null(n_esperado) && length(z) != n_esperado) stop("Se esperaban ", n_esperado, " valores en: ", x)
  z
}
softmax_masas <- function(eta0, eta1) {
  mm <- pmax(eta0, eta1, 0)
  den <- exp(-mm) + exp(eta0 - mm) + exp(eta1 - mm)
  list(p0 = exp(eta0 - mm) / den, p1 = exp(eta1 - mm) / den, pc = exp(-mm) / den)
}

## Enlaces: g^{-1} y log(-log theta) estable (de la simulación 1).
inv_enlace <- function(eta, enlace) {
  switch(enlace, logit = plogis(eta), probit = pnorm(eta), cloglog = 1 - exp(-exp(eta)))
}
g_enlace <- function(theta, enlace) {
  switch(enlace, logit = qlogis(theta), probit = qnorm(theta), cloglog = log(-log1p(-theta)))
}
lth_cloglog <- function(eta) {
  x <- exp(eta)
  ifelse(x < 700, log(-log1p(-exp(-pmin(x, 700)))), -x)
}
lth_enlace <- function(eta, enlace) {
  switch(
    enlace,
    logit = log(log1p(exp(-pmax(eta, -700)))),
    probit = log(-pnorm(eta, log.p = TRUE)),
    cloglog = {
      x <- exp(pmin(eta, 700))
      suppressWarnings(ifelse(x > 1, lth_cloglog(eta), log(-log(-expm1(-x)))))
    },
    stop("Enlace no reconocido: ", enlace)
  )
}

## ---------------------------------------------------------------------------------------------------
## Familias: desplazamiento c(u, phi) en escala log(-log y)
## ---------------------------------------------------------------------------------------------------
c_cuantil <- function(u, phi, familia) {
  switch(
    familia,
    ug = log(qgamma(u, shape = phi, rate = 1, lower.tail = FALSE)),       # q_Gamma(1 - u)
    uw = log(-log(u)) / phi,
    ku = log(-log(-expm1(log1p(-u) / phi))),
    stop("Familia no reconocida: ", familia)
  )
}
nombre_familia <- c(ug = "Unit-Gamma", uw = "Unit-Weibull", ku = "Kumaraswamy")

## DE de log(-log Y) dentro del sujeto implicada por phi (criterio de los niveles fB, fM y fA).
de_interior <- function(phi, familia, k = 20000L) {
  switch(familia,
         ug = sqrt(trigamma(phi)),
         uw = pi / (sqrt(6) * phi),
         ku = { u <- (seq_len(k) - 0.5) / k; sd(c_cuantil(u, phi, "ku")) },
         stop("Familia no reconocida: ", familia))
}
## Tabla de los niveles de phi: valor, DE implicada y percentil en la previa de la familia.
niveles_phi <- function(familias = c("ug", "uw", "ku")) {
  do.call(rbind, lapply(familias, function(f) {
    pp <- leer_numeros(P3[[paste0("previa_phi_", f)]], 2)
    do.call(rbind, lapply(c("fB", "fM", "fA"), function(l) {
      phi <- as.numeric(P3[[sprintf("phi_%s_%s", l, f)]])
      data.frame(Familia = f, Nivel = l, phi = phi, DE_loglog = de_interior(phi, f),
                 Percentil_previa = pgamma(phi, pp[1], pp[2]))
    }))
  }))
}

## ---------------------------------------------------------------------------------------------------
## Escenarios: inflación (3) x sigma_b (3) x phi (3) = 27, factorial completo
## ---------------------------------------------------------------------------------------------------
## Identificador: i<p0p1>_s<sigma x 10>_f<nivel>; el escenario base (el de la simulación 1) es i1008_s10_fM.
escenarios_sim2 <- function() {
  g <- expand.grid(inflacion = c("i1008", "i2030", "i3020"), sigma = c("s05", "s10", "s15"),
                   forma = c("fB", "fM", "fA"), stringsAsFactors = FALSE)
  g$id <- sprintf("%s_%s_%s", g$inflacion, g$sigma, g$forma)
  g$base <- g$id == "i1008_s10_fM"
  g
}

## Perfiles donde se mide el error de los cuantiles (b = 0, z = 0), leídos de P3$perfiles.
perfiles_sim2 <- function() {
  m <- do.call(rbind, lapply(strsplit(strsplit(P3$perfiles, ";")[[1]], ","), function(x) as.numeric(trimws(x))))
  d <- data.frame(sexo = m[, 1], edad_s = m[, 2])
  d$etiqueta <- sprintf("sexo=%g,edad_s=%+g", d$sexo, d$edad_s)
  d
}

definir_escenario <- function(id, familia) {
  g <- escenarios_sim2(); fila <- g[g$id == id, ]
  if (!nrow(fila)) stop("Escenario no reconocido: ", id)
  if (!familia %in% names(nombre_familia)) stop("Familia no reconocida: ", familia)
  masas <- leer_numeros(P3[[paste0("masas_", fila$inflacion)]], 2)
  pc <- 1 - sum(masas)
  if (pc <= 0) stop("p0 + p1 debe ser menor que 1 en masas_", fila$inflacion)
  list(
    id = id, familia = familia,
    beta = setNames(leer_numeros(P3$beta, 3), c("Intercepto", "Sexo", "Edad_s")),
    phi = as.numeric(P3[[sprintf("phi_%s_%s", fila$forma, familia)]]),
    sigma_b = as.numeric(P3[[paste0("sigma_", fila$sigma)]]),
    gamma0 = c(Intercepto = log(masas[1] / pc), z = as.numeric(P3$gamma_z)),
    gamma1 = c(Intercepto = log(masas[2] / pc), z = as.numeric(P3$delta_z)),
    inflacion = fila$inflacion, sigma = fila$sigma, forma = fila$forma
  )
}

## ---------------------------------------------------------------------------------------------------
## Generador bien especificado (simulación 1): un conjunto de datos por familia, enlace y tau
## ---------------------------------------------------------------------------------------------------
## Mismo orden de sorteos y semillas que la simulación 1 (ZOAUNIT_RE_simulacion1*.qmd): sexo, edad, z, b, U interior,
## U masa.
simular_dataset <- function(esc, tau, enlace, n_sujetos, replica, m = as.integer(P3$m_repeticiones),
                            seed_base = as.integer(P3$semilla)) {
  tau <- as.numeric(tau); n_sujetos <- as.integer(n_sujetos); replica <- as.integer(replica)
  seed <- seed_base + 100000L * replica + as.integer(round(10000 * tau)) + n_sujetos
  set.seed(seed)
  sexo <- rbinom(n_sujetos, 1, calibracion$p_mujer)
  edad <- rnorm(n_sujetos, calibracion$edad_media, calibracion$edad_sd)
  edad_s <- (edad - calibracion$edad_media) / calibracion$edad_sd
  z <- rnorm(n_sujetos)
  b <- rnorm(n_sujetos, 0, esc$sigma_b)
  id <- rep(seq_len(n_sujetos), each = m)
  X_suj <- cbind(Intercepto = 1, Sexo = sexo, Edad_s = edad_s)
  Z_suj <- cbind(Intercepto = 1, z = z)
  X <- X_suj[id, , drop = FALSE]; Z <- Z_suj[id, , drop = FALSE]
  eta <- as.vector(X %*% esc$beta) + b[id]
  ## Componente interior: log(-log Y) = log(-log theta) - c_d(tau, phi) + c_d(U, phi); Q_d(tau) = theta.
  u_int <- runif(length(id))
  lly <- lth_enlace(eta, enlace) - c_cuantil(tau, esc$phi, esc$familia) + c_cuantil(u_int, esc$phi, esc$familia)
  pr <- softmax_masas(as.vector(Z %*% esc$gamma0), as.vector(Z %*% esc$gamma1))
  u_masa <- runif(length(id))
  ## El grupo sale del sorteo de la masa: un interior que vale 1 en doble precisión sigue siendo interior.
  grp <- ifelse(u_masa < pr$p0, 1L, ifelse(u_masa < pr$p0 + pr$p1, 2L, 3L))
  y <- ifelse(grp == 1L, 0, ifelse(grp == 2L, 1, exp(-exp(lly))))
  interior <- which(grp == 3L)
  list(
    escenario = esc, tau = tau, enlace = enlace, replica = replica, seed = seed, n = n_sujetos,
    N = length(y), Nc = length(interior),
    id = id, y = y, grp = grp, interior = interior, Z = Z, X_subject = X_suj, Z_subject = Z_suj,
    idc = id[interior], llyc = lly[interior], lyc = -exp(lly[interior]),
    eta = eta, b = b
  )
}

## ---------------------------------------------------------------------------------------------------
## Valores verdaderos: todos los parámetros (el modelo ajustado es el generador)
## ---------------------------------------------------------------------------------------------------
parametros_verdaderos <- function(esc) {
  c(`beta[1]` = unname(esc$beta[1]), `beta[2]` = unname(esc$beta[2]), `beta[3]` = unname(esc$beta[3]),
    `gamma0[1]` = unname(esc$gamma0[1]), `gamma0[2]` = unname(esc$gamma0[2]),
    `gamma1[1]` = unname(esc$gamma1[1]), `gamma1[2]` = unname(esc$gamma1[2]),
    varphi = esc$phi, sigma.b = esc$sigma_b)
}
## Cuantiles verdaderos en los perfiles (b = 0, z = 0): theta_tau de la componente interior y Q_Y(tau) de
## la respuesta completa (cero o uno cuando la masa supera el orden).
cuantiles_verdaderos <- function(esc, tau, enlace, perfiles = perfiles_sim2()) {
  eta <- esc$beta[1] + esc$beta[2] * perfiles$sexo + esc$beta[3] * perfiles$edad_s
  pr <- softmax_masas(esc$gamma0[1], esc$gamma1[1])
  theta <- inv_enlace(eta, enlace)
  lam <- lth_enlace(eta, enlace) - c_cuantil(tau, esc$phi, esc$familia)
  u <- (tau - pr$p0) / pr$pc
  qy <- if (tau <= pr$p0) rep(0, length(eta)) else if (tau >= 1 - pr$p1) rep(1, length(eta)) else
    exp(-exp(lam + c_cuantil(u, esc$phi, esc$familia)))
  c(setNames(theta, paste0("theta[", perfiles$etiqueta, "]")),
    setNames(qy, paste0("QY[", perfiles$etiqueta, "]")))
}

## Cuantiles del modelo ajustado (tau, enlace) en los perfiles, para cada extracción posterior.
cuantiles_ajustados <- function(beta, g0, g1, phi, tau, enlace, familia, perfiles = perfiles_sim2()) {
  Xp <- cbind(1, perfiles$sexo, perfiles$edad_s)
  eta <- beta %*% t(Xp)
  pr <- softmax_masas(g0[, 1], g1[, 1])
  out <- list()
  for (k in seq_len(nrow(perfiles))) {
    theta <- inv_enlace(eta[, k], enlace)
    lam <- lth_enlace(eta[, k], enlace) - c_cuantil(tau, phi, familia)
    u <- pmin(pmax((tau - pr$p0) / pr$pc, 1e-12), 1 - 1e-12)
    qy <- ifelse(tau <= pr$p0, 0, ifelse(tau >= 1 - pr$p1, 1, exp(-exp(lam + c_cuantil(u, phi, familia)))))
    out[[paste0("theta[", perfiles$etiqueta[k], "]")]] <- theta
    out[[paste0("QY[", perfiles$etiqueta[k], "]")]] <- qy
  }
  out
}

## ---------------------------------------------------------------------------------------------------
## Modelo JAGS (parametrización cuantílica de la simulación 1; previas nuevas)
## ---------------------------------------------------------------------------------------------------
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
## Log-densidad interior de cada familia en escala llt (idéntica a la simulación 1).
densidad_jags <- list(
  ug = list(pre = "  lgq <- log(qgamma(1 - tau, varphi, 1))\n", cuerpo = "
    lb[m] <- lgq - llt[m]
    logdc[m] <- varphi * lb[m] - loggam(varphi) + (varphi - 1) * llyc[m] +
                (exp(min(lb[m], 700)) - 1) * lyc[m]
"),
  uw = list(pre = "", cuerpo = "
    lr[m] <- llyc[m] - llt[m]
    logdc[m] <- log(varphi) + log(-log(tau)) - llt[m] + (varphi - 1) * lr[m] -
                lyc[m] + log(tau) * exp(min(varphi * lr[m], 700))
"),
  ku = list(pre = "  lqb <- log(-log(min(max(1.0E-12, 1 - pow(1 - tau, 1 / varphi)), 1 - 1.0E-12)))\n", cuerpo = "
    la[m] <- lqb - llt[m]
    lt[m] <- la[m] + llyc[m]
    tk[m] <- exp(min(lt[m], 700))
    ts[m] <- min(tk[m], 1.0E-4)
    tg[m] <- max(tk[m], 1.0E-4)
    l1my[m] <- ifelse(step(tk[m] - 1.0E-4),
                      log(1 - exp(-tg[m])),
                      lt[m] + log(1 - ts[m] / 2 + pow(ts[m], 2) / 6))
    logdc[m] <- la[m] + log(varphi) + (exp(min(la[m], 700)) - 1) * lyc[m] + (varphi - 1) * l1my[m]
")
)

modelo_jags_sim2 <- function(familia, enlace) {
  prec <- 1 / as.numeric(P3$previa_sd_coef)^2
  pphi <- leer_numeros(P3[[paste0("previa_phi_", familia)]], 2)
  paste0("
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
", densidad_jags[[familia]]$pre, "  for (m in 1:Nc) {
    eta[m] <- mu[idc[m]]
", llt_jags[[enlace]], densidad_jags[[familia]]$cuerpo, "    zeros.c[m] ~ dpois(Cz - logdc[m])
  }
  ## Centrado jerárquico: mu[i] = X_i beta + b_i, b_i ~ N(0, sigma.b^2).
  for (i in 1:n) {
    mu[i] ~ dnorm(inprod(Xs[i, ], beta[]), tau.b)
  }
  ## Previas (comentario 32 del director): N(0, ", P3$previa_sd_coef, "^2), half-normal(0, ", P3$escala_sigma,
  "), Gamma(", pphi[1], ", ", pphi[2], ").
  for (j in 1:px) { beta[j] ~ dnorm(0, ", prec, ") }
  for (j in 1:pz) {
    gamma0[j] ~ dnorm(0, ", prec, ")
    gamma1[j] ~ dnorm(0, ", prec, ")
  }
  varphi ~ dgamma(", pphi[1], ", ", pphi[2], ")
  sigma.b ~ dnorm(0, ", 1 / as.numeric(P3$escala_sigma)^2, ") T(0, )
  tau.b <- pow(sigma.b, -2)
}
")
}

## ---------------------------------------------------------------------------------------------------
## Ajuste con caché
## ---------------------------------------------------------------------------------------------------
config_mcmc <- function() list(
  n_adapt = as.integer(P3$n_adapt), n_burn = as.integer(P3$n_burn),
  n_iter = as.integer(P3$n_iter), n_thin = as.integer(P3$n_thin),
  n_chains = as.integer(P3$n_chains), ext_max = as.integer(P3$ext_max)
)
huella_sim2 <- function(esc) {
  h <- paste(signif(c(esc$beta, esc$phi, esc$sigma_b, esc$gamma0, esc$gamma1,
                      as.numeric(P3$previa_sd_coef), as.numeric(P3$escala_sigma),
                      leer_numeros(P3[[paste0("previa_phi_", esc$familia)]], 2)), 10), collapse = "_")
  cod <- utf8ToInt(h)
  sprintf("%07d", sum(cod * seq_along(cod)) %% 9999991)
}
clave_sim2 <- function(familia, enlace, esc_id, tau, n, replica) {
  esc <- definir_escenario(esc_id, familia); cm <- config_mcmc()
  sprintf("s2%s_%s_%s_tau%03d_n%04d_rep%04d_%d-%d-%d-%d-%d-x%d_g%s",
          familia, enlace, esc_id, round(100 * tau), as.integer(n), as.integer(replica),
          cm$n_adapt, cm$n_burn, cm$n_iter, cm$n_thin, cm$n_chains, cm$ext_max, huella_sim2(esc))
}

ajustar_sim2 <- function(sim, tau = sim$tau, enlace = sim$enlace, usar_cache = TRUE) {
  ## El modelo ajustado es el generador: mismo tau y mismo enlace que los datos.
  if (!isTRUE(all.equal(tau, sim$tau)) || !identical(enlace, sim$enlace))
    stop("Los datos se generaron con tau = ", sim$tau, " y enlace ", sim$enlace, ".", call. = FALSE)
  esc <- sim$escenario; cm <- config_mcmc()
  key <- clave_sim2(esc$familia, enlace, esc$id, tau, sim$n, sim$replica)
  ruta <- file.path(P3$cache_dir, paste0(key, ".rds"))
  if (usar_cache && file.exists(ruta) && !isTRUE(P3$force_refit)) {
    out <- readRDS(ruta); out$cache_status <- "reutilizado"; return(out)
  }
  if (sim$Nc < 10) stop("Muy pocas observaciones interiores para ajustar el modelo.")
  dat <- list(N = sim$N, Nc = sim$Nc, n = sim$n, px = 3L, pz = 2L, Xs = sim$X_subject, Z = sim$Z,
              idc = sim$idc, lyc = sim$lyc, llyc = sim$llyc, grp = sim$grp, tau = tau, Cz = 1e6,
              zeros.c = rep(0, sim$Nc))
  semilla_aj <- sim$seed + as.integer(round(1000 * tau)) + 7L * match(enlace, c("logit", "probit", "cloglog"))
  inits <- lapply(seq_len(cm$n_chains), function(ch) {
    set.seed(semilla_aj + 5000L * ch)
    beta0 <- rnorm(3, 0, 0.10)
    list(beta = beta0,
         gamma0 = c(rnorm(1, -2, 0.20), rnorm(1, 0, 0.10)),
         gamma1 = c(rnorm(1, -2, 0.20), rnorm(1, 0, 0.10)),
         varphi = runif(1, 1, 6), sigma.b = runif(1, 0.3, 1.2),
         mu = as.vector(sim$X_subject %*% beta0) + rnorm(sim$n, 0, 0.10),
         .RNG.name = "base::Wichmann-Hill", .RNG.seed = semilla_aj + 10000L * ch)
  })
  dir.create(P3$cache_dir, showWarnings = FALSE, recursive = TRUE)
  inicio <- proc.time()[[3]]
  vars <- c("beta", "gamma0", "gamma1", "varphi", "sigma.b")
  ## Un error de JAGS (p. ej., una densidad no finita durante la adaptación) no detiene el bloque: la
  ## réplica queda guardada como inestable, con el mensaje, y se cuenta en las métricas de estabilidad.
  error_jags <- NULL
  extensiones <- 0L
  mcmc <- tryCatch({
    jm <- rjags::jags.model(textConnection(modelo_jags_sim2(esc$familia, enlace)), data = dat, inits = inits,
                            n.chains = cm$n_chains, n.adapt = cm$n_adapt, quiet = TRUE)
    update(jm, n.iter = cm$n_burn, progress.bar = "none")
    mc <- rjags::coda.samples(jm, variable.names = vars, n.iter = cm$n_iter, thin = cm$n_thin, progress.bar = "none")
    ## Extensión de las cadenas mientras no se cumpla el protocolo de convergencia.
    cumple <- function(mc) {
      d <- vapply(parametros_traza, function(v) {
        m <- sapply(mc, function(ch) as.numeric(ch[, v]))
        if (!all(is.finite(m))) return(c(Inf, 0, 0))
        c(posterior::rhat(m), posterior::ess_bulk(m), posterior::ess_tail(m))
      }, numeric(3))
      max(d[1, ]) <= as.numeric(P3$rhat_max) && min(d[2:3, ]) >= as.numeric(P3$ess_min)
    }
    while (extensiones < cm$ext_max && !cumple(mc)) {
      extra <- rjags::coda.samples(jm, variable.names = vars, n.iter = cm$n_iter, thin = cm$n_thin, progress.bar = "none")
      mc <- coda::as.mcmc.list(lapply(seq_along(mc), function(ch)
        coda::mcmc(rbind(as.matrix(mc[[ch]]), as.matrix(extra[[ch]])), start = start(mc[[ch]]), thin = cm$n_thin)))
      extensiones <- extensiones + 1L
    }
    mc
  }, error = function(e) { error_jags <<- conditionMessage(e); NULL })
  out <- list(
    simulacion = "sim2", familia = esc$familia, enlace = enlace, escenario = esc, tau = tau,
    n = sim$n, replica = sim$replica, Nc = sim$Nc, mcmc = mcmc, extensiones = extensiones,
    elapsed = proc.time()[[3]] - inicio, cache_status = "calculado", clave = key,
    error = error_jags,
    ## Comentario 22: cota de la log-densidad interior frente a Cz (intensidad de Poisson no negativa).
    ceros = if (is.null(mcmc)) NULL else chequeo_truco_ceros(esc$familia, as.vector(sapply(mcmc, function(ch)
      as.numeric(ch[, "varphi"]))), -sim$lyc, dat$Cz),
    ## Comentario 34: semillas, configuración y entorno de cómputo.
    reproducibilidad = list(semilla_datos = sim$seed, semilla_ajuste = semilla_aj, config_mcmc = cm,
                            previas = list(sd_coef = as.numeric(P3$previa_sd_coef),
                                           escala_sigma = as.numeric(P3$escala_sigma),
                                           phi = P3[[paste0("previa_phi_", esc$familia)]]),
                            entorno = info_entorno())
  )
  if (identical(as.character(P3$guardar), "liviano")) {
    out$resumen <- resumir_ajuste_sim2(out)
    ## Trazas solo en las primeras réplicas (diagnóstico visual); en las demás basta el resumen, que ya
    ## tiene R-hat y ESS calculados con todas las extracciones. El entorno de cómputo queda en registro/.
    if (!is.null(out$mcmc)) out$mcmc <- if (sim$replica <= as.integer(P3$replicas_traza))
      adelgazar_mcmc(out$mcmc, as.integer(P3$extracciones_traza)) else NULL
    out$reproducibilidad$entorno <- out$reproducibilidad$entorno[c("R", "JAGS", "commit", "corrida")]
    out$liviano <- TRUE
  }
  saveRDS(out, ruta)
  out
}

## ---------------------------------------------------------------------------------------------------
## Resumen de un ajuste (parámetros y cuantiles) y métricas Monte Carlo
## ---------------------------------------------------------------------------------------------------
## Entorno de cómputo (comentario 34): versiones de R, JAGS y paquetes, sistema y, en GitHub Actions,
## imagen del ejecutor, repositorio, commit y corrida. Se calcula una vez por proceso.
info_entorno <- local({
  guardado <- NULL
  function() {
    if (is.null(guardado)) {
      ver <- function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
      guardado <<- list(
        R = R.version.string, JAGS = tryCatch(as.character(rjags::jags.version()), error = function(e) NA_character_),
        paquetes = vapply(c("rjags", "coda", "posterior"), ver, character(1)),
        sistema = paste(Sys.info()[c("sysname", "release")], collapse = " "),
        imagen = paste(Sys.getenv("ImageOS"), Sys.getenv("ImageVersion")),
        repositorio = Sys.getenv("GITHUB_REPOSITORY"), commit = Sys.getenv("GITHUB_SHA"),
        corrida = Sys.getenv("GITHUB_RUN_ID"), inicio_proceso = format(Sys.time(), tz = "UTC", usetz = TRUE)
      )
    }
    guardado
  }
})

## Comprobación del truco de los ceros (comentario 22). JAGS usa zeros.c[m] ~ dpois(Cz - logdc[m]); la media
## debe ser no negativa en todo el muestreo. Condición suficiente: Cz mayor que el máximo de la log-densidad
## interior sobre TODO valor de theta (es decir, de eta y b_i), para los phi muestreados y los y interiores.
## Con t = -log y, ese máximo es h(phi) + t - log t:
##   Unit-Gamma:   h = phi log phi - phi - lgamma(phi)
##   Unit-Weibull: h = log phi - 1
##   Kumaraswamy:  h = log phi + max_w {log w - w + (phi - 1) log(1 - e^-w)}
chequeo_truco_ceros <- function(familia, phi, t, Cz) {
  phi <- phi[is.finite(phi) & phi > 0]
  if (!length(phi) || !length(t)) return(list(cota = NA_real_, margen = NA_real_, ok = NA))
  if (length(phi) > 200) phi <- c(range(phi), quantile(phi, seq(0, 1, length.out = 198), names = FALSE))
  h <- switch(
    familia,
    ug = max(phi * log(phi) - phi - lgamma(phi)),
    uw = max(log(phi) - 1),
    ku = max(vapply(phi, function(p) log(p) + optimize(function(lw)
      lw - exp(lw) + (p - 1) * log(-expm1(-exp(lw))), c(-30, 10), maximum = TRUE)$objective, numeric(1))),
    stop("Familia no reconocida: ", familia)
  )
  cota <- h + max(t - log(t))
  list(cota = cota, margen = Cz - cota, ok = cota <= Cz)
}

## Deja a lo sumo `k` extracciones por cadena (para guardar trazas sin todo el tamaño del ajuste).
adelgazar_mcmc <- function(mc, k) {
  coda::as.mcmc.list(lapply(mc, function(ch) {
    m <- as.matrix(ch); paso <- max(1L, ceiling(nrow(m) / k))
    coda::mcmc(m[seq(1, nrow(m), by = paso), , drop = FALSE], start = start(ch), thin = coda::thin(ch) * paso)
  }))
}

resumir_ajuste_sim2 <- function(fit, nivel = as.numeric(P3$icr), perfiles = perfiles_sim2()) {
  ## Ajustes livianos: el resumen ya se calculó con todas las extracciones.
  if (!is.null(fit$resumen)) return(fit$resumen)
  alpha <- (1 - nivel) / 2
  pars <- c("beta[1]", "beta[2]", "beta[3]", "gamma0[1]", "gamma0[2]", "gamma1[1]", "gamma1[2]",
            "varphi", "sigma.b")
  ceros_ok <- if (is.null(fit$ceros)) NA else fit$ceros$ok
  ## Ajuste con error de JAGS: filas sin estimaciones, contadas como inestables.
  if (!is.null(fit$error) || is.null(fit$mcmc)) {
    verdad <- c(parametros_verdaderos(fit$escenario), cuantiles_verdaderos(fit$escenario, fit$tau, fit$enlace, perfiles))
    obj <- union(pars, names(verdad))
    tipo <- ifelse(startsWith(obj, "theta["), "cuantil_interior", ifelse(startsWith(obj, "QY["), "cuantil_total",
             ifelse(obj %in% names(verdad), "parametro_exacto", "parametro_sin_verdadero")))
    return(data.frame(
      Familia = fit$familia, Enlace = fit$enlace, Escenario = fit$escenario$id, Tau = fit$tau, n = fit$n,
      Replica = fit$replica, Tipo = tipo, Objetivo = obj,
      Verdadero = unname(ifelse(obj %in% names(verdad), verdad[obj], NA_real_)),
      Estimacion = NA_real_, Mediana = NA_real_, DE_post = NA_real_, LI = NA_real_, LS = NA_real_,
      Rhat = NA_real_, ESS_bulk = NA_real_, ESS_tail = NA_real_, Finitos = FALSE, N_cadenas = NA_integer_,
      Tiempo_min = fit$elapsed / 60, Extensiones = fit$extensiones, Ceros_ok = ceros_ok,
      Error = if (is.null(fit$error)) "sin extracciones" else fit$error, stringsAsFactors = FALSE
    ))
  }
  arr <- lapply(pars, function(v) sapply(fit$mcmc, function(ch) as.numeric(ch[, v]))); names(arr) <- pars
  plano <- function(v) as.vector(arr[[v]])
  q_post <- cuantiles_ajustados(cbind(plano("beta[1]"), plano("beta[2]"), plano("beta[3]")),
                                cbind(plano("gamma0[1]"), plano("gamma0[2]")),
                                cbind(plano("gamma1[1]"), plano("gamma1[2]")),
                                plano("varphi"), fit$tau, fit$enlace, fit$familia, perfiles)
  ncad <- length(fit$mcmc)
  objetivos <- c(arr, lapply(q_post, function(x) matrix(x, ncol = ncad)))
  verdad <- c(parametros_verdaderos(fit$escenario), cuantiles_verdaderos(fit$escenario, fit$tau, fit$enlace, perfiles))
  filas <- lapply(names(objetivos), function(v) {
    m <- objetivos[[v]]; x <- as.vector(m); ok <- all(is.finite(x))
    tipo <- if (startsWith(v, "theta[")) "cuantil_interior" else if (startsWith(v, "QY[")) "cuantil_total" else
      if (v %in% names(verdad)) "parametro_exacto" else "parametro_sin_verdadero"
    data.frame(
      Familia = fit$familia, Enlace = fit$enlace, Escenario = fit$escenario$id, Tau = fit$tau, n = fit$n,
      Replica = fit$replica, Tipo = tipo, Objetivo = v,
      Verdadero = if (v %in% names(verdad)) unname(verdad[v]) else NA_real_,
      Estimacion = mean(x), Mediana = median(x), DE_post = sd(x),
      LI = unname(quantile(x, alpha, na.rm = TRUE)), LS = unname(quantile(x, 1 - alpha, na.rm = TRUE)),
      Rhat = if (ok) posterior::rhat(m) else NA_real_,
      ESS_bulk = if (ok) posterior::ess_bulk(m) else NA_real_,
      ESS_tail = if (ok) posterior::ess_tail(m) else NA_real_,
      Finitos = ok, N_cadenas = ncad, Tiempo_min = fit$elapsed / 60,
      Extensiones = if (is.null(fit$extensiones)) NA_integer_ else fit$extensiones,
      Ceros_ok = ceros_ok, Error = NA_character_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, filas)
}

metricas_celda <- function(d) {
  R_total <- nrow(d)
  d <- d[is.finite(d$Estimacion) & is.finite(d$LI) & is.finite(d$LS) & is.finite(d$Verdadero), ]
  R <- nrow(d)
  if (R < 2) {
    nombres <- c("Verdadero", "Sesgo", "MCSE_Sesgo", "RB", "MCSE_RB", "MSE", "MCSE_MSE", "RMSE", "MCSE_RMSE",
                 "CP", "MCSE_CP", "Longitud_ICr", "MCSE_Longitud", "Rhat_mediana", "ESS_bulk_mediana")
    return(cbind(data.frame(R = R, R_excluidas = R_total - R),
                 as.data.frame(setNames(as.list(rep(NA_real_, length(nombres))), nombres))))
  }
  e <- d$Estimacion - d$Verdadero; v <- d$Verdadero[1]
  cubre <- d$LI <= d$Verdadero & d$Verdadero <= d$LS; lon <- d$LS - d$LI
  sesgo <- mean(e); mse <- mean(e^2); rmse <- sqrt(mse); cp <- mean(cubre); rel <- abs(v) > 1e-8
  data.frame(
    R = R, R_excluidas = R_total - R, Verdadero = v,
    Sesgo = sesgo, MCSE_Sesgo = sd(e) / sqrt(R),
    RB = if (rel) sesgo / v else NA_real_, MCSE_RB = if (rel) sd(e) / sqrt(R) / abs(v) else NA_real_,
    MSE = mse, MCSE_MSE = sd(e^2) / sqrt(R),
    RMSE = rmse, MCSE_RMSE = if (rmse > 0) sd(e^2) / sqrt(R) / (2 * rmse) else NA_real_,
    CP = cp, MCSE_CP = sqrt(cp * (1 - cp) / R),
    Longitud_ICr = mean(lon), MCSE_Longitud = sd(lon) / sqrt(R),
    Rhat_mediana = median(d$Rhat, na.rm = TRUE), ESS_bulk_mediana = median(d$ESS_bulk, na.rm = TRUE)
  )
}

## Convergencia por réplica con el protocolo del director (R-hat <= 1.01; ESS bulk y tail >= 400).
## Una réplica se cuenta como numéricamente inestable si produjo extracciones no finitas o un error de JAGS.
## Métricas de estabilidad adicionales (diseño reformulado): error de JAGS, extensiones usadas, sigma_b en la
## frontera (mediana posterior < umbral_sigma), truco de los ceros cumplido y tiempo por ajuste.
convergencia_replica <- function(res, umbral_sigma = 0.05) {
  p <- res[res$Tipo %in% c("parametro_exacto", "parametro_sin_verdadero"), ]
  clave <- do.call(paste, c(p[, c("Familia", "Enlace", "Escenario", "Tau", "n", "Replica")], sep = "|"))
  agg <- function(x, f) tapply(x, clave, f)
  u <- !duplicated(clave)
  out <- p[u, c("Familia", "Enlace", "Escenario", "Tau", "n", "Replica")]
  k <- clave[u]
  out$Rhat_max <- as.numeric(agg(p$Rhat, max)[k])
  out$ESS_bulk_min <- as.numeric(agg(p$ESS_bulk, min)[k])
  out$ESS_tail_min <- as.numeric(agg(p$ESS_tail, min)[k])
  out$Finitos <- as.logical(agg(p$Finitos, all)[k])
  out$Extensiones <- if ("Extensiones" %in% names(p)) as.numeric(agg(p$Extensiones, max)[k]) else NA_real_
  out$Conv_director <- out$Finitos & out$Rhat_max <= as.numeric(P3$rhat_max) &
    out$ESS_bulk_min >= as.numeric(P3$ess_min) & out$ESS_tail_min >= as.numeric(P3$ess_min)
  out$Conv_director[is.na(out$Conv_director)] <- FALSE
  out$Error_JAGS <- if ("Error" %in% names(p)) as.logical(agg(!is.na(p$Error), any)[k]) else FALSE
  out$Inestable <- !out$Finitos | out$Error_JAGS
  out$Conv_director[out$Inestable] <- FALSE
  sb <- p[p$Objetivo == "sigma.b", ]
  med_sb <- if ("Mediana" %in% names(sb)) sb$Mediana else rep(NA_real_, nrow(sb))
  names(med_sb) <- do.call(paste, c(sb[, c("Familia", "Enlace", "Escenario", "Tau", "n", "Replica")], sep = "|"))
  out$Sigma_b_mediana <- unname(med_sb[k])
  out$Sigma_b_frontera <- !is.na(out$Sigma_b_mediana) & out$Sigma_b_mediana < umbral_sigma
  out$Ceros_ok <- if ("Ceros_ok" %in% names(p)) as.logical(agg(p$Ceros_ok, function(z) if (all(is.na(z))) NA else all(z, na.rm = TRUE))[k]) else NA
  out$Tiempo_min <- if ("Tiempo_min" %in% names(p)) as.numeric(agg(p$Tiempo_min, max)[k]) else NA_real_
  rownames(out) <- NULL
  out
}

## Trazas de un ajuste (comentario 23: "sin patrones de deriva o falta de mezcla"): una fila por
## extracción guardada, cadena y parámetro. La iteración es la posterior al calentamiento.
parametros_traza <- c("beta[1]", "beta[2]", "beta[3]", "gamma0[1]", "gamma0[2]", "gamma1[1]", "gamma1[2]",
                      "varphi", "sigma.b")
## Etiquetas (plotmath) con la notación del informe.
etiquetas_traza <- c(`beta[1]` = "beta[0]", `beta[2]` = "beta[1]~(sexo)", `beta[3]` = "beta[2]~(edad)",
                     `gamma0[1]` = "gamma[0]", `gamma0[2]` = "gamma[1]", `gamma1[1]` = "delta[0]",
                     `gamma1[2]` = "delta[1]", varphi = "phi", `sigma.b` = "sigma[b]")
trazas_ajuste <- function(fit, pars = parametros_traza) {
  thin <- max(1L, as.integer(coda::thin(fit$mcmc)))
  do.call(rbind, lapply(seq_along(fit$mcmc), function(ch) {
    m <- as.matrix(fit$mcmc[[ch]])[, pars, drop = FALSE]
    data.frame(Iteracion = rep(seq_len(nrow(m)) * thin, ncol(m)), Cadena = factor(ch),
               Parametro = rep(pars, each = nrow(m)), Valor = as.vector(m), stringsAsFactors = FALSE)
  }))
}

## ---------------------------------------------------------------------------------------------------
## Simulación predictiva a priori (comentario 32)
## ---------------------------------------------------------------------------------------------------
## "Antes de ajustar los modelos a los datos observados, se realizarán simulaciones predictivas a priori
## para identificar especificaciones que generen concentraciones extremas o comportamientos incompatibles
## con el fenómeno estudiado." Para cada extracción de las previas se generan respuestas con covariables
## simuladas como en el generador y se resumen las masas en 0 y 1 y la concentración de la componente
## interior cerca de 0 y de 1.
predictiva_a_priori <- function(familia, enlace, tau, n_extracciones = 1000, n_sujetos = 100, semilla = 1,
                                sd_coef = as.numeric(P3$previa_sd_coef), escala_sigma = as.numeric(P3$escala_sigma),
                                previa_phi = P3[[paste0("previa_phi_", familia)]]) {
  if (identical(as.character(previa_phi), "POR_DEFINIR")) stop("Defina la previa de phi de ", familia, ".", call. = FALSE)
  pphi <- leer_numeros(previa_phi, 2)
  set.seed(semilla)
  m <- as.integer(P3$m_repeticiones)
  id <- rep(seq_len(n_sujetos), each = m)
  sexo <- rbinom(n_sujetos, 1, calibracion$p_mujer)
  edad_s <- rnorm(n_sujetos)
  z <- rnorm(n_sujetos)
  X <- cbind(1, sexo, edad_s)[id, ]; Z <- cbind(1, z)[id, ]
  filas <- lapply(seq_len(n_extracciones), function(s) {
    beta <- rnorm(3, 0, sd_coef); g0 <- rnorm(2, 0, sd_coef); g1 <- rnorm(2, 0, sd_coef)
    sig <- abs(rnorm(1, 0, escala_sigma)); phi <- rgamma(1, pphi[1], pphi[2])
    b <- rnorm(n_sujetos, 0, sig)
    eta <- as.vector(X %*% beta) + b[id]
    lam <- lth_enlace(eta, enlace) - c_cuantil(tau, phi, familia)
    y <- exp(-exp(lam + c_cuantil(runif(length(eta)), phi, familia)))
    pr <- softmax_masas(as.vector(Z %*% g0), as.vector(Z %*% g1))
    data.frame(phi = phi, sigma_b = sig, p0_media = mean(pr$p0), p1_media = mean(pr$p1),
               mediana_interior = median(y), frac_interior_menor_001 = mean(y < 0.01),
               frac_interior_mayor_099 = mean(y > 0.99))
  })
  do.call(rbind, filas)
}
