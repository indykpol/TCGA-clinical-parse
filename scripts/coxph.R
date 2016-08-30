function (formula, data, weights, subset, na.action, init, control, 
  ties = c("efron", "breslow", "exact"), singular.ok = TRUE, 
  robust = FALSE, model = FALSE, x = FALSE, y = TRUE, tt, 
  method = ties, ...) 
{
  ties <- match.arg(ties)
  Call <- match.call()
  indx <- match(c("formula", "data", "weights", "subset", 
    "na.action"), names(Call), nomatch = 0)
  if (indx[1] == 0) 
    stop("A formula argument is required")
  temp <- Call[c(1, indx)]
  temp[[1L]] <- quote(stats::model.frame)
  special <- c("strata", "cluster", "tt")
  temp$formula <- if (missing(data)) 
    terms(formula, special)
  else terms(formula, special, data = data)
  if (!is.null(attr(temp$formula, "specials")$tt)) {
    coxenv <- new.env(parent = environment(formula))
    assign("tt", function(x) x, env = coxenv)
    environment(temp$formula) <- coxenv
  }
  mf <- eval(temp, parent.frame())
  if (nrow(mf) == 0) 
    stop("No (non-missing) observations")
  Terms <- terms(mf)
  extraArgs <- list(...)
  if (length(extraArgs)) {
    controlargs <- names(formals(coxph.control))
    indx <- pmatch(names(extraArgs), controlargs, nomatch = 0L)
    if (any(indx == 0L)) 
      stop(gettextf("Argument %s not matched", names(extraArgs)[indx == 
        0L]), domain = NA)
  }
  if (missing(control)) 
    control <- coxph.control(...)
  Y <- model.extract(mf, "response")
  if (!inherits(Y, "Surv")) 
    stop("Response must be a survival object")
  type <- attr(Y, "type")
  if (type != "right" && type != "counting") 
    stop(paste("Cox model doesn't support \"", type, "\" survival data", 
      sep = ""))
  data.n <- nrow(Y)
  if (length(attr(Terms, "variables")) > 2) {
    ytemp <- terms.inner(formula[1:2])
    xtemp <- terms.inner(formula[-2])
    if (any(!is.na(match(xtemp, ytemp)))) 
      warning("a variable appears on both the left and right sides of the formula")
  }
  strats <- attr(Terms, "specials")$strata
  if (length(strats)) {
    stemp <- untangle.specials(Terms, "strata", 1)
    if (length(stemp$vars) == 1) 
      strata.keep <- mf[[stemp$vars]]
    else strata.keep <- strata(mf[, stemp$vars], shortlabel = TRUE)
    strats <- as.numeric(strata.keep)
  }
  timetrans <- attr(Terms, "specials")$tt
  if (missing(tt)) 
    tt <- NULL
  if (length(timetrans)) {
    timetrans <- untangle.specials(Terms, "tt")
    ntrans <- length(timetrans$terms)
    if (is.null(tt)) {
      tt <- function(x, time, riskset, weights) {
        obrien <- function(x) {
          r <- rank(x)
          (r - 0.5)/(0.5 + length(r) - r)
        }
        unlist(tapply(x, riskset, obrien))
      }
    }
    if (is.function(tt)) 
      tt <- list(tt)
    if (is.list(tt)) {
      if (any(!sapply(tt, is.function))) 
        stop("The tt argument must contain function or list of functions")
      if (length(tt) != ntrans) {
        if (length(tt) == 1) {
          temp <- vector("list", ntrans)
          for (i in 1:ntrans) temp[[i]] <- tt[[1]]
          tt <- temp
        }
        else stop("Wrong length for tt argument")
      }
    }
    else stop("The tt argument must contain a function or list of functions")
    if (ncol(Y) == 2) {
      if (length(strats) == 0) {
        sorted <- order(-Y[, 1], Y[, 2])
        newstrat <- rep.int(0L, nrow(Y))
        newstrat[1] <- 1L
      }
      else {
        sorted <- order(strats, -Y[, 1], Y[, 2])
        newstrat <- as.integer(c(1, 1 * (diff(strats[sorted]) != 
          0)))
      }
      if (storage.mode(Y) != "double") 
        storage.mode(Y) <- "double"
      counts <- .Call(Ccoxcount1, Y[sorted, ], as.integer(newstrat))
      tindex <- sorted[counts$index]
    }
    else {
      if (length(strats) == 0) {
        sort.end <- order(-Y[, 2], Y[, 3])
        sort.start <- order(-Y[, 1])
        newstrat <- c(1L, rep(0, nrow(Y) - 1))
      }
      else {
        sort.end <- order(strats, -Y[, 2], Y[, 3])
        sort.start <- order(strats, -Y[, 1])
        newstrat <- c(1L, as.integer(diff(strats[sort.end]) != 
          0))
      }
      if (storage.mode(Y) != "double") 
        storage.mode(Y) <- "double"
      counts <- .Call(Ccoxcount2, Y, as.integer(sort.start - 
        1L), as.integer(sort.end - 1L), as.integer(newstrat))
      tindex <- counts$index
    }
    mf <- mf[tindex, ]
    Y <- Surv(rep(counts$time, counts$nrisk), counts$status)
    type <- "right"
    strats <- rep(1:length(counts$nrisk), counts$nrisk)
    weights <- model.weights(mf)
    if (!is.null(weights) && any(!is.finite(weights))) 
      stop("weights must be finite")
    tcall <- attr(Terms, "variables")[timetrans$terms + 
      2]
    pvars <- attr(Terms, "predvars")
    pmethod <- sub("makepredictcall.", "", as.vector(methods("makepredictcall")))
    for (i in 1:ntrans) {
      newtt <- (tt[[i]])(mf[[timetrans$var[i]]], Y[, 1], 
        strats, weights)
      mf[[timetrans$var[i]]] <- newtt
      nclass <- class(newtt)
      if (any(nclass %in% pmethod)) {
        dummy <- as.call(list(as.name(class(newtt)[1]), 
          tcall[[i]][[2]]))
        ptemp <- makepredictcall(newtt, dummy)
        pvars[[timetrans$terms[i] + 2]] <- ptemp
      }
    }
    attr(Terms, "predvars") <- pvars
  }
  cluster <- attr(Terms, "specials")$cluster
  if (length(cluster)) {
    robust <- TRUE
    tempc <- untangle.specials(Terms, "cluster", 1:10)
    ord <- attr(Terms, "order")[tempc$terms]
    if (any(ord > 1)) 
      stop("Cluster can not be used in an interaction")
    cluster <- strata(mf[, tempc$vars], shortlabel = TRUE)
    dropterms <- tempc$terms
    xlevels <- .getXlevels(Terms[-tempc$terms], mf)
  }
  else {
    dropterms <- NULL
    if (missing(robust)) 
      robust <- FALSE
    xlevels <- .getXlevels(Terms, mf)
  }
  contrast.arg <- NULL
  attr(Terms, "intercept") <- 1
  adrop <- 0
  stemp <- untangle.specials(Terms, "strata", 1)
  if (length(stemp$vars) > 0) {
    hasinteractions <- FALSE
    for (i in stemp$vars) {
      if (any(attr(Terms, "order")[attr(Terms, "factors")[i, 
        ] > 0] > 1)) 
        hasinteractions <- TRUE
    }
    if (!hasinteractions) 
      dropterms <- c(dropterms, stemp$terms)
    else adrop <- c(0, match(stemp$var, colnames(attr(Terms, 
      "factors"))))
  }
  if (length(dropterms)) {
    temppred <- attr(terms, "predvars")
    Terms2 <- Terms[-dropterms]
    if (!is.null(temppred)) {
      attr(Terms2, "predvars") <- temppred[-(1 + dropterms)]
    }
    X <- model.matrix(Terms2, mf, constrasts = contrast.arg)
    renumber <- match(colnames(attr(Terms2, "factors")), 
      colnames(attr(Terms, "factors")))
    attr(X, "assign") <- c(0, renumber)[1 + attr(X, "assign")]
  }
  else X <- model.matrix(Terms, mf, contrasts = contrast.arg)
  Xatt <- attributes(X)
  xdrop <- Xatt$assign %in% adrop
  X <- X[, !xdrop, drop = FALSE]
  attr(X, "assign") <- Xatt$assign[!xdrop]
  attr(X, "contrasts") <- Xatt$contrasts
  offset <- model.offset(mf)
  if (is.null(offset) | all(offset == 0)) 
    offset <- rep(0, nrow(mf))
  else if (any(!is.finite(offset))) 
    stop("offsets must be finite")
  weights <- model.weights(mf)
  if (!is.null(weights) && any(!is.finite(weights))) 
    stop("weights must be finite")
  assign <- attrassign(X, Terms)
  contr.save <- attr(X, "contrasts")
  if (missing(init)) 
    init <- NULL
  else {
    if (length(init) != ncol(X)) 
      stop("wrong length for init argument")
    temp <- X %*% init - sum(colMeans(X) * init)
    if (any(temp < .Machine$double.min.exp | temp > .Machine$double.max.exp)) 
      stop("initial values lead to overflow or underflow of the exp function")
  }
  pterms <- sapply(mf, inherits, "coxph.penalty")
  if (any(pterms)) {
    pattr <- lapply(mf[pterms], attributes)
    pname <- names(pterms)[pterms]
    ord <- attr(Terms, "order")[match(pname, attr(Terms, 
      "term.labels"))]
    if (any(ord > 1)) 
      stop("Penalty terms cannot be in an interaction")
    pcols <- assign[match(pname, names(assign))]
    fit <- coxpenal.fit(X, Y, strats, offset, init = init, 
      control, weights = weights, method = method, row.names(mf), 
      pcols, pattr, assign)
  }
  else {
    if (method == "breslow" || method == "efron") {
      if (type == "right") 
        fitter <- get("coxph.fit")
      else fitter <- get("agreg.fit")
    }
    else if (method == "exact") {
      if (type == "right") 
        fitter <- get("coxexact.fit")
      else fitter <- get("agexact.fit")
    }
    else stop(paste("Unknown method", method))
    fit <- fitter(X, Y, strats, offset, init, control, weights = weights, 
      method = method, row.names(mf))
  }
  if (is.character(fit)) {
    fit <- list(fail = fit)
    class(fit) <- "coxph"
  }
  else {
    if (!is.null(fit$coefficients) && any(is.na(fit$coefficients))) {
      vars <- (1:length(fit$coefficients))[is.na(fit$coefficients)]
      msg <- paste("X matrix deemed to be singular; variable", 
        paste(vars, collapse = " "))
      if (singular.ok) 
        warning(msg)
      else stop(msg)
    }
    fit$n <- data.n
    fit$nevent <- sum(Y[, ncol(Y)])
    fit$terms <- Terms
    fit$assign <- assign
    class(fit) <- fit$method
    if (robust) {
      fit$naive.var <- fit$var
      fit$method <- method
      fit2 <- c(fit, list(x = X, y = Y, weights = weights))
      if (length(strats)) 
        fit2$strata <- strats
      if (length(cluster)) {
        temp <- residuals.coxph(fit2, type = "dfbeta", 
          collapse = cluster, weighted = TRUE)
        if (is.null(init)) 
          fit2$linear.predictors <- 0 * fit$linear.predictors
        else fit2$linear.predictors <- c(X %*% init)
        temp0 <- residuals.coxph(fit2, type = "score", 
          collapse = cluster, weighted = TRUE)
      }
      else {
        temp <- residuals.coxph(fit2, type = "dfbeta", 
          weighted = TRUE)
        fit2$linear.predictors <- 0 * fit$linear.predictors
        temp0 <- residuals.coxph(fit2, type = "score", 
          weighted = TRUE)
      }
      fit$var <- t(temp) %*% temp
      u <- apply(as.matrix(temp0), 2, sum)
      fit$rscore <- coxph.wtest(t(temp0) %*% temp0, u, 
        control$toler.chol)$test
    }
    if (length(fit$coefficients) && is.null(fit$wald.test)) {
      nabeta <- !is.na(fit$coefficients)
      if (is.null(init)) 
        temp <- fit$coefficients[nabeta]
      else temp <- (fit$coefficients - init[1:length(fit$coefficients)])[nabeta]
      fit$wald.test <- coxph.wtest(fit$var[nabeta, nabeta], 
        temp, control$toler.chol)$test
    }
    na.action <- attr(mf, "na.action")
    if (length(na.action)) 
      fit$na.action <- na.action
    if (model) {
      if (length(timetrans)) {
        mf[[".surv."]] <- Y
        mf[[".strata."]] <- strats
        stop("Time transform + model frame: code incomplete")
      }
      fit$model <- mf
    }
    if (x) {
      fit$x <- X
      if (length(strats)) {
        if (length(timetrans)) 
          fit$strata <- strats
        else fit$strata <- strata.keep
      }
    }
    if (y) 
      fit$y <- Y
  }
  if (!is.null(weights) && any(weights != 1)) 
    fit$weights <- weights
  names(fit$means) <- names(fit$coefficients)
  fit$formula <- formula(Terms)
  if (length(xlevels) > 0) 
    fit$xlevels <- xlevels
  fit$contrasts <- contr.save
  if (any(offset != 0)) 
    fit$offset <- offset
  fit$call <- Call
  fit$method <- method
  fit
}
