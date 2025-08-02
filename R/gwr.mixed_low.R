## Mixed GWR
gwr.mixed.2 <- function(x1, x2, y, loc, out.loc, adaptive = F, bw = sqrt(var(loc[, 1]) + var(loc[, 2])),
                        kernel = "bisquare", p = 2, theta = 0, longlat = F, dMat)
{
  # gwr.fitted <- function(x,b) apply(x*b,1,sum)
  dp.n <- nrow(loc)

  ncols.2 <- dim(x2)[2]
  x3 <- NULL
  if (missing(out.loc)) {
    dMat1 <- dMat
  } else {
    dMat1 <- gw.dist(loc, p = p, theta = theta, longlat = longlat)
  }
  
  for (i in 1:ncols.2) {
    m.temp <- gwr.q(x1, x2[, i], loc,
      adaptive = adaptive, bw = bw,
      kernel = kernel, p = p, theta = theta, longlat = longlat, dMat = dMat1
    )
    x3 <- cbind(x3, x2[, i] - gw_fitted(x1, m.temp))
  }
  colnames(x3) <- colnames(x2)

  m.temp <- gwr.q(x1, y, loc,
    adaptive = adaptive, bw = bw,
    kernel = kernel, p = p, theta = theta, longlat = longlat, dMat = dMat1
  )
  y2 <- y - gw.fitted(x1, m.temp)

  model2 <- gwr.q(x3, y2, loc,
    adaptive = TRUE, bw = 1.0e6, kernel = "boxcar",
    p = p, theta = theta, longlat = longlat, dMat = dMat1
  )
  fit2 <- gw_fitted(x2, model2)

  model1 <- gwr.q(x1, y - fit2, loc,
    out.loc = out.loc, adaptive = adaptive, bw = bw,
    kernel = kernel, p = p, theta = theta, longlat = longlat, dMat = dMat
  )

  if (!missing(out.loc)) {
    model2 <- gwr.q(x3, y2, loc,
      out.loc = out.loc, adaptive = TRUE, bw = 1.0e6, kernel = "boxcar",
      p = p, theta = theta, longlat = longlat, dMat = dMat
    )
  }
  list(local = model1, global = model2)
}


# Fix the column names for x3
gwr.mixed.trace <- function(x1, x2, y, loc, out.loc, adaptive = F, bw = sqrt(var(loc[, 1]) + var(loc[, 2])),
                            kernel = "bisquare", p = 2, theta = 0, longlat = F, dMat)
{
  gw.fitted <- gwr.fitted <- function(x, b) apply(x * b, 1, sum)
  e.vec <- function(m, n) as.numeric(m == 1:n)
  dp.n <- nrow(loc)
  if (missing(dMat)) {
    DM.given <- F
  } else {
    DM.given <- T
  }

  ncols.2 <- dim(x2)[2]
  # n.items <- length(y)
  if (missing(out.loc)) {
    dMat1 <- dMat
  } else {
    dMat1 <- gw.dist(loc, p = p, theta = theta, longlat = longlat)
  }
  x3 <- NULL
  for (i in 1:ncols.2) {
    m.temp <- gwr.q(x1, x2[, i], loc,
      adaptive = adaptive, bw = bw,
      kernel = kernel, p = p, theta = theta, longlat = longlat, dMat = dMat1
    )
    x3 <- cbind(x3, x2[, i] - gw_fitted(x1, m.temp))
  }

  colnames(x3) <- colnames(x2)
  hii <- NULL

  for (i in 1:dp.n) {
    m.temp <- gwr.q(x1, e.vec(i, dp.n), loc,
      adaptive = adaptive, bw = bw,
      kernel = kernel, p = p, theta = theta, longlat = longlat, dMat = dMat1
    )
    y2 <- e.vec(i, dp.n) - gw_fitted(x1, m.temp)

    model2 <- gwr.q(x3, y2, loc,
      adaptive = TRUE, bw = 1.0e6, kernel = "boxcar",
      p = p, theta = theta, longlat = longlat, dMat = dMat1
    )
    fit2 <- gw_fitted(x2, model2)
    if (DM.given) {
      model1 <- gwr.q(x1, e.vec(i, dp.n) - fit2, loc,
        out.loc = matrix(loc[i, ], ncol = 2), adaptive = adaptive, bw = bw,
        kernel = kernel, dMat = matrix(dMat[, i], ncol = 1)
      )

      model2 <- gwr.q(x3, y2, loc,
        out.loc = matrix(loc[i, ], ncol = 2), adaptive = TRUE, bw = 1.0e6, kernel = "boxcar",
        dMat = matrix(dMat[, i], ncol = 1)
      )
    } else {
      model1 <- gwr.q(x1, e.vec(i, dp.n) - fit2, loc,
        out.loc = matrix(loc[i, ], ncol = 2), adaptive = adaptive, bw = bw,
        kernel = kernel, p = p, theta = theta, longlat = longlat
      )
      model2 <- gwr.q(x3, y2, loc,
        out.loc = matrix(loc[i, ], ncol = 2), adaptive = TRUE, bw = 1.0e6, kernel = "boxcar",
        p = p, theta = theta, longlat = longlat
      )
    }
    hii <- c(hii, gw_fitted(matrix(x1[i, ], nrow = 1), model1) + gw_fitted(matrix(x2[i, ], nrow = 1), model2))
  }
  sum(hii)
}


gwr.q <- function(
    x, y, loc, out.loc = loc, adaptive = F,
    bw = sqrt(var(loc[, 1]) + var(loc[, 2])),
    kernel, p, theta, longlat, dMat, wt2 = rep(1, nrow(loc)))
{
  if (missing(out.loc)) {
    rp.n <- nrow(loc)
  } else {
    rp.n <- nrow(out.loc)
  }
  if (missing(dMat)) {
    dMat <- gw.dist(loc, out.loc, p, theta, longlat)
  }
  var.n <- ncol(x)
  betas <- gwr_q(x, y, dMat, bw, kernel, adaptive)
  colnames(betas) <- colnames(x)
  betas
}
