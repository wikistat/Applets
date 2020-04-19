init <- function(n = 20, var1 = 2, var2 = 1, cor12 = 0.8, seed = 2){
  c12 <- cor12 * sqrt(var1) * sqrt(var2)
  Sigma <- matrix(c(var1, c12, c12, var2), nrow = 2)
  n <- n
  set.seed(seed)
  x <- mvrnorm(n = n, mu = rep(0, 2), Sigma = Sigma)
  x <- scale(x, center = TRUE, scale = FALSE)
  colnames(x) <- c("x1", "x2")
  return(x)
}


inertia <- function(x){
  sum(x^2 / nrow(x))
}

inertiaPlot <- function(x, ...){
  for (i in 1:nrow(x)){
    arrows(x0 = x[i, 1], x1 = 0,
           y0 = x[i, 2], y1 = 0, code = 3, length = 0.15, ...)
  }
}

projPlot <- function(x, xproj, ...){
  for (i in 1:nrow(x)){
    lines(x = c(x[i, 1], xproj[i, 1]),
          y = c(x[i, 2], xproj[i, 2]), lty = "dotted")
  }
  points(xproj, pch = 1, cex = 0.8)
}

proj1d <- function(x, u) {
  sum(x * u) * u / sum(u ^ 2)
}

proj1dPlot <- function(x, v, plot = TRUE, col = "black", 
                       cex.points = 1, ...){
  xproj <- apply(x, 1, proj1d, u = v)
  xproj <- t(xproj)
  if (plot){
    if (v[1] == 0) {
      abline(v = 0, col = "blue", lwd = 2)
    } else { 
      abline(a = 0, b = v[2]/v[1], col = "blue", lwd = 2)
    }
    points(xproj, pch = 19, col = col, cex = cex.points)
    for (i in 1:nrow(x)){
      lines(c(x[i, 1], xproj[i, 1]), c(x[i, 2], xproj[i, 2]), 
            lty = "dotted", col = col, ...)
    }
  }
  return(xproj)
}

PCAturningLine <- function(x, solFlag, angle){
  n <- nrow(x)
  
  if (solFlag) {
    SigmaHat <- cov(x) * (n-1) / n 
    E <- eigen(SigmaHat, symmetric = TRUE)
    angleRad <- atan(E$vectors[2, 1] / E$vectors[1, 1])
  } else {
    angleRad <- angle * pi / 180
  }

  v1 <- c(cos(angleRad), sin(angleRad))
  v2 <- c(- sin(angleRad), cos(angleRad))
  xproj1 <- proj1dPlot(x = x, v = c(cos(angleRad), sin(angleRad)), plot = FALSE)
  xproj2 <- proj1dPlot(x = x, v = c(-sin(angleRad), cos(angleRad)), plot = FALSE)
  projInertia1 <- inertia(xproj1)
  projInertia2 <- inertia(xproj2)
  totalInertia <- projInertia1 + projInertia2
  
  par(mfrow = c(1, 2))
  plot(x, asp = 1,
       main = paste0("Variance of projections (axis 1): ", 
                     round(projInertia1/totalInertia*100, 2), "%"))
  proj1dPlot(x = x, v = v1, plot = TRUE)
  plot(x, asp = 1, 
       main = paste0("Variance of projections (orthogonal of axis 1): ", 
                     round(projInertia2/totalInertia*100, 2), "%"))
  proj1dPlot(x = x, v = v2)
}