LDAplot <- function(nClass = 5, n = 20, rhoW = 0.8, varW = 0.1, 
                    pred = FALSE, seed = 1){

  set.seed(seed)
  cex.points <- 0.8
  cex.centroids <- 2

  rhoB <- 0
  B <- matrix(c(1, rhoB, rhoB, 1), nrow = 2)
  W <- matrix(varW * c(1, rhoW, rhoW, 1), nrow = 2)

  d <- 2
  g <- mvrnorm(n = nClass, mu = rep(0, d), Sigma = B)

  cloud <- centroids <- c()
  What <- matrix(0, d, d)
  
  for (ell in 1:nClass){
    currentCloud <- mvrnorm(n = n, mu = g[ell, ], Sigma = W) 
    currentCentroid <- colMeans(currentCloud)
    currentWhat <- cov(currentCloud) * (n - 1) / n
    cloud <- rbind(cloud, cbind(currentCloud, rep(ell, n)))
    centroids <- rbind(centroids, c(currentCentroid, ell))
    What <- What + currentWhat
  }
  What <- What / nClass
  Bhat <- cov(centroids[, 1:2]) * (nClass - 1)/nClass
  cloud[, 1:2] <- scale(cloud[, 1:2], scale = FALSE)
  centroids[, 1:2] <- scale(centroids[, 1:2], scale = FALSE)

  mydata <- data.frame(cloud)
  centroids <- data.frame(centroids)
  names(mydata) <- names(centroids) <- c("x1", "x2", "class")

  ## Sphere the data
  Weig <- eigen(What)
  Wroot <- Weig$vectors %*% diag(sqrt(Weig$values)) %*% t(Weig$vectors)
  WrootInv <- Weig$vectors %*% diag(1/sqrt(Weig$values)) %*% t(Weig$vectors)
  # verification: sum(abs(WrootInv %*% WrootInv - solve(What)))
  mydataSphered <- mydata
  mydataSphered[, 1:2] <- as.matrix(mydata[, 1:2]) %*% WrootInv
  centroidsSphered <- centroids
  centroidsSphered[, 1:2] <- as.matrix(centroids[, 1:2]) %*% WrootInv

  ## PCA
  S <- WrootInv %*% Bhat %*% WrootInv  # Bhat %*% solve(What)
  vSphered <- eigen(S)$vectors
  v <- Wroot %*% vSphered

  if (pred){
    limSphered <- c(min(mydataSphered[, 1:2]), max(mydataSphered[, 1:2]))
    tesselSphered <- deldir(x = centroidsSphered$x1, y = centroidsSphered$x2,
                     rw = c(limSphered, limSphered) * 100)
    tessel <- tesselSphered
    tessel$dirsgs[, 1:2] <- as.matrix(tessel$dirsgs[, 1:2]) %*% Wroot
    tessel$dirsgs[, 3:4] <- as.matrix(tessel$dirsgs[, 3:4]) %*% Wroot
  }
      
  # Plots
  t <- c(-100, 100) 
  par(mfrow = c(1, 2))
  plot(mydata[1:2], asp = 1, main = "Original data", 
       col = mydata$class, pch = 20, cex = cex.points)
  points(centroids, pch = 19, col = centroids$class, cex = cex.centroids)
  if (pred){
    plot(tessel, asp = 1, wlines = "tess", wpoints = "none", 
         lty = 1, cex = 2, add = TRUE, xlim = limSphered, ylim = limSphered)
  } else {
    lines(t * v[1, 1], t * v[2, 1], col = "blue")
    lines(t * v[1, 2], t * v[2, 2], col = "blue", lty = "dotted")
  }
  plot(mydataSphered[1:2], asp = 1, main = "Sphered data",
         col = mydataSphered$class, pch = 20, cex = cex.points)
  points(centroidsSphered, pch = 19, col = centroids$class, cex = cex.centroids)
  if (pred){
    plot(tesselSphered, asp = 1, wlines = "tess", wpoints = "none", 
         lty = 1, cex = 2, add = TRUE, xlim = limSphered, ylim = limSphered)
  } else {
    lines(t * vSphered[1, 1], t * vSphered[2, 1], col = "blue")
    lines(t * vSphered[1, 2], t * vSphered[2, 2], col = "blue", lty = "dotted")
  }
} # End function
