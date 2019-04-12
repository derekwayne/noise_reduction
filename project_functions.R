# Fast Walsh-Hadamard Transform
fwht2d <- function(x) {
  h <- 1
  length <- ncol(x)
  while (h < length) {
    for (i in seq(1, length, by=h*2)) {
      for (j in seq(i,i+h-1)) {
        a <- x[,j]
        b <- x[,j+h]
        x[,j] <- a + b
        x[,j+h] <- a - b
      }
    }
    h <- 2*h
  }
  h <- 1
  length <- nrow(x)
  while (h < length) {
    for (i in seq(1, length, by=h*2)) {
      for (j in seq(i,i+h-1)) {
        a <- x[j, ]
        b <- x[j+h, ]
        x[j, ] <- a + b
        x[j+h, ] <- a - b
      }
    }
    h <- 2*h
  }
  x
}

hard.thresh <- function(x, lambda) {ifelse(abs(x) > lambda, x, 0)}
soft.thresh <- function(x, lambda) {
  sign(x) * (abs(x) >= lambda) * (abs(x) - lambda)
}

sub.im.transform <- function(image, size, Th.type="soft", thresh) {
  # image size must be square
  N <- ncol(image)
  
  if ((N %% size) != 0) {
    stop("image size is not divisible by specified sub-image size")
  }
  # create a list of submatrices to apply the transform to
  flat <- as.vector(image)
  subimages <- list()
  A <- N^2/size^2 # index span for each subimage 
  for (j in 0:(A-1)) {
    subimages[[j+1]] <- matrix(flat[(j*(size^2)+1):((j+1)*(size^2))],ncol=size)
  }
  # apply transform to each sub-image and reconstruct the original matrix/image
  subimages.hat <- list()
  i <- 1
  for (im in subimages) {
    Zhat <- fwht2d(im)
    if (Th.type == "soft") {
      Zhat.star <- soft.thresh(Zhat, thresh)
    } else {
      Zhat.star <- hard.thresh(Zhat, thresh)
    }
    Zstar <- fwht2d(Zhat.star) / size
    subimages.hat[[i]] <- Zstar
    i <- i + 1
  }
  
  # reconstruct original image
  original <- c()
  for (i in 1:A) {
    original <- c(original, as.vector(subimages.hat[[i]]))
  }
  original <- matrix(original, ncol = 256)
  original
}