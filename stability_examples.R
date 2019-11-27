

## one predator, two prey

## interaction matrix: B
Bvals <- c( 0.7,  0.4,  0.2,
           -0.2,  0.6, -0.2,
           -0.2, -0.1,  0.5)
BB <- matrix(Bvals, 3, 3, byrow = TRUE)


## covariate effects

## temp has positive effects on growth
CC <- matrix(c(0.1, 0.2, 0.3), 3, 1)

## identity matrix
II <- diag(3)


tt <- 100

## process errors
QQ <- diag(c(0.1, 0.2, 0.2))
ww <- xx <- t(MASS::mvrnorm(tt, matrix(0,3,1), QQ))

## dummy covariate: discrete sine wave
cc <- matrix(sin(2*pi*seq(tt)/12), 1, tt)
cc[51:100] <- cc[51:100] + 2


## simulate MAR(1) process
for(t in 2:tt) {
  xx[,t] <- BB %*% xx[,t-1] + CC %*% cc[t] + ww[,t]
}

matplot(t(xx), type = "l", lty = "solid", col = c("blue", "orange", "red"))

x2 <- ww

## simulate MAR(1) process
for(t in 2:tt) {
  x2[,t] <- BB %*% x2[,t-1] + CC %*% (cc[t] + 1) + ww[,t]
}

matlines(t(x2), type = "l", lty = "dotted")




long_term_change <- function(BB, CC) {
  ltc <- vector("list", 2)
  ltc[[1]] <- solve(II - BB) %*% CC
  BB <- II - BB
  nn <- dim(BB)[1]
  mm <- dim(CC)[2]
  tmp <- matrix(NA, nn, mm)
  for(i in 1:nn) {
    for(j in 1:mm) {
      if(i == 1) {
        num <- cbind(CC[,j], BB[,-i])
      } else {
        if(i == nn) {
          num <- cbind(BB[,-i], CC[,j])
        } else {
          num <- cbind(BB[,1:(i-1)], CC[,j], BB[,(i+1):nn]) 
        }
      }
      tmp[i,j] <- det(num) / det(BB)
    }
  }
  ltc[[2]] <- tmp
  return(ltc)
}
  

## long-term change in X via Ives et al. 1999
long_term_change(BB, CC)

m1 <- apply(xx[,1:50], 1, mean)
m2 <- apply(xx[,51:100], 1, mean)


## change in X per unit change in covariate via Cottingham & Klug (2001)
solve(II - BB) %*% CC

- solve(BB) %*% CC


BB[,1]
BB[,2]
BB[,3]

BB[,1] %*% BB[,2]






x <- c(1,0,0)
y <- c(sqrt(0.5),sqrt(0.5),0)


xprod <- function(x, y) {
  i1 <- c(2,3,1)
  i2 <- c(3,1,2)
  return (x[i1]*y[i2] - x[i2]*y[i1])
}

xprod(x,y)

x %*% y



## compute cross product to get orthogonal vector O
ovec <- xprod(BB[,1], BB[,2])
## normalize to unit vector z
z <- ovec / sqrt(sum(ovec^2))
## get projection of B_i onto z
Bz <- z * as.vector(BB[,3] %*% z)
## get projection of C_j on z
Cz <- CC %*% z


ee <- c(0,0,0)

clr <- c("red", "orange", "blue")

x <- c(1,0,0)
y <- c(0,1,0)
z <- c(0,0,1)

arrows3D(ee, ee, ee, x1 = x, y1 = y, z1 = z,  
         phi = 0, theta = 45, lwd = 2,
         col = clr, bty = "n", type = "triangle")

text3D(x, y, z, c(expression(B[1]), expression(B[2]), expression(B[3])),
       phi = 0, theta = 45,
       col = NULL, NAcol = "white",  breaks = NULL,
       colkey = NULL, panel.first = NULL, 
       clim = NULL, clab = NULL, 
       bty = "n", add = TRUE, plot = TRUE)


ee <- rep(0,6)
aa <- cbind(diag(3), BB)

x <- aa[1,]
y <- aa[2,]
z <- aa[3,]

segments3D(ee, ee, ee, x1 = x, y1 = y, z1 = z,  
         phi = 0, theta = 45, lwd = 2,
         col = c(rep("gray",3), clr), bty = "n")

text3D(x, y, z, c("x", "y", "z",
                  expression(B[1]), expression(B[2]), expression(B[3])),
       phi = 0, theta = 45,
       col = NULL, NAcol = "white",  breaks = NULL,
       colkey = NULL, panel.first = NULL, 
       clim = NULL, clab = NULL, 
       bty = "n", add = TRUE, plot = TRUE)


segments3D(0, 0, 0, x1 = ovec[1], y1 = ovec[2], z1 = 1.2*ovec[3],  
         phi = 0, theta = 60, lwd = 2, add = TRUE,
         col = "black", bty = "n", type = "triangle")
segments3D(0, 0, 0, x1 = -ovec[1], y1 = -ovec[2], z1 = -1.2*ovec[3],  
         phi = 0, theta = 60, lwd = 2, add = TRUE,
         col = "darkgray", bty = "n", type = "triangle")

