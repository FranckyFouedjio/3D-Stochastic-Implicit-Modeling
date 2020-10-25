############################################################################################
#                                                                                          #
#                                FFTMA FUNCTIONS                                           #
#                                                                                          #
############################################################################################


grid.prep <- function(W, M, N, L, ext = 3) {
  cell.width <- diff(W$xrange)/M
  cell.height <- diff(W$yrange)/N
  cell.length <- diff(W$zrange)/L
  
  mgrid <- seq(W$xrange[1], W$xrange[2], by = cell.width)
  ngrid <- seq(W$yrange[1], W$yrange[2], by = cell.height)
  lgrid <- seq(W$zrange[1], W$zrange[2], by = cell.length)
  mcens <- (mgrid + 0.5 * cell.width)[-(M + 1)]
  ncens <- (ngrid + 0.5 * cell.height)[-(N + 1)]
  lcens <- (ngrid + 0.5 * cell.length)[-(L + 1)]
  
  if (ext <= 1) 
    mgrid.ext <- ngrid.ext <-lgrid.ext <- mcens.ext <- ncens.ext <- lcens.ext  <- M.ext <- N.ext <- L.ext <- NULL else {
      M.ext <- ext * M
      N.ext <- ext * N
      L.ext <- ext * L
      mgrid.ext <- seq(W$xrange[1], W$xrange[2] + (ext - 1) * diff(W$xrange), by = cell.width)
      ngrid.ext <- seq(W$yrange[1], W$yrange[2] + (ext - 1) * diff(W$yrange), by = cell.height)
      lgrid.ext <- seq(W$zrange[1], W$zrange[2] + (ext - 1) * diff(W$zrange), by = cell.length)
      mcens.ext <- (mgrid.ext + 0.5 * cell.width)[-(M.ext + 1)]
      ncens.ext <- (ngrid.ext + 0.5 * cell.height)[-(N.ext + 1)]
      lcens.ext <- (lgrid.ext + 0.5 * cell.length)[-(L.ext + 1)]
    }
  
  return(list(M = M, N = N, L=L, mgrid = mgrid, ngrid = ngrid, lgrid=lgrid, mcens = mcens, ncens = ncens, lcens=lcens,
              cell.width = cell.width, cell.height = cell.height, cell.length = cell.length, M.ext = M.ext, N.ext = N.ext, L.ext = L.ext, 
              mgrid.ext = mgrid.ext, ngrid.ext = ngrid.ext, lgrid.ext = lgrid.ext,mcens.ext = mcens.ext, ncens.ext = ncens.ext,lcens.ext = lcens.ext))
}


fft_ma_3d=function(mygrid,M,N,L,ext=2,sill,range,std)
{
  
  r=function(u,phi_r) (1-(1.5*u/phi_r)+0.5*(u/phi_r)^3)*as.numeric(u<=phi_r) + (as.numeric(u>phi_r))*0
  

  cent <- expand.grid(mygrid$mcens, mygrid$ncens,mygrid$lcens)
  Rx <- mygrid$M.ext * mygrid$cell.width
  Ry <- mygrid$N.ext * mygrid$cell.height
  Rz <- mygrid$L.ext * mygrid$cell.length
  m.abs.diff.row1 <- abs(mygrid$mcens.ext[1] - mygrid$mcens.ext)
  m.diff.row1 <- pmin(m.abs.diff.row1, Rx - m.abs.diff.row1)
  n.abs.diff.row1 <- abs(mygrid$ncens.ext[1] - mygrid$ncens.ext)
  n.diff.row1 <- pmin(n.abs.diff.row1, Ry - n.abs.diff.row1)
  l.abs.diff.row1 <- abs(mygrid$lcens.ext[1] - mygrid$lcens.ext)
  l.diff.row1 <- pmin(l.abs.diff.row1, Rz - l.abs.diff.row1)
  cent.ext.row1 <- expand.grid(m.diff.row1, n.diff.row1, l.diff.row1)
  
  D.ext.row1 <- array(sqrt(cent.ext.row1[, 1]^2 + cent.ext.row1[, 2]^2+ cent.ext.row1[, 3]^2), c(mygrid$M.ext,mygrid$N.ext,mygrid$L.ext))
  SIGMA.Y.ext.row1 <- sill * r(D.ext.row1, range)
  
  
  d <- dim(SIGMA.Y.ext.row1)
  dp <- prod(d)
  sdp <- sqrt(dp)
  prefix <- sqrt(Re(fft(SIGMA.Y.ext.row1, TRUE)))
  realz <- prefix * (fft(array(std, c(d[1], d[2], d[3])))/sdp)
  realz <- as.vector(Re(fft(realz, TRUE)/sdp)[1:M, 1:N, 1:L])
  return(realz)
}



fft_ma_3d_optimal=function(D.ext.row1,M,N,L,sill,range,std)
{
  
  r=function(u,phi_r) (1-(1.5*u/phi_r)+0.5*(u/phi_r)^3)*as.numeric(u<=phi_r) + (as.numeric(u>phi_r))*0
  
  SIGMA.Y.ext.row1 <- sill * r(D.ext.row1, range)
  
  
  d <- dim(SIGMA.Y.ext.row1)
  dp <- prod(d)
  sdp <- sqrt(dp)
  prefix <- sqrt(Re(fft(SIGMA.Y.ext.row1, TRUE)))
  realz <- prefix * (fft(array(std, c(d[1], d[2], d[3])))/sdp)
  realz <- as.vector(Re(fft(realz, TRUE)/sdp)[1:M, 1:N, 1:L])
  return(realz)
}

