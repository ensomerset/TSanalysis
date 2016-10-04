#Bartlett M-Test Function 

require('multitaper')

spec.block <- function(blocksize, data, nw, k,overlap)
{
N <- length(data)
incr <- ceiling(blocksize*(1-overlap))
sectInd <- seq(1,N-blocksize+1, by=incr)
numSect <- length(sectInd)
nFFT <- 2^(floor(log2(N)+4))
Spec <- matrix(NA, ncol=numSect, nrow = (nFFT/2+1))

for (i in 1:numSect)
{bloc <- data[sectInd[i]:(sectInd[i]+blocksize -1)]
Spec[,i] <- spec.mtm(bloc, nw=nw, k=k,deltat=1, nFFT=nFFT, plot=F)$spec}

freq <- spec.mtm(bloc, nw=nw, k=k, deltat=1, nFFT=nFFT, plot=F)$freq
S <- cbind(Spec, freq)
return(S)}

BMT <- function(blocksize, data, nw, k, overlap)
{
N <- length(data)
nu <- 2*k
S <- spec.block(blocksize, data, nw, k, overlap)
freq.ind <- dim(S)[2]
numblocks <- dim(S)[2]-1
index <- 1:(dim(S)[1])
M <- as.vector(NA)
for (i in index)
{
M[i] <- numblocks*nu*log((1/numblocks)*sum(S[i,1:numblocks])) -nu*sum(log(S[i,1:numblocks]))}
M <- M/(1+((numblocks+1)/(3*nu*numblocks)))
M[1] <- M[1]/2
M[length(M)] <- M[length(M)]/2
plot(S[,freq.ind], M, type="l", xlab="Frequency [Cycles/Month]", ylab="M(f)")
expectation <- numblocks-1
chi <- qchisq(0.999, df= expectation)
out <- list(M, expectation,chi)
return(out)}



