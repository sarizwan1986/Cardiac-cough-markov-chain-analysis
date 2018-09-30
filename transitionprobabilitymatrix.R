nstate = 3
N = matrix(c(18,13,2,8,63,21,2,26,50),nstate, nstate, byrow=T)
P = N/rowSums(N) 
P