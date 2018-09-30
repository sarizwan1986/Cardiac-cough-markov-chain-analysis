                          # PART I

nstate = 3
N = matrix(c(59, 37, 8, 38, 244, 67, 7, 91, 149),nstate, nstate, byrow=T)
P = N/rowSums(N) 
                          # PART II 


mpt=function(f) 
 {

                         # PART II A
for(i in 2:k) 
{
Q=Q%*%P 
term2 = 0
Q1 = diag(nstate) 
for(j in (i-1):1){
         Q1 = Q1%*%P
         term2 = term2 + f[j]*Q1[s,s] 
                }
f[i] = Q[r,s] - term2
}
k0 = sum(cumsum(f) <= 0.99999)
f0 = c(f[1:k0], 1-sum(f[1:k0])) 

                             # PART II B

m = 10000
frs = rep(0,m)
for(i in 1:m)
{
  frs[i] = mean(sample(1:(k0+1),1000,prob=f0,replace=TRUE))
          }
frs[i] 
frs.LCL = quantile(frs,0.025) 
frs.UCL = quantile(frs,0.975)
     ppt = list(frs.LCL, mean(frs), frs.UCL)
}

                            # PART III

k = 1000
mpt.mean = matrix(0,nstate, nstate)
mpt.LCL = matrix(0,nstate, nstate)
mpt.UCL = matrix(0,nstate, nstate)
for(r in 1:nstate){
for(s in 1:nstate){
f = rep(0,k+1) 
f[1] = P[r,s] 
Q=P
temp = mpt(f) 
mpt.LCL[r,s] = temp[[1]]

mpt.mean[r,s] = temp[[2]]
mpt.UCL[r,s] = temp[[3]]
}
}
                            ### Printing the output
mpt.LCL 
mpt.mean 
mpt.UCL 


