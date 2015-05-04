A <- matrix(rbinom(10*20,1,0.5),10,20)
B <- colSums(A)
C <- factor(B,levels=seq(1,10,1))
table(C)
