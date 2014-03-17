source('./scripts/spat_analytical_prob_funcs.R')

load_heap('./scripts')

A = 1
A0 = 4
n0 = 10
psi = .5

pSingle = sapply(0:10, function(n) single_prob(n, n0, psi, c=4))
pBisect = sapply(0:10, function(n) bisect_prob(n, A, n0, A0, psi))
pQuad = sapply(0:10, function(n) quad_prob(n, A, n0, A0, psi))
pMETE = c(0.2625532468917048, 0.19645684221237666, 0.14699986120521316, 0.1099934161162577, 0.08230314974266902, 0.06158376288999527, 0.046080373151563055, 0.034479880574693735, 0.02579975124191928, 0.019304799003087254, 0.014444916970520172)

plot(0:10, pSingle, type='l', ylim=range(pSingle, pBisect, pQuad, pMETE),
     ylab='Probability', xlab='Cell Abundance')
lines(0:10, pBisect, type='l', col='dodgerblue')
lines(0:10, pQuad, type='l', col='green3')
lines(0:10, pMETE, type='l', col='red')
legend('topright', c('single', 'bisection', 'quadsection', 'METE'),
       col=c('black', 'dodgerblue', 'green3', 'red'), lty=1,bty='n')

##
A0 = 16

lwd = 3
pSingle = sapply(0:10, function(n) single_prob(n, n0, psi, c=4))
pBisect = sapply(0:10, function(n) bisect_prob(n, A, n0, A0, psi))
pQuad = sapply(0:10, function(n) quad_prob(n, A, n0, A0, psi))
pMETE = c(0.6152875860366914, 0.236719120670997, 0.09107276558625134, 0.03503835519504151, 0.01348027949816897, 0.005186257583645668, 0.0019953048990991374, 0.000767652122201461, 0.00029533821171213786, 0.00011362524348552834, 4.371495270557634e-05)

plot(0:10, pSingle, type='l', ylim=range(pSingle, pBisect, pQuad, pMETE),
     ylab='Probability', xlab='Cell Abundance', lwd=lwd)
lines(0:10, pBisect, type='l', col='dodgerblue', lwd=lwd)
lines(0:10, pQuad, type='l', col='green3', lwd=lwd)
lines(0:10, pMETE, type='l', col='red', lwd=lwd)
legend('topright', c('single', 'bisection', 'quadsection', 'METE'),
       col=c('black', 'dodgerblue', 'green3', 'red'), lty=1,bty='n')
