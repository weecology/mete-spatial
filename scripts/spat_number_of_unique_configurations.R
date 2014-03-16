source('./scripts/spat_analytical_prob_funcs.R')
source('./scripts/spat_functions.R')

A = 1
A0 = 2^(1:7)

n0 = 2^(1:7)

inputs = expand.grid(n0, A0)

g = mapply(calc_g, inputs[,1], inputs[,2])

out = data.frame(n0=inputs[,1], A0=inputs[,2], g)

AVN = 6.0221413e+23

col = terrain.colors(length(A0) + 3)

png('./figs/number_of_configurations.png', width=480*2, height=480)
ylab = 'Number of Unique Configurations of\n Indistinguishable Individuals'
par(mfrow=c(1,2))
plot(g ~ A0, data=out, log='xy', type='n', ylab='')
i = 1
for(n in n0) {
  lines(g ~ A0, data=out, subset=n0 == n, col=col[i], lwd=3)
  i = i + 1
}
abline(h=AVN)
mtext(ylab, 2, padj=-1.25)
legend('topleft', paste('n0 =', n0), col=col, lty=1, lwd=3, bty='n')
##
plot(g ~ n0, data=out, log='xy', type='n', ylab='')
i = 1
for(A in A0) {
  lines(g ~ n0, data=out, subset=A0 == A, col=col[i], lwd=3)
  i = i + 1
}
abline(h=AVN)
mtext(ylab, 2, padj=-1.25)
legend('topleft', paste('A0 =', A0), col=col, lty=1, lwd=5, bty='n')
dev.off()