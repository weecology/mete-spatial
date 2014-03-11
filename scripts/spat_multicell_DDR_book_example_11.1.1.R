## Investigate multi-cell probabilies for different spatial
## configurations described in Harte 2011 11.1.1

source('./scripts/spat_functions.R')
source('./scripts/spat_analytical_prob_funcs.R')

load_heap('./scripts')

## possible spatial arragments for n0 = 2 and A0 = 4
## with a bisection scheme

2 0   0 0    0 0   0 2
0 0   2 0    0 2   0 0
abu1 = c(2, 0, 0, 0, 2, 0)

1 0    0 1  
1 0    0 1
abu2 = c(1, 1, 0, 0, 2, 0)


1 0   1 1    0 0   0 1
0 1   0 0    1 1   1 0
abu3 = c(1, 0, 1, 0, 1, 1)

grains = c(rep(1, 4), rep(2, 2))

## 1/9
multi_prob(get_bisection_history(grains, abu1), psi=0.5)

## 1/9
multi_prob(get_bisection_history(grains, abu2), psi=0.5)

## 1/12
multi_prob(get_bisection_history(grains, abu3), psi=0.5)

2 * multi_prob(get_bisection_history(grains, abu1), psi=0.5) + 
4 * multi_prob(get_bisection_history(grains, abu2), psi=0.5) + 
4 * multi_prob(get_bisection_history(grains, abu3), psi=0.5)

## number of shared species on average for distance-based classification
(4 * (1/9) * 0) + (4 * (1/9) * (1/6)) + (2 * (1/18) * (1/6))

## number of shared species on average for bisection-based classification
(4 * (1/9) * 0) + (2 * (1/9) * (1/6)) + (4 * (1/12) * (1/6))


## so they have an identical number of shared species overall 
## they must differ in how those are spatially likely
n = c(2, 4)
avgS = (4 * .25 + 6 * .5) / 10
sor = sor_heap(1, 2, 4, 'golden', sor_use_c=T, heap_use_c=T)$Sor
avgSor = sum(sor * n)/sum(n)
sor * avgS
avgSor * avgS ## average shared species

## Shared species as a function of distance or seperation order
## distance-based
## dist = 1 and sqrt(2)
sh_db_r = c((4 * (1/9) * (1/6)) ,
            (2 * (1/18) * (1/6)))

sh_db_nr = c((4 * .1161 * (1/6)),
             (2 * .1011 * (1/6)))
  
## bisection-based  
## seperation order, j = 2 distance = 1 and j = 1 distance = sqrt(2)
sh_bb = c((1/2) * (2/9), 
          (1/4) * (4/12))

## assume dimenstions of square are 1 x 1
d_db = c(1, sqrt(2))
d_bb = c(1, mean(c(1, 1, sqrt(2), sqrt(2))))

png('./figs/distance_vs_bisection_based_ddr.png')
plot(d_db, sh_db_nr, type='o', xlab='Distance', ylab='# shared species',
     ylim=range(c(sh_db_r, sh_db_nr, sh_bb)), xlim=range(c(d_db, d_bb)))
points(d_bb, sh_bb, type='o', col='red', pch=19)
points(d_db, sh_db_r, type='o', pch=19)
lgd = c('Bisection-based, Recursive', 'Distance-based, Recursive',
        'Distance-based, Non-recursive')
legend('topright', lgd,
       col=c('red', 'black', 'black'), pch=c(19, 19, 1), lty=1, bty='n')
dev.off()

dist(sh_db) / dist(d_db)
dist(sh_bb) / dist(d_bb)








