from mete import *

p = []
Svals = [10, 11, 13, 14, 16, 18, 21, 23, 26, 30, 34, 38, 43, 48, 55, 62, 70,
         78, 89, 100]

Nvals = [120, 186, 289, 447, 694, 1076, 1668, 2587, 4011, 6220, 9646, 14957,
         23193, 35965, 55769, 86479, 134099, 207941, 322444, 500000]
for s in Svals:
    for n in Nvals:
        p.append(exp(-get_beta(s, n)))


def gr_one(x): return x > 1

which(map(gr_one, p))


file = open('param_space_p_vals.csv','w')
file.write('S' + ',' + 'N' + ',' + 'p' + '\n')
i = 0
for s in Svals:
    for n in Nvals:
        file.write(str(s) + ',' + str(n) + ',' + str(p[i]) + '\n')
        i += 1

file.close()
