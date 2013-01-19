
conlisk = function(psi,l,r) ((psi * l) + (1-psi)) / (psi*(l + r) + 2*(1 - psi))
harte = function(the,l,r) ((the * l) + 1) / (the*(l+r) + 2)

l = 1:10
r = 0

## HEAP special case
plot(l,conlisk(.5,l,r),type='o',ylim=c(0,1))
points(l,harte(1,l,r),col='red',pch=19)

## what about more generally 
## does theta = psi/(1-psi) ?

psi = .25
plot(l,conlisk(psi,l,r),type='o',ylim=c(0,1))
points(l,harte(psi/(1-psi),l,r),col='red',pch=19)

psi = .8
plot(l,conlisk(psi,l,r),type='o',ylim=c(0,1))
points(l,harte(psi/(1-psi),l,r),col='red',pch=19)


## YES !
