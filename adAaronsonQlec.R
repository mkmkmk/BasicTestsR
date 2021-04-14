# ----------------------
# ----------     ad qlec
# ----------     https://www.scottaaronson.com/qclec.pdf
# ----------------------
rm(list=ls())

# ------ ad (2.8)
{
    r = runif(1)
    print(r)
    matrix(c(1/2,1/2,1/2,1/2), nrow=2) %*% t(t(c(r,1-r)))
}


# ------ ad fig. 3.2, error
U = matrix(c(1,1,-1,1), nrow=2) / sqrt(2)
U

qp = t(t(c(1,1))) / sqrt(2)
qp  # |+〉

round(U %*% qp, 10)

q0 = t(t(1:0))  # |0〉
q0
U %*% q0

q1 = t(t(0:1))  # |1〉
q1
U %*% q1

round(U %*% (U %*% q0), 10)

# ----- AD lecture 4 

# eq. (4.2) Hadamard
H = matrix(c(1,1,1,-1), nrow=2) / sqrt(2)
H
q0 = t(t(1:0))  # |0〉
q1 = t(t(0:1))  # |1〉

H %*% q0
round(H %*% (H %*% q0), 10)

H %*% q1
round(H %*% (H %*% q1), 10)
H %*% (H %*% (H %*% q1))


# --------- ad bomb (trochę w przód wybiegło)
q00 = kronecker(q0, q0)
q01 = kronecker(q0, q1)
q10 = kronecker(q1, q0)
q11 = kronecker(q1, q1)
not = (1 - diag(1, 2))
cnot = kronecker (q0 %*% t(q0), diag(1,2)) + kronecker (q1 %*% t(q1), not)
rotgate = function(th) matrix(c(cos(th), sin(th), -sin(th), cos(th)), 2, 2)
N = 10
step1_dud = kronecker(rotgate(pi/2/N) %*% q0, q0)
step1_bomb = cnot %*% step1_dud
round(step1_bomb, 3)
round(step1_dud, 3)
all(step1_dud == cnot %*% cnot %*% step1_dud) # $$$

# after step 1
#
# |C〉 - control bit, |P〉 - probe bit
#
#     | before meas          | after meas, if P=|0〉 eq. (5.1)
#     |----------------------|---------------------------------
#     | |CP〉| bomb  | dud   |   | bomb         | dud           
#     |------|-------|-------|---|--------------|-------------
#     | |00〉| 0.988 | 0.988 | * | 0.988/0.988  | 0.988/1     
#     | |01〉| 0.000 | 0.000 |   | 0            | 0         
#     | |10〉| 0.000 | 0.156 | * | 0.000/0.988  | 0.156/1 
#     | |11〉| 0.156 | 0.000 |   | 0            | 0         

# post partial measurement state
# measBasis - measurement basis-state list
postmeas = function(state, ...)
{
    measBasis = list(...)
    if(is.list(measBasis[[1]]))
        measBasis = unlist(measBasis, recursive = FALSE)
    res = 0
    for (basis in measBasis)
        res = res + (Conj(t(basis)) %*% state)[1,1] * basis
    res / sqrt(sum(Mod(res)^2)) # eq. (5.1)
}
round(postmeas(step1_bomb, q00, q10), 3)
round(postmeas(step1_dud, q00, q10), 3)
round(postmeas(step1_bomb, list(q00, q10)), 3)

# prawie symulacja bomby
bomb = cnot %*% kronecker(rotgate(pi/2/N), diag(1,2))
dud = kronecker(rotgate(pi/2/N), diag(1,2))
qcirc = dud
qcirc = bomb
{
    ctrl = q0
    probe = q0
    state = kronecker(ctrl, probe)
    for(i in 1:N)
    {
        state = qcirc %*% state
        # zakładamy że bomba nie wybuchła:
        state = postmeas(state, q00, q10)
    }
    round(state, 3)
}


# ----- AD lecture 5 

# eq. (5.5)
not = (1 - diag(1, 2))
noti = kronecker(not, diag(1, 2))
noti

q00 = kronecker(q0, q0)
q01 = kronecker(q0, q1)
q10 = kronecker(q1, q0)
q11 = kronecker(q1, q1)

q00
q01
q10
q11

noti %*% q00   # => q10
noti %*% q01   # => q11
noti %*% q10   # => q00
noti %*% q11   # => q01


# eq. (5.6)
H = (1 - 2 * (0:1) %*% t(0:1)) / 2^.5
H
H = matrix(c(rep(1, 3), -1), 2, 2) / 2^.5
H
H %*% 1:0
H %*% 0:1

H2 = kronecker(H, H)
H2
H2 %*% q00   # |++〉
H2 %*% q01   # |+-〉
H2 %*% q10   # |-+〉
H2 %*% q11   # |--〉


sum((H2 %*% q00)^2)
sum((H2 %*% q01)^2)
sum((H2 %*% q10)^2)
sum((H2 %*% q11)^2)

round(H2 %*% H2 %*% q00, 10)
round(H2 %*% H2 %*% q01, 10)
round(H2 %*% H2 %*% q10, 10)
round(H2 %*% H2 %*% q11, 10)


# ad order
kronecker(H, diag(1, 2)) * 2^.5
kronecker(diag(1, 2), H) * 2^.5

# eq. (5.7) (5.8) 
cnot = kronecker (q0 %*% t(q0), diag(1,2)) + kronecker (q1 %*% t(q1), not)
cnot
not
I = diag(1, 2)
HI = kronecker(H, I)
HI

# (H ⊗ I) (q0 ⊗ q0) == (H q0) ⊗ q0
kronecker(H, I) %*% kronecker(q0, q0) == kronecker(H %*% q0, q0)
HI %*% q00 == kronecker(H %*% q0, q0)

cnot %*% kronecker(H %*% q0 , q0) # eq. (5.7)
cnot %*% HI %*% q00

# ok, można niezależnie wpuścić sumę na CNOT, i wtedy zgadnąć wynik
cnot %*% HI %*% q00
cnot %*% kronecker(qp, q0)
cnot %*% (q00 + q10) / sqrt(2)
(cnot %*% q00 + cnot %*% q10) / sqrt(2)
(q00 + q11) / sqrt(2)

bell = cnot %*% HI %*% q00
bell

H %*% q0   #  |+〉
(q00 + q10) / 2^.5
kronecker(H %*% q0, q0)
(q00 + q11) / 2^.5

HI %*% bell

# If Alice sees |0〉
# wg. eq. (5.1)
alfa = .5 # == beta
scale = (alfa^2+alfa^2)^.5
alfa/scale
1/2^.5

# Bell pair, if (after meas) one of pair is |0〉then ...
postmeas(bell, q00, q10) # ... both are |0〉
postmeas(bell, q00, q01) # ... both are |0〉
# Bell pair, if one of pair is |1〉then ...
postmeas(bell, q11, q01) # ... both are |1〉
postmeas(bell, q10, q11) # ... both are |1〉

ifOneIs = q0 # then ...
postmeas(bell, kronecker(ifOneIs, q0), kronecker(ifOneIs, q1))
postmeas(bell, kronecker(q0, ifOneIs), kronecker(q1, ifOneIs))

ifOneIs = q1 # then ...
postmeas(bell, kronecker(q0, ifOneIs), kronecker(q1, ifOneIs))
postmeas(bell, kronecker(ifOneIs, q0), kronecker(ifOneIs, q1))


# ----------- lib
multikron = function(...)
{
    inp = list(...)
    if(is.list(inp[[1]]))
        inp = unlist(inp, recursive = FALSE)
    res = inp[[1]]
    if (length(inp) < 2)
        return(inp[[1]])
    for (ii in 2:length(inp))
        res = kronecker(res, inp[[ii]])
    res
}
IH = kronecker(I, H)
stopifnot(all(multikron(I) == I))
stopifnot(all(multikron(I,H) == IH))
stopifnot(all(multikron(H,I,I) == kronecker(HI, I)))
stopifnot(all(multikron(I, H, I, I) == kronecker(I, kronecker(HI, I))))



# ----- swap 
#
#     00 01 10 11    
# 00   1  0  0  0 
# 01   0  0  1  0    
# 10   0  1  0  0
# 11   0  0  0  1
swap = matrix(c(1,0,0,0, 0,0,1,0, 0,1,0,0, 0,0,0,1), 4)
swap %*% q00 == q00
swap %*% q10 == q01 
swap %*% q01 == q10
swap %*% q11 == q11


# ----- AD lecture 6

# mixed states eq. (6.2)

# outer product with itself
q0 %*% t(Conj(q0))
q1 %*% t(Conj(q1))

# mixed state eq. (6.3), ok!
(q0 %*% t(Conj(q0)) + q1 %*% t(Conj(q1))) / 2 


bell %*%  t(Conj(bell)) 

# eq. (6.4), ok!
qp = H %*% q0  # |+〉
qm = H %*% q1  # |-〉

(qp %*% t(Conj(qp)) + qm %*% t(Conj(qm))) / 2 

# ad order
qp %*% t(Conj(qp))
t(Conj(qp)) %*% qp  

# ad p. 6.1.2 〈+|M|+〉= 19/2
M = matrix(c(1/2,-10, -10,1/2 ), 2, 2)
t(Conj(qp)) %*% M %*% qp  #〈+|M|+〉
#〈+|M|+〉== (H ρ H^†)[1,1]
HMH = round(H %*% M %*% t(Conj(H)), 10)
HMH
t(Conj(qp)) %*% M %*% qp  #〈+|M|+〉
t(Conj(qm)) %*% M %*% qm  #〈-|M|-〉
t(Conj(qm)) %*% M %*% qp  #〈-|M|+〉
t(Conj(qp)) %*% M %*% qm  #〈+|M|-〉

# ad rank
qr(q1 %*% t(Conj(q1)))$rank
qr((q0 %*% t(Conj(q0)) + q1 %*% t(Conj(q1))) / 2 )$rank
qr(qp %*% t(Conj(qp)))$rank

# eq. (6.9) --> ok ale brak ⊗
(q00 + q01 + q10) / sqrt(3)
kronecker(q0, qp) * sqrt(2/3) + kronecker(q1, q0) * sqrt(1/3)  # ok!
eq69 = (q00 + q01 + q10) / sqrt(3)

# If Alice sees |0〉Bob sees 
bobA0 = sqrt(2/3) * qp
bobA0
# If Alice sees |1〉Bob sees 
bobA1 = sqrt(1/3) * q0
bobA1

# density mx
bobA0 %*% t(Conj(bobA0)) + bobA1 %*% t(Conj(bobA1))

# ad eq (6.10) tracing out - ok
# (ρB)jj′=∑i αij (αij′)*
aij = matrix(eq69, 2, 2)
rhoB = matrix(0, 2, 2)
for(j in 1+0:1)
    for(jp in 1+0:1)
        rhoB[j, jp] = aij[1, j] * Conj(aij[1, jp]) + aij[2, j] * Conj(aij[2, jp])
rhoB
# αij * αij^†
aij %*% t(Conj(aij))

# ----- AD lecture 7

# --------- ad Bloch sphere
if (FALSE)
{
    plot(exp(1i*pi/2), xlim = c(-1,1), ylim = c(-1,1))
    
    th = 0:99/100*pi
    plot(th, cos(th/2), type = 'l', col = "blue")
    lines(th, sin(th/2), type = 'l', col = "red")
}

toBloch = function (qbit)
{
    stopifnot(sum(Mod(qbit)^2) < 1 + 1e-10)
    # TODO if not, rotate |0〉
    stopifnot(abs(Im(qbit[1])) < 1e-10)
    th = 2 * acos(qbit[1])
    #sin = (1 - qbit[1]^2)^.5
    #stopifnot(abs(sin - sin(th/2)) < 1e-9)
    #phi = if (abs(th) > 1e-9) acos(Re(qbit[2]) / sin) else 0
    phi = Arg(qbit[2])
    round(c(th, phi) / pi, 8) 
}

toStateTb = function(th, phi)
    round(q0 %*% t(cos(th/2)) + q1 %*% t(exp(1i*phi)*sin(th/2)), 10)

toState = function(bloch)
    toStateTb(bloch[1],  bloch[2])

th  = c(0, pi, pi/2, pi/2, pi/2, pi/2)
phi = c(0, 0,  0,    pi,   pi/2, -pi/2) 

bloch = toStateTb(th, phi)
bloch
all(Mod(bloch[,1] - q0) < 1e-10)
all(Mod(bloch[,2] - q1) < 1e-10)
all(Mod(bloch[,3] - qp) < 1e-10)
all(Mod(bloch[,4] - qm) < 1e-10)

toBloch(q0)
toBloch(q1)
toBloch(qp)
toBloch(qm)
toBloch((q1 - q0) / sqrt(2))

pauliZ = matrix(c(1, 0, 0, -1), 2, 2)
pauliY = matrix(c(0, 1i, -1i, 0), 2, 2)
pauliX = not
sgate = matrix(c(1, 0, 0, 1i), 2, 2)
tgate = matrix(c(1, 0, 0, exp(1i*pi/4)), 2, 2)

all(Mod(sgate - (tgate%*%tgate)) < 1e-9)

toBloch(pauliZ %*% qp)
toBloch(qm)
all(Mod(toState(toBloch(pauliZ %*% qp) * pi) - qm) < 1e-9)

toBloch(pauliX %*% q0) == toBloch(q1)

toBloch(sgate %*% qp)
toBloch(sgate %*% qm)
sgate %*% sgate %*% qm == qp

toState(toBloch(sgate %*% qm))


# p.7.2 cloning
kronecker(.33*q0+.44*q1, q0)

cnot %*% kronecker(qp, q0)


# ----- AD lecture 8

# quantum money bomb attack
# ad eq. (8.1)
rotgate = function(th) matrix(c(cos(th), sin(th), -sin(th), cos(th)), 2, 2)

round(rotgate(pi/2), 10)
round(rotgate(pi/4), 10)
round(rotgate(0), 10)
N=10

kronecker(rotgate(pi/2/N) %*% q0, q0)

cnot %*% kronecker(rotgate(pi/2/N) %*% q0, q0)  # ok

#  |Cα〉
#  |00〉
#  |01〉
#  |10〉
#  |11〉
postmeas(cnot %*% kronecker(rotgate(pi/2/N) %*% q0, q0), q00, q10) # |00〉-> ok
postmeas(cnot %*% kronecker(rotgate(pi/2/N) %*% q0, q1), q01, q11) # |01〉-> ok

                       #  |Cα
qpp = kronecker(qp,qp) #  |++〉
qpm = kronecker(qp,qm) #  |+-〉
qmp = kronecker(qm,qp) #  |-+〉
qmm = kronecker(qm,qm) #  |--〉

postmeas(cnot %*% kronecker(rotgate(pi/2/N) %*% q0, qm), qpm, qmm)

# -----------------
qcirc = cnot %*% kronecker(rotgate(pi/2/N), diag(1,2))
hbase = kronecker(diag(1,2), H)
H
# --
alfa = q0
meas = list(q00, q10)
meas = list(kronecker(q0, alfa), kronecker(q1, alfa))
# --
alfa = q1
meas = list(q01, q11)
meas = list(kronecker(q0, alfa), kronecker(q1, alfa))
# --
# dla |+〉i |-〉obracam przed pomiarem bazę żeby pomiar był w bazie { |+〉, |-〉}
# wtedy |+〉<--> |0〉i |-〉<--> |1〉
qcirc = hbase %*% cnot %*% kronecker(rotgate(pi/2/N), diag(1,2))
# --
#alfa = H %*% qm
alfa = qm
#meas = list(qpm, qmm)
meas = list(q01, q11)
meas = list(kronecker(q0, H %*% alfa), kronecker(q1, H %*% alfa))
round(meas[[2]], 3)
# ok - jest nawet opisane to małe skakanie !!!
# --
#alfa = H %*% qp
alfa = qp
#meas = list(qpp, qmp)
meas = list(q00, q10)
meas = list(kronecker(q0, H %*% alfa), kronecker(q1, H %*% alfa))
# --
# pomiar |-〉 (obraca się)
-cnot
cminusnot = kronecker (q0 %*% t(q0), diag(1,2)) + kronecker (q1 %*% t(q1), -not)
qcirc = hbase %*% cminusnot %*% kronecker(rotgate(pi/2/N), diag(1,2))
alfa = qm
meas = list(q01, q11)
{
    ctrl = q0
    state = kronecker(ctrl, alfa)
    for(i in 1:N)
    {
        state = postmeas(qcirc %*% state, meas)
        # state = meas %*% state / sqrt(Conj(t(state)) %*% meas %*% state)[1,1]
        round(state, 3)
    }
    round(state, 3)
}
sum(Mod(state)^2)


# -----------------
# -- bomb - swaped bits - ok
#  |αC〉
#  |00〉
#  |01〉
#  |10〉
#  |11〉
# error: swcnot = kronecker (q0 %*% t(q0), not) + kronecker (q1 %*% t(q1),  diag(1,2))
swcnot = swap %*% cnot %*% swap
all(swcnot == kronecker (diag(1,2), q0 %*% t(q0)) + kronecker (not, q1 %*% t(q1)) )
qcirc = swcnot %*% kronecker(diag(1,2), rotgate(pi/2/N))
# --
alfa = q0
meas = list(q00, q01)
meas = list(kronecker(alfa, q0), kronecker(alfa, q1))
# --
alfa = q1
meas = list(q10, q11)
# --
{
    ctrl = q0
    state = kronecker(alfa, ctrl)
    for(i in 1:N)
        state = postmeas(qcirc %*% state, meas)
    state
    all(state == kronecker(alfa, q0))
}
sum(Mod(state)^2)

    
#------------------------
# Nimish Mishra, Understanding the basics of measurements in Quantum Computation
# https://towardsdatascience.com/understanding-basics-of-measurements-in-quantum-computation-4c885879eba0
M0 = q0 %*% t(Conj(q0))

M00 = q00 %*% t(Conj(q00))
M01 = q01 %*% t(Conj(q01))
M10 = q10 %*% t(Conj(q10))
M11 = q11 %*% t(Conj(q11))

state = step1_bomb 
state = step1_dud 
meas = M00 + M10
meas %*% state / sqrt(Conj(t(state)) %*% meas %*% state)[1,1]

postmeas2 = function(state, measBasis)
{
    meas = 0
    for (it in measBasis)
        meas = meas + it %*% t(Conj(it))
    meas %*% state / sqrt(Conj(t(state)) %*% meas %*% state)[1,1]
}

round(postmeas2(step1_bomb, list(q00, q10)), 3)
round(postmeas2(step1_dud, list(q00, q10)), 3)

qmm %*% t(Conj(qmm)) + qpm %*% t(Conj(qpm)) + qmp %*% t(Conj(qmp)) + qpp %*% t(Conj(qpp)) 
state = cnot %*% kronecker(rotgate(pi/2/N) %*% q0, qm)
postmeas(state , list(qmm, qpm))
meas = qmm %*% t(Conj(qmm)) +  qpm %*% t(Conj(qpm))
meas %*% state / sqrt(Conj(t(state)) %*% meas %*% state)[1,1]

postmeas(state , list(qmm, qpm))


# ad flip does nothing to |+〉
not %*% qp == qp
cnot %*% qpp == qpp 

# CNOT*CNOT == I, no tak --> crypto: c ⊕ k = m ⊕ k ⊕ k = m
all(cnot %*% cnot == diag(1,4)) 


# ----- AD lecture 9 
# superdense coding

# eq. (9.1), ok chyba trzeba zamienić kolejność w ostatnim 

XI = kronecker(not, diag(1,2)) 
ZI = kronecker(pauliZ, diag(1,2))
I = diag(1,2)
II = kronecker(I, I)

# wszędzie dla 2go bit-a jest I czyli nie ruszamy go
II %*% bell
XI %*% bell
ZI %*% bell
(ZI %*% XI) %*% bell

I %*% q1 

XI %*% ZI %*% bell # error


# eq. (9.2) - error!!
IH = kronecker(diag(1,2),H)
errBobT = matrix(c(1,1,0,0, 0,0,-1,1, 0,0,1,-1, 1,-1,0,0), 4) / sqrt(2)
bobT = swap %*% HI %*% cnot %*% swap 

bobT == errBobT

sqrt(2) * errBobT
sqrt(2) * bobT
sqrt(2) * IH %*% swap %*% cnot %*% swap  
sqrt(2) * IH %*% swcnot

#decode - ok
round(bobT %*% bell, 9)
round(bobT %*% ZI %*% bell, 9)
round(bobT %*% XI %*% bell, 9)
round(bobT %*% XI %*% ZI %*% bell, 9)

if (FALSE)
{
    round(errBobT %*% XI %*% bell, 9)
    round(errBobT %*% ZI %*% bell, 9)
    round(errBobT %*% ZI %*% XI %*% bell, 9)
}


# ----- AD lecture 10

# quantum teleportation

# eq. 10.2
alfa = .666
beta = (1-alfa^2)^.5
eq102 = kronecker(alfa*q0 + beta*q1, bell) 
eq102 * sqrt(2) # ok

# eq. (10.3)
eq103 = kronecker(cnot, I) %*% eq102 # ok
eq103 * sqrt(2) # ok

# eq. (10.4)
eq104 =  kronecker(HI, I) %*% eq103 
2 * eq104 # ok!

# Table 10.1
# if Alice sees:
ifAliceSees = q00
postm = postmeas(eq104, list(kronecker(ifAliceSees, q0), kronecker(ifAliceSees, q1)))
postm

ifAliceSees = q01
postm = postmeas(eq104, list(kronecker(ifAliceSees, q0), kronecker(ifAliceSees, q1)))
postm
postm = kronecker(II, pauliX) %*% postm
postm
all(Mod(postm - multikron(ifAliceSees, alfa*q0 + beta*q1)) < 1e-9)


ifAliceSees = q10
postm = postmeas(eq104, list(kronecker(ifAliceSees, q0), kronecker(ifAliceSees, q1)))
postm
kronecker(II, pauliZ) %*% postm

ifAliceSees = q11
postm = postmeas(eq104, list(kronecker(ifAliceSees, q0), kronecker(ifAliceSees, q1)))
postm
kronecker(II, pauliZ) %*% kronecker(II, pauliX) %*% postm
# ok!!!!

# -------
# teleportacja połówki innej pary Bell
# bit 0 to dodatkowy bit bity 1 2 3 to tj. Fig. 10.1
state = multikron(I, H, I, I) %*% multikron(I, cnot, I) %*% kronecker(bell, bell) 
ifAliceSees = q00
postm = postmeas(state, multikron(q0, ifAliceSees, q0), multikron(q0, ifAliceSees, q1), multikron(q1, ifAliceSees, q0), multikron(q1, ifAliceSees, q1))
postm
all(postm == multikron(swap, swap) %*% multikron(q0, bell, q0))
all(multikron(I, swap, I) %*% multikron(swap, I, I) %*% postm == multikron(ifAliceSees, bell))
## --
ifAliceSees = q01
postm = postmeas(state, multikron(q0, ifAliceSees, q0), multikron(q0, ifAliceSees, q1), multikron(q1, ifAliceSees, q0), multikron(q1, ifAliceSees, q1))
postm
postm = multikron(I,I,I, pauliX) %*% postm
postm
all(postm == multikron(swap, swap) %*% multikron(q0, bell, q1))
all(multikron(I, swap, I) %*% multikron(swap, I, I) %*% postm == multikron(ifAliceSees, bell))
## --
ifAliceSees = q10
postm = postmeas(state, multikron(q0, ifAliceSees, q0), multikron(q0, ifAliceSees, q1), multikron(q1, ifAliceSees, q0), multikron(q1, ifAliceSees, q1))
postm
postm = multikron(I,I,I, pauliZ) %*% postm
postm
all(postm == multikron(swap, swap) %*% multikron(q1, bell, q0))
all(multikron(I, swap, I) %*% multikron(swap, I, I) %*% postm == multikron(ifAliceSees, bell))
## --
ifAliceSees = q11
postm = postmeas(state, multikron(q0, ifAliceSees, q0), multikron(q0, ifAliceSees, q1), multikron(q1, ifAliceSees, q0), multikron(q1, ifAliceSees, q1))
postm
postm = multikron(I,I,I, pauliZ) %*% multikron(I,I,I, pauliX) %*% postm
postm
all(postm == multikron(swap, swap) %*% multikron(q1, bell, q1))
all(multikron(I, swap, I) %*% multikron(swap, I, I) %*% postm == multikron(ifAliceSees, bell))

#   ---x--- a
#      |
# a ---x---
#                  
# b ---x---
#      |
#   ---x--- b 

# a ---x-------
#      |          
#   ---x---x---
#          |
#   -------x--- a
#              
# b ----------- b

# -----------------
# ad GHZ
GHZ = (multikron(q0, q0, q0) + multikron(q1, q1, q1)) / sqrt(2)
GHZ

ifISee = q0
postmeas(GHZ, multikron(q0, q0, ifISee), multikron(q0, q1, ifISee), multikron(q1, q0, ifISee), multikron(q1, q1, ifISee))

ifISee = q1
postmeas(GHZ, multikron(q0, q0, ifISee), multikron(q0, q1, ifISee), multikron(q1, q0, ifISee), multikron(q1, q1, ifISee))


someEn = (multikron(q1, q0, q0) + multikron(q0, q1, q0) + multikron(q0, q0, q1)) / sqrt(3)
ifISee = q0
postmeas(someEn, multikron(q0, q0, ifISee), multikron(q0, q1, ifISee), multikron(q1, q0, ifISee), multikron(q1, q1, ifISee))
ifISee = q1
postmeas(someEn, multikron(q0, q0, ifISee), multikron(q0, q1, ifISee), multikron(q1, q0, ifISee), multikron(q1, q1, ifISee))



# ----- AD lecture 14
requal = function(a, b, dig) round(a, dig) == round(b, dig)



ifOneIs = q1 # then ...
postmeas(bell, kronecker(q0, ifOneIs), kronecker(q1, ifOneIs))

ifOneIs = q0 # then ...
postmeas(bell, kronecker(q0, ifOneIs), kronecker(q1, ifOneIs))

Mod(Conj(t(q00)) %*% bell)^2
Mod(Conj(t(q11)) %*% bell)^2

all(requal(rotgate(pi * 1/8) %*% q1, rotgate(pi * 5/8) %*% q0, 9))


# ad eq. (14.3)
Mod(Conj(t(q0)) %*% rotgate(pi/8) %*% q0)^2 == cos(pi/8)^2
Mod(Conj(t(q0)) %*% rotgate(pi/8) %*% q1)^2 == sin(pi/8)^2
# ad eq. (14.4)
Mod(Conj(t(q1)) %*% rotgate(pi/8) %*% q0)^2 == sin(pi/8)^2
Mod(Conj(t(q1)) %*% rotgate(pi/8) %*% q1)^2 == cos(pi/8)^2



# x=0 y=0 -> xy==0
Mod(Conj(t(q0)) %*% rotgate(pi/8) %*% q0)^2

meas = q0
meas = q1
Mod(Conj(t(q0)) %*% meas)^2
Mod(Conj(t(q0)) %*% rotgate(pi/8) %*% meas)^2


round(H2 %*% bell, 9)
round((qpp + qmm) / sqrt(2), 9)
round(H2 %*% bell, 9) == round((qpp + qmm) / sqrt(2), 9)

bell = cnot %*% HI %*% q00

swap %*% kronecker(H %*% q0, q0)
kronecker(H %*% q0, q0)

kronecker(H, rotgate(pi/8)) %*% bell 

# a + b = xy

# x=0 y=0  
# if x=0 Alice measures in {|0〉,|1〉}
# if y=0 Bob measures in {|π/8〉,|5π/8〉}
# sukces gdy a=0 b=0 or a=1 b=1 ---> Table 13.1
# P(a=0)P(b=0|a=0) + P(a=1)P(b=1|a=1)
# P(a=0) == P. Alice sees |0〉
pa0 = as.numeric(Mod(Conj(t(q00)) %*% bell)^2)
pa0
# P(b=0|a=0) == P. Bob sees |π/8〉when Alice sees |0〉
pb0 = Mod(Conj(t(q0)) %*% rotgate(pi/8) %*% q0)^2     # (14.3)
pb0
# P(a=1) == P. Alice sees |1〉
pa1 = as.numeric(Mod(Conj(t(q11)) %*% bell)^2)
pa1
# P(b=1|a=1) == P. Bob sees |5π/8〉when Alice sees |1〉
pb1 = as.numeric(Mod(Conj(t(q1)) %*% rotgate(pi * 5/8) %*% q0)^2)      # (14.4)
pb1

pa0*pb0 + pa1*pb1
cos(pi/8)^2

# x=1 y=0  
# sukces gdy a=0 b=0 or a=1 b=1

# if x=1 Alice measures in {|+〉,|−〉}
# if y=0 Bob measures in {|π/8〉,|5π/8〉}

# P(a=0) == P. Alice sees |+〉
pa0 = as.numeric(Mod(Conj(t(q0)) %*% qp)^2)
pa0 
# P(b=0|a=0) == P. Bob sees |π/8〉when Alice sees |+〉
pb0 = Mod(Conj(t(qp)) %*% rotgate(pi/8) %*% q0)^2    
pb0
# P(a=1) == P. Alice sees |-〉
pa1 = as.numeric(Mod(Conj(t(q1)) %*% qm)^2)
pa1 
# P(b=1|a=1) == P. Bob sees |5π/8〉when Alice sees |+〉
pb1 = as.numeric(Mod(Conj(t(q1)) %*% rotgate(pi * 5/8) %*% q0)^2)     
pb1

Mod(Conj(t(rotgate(pi/8) %*% q0)) %*% qp)^2    

requal((qpp + qmm) / sqrt(2), (q00 + q11) / sqrt(2), 9)
requal(H2 %*% bell, bell, 9)
requal(H2 %*% H2 %*% H2 %*% H2 %*% H2 %*% bell, bell, 9)
requal(H2, H2 %*% H2 %*% H2 %*% H2 %*% H2, 9)

# ---- H --- M
# ---- R --- M

postmeas(bell, q00, q10)

Mod(Conj(t(q00)) %*% postmeas(IH %*% bell, q00, q10))^2
Mod(Conj(t(q00)) %*% postmeas(HI %*% bell, q00, q10))^2



# ----- AD lecture 17
# Deutch
#---- (17.5) 
f = c(0, 1)
(-1)^(1-f)
(-1)^(1) * (-1)^(-f)
-(-1)^(-f)
-(-1)^f
(-1)^(1-f) == -(-1)^f


# ---- (17.9) Deutsch-Josza
multikron(H %*% q0, H %*% q0, H %*% q0, H %*% q0)
multikron(H %*% q0, H %*% q0, H %*% q0, H %*% q1)
multikron(H %*% q0, H %*% q0, H %*% q1, H %*% q0)
multikron(H %*% q0, H %*% q0, H %*% q1, H %*% q1)
multikron(H %*% q0, H %*% q1, H %*% q0, H %*% q0)
multikron(H %*% q0, H %*% q1, H %*% q0, H %*% q1)
multikron(H %*% q0, H %*% q1, H %*% q1, H %*% q0)
multikron(H %*% q0, H %*% q1, H %*% q1, H %*% q1)
multikron(H %*% q1, H %*% q0, H %*% q0, H %*% q0)
multikron(H %*% q1, H %*% q0, H %*% q0, H %*% q1)
multikron(H %*% q1, H %*% q0, H %*% q1, H %*% q0)
multikron(H %*% q1, H %*% q0, H %*% q1, H %*% q1)
multikron(H %*% q1, H %*% q1, H %*% q0, H %*% q0)
multikron(H %*% q1, H %*% q1, H %*% q0, H %*% q1)
multikron(H %*% q1, H %*% q1, H %*% q1, H %*% q0)
multikron(H %*% q1, H %*% q1, H %*% q1, H %*% q1)

# left i right to strony wzoru (17.9)
n = 10
n = 8
n = 6
n = 4
allbits = expand.grid(rep(list(0:1), n))
for(xrow in 1:2^n)
{
    x = unlist(allbits[xrow,], use.names = F)
    left = c()
    for(yrow in 1:2^n)
    {
        y = unlist(allbits[yrow,], use.names = F)
        left = c(left, (-1)^sum(rev(x)*y))
    }
    left = t(t(left)) / sqrt(2^n)
    
    right = list()
    for(i in 1:length(x))
        right[[i]] = if (x[i] == 0) H %*% q0 else H %*% q1
    right = multikron(right)
    stopifnot(max(Mod(left - right)) < 1e-12)
}




























