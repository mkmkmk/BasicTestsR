# ----------------------
# https://gigadom.in/2016/06/23/introducing-qcsimulator-a-5-qubit-quantum-computing-simulator-in-r/
# ----------------------

if (F)
    install.packages("QCSimulator")

library(QCSimulator)

rm(list=ls())
init()
ls()


Hadamard(q0_)

a = SGate(TGate(Hadamard(TGate(Hadamard(q0_)))))

measurement(a)
res = measurement(a)


a = TensorProd(Hadamard(I2),I2)
b = DotProduct(a,q00_)
measurement(b)

t(Mod(b)^2)


if (F)
    plotMeasurement(measurement(b))


# -----------------
a = TensorProd(Hadamard(I2),I2)
b = CNOT2_01(a)
c = Hadamard(TGate(Hadamard(SGate(I2))))
d = TensorProd(I2,c)
e = DotProduct(b,d)
f = DotProduct(e,q00_)
measurement(f)




detach("package:QCSimulator", unload=TRUE)
