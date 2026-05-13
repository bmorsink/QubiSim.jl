using QubiSim
using LinearAlgebra, PlotlyJS


print("\n")
print("1. Create entangled state of 2 bits, take partial trace on 2nd bit to obtain reduced density operator and evaluate purity.\n")

# rho_ab=(|00>+|11>)ab(<00|+<11|)ab
# rho_a=|0>a<0|a+|1>a<1|a
qc=createQuantumCircuit(2)
hGate!(qc, 1)
cnotGate!(qc, 1, 2)
qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
iqs=createInitialQubitState(density, [[1.00,[1. 0.;1. 0.]],[0.00,[0. 1.;0. 1.]]])
qo=runQuantumProgram(qp, iqs, 1)
rho=partialTrace(qo.output[2].rho, traceIndex=2)
(entropy, purity, purity_norm, purity_str)=calculateEntropyAndPurity(rho)
print("entropy=",entropy,", purity=",purity,", normalized purity=",purity_norm,", state identified as: ",purity_str,"\n")

print("\n")
print("2. Create non-entangled state of 2 bits, take partial trace on 2nd bit to obtain reduced density operator and evaluate purity.\n")

# rho_ab=(|0>+|1>)a*(|0>+|1>)b(<0|+<1|)a*(<0|+<1|)b
# rho_a=(|0>+|1>)a(<0|+<1|)a
qc=createQuantumCircuit(2)
hGate!(qc, 1)
hGate!(qc, 2)
qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
iqs=createInitialQubitState(density, [[1.00,[1. 0.;1. 0.]],[0.00,[0. 1.;0. 1.]]])
qo=runQuantumProgram(qp, iqs, 1)
rho=partialTrace(qo.output[2].rho, traceIndex=2)
(entropy, purity, purity_norm, purity_str)=calculateEntropyAndPurity(rho)
print("entropy=",entropy,", purity=",purity,", normalized purity=",purity_norm,", state identified as: ",purity_str,"\n")
