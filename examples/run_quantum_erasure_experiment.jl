using QubiSim
using LinearAlgebra, PlotlyJS

print("\n")
print("GHJW theorem on entangled pair |Φ>AB=(|00>+|11>)/sqrt(2) as purification of mixed state ρA=0.5*(|0><0|+|1><1|).\n")
print("\n")

function printEntropyAndPurity(rho)
    (entropy, _, purityNorm, purityString)=calculateEntropyAndPurity(rho)
    print("Entropy=",entropy,", PurityNorm=",purityNorm,", purityString=",purityString,"\n\n")
end


print("\n")
print("1. Quantum erasure experiment: measure and forget qubit B (along σZ), hence decohere qubit A.\n")

qc=createQuantumCircuit(2)
hGate!(qc, 1)
cnotGate!(qc, 1, 2)
measureGate!(qc, [2], [sigmaZ()], forgetOutcome=true)
qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=false)
iqs=createInitialQubitState(density, [[1.00,[1. 0.;1. 0.]],[0.00,[0. 1.;0. 1.]]])
qo=runQuantumProgram(qp, iqs, 1)

print("\n")
print("Density ρAB before ignorant measurement of qubit B=0.5*(|00><00|+|00><11|+|11><00|+|11><11|)\n")
print("ρAB=",round.(real(qo.output[3].rho),digits=2),"\n")
printEntropyAndPurity(qo.output[3].rho)
print("Density ρAB after ignorant measurement of qubit B=0.5*(|00><00|+|11><11|)\n")
print("ρAB=",round.(real(qo.output[4].rho),digits=2),"\n")
printEntropyAndPurity(qo.output[4].rho)
print("Reduced density ρA (partial trace over qubit B) before ignorant measurement of qubit B=0.5*(|0><0|+|1><1|)\n")
print("ρA=",round.(real(partialTrace(qo.output[3].rho, traceIndex=2)),digits=2),"\n")
printEntropyAndPurity(partialTrace(qo.output[4].rho, traceIndex=2))


print("\n")
print("2. Quantum erasure experiment: measure qubit B (along σZ), hence extract which path information of qubit A.\n")

qc=createQuantumCircuit(2)
hGate!(qc, 1)
cnotGate!(qc, 1, 2)
measureGate!(qc, [2], [sigmaZ()])
qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=false)
iqs=createInitialQubitState(density, [[1.00,[1. 0.;1. 0.]],[0.00,[0. 1.;0. 1.]]])
qo=runQuantumProgram(qp, iqs, 1)

print("\n")
print("Density ρAB before measurement of qubit B along σZ=0.5*(|00><00|+|00><11|+|11><00|+|11><11|)\n")
print("ρAB=",round.(real(qo.output[3].rho),digits=2),"\n")
printEntropyAndPurity(qo.output[3].rho)
print("Density ρAB after measurement of qubit B along σZ=|00><00| or |11><11|\n")
print("ρAB=",round.(real(qo.output[4].rho),digits=2),"\n")
printEntropyAndPurity(qo.output[4].rho)
print("Reduced density ρA (partial trace over qubit B) after measurement of qubit B along σZ=|0><0| or |1><1|\n")
print("ρA=",round.(real(partialTrace(qo.output[4].rho, traceIndex=2)),digits=2),"\n")
printEntropyAndPurity(partialTrace(qo.output[4].rho, traceIndex=2))


print("\n")
print("3. Quantum erasure experiment: measure and forget qubit B (along σX), hence decohere qubit A.\n")

qc=createQuantumCircuit(2)
hGate!(qc, 1)
cnotGate!(qc, 1, 2)
measureGate!(qc, [2], [sigmaX()], forgetOutcome=true)
qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=false)
iqs=createInitialQubitState(density, [[1.00,[1. 0.;1. 0.]],[0.00,[0. 1.;0. 1.]]])
qo=runQuantumProgram(qp, iqs, 1)

print("\n")
print("Density ρAB before ignorant measurement of qubit B along σX=0.5*(|00><00|+|00><11|+|11><00|+|11><11|)\n")
print("ρAB=",round.(real(qo.output[3].rho),digits=2),"\n")
printEntropyAndPurity(qo.output[3].rho)
print("Density ρAB after ignorant measurement of qubit B along σX=|++><++|+|--><--|\n")
print("ρAB=",round.(real(qo.output[4].rho),digits=2),"\n")
printEntropyAndPurity(qo.output[4].rho)
print("Reduced density ρA (partial trace over qubit B) after ignorant measurement of qubit B along σX=|+><+|+|-><-|\n")
print("ρA=",round.(real(partialTrace(qo.output[4].rho, traceIndex=2)),digits=2),"\n")
printEntropyAndPurity(partialTrace(qo.output[4].rho, traceIndex=2))


print("\n")
print("4. Quantum erasure experiment: measure qubit B (along σX), hence erase which path information and retrieve interference (prevent decoherence) of qubit A.\n")

qc=createQuantumCircuit(2)
hGate!(qc, 1)
cnotGate!(qc, 1, 2)
measureGate!(qc, [2], [sigmaX()])
qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=false)
iqs=createInitialQubitState(density, [[1.00,[1. 0.;1. 0.]],[0.00,[0. 1.;0. 1.]]])
qo=runQuantumProgram(qp, iqs, 1)

print("\n")
print("Density ρAB before measurement of qubit B along σX=0.5*(|00><00|+|00><11|+|11><00|+|11><11|)\n")
print("ρAB=",round.(real(qo.output[3].rho),digits=2),"\n")
printEntropyAndPurity(qo.output[3].rho)
print("Density ρAB after measurement of qubit B along σX=|++><++| or |--><--|\n")
print("ρAB=",round.(real(qo.output[4].rho),digits=2),"\n")
printEntropyAndPurity(qo.output[4].rho)
print("Reduced density ρA (partial trace over qubit B) after measurement of qubit B along σX=|+><+| or |-><-|\n")
print("ρA=",round.(real(partialTrace(qo.output[4].rho, traceIndex=2)),digits=2),"\n")
printEntropyAndPurity(partialTrace(qo.output[4].rho, traceIndex=2))


print("\n")
print("5. If qubits A and B are not entangled, measuring and forgetting qubit B (along σZ) does not affect qubit A.\n")

qc=createQuantumCircuit(2)
hGate!(qc, 1)
idGate!(qc, 1)
measureGate!(qc, [2], [sigmaZ()], forgetOutcome=true)
# measureGate!(qc, [2], [sigmaZ()])
qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=false)
iqs=createInitialQubitState(density, [[1.00,[1. 0.;1. 0.]],[0.00,[0. 1.;0. 1.]]])
qo=runQuantumProgram(qp, iqs, 1)

print("\n")
print("Density ρAB before measurement of qubit B along σZ=|+0><+0|\n")
print("ρAB=",round.(real(qo.output[3].rho),digits=2),"\n")
printEntropyAndPurity(qo.output[3].rho)
print("Density ρAB after measurement of qubit B along σZ=|+0><+0|\n")
print("ρAB=",round.(real(qo.output[4].rho),digits=2),"\n")
printEntropyAndPurity(qo.output[4].rho)
print("Reduced density ρA (partial trace over qubit B) after measurement of qubit B along σZ=|+><+|\n")
print("ρA=",round.(real(partialTrace(qo.output[4].rho, traceIndex=2)),digits=2),"\n")
printEntropyAndPurity(partialTrace(qo.output[4].rho, traceIndex=2))
