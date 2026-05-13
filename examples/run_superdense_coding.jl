using QubiSim
using LinearAlgebra, PlotlyJS

print("1. Superdense Coding: Transmit two classical bits using one qubit\n")
print("\n")

print("Quantum teleportation is a process by which the state of qubit can be\n")
print("transmitted from one location to another, using two bits of classical\n")
print("communication and a Bell pair. In other words, we can say it is a\n")
print("protocol that destroys the quantum state of a qubit in one location\n")
print("and recreates it on a qubit at a distant location, with the help of\n")
print("shared entanglement. Superdense coding is a procedure that allows\n")
print("someone to send two classical bits to another party using just a single\n")
print("qubit of communication. The teleportation protocol can be thought of as\n")
print("a flipped version of the superdense coding protocol, in the sense that\n")
print("Alice and Bob merely swap their equipment.\n")
print("\n")
print("Teleportation: Transmit one qubit using two classical bits\n")
print("Superdense Coding: Transmit two classical bits using one qubit\n")
print("\n")

code=2

qc=createQuantumCircuit(2)
# create entangled pair for Alice and Bob
hGate!(qc, 1)
cnotGate!(qc, 1, 2)
# Alice encodes the first qubit (00=I, 01=X, 10=Z, 11=ZX)
if code==0
    idGate!(qc, 1)
elseif code==1
    xGate!(qc, 1)
elseif code==2
    zGate!(qc, 1)
elseif code==3
    zGate!(qc, 1)
    xGate!(qc, 1)
else
    error("wrong code.")
end
# Bob decodes Alice's qubit and measures the qubits which now contain Alice's code
cnotGate!(qc, 1, 2)
hGate!(qc, 1)
measureGate!(qc, [1, 2])

qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
iqs=createInitialQubitState(vector, [1. 0.;1. 0.])
qo=runQuantumProgram(qp, iqs, 1)

probeMeasureOutcome(qo, 3, "1. Decoding 2-bits code="*string(code)*" over superdense coded 1-qubit")


print("2. Superdense Coding: Transmit two classical bits chosen at random using one qubit\n")
print("\n")

qc=createQuantumCircuit(3)
# create entangled pair for Alice and Bob
hGate!(qc, 2)
cnotGate!(qc, 2, 3)
# Alice encodes the first qubit at random (00=I, 01=X, 10=Z, 11=ZX)
hGate!(qc,1)
measureGate!(qc,[1]) # random generator used in next controlledUGate
controlledUGate!(qc, [1], [2], createSingleQubitOperationZ())
hGate!(qc,1)
measureGate!(qc,[1]) # random generator used in next controlledUGate
controlledUGate!(qc, [1], [2], createSingleQubitOperationX())
# Bob decodes Alice's qubit and measures the qubits which now contain Alice's code
cnotGate!(qc, 2, 3)
hGate!(qc, 2)
measureGate!(qc, [2, 3])

qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
iqs=createInitialQubitState(vector, [1. 0.;1. 0.;1. 0.])
qo=runQuantumProgram(qp, iqs, 1000)

probeMeasureOutcome(qo, 7, "2. Decoding random 2-bits codes over superdense coded 1-qubit")
