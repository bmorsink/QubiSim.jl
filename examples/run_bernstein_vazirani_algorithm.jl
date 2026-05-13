using QubiSim
using LinearAlgebra, PlotlyJS

print("1. Bernstein Vazirani algorithm: fixed 3 bits\n")
print("\n")

qc=createQuantumCircuit(4)
hGate!(qc, 1)
hGate!(qc, 2)
hGate!(qc, 3)
hGate!(qc, 4)

# 110: 1>CNOT(1,4), 1>CNOT(2,4), 0>I(3)
cnotGate!(qc, 1, 4)
cnotGate!(qc, 2, 4)

hGate!(qc, 1)
hGate!(qc, 2)
hGate!(qc, 3)
hGate!(qc, 4)

measureGate!(qc, [1, 2, 3])

qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
iqs=createInitialQubitState(vector, [1. 0.;1. 0.;1. 0.;0. 1.]) # |high=control,...,low=controlled>=|00...01>
qo=runQuantumProgram(qp, iqs, 1)

probeMeasureOutcome(qo, 3, "1. Oracle = 110")


print("2. Bernstein Vazirani algorithm: variable n=5 bits\n")
print("\n")

numberOfQubits = 5
bitString = [1,1,0,1,1]

qc=createQuantumCircuit(numberOfQubits+1)
for k in 1:(numberOfQubits+1)
    hGate!(qc, k)
end

for k in 1:numberOfQubits
    if bitString[k] == 1
        cnotGate!(qc, k, numberOfQubits+1)
    else
        idGate!(qc, k)
    end
end

for k in 1:(numberOfQubits+1)
    hGate!(qc, k)
end

measureGate!(qc, Vector(1:numberOfQubits))

qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
iqs=createInitialQubitState(vector, [repeat([1. 0.], outer=[numberOfQubits,1]);0. 1.]) # |high=control,...,low=controlled>=|00...01>
qo=runQuantumProgram(qp, iqs, 1)

probeMeasureOutcome(qo, 3, "2. Oracle = "*join(string.(bitString)))
