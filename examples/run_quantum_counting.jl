using QubiSim
using LinearAlgebra, PlotlyJS

print("\n")
print("1. Quantum counting: first test oracle with 5 marked states and verify oracle behavior classically.\n")

numberOfQubits = 4

qcOracle=createQuantumCircuit(numberOfQubits)
hGate!(qcOracle, 3)
hGate!(qcOracle, 4)
controlledUGate!(qcOracle, [1,2], [3], createSingleQubitOperationX())
hGate!(qcOracle, 3)
xGate!(qcOracle, 3)
controlledUGate!(qcOracle, [1,3], [4], createSingleQubitOperationX())
xGate!(qcOracle, 3)
hGate!(qcOracle, 4)
xGate!(qcOracle, 2)
xGate!(qcOracle, 4)
hGate!(qcOracle, 3)
controlledUGate!(qcOracle, [1,2,4], [3], createSingleQubitOperationX())
xGate!(qcOracle, 2)
xGate!(qcOracle, 4)
hGate!(qcOracle, 3)

qpOracle=compileQuantumCircuit(qcOracle; optimizeNumberOfSteps=true)

oracleResponse=zeros(ComplexF64,16,16)
for oracle in 0:(2^numberOfQubits)-1
    local binaryOracle = reverse(digits(oracle, base=2, pad=numberOfQubits)')
    value=zeros(numberOfQubits,2)
    for k in 1:numberOfQubits
        if binaryOracle[k]==1
            value[k,2]=1.
        else
            value[k,1]=1.
        end
    end
    local iqs=createInitialQubitState(vector, value)
    local qo=runQuantumProgram(qpOracle, iqs, 1)
    oracleResponse[:,oracle+1]=qo.output[2].q
end

trace=scatter(y=angle.(diag(oracleResponse)), mode="lines", name="oracle response")
layout=Layout(title_text="Oracle response given qubits input", xaxis_title_text="qubits input", yaxis_title_text="angle(diag(oracleResponse))")
p1=plot([trace], layout)
display(p1)


print("\n")
print("2. Quantum counting: first test oracle with 5 marked states and find them with grovers algorithm.\n")

numberOfQubits = 4

qcOracle=createQuantumCircuit(numberOfQubits)
hGate!(qcOracle, 3)
hGate!(qcOracle, 4)
controlledUGate!(qcOracle, [1,2], [3], createSingleQubitOperationX())
hGate!(qcOracle, 3)
xGate!(qcOracle, 3)
controlledUGate!(qcOracle, [1,3], [4], createSingleQubitOperationX())
xGate!(qcOracle, 3)
hGate!(qcOracle, 4)
xGate!(qcOracle, 2)
xGate!(qcOracle, 4)
hGate!(qcOracle, 3)
controlledUGate!(qcOracle, [1,2,4], [3], createSingleQubitOperationX())
xGate!(qcOracle, 2)
xGate!(qcOracle, 4)
hGate!(qcOracle, 3)

qcInitializeSuperposition=createQuantumCircuit(numberOfQubits)
for k in 1:numberOfQubits
    hGate!(qcInitializeSuperposition, k)
end

qcAmplification=createQuantumCircuit(numberOfQubits)
for k in 1:numberOfQubits
    hGate!(qcAmplification, k)
    xGate!(qcAmplification, k)
end
controlledUGate!(qcAmplification, Vector(1:(numberOfQubits-1)), [numberOfQubits], createSingleQubitOperationZ())
for k in 1:numberOfQubits
    xGate!(qcAmplification, k)
    hGate!(qcAmplification, k)
end

qcMeasure=createQuantumCircuit(numberOfQubits)
measureGate!(qcMeasure, Vector(1:numberOfQubits))

qcGrover=concatenateQuantumCircuits(qcInitializeSuperposition,
        concatenateQuantumCircuits(qcOracle,
        concatenateQuantumCircuits(qcAmplification,qcMeasure)))

qpGrover=compileQuantumCircuit(qcGrover; optimizeNumberOfSteps=true)
iqs=createInitialQubitState(vector, repeat([1. 0.], outer=[numberOfQubits,1])) # |high=control,...,low=controlled>=|00...0>
qo=runQuantumProgram(qpGrover, iqs, 1)

probeMeasureProbability(qo, 3, "Probabilities of grovers algorithm on test oracle with 5 marked states")


print("\n")
print("3. Quantum counting: quantum counting algorithm on oracle with 5 marked states.\n")

numberOfQubits = 4

qcOracle=createQuantumCircuit(numberOfQubits)
hGate!(qcOracle, 3)
hGate!(qcOracle, 4)
controlledUGate!(qcOracle, [1,2], [3], createSingleQubitOperationX())
hGate!(qcOracle, 3)
xGate!(qcOracle, 3)
controlledUGate!(qcOracle, [1,3], [4], createSingleQubitOperationX())
xGate!(qcOracle, 3)
hGate!(qcOracle, 4)
xGate!(qcOracle, 2)
xGate!(qcOracle, 4)
hGate!(qcOracle, 3)
controlledUGate!(qcOracle, [1,2,4], [3], createSingleQubitOperationX())
xGate!(qcOracle, 2)
xGate!(qcOracle, 4)
hGate!(qcOracle, 3)

qcAmplification=createQuantumCircuit(numberOfQubits)
for k in 1:numberOfQubits
    hGate!(qcAmplification, k)
    xGate!(qcAmplification, k)
end
controlledUGate!(qcAmplification, Vector(1:(numberOfQubits-1)), [numberOfQubits], createSingleQubitOperationZ())
for k in 1:numberOfQubits
    xGate!(qcAmplification, k)
    hGate!(qcAmplification, k)
end

qcGroverIteration=concatenateQuantumCircuits(qcOracle,qcAmplification)
qpGroverIteration=compileQuantumCircuit(qcGroverIteration; optimizeNumberOfSteps=true)

qcControlledGroverIteration=createQuantumCircuit(numberOfQubits+1)
controlledUGate!(qcControlledGroverIteration, [1], [2,3,4,5], qpGroverIteration.program[1])
qpControlledGroverIteration=compileQuantumCircuit(qcControlledGroverIteration; optimizeNumberOfSteps=true)
U=qpControlledGroverIteration.program[1]

numberOfCountingQubits = 4
numberOfSearchingQubits = 4
numberOfQubits = numberOfCountingQubits + numberOfSearchingQubits

qcQuantumCounting=createQuantumCircuit(numberOfQubits)
for k in 1:numberOfQubits
    hGate!(qcQuantumCounting, k)
end
for k in 1:numberOfCountingQubits
    for m in 1:(2^(k-1))
        unitaryUGate!(qcQuantumCounting, [5-k; numberOfCountingQubits .+ (1:numberOfSearchingQubits)], U)
    end
end
iqftGate!(qcQuantumCounting, Vector(1:4))
measureGate!(qcQuantumCounting, Vector(1:numberOfCountingQubits))
qpQuantumCounting=compileQuantumCircuit(qcQuantumCounting; optimizeNumberOfSteps=true)
iqs=createInitialQubitState(vector, repeat([1. 0.], outer=[numberOfQubits,1])) # |high=control,...,low=controlled>=|00...0>
qo=runQuantumProgram(qpQuantumCounting, iqs, 1)

probeMeasureProbability(qo, 3, "Probabilities of quantum counting algorithm on test oracle with 5 marked states")

measuredInt=5 # or 11 from plot, FIX THIS

theta=(measuredInt/(2^numberOfCountingQubits))*2*pi
N=2^numberOfSearchingQubits
M=N*sin(theta/2)^2
m=numberOfCountingQubits-1
err=(sqrt(2*M*N) + N/(2. ^(m+1)))*(2. ^(-m))
print("Estimated number of solutions = ",round(N-M, digits=3), ", error < ",round(err, digits=3))
