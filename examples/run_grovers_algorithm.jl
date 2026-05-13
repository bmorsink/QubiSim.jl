using QubiSim

print("\n")
print("1. grovers algorithm: arbitrary number of bits with simple oracle based on controlled_z.\n")

numberOfQubits = 4

oracle = 9 # in between 0 and (2^nr_qubits)-1
binaryOracle = reverse(digits(oracle, base=2, pad=numberOfQubits)')

qcInitializeSuperposition=createQuantumCircuit(numberOfQubits)
for k in 1:numberOfQubits
    hGate!(qcInitializeSuperposition, k)
end

qcOracle=createQuantumCircuit(numberOfQubits)
for k in 1:numberOfQubits
    if binaryOracle[k]==0
        xGate!(qcOracle, k)
    else
        idGate!(qcOracle, k)
    end
end
controlledUGate!(qcOracle, Vector(1:(numberOfQubits-1)), [numberOfQubits], createSingleQubitOperationZ())
for k in 1:numberOfQubits
    if binaryOracle[k]==0
        xGate!(qcOracle, k)
    else
        idGate!(qcOracle, k)
    end
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

#qcGrover=concatenateQuantumCircuits(qcInitializeSuperposition,
#        concatenateQuantumCircuits(qcOracle,
#        concatenateQuantumCircuits(qcAmplification,qcMeasure)))

qcGrover=concatenateQuantumCircuits(qcInitializeSuperposition,
          concatenateQuantumCircuits(qcOracle,
          concatenateQuantumCircuits(qcAmplification,
          concatenateQuantumCircuits(qcOracle,
          concatenateQuantumCircuits(qcAmplification,qcMeasure)))))

qpGrover=compileQuantumCircuit(qcGrover, optimizeNumberOfSteps=true)
iqs=createInitialQubitState(vector, repeat([1. 0.], outer=[numberOfQubits,1])) # |high=control,...,low=controlled>=|00...0>
qo=runQuantumProgram(qpGrover, iqs, 1)

probeMeasureProbability(qo, 3, "Probabilities of grovers algorithm with arbitrary bits with simple oracle (marked at " * string(oracle) * ") based on controlled_z")


print("\n")
print("2. grovers algorithm: arbitrary number of bits with simple oracle based on controlled_z, and 2 marked states\n")

numberOfQubits = 4

oracle1 = 7 # in between 0 and (2^nr_qubits)-1
binaryOracle1 = reverse(digits(oracle1, base=2, pad=numberOfQubits)')

oracle2 = 14 # in between 0 and (2^nr_qubits)-1
binaryOracle2 = reverse(digits(oracle2, base=2, pad=numberOfQubits)')

qcInitializeSuperposition=createQuantumCircuit(numberOfQubits)
for k in 1:numberOfQubits
    hGate!(qcInitializeSuperposition, k)
end

qcOracle1=createQuantumCircuit(numberOfQubits)
for k in 1:numberOfQubits
    if binaryOracle1[k]==0
        xGate!(qcOracle1, k)
    else
        idGate!(qcOracle1, k)
    end
end
controlledUGate!(qcOracle1, Vector(1:(numberOfQubits-1)), [numberOfQubits], createSingleQubitOperationZ())
for k in 1:numberOfQubits
    if binaryOracle1[k]==0
        xGate!(qcOracle1, k)
    else
        idGate!(qcOracle1, k)
    end
end

qcOracle2=createQuantumCircuit(numberOfQubits)
for k in 1:numberOfQubits
    if binaryOracle2[k]==0
        xGate!(qcOracle2, k)
    else
        idGate!(qcOracle2, k)
    end
end
controlledUGate!(qcOracle2, Vector(1:(numberOfQubits-1)), [numberOfQubits], createSingleQubitOperationZ())
for k in 1:numberOfQubits
    if binaryOracle2[k]==0
        xGate!(qcOracle2, k)
    else
        idGate!(qcOracle2, k)
    end
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

#qcGrover=concatenateQuantumCircuits(qcInitializeSuperposition,
#        concatenateQuantumCircuits(qcOracle1,
#        concatenateQuantumCircuits(qcOracle2,
#        concatenateQuantumCircuits(qcAmplification,qcMeasure))))

qcGrover=concatenateQuantumCircuits(qcInitializeSuperposition,
          concatenateQuantumCircuits(qcOracle1,
          concatenateQuantumCircuits(qcOracle2,
          concatenateQuantumCircuits(qcAmplification,
          concatenateQuantumCircuits(qcOracle1,
          concatenateQuantumCircuits(qcOracle2,
          concatenateQuantumCircuits(qcAmplification,qcMeasure)))))))

qpGrover=compileQuantumCircuit(qcGrover, optimizeNumberOfSteps=true)
iqs=createInitialQubitState(vector, repeat([1. 0.], outer=[numberOfQubits,1])) # |high=control,...,low=controlled>=|00...0>
qo=runQuantumProgram(qpGrover, iqs, 1)

probeMeasureProbability(qo, 3, "Probabilities of grovers algorithm with arbitrary bits with simple oracle (marked at " * string(oracle1) * " and " * string(oracle2) * ") based on controlled_z")


print("\n")
print("3. grovers algorithm: arbitrary number of bits with unitary oracle based on controlled_x and ancilla bit.\n")

numberOfQubits = 4

oracle = 9 # in between 0 and (2^nr_qubits)-1
binaryOracle = reverse(digits(oracle, base=2, pad=numberOfQubits)')

qcInitializeSuperposition=createQuantumCircuit(numberOfQubits+1)
for k in 1:numberOfQubits
    hGate!(qcInitializeSuperposition, k)
end
hGate!(qcInitializeSuperposition, numberOfQubits+1) # ancillary bit also in superposition

qcOracle=createQuantumCircuit(numberOfQubits+1)
for k in 1:numberOfQubits
    if binaryOracle[k]==0
        xGate!(qcOracle, k)
    else
        idGate!(qcOracle, k)
    end
end
controlledUGate!(qcOracle, Vector(1:(numberOfQubits)), [numberOfQubits+1], createSingleQubitOperationX())
for k in 1:numberOfQubits
    if binaryOracle[k]==0
        xGate!(qcOracle, k)
    else
        idGate!(qcOracle, k)
    end
end

qcAmplification=createQuantumCircuit(numberOfQubits+1)
for k in 1:numberOfQubits
    hGate!(qcAmplification, k)
    xGate!(qcAmplification, k)
end
controlledUGate!(qcAmplification, Vector(1:(numberOfQubits-1)), [numberOfQubits], createSingleQubitOperationZ())
for k in 1:numberOfQubits
    xGate!(qcAmplification, k)
    hGate!(qcAmplification, k)
end

qcMeasure=createQuantumCircuit(numberOfQubits+1)
measureGate!(qcMeasure, Vector(1:numberOfQubits))

#qcGrover=concatenateQuantumCircuits(qcInitializeSuperposition,
#        concatenateQuantumCircuits(qcOracle,
#        concatenateQuantumCircuits(qcAmplification,qcMeasure)))

qcGrover=concatenateQuantumCircuits(qcInitializeSuperposition,
          concatenateQuantumCircuits(qcOracle,
          concatenateQuantumCircuits(qcAmplification,
          concatenateQuantumCircuits(qcOracle,
          concatenateQuantumCircuits(qcAmplification,qcMeasure)))))

qpGrover=compileQuantumCircuit(qcGrover, optimizeNumberOfSteps=true)
iqs=createInitialQubitState(vector, [repeat([1. 0.], outer=[numberOfQubits,1]);0. 1.]) # |high=control,...,low=controlled>=|00...01>
qo=runQuantumProgram(qpGrover, iqs, 10)

probeMeasureOutcome(qo, 3, "Probabilities of grovers algorithm with arbitrary number of bits with unitary oracle (marked at " * string(oracle) * ") based on controlled_x and ancilla bit")


print("\n")
print("4. solving 2x2 sudoku using grovers algorithm.\n")

numberOfQubits = 9

qcInitializeSuperposition=createQuantumCircuit(numberOfQubits)
hGate!(qcInitializeSuperposition, 1)
hGate!(qcInitializeSuperposition, 2)
hGate!(qcInitializeSuperposition, 3)
hGate!(qcInitializeSuperposition, 4)
hGate!(qcInitializeSuperposition, 9)

qcOracle=createQuantumCircuit(numberOfQubits)
cnotGate!(qcOracle, 1,5)
cnotGate!(qcOracle, 2,5)
cnotGate!(qcOracle, 1,6)
cnotGate!(qcOracle, 3,6)
cnotGate!(qcOracle, 2,7)
cnotGate!(qcOracle, 4,7)
cnotGate!(qcOracle, 3,8)
cnotGate!(qcOracle, 4,8)
controlledUGate!(qcOracle, Vector(5:8), [9], createSingleQubitOperationX())
cnotGate!(qcOracle, 1,5)
cnotGate!(qcOracle, 2,5)
cnotGate!(qcOracle, 1,6)
cnotGate!(qcOracle, 3,6)
cnotGate!(qcOracle, 2,7)
cnotGate!(qcOracle, 4,7)
cnotGate!(qcOracle, 3,8)
cnotGate!(qcOracle, 4,8)

qcAmplification=createQuantumCircuit(numberOfQubits)
for k in 1:4
    hGate!(qcAmplification, k)
    xGate!(qcAmplification, k)
end
controlledUGate!(qcAmplification, Vector(1:3), [4], createSingleQubitOperationZ())
for k in 1:4
    xGate!(qcAmplification, k)
    hGate!(qcAmplification, k)
end

qcMeasure=createQuantumCircuit(numberOfQubits)
measureGate!(qcMeasure, Vector(1:4))

#qcGrover=concatenateQuantumCircuits(qcInitializeSuperposition,
#        concatenateQuantumCircuits(qcOracle,
#        concatenateQuantumCircuits(qcAmplification,qcMeasure)))

qcGrover=concatenateQuantumCircuits(qcInitializeSuperposition,
          concatenateQuantumCircuits(qcOracle,
          concatenateQuantumCircuits(qcAmplification,
          concatenateQuantumCircuits(qcOracle,
          concatenateQuantumCircuits(qcAmplification,qcMeasure)))))

qpGrover=compileQuantumCircuit(qcGrover, optimizeNumberOfSteps=true)
iqs=createInitialQubitState(vector, [repeat([1. 0.], outer=[8,1]);0. 1.]) # |high=control,...,low=controlled>=|00...01>
qo=runQuantumProgram(qpGrover, iqs, 1)

probeMeasureProbability(qo, 3, "Probabilities of solving 2x2 sudoku using grovers algorithm")
