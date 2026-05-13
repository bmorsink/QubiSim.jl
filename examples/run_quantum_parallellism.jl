using QubiSim

print("1. Quantum parallelism: |0yx> --> |(x+y)yx>\n")
print("\n")
print("Quantum computer can calculate many different instances simultaneously.\n")
print("|0yx> --> |(x+y)yx>\n")
print("|000> --> |000>\n")
print("|001> --> |101>\n")
print("|010> --> |110>\n")
print("|011> --> |011>\n")
print("We want to compute all these four cases at once.\n")
print("Bring qubits 1(y) and 2(z) in superposition. Assume zeroBasedNumbering.\n")

qc = createQuantumCircuit(3; zeroBasedNumbering=true)

hGate!(qc, 0)
hGate!(qc, 1)
cnotGate!(qc, 0, 1)
cnotGate!(qc, 1, 2)
cnotGate!(qc, 0, 1)
barrier!(qc)
measureGate!(qc, [0, 1, 2])

qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

iqs=createInitialQubitState(vector, [1. 0.; 1. 0.; 1. 0.])
qo=runQuantumProgram(qp, iqs, 1024)
probeMeasureOutcome(qo, 2, "Quantum parallelism: |0yx> --> |(x+y)yx> with both x & y in superposition")