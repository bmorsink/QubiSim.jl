using QubiSim

print("Quantum teleportation: cloning qubit 1 of Alice with a prepared state to qubit 3 of Bob.\n")
print("Due to no-cloning theorem qubit 1 of Alice will be destroyed along the cloning process.\n")


print("\n")
print("1. Probabilities\n")

alpha=1/sqrt(2)
beta=-1/sqrt(2)

qc=createQuantumCircuit(3)
# create entangled pair for Alice and Bob
hGate!(qc, 2)
cnotGate!(qc, 2, 3)
# Alice performs some steps
cnotGate!(qc, 1, 2)
hGate!(qc, 1)
# Alice performs measurements and transmits results to Bob who then performs some steps
measureGate!(qc, [1,2])
cnotGate!(qc, 2, 3)
controlledUGate!(qc, [1], [3], createSingleQubitOperationZ())
# third qubit of Bob now contains secret code
measureGate!(qc, [3])

qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

iqs=createInitialQubitState(vector, [alpha beta;1. 0.;1. 0.])
qo=runQuantumProgram(qp, iqs, 1024)

probeMeasureOutcome(qo,5,"1. Measure 3rd qubit (Bob's) which should be |->")


print("\n")
print("2. Probe state multibloch vectors\n")

qc=createQuantumCircuit(3)
# create entangled pair for Alice and Bob
hGate!(qc, 2)
cnotGate!(qc, 2, 3)
# Alice performs some steps
cnotGate!(qc, 1, 2)
hGate!(qc, 1)
# Alice performs measurements and transmits results to Bob who then performs some steps
measureGate!(qc, [1,2])
cnotGate!(qc, 2, 3)
controlledUGate!(qc, [1], [3], createSingleQubitOperationZ())

qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

iqs=createInitialQubitState(vector, [alpha beta;1. 0.;1. 0.])
qo=runQuantumProgram(qp, iqs, 1)

probeStateMultiBlochVector(qo, 1, [1], "Alice")
probeStateMultiBlochVector(qo, 4, [3], "Bob")


print("\n")
print("3. Probe state multibloch vectors\n")

theta=pi*rand()
phi=2*pi*rand()            

qc=createQuantumCircuit(3)
# Alice creates her secret state in qubit 0: |psi>=u3|0>
u3Gate!(qc, 1, theta, phi, pi - phi)
barrier!(qc)
# create entangled pair for Alice and Bob
hGate!(qc, 2)
cnotGate!(qc, 2, 3)
# Alice performs some steps
cnotGate!(qc, 1, 2)
hGate!(qc, 1)
# Alice performs measurements and transmits results to Bob who then performs some steps
measureGate!(qc, [1,2])
cnotGate!(qc, 2, 3)
controlledUGate!(qc, [1], [3], createSingleQubitOperationZ())

qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=false)

iqs=createInitialQubitState(vector, [1. 0.;1. 0.;1. 0.])
qo=runQuantumProgram(qp, iqs, 1)

probeStateMultiBlochVector(qo, 2, [1], "Alice")
probeStateMultiBlochVector(qo, 9, [3], "Bob")


print("\n")
print("4. Alice secret state, initialize |+>=H|0>, inverse initialize qubit 2: Bob |+> --> H|+>=|0>\n")

print("Quantum teleportation is designed to send qubits between two parties. We do not have the hardware\n")
print("to demonstrate this, but we can demonstrate that the gates perform the correct transformations on\n")
print("a single quantum chip. Here we use the QASM simulator to simulate how we might test our protocol.\n")
print("On a real quantum computer, we would not be able to sample the statevector, so if we wanted to\n")
print("check our teleportation circuit is working, we need to do things slightly differently.\n")

qc=createQuantumCircuit(3)
# Alice creates her secret state in qubit 1: |psi>=|+>=H|0>
hGate!(qc, 1)
# create entangled pair for Alice and Bob
hGate!(qc, 2)
cnotGate!(qc, 2, 3)
# Alice performs some steps
cnotGate!(qc, 1, 2)
hGate!(qc, 1)
# Alice performs measurements and transmits results to Bob who then performs some steps
measureGate!(qc, [1,2])
cnotGate!(qc, 2, 3)
controlledUGate!(qc, [1], [3], createSingleQubitOperationZ())
# after inverse initialize on Bob's qubit with H|+>=|0>, measure 3rd qubit which should be 100% |0>
hGate!(qc, 3)
measureGate!(qc, [3])

qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

iqs=createInitialQubitState(vector, [1. 0.;1. 0.;1. 0.])
qo=runQuantumProgram(qp, iqs, 1024)

probeMeasureOutcome(qo,5,"4. Measure 3rd qubit (Bob's) which should be 100% |0>")
# probeMeasureProbability(qo, 5, "Measure 3rd qubit (Bob's) which should be 100% |0>")


print("\n")
print("5. Alice secret state, initialize |psi>=u3|0>, inverse initialize qubit 2: Bob |psi> --> u3|psi>=|0>\n")

print("use unitary property: u3*u3=I so initialize and inverse initialize with same u3\n")

theta=rand()
phi=rand()            

qc=createQuantumCircuit(3)
# Alice creates her secret state in qubit 0: |psi>=u3|0>
u3Gate!(qc, 1, theta, phi, pi - phi)
# create entangled pair for Alice and Bob
hGate!(qc, 2)
cnotGate!(qc, 2, 3)
# Alice performs some steps
cnotGate!(qc, 1, 2)
hGate!(qc, 1)
# Alice performs measurements and transmits results to Bob who then performs some steps
measureGate!(qc, [1,2])
cnotGate!(qc, 2, 3)
controlledUGate!(qc, [1], [3], createSingleQubitOperationZ())
# after inverse initialize on Bob's qubit with u3|psi>=|0>, measure 3rd qubit which should be 100% |0>
u3Gate!(qc, 3, theta, phi, pi - phi)
measureGate!(qc, [3])

qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

iqs=createInitialQubitState(vector, [1. 0.;1. 0.;1. 0.])
qo=runQuantumProgram(qp, iqs, 1024)

probeMeasureOutcome(qo,5,"5. Measure 3rd qubit (Bob's) which should be 100% |0>")


print("\n")
print("6. Alice secret state, initialize |psi>=u3|0>, inverse initialize qubit 2: Bob |psi> --> u3|psi>=|0>\n")

print("use unitary property: u3*u3=I so initialize and inverse initialize with same u3\n")

print("The IBM quantum computers currently do not support instructions after measurements,\n")
print("meaning we cannot run the quantum teleportation in its current form on real hardware.\n")
print("Fortunately, this does not limit our ability to perform any computations due to the\n")
print("deferred measurement principle discussed in chapter 4.4 of [1]. The principle states\n")
print("that any measurement can be postponed until the end of the circuit, i.e. we can move\n")
print("all the measurements to the end, and we should see the same results.\n")

print("Any benefits of measuring early are hardware related: If we can measure early, we may\n")
print("be able to reuse qubits, or reduce the amount of time our qubits are in their fragile\n")
print("superposition. In this example, the early measurement in quantum teleportation would\n")
print("have allowed us to transmit a qubit state without a direct quantum communication channel.\n")

print("While moving the gates allows us to demonstrate the teleportation circuit on real\n")
print("hardware, it should be noted that the benefit of the teleportation process (transferring\n") 
print("quantum states via classical channels) is lost.\n")

theta=rand()
phi=rand()            

qc=createQuantumCircuit(3)
# Alice creates her secret state in qubit 0: |psi>=u3|0>
u3Gate!(qc, 1, theta, phi, pi - phi)
# create entangled pair for Alice and Bob
hGate!(qc, 2)
cnotGate!(qc, 2, 3)
# Alice performs some steps
cnotGate!(qc, 1, 2)
hGate!(qc, 1)
# Bob then performs some steps, third qubit of Bob now contains secret code
cnotGate!(qc, 2, 3)
controlledUGate!(qc, [1], [3], createSingleQubitOperationZ())
# after inverse initialize on Bob's qubit with u3|psi>=|0>, measure all qubits, 3rd qubit should be 100% |0>
u3Gate!(qc, 3, theta, phi, pi - phi)
measureGate!(qc, [1,2,3])

qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

iqs=createInitialQubitState(vector, [1. 0.;1. 0.;1. 0.])
qo=runQuantumProgram(qp, iqs, 1024)

probeMeasureOutcome(qo,3,"6. Measure all qubits, 3rd qubit (Bob's) should be 100% |0>")
