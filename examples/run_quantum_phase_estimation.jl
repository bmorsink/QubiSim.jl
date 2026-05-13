using QubiSim
using LinearAlgebra, PlotlyJS, StatsBase


print("\n")
print("1a. Big endian 4 qubit accurate quantum phase estimation of a phase gate\n")

acc=4
phasePart=11/16
theta=2*pi*phasePart
(2^acc)*phasePart
qc=createQuantumCircuit(acc+1)
xGate!(qc,acc+1) # to create eigenstate |1> with eigenphase theta from initial state |0>
qpeGate!(qc, Vector(1:acc), [acc+1], createSingleQubitOperationU1(theta))
measureGate!(qc, Vector(1:acc))
qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
iqs=createInitialQubitState(vector, repeat([1. 0.], outer=[acc+1,1]))
qo=runQuantumProgram(qp, iqs, 1000)
probeMeasureOutcome(qo, 3, "1a. Estimated phase (big endian) of a "*string(round(theta, digits=3))*"[rad] phase gate")

phaseEstimate = 2 * pi * findmax(countmap([shot.measured.outcome for shot in qo.output[3, :]]))[2] / (2^acc)
println("Estimated phase of a $(round(theta, digits=3)) [rad] phase gate: $(round(phaseEstimate, digits=3)) [rad]")
print("\n")


print("\n")
print("1b. Little endian 4 qubit accurate quantum phase estimation of a phase gate\n")

acc=4
phasePart=11/16 # 0.75
theta=2*pi*phasePart
(2^acc)*phasePart
qc=createQuantumCircuit(acc+1, indexType=circuitIndexLittleEndian)
xGate!(qc,acc+1) # to create eigenstate |1> with eigenphase theta from initial state |0>
qpeGate!(qc, Vector(1:acc), [acc+1], createSingleQubitOperationU1(theta))
measureGate!(qc, Vector(1:acc))
qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
iqs=createInitialQubitState(vector, repeat([1. 0.], outer=[acc+1,1]))
qo=runQuantumProgram(qp, iqs, 1000)
probeMeasureOutcome(qo, 3, "1b. Estimated phase (little endian) of a "*string(round(theta, digits=3))*"[rad] phase gate")

phaseEstimate = 2 * pi * findmax(countmap([shot.measured.outcome for shot in qo.output[3, :]]))[2] / (2^acc)
println("Estimated phase of a $(round(theta, digits=3)) [rad] phase gate: $(round(phaseEstimate, digits=3)) [rad]")
print("\n")


print("\n")
print("2a. Big endian 4 qubit accurate quantum phase estimation of a Td-gate\n")

acc=4
qc=createQuantumCircuit(acc+1)
xGate!(qc,acc+1) # to create eigenstate |1> with eigenphase theta from initial state |0>
qpeGate!(qc, Vector(1:acc), [acc+1], createSingleQubitOperationTd())
measureGate!(qc, Vector(1:acc))
qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
iqs=createInitialQubitState(vector, repeat([1. 0.], outer=[acc+1,1]))
qo=runQuantumProgram(qp, iqs, 1000)
probeMeasureOutcome(qo, 3, "2a. Estimated phase (big endian) of a Td-gate ("*string(round(-pi/4, digits=3))*"[rad])")

phaseEstimate = 2 * pi * findmax(countmap([shot.measured.outcome for shot in qo.output[3, :]]))[2] / (2^acc)
println("Estimated phase of a Td-gate ($(round(-pi/4, digits=3)) [rad]): $(round(phaseEstimate, digits=3)) [rad]")
print("\n")


print("\n")
print("2b. Little endian 4 qubit accurate quantum phase estimation of a Td-gate\n")

acc=4
qc=createQuantumCircuit(acc+1, indexType=circuitIndexLittleEndian)
xGate!(qc,acc+1) # to create eigenstate |1> with eigenphase theta from initial state |0>
qpeGate!(qc, Vector(1:acc), [acc+1], createSingleQubitOperationTd())
measureGate!(qc, Vector(1:acc))
qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
iqs=createInitialQubitState(vector, repeat([1. 0.], outer=[acc+1,1]))
qo=runQuantumProgram(qp, iqs, 1000)
probeMeasureOutcome(qo, 3, "2b. Estimated phase (little endian) of a Td-gate ("*string(round(-pi/4, digits=3))*"[rad])")

phaseEstimate = 2 * pi * findmax(countmap([shot.measured.outcome for shot in qo.output[3, :]]))[2] / (2^acc)
println("Estimated phase of a Td-gate ($(round(-pi/4, digits=3)) [rad]): $(round(phaseEstimate, digits=3)) [rad]")
print("\n")


print("\n")
print("3a. Big endian 3 qubit accurate quantum phase estimation of a phase gate\n")

acc=3
theta=.7*2*pi
qc=createQuantumCircuit(acc+1)
xGate!(qc,acc+1) # to create eigenstate |1> with eigenphase theta from initial state |0>
qpeGate!(qc, Vector(1:acc), [acc+1], createSingleQubitOperationU1(theta))
measureGate!(qc, Vector(1:acc))
qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
iqs=createInitialQubitState(vector, repeat([1. 0.], outer=[acc+1,1]))
qo=runQuantumProgram(qp, iqs, 1000)
probeMeasureOutcome(qo, 3, "3a. Estimated phase (big endian) of a "*string(round(theta, digits=3))*"[rad] phase gate")

phaseEstimate = 2 * pi * findmax(countmap([shot.measured.outcome for shot in qo.output[3, :]]))[2] / (2^acc)
println("Estimated phase of a $(round(theta, digits=3)) [rad] phase gate: $(round(phaseEstimate, digits=3)) [rad]")
print("\n")


print("\n")
print("3b. Little endian 3 qubit accurate quantum phase estimation of a phase gate\n")

acc=3
theta=.7*2*pi
qc=createQuantumCircuit(acc+1, indexType=circuitIndexLittleEndian)
xGate!(qc,acc+1) # to create eigenstate |1> with eigenphase theta from initial state |0>
qpeGate!(qc, Vector(1:acc), [acc+1], createSingleQubitOperationU1(theta))
measureGate!(qc, Vector(1:acc))
qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
iqs=createInitialQubitState(vector, repeat([1. 0.], outer=[acc+1,1]))
qo=runQuantumProgram(qp, iqs, 1000)
probeMeasureOutcome(qo, 3, "3b. Estimated phase (little endian) of a "*string(round(theta, digits=3))*"[rad] phase gate")

phaseEstimate = 2 * pi * findmax(countmap([shot.measured.outcome for shot in qo.output[3, :]]))[2] / (2^acc)
println("Estimated phase of a $(round(theta, digits=3)) [rad] phase gate: $(round(phaseEstimate, digits=3)) [rad]")
print("\n")


print("\n")
print("4. Big endian 4 qubit accurate quantum phase estimation of a 1-qubit-controlled phase gate\n")

acc=4
phasePart=12/16
theta=2*pi*phasePart
(2^acc)*phasePart
qc=createQuantumCircuit(acc+2)
xGate!(qc,acc+1) # to create eigenstate |11> with eigenphase theta from initial state |00>
xGate!(qc,acc+2)
qpeGate!(qc, Vector(1:acc), [acc+1,acc+2], createNQubitOperationControlledU(1, 1, createSingleQubitOperationU1(theta)))
measureGate!(qc, Vector(1:acc))
qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
iqs=createInitialQubitState(vector, repeat([1. 0.], outer=[acc+2,1]))
qo=runQuantumProgram(qp, iqs, 1000)
probeMeasureOutcome(qo, 3, "4. Estimated phase (big endian) of a 1-qubit-controlled "*string(round(theta, digits=3))*"[rad] phase gate")

phaseEstimate = 2 * pi * findmax(countmap([shot.measured.outcome for shot in qo.output[3, :]]))[2] / (2^acc)
println("Estimated phase of a 1-qubit-controlled $(round(theta, digits=3)) [rad] phase gate: $(round(phaseEstimate, digits=3)) [rad]")
print("\n")


print("\n")
print("5. Big endian 4 qubit accurate quantum phase estimation of a SWAP gate\n")

acc=4
qc=createQuantumCircuit(acc+2)
# eigvecs(createDoubleQubitOperationSWAP().U)
# 4×4 Matrix{ComplexF64}:
#        0.0+0.0im  1.0+0.0im       0.0+0.0im  0.0+0.0im
#   0.707107+0.0im  0.0+0.0im  0.707107-0.0im  0.0+0.0im
#  -0.707107-0.0im  0.0+0.0im  0.707107+0.0im  0.0+0.0im
#        0.0+0.0im  0.0+0.0im       0.0+0.0im  1.0+0.0im
# angle.(eigvals(createDoubleQubitOperationSWAP().U))
# 4-element Vector{Float64}:
#  3.141592653589793
#  0.0
#  0.0
#  0.0
eigenvec=1 # select eigenstate
if eigenvec==1 # (|01>-|10>)/sqrt(2)
	xGate!(qc,acc+1)
	xGate!(qc,acc+2)
	hGate!(qc,acc+1)
	cnotGate!(qc,acc+1,acc+2)
elseif eigenvec==2 # |00>
elseif eigenvec==3 # (|01>+|10>)/sqrt(2)
	xGate!(qc,acc+2)
	hGate!(qc,acc+1)
	cnotGate!(qc,acc+1,acc+2)
elseif eigenvec==4 # |11>
	xGate!(qc,acc+1)
	xGate!(qc,acc+2)
end
qpeGate!(qc, Vector(1:acc), [acc+1,acc+2], createDoubleQubitOperationSWAP())
measureGate!(qc, Vector(1:acc))
qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
iqs=createInitialQubitState(vector, repeat([1. 0.], outer=[acc+2,1]))
qo=runQuantumProgram(qp, iqs, 1000)
probeMeasureOutcome(qo, 3, "5. Estimated phase (big endian) of a SWAP gate")

phaseEstimate = 2 * pi * findmax(countmap([shot.measured.outcome for shot in qo.output[3, :]]))[2] / (2^acc)
println("Estimated phase of a SWAP gate: $(round(phaseEstimate, digits=3)) [rad]")
print("\n")


print("\n")
print("6. Big endian 4 qubit accurate quantum phase estimation of a phase gate for all eigenstates\n")

acc=4
phasePart=12/16
theta=2*pi*phasePart
(2^acc)*phasePart
qc=createQuantumCircuit(acc+1)
hGate!(qc,acc+1) # to create mix of eigenstates |0> and |1> with eigenphase 0 and theta from initial state |0>
qpeGate!(qc, Vector(1:acc), [acc+1], createSingleQubitOperationU1(theta))
measureGate!(qc, Vector(1:acc))
qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
iqs=createInitialQubitState(vector, repeat([1. 0.], outer=[acc+1,1]))
qo=runQuantumProgram(qp, iqs, 1000)
probeMeasureOutcome(qo, 3, "6. Estimated phase (big endian) of a "*string(round(theta, digits=3))*"[rad] phase gate")

phaseEstimate = 2 * pi * findmax(countmap([shot.measured.outcome for shot in qo.output[3, :]]))[2] / (2^acc)
println("Estimated phase of a $(round(theta, digits=3)) [rad] phase gate: $(round(phaseEstimate, digits=3)) [rad]")
print("\n")


print("\n")
print("7. Big endian 4 qubit accurate quantum phase estimation of a unitary gate with eigenphases 0, pi/2, pi, -pi/2 for all eigenstates\n")

acc=4
qc=createQuantumCircuit(acc+2)
hGate!(qc,acc+1) # to create mix of eigenstates |00>, |01>, |10>  and |11> with eigenphases 0, pi/2, pi, -pi/2 from initial state |00>
hGate!(qc,acc+2) 
qpeGate!(qc, Vector(1:acc), [acc+1, acc+2], UnitaryOperation([1. 0. 0. 0.; 0. -1. 0. 0.; 0. 0. 1.0im 0.; 0. 0. 0. -1.0im]))
measureGate!(qc, Vector(1:acc))
qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
iqs=createInitialQubitState(vector, repeat([1. 0.], outer=[acc+2,1]))
qo=runQuantumProgram(qp, iqs, 1000)
probeMeasureOutcome(qo, 3, "7. Estimated phase (big endian) of a phase gate")

phaseEstimate = 2 * pi * findmax(countmap([shot.measured.outcome for shot in qo.output[3, :]]))[2] / (2^acc)
println("Estimated phase of a unitary gate with eigenphases 0, pi/2, pi, -pi/2 [rad]: $(round(phaseEstimate, digits=3)) [rad]")
print("\n")