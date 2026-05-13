using QubiSim
using LinearAlgebra, PlotlyJS


print("1. Entanglement alone cannot transmit information: Alice’s reduced density matrix is independent of Bob’s measurement basis, and the basis can only be inferred after classical communication.\n")
print("\n")

# We create many copies of a two-qubit entangled state. The first qubit of each pair is sent to Alice and 
# the second to Bob. Bob now attempts to send a one-bit message to Alice using these entangled pairs. He 
# does this by measuring all of his qubits either in the z-direction or in the x-direction. By choosing 
# the measurement basis, he prepares Alice’s qubits in either the ensemble {|up,z>,|down,z>} or the ensemble 
# {|up,x>,|down,x>}.

# At first glance, it might seem that Alice could read Bob’s message by measuring all of her qubits along the 
# same direction. However, this does not work. Both ensembles on Alice’s side are described by the same density 
# matrix, so from her local measurements alone Alice cannot distinguish whether Bob measured in the z- or the 
# x-basis.

# However, the situation changes if Bob later communicates his measurement outcomes to Alice through a classical 
# channel, without revealing the direction in which he measured. Alice can then compare Bob’s outcomes with her 
# own measurement results. If the two sets of outcomes always match, this indicates that Alice and Bob measured 
# along the same direction. If they measured along different directions, their results will agree with probability 
# about 1/2 and disagree with probability about 1/2.

# By analyzing these correlations over many entangled pairs, Alice can determine Bob’s measurement basis with 
# increasing confidence. The larger the number of shared entangled states, the more reliably she can infer Bob’s 
# choice.

nrShots=100

print("Both Bob and Alice measure along the z-direction:\n")

qc=createQuantumCircuit(2)
hGate!(qc, 1)
cnotGate!(qc, 1, 2) # entangled
measureGate!(qc, [2], [sigmaZ()]) # Bob tries to send a 0=sigmaZ
measureGate!(qc, [1], [sigmaZ()]) # Alice measures in a fixed basis=sigmaZ
qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=false)
iqs=createInitialQubitState(density, [[1.00,[1. 0.;1. 0.]],[0.00,[0. 1.;0. 1.]]])
qo=runQuantumProgram(qp, iqs, nrShots)
probeMeasureOutcome(qo, 5, "Probabilities of Alice's measurement when both Bob and Alice measure along the z-direction")

Alice = [qo.output[5,k].measured.outcome for k in 1:nrShots]
Bob   = [qo.output[4,k].measured.outcome for k in 1:nrShots]

traces=GenericTrace[]
push!(traces,scatter(x=1:nrShots, y=Bob, mode="lines", name="B"))
push!(traces,scatter(x=1:nrShots, y=Alice, mode="lines", name="A"))
layout=Layout(title_text="Outcomes of both Bob and Alice measuring along the z-direction", xaxis_title_text="Shot", yaxis_title_text="Outcome")
p6=plot(traces, layout)
display(p6)

# probability matrix P[Alice, Bob]
P = zeros(2,2)
for (a,b) in zip(Alice, Bob)
    P[a+1, b+1] += 1
end
P ./= nrShots
println("Joint probability matrix P(Alice, Bob):", P)
println("P(Alice=0) = ", sum(P[1,:]))
println("P(Alice=1) = ", sum(P[2,:]))

qc=createQuantumCircuit(2)
hGate!(qc, 1)
cnotGate!(qc, 1, 2) # entangled
measureGate!(qc, [1, 2], [sigmaZ(), sigmaZ()])
qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=false)
iqs=createInitialQubitState(density, [[1.00,[1. 0.;1. 0.]],[0.00,[0. 1.;0. 1.]]])
qo=runQuantumProgram(qp, iqs, nrShots)
probeMeasureOutcome(qo, 4, "Probabilities of Alice's and Bob's measurements when both Bob and Alice measure along the z-direction")

print("\n")
print("Bob and Alice measure along the x- and z-direction, respectively:\n")

qc=createQuantumCircuit(2)
hGate!(qc, 1)
cnotGate!(qc, 1, 2) # entangled
measureGate!(qc, [2], [sigmaX()]) # Bob tries to send a 1=sigmaX
measureGate!(qc, [1], [sigmaZ()]) # Alice measures in a fixed basis=sigmaZ
qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=false)
iqs=createInitialQubitState(density, [[1.00,[1. 0.;1. 0.]],[0.00,[0. 1.;0. 1.]]])
qo=runQuantumProgram(qp, iqs, nrShots)
probeMeasureOutcome(qo, 5, "Probabilities of Alice's measurement when Bob and Alice measure along the x- and z-direction, respectively")

Alice = [qo.output[5,k].measured.outcome for k in 1:nrShots]
Bob   = [qo.output[4,k].measured.outcome for k in 1:nrShots]

traces=GenericTrace[]
push!(traces,scatter(x=1:nrShots, y=Bob, mode="lines", name="B"))
push!(traces,scatter(x=1:nrShots, y=Alice, mode="lines", name="A"))
layout=Layout(title_text="Outcomes of Bob and Alice measuring along the x- and z-direction, respectively", xaxis_title_text="Shot", yaxis_title_text="Outcome")
p6=plot(traces, layout)
display(p6)

# probability matrix P[Alice, Bob]
P = zeros(2,2)
for (a,b) in zip(Alice, Bob)
    P[a+1, b+1] += 1
end
P ./= nrShots
println("Joint probability matrix P(Alice, Bob):", P)
println("P(Alice=0) = ", sum(P[1,:]))
println("P(Alice=1) = ", sum(P[2,:]))

qc=createQuantumCircuit(2)
hGate!(qc, 1)
cnotGate!(qc, 1, 2) # entangled
measureGate!(qc, [1, 2], [sigmaX(), sigmaZ()])
qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=false)
iqs=createInitialQubitState(density, [[1.00,[1. 0.;1. 0.]],[0.00,[0. 1.;0. 1.]]])
qo=runQuantumProgram(qp, iqs, nrShots)
probeMeasureOutcome(qo, 4, "Probabilities of Alice's and Bob's measurements when Bob and Alice measuring along the x- and z-direction, respectively")


print("\n")
print("2. Check quantum no-signaling theorem: reduced operator of 1st qubit not influenced by ignorant measurement of 2nd qubit.\n")

# Create entangled state of 2 qubits, check quantum no-signalling theorem by checking that
# the reduced operator of the 1st qubit before and after an ignorant measurement of the 2nd 
# qubit is the same.

qc=createQuantumCircuit(2)
hGate!(qc, 1)
cnotGate!(qc, 1, 2) # entangled
measureGate!(qc, [2], forgetOutcome=true)
qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=false)
iqs=createInitialQubitState(density, [[1.00,[1. 0.;1. 0.]],[0.00,[0. 1.;0. 1.]]])
qo=runQuantumProgram(qp, iqs, 1)
reducedrho2before=partialTrace(qo.output[3].rho, traceIndex=2)
reducedrho2after=partialTrace(qo.output[4].rho, traceIndex=2)

print("Reduced rho of 1st qubit before ignorant measurement of 2nd qubit=",reducedrho2before,"\n")
print("Reduced rho of 1st qubit after ignorant measurement of 2nd qubit=",reducedrho2after,"\n")
