using QubiSim
using LinearAlgebra, PlotlyJS

print("\n")
print("1. Create 2-qubit entangled state and measure the violation of the CHSH inequality |S|<2\n")

theta=Vector(range(start=0, stop=2*pi, step=0.1))
corr1=zeros(Float64,length(theta))
corr2=zeros(Float64,length(theta))

function evaluateCorrelationInTwoQubitSystem(iqs, theta, listOfPauliOperators, useFromMeasurement::UseFromMeasurement, nrShots)
    qc=createQuantumCircuit(2)
    hGate!(qc, 1)
    cnotGate!(qc, 1, 2)
    ryGate!(qc, 1, theta)
    measureGate!(qc, [1,2], listOfPauliOperators)
    qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
    qo=runQuantumProgram(qp, iqs, nrShots)
    expAB=calculateExpectationValueOfProductOfMeasuredQubits(qo, 3, useFromMeasurement)
    return expAB
end

nrShots=1
iqs=createInitialQubitState(vector, [1. 0.;1. 0.])
for m in 1:length(theta)
    expZZ=evaluateCorrelationInTwoQubitSystem(iqs, theta[m], [sigmaZ(),sigmaZ()], useProbability, nrShots)
    expZX=evaluateCorrelationInTwoQubitSystem(iqs, theta[m], [sigmaZ(),sigmaX()], useProbability, nrShots)
    expXZ=evaluateCorrelationInTwoQubitSystem(iqs, theta[m], [sigmaX(),sigmaZ()], useProbability, nrShots)
    expXX=evaluateCorrelationInTwoQubitSystem(iqs, theta[m], [sigmaX(),sigmaX()], useProbability, nrShots)

    corr1[m]=(expZZ-expZX+expXZ+expXX)
    corr2[m]=(expZZ+expZX-expXZ+expXX)
end

traces=GenericTrace[]
push!(traces,scatter(x=theta, y=corr1, mode="lines", name="S=<ZZ>-<ZX>+<XZ>+<XX>"))
push!(traces,scatter(x=theta, y=corr2, mode="lines", name="S=<ZZ>+<ZX>-<XZ>+<XX>"))
layout=Layout(title_text="CHSH inequality |S|<2", xaxis_title_text="theta[rad]", yaxis_title_text="Expected correlation |S|")
p1=plot(traces, layout)
display(p1)


nrShots=100
iqs=createInitialQubitState(vector, [1. 0.;1. 0.])
for m in 1:length(theta)
    expZZ=evaluateCorrelationInTwoQubitSystem(iqs, theta[m], [sigmaZ(),sigmaZ()], useOutcome, nrShots)
    expZX=evaluateCorrelationInTwoQubitSystem(iqs, theta[m], [sigmaZ(),sigmaX()], useOutcome, nrShots)
    expXZ=evaluateCorrelationInTwoQubitSystem(iqs, theta[m], [sigmaX(),sigmaZ()], useOutcome, nrShots)
    expXX=evaluateCorrelationInTwoQubitSystem(iqs, theta[m], [sigmaX(),sigmaX()], useOutcome, nrShots)

    corr1[m]=(expZZ-expZX+expXZ+expXX)
    corr2[m]=(expZZ+expZX-expXZ+expXX)
end

traces=GenericTrace[]
push!(traces,scatter(x=theta, y=corr1, mode="lines", name="S=<ZZ>-<ZX>+<XZ>+<XX>"))
push!(traces,scatter(x=theta, y=corr2, mode="lines", name="S=<ZZ>+<ZX>-<XZ>+<XX>"))
layout=Layout(title_text="CHSH inequality |S|<2", xaxis_title_text="theta[rad]", yaxis_title_text="Measured correlation |S|")
p1=plot(traces, layout)
display(p1)


print("\n")
print("2. Play CHSH game as described in IBM quantum information course.\n")

function calculateWinningProbabilityCHSHGame(x,y)
    qc=createQuantumCircuit(4)
    hGate!(qc,2)
    cnotGate!(qc,2,3)

    controlledUGate!(qc, [1], [2], createSingleQubitOperationRy(-pi/2.))
    xGate!(qc,1)
    controlledUGate!(qc, [1], [2], createSingleQubitOperationRy(0.))

    controlledUGate!(qc, [4], [3], createSingleQubitOperationRy(pi/4.))
    xGate!(qc,4)
    controlledUGate!(qc, [4], [3], createSingleQubitOperationRy(-pi/4.))

    measureGate!(qc,[2,3])

    qp=compileQuantumCircuit(qc)

    iqsWithQuestions=[1. 0.;1. 0.;1. 0.;1. 0.]
    if x==1
        iqsWithQuestions[1,:]=[0. 1.]
    end
    if y==1
        iqsWithQuestions[4,:]=[0. 1.]
    end
    iqs=createInitialQubitState(vector, iqsWithQuestions)

    qo=runQuantumProgram(qp, iqs, 1)
    probeMeasureProbability(qo, 7, "CHSH Game for questions (x,y)=("*string(x)*","*string(y)*")")
    return qo.output[end].measured.probability
end

probabilityAnswers00=calculateWinningProbabilityCHSHGame(0,0)
winningProbability00=probabilityAnswers00[1]+probabilityAnswers00[4] # win if a==b so for 00 and 11

probabilityAnswers01=calculateWinningProbabilityCHSHGame(0,1)
winningProbability01=probabilityAnswers01[1]+probabilityAnswers01[4] # win if a==b so for 00 and 11

probabilityAnswers10=calculateWinningProbabilityCHSHGame(1,0)
winningProbability10=probabilityAnswers10[1]+probabilityAnswers10[4] # win if a==b so for 00 and 11

probabilityAnswers11=calculateWinningProbabilityCHSHGame(1,1)
winningProbability11=probabilityAnswers11[2]+probabilityAnswers11[3] # win if a!=b so for 01 and 10

print("\n")
print("Winning strategy for questions (x,y) is answers (a,b) with a⊕b=x∧y\n")
print("Winning probability for questions (x,y)=(0,0) is answers (a,b) with a==b so for (0,0) and (1,1): "*string(winningProbability00)*"\n")
print("Winning probability for questions (x,y)=(0,1) is answers (a,b) with a==b so for (0,0) and (1,1): "*string(winningProbability01)*"\n")
print("Winning probability for questions (x,y)=(1,0) is answers (a,b) with a==b so for (0,0) and (1,1): "*string(winningProbability10)*"\n")
print("Winning probability for questions (x,y)=(1,1) is answers (a,b) with a!=b so for (0,1) and (1,0): "*string(winningProbability11)*"\n")
