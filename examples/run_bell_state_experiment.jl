using QubiSim
using LinearAlgebra, PlotlyJS

print("\n")
print("1. Create entangled output of 2 bits using density operator and measure 2 bits along arbitrary axes\n")

p=1
# 1=pure entangled, 0=impure, Bell's inequality is satisfied if p<1/sqrt(2), PPT condition is satisfied if p<1/3
iqs=createDoubleQubitWernerDensityState(p)

theta=Vector(range(start=0, stop=pi, step=0.05))
phi=zeros(Float64,length(theta))

corr1=zeros(Float64,length(theta))
corr2=zeros(Float64,length(theta))

for m in 1:length(theta)
    theta_a=0.
    phi_a=0.
    
    theta_b=theta[m]
    phi_b=phi[m]
    
    local qc=createQuantumCircuit(2)
    u3Gate!(qc, 1, theta_a, phi_a, pi-phi_a)
    u3Gate!(qc, 2, theta_b, phi_b, pi-phi_b)
    measureGate!(qc, [1,2])
    local qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

    local qo=runQuantumProgram(qp, iqs, 1)
    corr1[m]=0.25*calculateExpectationValueOfProductOfMeasuredQubits(qo, 3, useProbability)
    
    local nrShots=1000
    qo=runQuantumProgram(qp, iqs, nrShots)
    corr2[m]=0.25*calculateExpectationValueOfProductOfMeasuredQubits(qo, 3, useOutcome)
end

traces=GenericTrace[]
push!(traces,scatter(x=theta, y=corr1, mode="lines", name="expectation value"))
push!(traces,scatter(x=theta, y=corr2, mode="lines", name="random experiment"))
layout=Layout(title_text="Bell correlation evaluation of entangled bit-pairs", xaxis_title_text="theta[rad]", yaxis_title_text="correlation <σ(A),σ(B)>")
p1=plot(traces, layout)
display(p1)


print("\n")
print("2. Create entangled output of 2 bits using density operator and evaluate Bells inequality\n")

function evaluateCorrelationInTwoQubitSystem(iqs, theta_1, phi_1, theta_2, phi_2)
    qc=createQuantumCircuit(2)
    u3Gate!(qc, 1, theta_1, phi_1, pi-phi_1)
    u3Gate!(qc, 2, theta_2, phi_2, pi-phi_2)
    measureGate!(qc, [1,2])
    qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

    qo=runQuantumProgram(qp, iqs, 1)
    corr1=0.25*calculateExpectationValueOfProductOfMeasuredQubits(qo, 3, useProbability)
    
    nrShots=1000
    qo=runQuantumProgram(qp, iqs, nrShots)
    corr2=0.25*calculateExpectationValueOfProductOfMeasuredQubits(qo, 3, useOutcome)

    return corr1, corr2
end

p=1
# 1=pure entangled, 0=impure, Bell's inequality is satisfied if p<1/sqrt(2), PPT condition is satisfied if p<1/3
iqs=createDoubleQubitWernerDensityState(p)

theta=Vector(range(start=0, stop=pi/2, step=0.05))

corr1_ab=zeros(Float64,length(theta))
corr1_aba=zeros(Float64,length(theta))
corr1_bba=zeros(Float64,length(theta))
corr2_ab=zeros(Float64,length(theta))
corr2_aba=zeros(Float64,length(theta))
corr2_bba=zeros(Float64,length(theta))

for m in 1:length(theta)
    theta_a=0.
    phi_a=0.

    theta_b=pi/2
    phi_b=0.

    theta_ba=pi/2-theta[m]
    phi_ba=0.
    
    (corr1_ab[m],corr2_ab[m])=evaluateCorrelationInTwoQubitSystem(iqs,theta_a,phi_a,theta_b,phi_b)
    (corr1_aba[m],corr2_aba[m])=evaluateCorrelationInTwoQubitSystem(iqs,theta_a,phi_a,theta_ba,phi_ba)    
    (corr1_bba[m],corr2_bba[m])=evaluateCorrelationInTwoQubitSystem(iqs,theta_b,phi_b,theta_ba,phi_ba)
end

traces=GenericTrace[]
push!(traces,scatter(x=theta, y=abs.(corr1_ab.-corr1_aba), mode="lines", name="LHS: |<σ(a),σ(b)>-<σ(a),σ(ba)>|", line_color="blue"))
push!(traces,scatter(x=theta, y=0.25.+corr1_bba, mode="lines", name="RHS: 1/4+<σ(b),σ(ba)>", line_color="red"))
push!(traces,scatter(x=theta, y=abs.(corr2_ab.-corr2_aba), mode="lines", line_color="blue", showlegend=false, line = attr(dash="dot")))
push!(traces,scatter(x=theta, y=0.25.+corr2_bba, mode="lines", line_color="red", showlegend=false, line = attr(dash="dot")))
layout=Layout(title_text="Bells inequality |<σ(a),σ(b)>-<σ(a),σ(ba)>| <= 1/4+<σ(b),σ(ba)>", xaxis_title_text="theta[rad]")
p2=plot(traces, layout)
display(p2)
