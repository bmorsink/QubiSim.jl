# Complete function and type documentation

## Quantum Circuits

### Functions and types for constructing circuits and adding gates

```@docs
IndexType
```

```@docs
Settings
```

```@docs
createSettings
```

```@docs
QuantumCircuit
```

```@docs
Circuit
```

```@docs
Gate
```

```@docs
Step
```

```@docs
Qubits
```

```@docs
KrausOperator
```

```@docs
KrausOperators
```

```@docs
Sigmas
```

```@docs
UnitaryGate
```

```@docs
qubits
```

```@docs
MeasureGateType
```

```@docs
MeasureGate
```

```@docs
QuantumChannelGate
```

```@docs
UnitaryStep
```

```@docs
MeasurementStep
```

```@docs
QuantumChannelStep
```

```@docs
getStep
```

```@docs
getGate
```

```@docs
createQuantumCircuit
```

```@docs
createQuantumCircuitWithSettings
```

```@docs
concatenateQuantumCircuits
```

```@docs
inverseQuantumCircuit
```

```@docs
barrier!
```

```@docs
unitaryUGate!
```

```@docs
measureGate!
```

```@docs
quantumChannelGate!
```

```@docs
generateKrausOperatorsForPVMMeasurement
```

```@docs
generateKrausOperatorsForDepolarizingChannel
```

```@docs
generateKrausOperatorsForPhaseDampingChannel
```

```@docs
generateKrausOperatorsForAmplitudeDampingChannel
```

```@docs
pauliX
```

```@docs
pauliY
```

```@docs
pauliZ
```

```@docs
sigmaX
```

```@docs
sigmaY
```

```@docs
sigmaZ
```

```@docs
sigmaN
```

```@docs
u1Gate!
```

```@docs
u2Gate!
```

```@docs
u3Gate!
```

```@docs
xGate!
```

```@docs
yGate!
```

```@docs
zGate!
```

```@docs
hGate!
```

```@docs
idGate!
```

```@docs
rxGate!
```

```@docs
ryGate!
```

```@docs
rzGate!
```

```@docs
rotationGate!
```

```@docs
projectionGate!
```

```@docs
reflectionGate!
```

```@docs
tGate!
```

```@docs
tdGate!
```

```@docs
sGate!
```

```@docs
sdGate!
```

```@docs
cnotGate!
```

```@docs
cnotReverseGate!
```

```@docs
swapGate!
```

```@docs
phaseGate!
```

```@docs
expHGate!
```

```@docs
qftGate!
```

```@docs
createNQubitQFTQuantumCircuit
```

```@docs
iqftGate!
```

```@docs
createNQubitIQFTQuantumCircuit
```

```@docs
controlledUGate!
```

```@docs
controlledHGate!
```

```@docs
controlledXGate!
```

```@docs
controlledYGate!
```

```@docs
controlledZGate!
```

```@docs
controlledTGate!
```

```@docs
controlledTdGate!
```

```@docs
controlledSGate!
```

```@docs
controlledSdGate!
```

```@docs
qpeGate!
```

```@docs
createNQubitQPEQuantumCircuit
```

## Compilation

### Tools for translating circuits into executable quantum programs

```@docs
Program
```

```@docs
Operation
```

```@docs
QuantumProgram
```

```@docs
UnitaryOperation
```

```@docs
MeasureOperation
```

```@docs
MeasureAndForgetOperation
```

```@docs
QuantumChannelOperation
```

```@docs
getOperation
```

```@docs
createSingleQubitOperationU1
```

```@docs
createSingleQubitOperationU2
```

```@docs
createSingleQubitOperationU3
```

```@docs
createSingleQubitOperationX
```

```@docs
createSingleQubitOperationY
```

```@docs
createSingleQubitOperationZ
```

```@docs
createSingleQubitOperationH
```

```@docs
createSingleQubitOperationId
```

```@docs
createSingleQubitOperationRx
```

```@docs
createSingleQubitOperationRy
```

```@docs
createSingleQubitOperationRz
```

```@docs
createSingleQubitOperationRotation
```

```@docs
createNQubitOperationProjection
```

```@docs
createNQubitOperationReflection
```

```@docs
createSingleQubitOperationT
```

```@docs
createSingleQubitOperationTd
```

```@docs
createSingleQubitOperationS
```

```@docs
createSingleQubitOperationSd
```

```@docs
createDoubleQubitOperationCNOT
```

```@docs
createDoubleQubitOperationCNOTReverse
```

```@docs
createDoubleQubitOperationSWAP
```

```@docs
createDoubleQubitOperationPhase
```

```@docs
createNQubitOperationExpH
```

```@docs
createNQubitOperationQFT
```

```@docs
createNQubitOperationIQFT
```

```@docs
createNQubitOperationControlledU
```

```@docs
createDoubleQubitOperationControlledH
```

```@docs
createDoubleQubitOperationControlledX
```

```@docs
createDoubleQubitOperationControlledY
```

```@docs
createDoubleQubitOperationControlledZ
```

```@docs
createDoubleQubitOperationControlledT
```

```@docs
createDoubleQubitOperationControlledTd
```

```@docs
createDoubleQubitOperationControlledS
```

```@docs
createDoubleQubitOperationControlledSd
```

```@docs
createNQubitOperationQPE
```

```@docs
compileQuantumCircuit
```

```@docs
compileToSingleGate
```


## Execution

### Functionality for running quantum programs on quantum states

```@docs
QuantumOutput
```

```@docs
Measured
```

```@docs
State
```

```@docs
VectorState
```

```@docs
DensityState
```

```@docs
convertVectorStateToDensityState
```

```@docs
getState
```

```@docs
getShot
```

```@docs
getStateShot
```

```@docs
createByteIndexVector
```

```@docs
StateSpace 
```

```@docs
createInitialQubitState
```

```@docs
createDoubleQubitWernerDensityState
```

```@docs
createSingleQubitBlochDensityState
```

```@docs
runQuantumProgram
```


## Results

### Tools for visualizing qubit states and measurement outcomes

```@docs
probeMeasureOutcome
```

```@docs
probeMeasureProbability
```

```@docs
probeStateProbability
```

```@docs
probeStateMultiBlochVector
```

```@docs
UseFromMeasurement
```

```@docs
calculateExpectationValueOfProductOfMeasuredQubits
```

## Miscellaneous

### Supporting and utility functions and types

```@docs
ComplexMatrix
```

```@docs
ComplexVector
```

```@docs
tensorProduct
```

```@docs
partialTrace
```

```@docs
fidelity
```

```@docs
calculateEntropyAndPurity
```

```@docs
calculateCorrelation
```

```@docs
createToggleSwapList
```
