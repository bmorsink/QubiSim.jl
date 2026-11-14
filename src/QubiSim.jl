module QubiSim

export ComplexVector, ComplexMatrix
export IndexType, circuitIndexBigEndian, circuitIndexLittleEndian, byteIndex
export StateSpace, vector, density
export UseFromMeasurement, useProbability, useOutcome

export Settings, createSettings, zeroBasedNumbering, indexType

export Qubits, Gate, Step, Circuit, Sigmas, Operation, Program, State, KrausOperator, KrausOperators
export QuantumCircuit, UnitaryStep, UnitaryGate, qubits, MeasurementStep, MeasureGateType, MeasureGate, QuantumChannelStep, QuantumChannelGate, getStep, getGate
export QuantumProgram, UnitaryOperation, MeasureOperation, MeasureAndForgetOperation, QuantumChannelOperation, getOperation
export QuantumOutput, Measured, VectorState, DensityState, getState, getShot, getStateShot, convertVectorStateToDensityState

export createQuantumCircuit, createQuantumCircuitWithSettings, concatenateQuantumCircuits, inverseQuantumCircuit
export compileQuantumCircuit, compileToSingleGate

export createByteIndexVector, createInitialQubitState, createDoubleQubitWernerDensityState, createSingleQubitBlochDensityState
export runQuantumProgram, calculateEntropyAndPurity
export calculateCorrelation, calculateExpectationValueOfProductOfMeasuredQubits

export createToggleSwapList, tensorProduct, partialTrace, fidelity

export barrier!
export unitaryUGate!
export measureGate!, sigmaX, sigmaY, sigmaZ, sigmaN, pauliX, pauliY, pauliZ, generateKrausOperatorsForPVMMeasurement
export generateKrausOperatorsForDepolarizingChannel, generateKrausOperatorsForPhaseDampingChannel, generateKrausOperatorsForAmplitudeDampingChannel, generateKrausOperatorsForBitFlipChannel
export quantumChannelGate!
export u1Gate!, createSingleQubitOperationU1
export u2Gate!, createSingleQubitOperationU2
export u3Gate!, createSingleQubitOperationU3
export xGate!, createSingleQubitOperationX
export yGate!, createSingleQubitOperationY
export zGate!, createSingleQubitOperationZ
export hGate!, createSingleQubitOperationH
export idGate!, createSingleQubitOperationId
export rxGate!, createSingleQubitOperationRx
export ryGate!, createSingleQubitOperationRy
export rzGate!, createSingleQubitOperationRz
export rotationGate!, createSingleQubitOperationRotation
export projectionGate!, createNQubitOperationProjection
export reflectionGate!, createNQubitOperationReflection
export tGate!, createSingleQubitOperationT
export tdGate!, createSingleQubitOperationTd
export sGate!, createSingleQubitOperationS
export sdGate!, createSingleQubitOperationSd
export cnotGate!, createDoubleQubitOperationCNOT
export cnotReverseGate!, createDoubleQubitOperationCNOTReverse
export swapGate!, createDoubleQubitOperationSWAP
export phaseGate!, createDoubleQubitOperationPhase
export expHGate!, createNQubitOperationExpH
export qftGate!, createNQubitOperationQFT, createNQubitQFTQuantumCircuit
export iqftGate!, createNQubitOperationIQFT, createNQubitIQFTQuantumCircuit
export controlledUGate!, createNQubitOperationControlledU
export controlledHGate!, createDoubleQubitOperationControlledH
export controlledXGate!, createDoubleQubitOperationControlledX
export controlledYGate!, createDoubleQubitOperationControlledY
export controlledZGate!, createDoubleQubitOperationControlledZ
export controlledTGate!, createDoubleQubitOperationControlledT
export controlledTdGate!, createDoubleQubitOperationControlledTd
export controlledSGate!, createDoubleQubitOperationControlledS
export controlledSdGate!, createDoubleQubitOperationControlledSd
export qpeGate!, createNQubitOperationQPE, createNQubitQPEQuantumCircuit

export BernoulliGate, BernoulliStep, bernoulliGate!, BernoulliOperation
export StochasticGate, StochasticStep, stochasticGate!, StochasticOperation

export probeMeasureOutcome, probeMeasureProbability, probeStateProbability, probeStateMultiBlochVector

using LinearAlgebra, PlotlyJS, StatsBase

import Base: broadcastable

include("MiscellaneousTypes.jl")

include("QuantumConfigurationDefinitions.jl")

include("QuantumCircuitDefinitions.jl")

include("QuantumProgramDefinitions.jl")

include("QuantumOutputDefinitions.jl")

include("SingleQubitUnitaryGates.jl")

include("NQubitUnitaryGates.jl")

include("DoubleQubitUnitaryGates.jl")

include("OtherNQubitGates.jl")

include("QuantumCircuitCreation.jl")

include("MeasurementDirectionOperators.jl")

include("QuantumCircuitCompilation.jl")

include("QuantumStateCreation.jl")

include("QuantumProgramExecution.jl")

include("KrausOperatorsCreation.jl")

include("MiscellaneousFunctions.jl")

include("PlottingFunctions.jl")

end