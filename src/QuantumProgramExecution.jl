############################################
# QuantumProgramExecution

mutable struct Cache
	cV1::ComplexMatrix
	cV2::ComplexMatrix
	cM1::ComplexMatrix
	cM2::ComplexMatrix
	outcome::Int
end

function applyOperationOnQubitState(unitaryOperation::UnitaryOperation, vectorState::VectorState, cache::Cache)
	return VectorState(unitaryOperation.U * vectorState.q)
end

function applyOperationOnQubitState(unitaryOperation::UnitaryOperation, densityState::DensityState, cache::Cache)
    mul!(cache.cM1, densityState.rho, (unitaryOperation.U'))
	return DensityState(unitaryOperation.U * cache.cM1)
end

function applyOperationOnQubitState(measureOperation::MeasureOperation, vectorState::VectorState, cache::Cache)
	prob = zeros(length(measureOperation.krausOperators))
	mul!(cache.cM1, vectorState.q, vectorState.q')
	for m in axes(measureOperation.krausOperators, 1)
		mul!(cache.cM2, measureOperation.measurementOperators[m], cache.cM1)
		prob[m] = real(tr(cache.cM2))
	end
	cache.outcome = sampleMeasurementOutcome(prob)
	mul!(cache.cV1, vectorState.q, 1 / sqrt(prob[cache.outcome+1]))
	return VectorState(measureOperation.measurementOperators[cache.outcome+1] * cache.cV1; measured = Measured(prob, cache.outcome, measureOperation.labels))
end

function applyOperationOnQubitState(measureOperation::MeasureOperation, densityState::DensityState, cache::Cache)
	prob = zeros(length(measureOperation.krausOperators))
	for m in axes(measureOperation.krausOperators, 1)
		mul!(cache.cM2, measureOperation.measurementOperators[m], densityState.rho)
		prob[m] = real(tr(cache.cM2))
	end
	cache.outcome = sampleMeasurementOutcome(prob)
	mul!(cache.cM1, measureOperation.krausOperators[cache.outcome+1].E', 1 / prob[cache.outcome+1])
	mul!(cache.cM2, densityState.rho, cache.cM1)
	return DensityState(measureOperation.krausOperators[cache.outcome+1].E * cache.cM2; measured = Measured(prob, cache.outcome, measureOperation.labels))
end

function applyOperationOnQubitState(measureAndForgetOperation::MeasureAndForgetOperation, vectorState::VectorState, cache::Cache)
	throw("Combination of a MeasureAndForgetOperation with a VectorState is not possible...")
end

function applyOperationOnQubitState(measureAndForgetOperation::MeasureAndForgetOperation, densityState::DensityState, cache::Cache)
	prob = zeros(length(measureAndForgetOperation.krausOperators))
	rhoa = zeros(ComplexF64, size(densityState.rho))
	for m in axes(measureAndForgetOperation.krausOperators, 1)
		mul!(cache.cM1, measureAndForgetOperation.measurementOperators[m], densityState.rho)
		prob[m] = real(tr(cache.cM1))
		mul!(cache.cM1, densityState.rho, measureAndForgetOperation.krausOperators[m].E')
		mul!(cache.cM2, measureAndForgetOperation.krausOperators[m].E, cache.cM1)
		rhoa .+= cache.cM2
	end
	return DensityState(rhoa; measured = Measured(prob, nothing, measureAndForgetOperation.labels))
end

function sampleMeasurementOutcome(probabilities::Vector{Float64})
	p = rand()
	cumulative = cumsum(probabilities, dims=1)
	indices = findmin(findall(x -> p <= x, cumulative))
	return first(indices)-1
end

function applyOperationOnQubitState(quantumChannelOperation::QuantumChannelOperation, vectorState::VectorState, cache::Cache)
	throw("Combination of a QuantumChannelOperation with a VectorState is not possible...")
end

function applyOperationOnQubitState(quantumChannelOperation::QuantumChannelOperation, densityState::DensityState, cache::Cache)
	rhoa = zeros(ComplexF64, size(densityState.rho))
	for krausOperator in quantumChannelOperation.krausOperators
		mul!(cache.cM1, densityState.rho, krausOperator.E')
		mul!(cache.cM2, krausOperator.E, cache.cM1)
		rhoa .+= cache.cM2
	end
	return DensityState(rhoa)
end

function applyOperationOnQubitState(bernoulliOperation::BernoulliOperation, vectorState::VectorState, cache::Cache)
	if rand() <= bernoulliOperation.probability
		return VectorState(bernoulliOperation.U * vectorState.q)
	else
		return vectorState
	end
end

function applyOperationOnQubitState(bernoulliOperation::BernoulliOperation, densityState::DensityState, cache::Cache)
	if rand() <= bernoulliOperation.probability
	    mul!(cache.cM1, densityState.rho, (bernoulliOperation.U'))
		return DensityState(bernoulliOperation.U * cache.cM1)
	else
		return densityState
	end
end

function applyOperationOnQubitState(stochasticOperation::StochasticOperation, vectorState::VectorState, cache::Cache)
	UStep = stochasticOperation.unitaryOperationFactory()
	if all(abs.(UStep.U'*UStep.U-ComplexMatrix(I,size(UStep.U))).<stochasticOperation.accuracyCheckForUnitarity)
    	UStep = moveNQubitOperationToListOfQubits(UStep, stochasticOperation.qubits, stochasticOperation.numberOfQubits)
		return VectorState(UStep.U * vectorState.q)
	else
		throw("U:$UStep is not unitary")
	end
end

function applyOperationOnQubitState(stochasticOperation::StochasticOperation, densityState::DensityState, cache::Cache)
	UStep = stochasticOperation.unitaryOperationFactory()
	if all(abs.(UStep.U'*UStep.U-ComplexMatrix(I,size(UStep.U))).<stochasticOperation.accuracyCheckForUnitarity)
    	UStep = moveNQubitOperationToListOfQubits(UStep, stochasticOperation.qubits, stochasticOperation.numberOfQubits)
		mul!(cache.cM1, densityState.rho, (UStep.U'))
		return DensityState(UStep.U * cache.cM1)
	else
		throw("U:$UStep is not unitary")
	end
end

""" 
    runQuantumProgram(
		quantumProgram::QuantumProgram, 
		initialQubitState::S, 
		numberOfShots::Int
	) -> QuantumOutput{Z, I, S} where {Z, I, S <: State}

Executes a quantum program multiple times using a specified initial quantum state.

This function runs the provided `quantumProgram` starting from the `initialQubitState` and simulates its execution `numberOfShots` times. The result includes the full quantum state evolution across all operations for each shot.

This function is **compile-time specialized** on type parameter `S <: State` (that forces the struct to hold type-stable quantum state objects), ensuring optimized and branch-free code generation.

# Type parameters
- `Z`: Boolean type parameter indicating whether **zero-based numbering** is used.
- `I`: `IndexType` parameter specifying the **index mapping convention**.
- `S <: State`: The concrete type of the quantum states being stored, e.g., `VectorState` or `DensityState`.

# Arguments
- `quantumProgram::QuantumProgram` - A compiled quantum program containing a sequence of quantum operations.
- `initialQubitState::S` - The starting state for all shots. This can be a `VectorState` (pure state) or `DensityState` (mixed state).
- `numberOfShots::Int` - Number of repetitions (shots) to run the program. Each shot evolves independently.

# Supported State Types
- `VectorState` - For pure states.
- `DensityState` - For mixed states.

# Supported Operations
- `UnitaryOperation` - Applies a reversible gate.
- `MeasureOperation` - Simulates measurement with revealing outcome by collapsing the state to one of the possible outcomes.
- `MeasureAndForgetOperation` - Simulates measurement without revealing outcome by creating a mixed state over all possible outcomes.
- `QuantumChannelOperation` - Applies a general (possibly non-unitary) quantum channel.

# Returns
- `QuantumOutput{Z, I, S}` - Contains all intermediate states across all operations and shots. The `.output` field is a `Matrix{S}` of size `(numberOfStates, numberOfShots)`.

# Output Layout
If the program has `n` operations, the output is an `(n+1) × numberOfShots` matrix, where:
- Row `1` contains the initial state.
- Row `2` through `n+1` contain the resulting states after each operation.
- Each column corresponds to a different shot.

# Notes
- Each shot is processed independently; the program is re-applied from scratch for each shot.
- The function applies the sequence of operations to the initial state and stores the intermediate results for analysis or visualization.
- Measurement outcomes and post-measurement states are captured in the final states.

# Quantum state evolution 

- Under a unitary operation
    - Let U be the unitary matrix defining the unitary operation.
        - **Quantum vector state**:  
            - A pure state ``|ψ⟩`` evolves as: ``|ψ⟩ → U·|ψ⟩``
        - **Quantum density state**:
            - A density operator ``ρ`` evolves as: ``ρ → U·ρ·U^†``

- Under a generalized measurement (POVM)
    - Let ``{Mₖ}`` be measurement operators satisfying ``∑ₖ Mₖ^†·Mₖ = I``. Each operator is given by ``Mₖ = Eₖ^†·Eₖ`` where ``{Eₖ}`` are the Kraus operators.
        - **Quantum vector state**:
            - Probability of outcome ``k``: ``pₖ = Tr(Mₖ·|ψ⟩⟨ψ|)``
            - Post-measurement state for outcome ``k``: ``|ψ⟩ → Mₖ·|ψ⟩/√pₖ``
        - **Quantum density state**:
            - Probability of outcome ``k``: ``pₖ = Tr(Mₖ·ρ)``
            - Post-measurement state (outcome ``k`` observed): ``ρ → Eₖ·ρ·Eₖ^†/pₖ``
            - Post-measurement state (outcome not observed): ``ρ → ∑ₖ Eₖ·ρ·Eₖ^†``

- Under a projective measurement (PVM)
    - A **PVM** is a special case of the **POVM** where each Kraus operator has the form: ``Eₖ = T^†·|k⟩⟨k|·T``. Here, ``|k⟩⟨k|`` denotes the rank-1 projection operator for outcome ``k``, while ``T`` is the tensor product, taken over the measured qubits, of the corresponding **measurement direction operators** ``M(nₖ)`` as defined in `Sigmas`. Then the corresponding measurement operator becomes: ``Mₖ = Eₖ^†·Eₖ = T^†·|k⟩⟨k|·T = Eₖ``.

- Under a quantum channel
    - Let ``{Eₖ}`` be the Kraus operators of the quantum channel.
        - **Quantum density state**:
            - A density operator ``ρ`` evolves as: ``ρ → ∑ₖ Eₖ·ρ·Eₖ^†``

# See Also
- [`QuantumProgram`](@ref)
- [`QuantumOutput`](@ref)
- [`State`](@ref)
- [`VectorState`](@ref)
- [`DensityState`](@ref)
- [`UnitaryOperation`](@ref)
- [`MeasureOperation`](@ref)
- [`MeasureAndForgetOperation`](@ref)
- [`QuantumChannelOperation`](@ref)
- [`Sigmas`](@ref)
- [`KrausOperators`](@ref)

# Examples

## Example 1: Single-qubit vector initial state ``|0⟩``, Hadamard gate, then measurement. Take 5 shots and extract the last state after the measurement.
```julia-repl
julia> qc=createQuantumCircuit(1)
julia> hGate!(qc,1)
julia> measureGate!(qc,[1])
julia> qp=compileQuantumCircuit(qc)
julia> iqs=createInitialQubitState(vector, [1. 0.;])
julia> qo=runQuantumProgram(qp, iqs, 5)
julia> getState(qo, 3)
5-element Vector{VectorState}:
 VectorState(ComplexF64[1.0 + 0.0im; 0.0 + 0.0im;;], Measured([0.5000000000000001, 0.4999999999999999], 0, ["0", "1"]))
 VectorState(ComplexF64[0.0 + 0.0im; 1.0 + 0.0im;;], Measured([0.5000000000000001, 0.4999999999999999], 1, ["0", "1"]))
 VectorState(ComplexF64[0.0 + 0.0im; 1.0 + 0.0im;;], Measured([0.5000000000000001, 0.4999999999999999], 1, ["0", "1"]))
 VectorState(ComplexF64[0.0 + 0.0im; 1.0 + 0.0im;;], Measured([0.5000000000000001, 0.4999999999999999], 1, ["0", "1"]))
 VectorState(ComplexF64[0.0 + 0.0im; 1.0 + 0.0im;;], Measured([0.5000000000000001, 0.4999999999999999], 1, ["0", "1"]))
```

## Example 2: Same circuit with a zero mixed (density) input state. Take 5 shots and extract the last state after the measurement.
```julia-repl
julia> qc=createQuantumCircuit(1)
julia> hGate!(qc,1)
julia> measureGate!(qc,[1])
julia> qp=compileQuantumCircuit(qc)
julia> iqs=createInitialQubitState(density, [[1.00,[1. 0.;]]])
julia> qo=runQuantumProgram(qp, iqs, 5)
julia> getState(qo, 3)
5-element Vector{DensityState}:
 DensityState(ComplexF64[1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im], 0.0, 1.0, 1.0, Measured([0.5000000000000001, 0.4999999999999999], 0, ["0", "1"]))
 DensityState(ComplexF64[1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im], 0.0, 1.0, 1.0, Measured([0.5000000000000001, 0.4999999999999999], 0, ["0", "1"]))
 DensityState(ComplexF64[1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im], 0.0, 1.0, 1.0, Measured([0.5000000000000001, 0.4999999999999999], 0, ["0", "1"]))
 DensityState(ComplexF64[1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im], 0.0, 1.0, 1.0, Measured([0.5000000000000001, 0.4999999999999999], 0, ["0", "1"]))
 DensityState(ComplexF64[1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im], 0.0, 1.0, 1.0, Measured([0.5000000000000001, 0.4999999999999999], 0, ["0", "1"]))
```
"""
function runQuantumProgram(quantumProgram::QuantumProgram, initialQubitState::S, numberOfShots::Int) where {S <: State}
	N = 2^quantumProgram.numberOfQubits
	cache = Cache(ComplexMatrix(undef, N, 1), ComplexMatrix(undef, N, 1), ComplexMatrix(undef, N, N), ComplexMatrix(undef, N, N), 0)
	numberOfStates = length(quantumProgram.program)+1

    # Create a type-stable matrix with concrete type S. We use `undef` because we will fill it in the loop.
    output = Matrix{S}(undef, numberOfStates, numberOfShots)

    # Initialize the first state for all shots (using broadcasting (.=) for VectorState and DensityState).
    output[1, :] .= initialQubitState

    # Loop in cache-friendly order (over columns inside loopOverShots!)
	for (m, operation) in enumerate(quantumProgram.program)
		# call _loopOverShots with type-stable operation::O where {O <: Operation}
		loopOverShots!(operation, output, m, numberOfShots, cache)
    end

	return QuantumOutput(quantumProgram.settings, quantumProgram.numberOfQubits, numberOfStates, numberOfShots, output)
end

function loopOverShots!(operation::O, output::Matrix{S}, m::Int, numberOfShots::Int, cache::Cache) where {O <: Operation, S <: State}
    @inbounds @simd for k in 1:numberOfShots
        output[m+1, k] = applyOperationOnQubitState(operation, output[m, k], cache)
    end
end
