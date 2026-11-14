############################################
# QuantumCircuitCompilation

function cascadeOperations(U1::UnitaryOperation, U2::UnitaryOperation)
	UnitaryOperation(U2.U * U1.U) # watch out for reverse order!!! first U1 and then U2, i.e., U2(U1(q))=U2*U1*q
end

function mergeTwoOperators(A::ComplexMatrix, B::ComplexMatrix, listA::Vector{Int64}, listB::Vector{Int64})
	nrBitsA = log2(size(A, 1))
	nrBitsB = log2(size(B, 1))
	nrBitsU = Int(nrBitsA + nrBitsB)
	N = 2^nrBitsU
	U = zeros(N, N)
	iA = ones(Int, N) # to force a Vector that can be used to index into Matrix
	if !isempty(listA)
		for k in 0:(N-1)
			for (m, element) in enumerate(listA)
				iA[k+1] += mod(floor(k / 2^(nrBitsU - element)), 2) * 2^(length(listA) - m);
			end
		end
	end
	iB = ones(Int, N) # to force a Vector that can be used to index into Matrix
	if !isempty(listB)
		for k in 0:(N-1)
			for (m, element) in enumerate(listB)
				iB[k+1] += mod(floor(k / 2^(nrBitsU - element)), 2) * 2^(length(listB) - m);
			end
		end
	end
	U = A[iA, iA] .* B[iB, iB]
end

function moveNQubitOperationToListOfQubits(U::UnitaryOperation, list::Vector{Int64}, n::Int)
	if size(U.U,2) == size(U.U,1) == 2^length(list)
		U = UnitaryOperation(mergeTwoOperators(U.U, createNQubitOperatorId(n-length(list)), list, setdiff(1:n, list)))
	else
		throw("size of UnitaryOperation does not match length of qubits list")
	end
end

function moveNQubitKrausOperatorToListOfQubits(E::KrausOperator, list::Vector{Int64}, n::Int)
	if size(E.E,2) == size(E.E,1) == 2^length(list)
		E = KrausOperator(mergeTwoOperators(E.E, createNQubitOperatorId(n-length(list)), list, setdiff(1:n, list)), E.label)
	else
		throw("size of KrausOperator does not match length of qubits list")
	end
end

@enum StepType unitary measure quantumChannel unknown

""" 
    compileQuantumCircuit(
		quantumCircuit::QuantumCircuit; 
		optimizeNumberOfSteps::Bool = false
	) -> QuantumProgram{Z, I}

Compiles a `quantumCircuit` into a `QuantumProgram` by translating each logical step into a corresponding quantum operation.

Each step in `quantumCircuit.circuit::Vector{Step}` is a vector of gates, where **all gates must be of the same gate type**. Supported gate types include:
- `UnitaryGate`
- `MeasureGate`
- `QuantumChannelGate`

### Step Compilation Behavior

- **Unitary Steps**:  
  For steps containing only `UnitaryGate`s, a tensor product of all gate unitary matrices in that step is constructed to form a single `UnitaryOperation`.

- **Measurement Steps**:  
  For steps containing `MeasureGate`s, a `MeasureOperation` or `MeasureAndForgetOperation` is created based on the qubits, measurement direction operators (`sigmas`) or Kraus operators (`krausOperators`), and measurement settings.

- **Quantum Channel Steps**:  
  For steps containing `QuantumChannelGate`s, the defined Kraus operators are used to form a `QuantumChannelOperation`.

### Optimization

If `optimizeNumberOfSteps = true`, all **consecutive steps** made up entirely of `UnitaryGate`s are merged into a **single** `UnitaryOperation` to reduce circuit depth and improve execution performance.

Otherwise (default behavior), each step is compiled into a **separate** `UnitaryOperation`.

# Type Parameters
- `Z`: Boolean type parameter indicating whether **zero-based numbering** is used.
- `I`: `IndexType` parameter specifying the **index mapping convention**.

# Arguments
- `quantumCircuit::QuantumCircuit` - The quantum circuit to be compiled.
- `optimizeNumberOfSteps::Bool` (optional) - Whether to combine consecutive unitary steps into one operation. Defaults to `false`.

# Returns
- `QuantumProgram{Z, I}`: The compiled quantum program ready for execution.

# See Also
- [`QuantumCircuit`](@ref)
- [`QuantumProgram`](@ref)
- [`Step`](@ref)
- [`UnitaryGate`](@ref)
- [`MeasureGate`](@ref)
- [`QuantumChannelGate`](@ref)
- [`UnitaryOperation`](@ref)
- [`MeasureOperation`](@ref)
- [`MeasureAndForgetOperation`](@ref)
- [`QuantumChannelOperation`](@ref)

# Examples

**Example 1: Compile a quantum circuit with separate Hadamard and Pauli-X gates**

```julia-repl
julia> qc = createQuantumCircuit(1)
julia> hGate!(qc, 1)
julia> xGate!(qc, 1)
julia> qc.circuit
2-element Vector{Step}:
 UnitaryStep(UnitaryGate[UnitaryGate{QubiSim.var"#11#12"}(QubiSim.var"#11#12"(), [1], "H")])
 UnitaryStep(UnitaryGate[UnitaryGate{QubiSim.var"#5#6"}(QubiSim.var"#5#6"(), [1], "X")])
julia> qp = compileQuantumCircuit(qc)
julia> qp.program
2-element Vector{Operation}:
 UnitaryOperation(ComplexF64[0.7071067811865476 + 0.0im 0.7071067811865475 - 8.659560562354932e-17im; 0.7071067811865475 + 0.0im -0.7071067811865476 + 8.659560562354934e-17im])
 UnitaryOperation(ComplexF64[6.123233995736766e-17 + 0.0im 1.0 - 1.2246467991473532e-16im; 1.0 + 0.0im -6.123233995736766e-17 + 7.498798913309288e-33im])
```

**Example 2: Compile the same circuit with step optimization**

```julia-repl
julia> qp = compileQuantumCircuit(qc; optimizeNumberOfSteps = true)
julia> qp.program
1-element Vector{Operation}:
 UnitaryOperation(ComplexF64[0.7071067811865475 + 0.0im -0.7071067811865476 + 0.0im; 0.7071067811865476 + 0.0im 0.7071067811865475 + 0.0im])
```
"""
function compileQuantumCircuit(quantumCircuit::QuantumCircuit; optimizeNumberOfSteps::Bool = false)
	quantumProgram = QuantumProgram(quantumCircuit.settings, quantumCircuit.numberOfQubits, Program())
    for step in quantumCircuit.circuit
        compileStep!(quantumProgram, step, quantumCircuit.numberOfQubits, optimizeNumberOfSteps)
    end
    return quantumProgram
end

function compileStep!(program::QuantumProgram, step::UnitaryStep, numberOfQubits::Int, optimize::Bool)
    UStep = createNQubitOperationId(numberOfQubits)
    for gate in step.gates
		UGate = gate.unitaryOperationFactory()
        UGate = moveNQubitOperationToListOfQubits(UGate, gate.qubits, numberOfQubits)
        UStep = cascadeOperations(UStep, UGate)
    end

	# addStepToProgram!(program, UStep, Val(optimize))
    if optimize && !isempty(program.program) && program.program[end] isa UnitaryOperation
        # Merge with previous unitary
        program.program[end] = cascadeOperations(program.program[end], UStep)
    else
        push!(program.program, UStep)
    end
end

# function addStepToProgram!(program::QuantumProgram, UStep::UnitaryOperation, ::Val{true})
# 	if !isempty(program.program) && program.program[end] isa UnitaryOperation
#         # Merge with previous unitary
#         program.program[end] = cascadeOperations(program.program[end], UStep)
# 	else
# 	    push!(program.program, UStep)
#     end
# end

# function addStepToProgram!(program::QuantumProgram, UStep::UnitaryOperation, ::Val{false})
#     push!(program.program, UStep)
# end

function compileStep!(program::QuantumProgram, step::MeasurementStep, numberOfQubits::Int, optimize::Bool)
    gate = step.gate
    qubits = gate.qubits
	if gate.measureGateType == PVM
		sigmas = isempty(gate.sigmas) ? [sigmaZ() for _ in qubits] : gate.sigmas
    	krausOperators = constructKrausOperatorsForPVMMeasurement(numberOfQubits, qubits, sigmas)
	elseif gate.measureGateType == POVM
		krausOperators = [moveNQubitKrausOperatorToListOfQubits(krausOperator, qubits, numberOfQubits) for krausOperator in gate.krausOperators]
	end

	measurementOperators = [krausOperators[m].E'*krausOperators[m].E for m in 1:length(krausOperators)]
	labels=[krausOperator.label for krausOperator in krausOperators]

	if gate.forgetOutcome
	    push!(program.program, MeasureAndForgetOperation(krausOperators, measurementOperators, labels))
	else
    	push!(program.program, MeasureOperation(krausOperators, measurementOperators, labels))
	end
end

function compileStep!(program::QuantumProgram, step::QuantumChannelStep, numberOfQubits::Int, optimize::Bool)
    gate = step.gate
    qubits = gate.qubits
	krausOperators = [moveNQubitKrausOperatorToListOfQubits(krausOperator, qubits, numberOfQubits) for krausOperator in gate.krausOperators]
    push!(program.program, QuantumChannelOperation(krausOperators))
end

function compileStep!(program::QuantumProgram, step::BernoulliStep, numberOfQubits::Int, optimize::Bool)
	UStep = step.gate.unitaryOperationFactory()
    UStep = moveNQubitOperationToListOfQubits(UStep, step.gate.qubits, numberOfQubits)

    push!(program.program, BernoulliOperation(UStep.U, step.gate.probability))
end

function compileStep!(program::QuantumProgram, step::StochasticStep, numberOfQubits::Int, optimize::Bool)
    push!(program.program, StochasticOperation(step.gate.unitaryOperationFactory, step.gate.qubits, numberOfQubits, step.gate.accuracyCheckForUnitarity))
end

function compileStep!(program::QuantumProgram, step::Step, numberOfQubits::Int, optimize::Bool)
    error("Unsupported quantum step type: $(typeof(step))")
end

# 

""" 
    compileToSingleGate(
		quantumCircuit::QuantumCircuit
	) -> UnitaryOperation

Compiles a given `quantumCircuit` that consist of only `UnitaryGate` types into a single `UnitaryOperation`.

# Arguments
- `quantumCircuit::QuantumCircuit` - The quantum circuit to compile.

# Returns
- `UnitaryOperation` - A single unitary representation of the quantum circuit.

# See Also
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`UnitaryOperation`](@ref)

# Throws
- A `String` exception if the quantum circuit does not contain only `UnitaryGate` types.

# Example

Step 1 - Compile a quantum circuit that applies a **Hadamard gate** followed by an **Pauli-X gate** into a single **HX unitary operation**:
```julia-repl
julia> qc = createQuantumCircuit(1)
julia> hGate!(qc, 1)
julia> xGate!(qc, 1)
julia> hxGate = compileToSingleGate(qc)
julia> hxGate.U
2×2 Matrix{ComplexF64}:
 0.707107+0.0im  -0.707107+0.0im
 0.707107+0.0im   0.707107+0.0im
```
Step 2 - Create and compile a quantum circuit containing the single **HX unitary operation* created above:
```julia-repl
julia> qc = createQuantumCircuit(1)
julia> unitaryUGate!(qc, [1], hxGate)
julia> qp = compileQuantumCircuit(qc)
julia> qp.program[1].U
2×2 Matrix{ComplexF64}:
 0.707107+0.0im  -0.707107+0.0im
 0.707107+0.0im   0.707107+0.0im
```
"""
function compileToSingleGate(quantumCircuit::QuantumCircuit)
	quantumProgram = compileQuantumCircuit(quantumCircuit; optimizeNumberOfSteps = true)
	if length(quantumProgram.program) == 1
		return quantumProgram.program[1]
	else
		throw("The quantumCircuit does not have only UnitaryGate types")
	end
end