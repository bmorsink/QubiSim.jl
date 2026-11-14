############################################
# MiscellaneousFunctions

""" 
    tensorProduct(
		a, 
		b
	) -> AbstractArray

Computes the Kronecker (tensor) product ``a ⊗ b`` of two quantum operations or states ``a`` and ``b``.

This operation is used to combine two quantum systems into a single composite system. It supports both vectors (quantum states) and matrices (quantum operations).

# Arguments
- `a` - A vector or matrix ``a`` representing a quantum state or operation.
- `b` - A vector or matrix ``b`` representing a quantum state or operation.

# Returns
- `AbstractArray` - A new array representing the tensor (Kronecker) product of ``a`` and ``b``. 

# Examples

**Example 1: Tensor product of the ``|0⟩`` and ``|+⟩`` vector states**
```julia-repl
julia> tensorProduct([1; 0], [1; 1]/sqrt(2))
4-element Vector{Float64}:
 0.7071067811865475
 0.7071067811865475
 0.0
 0.0
```
**Example 2: Tensor product of the identity and Hadamard unitary operations**
```julia-repl
julia> tensorProduct([1 0; 0 1], [1 1; 1 -1]/sqrt(2))
4×4 Matrix{Float64}:
 0.707107   0.707107  0.0        0.0
 0.707107  -0.707107  0.0       -0.0
 0.0        0.0       0.707107   0.707107
 0.0       -0.0       0.707107  -0.707107
```
**Example 3: Tensor product of the Hadamard and identity unitary operations**
```julia-repl
julia> tensorProduct([1 1; 1 -1]/sqrt(2), [1 0; 0 1])
4×4 Matrix{Float64}:
 0.707107  0.0        0.707107   0.0
 0.0       0.707107   0.0        0.707107
 0.707107  0.0       -0.707107  -0.0
 0.0       0.707107  -0.0       -0.707107
```
"""
function tensorProduct(a, b)
	kron(a, b)
end

"""
    createToggleSwapList(
		numberOfQubits::Int
	) -> Vector{Int}

Generates a list of indices representing bit-reversed ordering for a quantum system with `numberOfQubits` qubits.

# Arguments
- `numberOfQubits::Int` - The number of qubits ``M``.

# Returns
- `Vector{Int}` - A list of indices (1-based) of length ``2^M`` that perform bit-reversal reordering.

# Notes
- The returned list can be used to permute state vectors or matrices: `UnitaryOperation.U[swapList, swapList]`, `VectorState.q[swapList]` and `DensityState.rho[swapList, swapList]`.

# See Also
- [`UnitaryOperation`](@ref)
- [`VectorState`](@ref)
- [`DensityState`](@ref)

# Examples

**Example 1: Reordering for 2 qubits**

Binary representations:
- index `1` → `0` (00) → (00) `0` → index `1`
- index `2` → `1` (01) → (10) `2` → index `3`
- index `3` → `2` (10) → (01) `1` → index `2`
- index `4` → `3` (11) → (11) `3` → index `4`

```julia-repl
julia> createToggleSwapList(2)
4-element Vector{Int64}:
 1
 3
 2
 4
```

**Example 2: Reordering for 3 qubits**

Binary representations:
- index `1` → `0` (000) → (000) `0` → index `1`
- index `2` → `1` (001) → (100) `4` → index `5`
- index `3` → `2` (010) → (010) `2` → index `3`
- index `4` → `3` (011) → (110) `6` → index `7`
- index `5` → `4` (100) → (001) `1` → index `2`
- index `6` → `5` (101) → (101) `5` → index `6`
- index `7` → `6` (110) → (011) `3` → index `4`
- index `8` → `7` (111) → (111) `7` → index `8`

```julia-repl
julia> createToggleSwapList(3)
8-element Vector{Int64}:
 1
 5
 3
 7
 2
 6
 4
 8
```
"""
function createToggleSwapList(numberOfQubits::Int)
	n = 2^numberOfQubits
	swapList = ones(Int, n)
	for m in 1:n
        # Compute bit-reversed index (0-based), convert to 1-based
  		# bits = reverse(digits(m - 1, base = 2, pad = numberOfQubits)')
  		# swapList[m] = sum([bits[i] * 2^(i - 1) for i in 1:numberOfQubits]) + 1
  		bits = reverse(digits(m - 1, base = 2, pad = numberOfQubits))
  		swapList[m] = sum(bits[i] * 2^(i - 1) for i in 1:numberOfQubits) + 1
	end
	return swapList
end

""" 
    calculateEntropyAndPurity(
		rho; 
		accuracyCheckForPurity = 10 * eps(1.0)
	) -> (entropy, purity, purityNorm, purityString)

Computes the **von Neumann entropy** and **purity** of a quantum density matrix ``ρ``.

This function analyzes the degree of *mixedness* of a quantum state. The entropy quantifies the uncertainty (disorder), and the purity characterizes how close the state is to being a pure state. Results include a normalized purity and a human-readable description of the state’s purity level.

# Arguments
- `rho::Matrix{<:Complex}` - Density matrix ``ρ`` representing a quantum state (Hermitian, positive semi-definite, trace 1).
- `accuracyCheckForPurity::Float64` (optional) - Numerical tolerance to determine whether the purity is effectively 1 or `1/n`. Defaults to `10 * eps(1.0)`.

# Returns
- `entropy::Float64` - Von Neumann entropy, computed as ``−Tr(ρ·log(ρ))``. Zero for pure states.
- `purity::Float64` - Raw purity value: ``Tr(ρ²)``, where:
  - ``1/n`` indicates a maximally mixed state,
  - ``1`` indicates a pure state.
- `purityNorm::Float64` - Normalized purity in range ``[0, 1]``, where:
  - ``0`` indicates a maximally mixed state,
  - ``1`` indicates a pure state.
- `purityString::String` - Descriptive label of purity:
  - `"pure"` if ``ρ`` is approximately pure,
  - `"impure, totally mixed"` if ``ρ`` is close to ``(1/n)·I``,
  - `"impure, partially mixed"` otherwise.

# Purity Range
- Pure state: purity = ``1``
- Maximally mixed state: purity = ``1 / n``
- Normalized purity = ``(purity·n - 1) / (n - 1)``

# Notes
- Entropy uses eigenvalue decomposition: entropy = ``−∑ᵢ λᵢ·log(λᵢ)``, skipping ``λᵢ ≈ 0`` or ``1``.
- The function assumes ``ρ`` is a valid density matrix (Hermitian, positive semi-definite, unit trace).
- The purityString classification is tolerant to floating-point precision via accuracyCheckForPurity.

# See Also
- [`DensityState`](@ref)

# Examples

**Example 1**: Totally mixed density state
```julia-repl
julia> (entropy, purity, purityNorm, purityString)=calculateEntropyAndPurity([0.5 0.0; 0.0 0.5])
(0.6931471805599453, 0.5, 0.0, "impure, totally mixed")
```

**Example 2**: Pure density state
```julia-repl
julia> (entropy, purity, purityNorm, purityString)=calculateEntropyAndPurity([0.5 0.5; 0.5 0.5])
(0, 1.0, 1.0, "pure")
```
"""
function calculateEntropyAndPurity(rho; accuracyCheckForPurity = 10*eps(1.0))
	n = size(rho, 1)
	# accuracyCheckForPurity = 1e-6
	entropy = 0
	purity = 0
	L = eigvals(rho)
	for i in 1:n
		purity = purity + L[i]^2 # purity=Tr(rho^2)
		if ((L[i] != 0) && (L[i] != 1))
			entropy = entropy - L[i] * log(abs(L[i])) # entropy=-Tr(rho*ln(rho))
			# use limit x*log(x) --> 0 (x --> 0)
			# use limit x*log(x) --> 0 (x --> 1)
		end
	end
	# for finite system, entropy quantifies departure of system from a pure state. 
	# it codifies degree of mixing of given finite state. for pure state entropy==0
	# purity is in the range [1/n,...,1] where [1/n=impure totally mixed, impure partially mixed, 1=pure]
	purityNorm = (purity * n - 1) / (n - 1) # purityNorm(purity)=0+((1-0)/(1-1/n))*(purity-1/n)
	if abs(purity - (1 / n)) < accuracyCheckForPurity
		purityString = "impure, totally mixed"
	elseif abs(purity - 1) < accuracyCheckForPurity
		purityString = "pure"
	else
		purityString = "impure, partially mixed"
	end
	return real(entropy), real(purity), real(purityNorm), purityString
end

# """ 
#     calculateExpectationValueOfQubitState(q) -> Vector{Float64}

# Computes the expectation value (probability distribution) of a qubit state ``|q⟩``.

# This function returns a vector where each element represents the probability of measuring the qubit in the corresponding basis state. The probabilities are computed as the squared magnitudes of the complex amplitudes in ``|q⟩``.

# # Arguments
# - `q::AbstractVector` - A normalized qubit state vector ``|q⟩`` (typically of complex numbers).

# # Returns
# - `Vector{Float64}` - A vector of probabilities corresponding to each basis state.

# # Example
# ```julia-repl
# julia> calculateExpectationValueOfQubitState([1; 1]/sqrt(2))
# 2-element Vector{Float64}:
#  0.4999999999999999
#  0.4999999999999999
# ```
# """
function calculateExpectationValueOfQubitState(q)
	return abs.(q).^2
end

"""
    calculateCorrelation(
		x1::AbstractVector, 
		x2::AbstractVector
	) -> Float64

Compute the Pearson correlation coefficient between two vectors `x1` and `x2`.

This function measures the linear correlation between `x1` and `x2`. A value of ``1.0`` indicates perfect positive linear correlation, ``-1.0`` indicates perfect negative linear correlation, and ``0.0`` indicates no linear correlation.

# Arguments
- `x1::AbstractVector` - The first data vector.
- `x2::AbstractVector` - The second data vector. Must be the same length as `x1`.

# Returns
- `Float64` - The Pearson correlation coefficient between `x1` and `x2`.

# Examples
```julia-repl
julia> calculateCorrelation([1 2 3], [3 2 1])
-1.0
```
"""
function calculateCorrelation(x1, x2)
	xa1 = x1.-sum(x1)/length(x1)
	xa2 = x2.-sum(x2)/length(x2)
	sdxa12 = sqrt(sum(xa1.^2)*sum(xa2.^2))
	correlation = sum(xa1.*xa2)/sdxa12
	return correlation
end

"""
    @enum UseFromMeasurement useProbability useOutcome

Specify the source of data for computing expectation values from a quantum measurement.

# Variants
- `useProbability` - Use the theoretical probability distribution from the quantum measurement.
- `useOutcome` - Use the actual sampled outcomes (e.g., from repeated measurements or shots).

This enum is used to control whether to compute expectation values based on ideal probabilities or empirical results.
"""
@enum UseFromMeasurement useProbability useOutcome

"""
    calculateExpectationValueOfProductOfMeasuredQubits(
		quantumOutput::QuantumOutput, 
		stateId::Int, 
		useFromMeasurement::UseFromMeasurement
	) -> Float64

Computes the expectation value of the **product** of measured qubits from the given quantum simulation output.

The expectation value is computed by checking the parity of the measurement outcomes:
- If the number of `1`s is even, contribute `+1`.
- If the number of `1`s is odd, contribute `-1`.

# Parity Evaluation Example (2-qubit system)

| Bitstring | Number of 1s | Parity | Contribution |
| --------- | ------------ | ------ | ------------ |
| `00`      | 0            | Even   | `+1`         |
| `01`      | 1            | Odd    | `-1`         |
| `10`      | 1            | Odd    | `-1`         |
| `11`      | 2            | Even   | `+1`         |

# Arguments
- `quantumOutput::QuantumOutput` - The result of running a quantum program.
- `stateId::Int` - The index of the quantum state to extract results from.
- `useFromMeasurement::UseFromMeasurement` - Determines whether to use probabilities (`useProbability`) or actual measurement outcomes (`useOutcome`) for the computation.

# Returns
- `Float64` - The expectation value of the product of the measured qubit values.

# See Also
- [`QuantumOutput`](@ref)
- [`UseFromMeasurement`](@ref)

# Examples

**Uncorrelated 2-qubit superposition state**

Each qubit is put in superposition using a Hadamard gate: ``|ψ⟩ = (|0⟩ + |1⟩)/√2 ⊗ (|0⟩ + |1⟩)/√2 = (1/2)(|00⟩ + |01⟩ + |10⟩ + |11⟩)``.

All bitstrings occur with equal probability ``0.25``. Expectation value is:

| Bitstring | Number of 1s | Parity | Contribution | Probability | Weighted Contribution |
| --------- | ------------ | ------ | ------------ | ----------- | --------------------- |
| `00`      | 0            | Even   | `+1`         | 0.25        | `+0.25`               |
| `01`      | 1            | Odd    | `-1`         | 0.25        | `-0.25`               |
| `10`      | 1            | Odd    | `-1`         | 0.25        | `-0.25`               |
| `11`      | 2            | Even   | `+1`         | 0.25        | `+0.25`               |

Final Expectation Value: ``0.25 − 0.25 − 0.25 + 0.25 = 0.0``

```julia-repl
julia> qc = createQuantumCircuit(2)
julia> hGate!(qc, 1)
julia> hGate!(qc, 2)
julia> measureGate!(qc, [1, 2])
julia> qp = compileQuantumCircuit(qc; optimizeNumberOfSteps = true)
julia> iqs = createInitialQubitState(vector, [1. 0.; 1. 0.])
julia> qo = runQuantumProgram(qp, iqs, 1)
julia> corr = calculateExpectationValueOfProductOfMeasuredQubits(qo, 3, useProbability)
0.0
```

**Correlated 2-qubit Bell state**

We prepare an entangled Bell state: ``|ψ⟩ = (|00⟩ + |11⟩) / √2``

Only `00` and `11` occur, each with probability ``0.5``. Expectation value is:

| Bitstring | Number of 1s | Parity | Contribution | Probability | Weighted Contribution |
| --------- | ------------ | ------ | ------------ | ----------- | --------------------- |
| `00`      | 0            | Even   | `+1`         | 0.5         | `+0.5`                |
| `01`      | 1            | Odd    | `-1`         | 0.0         | ` 0.0`                |
| `10`      | 1            | Odd    | `-1`         | 0.0         | ` 0.0`                |
| `11`      | 2            | Even   | `+1`         | 0.5         | `+0.5`                |

Final Expectation Value: ``0.5 + 0.0 + 0.0 + 0.5 = 1.0``

```julia-repl
julia> qc = createQuantumCircuit(2)
julia> hGate!(qc, 1)
julia> cnotGate!(qc, 1, 2)
julia> measureGate!(qc, [1, 2])
julia> qp = compileQuantumCircuit(qc; optimizeNumberOfSteps = true)
julia> iqs = createInitialQubitState(vector, [1. 0.; 1. 0.])
julia> qo = runQuantumProgram(qp, iqs, 1)
julia> corr = calculateExpectationValueOfProductOfMeasuredQubits(qo, 3, useProbability)
1.0
```
"""
function calculateExpectationValueOfProductOfMeasuredQubits(quantumOutput::QuantumOutput, stateId::Int, useFromMeasurement::UseFromMeasurement)
    nrMeasuredBits = Int(log2(length(quantumOutput.output[stateId].measured.probability)))
    expVal = 0.
    if useFromMeasurement == useProbability
        for m in 1:2^nrMeasuredBits
            bin = reverse(digits(m-1, base = 2, pad = nrMeasuredBits)')
            if iseven(sum(bin))
                expVal += quantumOutput.output[stateId].measured.probability[m]
            else
                expVal += -quantumOutput.output[stateId].measured.probability[m]
            end
        end
    else # useOutcome
        nrShots = size(quantumOutput.output, 2)
        for k in 1:nrShots
            bin = reverse(digits(quantumOutput.output[stateId, k].measured.outcome, base = 2, pad = nrMeasuredBits)')
            if iseven(sum(bin))
                expVal += 1.
            else
                expVal += -1.
            end
        end
        expVal = expVal/nrShots
    end
    return expVal
end

""" 
    partialTrace(
		rho::Matrix{<:Complex}; 
		traceIndex::Int, 
		zeroBasedNumbering = false, 
		indexType = circuitIndexBigEndian
	) -> Matrix{<:Complex}

Computes the partial trace of a given density matrix ``ρ`` over a specified subsystem index `traceIndex`. 
This operation removes (traces out) one qubit from a larger system, resulting in a reduced density matrix for the remaining subsystem(s).

# Arguments
- `rho::Matrix{<:Complex}` - Density matrix ``ρ`` representing a quantum state (Hermitian, positive semi-definite, trace 1).
- `traceIndex::Int` — The qubit index to trace out, interpreted using the settings for `zeroBasedNumbering` and `indexType`.
- `zeroBasedNumbering::Bool` (optional) — Whether `traceIndex` uses zero-based indexing. Defaults to `false`.
- `indexType::IndexType` (optional) — Interpretation of `traceIndex` based on circuit/byte ordering. Defaults to `circuitIndexBigEndian`.

# Returns
- `Matrix{<:Complex}` — Reduced density matrix.

# See Also
- [`IndexType`](@ref)

# Example
Trace out the first qubit from the Bell state ``|ψ⟩ = (|00⟩ + |11⟩) / √2``:
```julia-repl
julia> partialTrace([1 0 0 1; 0 0 0 0; 0 0 0 0; 1 0 0 1]/2., traceIndex = 1)
2×2 Matrix{Float64}:
 0.5  0.0
 0.0  0.5
```
"""
function partialTrace(rho; traceIndex::Int, zeroBasedNumbering = false, indexType = circuitIndexBigEndian)
	numberOfQubits = Int(log2(size(rho, 1)))

	# assume traceIndex conforms to Settings
	traceIndex = convertToByteIndex(Settings(zeroBasedNumbering, indexType), numberOfQubits, traceIndex)

	list0 = zeros(Int, 2^(numberOfQubits-1), 1)
	list1 = zeros(Int, 2^(numberOfQubits-1), 1)
	p = 1
	q = 1
	for m in 1:2^numberOfQubits
		bin = reverse(digits(m-1, base = 2, pad = numberOfQubits)')
		if bin[traceIndex] == 0
			list0[p] = m
			p += 1
		else
			list1[q] = m
			q += 1
		end
	end

	T0 = zeros(Int, 2^(numberOfQubits-1), 2^numberOfQubits)
	T1 = zeros(Int, 2^(numberOfQubits-1), 2^numberOfQubits)
	for m in 1:2^(numberOfQubits-1)
		T0[m, list0[m]] = 1
		T1[m, list1[m]] = 1
	end

	return T0 * rho * T0' + T1 * rho * T1'
end

""" 
    fidelity(
		rho::Matrix{<:Complex}, 
		sigma::Matrix{<:Complex}
	) -> Float64

Compute the **quantum fidelity** between two density matrices ``ρ`` and ``σ``.

Fidelity is a measure of similarity between quantum states, defined as: ``F(ρ, σ) = Tr[√ρ·σ·√ρ]``

# Arguments
- `rho::Matrix{<:Complex}` - Density matrix ``ρ`` representing a quantum state (Hermitian, positive semi-definite, trace 1).
- `sigma::Matrix{<:Complex}` - Another density matrix ``σ`` representing a quantum state of the same size as ``ρ``.

# Returns
- `Float64` - The fidelity between ``ρ`` and ``σ``, a value between 0 and 1

# See Also
- [`DensityState`](@ref)

# Example
Compute the fidelity between the pure states ``ρ = |00⟩⟨00|`` and ``σ = |11⟩⟨11|``:
```julia-repl
julia> fidelity([1 0; 0 0], [0 0; 0 1])
0.0
```
"""
function fidelity(rho, sigma)
	return tr(sqrt(rho) * sigma * sqrt(rho))
end