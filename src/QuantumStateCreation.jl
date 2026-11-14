############################################
# QuantumStateCreation

""" 
    createSingleQubitBlochDensityState(
		n::Vector{Float64}
	) -> DensityState
	
	createSingleQubitBlochDensityState(
		theta::Float64, 
		phi::Float64
	) -> DensityState

Creates a single qubit **Bloch density state** in the direction specified either by a 3D unit vector ``n`` or by spherical polar angles ``θ`` and ``φ``.

The Bloch density state is given by: 
```math
ρ = 0.5 * (I + n ⋅ σ) = 0.5 * (I + n_x ⋅ σ_x + n_y ⋅ σ_y + n_z ⋅ σ_z)
```

where:
- ``I`` is the 2×2 identity matrix,
- ``n = [n_x, n_y, n_z]`` is a real 3D unit vector representing a point on the Bloch sphere,
- ``n = [sin(θ)⋅cos(φ), sin(θ)⋅sin(φ), cos(θ)]`` when parameterized by spherical polar angles,
- ``σ = [σ_x, σ_y, σ_z]`` is the vector of Pauli matrices associated with the X, Y, and Z axes, respectively.

# Arguments
- `n::Vector{Float64}` - A 3-dimensional unit vector ``n = [n_x, n_y, n_z]`` specifying a direction on the Bloch sphere.
- `theta::Float64` - Polar angle ``θ`` in radians (angle from the Z-axis).
- `phi::Float64` - Azimuthal angle ``φ`` in radians (angle from the X-axis in the XY-plane).

# Returns
- `DensityState` - A single-qubit **Bloch density state** ρ aligned with the specified Bloch vector direction.

# Notes
- The resulting matrix ``ρ`` is Hermitian, positive semi-definite, and has trace 1.
- This representation is useful for simulating qubit states in both pure and mixed form.
- When ``‖n‖ = 1``, the state is pure. When ``‖n‖ < 1``, the state is mixed.

# Calling structure
```julia-repl
createSingleQubitBlochDensityState(theta::Float64, phi::Float64)
    └── createSingleQubitBlochDensityState(n::Vector{Float64})
        └── DensityState(...)
```

# See Also
- [`pauliX`](@ref)
- [`pauliY`](@ref)
- [`pauliZ`](@ref)
- [`DensityState`](@ref)

# Examples

**Example 1: Using a normalized Bloch vector**

```julia-repl
julia> createSingleQubitBlochDensityState(vec([1 1 1])/sqrt(3)).rho
2×2 Matrix{ComplexF64}:
 0.788675+0.0im       0.288675-0.288675im
 0.288675+0.288675im  0.211325+0.0im
```

**Example 2: Using spherical coordinates ``θ = π/3``, ``φ = π/4``**

```julia-repl
julia> createSingleQubitBlochDensityState(pi/3., pi/4.).rho
2×2 Matrix{ComplexF64}:
     0.75+0.0im       0.306186-0.306186im
 0.306186+0.306186im      0.25+0.0im
```
"""
function createSingleQubitBlochDensityState(n::Vector{Float64})
	# rho = 0.5*(identity+n*sigma) with sigma the 3 pauli matrices
	# rho = 0.5*(eye(2)+n(1)*[[0,1];[1,0]]+n(2)*[[0,-i];[i,0]]+n(3)*[[1,0];[0,-1]]);
	return DensityState(0.5*[1+n[3] n[1]-im*n[2]; n[1]+im*n[2] 1-n[3]])
end

function createSingleQubitBlochDensityState(theta::Float64, phi::Float64)
	return createSingleQubitBlochDensityState(vec([sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)]))
end

"""
    @enum StateSpace vector density

Supported quantum mechanical state representations:

- `vector`: **State vector** representation ``|Ψ⟩``
  Represents pure quantum states as complex-valued vectors in Hilbert space.

- `density`: **Density matrix** representation ``ρ``
  Represents statistical mixtures of pure states and supports both pure and mixed quantum states via a positive semi-definite matrix with unit trace.
"""
@enum StateSpace vector density

"""
    createByteIndexVector(
		value::Int, 
		numberOfQubits::Int
	) -> Matrix{Float64}

Constructs a **multi-qubit product vector state** from an integer `value` for a system of `numberOfQubits` qubits.

The resulting quantum vector state has the form: ``|ψ⟩ = |ψ₁⟩ ⊗ |ψ₂⟩ ⊗ ⋯ ⊗ |ψₙ⟩`` where each individual qubit is in the basis vector state: ``|ψⱼ⟩ = αⱼ|0⟩ + βⱼ|1⟩``

The function returns a 2×`numberOfQubits` matrix of the form: ``[α_1 ... α_n; β_1 ... β_n]``.

Each column represents one qubit, and is initialized to:
- ``|0⟩ = [1; 0]`` if the corresponding bit in the binary representation of `value` is 0
- ``|1⟩ = [0; 1]`` if the bit is 1

The least significant bit of `value` corresponds to the last qubit (rightmost column), consistent with standard quantum computing conventions.

This representation is useful as input for `createInitialQubitState` when using `indexType = byteIndex`.

# Arguments
- `value::Integer` - Integer value to convert into binary qubit states.
- `numberOfQubits::Integer` - Number of qubits to encode.

# Returns
- `Matrix{Float64}` - A 2×`numberOfQubits` matrix representing the qubit product state.

# See Also
- [`IndexType`](@ref)
- [`createInitialQubitState`](@ref)

# Example
Create a 3-qubit state from the value 6 (binary 110), producing ``|ψ⟩ = |1⟩ ⊗ |1⟩ ⊗ |0⟩``:
```julia-repl
julia> createByteIndexVector(6, 3)
2×3 Matrix{Float64}:
 0.0  0.0  1.0
 1.0  1.0  0.0
```
"""
function createByteIndexVector(value::Int, numberOfQubits::Int)
    valueBinary = reverse(digits(value, base = 2, pad = numberOfQubits)')
    byteIndexVector = zeros(2, numberOfQubits)
    for k in 1:numberOfQubits
        if valueBinary[k] == 0
            byteIndexVector[1, k] = 1.
        else
            byteIndexVector[2, k] = 1.
        end
    end
    return byteIndexVector
end

"""
	createInitialQubitState(
		stateSpace::StateSpace, 
		value; 
		indexType = circuitIndexBigEndian, 
		blochRepresentation = false
	) -> State

Creates a quantum state (`State`) in the given `stateSpace` using the provided `value`. This function supports both pure states (`vector`)
and mixed states (`density`), and allows initialization using either amplitudes or Bloch sphere parameters.

# Arguments
- `stateSpace::StateSpace` - The target quantum state space. Must be either `vector` (for pure states) or `density` (for mixed states).
- `value` - Array or structured data describing the initial qubit state(s). Format depends on `stateSpace` and `blochRepresentation`.
- `indexType::IndexType` (optional) - Determines the qubit indexing convention. Defaults to `circuitIndexBigEndian`.
- `blochRepresentation::Bool` (optional) - If `true`, interprets values using Bloch sphere angles ``(θ, ϕ)``. Defaults to `false`.

# Returns
- `State`: A `VectorState` or `DensityState` depending on the chosen `stateSpace`.

# Pure States (`stateSpace` = `vector`)

## Representations

- **Amplitude representation** (default):
  A single-qubit state is defined as: ``|ψ⟩ = α|0⟩ + β|1⟩``
  For multiple qubits, the full state is constructed as: ``|ψ⟩ = |ψ₁⟩ ⊗ |ψ₂⟩ ⊗ ... ⊗ |ψₙ⟩``

- **Bloch representation** (`blochRepresentation = true`):
  Each qubit is defined using Bloch sphere angles: ``|ψ⟩ = cos(θ/2)|0⟩ + sin(θ/2)·exp(i·ϕ)|1⟩``

## Accepted Formats
```julia-repl
createInitialQubitState(vector, [α₁ β₁; ... ; αₙ βₙ])
createInitialQubitState(vector, [αₙ βₙ; ... ; α₁ β₁], indexType=circuitIndexLittleEndian)
createInitialQubitState(vector, [α₁ ... αₙ; β₁ ... βₙ], indexType=byteIndex)

createInitialQubitState(vector, [θ₁ ϕ₁; ... ; θₙ ϕₙ], blochRepresentation=true)
createInitialQubitState(vector, [θₙ ϕₙ; ... ; θ₁ ϕ₁], indexType=circuitIndexLittleEndian, blochRepresentation=true)
createInitialQubitState(vector, [θ₁ ... θₙ; ϕ₁ ... ϕₙ], indexType=byteIndex, blochRepresentation=true)
```

# Mixed States (`stateSpace` = `density`)

## Representations

- **Amplitude representation** (default):
  Mixed states are convex combinations of pure states: ``ρ = p₁|ψ₁⟩⟨ψ₁| + p₂|ψ₂⟩⟨ψ₂| + ... + pₘ|ψₘ⟩⟨ψₘ|``
  Each ``|ψₖ⟩`` is built from individual qubits: ``|ψₖ⟩ = |ψₖ₁⟩ ⊗ |ψₖ₂⟩ ⊗ ... ⊗ |ψₖₙ⟩``
  With each qubit: ``|ψₖⱼ⟩ = αₖⱼ|0⟩ + βₖⱼ|1⟩``

- **Bloch representation** (`blochRepresentation = true`):
  Each qubit in each component state is defined using: ``|ψₖⱼ⟩ = cos(θₖⱼ/2)|0⟩ + sin(θₖⱼ/2)·exp(i·ϕₖⱼ)|1⟩``

## Accepted Formats
```julia-repl
createInitialQubitState(density, [[p₁, [α₁₁ β₁₁; ... ; α₁ₙ β₁ₙ]], ..., [pₘ, [αₘ₁ βₘ₁; ... ; αₘₙ βₘₙ]]])
createInitialQubitState(density, [[p₁, [α₁ₙ β₁ₙ; ... ; α₁₁ β₁₁]], ..., [pₘ, [αₘₙ βₘₙ; ... ; αₘ₁ βₘ₁]]], indexType=circuitIndexLittleEndian)
createInitialQubitState(density, [[p₁, [α₁₁ ... α₁ₙ; β₁₁ ... β₁ₙ]], ..., [pₘ, [αₘ₁ ... αₘₙ; βₘ₁ ... βₘₙ]]], indexType=byteIndex)

createInitialQubitState(density, [[p₁, [θ₁₁ ϕ₁₁; ... ; θ₁ₙ ϕ₁ₙ]], ..., [pₘ, [θₘ₁ ϕₘ₁; ... ; θₘₙ ϕₘₙ]]], blochRepresentation=true)
createInitialQubitState(density, [[p₁, [θ₁ₙ ϕ₁ₙ; ... ; θ₁₁ ϕ₁₁]], ..., [pₘ, [θₘₙ ϕₘₙ; ... ; θₘ₁ ϕₘ₁]]], indexType=circuitIndexLittleEndian, blochRepresentation=true)
createInitialQubitState(density, [[p₁, [θ₁₁ ... θ₁ₙ; ϕ₁₁ ... ϕ₁ₙ]], ..., [pₘ, [θₘ₁ ... θₘₙ; ϕₘ₁ ... ϕₘₙ]]], indexType=byteIndex, blochRepresentation=true)
```

# See Also
- [`IndexType`](@ref)
- [`StateSpace`](@ref)
- [`State`](@ref)
- [`VectorState`](@ref)
- [`DensityState`](@ref)

# Examples

**Pure States (vector)**

**Example 1: Big-endian, ``|ψ⟩ = |0⟩⊗|1⟩ → [0; 1; 0; 0]``**
```julia-repl
julia> createInitialQubitState(vector, [1. 0.; 0. 1.])
```
**Example 2: Little-endian, same state**
```julia-repl
julia> createInitialQubitState(vector, [0. 1.; 1. 0.], indexType=circuitIndexLittleEndian)
```
**Example 3: Byte index, ``|ψ⟩ = |0⟩⊗|1⟩⊗|1⟩ → [0; 0; 0; 1; 0; 0; 0; 0]``**
```julia-repl
julia> createInitialQubitState(vector, [1. 0. 0.; 0. 1. 1.], indexType=byteIndex)
```
**Example 4: Byte index for 3 qubits, ``|ψ⟩ = |6⟩ = |1⟩⊗|1⟩⊗|0⟩ → [0; 0; 0; 0; 0; 0; 1; 0]``**
```julia-repl
julia> createInitialQubitState(vector, createByteIndexVector(6, 3), indexType=byteIndex)
```
**Example 5: Big-endian, Bloch representation: ``|ψ⟩ = |+⟩⊗|–⟩⊗|1⟩``**
```julia-repl
julia> createInitialQubitState(vector, [π/2 0.; π/2 π/2], blochRepresentation=true)
```

**Mixed States (density)**

**Example 7: Big-endian, ``ρ = 0.25|+⟩⟨+| + 0.75|–⟩⟨–| → [0.5, -0.25; -0.25, 0.5]``**
```julia-repl
julia> createInitialQubitState(density, [[0.25, [1. 1.]/sqrt(2)], [0.75, [1. -1.]/sqrt(2)]])
```
**Example 8: Same using Bloch representation**
```julia-repl
julia> createInitialQubitState(density, [[0.25, [π/2 0.]], [0.75, [π/2 π]]], blochRepresentation=true)
```

# Notes
The qubit order interpretation depends on indexType. Use this to align with simulation or hardware layout.
All input vectors/matrices must be properly normalized or interpreted as unnormalized (e.g., Bloch angles).
"""
createInitialQubitState(stateSpace::StateSpace, value; indexType = circuitIndexBigEndian, blochRepresentation = false) = dispatchCreateInitialQubitState(Val(stateSpace), value, Val(indexType), Val(blochRepresentation))

dispatchCreateInitialQubitState(::Val{vector}, ab, indexType, blochRepresentation) = VectorState(buildStateVector(ab, indexType, blochRepresentation))
dispatchCreateInitialQubitState(::Val{density}, pab, indexType, blochRepresentation) = DensityState(buildDensityMatrix(pab, indexType, blochRepresentation))

function buildStateVector(ab, indexType, blochRepresentation)
    n = getQubitCount(ab, indexType)
    q = ones(ComplexF64, 2^n, 1)
    for r in 1:(2^n)
        bits = getBitPattern(r, n, indexType)
        q[r] = computeAmplitude(bits, ab, indexType, blochRepresentation)
    end
    return q
end

function buildDensityMatrix(pab, indexType, blochRepresentation)
    n = getQubitCount(pab[1][2], indexType)
    rho = zeros(ComplexF64, 2^n, 2^n)
    for (weight, ab) in pab
        q = buildStateVector(ab, indexType, blochRepresentation)
        rho += weight * q * q'
    end
    return rho
end

function computeAmplitude(bits, ab, indexType, blochRepresentation)
    amp = one(ComplexF64)
    for c in eachindex(bits)
        a, b = getAmplitudePair(ab, c, indexType)
        amp *= bits[c] == 0 ? firstComponent(a, b, blochRepresentation) : secondComponent(a, b, blochRepresentation)
    end
    return amp
end

firstComponent(a, _, ::Val{false}) = a
firstComponent(a, _, ::Val{true}) = cos(a / 2)
secondComponent(_, b, ::Val{false}) = b
secondComponent(a, b, ::Val{true}) = exp(im * b) * sin(a / 2)

getAmplitudePair(ab, c, ::Val{byteIndex}) = (ab[1, c], ab[2, c])
getBitPattern(r, n, ::Val{byteIndex}) = reverse(digits(r - 1, base = 2, pad = n)')
getQubitCount(ab, ::Val{byteIndex}) = size(ab, 2)

getAmplitudePair(ab, c, ::Val{circuitIndexBigEndian}) = (ab[c, 1], ab[c, 2])
getBitPattern(r, n, ::Val{circuitIndexBigEndian}) = reverse(digits(r - 1, base = 2, pad = n)')
getQubitCount(ab, ::Val{circuitIndexBigEndian}) = size(ab, 1)

getAmplitudePair(ab, c, ::Val{circuitIndexLittleEndian}) = (ab[c, 1], ab[c, 2])
getBitPattern(r, n, ::Val{circuitIndexLittleEndian}) = digits(r - 1, base = 2, pad = n)'
getQubitCount(ab, ::Val{circuitIndexLittleEndian}) = size(ab, 1)

"""
    createDoubleQubitWernerDensityState(
		p::Real; 
		indexType::IndexType = circuitIndexBigEndian
	) -> DensityState

Creates a **Werner density state** for two qubits with mixing parameter ``p``.

The Werner density state is defined as a mixture of a maximally entangled singlet state and the maximally mixed state: ``ρ = p·|Ψ⁻⟩⟨Ψ⁻| + (1 - p)·I/4``

where ``|Ψ⁻⟩ = (|01⟩ - |10⟩) / √2`` is the **singlet Bell state**, and ``I`` is the 4×4 identity matrix.

# Arguments
- `p::Real` - Mixing probability ``p``. Must be in the range ``[0, 1]``.
  - ``p = 1``: Fully entangled (pure singlet state).
  - ``p = 0``: Completely mixed state (maximally impure).
- `indexType::IndexType` (optional) - Determines qubit indexing convention. Defaults to `circuitIndexBigEndian`.

# Returns
- `DensityState` - A two-qubit mixed quantum state ``ρ`` represented as a `DensityState` object.

# Properties
- **Bell inequality** is violated if ``p > 1/√2``.
- **Positive partial transpose (PPT)** criterion indicates entanglement if ``p > 1/3``.

# See Also
- [`IndexType`](@ref)
- [`DensityState`](@ref)

# Examples

**Example 1**: Completely mixed state ``p = 0``
```julia-repl
julia> createDoubleQubitWernerDensityState(0).rho
4×4 Matrix{ComplexF64}:
 0.25+0.0im   0.0+0.0im   0.0+0.0im   0.0+0.0im
  0.0+0.0im  0.25+0.0im   0.0+0.0im   0.0+0.0im
  0.0+0.0im   0.0+0.0im  0.25+0.0im   0.0+0.0im
  0.0+0.0im   0.0+0.0im   0.0+0.0im  0.25+0.0im
```

**Example 2**: Pure singlet state ``p = 1``
```julia-repl
julia> createDoubleQubitWernerDensityState(1).rho
4×4 Matrix{ComplexF64}:
 0.0+0.0im   0.0+0.0im   0.0+0.0im  0.0+0.0im
 0.0+0.0im   0.5+0.0im  -0.5+0.0im  0.0+0.0im
 0.0+0.0im  -0.5+0.0im   0.5+0.0im  0.0+0.0im
 0.0+0.0im   0.0+0.0im   0.0+0.0im  0.0+0.0im
```
"""
function createDoubleQubitWernerDensityState(p::Real; indexType::IndexType = circuitIndexBigEndian)
	# rho=p*|singlet><singlet|+0.25*(1-p)*identity, here we choose |singlet>=(|01>-|10>)/sqrt(2)
	# 1=pure entangled, 0=impure, Bell's inequality is satisfied if p<1/sqrt(2), 
	# PPT condition is satisfied if p<1/3
	return DensityState(p * [0 0 0 0; 0 1 -1 0; 0 -1 1 0; 0 0 0 0] / 2 + (1-p) * [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1] / 4)
end

