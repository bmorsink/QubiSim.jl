using QubiSim
using LinearAlgebra
using Test

function toggleQubitOrdeningFor(A::UnitaryOperation)
    swapList = createToggleSwapList(Int(log2(size(A.U, 1))))
	return UnitaryOperation(A.U[swapList, swapList])
end

# what to do with Measured part of the VectorState?
function toggleQubitOrdeningFor(A::VectorState)
    swapList = createToggleSwapList(Int(log2(size(A.q, 1))))
	# return VectorState(reshape(A.q[swapList],size(A.q, 1),1), Measured(Vector{Float64}(),-1,Matrix{Int64}(undef,0,0),Vector{String}())) # reshape needed to convert Vector to Matrix
	# return VectorState(reshape(A.q[swapList],size(A.q, 1),1), nothing) # reshape needed to convert Vector to Matrix
	return VectorState(reshape(A.q[swapList],size(A.q, 1),1); measured = nothing) # reshape needed to convert Vector to Matrix
	# return VectorState(reshape(A.q[swapList],size(A.q, 1),1), []) # reshape needed to convert Vector to Matrix
end

# what to do with Measured part of the DensityState?
function toggleQubitOrdeningFor(A::DensityState)
    swapList = createToggleSwapList(Int(log2(size(A.rho, 1))))
	# return DensityState(A.rho[swapList, swapList], [], [], [], Measured(Vector{Float64}(),-1,Matrix{Int64}(undef,0,0),Vector{String}()))
	# return DensityState(A.rho[swapList, swapList], [], [], [], nothing)
	return DensityState(A.rho[swapList, swapList]; measured = nothing)
	# return DensityState(A.rho[swapList, swapList], [], [], [], [])
end

@testset "testing" verbose = true begin

    @testset "Qubit index conversions" verbose = true begin
        @testset "convertToOneBasedNumbering" begin
            @test QubiSim.convertToOneBasedNumbering(Settings(true, circuitIndexBigEndian), 0) == 1
            @test QubiSim.convertToOneBasedNumbering(Settings(true, circuitIndexBigEndian), 3) == 4
            @test QubiSim.convertToOneBasedNumbering(Settings(false, circuitIndexBigEndian), 1) == 1
            @test QubiSim.convertToOneBasedNumbering(Settings(false, circuitIndexLittleEndian), 5) == 5
        end

        @testset "convertFromOneToSettingsBasedNumbering" begin
            @test QubiSim.convertFromOneToSettingsBasedNumbering(Settings(true, circuitIndexLittleEndian), 1) == 0
            @test QubiSim.convertFromOneToSettingsBasedNumbering(Settings(true, circuitIndexLittleEndian), 5) == 4
            @test QubiSim.convertFromOneToSettingsBasedNumbering(Settings(false, byteIndex), 1) == 1
            @test QubiSim.convertFromOneToSettingsBasedNumbering(Settings(false, byteIndex), 100) == 100
        end

        @testset "convertToByteIndex" begin
            # Big-endian: should match one-based index
            @test QubiSim.convertToByteIndex(Settings(false, circuitIndexBigEndian), 4, 2) == 2
            @test QubiSim.convertToByteIndex(Settings(true, byteIndex), 4, 1) == 2  # index=1 → one-based=2

            # Little-endian: index should flip
            @test QubiSim.convertToByteIndex(Settings(false, circuitIndexLittleEndian), 4, 1) == 4
            @test QubiSim.convertToByteIndex(Settings(false, circuitIndexLittleEndian), 4, 2) == 3
            @test QubiSim.convertToByteIndex(Settings(true, circuitIndexLittleEndian), 4, 0) == 4  # 0 → 1 → flipped: 4
            @test QubiSim.convertToByteIndex(Settings(true, circuitIndexLittleEndian), 4, 1) == 3  # 1 → 2 → flipped: 3
        end
    end

    @testset "getStep and getGate" verbose = true begin
        qc = createQuantumCircuit(2)
        hGate!(qc, 2)
        xGate!(qc, 1)
        cnotGate!(qc, 1, 2)
        measureGate!(qc,[1,2])
        quantumChannelGate!(qc, [1,2], KrausOperators())
        measureGate!(qc,[1,2],KrausOperators())

        @test getGate(qc, getStep(qc, 1), 1) isa UnitaryGate
        @test getGate(qc, getStep(qc, 1), 1).unitaryOperationFactory() isa UnitaryOperation
        @test getGate(qc, getStep(qc, 1), 1).unitaryOperationFactory().U ≈ [0.7071067811865476 + 0.0im 0.7071067811865475 - 8.659560562354932e-17im; 0.7071067811865475 + 0.0im -0.7071067811865476 + 8.659560562354934e-17im] atol=0.0001
        @test getGate(qc, getStep(qc, 1), 1).qubits == [2]
        @test getGate(qc, getStep(qc, 1), 1).name == "H"

        @test getGate(qc, getStep(qc, 1), 2) isa UnitaryGate
        @test getGate(qc, getStep(qc, 1), 2).unitaryOperationFactory() isa UnitaryOperation
        @test getGate(qc, getStep(qc, 1), 2).unitaryOperationFactory().U ≈ [6.123233995736766e-17 + 0.0im 1.0 - 1.2246467991473532e-16im; 1.0 + 0.0im -6.123233995736766e-17 + 7.498798913309288e-33im] atol=0.0001
        @test getGate(qc, getStep(qc, 1), 2).qubits == [1]
        @test getGate(qc, getStep(qc, 1), 2).name == "X"

        @test getGate(qc, getStep(qc, 2), 1) isa UnitaryGate
        @test getGate(qc, getStep(qc, 2), 1).unitaryOperationFactory() isa UnitaryOperation
        @test getGate(qc, getStep(qc, 2), 1).unitaryOperationFactory().U ≈ [1.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 1.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 1.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 1.0 + 0.0im 0.0 + 0.0im] atol=0.0001
        @test getGate(qc, getStep(qc, 2), 1).qubits == [1, 2]
        @test getGate(qc, getStep(qc, 2), 1).name == "CNOT"

        @test getGate(qc, getStep(qc, 3), 1) isa MeasureGate
        @test getGate(qc, getStep(qc, 3), 1).measureGateType == QubiSim.PVM
        @test getGate(qc, getStep(qc, 3), 1).qubits == [1, 2]
        @test getGate(qc, getStep(qc, 3), 1).sigmas == Matrix{ComplexF64}[]
        @test isnothing(getGate(qc, getStep(qc, 3), 1).krausOperators)
        @test getGate(qc, getStep(qc, 3), 1).forgetOutcome == false

        @test getGate(qc, getStep(qc, 4), 1) isa QuantumChannelGate
        @test getGate(qc, getStep(qc, 4), 1).qubits == [1, 2]
        @test getGate(qc, getStep(qc, 4), 1).krausOperators == KrausOperator[]

        @test getGate(qc, getStep(qc, 5), 1) isa MeasureGate
        @test getGate(qc, getStep(qc, 5), 1).measureGateType == QubiSim.POVM
        @test getGate(qc, getStep(qc, 5), 1).qubits == [1, 2]
        @test isnothing(getGate(qc, getStep(qc, 5), 1).sigmas)
        @test getGate(qc, getStep(qc, 5), 1).krausOperators == KrausOperator[]
        @test getGate(qc, getStep(qc, 5), 1).forgetOutcome == false
    end

    @testset "getOperation" verbose = true begin
        qc = createQuantumCircuit(1)
        hGate!(qc, 1)
        xGate!(qc, 1)
        measureGate!(qc,[1])
        quantumChannelGate!(qc, [1], KrausOperators())
        measureGate!(qc,[1], forgetOutcome=true)
        qp = compileQuantumCircuit(qc)

        @test getOperation(qp, 1) isa UnitaryOperation
        @test getOperation(qp, 1).U ≈ [0.707107+0.0im   0.707107+0.0im; 0.707107+0.0im  -0.707107+0.0im] atol=0.0001

        @test getOperation(qp, 2) isa UnitaryOperation
        @test getOperation(qp, 2).U ≈ [ 0.0+0.0im   1.0+0.0im; 1.0+0.0im  0.0+0.0im] atol=0.0001

        @test getOperation(qp, 3) isa MeasureOperation
        # @test getOperation(qp, 3).forgetOutcome == false
        @test getOperation(qp, 3).krausOperators isa KrausOperators

        @test getOperation(qp, 4) isa QuantumChannelOperation
        @test getOperation(qp, 4).krausOperators isa KrausOperators

        @test getOperation(qp, 5) isa MeasureAndForgetOperation
        @test getOperation(qp, 5).krausOperators isa KrausOperators
    end

    @testset "getState, getShot, getStateShot" verbose = true begin
        qc = createQuantumCircuit(1)
        hGate!(qc, 1)
        xGate!(qc, 1)
        qp = compileQuantumCircuit(qc)
        iqs = createInitialQubitState(vector, [1. 0.;])
        qo = runQuantumProgram(qp, iqs, 5)

        @test length(getState(qo, 3)) == 5
        @test getState(qo, 3)[1] isa VectorState
        @test getState(qo, 3)[1].q ≈ [0.7071067811865475 + 0.0im; 0.7071067811865476 + 0.0im;;] atol=0.0001

        @test length(getShot(qo, 4)) == 3
        @test getShot(qo, 4)[1] isa VectorState
        @test getShot(qo, 4)[1].q ≈ [1.0 + 0.0im; 0.0 + 0.0im;;] atol=0.0001
        @test getShot(qo, 4)[2] isa VectorState
        @test getShot(qo, 4)[2].q ≈ [0.7071067811865476 + 0.0im; 0.7071067811865475 + 0.0im;;] atol=0.0001
        @test getShot(qo, 4)[3] isa VectorState
        @test getShot(qo, 4)[3].q ≈ [0.7071067811865475 + 0.0im; 0.7071067811865476 + 0.0im;;] atol=0.0001

        @test getStateShot(qo, 3, 4) isa VectorState
        @test getStateShot(qo, 3, 4).q ≈ [0.7071067811865475 + 0.0im; 0.7071067811865476 + 0.0im;;] atol=0.0001
    end

    @testset "struct (Settings, QuantumCircuit, QuantumProgram and QuantumOutput) type checking" verbose = true begin
        qc=createQuantumCircuit(2)
        hGate!(qc, 1)
        cnotGate!(qc, 1, 2)
        measureGate!(qc, [1,2])
        @test typeof(qc.settings) == Settings{false, circuitIndexBigEndian}
        @test typeof(qc) == QuantumCircuit{false, circuitIndexBigEndian}

        qp=compileQuantumCircuit(qc)
        @test typeof(qp) == QuantumProgram{false, circuitIndexBigEndian}
        @test typeof(qp.settings) == Settings{false, circuitIndexBigEndian}

        iqs=createInitialQubitState(vector, [1. 0.;1. 0.])
        qo=runQuantumProgram(qp, iqs, 2)
        @test typeof(qo) == QuantumOutput{false, circuitIndexBigEndian, VectorState}
        @test typeof(qo.settings) == Settings{false, circuitIndexBigEndian}
        @test typeof(qo.output) == Matrix{VectorState}

        iqs=createInitialQubitState(density, [[1.00,[1. 0.;1. 0.]],[0.00,[0. 1.;0. 1.]]])
        qo=runQuantumProgram(qp, iqs, 2)
        @test typeof(qo) == QuantumOutput{false, circuitIndexBigEndian, DensityState}
        @test typeof(qp.settings) == Settings{false, circuitIndexBigEndian}
        @test typeof(qo.output) == Matrix{DensityState}


        qc=createQuantumCircuit(2,indexType=circuitIndexLittleEndian)
        hGate!(qc, 1)
        cnotGate!(qc, 1, 2)
        measureGate!(qc, [1,2])
        @test typeof(qc.settings) == Settings{false, circuitIndexLittleEndian}
        @test typeof(qc) == QuantumCircuit{false, circuitIndexLittleEndian}

        qp=compileQuantumCircuit(qc)
        @test typeof(qp) == QuantumProgram{false, circuitIndexLittleEndian}
        @test typeof(qp.settings) == Settings{false, circuitIndexLittleEndian}

        iqs=createInitialQubitState(vector, [1. 0.;1. 0.])
        qo=runQuantumProgram(qp, iqs, 2)
        @test typeof(qo) == QuantumOutput{false, circuitIndexLittleEndian, VectorState}
        @test typeof(qo.settings) == Settings{false, circuitIndexLittleEndian}
        @test typeof(qo.output) == Matrix{VectorState}

        iqs=createInitialQubitState(density, [[1.00,[1. 0.;1. 0.]],[0.00,[0. 1.;0. 1.]]])
        qo=runQuantumProgram(qp, iqs, 2)
        @test typeof(qo) == QuantumOutput{false, circuitIndexLittleEndian, DensityState}
        @test typeof(qp.settings) == Settings{false, circuitIndexLittleEndian}
        @test typeof(qo.output) == Matrix{DensityState}


        qc=createQuantumCircuit(2,zeroBasedNumbering=true)
        hGate!(qc, 0)
        cnotGate!(qc, 0, 1)
        measureGate!(qc, [0,1])
        @test typeof(qc.settings) == Settings{true, circuitIndexBigEndian}
        @test typeof(qc) == QuantumCircuit{true, circuitIndexBigEndian}

        qp=compileQuantumCircuit(qc)
        @test typeof(qp) == QuantumProgram{true, circuitIndexBigEndian}
        @test typeof(qp.settings) == Settings{true, circuitIndexBigEndian}

        iqs=createInitialQubitState(vector, [1. 0.;1. 0.])
        qo=runQuantumProgram(qp, iqs, 2)
        @test typeof(qo) == QuantumOutput{true, circuitIndexBigEndian, VectorState}
        @test typeof(qo.settings) == Settings{true, circuitIndexBigEndian}
        @test typeof(qo.output) == Matrix{VectorState}

        iqs=createInitialQubitState(density, [[1.00,[1. 0.;1. 0.]],[0.00,[0. 1.;0. 1.]]])
        qo=runQuantumProgram(qp, iqs, 2)
        @test typeof(qo) == QuantumOutput{true, circuitIndexBigEndian, DensityState}
        @test typeof(qp.settings) == Settings{true, circuitIndexBigEndian}
        @test typeof(qo.output) == Matrix{DensityState}
    end

    @testset "createInitialQubitState" verbose = true begin
        @testset "Output: Pure States (vector)" begin
            # **Example 1: Big-endian, ``|ψ⟩ = |0⟩⊗|1⟩ → [0; 1; 0; 0]``**
            iqs = createInitialQubitState(vector, [1. 0.; 0. 1.])
            @test iqs isa VectorState
            @test iqs.measured == nothing
            @test iqs.q ≈ ComplexF64[0.0 + 0.0im; 1.0 + 0.0im; 0.0 + 0.0im; 0.0 + 0.0im;;] atol=1e-10
            # **Example 2: Little-endian, same state**
            iqs = createInitialQubitState(vector, [0. 1.; 1. 0.], indexType=circuitIndexLittleEndian)
            @test iqs isa VectorState
            @test iqs.measured == nothing
            @test iqs.q ≈ ComplexF64[0.0 + 0.0im; 1.0 + 0.0im; 0.0 + 0.0im; 0.0 + 0.0im;;] atol=1e-10
            # **Example 3: Byte index, ``|ψ⟩ = |0⟩⊗|1⟩⊗|1⟩ → [0; 0; 0; 1; 0; 0; 0; 0]``**
            iqs = createInitialQubitState(vector, [1. 0. 0.; 0. 1. 1.], indexType=byteIndex)
            @test iqs isa VectorState
            @test iqs.measured == nothing
            @test iqs.q ≈ ComplexF64[0.0 + 0.0im; 0.0 + 0.0im; 0.0 + 0.0im; 1.0 + 0.0im; 0.0 + 0.0im; 0.0 + 0.0im; 0.0 + 0.0im; 0.0 + 0.0im;;] atol=1e-10
            # **Example 4: Byte index for 3 qubits, ``|ψ⟩ = |6⟩ = |1⟩⊗|1⟩⊗|0⟩ → [0; 0; 0; 0; 0; 0; 1; 0]``**
            iqs = createInitialQubitState(vector, createByteIndexVector(6, 3), indexType=byteIndex)
            @test iqs isa VectorState
            @test iqs.measured == nothing
            @test iqs.q ≈ ComplexF64[0.0 + 0.0im; 0.0 + 0.0im; 0.0 + 0.0im; 0.0 + 0.0im; 0.0 + 0.0im; 0.0 + 0.0im; 1.0 + 0.0im; 0.0 + 0.0im;;] atol=1e-10
            # **Example 5: Big-endian, Bloch representation: ``|ψ⟩ = |+⟩⊗|–⟩⊗|1⟩``**
            iqs = createInitialQubitState(vector, [π/2 0.; π/2 π/2], blochRepresentation=true)
            @test iqs isa VectorState
            @test iqs.measured == nothing
            @test iqs.q ≈ ComplexF64[0.5000000000000001 + 0.0im; 3.061616997868383e-17 + 0.5im; 0.5 + 0.0im; 3.0616169978683824e-17 + 0.4999999999999999im;;] atol=1e-10

            iqs = createInitialQubitState(vector, [1. 1.; 0. 1.; 1. 0.]/sqrt(2))
            @test iqs isa VectorState
            @test iqs.measured == nothing
            @test iqs.q ≈ ComplexF64[0.0 + 0.0im;
                0.0 + 0.0im;
 0.3535533905932737 + 0.0im;
                0.0 + 0.0im;
                0.0 + 0.0im;
                0.0 + 0.0im;
 0.3535533905932737 + 0.0im;
                0.0 + 0.0im] atol=1e-10

            iqs = createInitialQubitState(vector,[1. 1.; 0. 1.; 1. 0.]/sqrt(2), indexType=circuitIndexLittleEndian)
            @test iqs isa VectorState
            @test iqs.measured == nothing
            @test iqs.q ≈ ComplexF64[0.0 + 0.0im;
                0.0 + 0.0im;
 0.3535533905932737 + 0.0im;
 0.3535533905932737 + 0.0im;
                0.0 + 0.0im;
                0.0 + 0.0im;
                0.0 + 0.0im;
                0.0 + 0.0im] atol=1e-10

            iqs = createInitialQubitState(vector, [1. 1. 0.; -1. 1. 1.]/sqrt(2), indexType=byteIndex)
            @test iqs isa VectorState
            @test iqs.measured == nothing
            @test iqs.q ≈ ComplexF64[0.0 + 0.0im;
  0.3535533905932737 + 0.0im;
                 0.0 + 0.0im;
  0.3535533905932737 + 0.0im;
                -0.0 - 0.0im;
 -0.3535533905932737 - 0.0im;
                -0.0 - 0.0im;
 -0.3535533905932737 - 0.0im] atol=1e-10

                       iqs = createInitialQubitState(vector, [1. 1.; 0. 1.; 1. 0.]/sqrt(2), blochRepresentation=true)
            @test iqs isa VectorState
            @test iqs.measured == nothing
            @test iqs.q ≈ ComplexF64[0.8801222985378151 + 0.0im;
  0.3248184695400312 + 0.0im;
                 0.0 + 0.0im;
                 0.0 + 0.0im;
 0.24694148649818384 + 0.21101407630865632im;
 0.09113637484647143 + 0.07787698304184731im;
                 0.0 + 0.0im;
                 0.0 + 0.0im] atol=1e-10

            iqs = createInitialQubitState(vector,[1. 1.; 0. 1.; 1. 0.]/sqrt(2), indexType=circuitIndexLittleEndian, blochRepresentation=true)
            @test iqs isa VectorState
            @test iqs.measured == nothing
            @test iqs.q ≈ ComplexF64[0.8801222985378151 + 0.0im;
 0.24694148649818384 + 0.21101407630865632im;
                 0.0 + 0.0im;
                 0.0 + 0.0im;
  0.3248184695400312 + 0.0im;
 0.09113637484647143 + 0.07787698304184731im;
                 0.0 + 0.0im;
                 0.0 + 0.0im] atol=1e-10

            iqs = createInitialQubitState(vector, [1. 1. 0.; -1. 1. 1.]/sqrt(2), indexType=byteIndex, blochRepresentation=true)
            @test iqs isa VectorState
            @test iqs.measured == nothing
            @test iqs.q ≈ ComplexF64[0.8801222985378151 + 0.0im;
                 0.0 + 0.0im;
 0.24694148649818384 + 0.21101407630865635im;
                 0.0 + 0.0im;
 0.24694148649818384 - 0.21101407630865632im;
                 0.0 + 0.0im;
 0.11987770146218488 + 0.0im;
                 0.0 + 0.0im] atol=1e-10
        end

        @testset "Output: Mixed States (density)" begin
            # **Example 7: Big-endian, ``ρ = 0.25|+⟩⟨+| + 0.75|–⟩⟨–| → [0.5, -0.25; -0.25, 0.5]``**
            iqs = createInitialQubitState(density, [[0.25, [1. 1.]/sqrt(2)], [0.75, [1. -1.]/sqrt(2)]])
            @test iqs isa DensityState
            @test iqs.measured == nothing
            @test iqs.rho ≈ ComplexF64[0.4999999999999999 + 0.0im -0.24999999999999992 + 0.0im; -0.24999999999999992 + 0.0im 0.4999999999999999 + 0.0im] atol=1e-10
            # **Example 8: Same using Bloch representation**
            iqs = createInitialQubitState(density, [[0.25, [π/2 0.]], [0.75, [π/2 π]]], blochRepresentation=true)
            @test iqs isa DensityState
            @test iqs.measured == nothing
            @test iqs.rho ≈ ComplexF64[0.5000000000000001 + 0.0im -0.25 - 4.592425496802574e-17im; -0.25 + 4.592425496802574e-17im 0.4999999999999999 + 0.0im] atol=1e-10

            iqs = createInitialQubitState(density, [[0.25, [1. 1.; 0. 1.; 1. 0.]/sqrt(2)], [0.75, [1. -1.; 1. 0.; 1. -1.]/sqrt(2)]])
            @test iqs isa DensityState
            @test iqs.measured == nothing
            @test iqs.rho ≈ ComplexF64[  0.09375+0.0im  -0.09375+0.0im      0.0+0.0im  0.0+0.0im  -0.09375+0.0im   0.09375+0.0im      0.0+0.0im  0.0+0.0im;
 -0.09375+0.0im   0.09375+0.0im      0.0+0.0im  0.0+0.0im   0.09375+0.0im  -0.09375+0.0im      0.0+0.0im  0.0+0.0im;
      0.0+0.0im       0.0+0.0im  0.03125+0.0im  0.0+0.0im       0.0+0.0im       0.0+0.0im  0.03125+0.0im  0.0+0.0im;
      0.0+0.0im       0.0+0.0im      0.0+0.0im  0.0+0.0im       0.0+0.0im       0.0+0.0im      0.0+0.0im  0.0+0.0im;
 -0.09375+0.0im   0.09375+0.0im      0.0+0.0im  0.0+0.0im   0.09375+0.0im  -0.09375+0.0im      0.0+0.0im  0.0+0.0im;
  0.09375+0.0im  -0.09375+0.0im      0.0+0.0im  0.0+0.0im  -0.09375+0.0im   0.09375+0.0im      0.0+0.0im  0.0+0.0im;
      0.0+0.0im       0.0+0.0im  0.03125+0.0im  0.0+0.0im       0.0+0.0im       0.0+0.0im  0.03125+0.0im  0.0+0.0im;
      0.0+0.0im       0.0+0.0im      0.0+0.0im  0.0+0.0im       0.0+0.0im       0.0+0.0im      0.0+0.0im  0.0+0.0im] atol=1e-10

            iqs = createInitialQubitState(density, [[0.25, [1. 1.; 0. 1.; 1. 0.]/sqrt(2)], [0.75, [1. -1.; 1. 0.; 1. -1.]/sqrt(2)]], indexType=circuitIndexLittleEndian)
            @test iqs isa DensityState
            @test iqs.measured == nothing
            @test iqs.rho ≈ ComplexF64[  0.09375+0.0im  -0.09375+0.0im      0.0+0.0im      0.0+0.0im  -0.09375+0.0im   0.09375+0.0im  0.0+0.0im  0.0+0.0im;
 -0.09375+0.0im   0.09375+0.0im      0.0+0.0im      0.0+0.0im   0.09375+0.0im  -0.09375+0.0im  0.0+0.0im  0.0+0.0im;
      0.0+0.0im       0.0+0.0im  0.03125+0.0im  0.03125+0.0im       0.0+0.0im       0.0+0.0im  0.0+0.0im  0.0+0.0im;
      0.0+0.0im       0.0+0.0im  0.03125+0.0im  0.03125+0.0im       0.0+0.0im       0.0+0.0im  0.0+0.0im  0.0+0.0im;
 -0.09375+0.0im   0.09375+0.0im      0.0+0.0im      0.0+0.0im   0.09375+0.0im  -0.09375+0.0im  0.0+0.0im  0.0+0.0im;
  0.09375+0.0im  -0.09375+0.0im      0.0+0.0im      0.0+0.0im  -0.09375+0.0im   0.09375+0.0im  0.0+0.0im  0.0+0.0im;
      0.0+0.0im       0.0+0.0im      0.0+0.0im      0.0+0.0im       0.0+0.0im       0.0+0.0im  0.0+0.0im  0.0+0.0im;
      0.0+0.0im       0.0+0.0im      0.0+0.0im      0.0+0.0im       0.0+0.0im       0.0+0.0im  0.0+0.0im  0.0+0.0im] atol=1e-10

            iqs = createInitialQubitState(density, [[0.25, [1. 1. 0.; -1. 1. 1.]/sqrt(2)], [0.75, [1. -1. 1.; -1. 1. 0.]/sqrt(2)]], indexType=byteIndex)
            @test iqs isa DensityState
            @test iqs.measured == nothing
            @test iqs.rho ≈ ComplexF64[  0.09375+0.0im       0.0+0.0im  -0.09375+0.0im       0.0+0.0im  -0.09375+0.0im       0.0+0.0im   0.09375+0.0im       0.0+0.0im;
      0.0+0.0im   0.03125+0.0im       0.0+0.0im   0.03125+0.0im       0.0+0.0im  -0.03125+0.0im       0.0+0.0im  -0.03125+0.0im;
 -0.09375+0.0im       0.0+0.0im   0.09375+0.0im       0.0+0.0im   0.09375+0.0im       0.0+0.0im  -0.09375+0.0im       0.0+0.0im;
      0.0+0.0im   0.03125+0.0im       0.0+0.0im   0.03125+0.0im       0.0+0.0im  -0.03125+0.0im       0.0+0.0im  -0.03125+0.0im;
 -0.09375+0.0im       0.0+0.0im   0.09375+0.0im       0.0+0.0im   0.09375+0.0im       0.0+0.0im  -0.09375+0.0im       0.0+0.0im;
      0.0+0.0im  -0.03125+0.0im       0.0+0.0im  -0.03125+0.0im       0.0+0.0im   0.03125+0.0im       0.0+0.0im   0.03125+0.0im;
  0.09375+0.0im       0.0+0.0im  -0.09375+0.0im       0.0+0.0im  -0.09375+0.0im       0.0+0.0im   0.09375+0.0im       0.0+0.0im;
      0.0+0.0im  -0.03125+0.0im       0.0+0.0im  -0.03125+0.0im       0.0+0.0im   0.03125+0.0im       0.0+0.0im   0.03125+0.0im] atol=1e-10

            iqs = createInitialQubitState(density, [[0.25, [1. 1.; 0. 1.]/sqrt(2)], [0.75, [1. -1.; 1. 0.]/sqrt(2)]], blochRepresentation=true)
            @test iqs isa DensityState
            @test iqs.measured == nothing
            @test iqs.rho ≈ ComplexF64[0.800992+0.0im          0.21441+0.0im         0.224739+0.0865351im    0.0601584+0.051406im;
   0.21441+0.0im        0.0791303+0.0im        0.0601584+0.051406im     0.0222021+0.0189719im;
  0.224739-0.0865351im  0.0601584-0.051406im      0.1091+0.0im          0.0292039+3.37144e-18im;
 0.0601584-0.051406im   0.0222021-0.0189719im  0.0292039-3.37144e-18im   0.010778+0.0im] atol=1e-6

            iqs = createInitialQubitState(density, [[0.25, [1. 1.; 0. 1.]/sqrt(2)], [0.75, [1. -1.; 1. 0.]/sqrt(2)]], indexType=circuitIndexLittleEndian, blochRepresentation=true)
            @test iqs isa DensityState
            @test iqs.measured == nothing
            @test iqs.rho ≈ ComplexF64[0.800992+0.0im         0.224739+0.0865351im      0.21441+0.0im        0.0601584+0.051406im;
  0.224739-0.0865351im     0.1091+0.0im          0.0601584-0.051406im   0.0292039+3.37144e-18im;
   0.21441+0.0im        0.0601584+0.051406im     0.0791303+0.0im        0.0222021+0.0189719im;
 0.0601584-0.051406im   0.0292039-3.37144e-18im  0.0222021-0.0189719im   0.010778+0.0im] atol=1e-6

            iqs = createInitialQubitState(density, [[0.25, [1. 1.; -1. 1.]/sqrt(2)], [0.75, [1. -1.; -1. 1.]/sqrt(2)]], indexType=byteIndex, blochRepresentation=true)
            @test iqs isa DensityState
            @test iqs.measured == nothing
            @test iqs.rho ≈ ComplexF64[0.774615+0.0im          -0.108669+0.0928591im     0.217339+0.185718im   -0.0527535+0.0im;
  -0.108669-0.0928591im     0.105507+0.0im        -0.00822658-0.0521081im   0.0296028+0.0252959im;
   0.217339-0.185718im   -0.00822658+0.0521081im     0.105507+0.0im        -0.0148014+0.0126479im;
 -0.0527535+0.0im          0.0296028-0.0252959im   -0.0148014-0.0126479im   0.0143707+0.0im] atol=1e-6

#             # **Example 7: Big-endian, ``ρ = 0.25|+⟩⟨+| + 0.75|–⟩⟨–| → [0.5, -0.25; -0.25, 0.5]``**
#             iqs = createInitialQubitState(enrichedDensity, [[0.25, [1. 1.]/sqrt(2)], [0.75, [1. -1.]/sqrt(2)]])
#             @test iqs isa EnrichedDensityState
#             @test iqs.measured == nothing
#             @test iqs.rho ≈ ComplexF64[0.4999999999999999 + 0.0im -0.24999999999999992 + 0.0im; -0.24999999999999992 + 0.0im 0.4999999999999999 + 0.0im] atol=1e-10
#             @test iqs.entropy ≈ 0.5623351446188085 atol=1e-10
#             @test iqs.purityValue ≈ 0.6249999999999997 atol=1e-10
#             @test iqs.purityNorm ≈ 0.24999999999999933 atol=1e-10
#             # **Example 8: Same using Bloch representation**
#             iqs = createInitialQubitState(enrichedDensity, [[0.25, [π/2 0.]], [0.75, [π/2 π]]], blochRepresentation=true)
#             @test iqs isa EnrichedDensityState
#             @test iqs.measured == nothing
#             @test iqs.rho ≈ ComplexF64[0.5000000000000001 + 0.0im -0.25 - 4.592425496802574e-17im; -0.25 + 4.592425496802574e-17im 0.4999999999999999 + 0.0im] atol=1e-10
#             @test iqs.entropy ≈ 0.5623351446188083 atol=1e-10
#             @test iqs.purityValue ≈ 0.625 atol=1e-10
#             @test iqs.purityNorm ≈ 0.25 atol=1e-10

#             iqs = createInitialQubitState(enrichedDensity, [[0.25, [1. 1.; 0. 1.; 1. 0.]/sqrt(2)], [0.75, [1. -1.; 1. 0.; 1. -1.]/sqrt(2)]])
#             @test iqs isa EnrichedDensityState
#             @test iqs.measured == nothing
#             @test iqs.rho ≈ ComplexF64[  0.09375+0.0im  -0.09375+0.0im      0.0+0.0im  0.0+0.0im  -0.09375+0.0im   0.09375+0.0im      0.0+0.0im  0.0+0.0im;
#  -0.09375+0.0im   0.09375+0.0im      0.0+0.0im  0.0+0.0im   0.09375+0.0im  -0.09375+0.0im      0.0+0.0im  0.0+0.0im;
#       0.0+0.0im       0.0+0.0im  0.03125+0.0im  0.0+0.0im       0.0+0.0im       0.0+0.0im  0.03125+0.0im  0.0+0.0im;
#       0.0+0.0im       0.0+0.0im      0.0+0.0im  0.0+0.0im       0.0+0.0im       0.0+0.0im      0.0+0.0im  0.0+0.0im;
#  -0.09375+0.0im   0.09375+0.0im      0.0+0.0im  0.0+0.0im   0.09375+0.0im  -0.09375+0.0im      0.0+0.0im  0.0+0.0im;
#   0.09375+0.0im  -0.09375+0.0im      0.0+0.0im  0.0+0.0im  -0.09375+0.0im   0.09375+0.0im      0.0+0.0im  0.0+0.0im;
#       0.0+0.0im       0.0+0.0im  0.03125+0.0im  0.0+0.0im       0.0+0.0im       0.0+0.0im  0.03125+0.0im  0.0+0.0im;
#       0.0+0.0im       0.0+0.0im      0.0+0.0im  0.0+0.0im       0.0+0.0im       0.0+0.0im      0.0+0.0im  0.0+0.0im] atol=1e-10
#             @test iqs.entropy ≈ 0.5410977650193832 atol=1e-10
#             @test iqs.purityValue ≈ 0.14453124999999983 atol=1e-10
#             @test iqs.purityNorm ≈ 0.02232142857142838 atol=1e-10

#             iqs = createInitialQubitState(enrichedDensity, [[0.25, [1. 1.; 0. 1.; 1. 0.]/sqrt(2)], [0.75, [1. -1.; 1. 0.; 1. -1.]/sqrt(2)]], indexType=circuitIndexLittleEndian)
#             @test iqs isa EnrichedDensityState
#             @test iqs.measured == nothing
#             @test iqs.rho ≈ ComplexF64[  0.09375+0.0im  -0.09375+0.0im      0.0+0.0im      0.0+0.0im  -0.09375+0.0im   0.09375+0.0im  0.0+0.0im  0.0+0.0im;
#  -0.09375+0.0im   0.09375+0.0im      0.0+0.0im      0.0+0.0im   0.09375+0.0im  -0.09375+0.0im  0.0+0.0im  0.0+0.0im;
#       0.0+0.0im       0.0+0.0im  0.03125+0.0im  0.03125+0.0im       0.0+0.0im       0.0+0.0im  0.0+0.0im  0.0+0.0im;
#       0.0+0.0im       0.0+0.0im  0.03125+0.0im  0.03125+0.0im       0.0+0.0im       0.0+0.0im  0.0+0.0im  0.0+0.0im;
#  -0.09375+0.0im   0.09375+0.0im      0.0+0.0im      0.0+0.0im   0.09375+0.0im  -0.09375+0.0im  0.0+0.0im  0.0+0.0im;
#   0.09375+0.0im  -0.09375+0.0im      0.0+0.0im      0.0+0.0im  -0.09375+0.0im   0.09375+0.0im  0.0+0.0im  0.0+0.0im;
#       0.0+0.0im       0.0+0.0im      0.0+0.0im      0.0+0.0im       0.0+0.0im       0.0+0.0im  0.0+0.0im  0.0+0.0im;
#       0.0+0.0im       0.0+0.0im      0.0+0.0im      0.0+0.0im       0.0+0.0im       0.0+0.0im  0.0+0.0im  0.0+0.0im] atol=1e-10
#             @test iqs.entropy ≈ 0.5410977650193813 atol=1e-10
#             @test iqs.purityValue ≈ 0.1445312499999998 atol=1e-10
#             @test iqs.purityNorm ≈ 0.02232142857142835 atol=1e-10

#             iqs = createInitialQubitState(enrichedDensity, [[0.25, [1. 1. 0.; -1. 1. 1.]/sqrt(2)], [0.75, [1. -1. 1.; -1. 1. 0.]/sqrt(2)]], indexType=byteIndex)
#             @test iqs isa EnrichedDensityState
#             @test iqs.measured == nothing
#             @test iqs.rho ≈ ComplexF64[  0.09375+0.0im       0.0+0.0im  -0.09375+0.0im       0.0+0.0im  -0.09375+0.0im       0.0+0.0im   0.09375+0.0im       0.0+0.0im;
#       0.0+0.0im   0.03125+0.0im       0.0+0.0im   0.03125+0.0im       0.0+0.0im  -0.03125+0.0im       0.0+0.0im  -0.03125+0.0im;
#  -0.09375+0.0im       0.0+0.0im   0.09375+0.0im       0.0+0.0im   0.09375+0.0im       0.0+0.0im  -0.09375+0.0im       0.0+0.0im;
#       0.0+0.0im   0.03125+0.0im       0.0+0.0im   0.03125+0.0im       0.0+0.0im  -0.03125+0.0im       0.0+0.0im  -0.03125+0.0im;
#  -0.09375+0.0im       0.0+0.0im   0.09375+0.0im       0.0+0.0im   0.09375+0.0im       0.0+0.0im  -0.09375+0.0im       0.0+0.0im;
#       0.0+0.0im  -0.03125+0.0im       0.0+0.0im  -0.03125+0.0im       0.0+0.0im   0.03125+0.0im       0.0+0.0im   0.03125+0.0im;
#   0.09375+0.0im       0.0+0.0im  -0.09375+0.0im       0.0+0.0im  -0.09375+0.0im       0.0+0.0im   0.09375+0.0im       0.0+0.0im;
#       0.0+0.0im  -0.03125+0.0im       0.0+0.0im  -0.03125+0.0im       0.0+0.0im   0.03125+0.0im       0.0+0.0im   0.03125+0.0im] atol=1e-10
#             @test iqs.entropy ≈ 0.6277411625893783 atol=1e-10
#             @test iqs.purityValue ≈ 0.1562499999999999 atol=1e-10
#             @test iqs.purityNorm ≈ 0.03571428571428559 atol=1e-10

#             iqs = createInitialQubitState(enrichedDensity, [[0.25, [1. 1.; 0. 1.]/sqrt(2)], [0.75, [1. -1.; 1. 0.]/sqrt(2)]], blochRepresentation=true)
#             @test iqs isa EnrichedDensityState
#             @test iqs.measured == nothing
#             @test iqs.rho ≈ ComplexF64[0.800992+0.0im          0.21441+0.0im         0.224739+0.0865351im    0.0601584+0.051406im;
#    0.21441+0.0im        0.0791303+0.0im        0.0601584+0.051406im     0.0222021+0.0189719im;
#   0.224739-0.0865351im  0.0601584-0.051406im      0.1091+0.0im          0.0292039+3.37144e-18im;
#  0.0601584-0.051406im   0.0222021-0.0189719im  0.0292039-3.37144e-18im   0.010778+0.0im] atol=1e-6
#             @test iqs.entropy ≈ 0.2126431887606028 atol=1e-10
#             @test iqs.purityValue ≈ 0.8962621322536026 atol=1e-10
#             @test iqs.purityNorm ≈ 0.8616828430048035 atol=1e-10

#             iqs = createInitialQubitState(enrichedDensity, [[0.25, [1. 1.; 0. 1.]/sqrt(2)], [0.75, [1. -1.; 1. 0.]/sqrt(2)]], indexType=circuitIndexLittleEndian, blochRepresentation=true)
#             @test iqs isa EnrichedDensityState
#             @test iqs.measured == nothing
#             @test iqs.rho ≈ ComplexF64[0.800992+0.0im         0.224739+0.0865351im      0.21441+0.0im        0.0601584+0.051406im;
#   0.224739-0.0865351im     0.1091+0.0im          0.0601584-0.051406im   0.0292039+3.37144e-18im;
#    0.21441+0.0im        0.0601584+0.051406im     0.0791303+0.0im        0.0222021+0.0189719im;
#  0.0601584-0.051406im   0.0292039-3.37144e-18im  0.0222021-0.0189719im   0.010778+0.0im] atol=1e-6
#             @test iqs.entropy ≈ 0.2126431887606013 atol=1e-10
#             @test iqs.purityValue ≈ 0.8962621322536012 atol=1e-10
#             @test iqs.purityNorm ≈ 0.8616828430048016 atol=1e-10

#             iqs = createInitialQubitState(enrichedDensity, [[0.25, [1. 1.; -1. 1.]/sqrt(2)], [0.75, [1. -1.; -1. 1.]/sqrt(2)]], indexType=byteIndex, blochRepresentation=true)
#             @test iqs isa EnrichedDensityState
#             @test iqs.measured == nothing
#             @test iqs.rho ≈ ComplexF64[0.774615+0.0im          -0.108669+0.0928591im     0.217339+0.185718im   -0.0527535+0.0im;
#   -0.108669-0.0928591im     0.105507+0.0im        -0.00822658-0.0521081im   0.0296028+0.0252959im;
#    0.217339-0.185718im   -0.00822658+0.0521081im     0.105507+0.0im        -0.0148014+0.0126479im;
#  -0.0527535+0.0im          0.0296028-0.0252959im   -0.0148014-0.0126479im   0.0143707+0.0im] atol=1e-6
#             @test iqs.entropy ≈ 0.29468509128381004 atol=1e-10
#             @test iqs.purityValue ≈ 0.841739442768507 atol=1e-10
#             @test iqs.purityNorm ≈ 0.7889859236913427 atol=1e-10
        end
    end

    @testset "H gate" verbose = true begin
        @testset "Program: Big endian 1 qubit H-gate on qubit 1" begin
            qc=createQuantumCircuit(1)
            hGate!(qc, 1)
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

            @test qp.program[1].U ≈ [sqrt(2.)/2+0im sqrt(2.)/2+0im;sqrt(2.)/2+0im -sqrt(2.)/2+0im] atol=1e-10
        end

        @testset "Output: Little endian 2 qubits H-gate on qubit 1" begin
            # 1. Test 2 qubits, hGate on circuitIndex 1, littleEndian, initial state |01>, assert output=|0->
            qc=createQuantumCircuit(2, indexType=circuitIndexLittleEndian)
            hGate!(qc, 1)
            measureGate!(qc, [1,2])
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            iqs=createInitialQubitState(vector, [1. 0.;0. 1.], indexType=byteIndex)
            qo=runQuantumProgram(qp, iqs, 1)

            @test qo.output[2].q ≈ [0.7071;-0.7071;0;0;;] atol=0.0001
        end

        @testset "Output: Big endian 2 qubits H-gate on qubit 1" begin
            # 2. Test 2 qubits, hGate on circuitIndex 1, bigEndian, initial state |10>, assert output=|-0>
            qc=createQuantumCircuit(2, indexType=circuitIndexBigEndian)
            hGate!(qc, 1)
            measureGate!(qc, [1,2])
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            iqs=createInitialQubitState(vector, [0. 1.;1. 0.], indexType=byteIndex)
            qo=runQuantumProgram(qp, iqs, 1)

            @test qo.output[2].q ≈ [0.7071;0;-0.7071;0;;] atol=0.0001
        end        
        
        @testset "Output: Little endian 2 qubits H-gate on qubit 2" begin
            # 3. Test 2 qubits, hGate on circuitIndex 2, littleEndian, initial state |01>, assert output=|+1>
            qc=createQuantumCircuit(2, indexType=circuitIndexLittleEndian)
            hGate!(qc, 2)
            measureGate!(qc, [1,2])
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            iqs=createInitialQubitState(vector, [1. 0.;0. 1.], indexType=byteIndex)
            qo=runQuantumProgram(qp, iqs, 1)

            @test qo.output[2].q ≈ [0;0.7071;0;0.7071;;] atol=0.0001
        end        
        
        @testset "Output: Big endian 2 qubits H-gate on qubit 2" begin
            # 4. Test 2 qubits, hGate on circuitIndex 2, bigEndian, initial state |10>, assert output=|1+>
            qc=createQuantumCircuit(2, indexType=circuitIndexBigEndian)
            hGate!(qc, 2)
            measureGate!(qc, [1,2])
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            iqs=createInitialQubitState(vector, [0. 1.;1. 0.], indexType=byteIndex)
            qo=runQuantumProgram(qp, iqs, 1)

            @test qo.output[2].q ≈ [0;0;0.7071;0.7071;;] atol=0.0001
        end
    end

    @testset "X gate" verbose = true begin
        @testset "Program: Big endian 1 qubit X-gate on qubit 1" begin
            qc=createQuantumCircuit(1)
            xGate!(qc, 1)
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

            @test qp.program[1].U ≈ [0.0+0.0im 1.0+0.0im;1.0+0.0im 0.0+0.0im] atol=1e-10
        end
    end

    @testset "Y gate" verbose = true begin
        @testset "Program: Big endian 1 qubit Y-gate on qubit 1" begin
            qc=createQuantumCircuit(1)
            yGate!(qc, 1)
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

            @test qp.program[1].U ≈ [0.0+0.0im 0.0-1.0im;0.0+1.0im 0.0+0.0im] atol=1e-10
        end
    end

    @testset "Z gate" verbose = true begin
        @testset "Program: Big endian 1 qubit Z-gate on qubit 1" begin
            qc=createQuantumCircuit(1)
            zGate!(qc, 1)
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

            @test qp.program[1].U ≈ [1.0+0.0im 0.0+0.0im;0.0+0.0im -1.0+0.0im] atol=1e-10
        end
    end

    @testset "Id gate" verbose = true begin
        @testset "Program: Big endian 1 qubit Id-gate on qubit 1" begin
            qc=createQuantumCircuit(1)
            idGate!(qc, 1)
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

            @test qp.program[1].U ≈ [1.0+0.0im 0.0+0.0im;0.0+0.0im 1.0+0.0im] atol=1e-10
        end
    end

    @testset "T gate" verbose = true begin
        @testset "Program: Big endian 1 qubit T-gate on qubit 1" begin
            qc=createQuantumCircuit(1)
            tGate!(qc, 1)
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

            @test qp.program[1].U ≈ [1.0+0.0im 0.0+0.0im;0.0+0.0im (1.0+1.0im)*sqrt(2.)/2] atol=1e-10
        end
    end

    @testset "Td gate" verbose = true begin
        @testset "Program: Big endian 1 qubit Td-gate on qubit 1" begin
            qc=createQuantumCircuit(1)
            tdGate!(qc, 1)
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

            @test qp.program[1].U ≈ [1.0+0.0im 0.0+0.0im;0.0+0.0im (1.0-1.0im)*sqrt(2.)/2] atol=1e-10
        end
    end

    @testset "S gate" verbose = true begin
        @testset "Program: Big endian 1 qubit S-gate on qubit 1" begin
            qc=createQuantumCircuit(1)
            sGate!(qc, 1)
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

            @test qp.program[1].U ≈ [1.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+1.0im] atol=1e-10
        end
    end

    @testset "Sd gate" verbose = true begin
        @testset "Program: Big endian 1 qubit Sd-gate on qubit 1" begin
            qc=createQuantumCircuit(1)
            sdGate!(qc, 1)
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

            @test qp.program[1].U ≈ [1.0+0.0im 0.0+0.0im;0.0+0.0im 0.0-1.0im] atol=1e-10
        end
    end

    @testset "U1 gate" verbose = true begin
        @testset "Program: Big endian 1 qubit U1-gate on qubit 1" begin
            qc=createQuantumCircuit(1)
            u1Gate!(qc, 1, 1.2)
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

            @test qp.program[1].U ≈ [1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.3623577544766736 + 0.9320390859672263im] atol=1e-10
        end
    end

    @testset "U2 gate" verbose = true begin
        @testset "Program: Big endian 1 qubit U2-gate on qubit 1" begin
            qc=createQuantumCircuit(1)
            u2Gate!(qc, 1, 1.2, 1.3)
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

            @test qp.program[1].U ≈ [0.7071067811865476 + 0.0im -0.18915023567990386 - 0.6813385269763018im; 0.2562256254059859 + 0.6590511580183371im -0.5664940832575452 + 0.42318371144716027im] atol=1e-10
        end
    end

    @testset "U3 gate" verbose = true begin
        @testset "Program: Big endian 1 qubit U3-gate on qubit 1" begin
            qc=createQuantumCircuit(1)
            u3Gate!(qc, 1, 1.2, 1.3, 1.4)
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

            @test qp.program[1].U ≈ [0.8253356149096783 + 0.0im -0.09597066796307953 - 0.5564267729471539im; 0.1510412002248617 + 0.5440658770739959im -0.7461629372543612 + 0.35273183625281257im] atol=1e-10
        end
    end

    @testset "CNOT gate" verbose = true begin
        @testset "Program: Big endian 2 qubits CNOT-gate on qubits c=1 and t=2" begin
            qc=createQuantumCircuit(2)
            cnotGate!(qc, 1, 2)
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

            @test qp.program[1].U ≈ [1.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 1.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 1.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 1.0 + 0.0im 0.0 + 0.0im] atol=1e-10
        end

        @testset "Output: Little endian 3 qubits CNOT-gate on qubits c=1 and t=2" begin
            # 5. Test 3 qubits, cnotGate on circuitIndices c=1 and t=2, littleEndian, initial state |001>, assert output=|011>
            qc=createQuantumCircuit(3, indexType=circuitIndexLittleEndian)
            cnotGate!(qc, 1, 2)
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            iqs=createInitialQubitState(vector, [1. 1. 0.;0. 0. 1.], indexType=byteIndex)
            qo=runQuantumProgram(qp, iqs, 1)

            @test qo.output[2].q == [0;0;0;1;0;0;0;0;;]
        end
        
        @testset "Output: Little endian 3 qubits CNOT-gate on qubits c=2 and t=3" begin
            # 6. Test 3 qubits, cnotGate on circuitIndices c=2 and t=3, littleEndian, initial state |010>, assert output=|110>
            qc=createQuantumCircuit(3, indexType=circuitIndexLittleEndian)
            cnotGate!(qc, 2, 3)
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            iqs=createInitialQubitState(vector, [1. 0. 1.;0. 1. 0.], indexType=byteIndex)
            qo=runQuantumProgram(qp, iqs, 1)

            @test qo.output[2].q == [0;0;0;0;0;0;1;0;;]
        end
                
        @testset "Output: Little endian 3 qubits CNOT-gate on qubits c=1 and t=3" begin
            # 7. Test 3 qubits, cnotGate on circuitIndices c=1 and t=3, littleEndian, initial state |001>, assert output=|101>
            qc=createQuantumCircuit(3, indexType=circuitIndexLittleEndian)
            cnotGate!(qc, 1, 3)
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            iqs=createInitialQubitState(vector, [1. 1. 0.;0. 0. 1.], indexType=byteIndex)
            qo=runQuantumProgram(qp, iqs, 1)

            @test qo.output[2].q == [0;0;0;0;0;1;0;0;;]
        end        
        
        @testset "Output: Big endian 3 qubits CNOT-gate on qubits c=1 and t=2" begin
            # 8. Test 3 qubits, cnotGate on circuitIndices c=1 and t=2, bigEndian, initial state |100>, assert output=|110>
            qc=createQuantumCircuit(3, indexType=circuitIndexBigEndian)
            cnotGate!(qc, 1, 2)
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            iqs=createInitialQubitState(vector, [0. 1. 1.;1. 0. 0.], indexType=byteIndex)
            qo=runQuantumProgram(qp, iqs, 1)

            @test qo.output[2].q == [0;0;0;0;0;0;1;0;;]
        end        
        
        @testset "Output: Big endian 3 qubits CNOT-gate on qubits c=2 and t=3" begin
            # 9. Test 3 qubits, cnotGate on circuitIndices c=2 and t=3, bigEndian, initial state |010>, assert output=|011>
            qc=createQuantumCircuit(3, indexType=circuitIndexBigEndian)
            cnotGate!(qc, 2, 3)
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            iqs=createInitialQubitState(vector, [1. 0. 1.;0. 1. 0.], indexType=byteIndex)
            qo=runQuantumProgram(qp, iqs, 1)

            @test qo.output[2].q == [0;0;0;1;0;0;0;0;;]
        end        
        
        @testset "Output: Big endian 3 qubits CNOT-gate on qubits c=1 and t=3" begin
            # 10. Test 3 qubits, cnotGate on circuitIndices c=1 and t=3, bigEndian, initial state |100>, assert output=|101>
            qc=createQuantumCircuit(3, indexType=circuitIndexBigEndian)
            cnotGate!(qc, 1, 3)
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            iqs=createInitialQubitState(vector, [0. 1. 1.;1. 0. 0.], indexType=byteIndex)
            qo=runQuantumProgram(qp, iqs, 1)

            @test qo.output[2].q == [0;0;0;0;0;1;0;0;;]
        end
    end

    @testset "CNOTReverse gate" verbose = true begin
        @testset "Program: Big endian 2 qubits CNOTReverse-gate on qubit 1" begin
            qc=createQuantumCircuit(2)
            cnotReverseGate!(qc, 1, 2)
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

            @test qp.program[1].U ≈ [1.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 1.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 1.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im] atol=1e-10
        end
    end

    @testset "Swap gate" verbose = true begin
        @testset "Program: Big endian 2 qubits Swap-gate on qubit 1" begin
            qc=createQuantumCircuit(2)
            swapGate!(qc, 1, 2)
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

            @test qp.program[1].U ≈ [1.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 1.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 1.0 + 0.0im] atol=1e-10
        end
    end

    @testset "Phase gate" verbose = true begin
        @testset "Program: Big endian 2 qubits Phase-gate on qubit 1" begin
            qc=createQuantumCircuit(2)
            phaseGate!(qc, 1, 2, 1.5)
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

            @test qp.program[1].U ≈ [1.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 1.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0707372016677029 + 0.9974949866040544im] atol=1e-10
        end
    end

    @testset "UnitaryU gate" verbose = true begin
        @testset "Program: Big endian 1 qubit UnitaryU-gate on qubit 1" begin
            qc=createQuantumCircuit(1)
            unitaryUGate!(qc, [1], createSingleQubitOperationH())
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

            @test qp.program[1].U ≈ [0.7071067811865476 + 0.0im 0.7071067811865475 + 0.0im; 0.7071067811865475 + 0.0im -0.7071067811865476 + 0.0im] atol=1e-10
        end
    end

    @testset "compileToSingleGate" verbose = true begin
        @testset "Function: Big endian 2 qubit compile cnot circuit to single gate" begin
            cnotCircuit=createQuantumCircuit(2)
            hGate!(cnotCircuit, 2)
            controlledUGate!(cnotCircuit, [1], [2], createSingleQubitOperationZ())
            hGate!(cnotCircuit, 2)
            cnotGate=compileToSingleGate(cnotCircuit)

            @test cnotGate.U ≈ [1. 0. 0. 0.; 0. 1. 0. 0.; 0. 0. 0. 1.; 0. 0. 1. 0.] atol=1e-10
        end
        @testset "Program: Big endian 2 qubit compile cnot circuit to single gate" begin
            cnotCircuit=createQuantumCircuit(2)
            hGate!(cnotCircuit, 2)
            controlledUGate!(cnotCircuit, [1], [2], createSingleQubitOperationZ())
            hGate!(cnotCircuit, 2)
            cnotGate=compileToSingleGate(cnotCircuit)

            qc=createQuantumCircuit(2)
            unitaryUGate!(qc, [1,2], cnotGate)
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

            @test qp.program[1].U ≈ [1. 0. 0. 0.; 0. 1. 0. 0.; 0. 0. 0. 1.; 0. 0. 1. 0.] atol=1e-10
        end
    end

    @testset "ExpH gate" verbose = true begin
        @testset "Program: Big endian 1 qubit ExpH-gate on qubit 1" begin
            qc=createQuantumCircuit(1)
            expHGate!(qc, [1], createSingleQubitOperationH().U, 0.3)
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

            @test qp.program[1].U ≈ [0.9553364891256061 - 0.20896434210788312im 0.0 - 0.20896434210788312im; 0.0 - 0.2089643421078831im 0.9553364891256061 + 0.20896434210788312im] atol=1e-10
        end
    end

    @testset "Rx gate" verbose = true begin
        @testset "Program: Big endian 1 qubit Rx-gate on qubit 1" begin
            qc=createQuantumCircuit(1)
            rxGate!(qc, 1, pi/2.)
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

            @test qp.program[1].U ≈ [sqrt(2.)/2+0im -im*sqrt(2.)/2;-im*sqrt(2.)/2 sqrt(2.)/2+0im] atol=1e-10
        end
        @testset "Program: Big endian 1 qubit Rotation-gate on qubit 1" begin
            qc=createQuantumCircuit(1)
            rotationGate!(qc, 1, pi/2., 0., pi/2.)
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

            @test qp.program[1].U ≈ [sqrt(2.)/2+0im -im*sqrt(2.)/2;-im*sqrt(2.)/2 sqrt(2.)/2+0im] atol=1e-10
        end
    end

    @testset "Ry gate" verbose = true begin
        @testset "Program: Big endian 1 qubit Ry-gate on qubit 1" begin
            qc=createQuantumCircuit(1)
            ryGate!(qc, 1, pi/2.)
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

            @test qp.program[1].U ≈ [sqrt(2.)/2+0im -sqrt(2.)/2+0im;sqrt(2.)/2+0im sqrt(2.)/2+0im] atol=1e-10
        end
        @testset "Program: Big endian 1 qubit Rotation-gate on qubit 1" begin
            qc=createQuantumCircuit(1)
            rotationGate!(qc, 1, pi/2., pi/2., pi/2.)
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

            @test qp.program[1].U ≈ [sqrt(2.)/2+0im -sqrt(2.)/2+0im;sqrt(2.)/2+0im sqrt(2.)/2+0im] atol=1e-10
        end
    end

    @testset "Rz gate" verbose = true begin
        @testset "Program: Big endian 1 qubit Rz-gate on qubit 1" begin
            qc=createQuantumCircuit(1)
            rzGate!(qc, 1, pi/2.)
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

            @test qp.program[1].U ≈ [(1.0-im)/sqrt(2) 0.0;0.0 (1.0+im)/sqrt(2)] atol=1e-10
        end
        @testset "Program: Big endian 1 qubit Rotation-gate on qubit 1" begin
            qc=createQuantumCircuit(1)
            rotationGate!(qc, 1, 0., 0., pi/2.)
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

            @test qp.program[1].U ≈ [(1.0-im)/sqrt(2) 0.0;0.0 (1.0+im)/sqrt(2)] atol=1e-10
        end
    end

    @testset "Rotation gate" verbose = true begin
        @testset "Program: Big endian 1 qubit Rotation-gate on qubit 1" begin
            qc=createQuantumCircuit(1)
            rotationGate!(qc, 1, 0.2, 0.4, 0.6)
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

            @test qp.program[1].U ≈ [0.9553364891256061 - 0.2896294776255155im -0.022863063071221844 - 0.054076229366822125im; 0.02286306307122184 - 0.05407622936682213im 0.9553364891256064 + 0.28962947762551555im] atol=1e-10
        end
    end

    @testset "ControlledU gate" verbose = true begin
        @testset "Output: Little endian 3 qubits controlled-X-gate on qubits c=1 and t=2" begin
            # 11. Test 3 qubits, controlledUGate(xGate) on circuitIndices c=1 and t=2, littleEndian, initial state |001>, assert output=|011>
            qc=createQuantumCircuit(3, indexType=circuitIndexLittleEndian)
            controlledUGate!(qc, [1], [2], createSingleQubitOperationX())
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            iqs=createInitialQubitState(vector, [1. 1. 0.;0. 0. 1.], indexType=byteIndex)
            qo=runQuantumProgram(qp, iqs, 1)

            @test qo.output[2].q ≈ [0+0im;0+0im;0+0im;1+0im;0+0im;0+0im;0+0im;0+0im;;] atol=0.1
        end
        
        @testset "Output: Big endian 3 qubits controlled-X-gate on qubits c=1 and t=2" begin
            # 12. Test 3 qubits, controlledUGate(xGate) on circuitIndices c=1 and t=2, bigEndian, initial state |100>, assert output=|110>
            qc=createQuantumCircuit(3, indexType=circuitIndexBigEndian)
            controlledUGate!(qc, [1], [2], createSingleQubitOperationX())
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            iqs=createInitialQubitState(vector, [0. 1. 1.;1. 0. 0.], indexType=byteIndex)
            qo=runQuantumProgram(qp, iqs, 1)

            @test qo.output[2].q ≈ [0+0im;0+0im;0+0im;0+0im;0+0im;0+0im;1+0im;0+0im;;] atol=0.1
        end    

        @testset "Program: Big endian 2 qubits controlled-H-gate on qubits c=1 and t=2" begin
            qc=createQuantumCircuit(2, indexType=circuitIndexBigEndian)
            controlledUGate!(qc, [1], [2], createSingleQubitOperationH())
            controlledHGate!(qc, 1, 2)
            qp=compileQuantumCircuit(qc)

            @test qp.program[1].U ≈ qp.program[2].U atol=0.1
        end

        @testset "Program: Big endian 2 qubits controlled-X-gate on qubits c=1 and t=2" begin
            qc=createQuantumCircuit(2, indexType=circuitIndexBigEndian)
            controlledUGate!(qc, [1], [2], createSingleQubitOperationX())
            controlledXGate!(qc, 1, 2)
            qp=compileQuantumCircuit(qc)

            @test qp.program[1].U ≈ qp.program[2].U atol=0.1
        end

        @testset "Program: Big endian 2 qubits controlled-Y-gate on qubits c=1 and t=2" begin
            qc=createQuantumCircuit(2, indexType=circuitIndexBigEndian)
            controlledUGate!(qc, [1], [2], createSingleQubitOperationY())
            controlledYGate!(qc, 1, 2)
            qp=compileQuantumCircuit(qc)

            @test qp.program[1].U ≈ qp.program[2].U atol=0.1
        end

        @testset "Program: Big endian 2 qubits controlled-Z-gate on qubits c=1 and t=2" begin
            qc=createQuantumCircuit(2, indexType=circuitIndexBigEndian)
            controlledUGate!(qc, [1], [2], createSingleQubitOperationZ())
            controlledZGate!(qc, 1, 2)
            qp=compileQuantumCircuit(qc)

            @test qp.program[1].U ≈ qp.program[2].U atol=0.1
        end

        @testset "Program: Big endian 2 qubits controlled-T-gate on qubits c=1 and t=2" begin
            qc=createQuantumCircuit(2, indexType=circuitIndexBigEndian)
            controlledUGate!(qc, [1], [2], createSingleQubitOperationT())
            controlledTGate!(qc, 1, 2)
            qp=compileQuantumCircuit(qc)

            @test qp.program[1].U ≈ qp.program[2].U atol=0.1
        end

        @testset "Program: Big endian 2 qubits controlled-Td-gate on qubits c=1 and t=2" begin
            qc=createQuantumCircuit(2, indexType=circuitIndexBigEndian)
            controlledUGate!(qc, [1], [2], createSingleQubitOperationTd())
            controlledTdGate!(qc, 1, 2)
            qp=compileQuantumCircuit(qc)

            @test qp.program[1].U ≈ qp.program[2].U atol=0.1
        end

        @testset "Program: Big endian 2 qubits controlled-S-gate on qubits c=1 and t=2" begin
            qc=createQuantumCircuit(2, indexType=circuitIndexBigEndian)
            controlledUGate!(qc, [1], [2], createSingleQubitOperationS())
            controlledSGate!(qc, 1, 2)
            qp=compileQuantumCircuit(qc)

            @test qp.program[1].U ≈ qp.program[2].U atol=0.1
        end

        @testset "Program: Big endian 2 qubits controlled-Sd-gate on qubits c=1 and t=2" begin
            qc=createQuantumCircuit(2, indexType=circuitIndexBigEndian)
            controlledUGate!(qc, [1], [2], createSingleQubitOperationSd())
            controlledSdGate!(qc, 1, 2)
            qp=compileQuantumCircuit(qc)

            @test qp.program[1].U ≈ qp.program[2].U atol=0.1
        end
    end

    @testset "Reflection gate" verbose = true begin
        @testset "Program: Big endian 2 qubits Reflection-gate on qubits 1 and 2" begin
            qc=createQuantumCircuit(2)
            reflectionGate!(qc, [1,2], [1;0;0;0], 1.0*pi)
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

            @test qp.program[1].U ≈ [-1. 0. 0. 0.; 0. 1. 0. 0.; 0. 0. 1. 0.; 0. 0. 0. 1.] atol=1e-10
        end
    end

    @testset "Projection gate" verbose = true begin
        @testset "Program: Big endian 2 qubit Projection-gate on qubits 1 and 2" begin
            qc=createQuantumCircuit(2)
            projectionGate!(qc, [1,2], [1;0;0;0])
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

            @test qp.program[1].U ≈ [1. 0. 0. 0.; 0. 0. 0. 0.; 0. 0. 0. 0.; 0. 0. 0. 0.] atol=1e-10
        end
    end

    @testset "Measure" verbose = true begin
        @testset "Program & Output: Little endian 2 qubits measure on qubit 2" begin
            # 13. Test 2 qubits, measureGate on circuitIndex 2, littleEndian, initial state |01>, assert output=0*
            qc=createQuantumCircuit(2, indexType=circuitIndexLittleEndian)
            measureGate!(qc, [2])
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            iqs=createInitialQubitState(vector, [1. 0.;0. 1.], indexType=byteIndex)
            qo=runQuantumProgram(qp, iqs, 1)

            @test qo.output[2].measured.probability == [1.0;0.0]
            # @test qo.output[2].measured.list == [1 2;3 4]
            # @test qp.program[1].list == [1 2;3 4]
            @test qp.program[1].krausOperators[1].E == [1.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 1.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im]
            @test qp.program[1].krausOperators[2].E == [0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 1.0 + 0.0im]

        end        
        
        @testset "Program & Output: Little endian 2 qubits measure on qubit 1" begin
            # 14. Test 2 qubits, measureGate on circuitIndex 1, littleEndian, initial state |01>, assert output=*1
            qc=createQuantumCircuit(2, indexType=circuitIndexLittleEndian)
            measureGate!(qc, [1])
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            iqs=createInitialQubitState(vector, [1. 0.;0. 1.], indexType=byteIndex)
            qo=runQuantumProgram(qp, iqs, 1)

            @test qo.output[2].measured.probability == [0.0;1.0]
            # @test qo.output[2].measured.list == [1 3;2 4]
            # @test qp.program[1].list == [1 3;2 4]
            @test qp.program[1].krausOperators[1].E == [1.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im]
            @test qp.program[1].krausOperators[2].E == [0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 1.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 1.0 + 0.0im]

        end        
        
        @testset "Program & Output: Big endian 2 qubits measure on qubit 1" begin
            # 15. Test 2 qubits, measureGate on circuitIndex 1, bigEndian, initial state |10>, assert output=1*
            qc=createQuantumCircuit(2, indexType=circuitIndexBigEndian)
            measureGate!(qc, [1])
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            iqs=createInitialQubitState(vector, [0. 1.;1. 0.], indexType=byteIndex)
            qo=runQuantumProgram(qp, iqs, 1)

            @test qo.output[2].measured.probability == [0.0;1.0]
            # @test qo.output[2].measured.list == [1 2;3 4]
            # @test qp.program[1].list == [1 2;3 4]
            @test qp.program[1].krausOperators[1].E == [1.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 1.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im]
            @test qp.program[1].krausOperators[2].E == [0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 1.0 + 0.0im]
        end        
        
        @testset "Program & Output: Big endian 2 qubits measure on qubit 2" begin
            # 16. Test 2 qubits, measureGate on circuitIndex 2, bigEndian, initial state |10>, assert output=*0
            qc=createQuantumCircuit(2, indexType=circuitIndexBigEndian)
            measureGate!(qc, [2])
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            iqs=createInitialQubitState(vector, [0. 1.;1. 0.], indexType=byteIndex)
            qo=runQuantumProgram(qp, iqs, 1)

            @test qo.output[2].measured.probability == [1.0;0.0]
            # @test qo.output[2].measured.list == [1 3;2 4]
            # @test qp.program[1].list == [1 3;2 4]
            @test qp.program[1].krausOperators[1].E == [1.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im]
            @test qp.program[1].krausOperators[2].E == [0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 1.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 1.0 + 0.0im]
        end        
        
        @testset "Program & Output: Little endian 2 qubits measure on qubits 1 and 2" begin
            # 17. Test 2 qubits, measureGate on circuitIndices 1 and 2, littleEndian, initial state |01>, assert output=01
            qc=createQuantumCircuit(2, indexType=circuitIndexLittleEndian)
            measureGate!(qc, [1,2])
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            iqs=createInitialQubitState(vector, [1. 0.;0. 1.], indexType=byteIndex)
            qo=runQuantumProgram(qp, iqs, 1)

            @test qo.output[2].measured.probability == [0.0;1.0;0.0;0.0]
            # @test qo.output[2].measured.list == [1; 2; 3; 4;;]
            # @test qp.program[1].list == [1; 2; 3; 4;;]
            @test qp.program[1].krausOperators[1].E == [1.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im]
            @test qp.program[1].krausOperators[2].E == [0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 1.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im]
            @test qp.program[1].krausOperators[3].E == [0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im]
            @test qp.program[1].krausOperators[4].E == [0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 1.0 + 0.0im]
        end
        
        @testset "Program & Output: Little endian 3 qubits measure on qubits 1 and 2" begin
            # 18. Test 3 qubits, measureGate on circuitIndices 1 and 2, littleEndian, initial state |001>, assert output=*01
            qc=createQuantumCircuit(3, indexType=circuitIndexLittleEndian)
            measureGate!(qc, [1,2])
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            iqs=createInitialQubitState(vector, [1. 1. 0.;0. 0. 1.], indexType=byteIndex)
            qo=runQuantumProgram(qp, iqs, 1)

            @test qo.output[2].measured.probability == [0.0;1.0;0.0;0.0]
            # @test qo.output[2].measured.list == [1 5; 2 6; 3 7; 4 8]
            # @test qp.program[1].list == [1 5; 2 6; 3 7; 4 8]
            @test qp.program[1].krausOperators[1].E == [ 1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im]
            @test qp.program[1].krausOperators[2].E == [ 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im]
            @test qp.program[1].krausOperators[3].E == [ 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im]
            @test qp.program[1].krausOperators[4].E == [ 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im]
        end        
        
        @testset "Program & Output: Little endian 3 qubits measure on qubits 2 and 3" begin
            # 19. Test 3 qubits, measureGate on circuitIndices 2 and 3, littleEndian, initial state |101>, assert output=10*
            qc=createQuantumCircuit(3, indexType=circuitIndexLittleEndian)
            measureGate!(qc, [2,3])
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            iqs=createInitialQubitState(vector, [0. 1. 0.;1. 0. 1.], indexType=byteIndex)
            qo=runQuantumProgram(qp, iqs, 1)

            @test qo.output[2].measured.probability == [0.0;0.0;1.0;0.0]
            # @test qo.output[2].measured.list == [1 2; 3 4; 5 6; 7 8]
            # @test qp.program[1].list == [1 2; 3 4; 5 6; 7 8]
            @test qp.program[1].krausOperators[1].E == [ 1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im]
            @test qp.program[1].krausOperators[2].E == [ 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im]
            @test qp.program[1].krausOperators[3].E == [ 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im]
            @test qp.program[1].krausOperators[4].E == [ 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im]
        end        
        
        @testset "Program & Output: Big endian 2 qubits measure on qubits 1 and 2" begin
            # 20. Test 2 qubits, measureGate on circuitIndices 1 and 2, bigEndian, initial state |01>, assert output=01
            qc=createQuantumCircuit(2, indexType=circuitIndexBigEndian)
            measureGate!(qc, [1,2])
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            iqs=createInitialQubitState(vector, [1. 0.;0. 1.], indexType=byteIndex)
            qo=runQuantumProgram(qp, iqs, 1)

            @test qo.output[2].measured.probability == [0.0;1.0;0.0;0.0]
            # @test qo.output[2].measured.list == [1; 2; 3; 4;;]
            # @test qp.program[1].list == [1; 2; 3; 4;;]
            @test qp.program[1].krausOperators[1].E == [1.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im]
            @test qp.program[1].krausOperators[2].E == [0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 1.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im]
            @test qp.program[1].krausOperators[3].E == [0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im]
            @test qp.program[1].krausOperators[4].E == [0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 1.0 + 0.0im]
        end        

        @testset "Program & Output: Big endian 3 qubits measure on qubits 1 and 2" begin
            # 21. Test 3 qubits, measureGate on circuitIndices 1 and 2, bigEndian, initial state |001>, assert output=00*
            qc=createQuantumCircuit(3, indexType=circuitIndexBigEndian)
            measureGate!(qc, [1,2])
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            iqs=createInitialQubitState(vector, [1. 1. 0.;0. 0. 1.], indexType=byteIndex)
            qo=runQuantumProgram(qp, iqs, 1)

            @test qo.output[2].measured.probability == [1.0;0.0;0.0;0.0]
            # @test qo.output[2].measured.list == [1 2; 3 4; 5 6; 7 8]
            # @test qp.program[1].list == [1 2; 3 4; 5 6; 7 8]
            @test qp.program[1].krausOperators[1].E == [ 1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im]
            @test qp.program[1].krausOperators[2].E == [ 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im]
            @test qp.program[1].krausOperators[3].E == [ 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im]
            @test qp.program[1].krausOperators[4].E == [ 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im]
        end

        @testset "Program & Output: Big endian 3 qubits measure on qubits 2 and 3" begin
            # 22. Test 3 qubits, measureGate on circuitIndices 2 and 3, bigEndian, initial state |101>, assert output=*01        
            qc=createQuantumCircuit(3, indexType=circuitIndexBigEndian)
            measureGate!(qc, [2,3])
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            iqs=createInitialQubitState(vector, [0. 1. 0.;1. 0. 1.], indexType=byteIndex)
            qo=runQuantumProgram(qp, iqs, 1)

            @test qo.output[2].measured.probability == [0.0;1.0;0.0;0.0]
            # @test qo.output[2].measured.list == [1 5; 2 6; 3 7; 4 8]
            # @test qp.program[1].list == [1 5; 2 6; 3 7; 4 8]
            @test qp.program[1].krausOperators[1].E == [ 1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im]
            @test qp.program[1].krausOperators[2].E == [ 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im]
            @test qp.program[1].krausOperators[3].E == [ 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im]
            @test qp.program[1].krausOperators[4].E == [ 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
                                                       0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im]
        end
        
        @testset "Output: Measure 1 qubit along arbitrary axis" begin
            theta=pi/4
            phi=pi/8
            
            q=[1;1]
            q=q/norm(q)
            rho=q*q'
            iqs=DensityState(rho)
            
            qc=createQuantumCircuit(1)
            measureGate!(qc, [1],[sigmaN(theta,phi)])
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            qo1=runQuantumProgram(qp, iqs, 1)
            
            qc=createQuantumCircuit(1)
            u3Gate!(qc, 1, theta, phi, pi-phi)
            measureGate!(qc, [1])
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            qo2=runQuantumProgram(qp, iqs, 1)
            
            Ma=createSingleQubitBlochDensityState(theta, phi).rho
            prob_a=tr(Ma*rho)

            @test qo1.output[2].measured.probability ≈ qo2.output[3].measured.probability atol=1e-10
            @test qo1.output[2].measured.probability[1] ≈ prob_a atol=1e-10
            @test qo2.output[3].measured.probability[1] ≈ prob_a atol=1e-10
        end

        @testset "Output: Measure 2 qubits along arbitrary axes" begin
            q=[1;0;0;2]
            q=q/norm(q)
            rho=q*q'
            iqs=DensityState(rho)
            
            n_theta=0.6
            n_phi=0.2
            m_theta=0.3
            m_phi=0.
            
            qc=createQuantumCircuit(2)
            measureGate!(qc, [1,2],[sigmaN(n_theta, n_phi), sigmaN(m_theta, m_phi)])
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            qo1=runQuantumProgram(qp, iqs, 1)
            
            qc=createQuantumCircuit(2)
            u3Gate!(qc, 1, n_theta, n_phi, pi-n_phi)
            u3Gate!(qc, 2, m_theta, m_phi, pi-m_phi)
            measureGate!(qc, [1,2])
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            qo2=runQuantumProgram(qp, iqs, 1)
            
            Ma=createSingleQubitBlochDensityState(n_theta, n_phi).rho
            Mb=createSingleQubitBlochDensityState(m_theta, m_phi).rho
            Ma_Mb=tensorProduct(Ma,Mb)
            prob_a_and_b=tr(Ma_Mb*rho)
            
            @test qo1.output[2].measured.probability ≈ qo2.output[3].measured.probability atol=1e-10
            @test qo1.output[2].measured.probability[1] ≈ prob_a_and_b atol=1e-10
            @test qo2.output[3].measured.probability[1] ≈ prob_a_and_b atol=1e-10
        end

        @testset "Output: Measure 2 qubits using conditional probabilities along arbitrary axes" begin
            q=[1;0;0;2]
            q=q/norm(q)
            rho=q*q'
            iqs=DensityState(rho)
            
            n_theta=0.6
            n_phi=0.2
            m_theta=0.3
            m_phi=0.
            
            Ma=createSingleQubitBlochDensityState(n_theta, n_phi).rho
            Mb=createSingleQubitBlochDensityState(m_theta, m_phi).rho
            
            Ma_Mb=tensorProduct(Ma,Mb)
            prob_a_and_b=real(tr(Ma_Mb*rho))
            
            @test prob_a_and_b ≈ 0.21272491249215775 atol=1e-10

            # repeat multiple time to get 1st qubit to 0 (small prob)
            qc=createQuantumCircuit(2)
            # u3Gate!(qc, 2, m_theta, m_phi, pi-m_phi)
            # measureGate!(qc, [2])
            measureGate!(qc, [2], [sigmaN(m_theta, m_phi)])
            # u3Gate!(qc, 1, n_theta, n_phi, pi-n_phi)
            # measureGate!(qc, [1])
            measureGate!(qc, [1], [sigmaN(n_theta, n_phi)])
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            found0=false
            nrTrials=0
            qo=runQuantumProgram(qp, iqs, 1)
            while (!found0 && nrTrials<100)
                qo=runQuantumProgram(qp, iqs, 1)
                # found0=(qo.output[3].measured.outcome==0)
                found0=(qo.output[2].measured.outcome==0)
                nrTrials+=1
                # print(nrTrials)
            end
            # @test qo.output[5].measured.probability[1]*qo.output[3].measured.probability[1] ≈ 0.21272491249215775 atol=1e-10
            @test qo.output[3].measured.probability[1]*qo.output[2].measured.probability[1] ≈ 0.21272491249215775 atol=1e-10
            
            I_Mb=tensorProduct(Matrix{Float64}(I, 2, 2), Mb)
            prob_b=tr(I_Mb*rho)
            rho_collapsed=I_Mb*rho*I_Mb/prob_b
            rho_a=partialTrace(rho_collapsed, traceIndex=2)
            prob_a_cond_b=tr(Ma*rho_a)
            prob_a_and_b=real(prob_a_cond_b*prob_b)
            @test prob_a_and_b ≈ 0.21272491249215775 atol=1e-10
            
            # repeat multiple time to get 1st qubit to 0 (small prob)
            qc=createQuantumCircuit(2)
            # u3Gate!(qc, 1, n_theta, n_phi, pi-n_phi)
            # measureGate!(qc, [1])
            measureGate!(qc, [1], [sigmaN(n_theta, n_phi)])
            # u3Gate!(qc, 2, m_theta, m_phi, pi-m_phi)
            # measureGate!(qc, [2])
            measureGate!(qc, [2], [sigmaN(m_theta, m_phi)])
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            found0=false
            nrTrials=0
            while (!found0 && nrTrials<100)
                qo=runQuantumProgram(qp, iqs, 1)
                # found0=(qo.output[3].measured.outcome==0)
                found0=(qo.output[2].measured.outcome==0)
                nrTrials+=1
                # print(nrTrials)
            end
            # @test qo.output[5].measured.probability[1]*qo.output[3].measured.probability[1] ≈ 0.21272491249215775 atol=1e-10
            @test qo.output[3].measured.probability[1]*qo.output[2].measured.probability[1] ≈ 0.21272491249215775 atol=1e-10
            
            Ma_I=tensorProduct(Ma, Matrix{Float64}(I, 2, 2))
            prob_a=tr(Ma_I*rho)
            rho_collapsed=Ma_I*rho*Ma_I/prob_a
            rho_b=partialTrace(rho_collapsed, traceIndex=1)
            prob_b_cond_a=tr(Mb*rho_b)
            prob_a_and_b=real(prob_b_cond_a*prob_a)
            @test prob_a_and_b ≈ 0.21272491249215775 atol=1e-10
        end

        @testset "Output: Compare PVM and POVM measurement in 1-qubit circuit" begin
            # along arbitrary axes
            theta=pi/4
            phi=pi/8
            
            q=[1;1]
            q=q/norm(q)
            rho=q*q'
            iqs=DensityState(rho)
            
            qc=createQuantumCircuit(1)
            measureGate!(qc, [1],[sigmaN(theta,phi)])
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            qo1=runQuantumProgram(qp, iqs, 1)
            
            qc=createQuantumCircuit(1)
            measureGate!(qc, [1], generateKrausOperatorsForPVMMeasurement(1, [1], [sigmaN(theta,phi)]))
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            qo2=runQuantumProgram(qp, iqs, 1)

            Ma=createSingleQubitBlochDensityState(theta, phi).rho
            prob_a=tr(Ma*rho)

            @test qo1.output[2].measured.probability ≈ qo2.output[2].measured.probability atol=1e-10
            @test qo1.output[2].measured.probability[1] ≈ prob_a atol=1e-10
            @test qo2.output[2].measured.probability[1] ≈ prob_a atol=1e-10
        end

        @testset "Output: Compare PVM and POVM measurement of 1st qubit in 2-qubit circuit" begin
            # along arbitrary axes
            theta=pi/4
            phi=pi/8
            
            q=[1;0;0;2]
            q=q/norm(q)
            rho=q*q'
            iqs=DensityState(rho)
            
            qc=createQuantumCircuit(2)
            measureGate!(qc, [1],[sigmaN(theta,phi)])
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            qo1=runQuantumProgram(qp, iqs, 1)
            
            qc=createQuantumCircuit(2)
            measureGate!(qc, [1], generateKrausOperatorsForPVMMeasurement(1, [1], [sigmaN(theta,phi)]))
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            qo2=runQuantumProgram(qp, iqs, 1)

            Ma=createSingleQubitBlochDensityState(theta, phi).rho
            prob_a=real(tr(kron(Ma, Matrix{ComplexF64}(I, 2, 2))*rho))

            @test qo1.output[2].measured.probability ≈ qo2.output[2].measured.probability atol=1e-10
            @test qo1.output[2].measured.probability[1] ≈ prob_a atol=1e-10
            @test qo2.output[2].measured.probability[1] ≈ prob_a atol=1e-10
        end

        @testset "Output: Compare PVM and POVM measurement of 2nd qubit in 2-qubit circuit" begin
            # along arbitrary axes
            theta=pi/4
            phi=pi/8
            
            q=[1;0;0;2]
            q=q/norm(q)
            rho=q*q'
            iqs=DensityState(rho)
            
            qc=createQuantumCircuit(2)
            measureGate!(qc, [2],[sigmaN(theta,phi)])
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            qo1=runQuantumProgram(qp, iqs, 1)
            
            qc=createQuantumCircuit(2)
            measureGate!(qc, [2], generateKrausOperatorsForPVMMeasurement(1, [1], [sigmaN(theta,phi)]))
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            qo2=runQuantumProgram(qp, iqs, 1)

            Ma=createSingleQubitBlochDensityState(theta, phi).rho
            prob_a=real(tr(kron(Matrix{ComplexF64}(I, 2, 2), Ma)*rho))

            @test qo1.output[2].measured.probability ≈ qo2.output[2].measured.probability atol=1e-10
            @test qo1.output[2].measured.probability[1] ≈ prob_a atol=1e-10
            @test qo2.output[2].measured.probability[1] ≈ prob_a atol=1e-10
        end

        @testset "Output: Compare PVM and POVM measurement of 2 qubits" begin
            q=[1;0;0;2]
            q=q/norm(q)
            rho=q*q'
            iqs=DensityState(rho)
            
            # along arbitrary axes
            n_theta=0.6
            n_phi=0.2
            m_theta=0.3
            m_phi=0.
            
            qc=createQuantumCircuit(2)
            measureGate!(qc, [1,2],[sigmaN(n_theta, n_phi), sigmaN(m_theta, m_phi)])
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            qo1=runQuantumProgram(qp, iqs, 1)

            qc=createQuantumCircuit(2)
            measureGate!(qc, [1,2], generateKrausOperatorsForPVMMeasurement(2, [1,2], [sigmaN(n_theta, n_phi), sigmaN(m_theta, m_phi)]))
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            qo2=runQuantumProgram(qp, iqs, 1)

            Ma=createSingleQubitBlochDensityState(n_theta, n_phi).rho
            Mb=createSingleQubitBlochDensityState(m_theta, m_phi).rho
            Ma_Mb=tensorProduct(Ma,Mb)
            prob_a_and_b=tr(Ma_Mb*rho)
            
            @test qo1.output[2].measured.probability ≈ qo2.output[2].measured.probability atol=1e-10
            @test qo1.output[2].measured.probability[1] ≈ prob_a_and_b atol=1e-10
            @test qo2.output[2].measured.probability[1] ≈ prob_a_and_b atol=1e-10
        end

        @testset "Output: non-orthogonal (Trine) POVM measurement of 1 qubit with 3 Kraus operators" begin
            q0=createInitialQubitState(vector, [1. 0.]).q
            E0=sqrt(2/3)*q0*q0'
            q1=createInitialQubitState(vector, [-0.5 sqrt(3)/2]).q
            E1=sqrt(2/3)*q1*q1'
            q2=createInitialQubitState(vector, [-0.5 -sqrt(3)/2]).q
            E2=sqrt(2/3)*q2*q2'

            qc=createQuantumCircuit(1)
            measureGate!(qc, [1], Vector([KrausOperator(E0,"0"), KrausOperator(E1,"1"), KrausOperator(E2,"2")]))
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

            theta=0
            phi=0

            iqs = createInitialQubitState(vector, [theta phi], blochRepresentation=true)
            qo1=runQuantumProgram(qp, iqs, 1)

            iqs = createInitialQubitState(density, [[1., [theta phi]]], blochRepresentation=true)
            qo2=runQuantumProgram(qp, iqs, 1)

            @test qo1.output[2].measured.probability ≈ [0.6666666666666666, 0.16666666666666663, 0.16666666666666663] atol=1e-10
            @test qo2.output[2].measured.probability ≈ [0.6666666666666666, 0.16666666666666663, 0.16666666666666663] atol=1e-10
        end

        @testset "Output: non-orthogonal POVM measurement of 1 qubit with 2 Kraus operators" begin
            lambda=0.6
            q0=createInitialQubitState(vector, [1. 0.]).q
            q1=createInitialQubitState(vector, [0. 1.]).q
            E0=sqrt(lambda)*q0*q0'
            E1=sqrt(1-lambda)*q0*q0'+q1*q1'

            qc=createQuantumCircuit(1)
            measureGate!(qc, [1], Vector([KrausOperator(E0,"0"), KrausOperator(E1,"1")]))
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)

            theta=pi/2
            phi=0

            iqs = createInitialQubitState(vector, [theta phi], blochRepresentation=true)
            qo1=runQuantumProgram(qp, iqs, 1)

            iqs = createInitialQubitState(density, [[1., [theta phi]]], blochRepresentation=true)
            qo2=runQuantumProgram(qp, iqs, 1)

            @test qo1.output[2].measured.probability ≈ [lambda/2, 1-lambda/2] atol=1e-10
            @test qo2.output[2].measured.probability ≈ [lambda/2, 1-lambda/2] atol=1e-10
        end

        @testset "Output: Measurement on superposition state yields statistical distribution" begin
            qc=createQuantumCircuit(1)
            hGate!(qc, 1)
            measureGate!(qc, [1])
            qp=compileQuantumCircuit(qc)

            nrOfShots=100
            p=0.5
            E=nrOfShots*p
            SD=sqrt(nrOfShots*p*(1-p))

            iqs=createInitialQubitState(vector, [1. 0.;])
            qo=runQuantumProgram(qp, iqs, nrOfShots)
            sumOutcome=sum([qo.output[3,k].measured.outcome for k in 1:nrOfShots])
            # test that 99% of the outcomes are within 3-sigma of the expectation
            @test sumOutcome > (E-3*SD)
            @test sumOutcome < (E+3*SD)

            iqs=createInitialQubitState(density, [[1.,[1. 0.;]]])
            qo=runQuantumProgram(qp, iqs, nrOfShots)
            sumOutcome=sum([qo.output[3,k].measured.outcome for k in 1:nrOfShots])
            # test that 99% of the outcomes are within 3-sigma of the expectation
            @test sumOutcome > (E-3*SD)
            @test sumOutcome < (E+3*SD)
        end
    end

    @testset "QFT" verbose = true begin
        @testset "Output: Big endian and toggled little endian 2 qubit QFT" begin
            # 23. Test 2 qubits qftGate, initial state endian dependent velocity, assert output big endian and toggled little endian is equal
            numberOfQubits=2
            
            velocity=1
            qcb=createQuantumCircuit(numberOfQubits)
            qftGate!(qcb, Vector(1:numberOfQubits))
            qpb=compileQuantumCircuit(qcb)
            iqs=createInitialQubitState(vector, createByteIndexVector(velocity, numberOfQubits), indexType=byteIndex)
            qob=runQuantumProgram(qpb, iqs, 1)
            
            velocity=2
            qcl=createQuantumCircuit(numberOfQubits, indexType=circuitIndexLittleEndian)
            qftGate!(qcl, Vector(1:numberOfQubits))
            qpl=compileQuantumCircuit(qcl)
            iqs=createInitialQubitState(vector, createByteIndexVector(velocity, numberOfQubits), indexType=byteIndex)
            qol=runQuantumProgram(qpl, iqs, 1)
            
            @test vec(qob.output[2].q) ≈ [0.5, 0.5im, -0.5, -0.5im] atol=0.1
            @test (qob.output[2].q) == toggleQubitOrdeningFor(qol.output[2]).q
            @test qpb.program[1].U == toggleQubitOrdeningFor(qpl.program[1]).U
        end

        @testset "Output: Big endian and toggled little endian 3 qubit QFT" begin
            # 24. Test 3 qubits qftGate, initial state endian dependent velocity, assert output big endian and toggled little endian is equal
            numberOfQubits=3
            
            velocity=3
            qcb=createQuantumCircuit(numberOfQubits)
            qftGate!(qcb, Vector(1:numberOfQubits))
            qpb=compileQuantumCircuit(qcb)
            iqs=createInitialQubitState(vector, createByteIndexVector(velocity, numberOfQubits), indexType=byteIndex)
            qob=runQuantumProgram(qpb, iqs, 1)
            
            velocity=6
            qcl=createQuantumCircuit(numberOfQubits, indexType=circuitIndexLittleEndian)
            qftGate!(qcl, Vector(1:numberOfQubits))
            qpl=compileQuantumCircuit(qcl)
            iqs=createInitialQubitState(vector, createByteIndexVector(velocity, numberOfQubits), indexType=byteIndex)
            qol=runQuantumProgram(qpl, iqs, 1)
            
            @test vec(qob.output[2].q) ≈ [0.3536,-0.25+0.25im,-0.3536im,0.25+0.25im,-0.3536,0.25-0.25im,0.3536im,-0.25-0.25im] atol=0.0001
            @test (qob.output[2].q) == toggleQubitOrdeningFor(qol.output[2]).q
            @test qpb.program[1].U == toggleQubitOrdeningFor(qpl.program[1]).U
        end
    end

    @testset "IQFT" verbose = true begin
        @testset "Output: Big endian and toggled little endian 2 qubit IQFT" begin
            # 25. Test 2 qubits iqftGate, initial state endian dependent velocity, assert output big endian and toggled little endian is equal        
            numberOfQubits=2
        
            velocity=1 # 01=1
            qcb=createQuantumCircuit(numberOfQubits)
            iqftGate!(qcb, Vector(1:numberOfQubits))
            qpb=compileQuantumCircuit(qcb)
            iqs=createInitialQubitState(vector, createByteIndexVector(velocity, numberOfQubits), indexType=byteIndex)
            qob=runQuantumProgram(qpb, iqs, 1)
            
            velocity=2 # toggled(01=1)=10=2
            qcl=createQuantumCircuit(numberOfQubits, indexType=circuitIndexLittleEndian)
            iqftGate!(qcl, Vector(1:numberOfQubits))
            qpl=compileQuantumCircuit(qcl)
            iqs=createInitialQubitState(vector, createByteIndexVector(velocity, numberOfQubits), indexType=byteIndex)
            qol=runQuantumProgram(qpl, iqs, 1)
            
            @test vec(qob.output[2].q) ≈ [0.5,-0.5im,-0.5,0.5im] atol=0.1
            @test (qob.output[2].q) == toggleQubitOrdeningFor(qol.output[2]).q
            @test qpb.program[1].U == toggleQubitOrdeningFor(qpl.program[1]).U
        end

        @testset "Output: Big endian and toggled little endian 3 qubit IQFT" begin
            # 26. Test 3 qubits iqftGate, initial state endian dependent velocity, assert output big endian and toggled little endian is equal
            numberOfQubits=3
            
            velocity=3 # 011=3
            qcb=createQuantumCircuit(numberOfQubits)
            iqftGate!(qcb, Vector(1:numberOfQubits))
            qpb=compileQuantumCircuit(qcb)
            iqs=createInitialQubitState(vector, createByteIndexVector(velocity, numberOfQubits), indexType=byteIndex)
            qob=runQuantumProgram(qpb, iqs, 1)
            
            velocity=6 # toggled(011=3)=110=6
            qcl=createQuantumCircuit(numberOfQubits, indexType=circuitIndexLittleEndian)
            iqftGate!(qcl, Vector(1:numberOfQubits))
            qpl=compileQuantumCircuit(qcl)
            iqs=createInitialQubitState(vector, createByteIndexVector(velocity, numberOfQubits), indexType=byteIndex)
            qol=runQuantumProgram(qpl, iqs, 1)

            @test vec(qob.output[2].q) ≈ [0.3536,-0.25-0.25im,0.3536im,0.25-0.25im,-0.3536,0.25+0.25im,-0.3536im,-0.25+0.25im] atol=0.0001
            @test (qob.output[2].q) == toggleQubitOrdeningFor(qol.output[2]).q
            @test qpb.program[1].U == toggleQubitOrdeningFor(qpl.program[1]).U
        end
    end

    @testset "Quantum phase estimation" verbose = true begin
        @testset "Output: Big endian 3 qubit quantum phase estimation of phase gate" begin
            # 27. Test big endian 3 qubit accurate quantum phase estimation of controlled phase gate with theta=0.7*2*pi
            acc=3
            theta=.7*2*pi
            qcb=createQuantumCircuit(acc+1)
            xGate!(qcb,acc+1) # to create eigenstate |1> with eigenvalue 0.7*2pi from initial state |0>
            qpeGate!(qcb, Vector(1:acc), [acc+1], createSingleQubitOperationU1(theta))
            measureGate!(qcb, Vector(1:acc))
            qpb=compileQuantumCircuit(qcb; optimizeNumberOfSteps=true)
            iqs=createInitialQubitState(vector, repeat([1. 0.], outer=[acc+1,1]), indexType=circuitIndexBigEndian)
            qob=runQuantumProgram(qpb, iqs, 1000)
            # probeMeasureOutcome(qob, 3, "estimated phase")
            outcome=[qob.output[3,k].measured.outcome for k in 1:qob.numberOfShots]
            orderedOutcome=[count(outcome.==k) for k=0:length(qob.output[3,1].measured.labels)-1]
        
            @test (findmax(orderedOutcome)[2]-1)/(2^acc) == 0.75
        end

        @testset "Output: Little endian 3 qubit quantum phase estimation of phase gate" begin    
            # 28. Test little endian 3 qubit accurate quantum phase estimation of controlled phase gate with theta=0.7*2*pi
            acc=3
            theta=.7*2*pi
            qcl=createQuantumCircuit(acc+1, indexType=circuitIndexLittleEndian)
            xGate!(qcl,acc+1) # to create eigenstate |1> with eigenvalue 0.7*2pi from initial state |0>
            qpeGate!(qcl, Vector(1:acc), [acc+1], createSingleQubitOperationU1(theta))
            measureGate!(qcl, Vector(1:acc))
            qpl=compileQuantumCircuit(qcl; optimizeNumberOfSteps=true)
            iqs=createInitialQubitState(vector, repeat([1. 0.], outer=[acc+1,1]), indexType=circuitIndexBigEndian)
            qol=runQuantumProgram(qpl, iqs, 1000)
            # probeMeasureOutcome(qol, 3, "estimated phase")
            outcome=[qol.output[3,k].measured.outcome for k in 1:qol.numberOfShots]
            orderedOutcome=[count(outcome.==k) for k=0:length(qol.output[3,1].measured.labels)-1]
        
            @test (findmax(orderedOutcome)[2]-1)/(2^acc) == 0.75
        end

        @testset "Output: Big endian 4 qubit quantum phase estimation of Td gate" begin
            # 29. Test big endian 4 qubit accurate quantum phase estimation of controlled Td gate
            acc=4
            qcb=createQuantumCircuit(acc+1)
            xGate!(qcb,acc+1) # to create eigenstate |1> with eigenvalue 0.875*2pi from initial state |0>
            qpeGate!(qcb, Vector(1:acc), [acc+1], createSingleQubitOperationTd())
            measureGate!(qcb, Vector(1:acc))
            qpb=compileQuantumCircuit(qcb; optimizeNumberOfSteps=true)
            iqs=createInitialQubitState(vector, repeat([1. 0.], outer=[acc+1,1]), indexType=circuitIndexBigEndian)
            qob=runQuantumProgram(qpb, iqs, 1000)
            # probeMeasureOutcome(qob, 3, "estimated phase")
            outcome=[qob.output[3,k].measured.outcome for k in 1:qob.numberOfShots]
            orderedOutcome=[count(outcome.==k) for k=0:length(qob.output[3,1].measured.labels)-1]
        
            @test (findmax(orderedOutcome)[2]-1)/(2^acc) == 0.875
        end

        @testset "Output: Little endian 4 qubit quantum phase estimation of Td gate" begin
            # 30. Test little endian 4 qubit accurate quantum phase estimation of controlled Td gate
            acc=4
            qcl=createQuantumCircuit(acc+1, indexType=circuitIndexLittleEndian)
            xGate!(qcl,acc+1) # to create eigenstate |1> with eigenvalue 0.875*2pi from initial state |0>
            qpeGate!(qcl, Vector(1:acc), [acc+1], createSingleQubitOperationTd())
            measureGate!(qcl, Vector(1:acc))
            qpl=compileQuantumCircuit(qcl; optimizeNumberOfSteps=true)
            iqs=createInitialQubitState(vector, repeat([1. 0.], outer=[acc+1,1]), indexType=circuitIndexBigEndian)
            qol=runQuantumProgram(qpl, iqs, 1000)
            # probeMeasureOutcome(qol, 3, "estimated phase")
            outcome=[qol.output[3,k].measured.outcome for k in 1:qol.numberOfShots]
            orderedOutcome=[count(outcome.==k) for k=0:length(qol.output[3,1].measured.labels)-1]

            @test (findmax(orderedOutcome)[2]-1)/(2^acc) == 0.875
        end
    end

    @testset "Quantum no-signaling theorem" verbose = true begin
        @testset "Output: Reduced operator of qubit 2 independent of ignorant measurement of qubit 1" begin
            # 31. Check quantum no-signaling theorem: reduced operator of 2nd qubit not influenced by ignorant measurement of 1st qubit.
            qc=createQuantumCircuit(2)
            hGate!(qc, 1)
            cnotGate!(qc, 1, 2) # entangled
            measureGate!(qc, [1], forgetOutcome=true)
            qp=compileQuantumCircuit(qc)
            iqs=createInitialQubitState(density, [[1.00,[1. 0.;1. 0.]],[0.00,[0. 1.;0. 1.]]])
            qo=runQuantumProgram(qp, iqs, 1)
            reducedrho2before=partialTrace(qo.output[3].rho, traceIndex=1)
            reducedrho2after=partialTrace(qo.output[4].rho, traceIndex=1)
            
            @test reducedrho2before ≈ reducedrho2after atol=0.1
        end
        
        @testset "Output: Reduced operator of qubit 2 independent of quantum channel on qubits 1 and 2(I)" begin
            # 32. Check quantum no-signaling theorem: reduced operator of 2nd qubit not influenced by quantum channel gate operating on 1st qubit.
            qc=createQuantumCircuit(2)
            hGate!(qc, 1)
            cnotGate!(qc, 1, 2) # entangled
            quantumChannelGate!(qc, [1,2], generateKrausOperatorsForPVMMeasurement(2, [1], Sigmas()))
            qp=compileQuantumCircuit(qc)
            iqs=createInitialQubitState(density, [[1.00,[1. 0.;1. 0.]],[0.00,[0. 1.;0. 1.]]])
            qo=runQuantumProgram(qp, iqs, 1)
            reducedrho2before=partialTrace(qo.output[3].rho, traceIndex=1)
            reducedrho2after=partialTrace(qo.output[4].rho, traceIndex=1)
            
            @test reducedrho2before ≈ reducedrho2after atol=0.1
        end    

        @testset "Output: Reduced operator of qubit 2 independent of quantum channel on qubit 1" begin
            # 32. Check quantum no-signaling theorem: reduced operator of 2nd qubit not influenced by quantum channel gate operating on 1st qubit.
            qc=createQuantumCircuit(2)
            hGate!(qc, 1)
            cnotGate!(qc, 1, 2) # entangled
            quantumChannelGate!(qc, [1], generateKrausOperatorsForPVMMeasurement(1, [1], Sigmas()))
            qp=compileQuantumCircuit(qc)
            iqs=createInitialQubitState(density, [[1.00,[1. 0.;1. 0.]],[0.00,[0. 1.;0. 1.]]])
            qo=runQuantumProgram(qp, iqs, 1)
            reducedrho2before=partialTrace(qo.output[3].rho, traceIndex=1)
            reducedrho2after=partialTrace(qo.output[4].rho, traceIndex=1)
            
            @test reducedrho2before ≈ reducedrho2after atol=0.1
        end    

        @testset "Function: Prepare Kraus operators for measure and forget outcome" begin
            @test generateKrausOperatorsForPVMMeasurement(2, [2], Sigmas())[1].E == [1. 0. 0. 0.;0. 0. 0. 0.;0. 0. 1. 0.;0. 0. 0. 0.]
            @test generateKrausOperatorsForPVMMeasurement(2, [2], Sigmas())[2].E == [0. 0. 0. 0.;0. 1. 0. 0.;0. 0. 0. 0.;0. 0. 0. 1.]

            @test generateKrausOperatorsForPVMMeasurement(2, [1], Sigmas(), indexType = circuitIndexLittleEndian)[1].E == [1. 0. 0. 0.;0. 0. 0. 0.;0. 0. 1. 0.;0. 0. 0. 0.]
            @test generateKrausOperatorsForPVMMeasurement(2, [1], Sigmas(), indexType = circuitIndexLittleEndian)[2].E == [0. 0. 0. 0.;0. 1. 0. 0.;0. 0. 0. 0.;0. 0. 0. 1.]

            @test generateKrausOperatorsForPVMMeasurement(2, [0], Sigmas(), zeroBasedNumbering = true, indexType = circuitIndexLittleEndian)[1].E == [1. 0. 0. 0.;0. 0. 0. 0.;0. 0. 1. 0.;0. 0. 0. 0.]
            @test generateKrausOperatorsForPVMMeasurement(2, [0], Sigmas(), zeroBasedNumbering = true, indexType = circuitIndexLittleEndian)[2].E == [0. 0. 0. 0.;0. 1. 0. 0.;0. 0. 0. 0.;0. 0. 0. 1.]
        end
    end

    @testset "Reduced density operators" verbose = true begin
        @testset "Output: Purity of reduced density operator of qubit 1 of entangled qubit pair" begin
            # 33. Create entangled state of 2 bits, take partial trace on 2nd bit to obtain reduced density operator and evaluate purity.
            qc=createQuantumCircuit(2)
            hGate!(qc, 1)
            cnotGate!(qc, 1, 2)
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            iqs=createInitialQubitState(density, [[1.00,[1. 0.;1. 0.]],[0.00,[0. 1.;0. 1.]]])
            qo=runQuantumProgram(qp, iqs, 1)
            rho=partialTrace(qo.output[2].rho, traceIndex=2)
            (entropy, purity, purity_norm, purity_str)=calculateEntropyAndPurity(rho)
            
            @test entropy ≈ 0.6931471805599453 atol=1e-10
            @test purity ≈ 0.5 atol=1e-10
            @test purity_norm ≈ 0.0 atol=1e-10
            @test purity_str == "impure, totally mixed"
        end
        
        @testset "Output: Purity of reduced density operator of qubit 1 of non-entangled qubit pair" begin
            # 34. Create non-entangled state of 2 bits, take partial trace on 2nd bit to obtain reduced density operator and evaluate purity.    
            qc=createQuantumCircuit(2)
            hGate!(qc, 1)
            hGate!(qc, 2)
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            iqs=createInitialQubitState(density, [[1.00,[1. 0.;1. 0.]],[0.00,[0. 1.;0. 1.]]])
            qo=runQuantumProgram(qp, iqs, 1)
            rho=partialTrace(qo.output[2].rho, traceIndex=2)
            (entropy, purity, purity_norm, purity_str)=calculateEntropyAndPurity(rho)
            
            @test entropy ≈ 0.0 atol=1e-10
            @test purity ≈ 1.0 atol=1e-10
            @test purity_norm ≈ 1.0 atol=1e-10
            @test purity_str == "pure"
        end

        @testset "Output: Partial trace on 3 qubits where qubit 1 & 2 entangled and 3 in superposition" begin
            qc=createQuantumCircuit(3)
            hGate!(qc, 1)
            cnotGate!(qc, 1, 2)
            hGate!(qc, 3)
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            iqs=createInitialQubitState(density, [[1.00,[1. 0.;1. 0.;1. 0.]]])
            qo=runQuantumProgram(qp, iqs, 1)
            rho=partialTrace(qo.output[2].rho, traceIndex=1)

            @test partialTrace(qo.output[2].rho, traceIndex=3) ≈ [.5 0. 0. .5;0. 0. 0. 0.;0. 0. 0. 0.;.5 0. 0. .5] atol=1e-10
            @test partialTrace(qo.output[2].rho, traceIndex=1) ≈ [.25 .25 0. 0.;.25 .25 0. 0.;0. 0. .25 .25;0. 0. .25 .25] atol=1e-10
            @test partialTrace(qo.output[2].rho, traceIndex=2) ≈ [.25 .25 0. 0.;.25 .25 0. 0.;0. 0. .25 .25;0. 0. .25 .25] atol=1e-10
        end

        @testset "Function: Partial trace" begin
            rho=[1 2 3 4;5 6 7 8;9 10 11 12;13 14 15 16]
            
            @test partialTrace(rho, traceIndex=1) ≈ [12 14;20 22] atol=1e-10
            @test partialTrace(rho, traceIndex=2) ≈ [7 11;23 27] atol=1e-10

            @test partialTrace(rho, traceIndex=2, indexType = circuitIndexLittleEndian) ≈ [12 14;20 22] atol=1e-10
            @test partialTrace(rho, traceIndex=1, indexType = circuitIndexLittleEndian) ≈ [7 11;23 27] atol=1e-10

            @test partialTrace(rho, traceIndex=1, zeroBasedNumbering = true, indexType = circuitIndexLittleEndian) ≈ [12 14;20 22] atol=1e-10
            @test partialTrace(rho, traceIndex=0, zeroBasedNumbering = true, indexType = circuitIndexLittleEndian) ≈ [7 11;23 27] atol=1e-10
        end
    end

    @testset "CHSH inequalities" verbose = true begin
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
        theta=0.8
        iqs=createInitialQubitState(vector, [1. 0.;1. 0.])

        expZZ=evaluateCorrelationInTwoQubitSystem(iqs, theta, [sigmaZ(),sigmaZ()], useProbability, nrShots)
        expZX=evaluateCorrelationInTwoQubitSystem(iqs, theta, [sigmaZ(),sigmaX()], useProbability, nrShots)
        expXZ=evaluateCorrelationInTwoQubitSystem(iqs, theta, [sigmaX(),sigmaZ()], useProbability, nrShots)
        expXX=evaluateCorrelationInTwoQubitSystem(iqs, theta, [sigmaX(),sigmaX()], useProbability, nrShots)

        chsh1=(expZZ-expZX+expXZ+expXX)
        chsh2=(expZZ+expZX-expXZ+expXX)
        @test chsh1 ≈ 2.828125600493375 atol=1e-10
        @test chsh2 ≈ -0.04129876310471481 atol=1e-10
    end

    @testset "Bell inequalities" verbose = true begin
        function evaluateCorrelationInTwoQubitSystem(iqs, theta_1, phi_1, theta_2, phi_2)
            qc=createQuantumCircuit(2)
            u3Gate!(qc, 1, theta_1, phi_1, pi-phi_1)
            u3Gate!(qc, 2, theta_2, phi_2, pi-phi_2)
            measureGate!(qc, [1,2])
            qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            qo=runQuantumProgram(qp, iqs, 1)
            corr=0.25*calculateExpectationValueOfProductOfMeasuredQubits(qo, 3, useProbability)
            return corr
        end
        
        theta = 0.9
        iqs=createDoubleQubitWernerDensityState(1)
        corr1_ab=evaluateCorrelationInTwoQubitSystem(iqs,0.,0.,pi/2,0.)
        corr1_aba=evaluateCorrelationInTwoQubitSystem(iqs,0.,0.,pi/2-theta,0.)    
        corr1_bba=evaluateCorrelationInTwoQubitSystem(iqs,pi/2,0.,pi/2-theta,0.)
        LHS=abs.(corr1_ab.-corr1_aba)
        RHS=0.25.+corr1_bba
        @test LHS ≈ 0.19583172740687077 atol=1e-10
        @test RHS ≈ 0.09459750793233385 atol=1e-10
    end

    @testset "Quantum teleportation" verbose = true begin
        # Alice secret state, initialize |psi>=u3|0>, inverse initialize qubit 2: Bob |psi> --> u3|psi>=|0>
        # use unitary property: u3*u3=I so initialize and inverse initialize with same u3
        
        theta=rand()
        phi=rand()            
        
        qc=createQuantumCircuit(3)
        # Alice creates her secret state in qubit 0: |psi>=u3|0>
        u3Gate!(qc, 1, theta, phi, pi - phi)
        # create entangled pair for Alice and Bob
        hGate!(qc, 2)
        cnotGate!(qc, 2, 3)
        # Alice performs some steps
        cnotGate!(qc, 1, 2)
        hGate!(qc, 1)
        # Alice performs measurements and transmits results to Bob who then performs some steps
        measureGate!(qc, [1,2])
        cnotGate!(qc, 2, 3)
        controlledUGate!(qc, [1], [3], createSingleQubitOperationZ())
        # after inverse initialize on Bob's qubit with u3|psi>=|0>, measure 3rd qubit which should be 100% |0>
        u3Gate!(qc, 3, theta, phi, pi - phi)
        measureGate!(qc, [3])
        
        qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)        
        iqs=createInitialQubitState(vector, [1. 0.;1. 0.;1. 0.])
        qo=runQuantumProgram(qp, iqs, 1)
        
        @test qo.output[5].measured.probability[1] ≈ 1.0 atol=0.0001            
    end

    @testset "Grovers search algorithm" verbose = true begin
        numberOfQubits = 4
        oracle = 9 # in between 0 and (2^nr_qubits)-1
        
        qcInitializeSuperposition=createQuantumCircuit(numberOfQubits)
        for k in 1:numberOfQubits
            hGate!(qcInitializeSuperposition, k)
        end
        
        qcOracle=createQuantumCircuit(numberOfQubits)
        binaryOracle = reverse(digits(oracle, base=2, pad=numberOfQubits)')
        for k in 1:numberOfQubits
            if binaryOracle[k]==0
                xGate!(qcOracle, k)
            else
                idGate!(qcOracle, k)
            end
        end
        controlledUGate!(qcOracle, Vector(1:(numberOfQubits-1)), [numberOfQubits], createSingleQubitOperationZ())
        for k in 1:numberOfQubits
            if binaryOracle[k]==0
                xGate!(qcOracle, k)
            else
                idGate!(qcOracle, k)
            end
        end
        
        qcAmplification=createQuantumCircuit(numberOfQubits)
        for k in 1:numberOfQubits
            hGate!(qcAmplification, k)
            xGate!(qcAmplification, k)
        end
        controlledUGate!(qcAmplification, Vector(1:(numberOfQubits-1)), [numberOfQubits], createSingleQubitOperationZ())
        for k in 1:numberOfQubits
            xGate!(qcAmplification, k)
            hGate!(qcAmplification, k)
        end
        
        qcMeasure=createQuantumCircuit(numberOfQubits)
        measureGate!(qcMeasure, Vector(1:numberOfQubits))
        
        qcGrover=concatenateQuantumCircuits(qcInitializeSuperposition,
                  concatenateQuantumCircuits(qcOracle,
                  concatenateQuantumCircuits(qcAmplification,
                  concatenateQuantumCircuits(qcOracle,
                  concatenateQuantumCircuits(qcAmplification,qcMeasure)))))
        
        qpGrover=compileQuantumCircuit(qcGrover; optimizeNumberOfSteps=true)
        iqs=createInitialQubitState(vector, repeat([1. 0.], outer=[numberOfQubits,1])) # |high=control,...,low=controlled>=|00...0>
        qo=runQuantumProgram(qpGrover, iqs, 1)
                
        @test qo.output[3].measured.probability ≈ [0.006103515625, 0.006103515625, 0.006103515625, 0.006103515625, 0.006103515625, 0.006103515625, 0.006103515625, 0.006103515625, 0.006103515625, 0.908447265625, 0.006103515625, 0.006103515625, 0.006103515625, 0.006103515625, 0.006103515625, 0.006103515625] atol=0.00000001
    end

    @testset "Phase kickback" verbose = true begin
        function generalizedPhaseKickbackTestOnControlledU(U)
            V=eigvecs(U.U)
            D=eigvals(U.U)
            ind=sortperm(real.(D),rev=true)
            D=D[ind]
            V=V[:,ind]
            
            qc=createQuantumCircuit(2)
            controlledUGate!(qc, [1], [2], U)
            qp=compileQuantumCircuit(qc)
            
            iqs=createInitialQubitState(vector, [1. 0.; V[1,1] V[2,1]])
            qo=runQuantumProgram(qp, iqs, 1)
            oqs=createInitialQubitState(vector, [1. 0.; V[1,1] V[2,1]])
            @test qo.output[2].q ≈ oqs.q atol=0.0001
            
            iqs=createInitialQubitState(vector, [0. 1.; V[1,1] V[2,1]])
            qo=runQuantumProgram(qp, iqs, 1)
            oqs=createInitialQubitState(vector, [0. D[1]*1.; V[1,1] V[2,1]])
            @test qo.output[2].q ≈ oqs.q atol=0.0001
            
            iqs=createInitialQubitState(vector, [1. 0.; V[1,2] V[2,2]])
            qo=runQuantumProgram(qp, iqs, 1)
            oqs=createInitialQubitState(vector, [1. 0.; V[1,2] V[2,2]])
            @test qo.output[2].q ≈ oqs.q atol=0.0001
            
            iqs=createInitialQubitState(vector, [0. 1.; V[1,2] V[2,2]])
            qo=runQuantumProgram(qp, iqs, 1)
            oqs=createInitialQubitState(vector, [0. D[2]*1.; V[1,2] V[2,2]])
            @test qo.output[2].q ≈ oqs.q atol=0.0001            
        end

        @testset "Output: ControlledX gate" begin
            generalizedPhaseKickbackTestOnControlledU(createSingleQubitOperationX())
        end
        @testset "Output: ControlledY gate" begin
            generalizedPhaseKickbackTestOnControlledU(createSingleQubitOperationY())
        end
        @testset "Output: ControlledZ gate" begin
            generalizedPhaseKickbackTestOnControlledU(createSingleQubitOperationZ())
        end

        @testset "Output: Sandwiched H-controlledX-H equals reversed-controlledX" begin
            qc1=createQuantumCircuit(2)
            cnotGate!(qc1, 1, 2)
            qp1=compileQuantumCircuit(qc1; optimizeNumberOfSteps=true)
            qp1.program[1].U

            qc2=createQuantumCircuit(2)
            hGate!(qc2, 1)
            hGate!(qc2, 2)
            cnotGate!(qc2, 2, 1)
            hGate!(qc2, 1)
            hGate!(qc2, 2)
            qp2=compileQuantumCircuit(qc2; optimizeNumberOfSteps=true)
            qp2.program[1].U

            @test qp1.program[1].U ≈ qp2.program[1].U atol=0.0001
        end
    end

    @testset "Quantum erasure with entangled pair" verbose = true begin
        @testset "Output: Partial trace is equivalent to ignorant measurement of 2nd qubit" begin
            qc=createQuantumCircuit(2)
            hGate!(qc, 1)
            cnotGate!(qc, 1, 2)
            measureGate!(qc, [2], forgetOutcome=true)
            qp=compileQuantumCircuit(qc)
            iqs=createInitialQubitState(density, [[1.00,[1. 0.;1. 0.]],[0.00,[0. 1.;0. 1.]]])
            qo=runQuantumProgram(qp, iqs, 1)
            
            # Density ρAB before ignorant measurement of qubit B=0.5*(|00><00|+|00><11|+|11><00|+|11><11|)
            @test real(qo.output[3].rho) ≈ [0.5 0.0 0.0 0.5; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.5 0.0 0.0 0.5] atol=0.0001
            # Density ρAB after ignorant measurement of qubit B=0.5*(|00><00|+|11><11|)
            @test real(qo.output[4].rho) ≈ [0.5 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.5] atol=0.0001
            # Reduced density ρA (partial trace over qubit B) before ignorant measurement of qubit B=0.5*(|0><0|+|1><1|)
            @test real(partialTrace(qo.output[3].rho, traceIndex=2)) ≈ [0.5 0.0; 0.0 0.5] atol=0.0001
            # Taking the partial trace over qubit B is equivalent to an ignorant measurement of qubit B
        end

        @testset "Output: Measure 2nd qubit along σZ and extract which path information" begin
            qc=createQuantumCircuit(2)
            hGate!(qc, 1)
            cnotGate!(qc, 1, 2)
            measureGate!(qc, [2], [sigmaZ()])
            qp=compileQuantumCircuit(qc)
            iqs=createInitialQubitState(density, [[1.00,[1. 0.;1. 0.]],[0.00,[0. 1.;0. 1.]]])
            qo=runQuantumProgram(qp, iqs, 1)
            
            # Density ρAB before measurement of qubit B along σZ=0.5*(|00><00|+|00><11|+|11><00|+|11><11|)
            @test real(qo.output[3].rho) ≈ [0.5 0.0 0.0 0.5; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.5 0.0 0.0 0.5] atol=0.0001
            # Density ρAB after measurement of qubit B along σZ=|00><00| or |11><11|
            @test (round.(real(qo.output[4].rho),digits=2) == [1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0] || round.(real(qo.output[4].rho),digits=2) == [0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0])
            # Reduced density ρA (partial trace over qubit B) after measurement of qubit B along σZ=|0><0| or |1><1|
            @test (round.(real(partialTrace(qo.output[4].rho, traceIndex=2)),digits=2) == [1.0 0.0; 0.0 0.0] || round.(real(partialTrace(qo.output[4].rho, traceIndex=2)),digits=2) == [0.0 0.0; 0.0 1.0])
        end

        @testset "Output: Measure 2nd qubit along σX and erase which path information" begin
            qc=createQuantumCircuit(2)
            hGate!(qc, 1)
            cnotGate!(qc, 1, 2)
            measureGate!(qc, [2], [sigmaX()])
            qp=compileQuantumCircuit(qc)
            iqs=createInitialQubitState(density, [[1.00,[1. 0.;1. 0.]],[0.00,[0. 1.;0. 1.]]])
            qo=runQuantumProgram(qp, iqs, 1)
            
            # Density ρAB before measurement of qubit B along σX=0.5*(|00><00|+|00><11|+|11><00|+|11><11|)
            @test real(qo.output[3].rho) ≈ [0.5 0.0 0.0 0.5; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.5 0.0 0.0 0.5] atol=0.0001
            # Density ρAB after measurement of qubit B along σX=|++><++| or |--><--|
            @test (round.(real(qo.output[4].rho),digits=2) == [0.25 0.25 0.25 0.25; 0.25 0.25 0.25 0.25; 0.25 0.25 0.25 0.25; 0.25 0.25 0.25 0.25] || round.(real(qo.output[4].rho),digits=2) == [0.25 -0.25 -0.25 0.25; -0.25 0.25 0.25 -0.25; -0.25 0.25 0.25 -0.25; 0.25 -0.25 -0.25 0.25])
            # Reduced density ρA (partial trace over qubit B) after measurement of qubit B along σX=|+><+| or |-><-|
            @test (round.(real(partialTrace(qo.output[4].rho, traceIndex=2)),digits=2) == [0.5 0.5; 0.5 0.5] || round.(real(partialTrace(qo.output[4].rho, traceIndex=2)),digits=2) == [0.5 -0.5; -0.5 0.5])
        end
    end

    @testset "Fidelity" verbose = true begin
        @testset "Function: perfect fidelity" begin
            @test fidelity([1 0;0 0],[1 0;0 0]) ≈ 1.0 atol=1e-10
        end
        @testset "Function: zero fidelity" begin
            @test fidelity([1 0;0 0],[0 0;0 1]) ≈ 0.0 atol=1e-10
        end
    end

    @testset "Vector state to density state conversion" verbose = true begin
        vectorState=VectorState(reshape([1.; 1.]/sqrt(2).+0im, 2, 1))
        @test vectorState.q ≈ reshape([1.; 1.]/sqrt(2).+0im, 2, 1) atol=1e-10
        densityState=convertVectorStateToDensityState(vectorState)
        @test densityState.rho ≈ [0.5 0.5; 0.5 0.5].+0im atol=1e-10
        # @test densityState.entropy ≈ 0. atol=1e-10
        # @test densityState.purityValue ≈ 1. atol=1e-10
        # @test densityState.purityNorm ≈ 1. atol=1e-10
    end

    @testset "Depolarizing quantum channel" verbose = true begin
        @testset "Output: Acting on single qubit in superposition" begin
            p=0.75
            qc = createQuantumCircuit(1)
            quantumChannelGate!(qc, [1], generateKrausOperatorsForDepolarizingChannel(p))
            qp = compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            theta = pi/2
            phi = 0
            iqs = createInitialQubitState(density, [[1., [theta phi]]], blochRepresentation=true)
            qo1 = runQuantumProgram(qp, iqs, 1)

            @test qo1.output[2].rho ≈ [(1-4*p/3)*qo1.output[1].rho[1,1]+2*p/3 (1-4*p/3)*qo1.output[1].rho[1,2]; (1-4*p/3)*qo1.output[1].rho[2,1] (1-4*p/3)*qo1.output[1].rho[2,2]+2*p/3] atol=1e-6
        end

        @testset "Output: Acting on one of two entangled qubits" begin
            p=0.75
            # p=0.35
            qc = createQuantumCircuit(2)
            hGate!(qc, 1)
            cnotGate!(qc, 1, 2)
            quantumChannelGate!(qc, [1], generateKrausOperatorsForDepolarizingChannel(p))
            qp = compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            iqs = createInitialQubitState(density, [[1.00, [1. 0.; 1. 0.]]])
            qo2 = runQuantumProgram(qp, iqs, 1)

            # @test qo2.output[1].rho ≈ [1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im; 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im; 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im; 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im] atol=1e-10
            # @test qo2.output[1].entropy ≈ 0. atol=1e-10
            # @test qo2.output[1].purityValue ≈ 1.0 atol=1e-10
            # @test qo2.output[1].purityNorm ≈ 1.0 atol=1e-10

            # @test qo2.output[2].rho ≈ [0.5+0.0im  0.0+0.0im  0.0+0.0im  0.5+0.0im; 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im; 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im; 0.5+0.0im  0.0+0.0im  0.0+0.0im  0.5+0.0im] atol=1e-10
            # @test qo2.output[2].entropy ≈ 0. atol=1e-10
            # @test qo2.output[2].purityValue ≈ 1.0 atol=1e-10
            # @test qo2.output[2].purityNorm ≈ 1.0 atol=1e-10

            @test qo2.output[3].rho ≈ [0.25+0.0im   0.0+0.0im   0.0+0.0im   0.0+0.0im; 0.0+0.0im  0.25+0.0im   0.0+0.0im   0.0+0.0im; 0.0+0.0im   0.0+0.0im  0.25+0.0im   0.0+0.0im; 0.0+0.0im   0.0+0.0im   0.0+0.0im  0.25+0.0im] atol=1e-10
            # @test qo2.output[3].entropy ≈ 1.3862943611198906 atol=1e-10
            # @test qo2.output[3].purityValue ≈ 0.25 atol=1e-10
            # @test qo2.output[3].purityNorm ≈ 0. atol=1e-10
        end
    end

    @testset "Phase damping quantum channel" verbose = true begin
        @testset "Output: Acting on single qubit in superposition" begin
            p=0.75
            qc = createQuantumCircuit(1)
            quantumChannelGate!(qc, [1], generateKrausOperatorsForPhaseDampingChannel(p))
            qp = compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            theta = pi/4
            phi = 0
            iqs = createInitialQubitState(density, [[1., [theta phi]]], blochRepresentation=true)
            qo1 = runQuantumProgram(qp, iqs, 1)

            @test qo1.output[2].rho ≈ [qo1.output[1].rho[1,1] (1-p)*qo1.output[1].rho[1,2]; (1-p)*qo1.output[1].rho[2,1] qo1.output[1].rho[2,2]] atol=1e-6
        end
    end

    @testset "Amplitude damping quantum channel" verbose = true begin
        @testset "Output: Acting on single qubit in superposition" begin
            p=0.65
            qc = createQuantumCircuit(1)
            quantumChannelGate!(qc, [1], generateKrausOperatorsForAmplitudeDampingChannel(p))
            qp = compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            theta = pi/2
            phi = 0
            iqs = createInitialQubitState(density, [[1., [theta phi]]], blochRepresentation=true)
            qo1 = runQuantumProgram(qp, iqs, 1)

            @test qo1.output[2].rho ≈ [qo1.output[1].rho[1,1]+p*qo1.output[1].rho[2,2] sqrt(1-p)*qo1.output[1].rho[1,2]; sqrt(1-p)*qo1.output[1].rho[2,1] (1-p)*qo1.output[1].rho[2,2]] atol=1e-6
        end

        @testset "Output: Acting on single qubit in maximally mixed state" begin
            p=0.45
            qc = createQuantumCircuit(1)
            quantumChannelGate!(qc, [1], generateKrausOperatorsForAmplitudeDampingChannel(p))
            qp = compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
            iqs = createInitialQubitState(density, [[0.5, [1.0 0.0]],[0.5, [0.0 1.0]]])
            qo1 = runQuantumProgram(qp, iqs, 1)

            @test qo1.output[2].rho ≈ [qo1.output[1].rho[1,1]+p*qo1.output[1].rho[2,2] sqrt(1-p)*qo1.output[1].rho[1,2]; sqrt(1-p)*qo1.output[1].rho[2,1] (1-p)*qo1.output[1].rho[2,2]] atol=1e-6
        end
    end
end;
