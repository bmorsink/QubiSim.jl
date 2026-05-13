using QubiSim
using LinearAlgebra, PlotlyJS

function extractNullSpaceOverGF2(Y::Matrix{Int})
    # Solve Y * s = 0 over GF(2)
    # Y is a matrix with rows y_i (0/1 integers)
    Y = copy(Y) .% 2
    m, n = size(Y)

    row = 1
    pivots = []

    # Forward elimination
    for col in 1:n
        # Find pivot
        pivotRow = findfirst(r -> Y[r, col] == 1, row:m)
        if pivotRow === nothing
            continue
        end

        pivotRow += row - 1

        # Swap rows
        Y[row, :], Y[pivotRow, :] = Y[pivotRow, :], Y[row, :]

        push!(pivots, col)

        # Eliminate below
        for r in row+1:m
            if Y[r, col] == 1
                Y[r, :] .= xor.(Y[r, :], Y[row, :])
            end
        end

        row += 1
        if row > m
            break
        end
    end

    # Back substitution to find nontrivial nullspace vector
    s = zeros(Int, n)

    # Free variables: pick last one = 1 (typical Simon case gives 1D nullspace)
    freeVars = setdiff(1:n, pivots)
    if isempty(freeVars)
        error("Only trivial solution found")
    end

    s[freeVars[1]] = 1

    # Back-substitute
    for i in reverse(1:length(pivots))
        col = pivots[i]
        val = 0
        for j in col+1:n
            val ⊻= (Y[i, j] & s[j])
        end
        s[col] = val
    end

    return s
end

print("1. Simons algorithm: arbitrary number of bits\n")
print("\n")

function executeSimonsAlgorithm(s)
    n = length(s)

    qc=createQuantumCircuit(2*n)

    for i in 1:n
        hGate!(qc, i)
    end

    # Step 1: copy input → output
    for i in 1:n
        cnotGate!(qc, i, n+i)
    end

    # Step 2: enforce hidden string s
    onesIndices = findall(x -> x == 1, s)
    if length(onesIndices) > 1
        pivot = onesIndices[1]
        for i in onesIndices[1:end]
            cnotGate!(qc, pivot, n+i)
        end
    end
    # Pick an index j such that: s[j] == 1
    # Then for every i (including j) with s[i] == 1, add: cnotGate!(qc, j, n+i)
    # We are enforcing: output[i] = x[i] ⊕ x[j] (for i where s[i]=1, i≠j)
    # So flipping by s: flips both x[i] and x[j] and cancels out → output unchanged
    # Hence: f(x)=f(x⊕s)

    for i in 1:n
        hGate!(qc, i)
    end

    measureGate!(qc, Vector(1:n))

    numberOfShots = n+2
    qp=compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
    iqs=createInitialQubitState(vector, repeat([1. 0.], outer=[2*n, 1]))
    qo=runQuantumProgram(qp, iqs, 1)
    probeMeasureProbability(qo, 3, "Simon's oracle = " * string(s))
    qo=runQuantumProgram(qp, iqs, numberOfShots)

    Y = reduce(vcat, [reverse(digits(qo.output[3,m].measured.outcome, base=2, pad=n))' for m in 1:numberOfShots])

    println("True s:   ", s)
    println("Found s:  ", extractNullSpaceOverGF2(Y))
    println("")
end

executeSimonsAlgorithm([1,1])
executeSimonsAlgorithm([1,0,1])
executeSimonsAlgorithm([1,1,0,1])

# # Example usage
# n = 5
# # Suppose hidden s is:
# s_true = [1, 0, 1, 1, 0]

# # Generate random y such that y⋅s = 0 mod 2
# function random_y(s)
#     n = length(s)
#     while true
#         y = rand(0:1, n)
#         if sum(y .* s) % 2 == 0
#             return y
#         end
#     end
# end

# # Collect equations
# Y = reduce(vcat, [random_y(s_true)' for _ in 1:n+2])

# # Solve
# s_found = extractNullSpaceOverGF2(Y)

# println("True s:   ", s_true)
# println("Found s:  ", s_found)
