import Base: \

function \(F::LBL{T}, b::AbstractVector{T}) where T
    L, B, vec_P = F.L, F.B, F.vec_P
    n = size(vec_P, 1)

    # Permuter b : b_perm = P * b
    b_perm = b[vec_P]

    # Résoudre L * y = b_perm <=> y = L \ b_perm
    y = similar(b)
    for i in 1:n
        y[i] = b_perm[i]
        for j in 1:i-1
            y[i] -= L[i, j]*y[j]
        end
    end

    # Résoudre B * z = y <=> z = B \ y
    z = similar(b)
    i = 1
    while i < n-1
        if B[i+1, i] == 0
            z[i] = y[i]/B[i, i]

            i += 1
        else
            B_11, B_22, B_21 = B[i, i], B[i+1, i+1], B[i+1, i]
            B_12 = conj(B_21)
            det = (B_11*B_22 - B_21*B_12)
            inv_det = conj(det)/abs(det)^2
            z[i] = inv_det*(B_22*y[i] - B_12*y[i+1])
            z[i+1] = inv_det*(B_11*y[i+1] - B_21*y[i])

            i += 2
        end
    end
    if B[n, n-1] == 0
        z[n-1] = y[n-1]/B[n-1, n-1]
        z[n] = y[n]/B[n, n]
    else
        B_11, B_22, B_21 = B[n-1, n-1], B[n, n], B[n, n-1]
        B_12 = conj(B_21)
        det = (B_11*B_22 - B_21*B_12)
        inv_det = conj(det)/abs(det)^2
        z[n-1] = inv_det*(B_22*y[n-1] - B_12*y[n])
        z[n] = inv_det*(B_11*y[n] - B_21*y[n-1])
    end

    # Résoudre L'* x_perm = z <=> x_perm = L' \ z
    x_perm = similar(b)
    for i in n:-1:1
        x_perm[i] = z[i]
        for j in n:-1:i+1
            x_perm[i] -= conj(L[j, i])*x_perm[j]
        end
    end

    # Permuter x_perm : x = P' * x_perm
    x = similar(b)
    x[vec_P] = x_perm  # Inverse de P est P' car P est une permutation

    return x
end