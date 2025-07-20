using LinearAlgebra

function lbl_complete(A; tol=1e-14)
    n = size(A, 1)
    A = copy(A)  # Work on a copy
    P = collect(1:n)
    L = Matrix{Float64}(I, n, n)
    B = zeros(n, n)
    k = 1
    α = (1 + sqrt(17)) / 8  

    while k ≤ n
        # Identifie μ₀ et μ₁
        μ₁, r = findmax(abs.(diag(A)[k:end]))   # Extrait la diagonale et trouve le maximum
        r += k - 1

        # Loop sur tous les éléments hors diagonale et met à jour quand une valeur plus élevée est trouvée
        p, q = k, k
        μ₀ = zero(eltype(A))
        for i in k:n
            for j in i:n
                if i != j && abs(A[i, j]) > μ₀
                    μ₀ = abs(A[i, j])
                    p, q = i, j
                end
            end
        end

        # Choix du pivot selon le critère et mise à jour des matrices
        if μ₁ ≥ α * μ₀ || k == n
            # Pivot de taille 1×1

            # Permute les lignes r et k, si r≠k
            if r != k
                A[:, [k, r]] = A[:, [r, k]]
                A[[k, r], :] = A[[r, k], :]
                P[k], P[r] = P[r], P[k]         # Échange les valeurs dans le vecteur de permutation
            end

            # Met B à jour (directement, car E est un scalaire)
            B[k, k] = A[k, k]

            # Met L et A à jour
            if abs(B[k, k]) > tol && k < n
              # Met à jour les valeurs sous la diagonale pour L
                for i in k+1:n
                    L[i, k] = A[i, k] / B[k, k]
                    # Met à jour le bloc inférieur droit pour A
                    for j in k+1:i
                        A[i, j] -= L[i, k] * B[k, k] * L[j, k]
                        A[j, i] = A[i, j] 
                    end
                end
            end
            k += 1      # Mise à jour de k et retour à la while loop
        else

            # Pivot de taille 2×2
            
            # Permute les lignes p et k, si p≠k
            if p != k
                A[:, [k, p]] = A[:, [p, k]]
                A[[k, p], :] = A[[p, k], :]
                P[k], P[p] = P[p], P[k]
            end
            # Permute les lignes q et k+1, si q≠k+1
            if q != k+1
                A[:, [k+1, q]] = A[:, [q, k+1]]
                A[[k+1, q], :] = A[[q, k+1], :]
                P[k+1], P[q] = P[q], P[k+1]
            end

            # Extrait le bloc E (car pas un scalaire)
            E = A[k:k+1, k:k+1]
            # Met B à jourS
            B[k:k+1, k:k+1] = E

            # Met L et A à jour
            if abs(det(E)) > tol && k+1 < n
                E_inv = inv(E)      # Calcule une seule fois l'inverse de E
                # Met à jour les valeurs sous la diagonale pour L
                for i in k+2:n
                    L[i, k] = dot(A[i, k:k+1], E_inv[:, 1])
                    L[i, k+1] = dot(A[i, k:k+1], E_inv[:, 2])
                end
                # Met à jour le bloc inférieur droit pour A
                for i in k+2:n
                    for j in k+2:i
                        A[i, j] -= dot(L[i, k:k+1], E * L[j, k:k+1])
                        A[j, i] = A[i, j]
                    end
                end
            end
            k += 2        # Mise à jour de k et retour à la while loop
        end
    end

    return P, L, B
end



A = [0.0 1.0 2.0; 1.0 0.0 3.0; 2.0 3.0 0.0]
p, L, B = lbl_complete(A)

# Verify the factorization
A_perm = A[p, p]
A_fact = L * B * L'
A_all = A_fact[p, p]
residual = norm(A_perm - A_fact)

