function choose_pivot_bunch_parlett(A_, n_, alpha, T)    
    mu1 = -abs(one(T))
    r = 1
    # Diagonal
    for i in 1:n_
        magnitude_aii = abs(A_[i, i])
        if magnitude_aii > mu1
            mu1 = magnitude_aii
            r = i
        end
    end

    mu0 = -abs(one(T))
    p, q = 1, 1
    # Subdiagonal
    for i in 1:n_
        for j in i+1:n_
            magnitude_aij = abs(A_[i, j])
            if magnitude_aij > mu0
                mu0 = magnitude_aij
                p, q = i, j
            end
        end
    end

    if mu1 >= alpha * mu0
        return (1, r), nothing  # E de taille 1, P permute 1 et r
    else
        return (1, p), (2, q)  # E de taille 2, P permute 1 et p, puis 2 et q
    end
end

function choose_pivot_bunch_kaufman(A_, n_, alpha, T)    
    w_1 = -abs(one(T))
    r = 1
    # First column
    for i in 2:n_
        magnitude_a_i1 = abs(A_[i, 1])
        if magnitude_a_i1 > w_1
            w_1 = magnitude_a_i1
            r = i
        end
    end

    magnitude_a_11 = abs(A_[1, 1])
    if magnitude_a_11 >= alpha * w_1
        return nothing, nothing  # E de taille 1, P ne permute rien
    else
        w_r = -abs(one(T))
        # Subdiagonal
        for i in r+1:n_
            magnitude_a_ir = abs(A_[i, r])
            if magnitude_a_ir > w_r
                w_r = magnitude_a_ir
                r_new = i
            end
        end
        # Superdiagonal
        for j in 1:r-1
            magnitude_a_ir = abs(A_[r, j])
            if magnitude_a_ir > w_r
                w_r = magnitude_a_ir
            end
        end

        if magnitude_a_11 * w_r >= alpha * w_1^2
            return nothing, nothing  # E de taille 1, P ne permute rien
        else
            magnitude_a_rr = abs(A_[r, r])
            if magnitude_a_rr >= alpha * w_r
                return (1, r), nothing  # E de taille 1, P permute 1 et r
            else             
                return (2, r), nothing  # E de taille 2, P permute 2 et r       
            end
        end
    end
end

function choose_pivot_bunch_kaufman_rook(A_, n_, alpha, T)    
    w_1 = -abs(one(T))
    r = 1
    # First column
    for i in 2:n_
        magnitude_a_i1 = abs(A_[i, 1])
        if magnitude_a_i1 > w_1
            w_1 = magnitude_a_i1
            r = i
        end
    end

    magnitude_a_11 = abs(A_[1, 1])
    if magnitude_a_11 >= alpha * w_1
        return nothing, nothing  # E de taille 1, P ne permute rien
    else
        index = 1
        w_index = w_1
        r_new = r
        while true
            w_r = -abs(one(T))
            # Subdiagonal
            for i in r+1:n_
                magnitude_a_ir = abs(A_[i, r])
                if magnitude_a_ir > w_r
                    w_r = magnitude_a_ir
                    r_new = i
                end
            end
            # Superdiagonal
            for j in 1:r-1
                magnitude_a_ir = abs(A_[r, j])
                if magnitude_a_ir > w_r
                    w_r = magnitude_a_ir
                end
            end

            magnitude_a_rr = abs(A_[r, r])
            if magnitude_a_rr >= alpha * w_r  # E de taille 1, P permute 1 et r
                return (1, r), nothing
            elseif w_index == w_r  # E de taille 2, P permute 1 et index, puis 2 et r     
                return (1, index), (2, r)
            else
                index = r
                w_index = w_r
                r = r_new
            end
        end
    end
end

function facto_LBL(A::LowerTriangular, pivoting_strategy::Int=1)
    """
    Factorise A selon la factorisation LBL'
        Entrée : - Matrice A carrée, hermitienne et indéfinie
                 - Stratégie de pivotage souhaitée : 
                    - 1 pour pivotage complet (Bunch-Parlett)
                    - 2 pour pivotage partiel (Bunch-Kaufman)
                    - 3 pour pivotage de la tour (Rook pivoting)
        Sortie : - Structure qui contient :
                    - Matrice L triangulaire inférieure
                    - Matrice B tridiagonale hermitienne
                    - Vecteur de permutation vec_P
    """
    n = size(A, 1)
    T = eltype(A)
    subdiagonal = zeros(T, n-1)
    diagonal = zeros(T, n)
    L = UnitLowerTriangular(Matrix{T}(I, n, n))
    vec_P = collect(1:n)

    alpha = abs(T((1 + sqrt(17))/8))
    i_diagonal = 1
    while i_diagonal < n - 1
        ### Initialisation de la partie de A traitée
        A_ = view(A, i_diagonal:n, i_diagonal:n)
        n_ = size(A_, 1)

        ### Choix du pivot
        if pivoting_strategy === 1
            pivot1, pivot2 = choose_pivot_bunch_parlett(A_, n_, alpha, T)
        elseif pivoting_strategy === 2
            pivot1, pivot2 = choose_pivot_bunch_kaufman(A_, n_, alpha, T)
        elseif pivoting_strategy === 3
            pivot1, pivot2 = choose_pivot_bunch_kaufman_rook(A_, n_, alpha, T)
        end

        ### Pivotage
        if pivot2 !== nothing
            # Permutation dans A_
            r1, r2 = pivot1[1], pivot1[2]
            r3, r4 = pivot2[1], pivot2[2]
            perm_r1_et_r2!(A_, r1, r2)
            perm_r1_et_r2!(A_, r3, r4)

            # Permutation dans L
            for j in 1:i_diagonal-1
                L[i_diagonal+r1-1, j], L[i_diagonal+r2-1, j] = L[i_diagonal+r2-1, j], L[i_diagonal+r1-1, j]
            end
            for j in 1:i_diagonal-1
                L[i_diagonal+r3-1, j], L[i_diagonal+r4-1, j] = L[i_diagonal+r4-1, j], L[i_diagonal+r3-1, j]
            end
            
            # Permutation dans P
            vec_P[i_diagonal+r1-1], vec_P[i_diagonal+r2-1] = vec_P[i_diagonal+r2-1], vec_P[i_diagonal+r1-1]
            vec_P[i_diagonal+r3-1], vec_P[i_diagonal+r4-1] = vec_P[i_diagonal+r4-1], vec_P[i_diagonal+r3-1]
        elseif pivot1 !== nothing
            # Permutation dans A_
            r1, r2 = pivot1[1], pivot1[2]
            perm_r1_et_r2!(A_, r1, r2)

            # Permutation dans L
            for j in 1:i_diagonal-1
                L[i_diagonal+r1-1, j], L[i_diagonal+r2-1, j] = L[i_diagonal+r2-1, j], L[i_diagonal+r1-1, j]
            end
            
            # Permutation dans P
            vec_P[i_diagonal+r1-1], vec_P[i_diagonal+r2-1] = vec_P[i_diagonal+r2-1], vec_P[i_diagonal+r1-1]
        end

        ### Calcul propre à la factorisation
        if pivot2 !== nothing
            # Déterminant de E
            e_11, e_22, e_21 = A_[1, 1], A_[2, 2], A_[2, 1]
            e_12 = conj(e_21)
            det_E = (e_11*e_22 - e_21*e_12)
            det_E_conj = conj(det_E)
            magnitude_det_E = abs(det_E)
            inv_det_E = det_E_conj/magnitude_det_E^2

            # Complément de Schur
            for i in 3:n_
                a_i1, a_i2 = A_[i, 1], A_[i, 2]
                for j in 3:i
                    a_j1_conj, a_j2_conj = conj(A_[j, 1]), conj(A_[j, 2])
                    A_[i, j] -= inv_det_E*((a_i1*e_22 - a_i2*e_21)*a_j1_conj + (a_i2*e_11 - a_i1*e_12)*a_j2_conj)
                end
            end
            
            # Calcul de B
            subdiagonal[i_diagonal], subdiagonal[i_diagonal+1] = A_[2, 1], 0
            diagonal[i_diagonal], diagonal[i_diagonal + 1] = A_[1, 1], A_[2, 2]

            # Calcul de L
            for i in i_diagonal+2:n
                a_i1, a_i2 = A[i, i_diagonal], A[i, i_diagonal+1]
                L[i, i_diagonal] = inv_det_E*(a_i1*e_22 - a_i2*e_21)
                L[i, i_diagonal+1] = inv_det_E*(a_i2*e_11 - a_i1*e_12)
            end

            # Prochain bloc
            i_diagonal += 2
        else
            # Déterminant de E
            magnitude_a_11 = abs(A_[1, 1])
            a_11_conj = conj(A_[1, 1])
            inv_det_E = a_11_conj/magnitude_a_11^2

            # Complément de Schur
            for i in 2:n_
                a_i1 = A_[i, 1]
                for j in 2:i
                    a_j1_conj = conj(A_[j, 1])
                    A_[i, j] -= inv_det_E * a_i1 * a_j1_conj
                end
            end
            
            # Calcul de B
            diagonal[i_diagonal] = A_[1, 1]
            subdiagonal[i_diagonal] = 0

            # Calcul de L
            for i in 2:n_
                L[i_diagonal+i-1, i_diagonal] = inv_det_E*A_[i, 1]
            end

            # Prochain bloc
            i_diagonal += 1
        end
    end

    # Dernier bloc
    if i_diagonal == n - 1
        subdiagonal[i_diagonal] = A[i_diagonal+1, i_diagonal]
        diagonal[i_diagonal], diagonal[i_diagonal+1]  = A[i_diagonal, i_diagonal], A[i_diagonal+1, i_diagonal+1]
    else
        diagonal[i_diagonal] = A[i_diagonal, i_diagonal]
    end

    B = Hermitian(Tridiagonal(subdiagonal, diagonal, conj(subdiagonal)))
    result = LBL(L, B, vec_P)
    return result
end

function facto_LBL!(A::LowerTriangular, pivoting_strategy::Int=1)
    """
    Modifie A sous la forme compacte de la factorisation LBL'
    L et B sont stockés en écrasant A, où L se situe en dessous des éléments blocs diagonaux de B.
        Entrée : - Matrice A carrée, hermitienne et indéfinie
                 - Stratégie de pivotage souhaitée : 
                    - 1 pour pivotage complet (Bunch-Parlett)
                    - 2 pour pivotage partiel (Bunch-Kaufman)
                    - 3 pour pivotage de la tour (Rook pivoting)
        Sortie : - Vecteur de permutation vec_P
                 - Vecteur de position des blocs diagonaux 2x2 vec_2by2
    """
    n = size(A, 1)
    T = eltype(A)
    vec_P = collect(1:n)
    vec_2by2 = Int[]

    alpha = abs(T((1 + sqrt(17))/8))
    i_diagonal = 1
    while i_diagonal < n - 1
        ### Initialisation de la partie de A traitée
        A_ = view(A, i_diagonal:n, i_diagonal:n)
        n_ = size(A_, 1)

        ### Choix du pivot
        if pivoting_strategy === 1
            pivot1, pivot2 = choose_pivot_bunch_parlett(A_, n_, alpha, T)
        elseif pivoting_strategy === 2
            pivot1, pivot2 = choose_pivot_bunch_kaufman(A_, n_, alpha, T)
        elseif pivoting_strategy === 3
            pivot1, pivot2 = choose_pivot_bunch_kaufman_rook(A_, n_, alpha, T)
        end

        ### Pivotage
        if pivot2 !== nothing
            # Permutation dans A_
            r1, r2 = pivot1[1], pivot1[2]
            r3, r4 = pivot2[1], pivot2[2]
            perm_r1_et_r2!(A_, r1, r2)
            perm_r1_et_r2!(A_, r3, r4)

            # Permutation dans L
            for j in 1:i_diagonal-1
                A[i_diagonal+r1-1, j], A[i_diagonal+r2-1, j] = A[i_diagonal+r2-1, j], A[i_diagonal+r1-1, j]
            end
            for j in 1:i_diagonal-1
                A[i_diagonal+r3-1, j], A[i_diagonal+r4-1, j] = A[i_diagonal+r4-1, j], A[i_diagonal+r3-1, j]
            end
            
            # Permutation dans P
            vec_P[i_diagonal+r1-1], vec_P[i_diagonal+r2-1] = vec_P[i_diagonal+r2-1], vec_P[i_diagonal+r1-1]
            vec_P[i_diagonal+r3-1], vec_P[i_diagonal+r4-1] = vec_P[i_diagonal+r4-1], vec_P[i_diagonal+r3-1]
        elseif pivot1 !== nothing
            # Permutation dans A_
            r1, r2 = pivot1[1], pivot1[2]
            perm_r1_et_r2!(A_, r1, r2)

            # Permutation dans L
            for j in 1:i_diagonal-1
                A[i_diagonal+r1-1, j], A[i_diagonal+r2-1, j] = A[i_diagonal+r2-1, j], A[i_diagonal+r1-1, j]
            end
            
            # Permutation dans P
            vec_P[i_diagonal+r1-1], vec_P[i_diagonal+r2-1] = vec_P[i_diagonal+r2-1], vec_P[i_diagonal+r1-1]
        end

        ### Calcul propre à la factorisation
        if pivot2 !== nothing
            # Position du bloc 2 x 2
            push!(vec_2by2, i_diagonal)

            # Déterminant de E
            e_11, e_22, e_21 = A_[1, 1], A_[2, 2], A_[2, 1]
            e_12 = conj(e_21)
            det_E = (e_11*e_22 - e_21*e_12)
            det_E_conj = conj(det_E)
            magnitude_det_E = abs(det_E)
            inv_det_E = det_E_conj/magnitude_det_E^2

            # Complément de Schur
            for i in 3:n_
                a_i1, a_i2 = A_[i, 1], A_[i, 2]
                for j in 3:i
                    a_j1_conj, a_j2_conj = conj(A_[j, 1]), conj(A_[j, 2])
                    A_[i, j] -= inv_det_E*((a_i1*e_22 - a_i2*e_21)*a_j1_conj + (a_i2*e_11 - a_i1*e_12)*a_j2_conj)
                end
            end

            # Calcul de L
            for i in 3:n_
                a_i1, a_i2 = A_[i, 1], A_[i, 2]
                A_[i, 1] = inv_det_E*(a_i1*e_22 - a_i2*e_21)
                A_[i, 2] = inv_det_E*(a_i2*e_11 - a_i1*e_12)
            end

            # Prochain bloc
            i_diagonal += 2
        else
            # Déterminant de E
            magnitude_a_11 = abs(A_[1, 1])
            a_11_conj = conj(A_[1, 1])
            inv_det_E = a_11_conj/magnitude_a_11^2

            # Complément de Schur
            for i in 2:n_
                a_i1 = A_[i, 1]
                for j in 2:i
                    a_j1_conj = conj(A_[j, 1])
                    A_[i, j] -= inv_det_E * a_i1 * a_j1_conj
                end
            end
            
            # Calcul de L
            for i in 2:n_
                A_[i, 1] *= inv_det_E
            end

            # Prochain bloc
            i_diagonal += 1
        end
    end

    return vec_P, vec_2by2
end