function bunch_parlett(A::LowerTriangular)
  """
  Factorise A = LBL' selon la factorisation de Bunch-Parlett, utilisant la technique du pivotage complet.
  Entrée : Matrice A carrée, hermitienne et indéfinie
  Sortie : Matrice L triangulaire inférieure
           Matrice B tridiagonale hermitienne
           Vecteur de permutation P
  """
  n = size(A, 1)
  T = eltype(A)
  subdiagonal = zeros(T, n-1)
  diagonal = zeros(T, n)
  L = LowerTriangular(Matrix{T}(I, n, n))
  P = collect(1:n)

  alpha = (1 + sqrt(17))/8
  i_diagonal = 1
  while i_diagonal < n - 1
    ### Initialisation de la partie de A traitée
    A_ = view(A, i_diagonal:n, i_diagonal:n)
    n_ = size(A_, 1)

    ### Choix du pivot et pivotage
    mu1, r = findmax(abs.(diag(A_)))

    p, q = i_diagonal, i_diagonal
    mu0 = zero(Float64)
    for i in i_diagonal:n
      for j in (i+1):n
        if abs(A[i, j]) > mu0
          μ₀ = abs(A[i, j])
          p, q = i, j
        end
      end
    end

    if mu1 ≥ alpha * mu0 || i_diagonal == n # E de taille 1, permute les lignes r et 1

      # Permutation dans A_
      perm_r1_et_r2!(A_, 1, r)

      # Permutation dans L
      for j in 1:i_diagonal-1
        L[i_diagonal, j], L[i_diagonal+r-1, j] = L[i_diagonal+r-1, j], L[i_diagonal, j]
      end
                
      # Permutation dans P
      P[i_diagonal], P[i_diagonal+r-1] = P[i_diagonal+r-1], P[i_diagonal]

      # Déterminant de E
      a_11_conj = conj(A_[1, 1])
      inv_det_E = a_11_conj/mu1^2

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

    else    # E de taille 2, P permute 1 et p et 2 et q
              
      # Permutation dans A_
      perm_r1_et_r2!(A_, 1, p)
      perm_r1_et_r2!(A_, 2, q)

      # Permutation dans L
      for j in 1:i_diagonal-1
        L[i_diagonal, j], L[i_diagonal+p-1, j] = L[i_diagonal+p-1, j], L[i_diagonal, j]
      end
      for j in 1:i_diagonal-1
        L[i_diagonal+1, j], L[i_diagonal+q-1, j] = L[i_diagonal+q-1, j], L[i_diagonal+1, j]
      end
                  
      # Permutation dans P
      P[i_diagonal], P[i_diagonal+p-1] = P[i_diagonal+p-1], P[i_diagonal]
      P[i_diagonal+1], P[i_diagonal+q-1] = P[i_diagonal+q-1], P[i_diagonal+1]

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
      for i in i_diagonal + 2:n
        a_i1, a_i2 = A[i, i_diagonal], A[i, i_diagonal+1]
        L[i, i_diagonal] = inv_det_E*(a_i1*e_22 - a_i2*e_21)
        L[i, i_diagonal+1] = inv_det_E*(a_i2*e_11 - a_i1*e_12)
      end

      # Prochain bloc
      i_diagonal += 2
    end

  end

  if i_diagonal == n - 1
    subdiagonal[i_diagonal] = A[i_diagonal+1, i_diagonal]
    diagonal[i_diagonal], diagonal[i_diagonal+1]  = A[i_diagonal, i_diagonal], A[i_diagonal+1, i_diagonal+1]
  else
    diagonal[i_diagonal] = A[i_diagonal, i_diagonal]
  end

  B = Hermitian(Tridiagonal(subdiagonal, diagonal, conj(subdiagonal)))
  result = LBL(L, B, P)

  return result
end

function bunch_parlett!(A::LowerTriangular)
  """
  Factorise A = LBL' selon la factorisation de Bunch-Parlett, utilisant la technique du pivotage complet.
  Entrée : Matrice A carrée, hermitienne et indéfinie
  Sortie : Matrice A carrée et hermitienne contenant à la fois L et B :
              - L sur les éléments hors diagonaux
              - B sur les éléments blocs diagonaux ne faisant pas partie de L
            Vecteur de permutation P
            Vecteur de position des blocs 2x2
  """
  n = size(A, 1)
  P = collect(1:n)
  B2x2 = Int[]

  alpha = (1 + sqrt(17))/8
  i_diagonal = 1
  while i_diagonal < n - 1
    ### Initialisation de la partie de A traitée
    A_ = view(A, i_diagonal:n, i_diagonal:n)
    n_ = size(A_, 1)

    ### Choix du pivot et pivotage
    mu1, r = findmax(abs.(diag(A_)))

    p, q = i_diagonal, i_diagonal
    mu0 = zero(Float64)
    for i in i_diagonal:n
      for j in (i+1):n
        if abs(A[i, j]) > mu0
          μ₀ = abs(A[i, j])
          p, q = i, j
        end
      end
    end

    if mu1 ≥ alpha * mu0 || i_diagonal == n # E de taille 1, permute les lignes 1 et r
      # Permutation dans A_
      perm_r1_et_r2!(A_, 1, r)

      # Permutation dans L
      for j in 1:i_diagonal-1
        A[i_diagonal, j], A[i_diagonal+r-1, j] = A[i_diagonal+r-1, j], A[i_diagonal, j]
      end
      
      # Permutation dans P
      P[i_diagonal], P[i_diagonal+r-1] = P[i_diagonal+r-1], P[i_diagonal]

      # Déterminant de E
      a_11_conj = conj(A_[1, 1])
      inv_det_E = a_11_conj/mu1^2

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
        A_[i, 1] *= inv_det_E      # À VALIDER
      end
      
      # Prochain bloc
      i_diagonal += 1

    else

      push!(B2x2, i_diagonal)
                
      # Permutation dans A_
      perm_r1_et_r2!(A_, 1, p)
      perm_r1_et_r2!(A_, 2, q)

      # Permutation dans L
      for j in 1:i_diagonal-1
        A[i_diagonal, j], A[i_diagonal+p-1, j] = A[i_diagonal+p-1, j], A[i_diagonal, j]
      end
      for j in 1:i_diagonal-1
        A[i_diagonal+1, j], A[i_diagonal+q-1, j] = A[i_diagonal+q-1, j], A[i_diagonal+1, j]
      end
      
      # Permutation dans P
      P[i_diagonal], P[i_diagonal+p-1] = P[i_diagonal+p-1], P[i_diagonal]
      P[i_diagonal+1], P[i_diagonal+q-1] = P[i_diagonal+q-1], P[i_diagonal+1]

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
    end

  end

    return P, B2x2
end