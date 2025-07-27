using LinearAlgebra, SparseArrays
using Test

# Implémentation Bunch-Parlett
function perm_A!(M, r₁::Int, r₂::Int)
    """
    Permute les lignes et les colonnes r₁ et r₂ de M hermitienne, où r₂ > r₁
    """
    if r₂ < r₁
      error("Index r₂ must be greater than index r₁")
    end

    n = size(M, 1)

    M[r₁, r₁], M[r₂, r₂] = M[r₂, r₂], M[r₁, r₁]
    M[r₂, r₁] = conj(M[r₂, r₁])
    for k in r₁+1:r₂-1
      M[k, r₁], M[r₂, k] = conj(M[r₂, k]), conj(M[k, r₁])
    end
    for i in r₂+1:n
      M[i, r₁], M[i, r₂] = M[i, r₂], M[i, r₁]
    end
    for j in 1:r₁-1
      M[r₁, j], M[r₂, j] = M[r₂, j], M[r₁, j]
    end
end

function bunch_parlett(A::Hermitian)
  """
  Factorise A = LBL' selon la factorisation de Bunch-Parlett, utilisant la technique du pivotage complet.
  Entrée : Matrice A carrée, hermitienne et indéfinie
  Sortie : Matrice L triangulaire inférieure
           Matrice B tridiagonale hermitienne
           Vecteur de permutation P
  """
  A = LowerTriangular(Matrix(A))

  α = (1 + sqrt(17))/8
  n = size(A, 1)
  T = eltype(A)
  subdiagonal = zeros(T, n-1)
  diagonal = zeros(T, n)
  L = LowerTriangular(Matrix{T}(I, n, n))
  P = collect(1:n)

  k = 1

  while k < n - 1
    ### Initialisation de la partie de A traitée
    A_ = view(A, k:n, k:n)
    n_ = size(A_, 1)

    ### Choix du pivot et pivotage

    μ₁, r = findmax(abs.(diag(A_)))

    p, q = k, k
    μ₀ = zero(Float64)
    for i in k:n
      for j in (i+1):n
        if abs(A[i, j]) > μ₀
          μ₀ = abs(A[i, j])
          p, q = i, j
        end
      end
    end

    if μ₁ ≥ α * μ₀ || k == n # E de taille 1, permute les lignes r et 1

      # Permutation dans A_
      perm_A!(A_, 1, r)

      # Permutation dans L
      for j in 1:k-1
        L[k, j], L[k+r-1, j] = L[k+r-1, j], L[k, j]
      end
                
      # Permutation dans P
      P[k], P[k+r-1] = P[k+r-1], P[k]

      # Déterminant de E
      inv_det_E = conj(A_[1, 1])/μ₁^2

      # Complément de Schur
      for i in 2:n_
        for j in 2:i
          A_[i, j] -= inv_det_E*(A_[i, 1] * conj(A_[j, 1]))
        end
      end

      # Calcul de B
      diagonal[k] = A_[1, 1]
      subdiagonal[k] = 0

      # Calcul de L
      L[k+1:n, k] .= inv_det_E*A_[2:n_, 1]
      
      # Prochain bloc
      k += 1

    else    # E de taille 2, P permute 1 et p et 2 et q
              
      # Permutation dans A_
      perm_A!(A_, 1, p)
      perm_A!(A_, 2, q)

      # Permutation dans L
      for j in 1:k-1
        L[k, j], L[k+p-1, j] = L[k+p-1, j], L[k, j]
      end
      for j in 1:k-1
        L[k+1, j], L[k+q-1, j] = L[k+q-1, j], L[k+1, j]
      end
                  
      # Permutation dans P
      P[k], P[k+p-1] = P[k+p-1], P[k]
      P[k+1], P[k+q-1] = P[k+q-1], P[k+1]

      # Déterminant de E
      e_11, e_22, e_21 = A_[1, 1], A_[2, 2], A_[2, 1]
      e_12 = conj(e_21)
      det_E = (e_11*e_22 - e_21*e_12)
      inv_det_E = conj(det_E)/abs(det_E)^2

      # Complément de Schur
      for i in 3:n_
        for j in 3:i
          A_[i, j] -= inv_det_E*((A_[i, 1]*e_22 - A_[i, 2]*e_21)*conj(A_[j, 1]) + (A_[i, 2]*e_11 - A_[i, 1]*e_12)*conj(A_[j, 2]))
        end
      end
                  
      # Calcul de B
      subdiagonal[k], subdiagonal[k+1] = A_[2, 1], 0
      diagonal[k], diagonal[k + 1] = A_[1, 1], A_[2, 2]

      # Calcul de L
      for i in k + 2:n
        a_i1, a_i2 = A[i, k], A[i, k+1]
        L[i, k] = inv_det_E*(a_i1*e_22 - a_i2*e_21)
        L[i, k+1] = inv_det_E*(a_i2*e_11 - a_i1*e_12)
      end

      # Prochain bloc
      k += 2
    end

  end

  if k == n - 1
    subdiagonal[k] = A[k+1, k]
    diagonal[k], diagonal[k+1]  = A[k, k], A[k+1, k+1]
  else
    diagonal[k] = A[k, k]
  end

  B = Hermitian(Tridiagonal(subdiagonal, diagonal, conj(subdiagonal)))

  return L, B, P
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

  α = (1 + sqrt(17))/8
  n = size(A, 1)
  P = collect(1:n)
  B2x2 = Int[]

  k = 1
  while k < n - 1
    ### Initialisation de la partie de A traitée
    A_ = view(A, k:n, k:n)
    n_ = size(A_, 1)

    ### Choix du pivot et pivotage
    μ₁, r = findmax(abs.(diag(A_)))

    p, q = k, k
    μ₀ = zero(Float64)
    for i in k:n
      for j in (i+1):n
        if abs(A[i, j]) > μ₀
          μ₀ = abs(A[i, j])
          p, q = i, j
        end
      end
    end

    if μ₁ ≥ α * μ₀ || k == n # E de taille 1, permute les lignes r et 1
      # Permutation dans A_
      perm_A!(A_, 1, r)

      # Permutation dans L
      for j in 1:k-1
        A[k, j], A[k+r-1, j] = A[k+r-1, j], A[k, j]
      end
      
      # Permutation dans P
      P[k], P[k+r-1] = P[k+r-1], P[k]

      # Déterminant de E
      inv_det_E = conj(A_[1, 1])/μ₁^2

      # Complément de Schur
      for i in 2:n_
        for j in 2:i
          A_[i, j] -= inv_det_E*(A_[i, 1] * conj(A_[j, 1]))
        end
      end

      # Calcul de L
      A_[2:n_, 1] .*= inv_det_E
      
      # Prochain bloc
      k += 1

    else

      push!(B2x2, k)
                
      # Permutation dans A_
      perm_A!(A_, 1, p)
      perm_A!(A_, 2, q)

      # Permutation dans L
      for j in 1:k-1
        A[k, j], A[k+p-1, j] = A[k+p-1, j], A[k, j]
      end
      for j in 1:k-1
        A[k+1, j], A[k+q-1, j] = A[k+q-1, j], A[k+1, j]
      end
      
      # Permutation dans P
      P[k], P[k+p-1] = P[k+p-1], P[k]
      P[k+1], P[k+q-1] = P[k+q-1], P[k+1]

      # Déterminant de E
      e_11, e_22, e_21 = A_[1, 1], A_[2, 2], A_[2, 1]
      e_12 = conj(e_21)
      det_E = (e_11*e_22 - e_21*e_12)
      inv_det_E = conj(det_E)/abs(det_E)^2

      # Complément de Schur
      for i in 3:n_
        for j in 3:i
          A_[i, j] -= inv_det_E*((A_[i, 1]*e_22 - A_[i, 2]*e_21)*conj(A_[j, 1]) + (A_[i, 2]*e_11 - A_[i, 1]*e_12)*conj(A_[j, 2]))
        end
      end

      # Calcul de L
      for i in 3:n_
        a_i1, a_i2 = A_[i, 1], A_[i, 2]
        A_[i, 1] = inv_det_E*(a_i1*e_22 - a_i2*e_21)
        A_[i, 2] = inv_det_E*(a_i2*e_11 - a_i1*e_12)
      end

      # Prochain bloc
      k += 2
    end

  end

    return P, B2x2
end

function extract_L_and_B(A, B2x2)
  """
  Récupère la matrice L et B à partir de la forme compacte de la factorisation de Bunch-Parlett
  Entrée : Matrice A compacte factorisée selon Bunch-Parlett
  Sortie : Matrice L
            Matrice B
  """
  n = size(A, 1)
  T = eltype(A)

  subdiagonal = zeros(T, n-1)
  diagonal = zeros(T, n)
  L = LowerTriangular(Matrix{T}(I, n, n))

  k = 1
  while k < n - 1
      if k in B2x2
          # Matrice B
          subdiagonal[k], subdiagonal[k+1] = A[k+1, k], 0
          diagonal[k], diagonal[k + 1] = A[k, k], A[k+1, k+1]

          # Matrice L
          L[k + 2:n, k] .= A[k + 2:n, k]
          L[k + 2:n, k + 1] .= A[k + 2:n, k + 1]

          # Prochain bloc
          k += 2
      else
          # Matrice B
          subdiagonal[k] = 0
          diagonal[k] = A[k, k]

          # Matrice L
          L[k + 1:n, k] .= A[k + 1:n, k]

          # Prochain bloc
          k += 1
      end
  end

  if k == n - 1
      subdiagonal[k] = A[k+1, k]
      diagonal[k], diagonal[k+1]  = A[k, k], A[k+1, k+1]
  else
      diagonal[k] = A[k, k]
  end

  B = Hermitian(Tridiagonal(subdiagonal, diagonal, conj(subdiagonal)))

  return L, B
end

