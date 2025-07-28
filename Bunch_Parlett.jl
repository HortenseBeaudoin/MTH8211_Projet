using Pkg
Pkg.activate("Projet_env")
Pkg.add("BenchmarkTools")
Pkg.add("SuiteSparseMatrixCollection")
Pkg.add("MatrixMarket")

using LinearAlgebra, SparseArrays
using MatrixMarket
using SuiteSparseMatrixCollection
using Test
using BenchmarkTools

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

  α = (1 + sqrt(17))/8
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
      perm_r1_et_r2!(A_, 1, r)

      # Permutation dans L
      for j in 1:k-1
        L[k, j], L[k+r-1, j] = L[k+r-1, j], L[k, j]
      end
                
      # Permutation dans P
      P[k], P[k+r-1] = P[k+r-1], P[k]

      # Déterminant de E
      a_11_conj = conj(A_[1, 1])
      inv_det_E = a_11_conj/μ₁^2

      # Complément de Schur
      for i in 2:n_
        a_i1 = A_[i, 1]
        for j in 2:i
            a_j1_conj = conj(A_[j, 1])
            A_[i, j] -= inv_det_E * a_i1 * a_j1_conj
        end
      end

      # Calcul de B
      diagonal[k] = A_[1, 1]
      subdiagonal[k] = 0

      # Calcul de L
      for i in 2:n_
        L[k+i-1, k] = inv_det_E*A_[i, 1]
      end
      
      # Prochain bloc
      k += 1

    else    # E de taille 2, P permute 1 et p et 2 et q
              
      # Permutation dans A_
      perm_r1_et_r2!(A_, 1, p)
      perm_r1_et_r2!(A_, 2, q)

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

  α = (1 + sqrt(17))/8
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

    if μ₁ ≥ α * μ₀ || k == n # E de taille 1, permute les lignes 1 et r
      # Permutation dans A_
      perm_r1_et_r2!(A_, 1, r)

      # Permutation dans L
      for j in 1:k-1
        A[k, j], A[k+r-1, j] = A[k+r-1, j], A[k, j]
      end
      
      # Permutation dans P
      P[k], P[k+r-1] = P[k+r-1], P[k]

      # Déterminant de E
      a_11_conj = conj(A_[1, 1])
      inv_det_E = a_11_conj/μ₁^2

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
      k += 1

    else

      push!(B2x2, k)
                
      # Permutation dans A_
      perm_r1_et_r2!(A_, 1, p)
      perm_r1_et_r2!(A_, 2, q)

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
      k += 2
    end

  end

    return P, B2x2
end


struct LBL{T}
    """
    Structure pour une matrice factorisée selon PAP' = LBL' (Bunch-Parlett et Bunch-Kaufman)
    """
    L::LowerTriangular{T, Matrix{T}}
    B::Hermitian{T, Tridiagonal{T, Vector{T}}}
    vec_P::Vector{Int}
end

function extract_L_and_B_from_LBL_inplace(A, vec_2by2)
    """
    Récupère la matrice L et B à partir de la forme compacte de la factorisation LBL'
        Entrée : - Matrice A sous la forme compacte de la factorisation LBL'
                 - Vecteur de position des blocs diagonaux 2x2 vec_2by2
        Sortie : - Matrice L triangulaire inférieure
                 - Matrice B tridiagonale hermitienne
    """
    n = size(A, 1)
    T = eltype(A)
    subdiagonal = zeros(T, n-1)
    diagonal = zeros(T, n)
    L = LowerTriangular(Matrix{T}(I, n, n))

    i_diagonal = 1
    while i_diagonal < n - 1
        if i_diagonal in vec_2by2
            # Matrice B
            subdiagonal[i_diagonal], subdiagonal[i_diagonal+1] = A[i_diagonal+1, i_diagonal], 0
            diagonal[i_diagonal], diagonal[i_diagonal + 1] = A[i_diagonal, i_diagonal], A[i_diagonal+1, i_diagonal+1]

            # Matrice L
            L[i_diagonal + 2:n, i_diagonal] .= A[i_diagonal + 2:n, i_diagonal]
            L[i_diagonal + 2:n, i_diagonal + 1] .= A[i_diagonal + 2:n, i_diagonal + 1]

            # Prochain bloc
            i_diagonal += 2
        else
            # Matrice B
            subdiagonal[i_diagonal] = 0
            diagonal[i_diagonal] = A[i_diagonal, i_diagonal]

            # Matrice L
            L[i_diagonal + 1:n, i_diagonal] .= A[i_diagonal + 1:n, i_diagonal]

            # Prochain bloc
            i_diagonal += 1
        end
    end

    if i_diagonal == n - 1
        subdiagonal[i_diagonal] = A[i_diagonal+1, i_diagonal]
        diagonal[i_diagonal], diagonal[i_diagonal+1]  = A[i_diagonal, i_diagonal], A[i_diagonal+1, i_diagonal+1]
    else
        diagonal[i_diagonal] = A[i_diagonal, i_diagonal]
    end

    B = Hermitian(Tridiagonal(subdiagonal, diagonal, conj(subdiagonal)))

    return L, B
end

function extract_P_from_vec_P(vec_P)
    """
    Retourne la matrice de permutation P correspondant à un vecteur de permutation vec_P
        Entrée : - Vecteur de permutation vec_P
        Sortie : - Matrice P
    """
    P = Matrix{Int}(I, n, n)
    P = P[vec_P, :]
    return P
end

function lmult_vec_P(M, vec_P)
    """
    Permute les lignes de la matrice M selon vec_P
        Entrée : - Matrice M
                 - Vecteur de permutation vec_P
        Sortie : - Matrice M permutée
    """
    return M[vec_P, :]
end

function rmult_vec_Pt(M, vec_P)
    """
    Permute les colonnes de la matrice M selon un vecteur de permutation vec_P
        Entrée : - Matrice M
                 - Vecteur de permutation vec_P
        Sortie : - Matrice M permutée
    """
    return M[:, vec_P]
end

function transpose_vec_P(vec_P)
    """
    Transpose un vecteur de permutation vec_P
        Entrée : - Vecteur de permutation vec_P
        Sortie : - Vecteur de permutation transposé vec_P_t
    """
    n = size(vec_P, 1)
    vec_P_t = zeros(Int, n)
    for i in 1:n
        vec_P_t[vec_P[i]] = i
    end
    return vec_P_t
end

function perm_r1_et_r2!(M, r1::Int, r2::Int)
    """
    Modifie une matrice M hermitienne en permutant les lignes et les colonnes d'indices r1 et r2, où r2 >= r1
        Entrée : - Matrice M hermitienne
                 - Indice r1
                 - Indice r2
    """
    if r2 < r1
        error("r2 doit être plus grand ou égal à r1")
    end

    n = size(M, 1)

    M[r1, r1], M[r2, r2] = M[r2, r2], M[r1, r1]
    M[r2, r1] = conj(M[r2, r1])
    for k in r1+1:r2-1
        M[k, r1], M[r2, k] = conj(M[r2, k]), conj(M[k, r1])
    end
    for i in r2+1:n
        M[i, r1], M[i, r2] = M[i, r2], M[i, r1]
    end
    for j in 1:r1-1
        M[r1, j], M[r2, j] = M[r2, j], M[r1, j]
    end
end


# Fichier de tests pour Bunch-Parlett

N = 25
m, n = 50, 25
for i in 1:N
    A = randn(m, n)
    T = eltype(A)
    K_0 = LowerTriangular([sparse(Matrix{T}(I, m, m)) A; A' spzeros(n,n)])
    K = copy(K_0)

    result = bunch_parlett(K)
    L, B, vec_P = result.L, result.B, result.vec_P
    vec_Pt = transpose_vec_P(vec_P)
    K = L*B*L'
    K = lmult_vec_P(K, vec_Pt)
    K = rmult_vec_Pt(K, vec_Pt)
    K = LowerTriangular(K)
    @test norm(K - K_0) < 1e-10
    @test eltype(L) == T && eltype(B) == T
end

for i in 1:N
    A = randn(ComplexF64, m, n)
    T = eltype(A)
    K_0 = LowerTriangular([sparse(Matrix{T}(I, m, m)) A; A' spzeros(n,n)])
    K = copy(K_0)

    result = bunch_parlett(K)
    L, B, vec_P = result.L, result.B, result.vec_P
    vec_Pt = transpose_vec_P(vec_P)
    K = L*B*L'
    K = lmult_vec_P(K, vec_Pt)
    K = rmult_vec_Pt(K, vec_Pt)
    K = LowerTriangular(K)
    @test norm(K - K_0) < 1e-10
    @test eltype(L) == T && eltype(B) == T
end

for i in 1:N
    A = randn(BigFloat, m, n) + im * randn(BigFloat, m, n)
    T = eltype(A)
    K_0 = LowerTriangular([sparse(Matrix{T}(I, m, m)) A; A' spzeros(n,n)])
    K = copy(K_0)

    result = bunch_parlett(K)
    L, B, vec_P = result.L, result.B, result.vec_P
    vec_Pt = transpose_vec_P(vec_P)
    K = L*B*L'
    K = lmult_vec_P(K, vec_Pt)
    K = rmult_vec_Pt(K, vec_Pt)
    K = LowerTriangular(K)
    @test norm(K - K_0) < 1e-10
    @test eltype(L) == T && eltype(B) == T
end

N = 25
m, n = 50, 25
for i in 1:N
    A = randn(m, n)
    T = eltype(A)
    K_0 = LowerTriangular([sparse(Matrix{T}(I, m, m)) A; A' spzeros(n,n)])
    K = copy(K_0)

    vec_P, vec_2by2 = bunch_parlett!(K)
    L, B = extract_L_and_B_from_LBL_inplace(K, vec_2by2)
    vec_Pt = transpose_vec_P(vec_P)
    K = L*B*L'
    K = lmult_vec_P(K, vec_Pt)
    K = rmult_vec_Pt(K, vec_Pt)
    K = LowerTriangular(K)
    @test norm(K - K_0) < 1e-10
    @test eltype(K) == T && eltype(L) == T && eltype(B) == T
end

for i in 1:N
    A = randn(ComplexF64, m, n)
    T = eltype(A)
    K_0 = LowerTriangular([sparse(Matrix{T}(I, m, m)) A; A' spzeros(n,n)])
    K = copy(K_0)

    vec_P, vec_2by2 = bunch_parlett!(K)
    L, B = extract_L_and_B_from_LBL_inplace(K, vec_2by2)
    vec_Pt = transpose_vec_P(vec_P)
    K = L*B*L'
    K = lmult_vec_P(K, vec_Pt)
    K = rmult_vec_Pt(K, vec_Pt)
    K = LowerTriangular(K)
    @test norm(K - K_0) < 1e-10
    @test eltype(K) == T && eltype(L) == T && eltype(B) == T
end

for i in 1:N
    A = randn(BigFloat, m, n) + im * randn(BigFloat, m, n)
    T = eltype(A)
    K_0 = LowerTriangular([sparse(Matrix{T}(I, m, m)) A; A' spzeros(n,n)])
    K = copy(K_0)

    vec_P, vec_2by2 = bunch_parlett!(K)
    L, B = extract_L_and_B_from_LBL_inplace(K, vec_2by2)
    vec_Pt = transpose_vec_P(vec_P)
    K = L*B*L'
    K = lmult_vec_P(K, vec_Pt)
    K = rmult_vec_Pt(K, vec_Pt)
    K = LowerTriangular(K)
    @test norm(K - K_0) < 1e-10
    @test eltype(K) == T && eltype(L) == T && eltype(B) == T
end