# Fichier de tests pour Bunch-Parlett

N = 25
m, n = 50, 25
for i in 1:N
    A = randn(m, n)
    T = eltype(A)
    K_0 = Hermitian([sparse(Matrix{T}(I, m, m)) A; A' spzeros(n,n)])
    K = copy(K_0)

    L, B, vec_P = bunch_parlett(K)
    vec_P_transposed = transpose_vec_P(vec_P)
    K = L*B*L'
    K = perm_line_vec_P(K, vec_P_transposed)
    K = perm_column_vec_P(K, vec_P_transposed)
    @test norm(K - K_0) < 1e-10
end

for i in 1:N
    A = randn(ComplexF64, m, n)
    T = eltype(A)
    K_0 = Hermitian([sparse(Matrix{T}(I, m, m)) A; A' spzeros(n,n)])
    K = copy(K_0)

    L, B, vec_P = bunch_parlett(K)
    vec_P_transposed = transpose_vec_P(vec_P)
    K = L*B*L'
    K = perm_line_vec_P(K, vec_P_transposed)
    K = perm_column_vec_P(K, vec_P_transposed)
    @test norm(K - K_0) < 1e-10
    @test eltype(L) == T && eltype(B) == T
end


for i in 1:N
    A = randn(BigFloat, m, n) + im * randn(BigFloat, m, n)
    T = eltype(A)
    K_0 = Hermitian([sparse(Matrix{T}(I, m, m)) A; A' spzeros(n,n)])
    K = copy(K_0)

    L, B, vec_P = bunch_parlett(K)
    vec_P_transposed = transpose_vec_P(vec_P)
    K = L*B*L'
    K = perm_line_vec_P(K, vec_P_transposed)
    K = perm_column_vec_P(K, vec_P_transposed)
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
    L, B = extract_L_and_B(K, vec_2by2)
    vec_P_transposed = transpose_vec_P(vec_P)
    K = L*B*L'
    K = perm_line_vec_P(K, vec_P_transposed)
    K = perm_column_vec_P(K, vec_P_transposed)
    K = LowerTriangular(K)
    @test norm(K - K_0) < 1e-10
end

for i in 1:N
    A = randn(ComplexF64, m, n)
    T = eltype(A)
    K_0 = LowerTriangular([sparse(Matrix{T}(I, m, m)) A; A' spzeros(n,n)])
    K = copy(K_0)

    vec_P, vec_2by2 = bunch_parlett!(K)
    L, B = extract_L_and_B(K, vec_2by2)
    vec_P_transposed = transpose_vec_P(vec_P)
    K = L*B*L'
    K = perm_line_vec_P(K, vec_P_transposed)
    K = perm_column_vec_P(K, vec_P_transposed)
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
    L, B = extract_L_and_B(K, vec_2by2)
    vec_P_transposed = transpose_vec_P(vec_P)
    K = L*B*L'
    K = perm_line_vec_P(K, vec_P_transposed)
    K = perm_column_vec_P(K, vec_P_transposed)
    K = LowerTriangular(K)
    @test norm(K - K_0) < 1e-10
    @test eltype(K) == T && eltype(L) == T && eltype(B) == T
end
