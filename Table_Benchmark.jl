using PrettyTables

algos = [
  "bunch_parlett",
  "bunch_parlett!",
  "bunch_kaufman",
  "bunch_kaufman!",
  "bunch_kaufman_rook",
  "bunch_kaufman_rook!",
  "bunchkaufman LAPACK",
  "bunchkaufman! LAPACK"
]

headers = ["Algorithme", "75x75", "100x100", "250x250", "500x500", "1000x1000"]

allocs = Dict(
  "bunch_parlett" => [18, 18, 18, 28, 28],
  "bunch_parlett!" => [9, 9, 9, 16, 16],
  "bunch_kaufman" => [18, 18, 18, 28, 28],
  "bunch_kaufman!" => [10, 10, 11, 18, 18],
  "bunch_kaufman_rook" => [18, 18, 18, 26, 26],
  "bunch_kaufman_rook!" => [11, 10, 10, 16, 16],
  "bunchkaufman LAPACK" => [69, 98, 222, 472, 954],
  "bunchkaufman! LAPACK" => [63, 92, 216, 471, 953]
)

memory = Dict(
  "bunch_parlett" => round.(Int, [162624, 306912, 1585472, 10957440, 43737728] ./ 1024),
  "bunch_parlett!" => round.(Int, [115504, 222272, 1075376, 8941920, 35709280] ./ 1024),
  "bunch_kaufman" => round.(Int, [162624, 306912, 1585472, 10957440, 43737728] ./ 1024),
  "bunch_kaufman!" => round.(Int, [115600, 222368, 1075840, 8942384, 35709744] ./ 1024),
  "bunch_kaufman_rook" => round.(Int, [162624, 306912, 1585472, 10957376, 43737664] ./ 1024),
  "bunch_kaufman_rook!" => round.(Int, [115968, 222368, 1075472, 8942320, 35709680] ./ 1024),
  "bunchkaufman LAPACK" => round.(Int, [159840, 309184, 1570784, 5174976, 20607264] ./ 1024),
  "bunchkaufman! LAPACK" => round.(Int, [118112, 226272, 1085328, 8963776, 35754272] ./ 1024)
)

time = Dict(
  "bunch_parlett" => round.([5.465413e6, 1.4511723e7, 2.79721081e8, 2.959162298e9, 3.7112352078e10], digits=2),
  "bunch_parlett!" => round.([5.925377e6, 1.5993916e7, 2.81235827e8, 3.656638106e9, 3.6181912817e10], digits=2),
  "bunch_kaufman" => round.([4.781777e6, 1.3071789e7, 2.48297187e8, 2.915699917e9, 3.3871329546e10], digits=2),
  "bunch_kaufman!" => round.([4.80889e6, 1.3619278e7, 2.43764668e8, 2.966110703e9, 3.5017727643e10], digits=2),
  "bunch_kaufman_rook" => round.([4.854916e6, 1.3245779e7, 2.68957728e8, 3.354467484e9, 3.4739666518e10], digits=2),
  "bunch_kaufman_rook!" => round.([4.656729e6, 1.2921323e7, 2.64627048e8, 3.224862225e9, 3.5328886991e10], digits=2),
  "bunchkaufman LAPACK" => round.([3.426366e6, 9.644988e6, 2.22131215e8, 2.953089352e9, 3.3963627936e10], digits=2),
  "bunchkaufman! LAPACK" => round.([3.245781e6, 1.0127854e7, 2.20212231e8, 2.978141356e9, 3.5167386841e10], digits=2)
)

function build_table(metric)
    nrows = length(algos)
    ncols = length(headers)
    table = Matrix{Any}(undef, nrows, ncols)
    for (i, algo) in enumerate(algos)
        table[i, 1] = algo
        for j in 1:5
            table[i, j+1] = metric[algo][j]
        end
    end
    return table
end

println("Nombre d'allocations :")
pretty_table(build_table(allocs), header=headers)

println("Mémoire allouée (KiB) :")
pretty_table(build_table(memory), header=headers)

println("Temps d'exécution minimum (ns) :")
pretty_table(build_table(time), header=headers, formatters = ft_printf("%0.2e", 2:6))