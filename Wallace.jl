using FreqTables: freqtable

function simpsons(observations::Vector)::Float32

    n = length(observations)

    classes = unique(observations)

    species_counts = map(species -> sum(observations .== species), classes)

    proportions = map(species_count -> (species_count / n)^2, species_counts)

    1.0 - sum(proportions)

end

function wallace(A, B)::Float32

    mm = mismatch_matrix(A, B)

    mm.a / (mm.a + mm.b)

end

function adjwallace(A, B)::Float32

    wi = wallacei(A, B)

    w = wallace(A, B)

    (w - wi) / (1 - wi)
end

function wallacei(A, B)::Float32
    1.0 - simpsons(B)
end

function adjrand(A, B)::Float32

    mm = mismatch_matrix(A, B)

    nis = sum(mm.rowsums .^ 2)
    njs = sum(mm.colsums .^ 2)

    nc =
        (
            mm.n * (mm.n^2 + 1) - (mm.n + 1) * nis - (mm.n + 1) * njs +
            2 * (nis * njs) / mm.n
        ) / (2 * (mm.n - 1))

    max((mm.a + mm.d - nc) / (mm.a + mm.b + mm.c + mm.d - nc), 0.0)

end

function mismatch_matrix(
    A,
    B,
)

    # contingency table
    # needs the argument reversal to match Comparing Partitions
    ct = freqtable(B, A)

    colsums = map(sum, eachcol(ct))
    rowsums = map(sum, eachrow(ct))

    n = sum(ct)

    a::Int = sum((ct .* (ct .- 1) / 2))
    a′::Int = sum((colsums .* (colsums .- 1) / 2))
    b::Int = a′ - a
    c::Int = sum((rowsums .* (rowsums .- 1) / 2)) - a
    d::Int = (n * (n - 1) / 2) - a′ - c

    (a = a, b = b, c = c, d = d, n = n, colsums = colsums, rowsums = rowsums)
end
