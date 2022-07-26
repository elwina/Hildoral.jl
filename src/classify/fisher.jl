using Statistics
using LinearAlgebra
using Plots
using LaTeXStrings
gr()

function Fisher(train::Matrix{<:Real}, label::Vector{Int})
    μ_0 = mean(train[label.==0, :], dims=1)'
    μ_1 = mean(train[label.==1, :], dims=1)'
    Σ_0 = zeros(Float64, size(train, 2), size(train, 2))
    Σ_1 = zeros(Float64, size(train, 2), size(train, 2))
    for row in eachrow(Σ_0)
        x_i = row'
        Σ_0 += (x_i .- μ_0) * ((x_i .- μ_0)')
        Σ_1 += (x_i .- μ_1) * ((x_i .- μ_1)')
    end
    S_w = Σ_0 + Σ_1
    U, S, V = svd(S_w)
    S_w_1 = V * Diagonal(S) .* (-1) * U'
    ω = (S_w_1*(μ_0-μ_1))[:, 1]

    distance = Vector{Float64}(undef, 0)
    for row in eachrow(train)
        t = ω' * row
        distance = [distance; t]
    end

    m0 = mean(distance[label.==0, :])
    m1 = mean(distance[label.==1, :])
    middle = mean([m0, m1])

    jud(X) = (ω' * X) < middle ? 0 : 1

    predre = Vector{Bool}(undef, 0)
    errnum = 0
    for i in 1:(size(train)[1])
        row = train[i, :]
        println(ω' * row)
        if jud(row) == label[i]
            predre = [predre; true]
        else
            predre = [predre; false]
            errnum += 1
        end
    end

    err = errnum / length(predre)
    return (ω, middle, jud, predre, err)
end

function Fisher2DPlot(train, label, ω, sample=undef)
    figure = plot()
    xlabel!(figure, L"x_1")
    ylabel!(figure, L"x_2")
    scatter!(figure, train[label.==0, 1], train[label.==0, 2], label="group 0", color="blue")
    scatter!(figure, train[label.==1, 1], train[label.==1, 2], label="group 1", color="green")

    p = (ω * ω') ./ (ω' * ω)

    projection = Matrix{Float64}(undef, 0, 2)
    for row in eachrow(train)
        projection = [projection; (p*row)[1] (p*row)[2]]
    end
    scatter!(figure, projection[label.==0, 1], projection[label.==0, 2], label=false, color="blue")
    scatter!(figure, projection[label.==1, 1], projection[label.==1, 2], label=false, color="green")

    k = (ω[2] / ω[1])
    linex = range(-5, 12, length=100)
    plot!(figure, linex, k * linex, label="line", color="black")

    middle = mean([mean(projection[label.==0, 1]); mean(projection[label.==1, 1])])
    scatter!(figure, [middle], [k * middle], label="point", color="red")

    if sample != undef
        scatter!(figure, sample[:, 1], train[:, 2], label=false, color="cyan")
        tprojection = Matrix{Float64}(undef, 0, 2)
        for row in eachrow(sample)
            tprojection = [tprojection; (p*row)[1] (p*row)[2]]
        end
        scatter!(figure, tprojection[:, 1], tprojection[:, 2], label="sample", color="cyan")
    end

end