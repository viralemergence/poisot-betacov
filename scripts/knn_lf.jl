using CSV
using DataFrames

using EcologicalNetworks

using StatsBase

## Load the data and aggregate everything
csv_file = download("https://raw.githubusercontent.com/ViromeNet/cleanbats_betacov/master/clean%20data/BatCoV-assoc_clean.csv?token=AAENOAM6NGFDAOTL5IK2UJS6WV5LY")
virion = CSV.read(csv_file)
select!(virion, Not(:origin))
select!(virion, Not(:Column1))
select!(virion, Not(:host_species))
rename!(virion, :clean_hostnames => :host_species)
virion = unique(virion)

## Prepare the network
hosts = unique(virion.host_species)
viruses = unique(virion.virus_genus)
A = zeros(Bool, (length(viruses), length(hosts)))
U = BipartiteNetwork(A, viruses, hosts)

for interaction in eachrow(virion)
    U[interaction.virus_genus, interaction.host_species] = true
end

## kNN preparation
tanimoto(x::Set{T}, y::Set{T}) where {T} = length(x∩y)/length(x∪y)

## Main loop?
function knn_virus(train::T, predict::T; k::Integer=5, cutoff::Integer=1) where {T <: BipartiteNetwork}
    predictions = DataFrame(virus = String[], host = String[], match = Float64[])
    for s in species(predict; dims=1)
        @assert s in species(train) "The species $s is absent from the predicted network"
        hosts = train[s,:]
        neighbors = Dict([neighbor => tanimoto(hosts, train[neighbor,:]) for neighbor in filter(x -> x != s, species(train; dims=1))])
        top_k = sort(collect(neighbors), by=x->x[2], rev=true)[1:k]
        hosts_count = StatsBase.countmap(vcat(collect.([predict[n.first,:] for n in top_k])...))
        likely = filter(p -> p.second >= cutoff, sort(collect(hosts_count), by=x->x[2], rev=true))
        for l in likely
            l.first ∈ predict[s, :] && continue
            push!(predictions,
                (s, l.first, l.second/k)
            )
        end
    end
    return predictions
end


## Write this shit
predict_path = joinpath(pwd(), "predictions", "knn")
ispath(predict_path) || mkpath(predict_path)

## Predict and write
knn = knn_virus(U, U)

CSV.write(
    joinpath(predict_path, "PoisotTanimoto.csv"),
    knn;
    writeheader=false
)

## Linear filtering path
lf_path = joinpath(pwd(), "predictions", "linearfilter")
ispath(lf_path) || mkpath(lf_path)

## Linear filtering
predictions_lf = DataFrame(species=String[], score=Float64[])

α = [0.0, 1.0, 1.0, 1.0]

for i in interactions(linearfilter(U; α=α))
    U[i.from, i.to] && continue
    i.to ∈ species(U; dims=2) || continue
    if i.from == "Betacoronavirus"
        push!(predictions_lf, 
            (replace(i.to, " "=>"_"), i.probability)
        )
    end
end

sort!(predictions_lf, :score, rev=true)

CSV.write(
    joinpath(lf_path, "PoisotLinearFilter.csv"),
    predictions_lf;
    writeheader=false
)


## Do some LOO just for fun

#=
success = 0
attempts = 0
k = 8
for i in interactions(M)
    global success
    global attempts
    K = copy(M)
    K[i.from, i.to] = false
    simplify!(K)
    if richness(K) != richness(M)
        continue
    end
    neighbors = filter(x -> x != i.from, species(K; dims=1))
    scores = [tanimoto(K[i.from,:], K[neighbor,:]) for neighbor in neighbors]
    nearest_neighbors = neighbors[StatsBase.partialsortperm(scores, 1:k)]
    if i.to in keys(degree(simplify(K[nearest_neighbors, :]); dims=2))
        success += 1
    end
    attempts += 1
end
@info "$k \t $(success/attempts)"
=#