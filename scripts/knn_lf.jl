using CSV
using DataFrames
using EcologicalNetworks
using StatsBase

#--- Load the data and aggregate everything

virion = CSV.read(joinpath("data", "virionette.csv"))

#--- Function to make a network

function mknet(df)
    hosts = unique(df.host_species)
    viruses = unique(df.virus_genus)
    A = zeros(Bool, (length(viruses), length(hosts)))
    U = BipartiteNetwork(A, viruses, hosts)

    for interaction in eachrow(df)
        U[interaction.virus_genus, interaction.host_species] = true
    end
    return simplify(U)
end

#--- Prepare the network

ALL = mknet(virion)

BATS = mknet(virion[virion.host_order.=="Chiroptera", :])

#--- kNN preparation

function tanimoto(x::Set{T}, y::Set{T}) where {T}
    length(x∩y)/length(x∪y)
end

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

#--- Make the folders

knn_path = joinpath(pwd(), "predictions", "knn")
ispath(knn_path) || mkpath(knn_path)
lf_path = joinpath(pwd(), "predictions", "linearfilter")
ispath(lf_path) || mkpath(lf_path)

#--- Linear filtering

lf_bats = DataFrame(species=String[], score=Float64[])
lf_all = DataFrame(species=String[], score=Float64[])

α = [0.0, 1.0, 1.0, 1.0]

for i in interactions(linearfilter(BATS; α=α))
    BATS[i.from, i.to] && continue
    i.to ∈ species(BATS; dims=2) || continue
    if i.from == "Betacoronavirus"
        push!(lf_bats,
            (replace(i.to, " "=>"_"), i.probability)
        )
    end
end

for i in interactions(linearfilter(ALL; α=α))
    ALL[i.from, i.to] && continue
    i.to ∈ species(ALL; dims=2) || continue
    if i.from == "Betacoronavirus"
        push!(lf_all,
            (replace(i.to, " "=>"_"), i.probability)
        )
    end
end

sort!(lf_all, :score, rev=true)
sort!(lf_bats, :score, rev=true)

CSV.write(
    joinpath(lf_path, "PoisotLinearFilter-Bats.csv"),
    lf_bats;
    writeheader=false
)

CSV.write(
    joinpath(lf_path, "PoisotLinearFilter-Mammals.csv"),
    lf_all;
    writeheader=false
)

#--- kNN for hosts AND viruses

function Base.transpose(N::T) where {T <: EcologicalNetworks.AbstractBipartiteNetwork}
    B = zeros(eltype(N.A), reverse(size(N)))
    Y = typeof(N)(B, species(N; dims=2), species(N; dims=1))
    for i in interactions(N)
        Y[i.to, i.from] = N[i.from, i.to]
    end
    return Y
end

host_knn_bats = knn_virus(transpose(BATS), transpose(BATS))

rename!(host_knn_bats, :host => :virusname)
rename!(host_knn_bats, :virus => :hostname)

host_knn_bats = host_knn_bats[host_knn_bats.virusname .== "Betacoronavirus", :]
sort!(host_knn_bats, :match, rev=true)



host_knn_all = knn_virus(transpose(ALL), transpose(ALL))

rename!(host_knn_all, :host => :virusname)
rename!(host_knn_all, :virus => :hostname)

host_knn_all = host_knn_all[host_knn_all.virusname .== "Betacoronavirus", :]
sort!(host_knn_all, :match, rev=true)

CSV.write(
    joinpath(knn_path, "PoisotTanimotoHosts-Bats.csv"),
    host_knn_bats;
    writeheader=false
)

CSV.write(
    joinpath(knn_path, "PoisotTanimotoHosts-Mammals.csv"),
    host_knn_all;
    writeheader=false
)

knn_bats = knn_virus(BATS, BATS)
knn_bats = knn_bats[knn_bats.virus .== "Betacoronavirus", :]
sort!(knn_bats, :match, rev=true)

knn_all = knn_virus(ALL, ALL)
knn_all = knn_all[knn_all.virus .== "Betacoronavirus", :]
sort!(knn_all, :match, rev=true)

CSV.write(
    joinpath(knn_path, "PoisotTanimotoVirus-Bats.csv"),
    knn_bats;
    writeheader=false
)

CSV.write(
    joinpath(knn_path, "PoisotTanimotoVirus-Mammals.csv"),
    knn_all;
    writeheader=false
)

## Do some LOO just for fun
# UPDATE This was not fun

#=
success = 0
attempts = 0
k = 11
M = transpose(U)
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
