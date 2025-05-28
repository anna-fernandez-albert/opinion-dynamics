module LFRAnalysis

include("./model_simulation.jl")
include("../utils/lfr_generator.jl")
include("../utils/network_utils.jl")
using PyCall
using Statistics
using .OpinionSimulation
using .LFRGenerator
using .NetworkUtils

export run_lfr_analysis
#-----------------------------------------------------------------------------------------------------------------------------------
function run_lfr_analysis()
    
    N = [100, 200, 300, 400, 500]
    μ = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    k_avg = [5, 10, 15, 20, 25]
    
    # Create networks
    networks = LFRGenerator.generate_lfr_networks(N, μ, k_avg)

    results = Dict()
    for (name, network) in networks
        println("Executant simulació i anàlisi per a '$name'...")
        results[name] = OpinionSimulation.run_model_analysis(
            network["network"], 
            network["communities"], 
            "lfr_$(name)"
        )
    end
    
    return networks, results
end
end