include("./src/constants.jl")
include("./src/utils/network_utils.jl")
include("./src/utils/lfr_generator.jl")
include("./src/simulations/model_parameters_analysis.jl")
include("./src/simulations/lfr_analysis.jl")

using .NetworkUtils
using .LFRGenerator
using .ModelAnalysis
using .LFRAnalysis
#-----------------------------------------------------------------------------------------------------------------------------------
# Create directories if they don't exist
mkpath(PATH_TO_NETWORKS)
mkpath(PATH_TO_PLOTS)
mkpath(PATH_TO_RESULTS)

function main()
    # Configuration of the study
    run_sensitivity_analysis = false
    run_lfr_analysis = true
    
    # 1. Model analysis
    #    - Study of the results due to the model parameters
    if run_sensitivity_analysis
        println("\n=== MODEL PARAMETERS SENSIBILITY ANALYSIS ===")

        model, parameters = "LFR", [5000, 0.2, 100]
        name = "LFR_N_5000_mu_$(parameters[2])_k_avg_100"
        println("Network name: $name")
        mkpath("$PATH_TO_PLOTS/$(SENSITIVITY_ANALYSIS)_$(name)")
        if isfile("$PATH_TO_NETWORKS/$(name).lgz")
            graph, communities = NetworkUtils.load_network(name)
        else
            println("Generating network '$name'...")
            graph, communities = NetworkUtils.generate_network(model, parameters, name)
        end

        ModelAnalysis.run_model_sensibility_analysis(graph, communities, name)
        println("Sensibility analysis finished")
    end
    
    # 2. Opinion dynamics with LFR networks
    #    - Study of the opinion dynamics due to the network structure
    if run_lfr_analysis
        println("\n=== OPINION DYNAMICS WITH LFR NETWORKS ===")
        mkpath("$PATH_TO_PLOTS/$(LFR_ANALYSIS)")
        
        LFRAnalysis.run_lfr_analysis()
        println("LFR analysis finished")
    end
    
    println("\nFINISHING SIMULATION")
end

#-----------------------------------------------------------------------------------------------------------------------------------
main()