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
    run_sensitivity_analysis = true
    run_lfr_analysis = false
    
    # 1. Model analysis
    #    - Study of the results due to the model parameters
    if run_sensitivity_analysis
        println("\n=== MODEL PARAMETERS SENSIBILITY ANALYSIS ===")

        name = "N_100_mu_0.9_k_avg_25"
        mkpath("$PATH_TO_PLOTS/$(SENSITIVITY_ANALYSIS)_$(name)")
        if isfile("$PATH_TO_NETWORKS/lfr_$(name).net")
            nx_graph, communities = LFRGenerator.load_lfr_network(name, PATH_TO_NETWORKS)
        else
            println("Generating network '$name'...")
            nx_graph, communities = LFRGenerator.generate_lfr_networks([100], [0.5], [20])
        end

        ModelAnalysis.run_model_sensibility_analysis(nx_graph, communities, name)
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