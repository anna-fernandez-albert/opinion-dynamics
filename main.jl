include("./src/constants.jl")
include("./src/utils/network_utils.jl")
include("./src/simulations/model_simulation.jl")
include("./src/simulations/model_parameters_analysis.jl")
include("./src/simulations/lfr_analysis.jl")

using .NetworkUtils
using .OpinionSimulation
using .ModelAnalysis
using .LFRAnalysis
#-----------------------------------------------------------------------------------------------------------------------------------
# Create directories if they don't exist
mkpath(PATH_TO_NETWORKS)
mkpath(PATH_TO_PLOTS)
mkpath(PATH_TO_RESULTS)

function main()
    # Configuration of the study
    run_unique_network = false
    run_sensitivity_analysis = false
    run_lfr_analysis = true
    
    # 1. Network analysis
    #   - Study of the opinion convergence/polarization evolution
    #   - Study of the polarization/convergence as a function of the model parameters
    #   - Study the opnion dynamics for bridge nodes
    if run_unique_network
        println("\n=== NETWORK ANALYSIS ===")
        
        name = "synthetic_network_N_300_blocks_5_prr_0.08_prs_0.02"
        nx_graph, communities = NetworkUtils.load_network(name)
        results = OpinionSimulation.run_model_analysis(nx_graph, communities, name)

        NetworkUtils.save_network_analysis_results(results, name)
        println("Network analysis finished")
    end
    
    # 2. Model analysis
    #    - Study of the results due to the model parameters
    if run_sensitivity_analysis
        println("\n=== MODEL PARAMETERS SENSITIBITY ANALYSIS ===")

        name = "synthetic_network_N_300_blocks_5_prr_0.08_prs_0.02"
        nx_graph, communities = NetworkUtils.load_network(name)
        sensitivity_results = ModelAnalysis.run_model_sensitibily_analysis(nx_graph, communities)

        NetworkUtils.save_sensitivity_analysis_results(sensitivity_results, name)
        println("Sensibility analysis finished")
    end
    
    # 3. Opinion dynamics with LFR networks
    #    - Study of the opinion dynamics due to the network structure
    if run_lfr_analysis
        println("\n=== OPINION DYNAMICS WITH LFR NETWORKS ===")
        
        networks, lfr_results = LFRAnalysis.run_lfr_analysis()

        NetworkUtils.save_lfr_analysis(networks, lfr_results)
        println("LFR analysis finished")
    end
    
    println("\nFINISHING SIMULATION")
end

#-----------------------------------------------------------------------------------------------------------------------------------
main()