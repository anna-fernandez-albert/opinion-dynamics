module ModelAnalysis
include("../constants.jl")
include("../utils/network_utils.jl")
include("../model/opinion_dynamics.jl")
include("../visualization/opinion_plots.jl")

using Statistics
using Graphs
using CSV, DataFrames

using .NetworkUtils
using .OpinionDynamics
using .Visualization

export run_model_sensibility_analysis
#-----------------------------------------------------------------------------------------------------------------------------------
function run_model_sensibility_analysis(graph, communities, name)
    size = nv(graph)
    repetitions = size < 5000 ? 100 : 20
    
    for pmo_val in PROB_MAJORITY_OPINION_VALUES
        for trust_val in LOCAL_TRUST_VALUES
            for λ_val in λ_VALUES
                println("\n\nRunning simulation with pmo=$pmo_val, trust=$trust_val, λ=$λ_val\n")

                for i in 1:repetitions
                    println("Repetition $i of $repetitions")
                    if isfile("$PATH_TO_RESULTS/$(SENSITIVITY_ANALYSIS)_$(name)_parameters_$(pmo_val)_$(trust_val)_$(λ_val)_repetition_$(i).csv")
                        println("File already exists")
                    else
                        current_time = time()
                        initial_opinions = OpinionDynamics.initialize_opinions(size, communities, OPINION_VALUES, pmo_val)

                        opinion_history, t_execution = OpinionDynamics.run_simulation(graph, initial_opinions, communities, NUM_STEPS, λ_val, trust_val, TOLERANCE_STEPS)
                        print("Execution time: $t_execution steps")
                        df = DataFrame(opinion_history, :auto)
                        CSV.write("$PATH_TO_RESULTS/$(SENSITIVITY_ANALYSIS)_$(name)_parameters_$(pmo_val)_$(trust_val)_$(λ_val)_repetition_$(i).csv", df)
                        print("Time taken: $(time() - current_time) seconds")
                    end
                end
            end
        end
    end
end
#-----------------------------------------------------------------------------------------------------------------------------------
end