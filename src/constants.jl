# Paths
const PATH_TO_NETWORKS = "./data/networks"
const PATH_TO_RESULTS = "./data/results"
const PATH_TO_PLOTS = "./plots"
const SENSITIVITY_ANALYSIS = "sensitivity_analysis"
const LFR_ANALYSIS = "lfr_analysis"


# Model Parameters
const OPINION_VALUES = [-1, 1]
const GLOBAL_INFLUENCE = 1.0

const PROBABILITY_MAJORITY_OPINION = 0.7
const TRUST_LEVEL = 0.4
const λ_DEFAULT = 0.2

# Simulation Parameters - SENSIBILITY ANALYSIS
const λ_VALUES = collect(0.05:0.01:1.0)                               # Sensitivity to opinion distance influence in the opinion interaction
const LOCAL_TRUST_VALUES = collect(0.1:0.05:1.0)                      # Trust level for local interactions. From 0.1 to 1 in steps of 0.1. If the trust level is low, the node is more resilient to change.
const GLOBAL_INFLUENCE_PERIODS = [500, 1500, 3000, 5000]              # Global influence period. The global influence is applied every X steps
const PROB_MAJORITY_OPINION_VALUES = [0.1, 0.3, 0.5, 0.7, 0.9, 1.0]   # Probability of a node adopting the majority opinion of its community

# Simulation Parameters - LFR ANALYSIS
const N_VALUES = [100, 200, 300, 400]
const μ_VALUES = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
const k_avg_VALUES = [5, 10, 15, 20, 25]