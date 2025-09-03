# Paths
const PATH_TO_NETWORKS = "./data/networks"
const PATH_TO_RESULTS = "./data/results"
const PATH_TO_PLOTS = "./plots"
const SENSITIVITY_ANALYSIS = "sensitivity_analysis"
const LFR_ANALYSIS = "lfr_analysis"


# Model Parameters
const OPINION_VALUES = [-1, 1]
const NUM_STEPS = 5000
const TOLERANCE_STEPS = 50

const PROBABILITY_MAJORITY_OPINION = 0.7
const TRUST_LEVEL = 0.8
const λ_DEFAULT = 0.5

# Simulation Parameters - SENSIBILITY ANALYSIS
const λ_VALUES = collect(0.1:0.1:1.0)              # Sensitivity to opinion distance influence in the opinion interaction
const LOCAL_TRUST_VALUES = collect(0.1:0.1:1.0)    # Trust level for local interactions. From 0.1 to 1 in steps of 0.1. If the trust level is low, the node is more resilient to change.
const PROB_MAJORITY_OPINION_VALUES = [0.7]         # Probability of a node adopting the majority opinion of its community

# Simulation Parameters - LFR ANALYSIS
const N_VALUES = [5000]
const μ_VALUES = collect(0.1:0.1:1.0)
const k_avg_VALUES = [100]