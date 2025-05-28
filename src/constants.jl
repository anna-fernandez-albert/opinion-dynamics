# Paths
const PATH_TO_NETWORKS = "./data/networks"
const PATH_TO_RESULTS = "./data/results"
const PATH_TO_PLOTS = "./plots"
const NETWORK_ANALYSIS_RESULTS_FILE_NAME = "network_analysis_results"
const NETWORK_ANALYSIS_PLOTS = "network_analysis_plots"
const SENSITIVITY_ANALYSIS_RESULTS_FILE_NAME = "sensitivity_analysis_results"

# Model Parameters
const OPINION_RANGE = 0:2π
const OPINION_VALUE = [-1, 1]
const GLOBAL_INFLUENCE = 1.0
const PROBABILITY_MAJORITY_OPINION = 0.5
const TRUST_LEVEL = 0.5
const NUM_STEPS = 10000
const GLOBAL_PERIOD = 1000
const λ = 0.5

# Simulation Parameters
const MU_VALUES = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]