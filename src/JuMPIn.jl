module JuMPIn
# TODO: Explicitly export functions?

# Methods to identify variables
include("identify_variables.jl") # IdentifyVariables
identify_unique_variables = IdentifyVariables.identify_unique_variables

# Methods to identify equality constraints
include("get_equality.jl") # GetEquality
get_equality_constraints = GetEquality.get_equality_constraints

# Methods to get the incidence graph of constraints and variables
include("incidence_graph.jl") # IncidenceGraph
get_bipartite_incidence_graph = IncidenceGraph.get_bipartite_incidence_graph

# Methods to apply graph algorithms to JuMP models
include("interface.jl") # Interface
IncidenceGraphInterface = Interface.IncidenceGraphInterface
get_adjacent = Interface.get_adjacent
maximum_matching = Interface.maximum_matching
dulmage_mendelsohn = Interface.dulmage_mendelsohn

end
