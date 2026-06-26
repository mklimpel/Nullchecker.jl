module Nullchecker

# Exports
# Types
export Vertexindex # Vertexindex.jl
export Predecessor, Visited, Startzero, Forward, Commutation, Backward, SumShort, SumLong, Backwardedge # Predecessor.jl
export Monomial, Vertex, Monomialgraph # Monomialgraph.jl
export Tieredqueue, Vertexqueue # Vertexqueue.jl

# Functions
export getvertex, createk1monomialgraph, createnextmonomialgraph, createmonomialgraph, computemonomialindex, computevertexindex, getcommute, setcommute!, checkhypervertex # Monomialgraph.jl
export createcompletegraph, createsixpointgraph, createpetersengraph, createfourpointgraphs # Graphgenerators.jl
export distancematrix, startzeros! # Startzeros.jl
export traverse! # Traverse.jl
export postprocessing!, abequalaba # Postprocessing.jl
export resetvisit!, iteratetraverse! # Iteratetraverse.jl

using LinearAlgebra
using MaxPlus
using TimerOutputs

include("NullcheckerIncludes/Vertexindex.jl")
include("NullcheckerIncludes/Predecessor.jl")
include("NullcheckerIncludes/Monomialgraph.jl")
include("NullcheckerIncludes/Vertexqueue.jl")
include("NullcheckerIncludes/Graphgenerators.jl")
include("NullcheckerIncludes/Startzeros.jl")
include("NullcheckerIncludes/Traverse.jl")
include("NullcheckerIncludes/Postprocessing.jl")
include("NullcheckerIncludes/Iteratetraverse.jl")

function setupsixpoint()
	g = createsixpointgraph()
	k1 = createk1monomialgraph(g)
	k2 = createnextmonomialgraph(k1)
	k3 = createnextmonomialgraph(k2)
	q = Vertex[]
	return g, k3, q
end

function setuppetersen()
	g = createpetersengraph()
	k1 = createk1monomialgraph(g)
	k2 = createnextmonomialgraph(k1)
	k3 = createnextmonomialgraph(k2)
	q = Vertex[]
	return g, k3, q
end

end
