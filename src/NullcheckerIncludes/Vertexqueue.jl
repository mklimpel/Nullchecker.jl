"""
	Tieredqueue(k)

Data structure for a queue of vertices, where each vertex is put into an individual queue that depends on the length of its associated monomial.
"""
mutable struct Tieredqueue
	const tiers::Vector{Vector{Vertex}}
	const k::Integer
	length::Integer
end

function Tieredqueue(k::Integer)
	tiers = Vector{Vector{Vertex}}(undef, k)
	for i = 1:k
		tiers[i] = Vertex[]
	end
	return Tieredqueue(tiers, k, 0)
end

function Base.show(io::IO, q::Tieredqueue)
	for i = 1:q.k-1
		println(io, q.tiers[i])
	end
	print(io, q.tiers[q.k])
end

"""
	push!(q::Tieredqueue, v::Vertex)

Like push! for regular containers. Insert the vertex v into the correct bucket in q.
"""
function Base.push!(q::Tieredqueue, v::Vertex)
	k = v.monomial.length
	q.length = q.length + 1
	push!(q.tiers[k], v)
end

"""
	pop!(q::Tieredqueue)

Like pop! for regular containers. Return the shortest possible element in the queue.
"""
function Base.pop!(q::Tieredqueue)
	for i = 1:q.k
		if isempty(q.tiers[i]) == false
			q.length = q.length - 1
			return pop!(q.tiers[i])
		end
	end
	throw(ArgumentError("Tieredqueue must be non-empty"))
end

Base.length(q::Tieredqueue) = q.length

Base.isempty(q::Tieredqueue) = (length(q) == 0)

"""
Union type for Vector{Vertex} and Tieredqueue. Used in the definitions of functions that allow for both types.
"""
Vertexqueue = Union{Tieredqueue, Vector{Vertex}}
