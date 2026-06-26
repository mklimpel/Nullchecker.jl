"""
	distancematrix(g::Matrix{<:Integer})

Compute the distance matrix of the graph given by adjacency matrix g. Uses the MaxPlus.jl package.
"""
function distancematrix(g::Matrix{<:Integer})
	distance = [(x==1 ? x : Inf) for x in g]
	distance[CartesianIndex.(axes(distance,1),axes(distance,2))] .= 0
	distance = plustimes(MI(distance)^size(distance,1))
end

"""
	startzeros!(q::Vertexqueue, m::Monomialgraph, g::Matrix{<:Integer})

Compute starting zeros for the traversal algorithm. When a monomial is determined to be zero through one of the following rules, the corresponding vertex in k is added to the queue q.

Rules:
(1) u_{i,j} u_{k,l} = 0 if i~k and !k~l
(2) u_{i,j} u_{k,l} = 0 if !i~k and k~l
(3) u_{i,j} u_{i,l} = 0 if j≠l
(4) u_{i,j} u_{k,j} = 0 if i≠k
(5) u_{i,j} u_{k,l} = 0 if d(i,k)≠d(j,l), where d(i,k) is the distance of vertex i to vertex j
(6) u_{i,j} = 0 		if deg(i)≠deg(j) where deg(i) is the degree of vertex i

Rule (5) is stronger than rules (1)-(4) and therefore contains all these cases.
"""
function startzeros!(q::Vertexqueue, m::Monomialgraph, g::Matrix{<:Integer})
	n = m.n
	d = distancematrix(g)
	
	# Rule (5)
	if m.k >= 2
		for i=1:n, j=1:n, k=1:n, l=1:n
			if d[i,k] != d[j,l]
				w = m.vertexlist[2][computemonomialindex([(i,j),(k,l)], n)]
				if typeof(w.visited) == Nothing
					w.visited = Startzero(:distance)
					push!(q, w)
				end
			end
		end
	end

	# Rule (6)
	degrees = sum(g, dims=2)
	for i=1:n, j=1:n
		if degrees[i] != degrees[j]
			w = m.vertexlist[1][computemonomialindex([(i,j)], n)]
			if typeof(w.visited) == Nothing
				w.visited = Startzero(:degree)
				push!(q,w)
			end
		end
	end
end

