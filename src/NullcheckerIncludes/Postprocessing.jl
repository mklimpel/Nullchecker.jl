"""
	postprocessing!(m::Monomialgraph)

Run post processing procedures after a single iteration of the traverse! algorithm. The post processing consists of the following:

(1) Check if ab == aba for generators a and b and if so set a and b as commuting in the commutation table.
(2) Check if u_{i,j} == 1 by checking if all u_{i,k} == 0 except for k == j or if all u_{k,j} == 0 except for k == 1
"""
function postprocessing!(mgraph::Monomialgraph)
	k = mgraph.k
	n = mgraph.n
	if k > 2
		# Process (1)
		filterfunction = generatepostfilter(mgraph)
		# Filter out monomials m = ab that are 0 and ones where ab == ba is known already
		for v in filter(filterfunction, mgraph.vertexlist[2]) 
			if abequalaba(v, mgraph)
				#println(v, " abequalaba")
				setcommute!(v.monomial, mgraph, true)
			end
		end
	end

	# Process (2)
	for v in mgraph.vertexlist[1]
		if typeof(v.visited) != Nothing
			continue
		else
			for position in ['i', 'j']
				w = v
				nonzerocount = 0
				for k = 1:n
					w = mgraph.vertexlist[1][computenextindex(w.monomial, position, 1, n)]
					nonzerocount += typeof(w.visited) != Nothing ? 0 : 1
				end
				if nonzerocount == 1
					# v is 1, set corresponding row and column in commutationtable to be true
					#println(v, " is 1")
					index = computemonomialindex(v.monomial, n)
					mgraph.commutationtable[index,:] .= mgraph.commutationtable[:,index] .= true
				end
			end
		end
	end
end

"""
	generatepostfilter(mgraph::Monomialgraph)

Create a function that takes a vertex corresponding to a monomial m = ab of length 2 and returns false if the vertex has been visited or if the generators a and b are known to commute in mgraph.
"""
function generatepostfilter(mgraph::Monomialgraph)
	x -> !(typeof(x.visited) != Nothing || getcommute(x.monomial, mgraph))
end

"""
	abequalaba(v::Vertex, m::Monomialgraph)

Check if ab == aba, where v is the vertex corresponding to the monomial ab. To this end, let a = u_{i,j}. Then check whether either abu_{i,l} = 0 for all l≠j or abu_{l,j} = 0 for all l≠i. If so, ab = aba
"""
function abequalaba(v::Vertex, m::Monomialgraph)
	vcontent = v.monomial.content
	isum = true
	jsum = true
	for index in v.forwardedges
		w = getvertex(index, m)
		wcontent = w.monomial.content
		if vcontent != wcontent[1:2]
			# This is a forward edge cab rather than abc
			continue
		end
		if wcontent[3] == vcontent[1]
			if typeof(w.visited) != Nothing
				# aba = 0, assumed that ab≠0 so return false
				return false
			end
		elseif wcontent[3][1] == vcontent[1][1]
			if typeof(w.visited) == Nothing
				# w is a vertex with monomial u_{i,j}bu_{i,l} ≠ 0, l≠j.
				jsum = false
			end
		elseif wcontent[3][2] == vcontent[1][2]
			if typeof(w.visited) == Nothing
				# w is a vertex with monomial u_{i,j}bu_{l,j} ≠ 0, l≠i.
				isum = false
			end
		end
	end

	return (isum || jsum)
end

