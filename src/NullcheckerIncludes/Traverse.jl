"""
	traverse!(q::Vertexqueue, mgraph::Monomialgraph, verbose::Integer=0, outfile=0)

Traverse the monomialgraph mgraph given the starting zeros in q. The traversal will follow forward edges from vertices visited, backward edges from hypervertices visited and sideways edges from vertices visited using commutation relations. Here, forward edges are edges going from a vertex to one with a longer monomial, backward edges are edges going from a hypervertex (a set of vertices with same length monomials) to a vertex with a shorter monomial and sideways edges are edges going from a vertex to another one with a monomial of the same length. Verbosity levels from 0 (default, no outputs printed) to 3 (most detailed outputs printed) can be set. Outputs, if produced, may be written to a file if a filename is specified in outfile. If outfile is not a String, stdout will be used by default.
Returns the number of visited vertices.
"""
function traverse!(q::Vertexqueue, mgraph::Monomialgraph, verbose::Integer=0, outfile=0)
	n = mgraph.n
	k = mgraph.k
	visited = 0
	nstarting = length(q)
	nforward = 0
	nbackward = 0
	nsideway = 0
	nnewrule = 0

	# Output stream
	output = stdout
	if typeof(outfile) == String
		output = open(outfile, "a")
	end

	# Verbosity outputs
	if verbose != 0
		println(output, "Traverse Monomialgraph with n = ", n, " and k = ", k)
		println(output, "Number of starting zeros: ", nstarting)
	elseif verbose == 3
		println(output, "Traverse Monomialgraph")
		println(output, mgraph)
		println(output, "Starting queue:")
		println(output, q)
	end

	while !isempty(q)
		v = pop!(q)
		vi = computevertexindex(v, n)
		
		# Verbosity outputs
		if verbose >= 2
			println(output, "Visiting ", v)
		end
	
		visited = visited+1

		# Mark commutation of two generators
		if v.monomial.length == 2
			setcommute!(v.monomial, mgraph, true)
		end

		# Forward edges
		for index in v.forwardedges	
			w = getvertex(index, mgraph)
			if typeof(w.visited) == Nothing

				# Verbosity outputs
				if verbose == 3
					println(output, "Adding forward edge ", w)
				end

				w.visited = Forward(vi)
				push!(q, w)
				nforward = nforward+1
			end
		end

		# Commutation edges
		for l = 1:(v.monomial.length-1)
			if getcommute(Monomial([v.monomial.content[l], v.monomial.content[l+1]]), mgraph)
				w = getvertex(v.commutationedges[l], mgraph)
				if typeof(w.visited) == Nothing

					# Verbosity outputs
					if verbose == 3
						println(output, "Adding commutation edge ", w)
					end

					w.visited = Commutation(vi, l)
					push!(q, w)
					nsideway = nsideway+1
				end
			end
		end

		
		for (it, be) in v.backwardedges	
			hvnumzeros, hvnonvisited = checkhypervertex(vi, it[1], it[2], mgraph, 1)
			target = getvertex(be.target, mgraph)

			if (hvnumzeros == n) && typeof(target.visited) == Nothing
				# The entire hyperedge has been visited
				# Backward edge

				# Verbosity outputs
				if verbose == 3
					println(output, "Adding backwards edge ", target)
				end
	
				target.visited = Backward(vi, it[1], it[2])
				push!(q, target)
				nbackward = nbackward+1
			elseif (hvnumzeros == n-1) && typeof(target.visited) != Nothing
				# One vertex of the hypervertex is nonzero
				# Sum Rule

				w = getvertex(hvnonvisited[1], mgraph)
				w.visited = SumLong(vi, it[1], it[2])
				push!(q, w)
				nnewrule = nnewrule+1
				
				# Verbosity outputs
				if verbose == 3
					println(output, "Adding vertex ", w, " due to sum rule")
				end
			end
		end

		# Sum edges when checking from the shorter monomial
		for (it, se) in v.sumedges
			hvnumzeros, hvnonvisited = checkhypervertex(se, it[1], it[2], mgraph, 1)
			if (hvnumzeros == n-1)
				w = getvertex(hvnonvisited[1], mgraph)
				w.visited = SumShort(vi, it[1], it[2])
				push!(q,w)
				nnewrule = nnewrule+1

				# Verbosity outputs
				if verbose == 3
					println(output, "Adding vertex ", w, " due to sum rule")
				end
			end
		end
	end

	if verbose != 0
		println(output, "Visited $(visited) vertices of which")
		println(output, "    $(nstarting) are starting zeros,")
		println(output, "    $(nforward) were added through forward edges,")
		println(output, "    $(nbackward) were added through backward edges,")
		println(output, "    $(nsideway) were added through commutation edges and")
		println(output, "    $(nnewrule) were added due to sum rule.")
	end

	if typeof(output) == IOStream
		close(output)
	end
	
	return visited
end

