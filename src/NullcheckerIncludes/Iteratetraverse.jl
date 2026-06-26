"""
	resetvisit!(m::Monomialgraph)

Set all vertices in the monomialgraph as unvisited.
"""
function resetvisit!(m::Monomialgraph)
	for samelengthlist in m.vertexlist
		for v in samelengthlist
			v.visited = nothing
		end
	end
end

"""
	iteratetraverse!(q::Vertexqueue, mgraph::Monomialgraph, g::Matrix{<:Integer}, verbose::Integer=0, outfile=0, graphname="")

Iterate the traversal of the monomialgraph mgraph and the post processing until the commutation table and number of visited vertices does not change anymore, i.e. additional iterations do not add any new commutation relations or find new zeros. Outputs, if produced, may be written to a file if a filename is specified in outfile. If outfile is not a String, stdout will be used by default. Optionally, a name for the graph may be specified. This is only used for the verbosity outputs.

Returns true if the quantum automorphism group has been found to be commutative.
"""
function iteratetraverse!(q::Vertexqueue, mgraph::Monomialgraph, g::Matrix{<:Integer}, verbose::Integer=0, outfile=0; graphname::String="")
	commutative_flag = false
	output = stdout
	if typeof(outfile) == String
		output = open(outfile, "a")
	end

	if verbose != 0
		println(output)
		println(output, "----------------------------")
		print(output, "Start iterated traversal")
		if graphname != ""
			print(output, " of graph ", graphname)
		end
		println(output)
	end

	if typeof(output) == IOStream
		close(output)
	end

	it = 1
	maxcommutations = length(mgraph.commutationtable)
	commutationtables = [copy(mgraph.commutationtable)]
	@timeit "start zeros" startzeros!(q, mgraph, g)
	zeros = [length(q)]


	while it == 1 || zeros[it] - zeros[it - 1] != 0 || sum(xor.(commutationtables[it], commutationtables[it-1])) != 0

		if verbose != 0
			if typeof(outfile) == String
				output = open(outfile, "a")
			end
			print(output, "There are currently ", zeros[it], " zeros and ")
			println(output, sum(commutationtables[it]), " commutations known.")
			println(output, "")
			println(output, "Iteration ", it, ":")
			if typeof(output) == IOStream
				close(output)
			end
		end
		
		if it != 1
			resetvisit!(mgraph)
			startzeros!(q, mgraph, g)	
		end

		@timeit "traversal ("*string(it)*")" push!(zeros, traverse!(q, mgraph, verbose, outfile)) # traverse! is called
		newcommutations = sum(xor.(commutationtables[it], mgraph.commutationtable))
		@timeit "post processing ("*string(it)*")" postprocessing!(mgraph)
		push!(commutationtables, copy(mgraph.commutationtable))

		it = it + 1

		if verbose != 0
			if typeof(outfile) == String
				output = open(outfile, "a")
			end
			newcommutationspost = sum(xor.(commutationtables[it], commutationtables[it - 1]))
			print(output, "Found ", zeros[it] - zeros[it - 1], " new zeros and ")
			println(output, newcommutationspost, " new commutations, of which ", newcommutationspost - newcommutations, " were added by post processing.")
			if typeof(output) == IOStream 
				close(output)
			end
		end

		if mgraph.commutationtable == ones(Bool, size(mgraph.commutationtable))
			if verbose != 0
				if typeof(outfile) == String
					output = open(outfile, "a")
				end
				println(output, "The quantum automorphism group is commutative! Stopping iteration.")
				if typeof(output) == IOStream
					close(output)
				end
			end
			commutative_flag = true
			break
		end
	end

	if verbose != 0
		if typeof(outfile) == String
			output = open(outfile, "a")
		end
		println(output, "Finished iterated traversal after ", it-1, " iterations, finding ", zeros[it], " zeros and ", sum(commutationtables[it]), " commutations.")
		if commutative_flag == true
			println(output)
			println(output, "The quantum automorphism group is commutative!")
		end
		if typeof(output) == IOStream
			close(output)
		end
	end

	return commutative_flag
end

