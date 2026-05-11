using DelimitedFiles
using TimerOutputs
using Oscar

include("nullchecker.jl")

function force_compile()
	g = createfourpointgraphs()[1][2]
	k = createk1monomialgraph(g)
	m = createnextmonomialgraph(k)
	q = Vertex[]
	temp_file = "temp.log"
	iteratetraverse!(q, m, g, 1, temp_file)
	rm(temp_file)
end

function write_header(graphpath, g, n, k, name=""; include_k=false)
	logname = graphpath * "/out"
	ctname = graphpath * "/commutation_table"
	if include_k == true
		logname = logname * "_" * string(k)
		ctname = ctname * "_" * string(k)
	end
	logname = logname * ".log"
	ctname = ctname * ".log"

	graphstring = split(repr("text/plain", g), "\n", limit=2)[2]
	f = open(logname, "a")
	if name != ""
		print(f, "Graph: ", name, "\n\n")
	end
	println(f, "Adjacency Matrix:\n", graphstring, "\n\nParameters: n = ", n, ", k = ", k)
	close(f)

	return logname, ctname
end

function write_results(logname, ctname, m)
	writedlm(ctname, Int.(m.commutationtable), ',')
	f = open(logname, "a")
	println(f, "\n#----------TIME-----------#\n")
	print_timer(f, sortby=:firstexec)
	println(f, "\n\n#----------ZEROS----------#\n")
	for i = 1:m.k
		println(f, "length ", i)
		for v in m.vertexlist[i]
			if typeof(v.visited) != Nothing
				vi = computevertexindex(v, m.n)
				println(f, "VI($(vi.k),$(vi.i)): ", v, " due to ", v.visited)
			end
		end
	end
	close(f)
end

function write_results_graph(logname, ctname, m, k, j)
	writedlm(ctname, Int.(m.commutationtable), ',')
	f = open(logname, "a")
	println(f, "\n#----------TIME-----------#\n")
	graphs_timer = TimerOutputs.get_defaulttimer()["monomial graphs ("*string(k)*")"]
	traverse_timer = TimerOutputs.get_defaulttimer()["gt j = "*string(j)*", k = "*string(k)]
	print_timer(f, merge(graphs_timer, traverse_timer), sortby=:firstexec)
	println(f, "\n\n#----------ZEROS----------#\n")
	for i = 1:m.k
		println(f, "length ", i)
		for v in m.vertexlist[i]
			if typeof(v.visited) != Nothing
				vi = computevertexindex(v, m.n)
				println(f, "VI($(vi.k),$(vi.i)): ", v, " due to ", v.visited)	
			end
		end
	end
	close(f)
end

"""
returns a vector [A, B, C, D] where
A = number of graphs with disjoint automorphisms
B = number of graphs with disjoint automorphisms and nullchecker output true
C = number of graphs without disjoint automorphisms and nullchecker output false
D = number of graphs with nullchecker output false and automorphism group z2

A is the number of graphs we know to have quantum symmetries
B is the number of errors!
C is the number of "interesting" graphs (i.e. those of which we neither know that they have no quantum symmetries nor that they don't have quantum symmetries)
D is the number of graphs with automorphism group z2 and of which we don't whether they have quantum symmetries
"""
function get_automorphism_statistics(graph_table, automorphism_data)
	disjoint_automorphisms = automorphism_data[:, 3]
	z2 = 2 .== automorphism_data[:, 1]
	qaut_commutative = graph_table[:, 2]

	num_disj_auts = sum(disjoint_automorphisms)
	num_disj_auts_and_qaut_comm = sum(qaut_commutative .&& disjoint_automorphisms)
	num_no_disj_auts_and_no_qaut_comm = size(graph_table, 1) - sum(qaut_commutative .|| disjoint_automorphisms)
	num_z2_and_not_qaut_comm = sum(z2 .&& .!(qaut_commutative))

	return [num_disj_auts, num_disj_auts_and_qaut_comm, num_no_disj_auts_and_no_qaut_comm, num_z2_and_not_qaut_comm]
end

function write_results_graph_list(outfile, n, max_k, l, graph_table; graph_names=0, class=0, automorphism_data=0)
	# get statistics
	commutation_statistic = [sum(x -> x == i, graph_table[:,3]) for i in 0:max_k]
	num_com = sum(commutation_statistic[2:end])
	if typeof(automorphism_data) != Int64
		automorphism_statistics = get_automorphism_statistics(graph_table, automorphism_data)
	end

	# write data
	# Header
	f = open(outfile, "a")
	println(f, "Results of running iterated traversal on graphs of same size")
	if typeof(class) != Int64
		println(f, "Class of graphs: ", class)
	end
	println(f, "Parameters: n = ", n, ", max_k = ", max_k, ", number of graphs = ", l)

	# Overview
	println(f, "\n#--------OVERVIEW---------#\n")
	println(f, "(a) k, the maximal length of monomials considered. k = 0 means that Qaut has not been found to be commutative.\n(b) Number of Qaut found to be commutative at that k. If a Qaut is found to be commutative for k, it would also be found to be commutative for all k' >= k. However, it is only counted for k here.\n")
	writedlm(f, ["(a)" "(b)"; 0:max_k commutation_statistic])

	print(f, "\nFound\n\t", num_com, " graphs with commutative quantum automorphism group")
	if typeof(automorphism_data) != Int64
		println(f, ",\n\t", automorphism_statistics[1], " graphs with disjoint automorphisms,")
		println(f, "\t", automorphism_statistics[2], " graphs with disjoint automorphisms and commutative quantum automorphism group (ERRORS),")
		println(f, "\t", automorphism_statistics[3], " graphs without disjoint automorphisms and with unknown quantum automorphism group,")
		println(f, "\t", automorphism_statistics[4], " graphs with automorphism group Z2 and unknown quantum automorphism group.")
	else
		println(f, ".")
	end

	# Detail
	println(f, "\n#--------DETAILED---------#\n")
	println(f, "(1) Graph number\n(2) Qaut commutative\n(3) Min k, the minimal k such that Qaut is found to be commutative or 0 if Qaut is not found to be commutative.")
	if typeof(graph_names) != Int64
		println(f, "(4) Graph identifier")
		graph_table = [graph_table graph_names]
		if typeof(automorphism_data) != Int64
			println(f, "(5) Order of automorphism group\n(6) Automorphism group\n(7) Disjoint automorphisms\n")
			graph_table = [graph_table automorphism_data]
		end
	else
		println(f, "\n")
	end
	header = reshape(["(" * string(x) * ")" for x in 1:size(graph_table, 2)], 1, :)
	writedlm(f, [header; graph_table])
	
	# Time
	println(f, "\n#----------TIME-----------#\n")
	print_timer(f, sortby=:firstexec)
	close(f)
end

function create_graph_table(l::Integer)
	gt = Matrix{Any}(undef, (l, 3))
	gt[:,1] = 1:l
	gt[:,2] .= false
	gt[:,3] .= 0
	return gt
end

function reset_commutation_table(m::Monomialgraph)
	m.commutationtable .= false
	m.commutationtable[CartesianIndex.(axes(m.commutationtable,1),axes(m.commutationtable,2))] .= true
end

function read_graphs(filename)
	g6 = readlines(filename * ".g6")
	graphs = readlines(filename * ".graphs")
	splitby = graphs[1]
	n = parse(Int, splitby)
	starts = vcat(findall(==(splitby), graphs), length(graphs)+1)
	ranges = (:).(starts[1:end-1].+1, starts[2:end].-1)
	adjacency_strings = getindex.(Ref(graphs), ranges)
	adjacency_matrices = reshape.((x -> parse.(Int, x)).(collect.(join.(adjacency_strings))), n, n)
	return [(g6[i], adjacency_matrices[i]) for i in 1:length(g6)]
end

function get_automorphism_group(graph::Matrix{Int64})
	oscar_graph = graph_from_adjacency_matrix(Undirected, graph)
	return automorphism_group(oscar_graph)
end

function disjoint_automorphisms(u::PermGroupElem, v::PermGroupElem)
	return isdisjoint(moved_points(u), moved_points(v))
end

function graph_has_disjoint_automorphisms(aut::PermGroup)
	for f in aut
		if isone(f) == true
			continue
		end
		for g in aut
			if isone(g) == true
				continue
			end
			if disjoint_automorphisms(f, g) == true
				return true
			end
		end
	end
	return false
end

function graph_has_disjoint_automorphisms(graph::Matrix{Int64})
	aut = get_automorphism_group(graph)
	return graph_has_disjoint_automorphisms(aut)
end

function compute_automorphism_data(read_graphs_output)
	adjacency_matrices = last.(read_graphs_output)
	automorphism_groups = get_automorphism_group.(adjacency_matrices)
	return [order.(automorphism_groups) describe.(automorphism_groups) graph_has_disjoint_automorphisms.(automorphism_groups)]
end
