module Nullchecker

using LinearAlgebra
using MaxPlus
using TimerOutputs

"""
	Vertexindex(k::Integer,i::Integer)

Type encoding the index of a monomial, consisting of k, the length of the monomial and i, the index in vertexlist[k].
"""
struct Vertexindex
	k::Integer
	i::Integer
end

"""
Type encoding the predecessor of a visited vertex in the monomial graph.
"""
abstract type Predecessor end

"""
Type encoding the predecessor of a vertex if it has been visited or nothing if it has not been visited in the monomial graph.
"""
Visited = Union{Predecessor, Nothing}

"""
	Startzero(rule::Symbol)

Type encoding that the reason a vertex has been visited is that it is a startzero. The reason that it is a startzero is encoded in the field rule that is a symbol. Currently included are the symbols :degree and :distance
"""
struct Startzero <: Predecessor
	rule::Symbol
end

Base.show(io::IO, ::MIME"text/plain", sz::Startzero) = print(io, "Startzero with rule ", sz.rule)
Base.show(io::IO, sz::Visited) = print(io, "Startzero with rule ", sz.rule)

"""
	Forward(pre::Vertexindex)

Type encoding that a vertex has been visited due to applying the forward rule from its predecessor pre.
"""
struct Forward <: Predecessor
	pre::Vertexindex
end

Base.show(io::IO, ::MIME"text/plain", fw::Forward) = print(io, "Forward edge from ", fw.pre)
Base.show(io::IO, fw::Forward) = print(io, "Forward edge from ", fw.pre)

"""
	Commutation(pre::Vertexindex, pos::Integer)

Type encoding that a vertex has been visited due to applying the commutation rule from its predecessor pre where the commutation happens at position pos.
"""
struct Commutation <: Predecessor
	pre::Vertexindex
	pos::Integer
end

Base.show(io::IO, ::MIME"text/plain", com::Commutation) = print(io, "Commutation edge from ", com.pre, ", pos=" , com.pos)
Base.show(io::IO, com::Commutation) = print(io, "Commutation edge from ", com.pre, ", pos=" , com.pos)

"""
	Backward(pre::Vertexindex, pos::Integer, index::Char)

Type encoding that a vertex has been visited due to applying the backward rule from its predecessor pre with parameters pos and index.
"""
struct Backward <: Predecessor
	pre::Vertexindex
	pos::Integer
	index::Char
end

Base.show(io::IO, ::MIME"text/plain", bw::Backward) = print(io, "Backward edge from ", bw.pre, ", pos=" , bw.pos, ", index=", bw.index)
Base.show(io::IO, bw::Backward) = print(io, "Backward edge from ", bw.pre, ", pos=" , bw.pos, ", index=", bw.index)

"""
	SumShort(pre::Vertexindex, pos::Integer, index::Char)

Type encoding that a vertex has been visited due to applying the sum rule from its predecessor pre, which is the short monomial, and with parameters pos and index.
"""
struct SumShort <: Predecessor
	pre::Vertexindex
	pos::Integer
	index::Char
end

Base.show(io::IO, ::MIME"text/plain", sum::SumShort) = print(io, "Sum edge from short monomial ", sum.pre, ", pos=" , sum.pos, ", index=", sum.index)
Base.show(io::IO, sum::SumShort) = print(io, "Sum edge from short monomial ", sum.pre, ", pos=" , sum.pos, ", index=", sum.index)

"""
	SumLong(pre::Vertexindex, pos::Integer, index::Char)

Type encoding that a vertex has been visited due to applying the sum rule from its predecessor pre, which is the long monomial, and with parameters pos and index.
"""
struct SumLong <: Predecessor
	pre::Vertexindex
	pos::Integer
	index::Char
end

Base.show(io::IO, ::MIME"text/plain", sum::SumLong) = print(io, "Sum edge from long monomial ", sum.pre, ", pos=" , sum.pos, ", index=", sum.index)
Base.show(io::IO, sum::SumLong) = print(io, "Sum edge from long monomial ", sum.pre, ", pos=" , sum.pos, ", index=", sum.index)

"""
	Backwardedge(position::Integer, ij::Char, target::Vertexindex, next::Vertexindex)

Type for backwardedges. Each backwardedge corresponds to a hyperedge where the range is given by target (a single Vertex encoded as a Vertexindex), and the source consists of all vertices ``vu_{i,j}w`` for fixed monomials ``v`` and ``w`` and either fixed ``i`` or ``j``. The field ``ij`` is ``'i'`` iff ``j`` is fixed (i.e. one iterates over ``i``) and ``'j'`` otherwise. Position determines the position of the changing ``u_{i,j}``, i.e. ``position = length(v) + 1``.
"""
struct Backwardedge
	position::Integer
	ij::Char
	target::Vertexindex
	next::Vertexindex
	Backwardedge(position, ij, target, next) = (ij == 'i' || ij == 'j') ? new(position, ij, target, next) : throw(ArgumentError("ij should be the character 'i' or 'j' but is $(ij)"))
end

"""
	Monomial(content)

Type for monomials. The content is a vector of integer pairs, where the content ``[(i_0, j_0), ... , (i_k, j_k)]`` corresponds to the monomial ``u_{i_0, j_0} ... u_{i_k, j_k}``. ``content`` is a ``Vector{Tuple{Integer, Integer}}``.
"""
struct Monomial
	length::Integer
	content::Vector{Tuple{Integer, Integer}}
	Monomial(content) = new(length(content), content)
end

Base.show(io::IO, ::MIME"text/plain", m::Monomial) = print(io, "Monomial ", [i for i in m.content])
Base.show(io::IO, m::Monomial) = print(io, [i for i in m.content])

"""
	Vertex(monomial)

Type for vertices in monomial graph. It is mutable to allow for easy marking during monomial graph traversal. The constructor only takes the monomial and creates a vertex without forward or backward edges. They need to be added in a later step. The vertices that are pointed to in the forwardedges and backwardedges fields are encoded through the Vertexindex corresponding to their position in the vertex list of the parent monomial graph and Backwardedge respectively. 

The field commutationedges is a list of length l-1 where l is the length of the associated monomial u. At position r, the entry is a Vertexindex pointing to the vertex associated to the monomial u^(C)_r.

The field sumedges of the monomial v associated to a monomial u_{i_1,j_1}...u_{i_l,j_l} is a dictionary linking a key (r,c,t):>(Integer,Char,Integer) to a Vertexindex U of the vertex u associated to a monomial u_{i_1,j_1}...u_{i_{r-1},j_{r-1}} u* u_{i_r,j_r}... with u* = u_{i*,t}, i* any value, if c='i' and u* = u_{t,j*}, j* any value, if c='j'. Thus, check_hypervertex(U,r,c,MG,1) iterates over the set B^(c)_{u,r}. If the first return value is n-1, the second one is a vertex v such that (S^(c)_{v,r},v) is a sum edge.
"""
mutable struct Vertex
	visited::Visited
	const monomial::Monomial
	const forwardedges::Vector{Vertexindex}
	const backwardedges::Dict{Tuple{Integer,Char}, Backwardedge}
	const sumedges::Dict{Tuple{Integer,Char,Integer}, Vertexindex}
	const commutationedges::Vector{Vertexindex}
	Vertex(monomial) = new(nothing, monomial, [], Dict(), Dict(), [])
end

function Base.show(io::IO, ::MIME"text/plain", v::Vertex)
	println(io, "Vertex with content ", v.monomial)
	typeof(v.visited) == Nothing ? println(io, "Not visited") : println(v.visited)
	println(io, "Forward edges:")
	for fw in v.forwardedges
		println(io, fw)
	end
	println(io, "Backward edges:")
	for bw in values(v.backwardedges)
		println(io, bw.target)
	end
	println(io, "Sum edges:")
	for sum in v.sumedges
		println(io, sum[1], " => ", sum[2])
	end
	println(io, "Commutation edges:")
	for com in v.commutationedges
		println(io, com)
	end
end

Base.show(io::IO, v::Vertex) = print(io, v.monomial)

"""
	Monomialgraph(vertexlist, commutationtable)

Type for monomial graphs. Access to the vertices is through accessing vertexlist. Note that it is a list of lists where the k-th list contains all vertices corresponding to monomials of length k. Within each list, the monomials appear in lexicographical order.
"""
struct Monomialgraph
	n::Integer
	k::Integer
	vertexlist::Vector{Vector{Vertex}}
	commutationtable::Matrix{Bool}
	Monomialgraph(vertexlist, commutationtable) = new(sqrt(length(vertexlist[1])), length(vertexlist), vertexlist, commutationtable)
end

"""
	getvertex(vi::Vertexindex, m::Monomialgraph)

Getter for the vertex in Monomialgraph m encoded by Vertexindex vi.
"""
function getvertex(vi::Vertexindex, m::Monomialgraph)
	return m.vertexlist[vi.k][vi.i]
end

"""
	createcommutationtable(n::Integer)

Create the default commutation table given the size of the input graph n.
"""
createcommutationtable(n::Integer) = Matrix{Bool}(I, n^2, n^2)

"""
	createvertexlist(n::Integer, k::Integer)

Create the empty vertex list for the size of the input graph n and maximal monomial length k.
"""
function createvertexlist(n::Integer, k::Integer)
	vertexlist = Vector{Vector{Vertex}}(undef, k)
	for i = 1:k
		vertexlist[i] = Vector{Vertex}(undef, n^(2*i))
	end
	vertexlist
end

"""
	createcompletegraph(n::Integer)

Create the adjacency matrix for the complete graph on n points.
"""
createcompletegraph(n::Integer) = ones(Bool, n, n) - Matrix{Bool}(I, n, n)

"""
	createk1monomialgraph(Γ)

Create the monomial graph for Γ and k = 1.
"""
function createk1monomialgraph(Γ)
	n = size(Γ,1)
	commutationtable = createcommutationtable(n)
	vertexlist = createvertexlist(n, 1)
	for i = 1:n
		for j = 1:n
			m = Monomial([(i,j)])
			vertexlist[1][computemonomialindex(m,n)] = Vertex(m)
		end
	end
	return Monomialgraph(vertexlist, commutationtable)
end

"""
	createnextmonomialgraph(G::Monomialgraph)

Create the (n, k+1)-monomial graph given as input the (n, k)-monomial graph.
"""
function createnextmonomialgraph(G::Monomialgraph)
	n, k = G.n, G.k
	newvertexlist = deepcopy(G.vertexlist)
	leaflist = newvertexlist[k]
	newleaflist = Vector{Vertex}(undef, n^(2*(k+1)))

	# New vertices and forward edges
	for v in leaflist
		for i = 1:n
			for j = 1:n
				newmonomial = Monomial(push!(copy(v.monomial.content), (i,j)))
				newvertex = Vertex(newmonomial)
				index = computemonomialindex(newmonomial, n)
				vindex = Vertexindex(k+1, index)
				try
					newvertex = newleaflist[index]
				catch
					newleaflist[index] = newvertex
				end
				push!(v.forwardedges, vindex)
			end
		end

		for i = 1:n
			for j = 1:n
				newmonomial = Monomial(pushfirst!(copy(v.monomial.content), (i,j)))
				newvertex = Vertex(newmonomial)
				index = computemonomialindex(newmonomial, n)
				vindex = Vertexindex(k+1, index)
				try
					newvertex = newleaflist[index]
				catch
					newleaflist[index] = newvertex
				end
				push!(v.forwardedges, vindex)
			end
		end
	end

	push!(newvertexlist, newleaflist)

	# Backward edges, sum edges and commutation edges
	for v in newleaflist
		for l = 1:(k+1)
			target = computetargetindex(v.monomial, l, n)
			nexti = computenextindex(v.monomial, 'i', l, n)
			nextj = computenextindex(v.monomial, 'j', l, n)

			vtarget = Vertexindex(k, target)
			vnexti = Vertexindex(k+1, nexti)
			vnextj = Vertexindex(k+1, nextj)
			backi = Backwardedge(l, 'i', vtarget, vnexti)
			backj = Backwardedge(l, 'j', vtarget, vnextj)

			v.backwardedges[(l, 'i')] =  backi
			v.backwardedges[(l, 'j')] =  backj
			
			leaflist[target].sumedges[(l,'i',v.monomial.content[l][2])] = vnexti
			leaflist[target].sumedges[(l,'j',v.monomial.content[l][1])] = vnextj
		end
	
		for r = 1:k
			vcom = computeswitchindex(v.monomial, r, n)
			push!(v.commutationedges, Vertexindex(k+1, vcom))
		end
	end

	return Monomialgraph(newvertexlist, copy(G.commutationtable))
end

"""
	createmonomialgraph(Γ, k)

Create the (n,l)-monomial graphs for n = #(vertices in Γ) and l = 1,...,k.
The return value is a vector graphs of length k such that the (n,l)-monomial graph is found as graphs[l].
"""
function createmonomialgraph(Γ, k)
	graphs = Vector{Monomialgraph}(undef, k)
	graphs[1] = @timeit "create (k1)" createk1monomialgraph(Γ)
	for i = 2:k
		graphs[i] = @timeit "create (k"*string(i)*")" createnextmonomialgraph(graphs[i-1])
	end
	return graphs
end

"""
	computemonomialindex(c::Vector{Tuple{Integer, Integer}}, n::Integer)

Compute the index of monomial content c given n.
"""
function computemonomialindex(c::Vector{<:Tuple{Integer, Integer}}, n::Integer)
	k = length(c)	
	t = 1

	for i = 1:k
		t = (t-1)*n^2 + ((c[i][1]-1)*n + c[i][2])
	end
	return t
end

"""
	computemonomialindex(m::Monomial, n::Integer)

Compute the index of monomial m given n.
"""
computemonomialindex(m::Monomial, n::Integer) = computemonomialindex(m.content, n)

"""
	computevertexindex(v::Vertex, n::Integer)

Compute the Vertexindex of vertex v, i.e. a tuple consisting of the length of v and its monomial index.
"""
function computevertexindex(v::Vertex, n::Integer)
	k = v.monomial.length
	i = computemonomialindex(v.monomial, n)
	return Vertexindex(k,i)
end

"""
	computetargetindex(m::Monomial, cut::Integer, n::Integer)

Compute the index of monomial m without the generator at position cut.
"""
computetargetindex(m::Monomial, cut::Integer, n::Integer) = computemonomialindex(m.content[1:end .!= cut], n)

"""
	computenextindex(m::Monomial, ij::Char, position::Integer, n::Integer)

Compute the index of monomial m' where m' is monomial m with the ij component of the generator at position incremented by 1 modulo n.
"""
function computenextindex(m::Monomial, ij::Char, position::Integer, n::Integer)
	if ij != 'i' && ij != 'j'
		throw(ArgumentError("ij should be the character 'i' or 'j' but is $(ij)"))
	end
	content = m.content
	length = m.length
	IJ = ij == 'i' ? 1 : 2
	currentIndex = computemonomialindex(m,n)

	return content[position][IJ] != n ? currentIndex + n^(2*(length - position) + (IJ % 2)) : currentIndex - (n - 1) * n^(2*(length - position) + (IJ % 2))
end

function computenextindex2(m::Monomial, ij::Char, position::Integer, n::Integer)
	content = m.content
	length = m.length
	if ij == 'i'
		computemonomialindex([i != position ? content[i] : (content[i][1] % n + 1, content[i][2]) for i in 1:length], n)
	elseif ij == 'j'
		computemonomialindex([i != position ? content[i] : (content[i][1], content[i][2] % n + 1) for i in 1:length], n)
	else
		throw(ArgumentError("ij should be the character 'i' or 'j' but is $(ij)"))
	end
end

"""
	computeswitchindex(m::Monomial, l::Integer, n::Integer)

Compute the index of monomial m' where m' is monomial m with the generators at positions l and l+1 switched.
"""
function computeswitchindex(m::Monomial, l::Integer, n::Integer)
	content = m.content
	length = m.length
	computemonomialindex([i == l ? content[l+1] : (i == l+1 ? content[l] : content[i]) for i in 1:length], n)
end

createsixpointgraph() = [0 1 0 0 0 0;
						 1 0 1 1 0 0;
						 0 1 0 1 0 0;
						 0 1 1 0 1 0;
						 0 0 0 1 0 1;
						 0 0 0 0 1 0]

createpetersengraph() = [0 1 0 0 1 1 0 0 0 0;
						 1 0 1 0 0 0 1 0 0 0;
						 0 1 0 1 0 0 0 1 0 0;
						 0 0 1 0 1 0 0 0 1 0;
						 1 0 0 1 0 0 0 0 0 1;
						 1 0 0 0 0 0 0 1 1 0;
						 0 1 0 0 0 0 0 0 1 1;
						 0 0 1 0 0 1 0 0 0 1;
						 0 0 0 1 0 1 1 0 0 0;
						 0 0 0 0 1 0 1 1 0 0]

createfourpointgraphs() = (
							("diamond",
							[0 0 1 1;
							0 0 1 1;
							1 1 0 1;
							1 1 1 0]),

							("K4^c", 
							[0 0 0 0;
							0 0 0 0;
							0 0 0 0;
							0 0 0 0]),

							("K4", 
							[0 1 1 1;
							1 0 1 1;
							1 1 0 1;
							1 1 1 0]),

							("diamond^c", 
							[0 0 0 0;
							0 0 0 0;
							0 0 0 1;
							0 0 1 0]),

							("K1,3 (claw)", 
							[0 0 0 1;
							0 0 0 1;
							0 0 0 1;
							1 1 1 0]),

							("C4^c = K2,2^c", 
							[0 1 0 0;
							1 0 0 0;
							0 0 0 1;
							0 0 1 0]),

							("P4 (4-line)", 
							[0 0 1 0;
							0 0 0 1;
							1 0 0 1;
							0 1 1 0]),

							("K1,3^c (claw^c)", 
							[0 0 0 0;
							0 0 1 1;
							0 1 0 1;
							0 1 1 0]),

							("paw", 
							[0 0 0 1;
							0 0 1 1;
							0 1 0 1;
							1 1 1 0]),

							("C4 = K2,2", 
							[0 1 1 0;
							1 0 0 1;
							1 0 0 1;
							0 1 1 0]),

							("paw^c (3-line)", 
							[0 0 0 0;
							0 0 0 1;
							0 0 0 1;
							0 1 1 0])
						  )

"""
	checkhypervertex(vertexindex::Vertexindex, position::Integer, ij::Char, m::Monomialgraph, maxnonvisited::Integer=-1)

Check how many vertices of a hypervertex have been visited. The hypervertex is given by the vertices that share the same content except at position. At position, they all share either the i or j component of the generator. Thus, there are n vertices to each hypervertex. Thus, the hypervertex is determined by the arguments vertexindex, position and ij. 
The optional parameter maxnonvisited allows for an early stop of the check that immediately returns the current number of visited vertices as soon as the number of there have been more non-visited vertices than the parameter allows. If unset the entire hypervertex is checked. 
The return type is a pair (visited, nvvertices) where visited is the number of visited vertices (until maxnonvisited is exceeded) and nvvertices is an array of Vertexindices denoting the nonvisited vertices.
"""
function checkhypervertex(vertexindex::Vertexindex, position::Integer, ij::Char, m::Monomialgraph, maxnonvisited::Integer=-1)
	visited = 0
	nonvisited = 0
	nvvertices = Vertexindex[]
	stop = maxnonvisited < 0 ? m.n : maxnonvisited

	current = getvertex(vertexindex, m)

	for i = 1:m.n
		if typeof(current.visited) != Nothing
			visited = visited + 1
		else
			nonvisited = nonvisited + 1
			push!(nvvertices, computevertexindex(current, m.n))
		end
		if nonvisited > stop
			break
		end
		current = getvertex(current.backwardedges[(position, ij)].next, m)
	end

	return (visited, nvvertices)
end

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
	
	#=
	for i = 1:n, j = 1:n
		# Rules (1) and (2)
		for k = 1:n, l = 1:n
			if k == i || j == l
				# Skip cases that are 0 by rules (3) or (4)
				continue
			end

			if isone(g[i,k]) != isone(g[j,l])
				w = m.vertexlist[2][computemonomialindex([(i,j),(k,l)], n)]
				w.visited = true
				push!(q, w)
			end
		end

		# Rules (3) and (4)
		for k = 1:n
			if k != j
				w = m.vertexlist[2][computemonomialindex([(i,j),(i,k)], n)]
				w.visited = true
				push!(q, w)
			end
			if i != k
				w = m.vertexlist[2][computemonomialindex([(i,j),(k,j)], n)]
				w.visited = true
				push!(q, w)
			end
		end
	end
	=#

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

Create a function that takes a vertex of corresponding to a monomial m = ab of length 2 and returns false if the vertex has been visited or if the generators a and b are known to commute in mgraph.
"""
function generatepostfilter(mgraph::Monomialgraph)
	x -> !(typeof(x.visited) != Nothing || getcommute(x.monomial, mgraph))
end

"""
	getcommute(m::Monomial, mgraph::Monomialgraph)

Check if the monomial m = ab of length 2 is made up of generators a and b that are known to commute in mgraph.
"""
function getcommute(m::Monomial, mgraph::Monomialgraph)
	i1 = computemonomialindex([m.content[1]], mgraph.n)
	i2 = computemonomialindex([m.content[2]], mgraph.n)
	return mgraph.commutationtable[i1,i2]
end

"""
	setcommute!(m::Monomial, mgraph::Monomialgraph, commute::Bool)

Set the commutation relation of generators a and b, where m = ab, to commute in mgraph.
"""
function setcommute!(m::Monomial, mgraph::Monomialgraph, commute::Bool)
	i1 = computemonomialindex([m.content[1]], mgraph.n)
	i2 = computemonomialindex([m.content[2]], mgraph.n)
	mgraph.commutationtable[i1,i2] = mgraph.commutationtable[i2,i1] = commute
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

end
