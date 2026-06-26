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

