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
