"""
	Vertexindex(k::Integer,i::Integer)

Type encoding the index of a monomial, consisting of k, the length of the monomial and i, the index in vertexlist[k].
"""
struct Vertexindex
	k::Integer
	i::Integer
end
