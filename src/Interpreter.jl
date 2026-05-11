include("nullchecker.jl")

abstract type Line end

struct Startzeroline <: Line
	vi::Vertexindex
	mon::Monomial
end

struct Forwardline <: Line
	vi::Vertexindex
	mon::Monomial
	from::Vertexindex
end

struct Commutationline <: Line
	vi::Vertexindex
	mon::Monomial
	from::Vertexindex
	pos::Integer
end

struct Backwardline <: Line
	vi::Vertexindex
	mon::Monomial
	from::Vertexindex
	pos::Integer
	index::Symbol
end

struct Sumshortline <: Line
	vi::Vertexindex
	mon::Monomial
	from::Vertexindex
	pos::Integer
	index::Symbol
end

struct Sumlongline <: Line
	vi::Vertexindex
	mon::Monomial
	from::Vertexindex
	pos::Integer
	index::Symbol
end

function parseline(line)
	vistring, reststring = split(line, ": ")

	# parse vistring
	k, i = parse.(Int, split(vistring, r"[\(\,\)]")[2:3])
	vi = Vertexindex(k, i)

	# parse restring
	# (1) parse monomial
	mon, reason = split(reststring, " due to ")
	mon = eval(Meta.parse("Monomial($mon)"))
	# (2) parse reason
	if contains(reason, "Startzero")
		return Startzeroline(vi, mon)
	else
		type, params = split(reason, "Vertexindex")

		if occursin("Forward", type)
			l, j = parse.(Int, split(params, r"\(|, |\)", keepempty=false))
			from = Vertexindex(l,j)
			return Forwardline(vi, mon, from)
		else
			l, j, params = split(params, r"\), |[\(\, ]", keepempty=false, limit=3)
			l, j = parse.(Int, [l,j])
			from = Vertexindex(l,j)
			params = split(params, r"[=, ]", keepempty=false)
			pos = parse(Int, params[2])

			if occursin("Commutation", type)
				return Commutationline(vi, mon, from, pos)
			else
				index = Meta.parse(params[4])
				if occursin("Backward", type)	
					return Backwardline(vi, mon, from, pos, index)
				elseif occursin("Sum edge from short", type)
					return Sumshortline(vi, mon, from, pos, index)
				elseif occursin("Sum edge from long", type)
					return Sumlongline(vi, mon, from, pos, index)
				end	
			end
		end
	end
end

function loadzeros(filename)
	lines = readlines(filename)
	lines = filter(s -> occursin("VI", s), lines)
	lineobjects = parseline.(lines)
	return Dict(Pair.([l.vi for l in lineobjects], lineobjects))
end

function gethypervertex(line::Line, zerodict, n::Integer)
	hypervertex = Vector{Vertexindex}(undef, n)
	index = line.index == :i ? 'i' : 'j'
	position = line.pos
	currentvi = 0

	if typeof(line) == Backwardline
		currentvi = line.from
		hypervertex[1] = currentvi
	elseif typeof(line) == Sumshortline
		currentvi = line.vi
		hypervertex[1] = line.from	
	elseif typeof(line) == Sumlongline
		currentvi = line.vi
		monomial = line.mon
		hypervertex[1] = Vertexindex(currentvi.k - 1, computetargetindex(monomial, position, n))
	end

	k = currentvi.k
	for i = 2:n
		monomial = zerodict[currentvi].mon
		currentvi = Vertexindex(k, computenextindex(monomial, index, position, n))
		hypervertex[i] = currentvi
	end
	return hypervertex
end

function fullbackward(line::Line, zerodict, n::Integer)
	return gethypervertex(line, zerodict, n)
end

function fullsumshort(line::Line, zerodict, n::Integer)
	return gethypervertex(line, zerodict, n)
end

function fullsumlong(line::Line, zerodict, n::Integer)
	return gethypervertex(line, zerodict, n)
end

abstract type Node end

mutable struct Midnode <: Node
	const line::Line
	const children::Vector{Node}
	index::Integer
	Midnode(line, children) = new(line, children, 0)
end

function Base.show(io::IO, ::MIME"text/plain", mn::Midnode) 
	line = mn.line
	println(io, "$(line.vi): $(line.mon)")
	if typeof(line) == Forwardline
		print(io, "Forward edge, ")
	elseif typeof(line) == Commutationline
		print(io, "Commutation edge, ")
	elseif typeof(line) == Backwardline
		print(io, "Backward edge, ")
	elseif typeof(line) == Sumshortline || typeof(line) == Sumlongline
		print(io, "Sum edge, ")
	end
	println(io, "$(length(mn.children)) predecessors:")
	for child in mn.children
		Base.show(io, "text/plain", child)
	end
end

Base.show(io::IO, mn::Midnode) = print(io, "$(mn.line.vi): $(mn.line.mon)")

mutable struct Leaf <: Node
	const line::Line
	index::Integer
	Leaf(line) = new(line, 0)
end

Base.show(io::IO, ::MIME"text/plain", l::Leaf) = println(io, "Leaf $(l.line.vi): $(l.line.mon)")
Base.show(io::IO, l::Leaf) = print(io, "Leaf $(l.line.vi): $(l.line.mon)")

function propagationtree(vi::Vertexindex, zerodict, n::Integer, maxdepth::Integer)
	line = zerodict[vi]

	if maxdepth == 0 
		return Midnode(line, [])
	end

	if typeof(line) == Startzeroline
		return Leaf(line)
	else
		if typeof(line) == Forwardline
			childnodes = [propagationtree(line.from, zerodict, n, maxdepth-1)]
			return Midnode(line, childnodes)
		elseif typeof(line) == Commutationline
			childnodes = [propagationtree(line.from, zerodict, n, maxdepth-1)]
			return Midnode(line, childnodes)
		elseif typeof(line) == Backwardline
			hypervertex = fullbackward(line, zerodict, n)
			childnodes = propagationtree.(hypervertex, Ref(zerodict), Ref(n), Ref(maxdepth-1))
			return Midnode(line, childnodes)
		elseif typeof(line) == Sumshortline
			hypervertex = fullsumshort(line, zerodict, n)
			childnodes = propagationtree.(hypervertex, Ref(zerodict), Ref(n), Ref(maxdepth-1))
			return Midnode(line, childnodes)
		elseif typeof(line) == Sumlongline
			hypervertex = fullsumlong(line, zerodict, n)
			childnodes = propagationtree.(hypervertex, Ref(zerodict), Ref(n), Ref(maxdepth-1))
			return Midnode(line, childnodes)
		end
	end
end

function enumeratetree(node::Node, current::Integer)
	node.index = current
	current = current+1
	if typeof(node) != Leaf
		for child in node.children
			current = enumeratetree(child, current)
		end
	end
	return current
end

function prettyprintplain(io::IO, r::Node, depth::Integer, star::Bool=false)
	line = r.line
	indent = "| "^depth
	print(io, indent)
	if star == true
		print(io, "*")
	end
	print(io, "$(line.vi): $(line.mon)")
	if typeof(line) == Startzeroline
		println(io, " Startzero")
	else
		println(io, "")
		if typeof(line) == Forwardline
			print(io, indent, "Forward edge, ")
		elseif typeof(line) == Commutationline
			print(io, indent, "Commutation edge, ")
		elseif typeof(line) == Backwardline
			print(io, indent, "Backward edge, ")
		elseif typeof(line) == Sumshortline || typeof(line) == Sumlongline
			print(io, indent, "Sum edge, ")
		end
		println(io, "$(length(r.children)) predecessors:")
		for child in r.children
			st = child.line.vi == line.from
			prettyprintplain(io, child, depth+1, st)
		end
	end
end

printtab(io, text, n) = println(io, "\t"^n, text)

typesetmonomial(m::Monomial) = join((t -> "u_{$(t[1]),$(t[2])}").(m.content))

function tikznodes(io::IO, r::Node, tab::Integer, color::String="black", parentvi::Vertexindex=Vertexindex(0, 0))
	line = r.line
	fill = "white"
	if typeof(r) == Leaf
		fill = "green!20"
	elseif length(r.children) == 0
		fill = "red!20"
	end
	printtab(io, "\\node ($(r.index)) [color=$color, fill=$fill] {\$$(typesetmonomial(r.line.mon))\$};", tab)
	if typeof(r) != Leaf
		for child in r.children
			childcolor = line.from == child.line.vi ? "blue" : "black"
			tikznodes(io, child, tab, childcolor, line.vi)
		end
	end
end

function tikzedges(io::IO, r::Node, tab::Integer, parentvi::Vertexindex=Vertexindex(0, 0))
	if typeof(r) != Leaf
		vi = r.line.vi
		path = "\\path[<-] ($(r.index))"
		for child in r.children
			path = path * " edge ($(child.index))"
		end
		path = path * ";"
		printtab(io, path, tab)

		for child in r.children
			tikzedges(io, child, tab, vi)
		end
	end
end

function prettyprintlatex(io::IO, r::Node)
	tab = 0
	printtab(io, "\\documentclass[tikz, border=1cm]{standalone}", tab)
	printtab(io, "\\usetikzlibrary{graphdrawing}", tab)
	printtab(io, "\\usegdlibrary{trees}", tab)
	printtab(io, "\\begin{document}", tab)

	tab = 1
	printtab(io, "\\begin{tikzpicture}[tree layout, nodes={draw, rounded corners}]", tab)

	tab = 2
	tikznodes(io, r, tab)
	printtab(io, "", tab)
	tikzedges(io, r, tab)

	tab =1
	printtab(io, "\\end{tikzpicture}", tab)

	tab = 0
	printtab(io, "\\end{document}", tab)
end

prettyprint(io::IO, r::Node, ::MIME"text/plain") = prettyprintplain(io, r, 0)

prettyprint(io::IO, r::Node, ::MIME"text/latex") = prettyprintlatex(io, r)
