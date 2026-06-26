"""
	createcompletegraph(n::Integer)

Create the adjacency matrix for the complete graph on n points.
"""
createcompletegraph(n::Integer) = ones(Bool, n, n) - Matrix{Bool}(I, n, n)

"""
	createsixpointgraph()

Create the adjacency matrix for the smallest non-trivial graph with trivial automorphism group.
"""
createsixpointgraph() = [0 1 0 0 0 0;
						 1 0 1 1 0 0;
						 0 1 0 1 0 0;
						 0 1 1 0 1 0;
						 0 0 0 1 0 1;
						 0 0 0 0 1 0]

"""
	createpetersengraph()

Create the adjacency matrix for the Petersen graph.
"""
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

"""
	createfourpointgraphs()

Creates a list of the adjacency matrices of all graphs on 4 vertices with their names. Each entry is a tuple of the form (Name, AdjacencyMatrix). The graphs in question are in order
	(1) diamond
	(2) K4^c
	(3) K4
	(4) diamond^c
	(5) K1,3 (claw)
	(6) C4^c = K2,2^c
	(7) P4 (4-line)
	(8) K1,3^c (claw^c)
	(9) paw
	(10) C4 = K2,2
	(11) paw^c (3-line)
"""
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

