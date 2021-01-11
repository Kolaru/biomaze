### A Pluto.jl notebook ###
# v0.12.18

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ d597cbe0-512f-11eb-09d7-c337a535978e
begin
	using Parameters
	using LightGraphs
	using MetaGraphs
	using DelimitedFiles
	using ColorSchemes
	using PlutoUI
	using Plots
end

# ╔═╡ c4ae48a0-51ec-11eb-10d4-b37e6e3d9613
auxin_cmap = ColorSchemes.davos

# ╔═╡ 9d679400-5367-11eb-10db-4d7dd001f20b
pins_cmap = ColorSchemes.bilbao

# ╔═╡ 45ae1c00-51de-11eb-1177-152d8608f121
function show_maze(maze)
	heatmap(maze, aspect_ratio=1)
end

# ╔═╡ f22cb0e0-512f-11eb-2ede-bb187108b83e
begin  # Maze definition
	maze = readdlm("maze.csv", ',', Int)
	show_maze(maze)
end

# ╔═╡ 6c924060-5136-11eb-3c87-e1673c1aa5c4
begin
	is_in_grid(maze, I::CartesianIndex) = (I in CartesianIndices(maze))
	
	dirs = [
		CartesianIndex(1, 0),
		CartesianIndex(0, 1),
		CartesianIndex(-1, 0),
		CartesianIndex(0, -1)
	]
	
	function neighboring_indices(maze, I)
		indices = []
		for D in dirs
			J = I + D
			if is_in_grid(maze, J) && !is_wall(maze[J])
				push!(indices, J)
			end
		end
		
		return indices
	end
	
	is_wall(tile_id) = (tile_id == 0)  # Walls have tile id 0, others tile id -1
end

# ╔═╡ 670adc60-51a9-11eb-335c-9b76117d94de
begin  # Definition of cell indices
	
	cell_ids = zeros(Int, size(maze))  # Define cell id for each non wall cell
	cid = 1  # Current cell id
	positions = []  # Store position of each cell for plotting later
	for I in CartesianIndices(maze)
		if !is_wall(maze[I])
			cell_ids[I] = cid
			push!(positions, [Tuple(I)...])
			cid += 1
		end
	end
	
	show_maze(cell_ids)
end

# ╔═╡ 42b903b0-5135-11eb-2d31-411c4d8d886a
begin  # Definition of the undirected graph
	g = Graph()
	n_vertices = maximum(cell_ids)
	add_vertices!(g, n_vertices)
	v = 1  # Current vertex

	for I in CartesianIndices(maze)
		if !is_wall(maze[I])
			for J in neighboring_indices(maze, I)
				add_edge!(g, cell_ids[I], cell_ids[J])
			end
			v += 1
		end
	end
	
	g
end

# ╔═╡ 89581372-51c5-11eb-178e-3b96ab51a5cb
@with_kw struct CellParameters
	αa::Float64 = 0  # Auxin production
	βa::Float64 = 0  # Auxin degradation
	γD::Float64 = 1  # Auxin diffusion rate
	γA::Float64 = 1  # Auxin active transport rate
	αp::Float64 = 0  # PINS production
	βp::Float64 = 0  # PINS degradation
	V::Float64 = 1  # Cell volume
	μ::Float64 = 0.1  # PINS removal rate
	λ::Float64 = 0.1  # PINS insertion rate
end

# ╔═╡ 7f6c1b80-5135-11eb-3ced-07c6f5301780
begin  # Define graph with auxin and pins
	system = MetaDiGraph(g)
	
	set_prop!(system, :size, size(maze))
	
	for v in vertices(system)
		set_prop!(system, v, :auxin, 0.0)
		set_prop!(system, v, :pins, 1.0)
		set_prop!(system, v, :position, positions[v])
		set_prop!(system, v, :params, CellParameters())
	end
	
	# Set source
	source = 1
	set_prop!(system, source, :params, CellParameters(αa=50))
	set_prop!(system, source, :source, true)
	
	# Set sink
	sink = 32
	set_prop!(system, sink, :params, CellParameters(αa=0, βa=10))
	set_prop!(system, sink, :sink, true)
	
	# Edges are directed so we get difference between ij and ji
	for edge in edges(system)
		set_prop!(system, edge, :pins, 0.0)
		set_prop!(system, edge, :S, 1.0)
	end
	
	system
end

# ╔═╡ 915c6cb0-5373-11eb-3425-33201835c539
corners = [
	[+1, +1],
	[+1, -1],
	[-1, -1],
	[-1, +1]
]

# ╔═╡ 5ae66c70-51d0-11eb-2bce-61a3ca3dd780
function cell_corners(system, i, pad)
	pos = get_prop(system, i, :position)
	outer_corners = [corner .+ pos for corner in corners .* pad]
	return Shape(Tuple.(outer_corners))
end

# ╔═╡ 287fd910-5392-11eb-3389-d387f969a34d
function cell_interface(system, i, j, inner_pad, outer_pad)
	posi = get_prop(system, i, :position)
	posj = get_prop(system, j, :position)
	
	center = (posi + posj) ./ 2
	dv = 0.5.*(posj - posi)  # Vector from i to j
	v_perp = [dv[2], dv[1]]  # Perpendicular vector to dv
	
	corners = [center + v_perp, center - v_perp]

	out_corners = [(c - posi) .* 2outer_pad + posi for c in corners]
	in_corners = [(c - posi) .* 2inner_pad + posi for c in corners]
	
	return Shape(Tuple.([out_corners[1], out_corners[2], in_corners[2], in_corners[1]]))
end

# ╔═╡ bb9cf8d0-51d6-11eb-28f3-17459599a99f
function draw_cell(plt, system, i, ; inner_pad = 0.3, outer_pad = 0.45, max_auxin = 10.0, max_pins = 0.5)
	ai = get_prop(system, i, :auxin)
	c = 1 - ai / max_auxin
	inner_cell = cell_corners(system, i, inner_pad)
	plot!(plt, inner_cell, fillcolor=get(auxin_cmap, c), linewidth=0)
	
	for j in neighbors(system, i)
		interface = cell_interface(system, i, j, inner_pad, outer_pad)
		pij = get_prop(system, i, j, :pins)
		c = pij / max_pins
		plot!(plt, interface, fillcolor=get(pins_cmap, c), linewidth=0)
	end
	
	outer_cell = cell_corners(system, i, outer_pad)
	plot!(plt, outer_cell, fillcolor=nothing)
	
	x, y = get_prop(system, i, :position)
	if has_prop(system, i, :source)
		plot!(plt, annotations=(x, y, "A"))
	end
	if has_prop(system, i, :sink)
		plot!(plt, annotations=(x, y, "Ω"))
	end
end

# ╔═╡ 7e0e0b80-51cc-11eb-32a1-23e61e9e7703
function draw_system(system ; size=(1200, 1200))
	w, h = get_prop(system, :size)
	plt = plot(
		bg = :white,
		framestyle = :none,
		size = size,
		legend = false,
		xlim=(0, w + 1),
		ylim=(0, h + 1)
	)
	cells = draw_cell.(Ref(plt), Ref(system), vertices(system))
	return plt
end

# ╔═╡ 5f9fa750-5297-11eb-1203-31f0a9e09f31
draw_system(system)

# ╔═╡ ea5ff850-51c4-11eb-27bc-53f1c34ef0bd
function J(system, i, j)
	@unpack γD, γA = get_prop(system, i, :params)
	pij = get_prop(system, i, j, :pins)
	pji = get_prop(system, j, i, :pins)
	ai = get_prop(system, i, :auxin)
	aj = get_prop(system, j, :auxin)

	return γD * (aj - ai) + γA * (ai*pij - aj*pji)
end

# ╔═╡ c4205290-51c2-11eb-066b-3fbd0dfd0f49
function da(system, i)
	@unpack αa, βa, V = get_prop(system, i, :params)
	ai = get_prop(system, i, :auxin)
	transfer = 0.0
	
	for j in neighbors(system, i)
		Sij = get_prop(system, i, j, :S)
		
		transfer += Sij * J(system, i, j)
	end
	
	return αa - βa*ai + transfer/V
end

# ╔═╡ d07ac5d0-51cb-11eb-062b-65b7191a028f
function h(J)
	return J^2 / (1 + J^2)
end

# ╔═╡ 90cf32b0-51c4-11eb-0071-41de58516cb9
function dp(system, i, j)
	@unpack λ, μ = get_prop(system, i, :params)
	Pi = get_prop(system, i, :pins)
	pij = get_prop(system, i, j, :pins)
	return λ * Pi * h(J(system, i, j)) - μ * pij
end

# ╔═╡ 0bb4f2e0-51c4-11eb-011c-3daa9ef0d7e3
function dP(system, i)
	@unpack αp, βp, V = get_prop(system, i, :params)
	Pi = get_prop(system, i, :pins)
	
	transfer = 0.0
	
	for j in neighbors(system, i)
		Sij = get_prop(system, i, j, :S)
		
		transfer -= Sij * dp(system, i, j)
	end
	
	return αp - βp * Pi + transfer / V
end

# ╔═╡ 80af7c3e-51ca-11eb-222a-0b0c6a9c4551
dt = 0.1

# ╔═╡ 85ee1b30-51ca-11eb-3d97-2b90f8d7c249
begin  # Simulate
	n_frames = 10000
	save_each = 100
	
	systems = [system]
	frames = [1]
	
	for frame in 2:n_frames
		new_system = copy(system)

		for i in vertices(system)
			auxin = get_prop(system, i, :auxin) + dt * da(system, i)
			pins = get_prop(system, i, :pins) + dt * dP(system, i)
			
			set_prop!(new_system, i, :auxin, auxin)
			set_prop!(new_system, i, :pins, pins)
			
			for j in neighbors(system, i)
				pins = get_prop(system, i, j, :pins) + dt * dp(system, i, j)
				
				set_prop!(new_system, i, j, :pins, pins)
			end
		end
		
		if frame % save_each == 0
			push!(systems, new_system)
			push!(frames, frame)
		end
		system = new_system
	end
end			

# ╔═╡ a623e0f0-5207-11eb-3b96-b3f0d076fde8
md"Diplayed frame $(@bind F PlutoUI.Slider(1:length(systems), show_value=true))"

# ╔═╡ 9a130ac0-5207-11eb-0a5d-2ff8421fba98
draw_system(systems[F])

# ╔═╡ Cell order:
# ╠═d597cbe0-512f-11eb-09d7-c337a535978e
# ╟─c4ae48a0-51ec-11eb-10d4-b37e6e3d9613
# ╟─9d679400-5367-11eb-10db-4d7dd001f20b
# ╠═45ae1c00-51de-11eb-1177-152d8608f121
# ╠═f22cb0e0-512f-11eb-2ede-bb187108b83e
# ╠═6c924060-5136-11eb-3c87-e1673c1aa5c4
# ╠═670adc60-51a9-11eb-335c-9b76117d94de
# ╠═42b903b0-5135-11eb-2d31-411c4d8d886a
# ╠═89581372-51c5-11eb-178e-3b96ab51a5cb
# ╠═7f6c1b80-5135-11eb-3ced-07c6f5301780
# ╠═5f9fa750-5297-11eb-1203-31f0a9e09f31
# ╠═915c6cb0-5373-11eb-3425-33201835c539
# ╠═5ae66c70-51d0-11eb-2bce-61a3ca3dd780
# ╠═287fd910-5392-11eb-3389-d387f969a34d
# ╠═bb9cf8d0-51d6-11eb-28f3-17459599a99f
# ╠═7e0e0b80-51cc-11eb-32a1-23e61e9e7703
# ╠═c4205290-51c2-11eb-066b-3fbd0dfd0f49
# ╠═ea5ff850-51c4-11eb-27bc-53f1c34ef0bd
# ╠═0bb4f2e0-51c4-11eb-011c-3daa9ef0d7e3
# ╠═90cf32b0-51c4-11eb-0071-41de58516cb9
# ╠═d07ac5d0-51cb-11eb-062b-65b7191a028f
# ╠═80af7c3e-51ca-11eb-222a-0b0c6a9c4551
# ╠═85ee1b30-51ca-11eb-3d97-2b90f8d7c249
# ╟─a623e0f0-5207-11eb-3b96-b3f0d076fde8
# ╠═9a130ac0-5207-11eb-0a5d-2ff8421fba98
