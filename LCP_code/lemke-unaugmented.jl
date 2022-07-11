using LinearAlgebra, PrettyTables, Printf
global const tol = eps(1e10)

"""
lcp structure
"""
mutable struct LCP
  q  ::  AbstractArray{Float64,1}
  M  ::  AbstractArray{Float64,2}
end
function Base.size(lcp::LCP); return length(lcp.q) ; end

"""
tableau array structure
"""
mutable struct TableauMatrix
  B    ::  Array{Int64,1}    # basis variable column indices
  arr  ::  Array{Float64,2}  # concatenated vectors
  #=
        ┌───────────────────────────────────────────────┐
        │ w_1 │ w_2 │...│ w_n │ z_1 │ z_2 │...│ z_n │ q │
        ├───────────────────────────────────────────────┤
  arr = │  :  │  :  │   │  :  │  :  │  :  │   │  :  │ : │
        │ B_1 │ B_2 │...│ B_n │ N_1 │ N_2 │...│ N_n │ q │
        │  :  │  :  │   │  :  │  :  │  :  │   │  :  │ : │
        └───────────────────────────────────────────────┘
  =#
end
function widx(tm::TableauMatrix);  n = size(tm); return collect(1:n); end
function zidx(tm::TableauMatrix);  n = size(tm); return collect((n+1):2n); end
function qidx(tm::TableauMatrix);  n = size(tm); return 2n+1; end

"""
convert LCP to tableau
"""
function TableauMatrix(lcp::LCP)
  n = size(lcp)
  B0 = collect(1:n)
  arr0 = hcat(I(n), -lcp.M, lcp.q)
  return TableauMatrix(B0, arr0)
end
function Base.size(tm::TableauMatrix); return length(tm.B); end
function Base.print(tm::TableauMatrix)
  n = size(tm)
  labs = hcat(["w$i" for i in 1:n]..., ["z$i" for i in 1:n]...)
  header = hcat("basis", [(i ∈ tm.B ? "*w$i" : "w$i") for i in 1:n]..., [(n+i ∈ tm.B ? "*z$i" : "z$i")  for i in 1:n]..., "q")
  # pretty_table(hcat(labs[tm.B], tm.arr), header)
  pretty_table(hcat(labs[tm.B], map(x -> (abs(x) >= tol ? @sprintf("%0.2f",x) : "⋅"), tm.arr)), header)
end

"""
pivot structure
"""
mutable struct Pivot
  l    ::  Int64  # leaving variable (index in array)
  e    ::  Int64  # entering variable (index in array)
  c    ::  Int64  # complement of leaving variable (index in array)
end
function pivot_to_name(pivot::Pivot, n::Int64)
  labs = hcat(["w$i" for i in 1:n]..., ["z$i" for i in 1:n]...)
  enter = labs[pivot.e]
  leave = labs[pivot.l]
  return enter, leave
end

"""
get complement
"""
function get_complement(vidx::Int64, tm::TableauMatrix)
  n = size(tm)
  vcidx = mod1(vidx + n, 2n)
  return vcidx
end

"""
get initial pivot
"""
function get_initial_pivot(tm::TableauMatrix, s::Int64)
  n = size(tm)
  #=
  find leaving variable via mrt: q / M_{s}
  =#
  mrt = tm.arr[:,qidx(tm)] ./ -tm.arr[:,n+s]
  mrt[-tm.arr[:,n+s] .<= tol] .= Inf
  l = tm.B[argmin(mrt)]

  #=
  find entering variable
  =#
  e = n+s

  #=
  find next pivot entering variable
  =#
  c = get_complement(l, tm)
  return Pivot(l, e, c)
end

"""
get new pivot
"""
function update_pivot!(pivot::Pivot, tm::TableauMatrix)
  n = size(tm)
  #=
  find entering variable
  =#
  pivot.e = pivot.c

  #=
  find leaving variable
  =#
  mrt = tm.arr[:,qidx(tm)] ./ tm.arr[:,pivot.e]
  mrt[tm.arr[:,pivot.e] .<= tol] .= Inf
  if all(mrt .== Inf); pivot.l = -1; end
  pivot.l = tm.B[argmin(mrt)]

  #=
  find next pivot entering variable
  =#
  pivot.c = get_complement(pivot.l, tm)
  nothing
end

"""
update basis
"""
function update_basis!(tm::TableauMatrix, pivot::Pivot, verb::Bool=true)
  if verb
    enter, leave = pivot_to_name(pivot, size(tm))
    println("(entering, leaving): ($enter, $leave)")
  end
  splice!(tm.B, findfirst(pivot.l .== tm.B), pivot.e)
  nothing
end

"""
perform pivot elimination
"""
function pivot_eliminate!(tm::TableauMatrix, pivot::Pivot)
  n = size(tm)
  Blidx = findfirst(pivot.l .== tm.B)
  λ = tm.arr[Blidx, pivot.e]
  r = tm.arr[Blidx,:]
  for i in 1:n
    if i == Blidx; continue; end
    c = tm.arr[i,pivot.e]
    tm.arr[i,:] .-= (c/λ) .* r
  end
  tm.arr[Blidx,:] .= r ./ λ
  nothing
end

"""
solve lcp
  ┌────────────┐
  │ w - Mz = q │
  │      w ≧ 0 │
  │      z ≧ 0 │
  │  z ∘ w = 0 │
  └────────────┘
via complementary lemke
"""
function process_lcp(lcp::LCP, s::Int64, verb::Bool=true)
  tm = TableauMatrix(lcp)
  n = size(tm)
  keepgoing = true

  ## initial basis and pivot
  if verb; print(tm); end
  pivot = get_initial_pivot(tm, s)
  pivot_eliminate!(tm, pivot)
  update_basis!(tm, pivot)
  if verb; print(tm); end

  ## complementary pivots until z0 drops out of basis
  while keepgoing
    update_pivot!(pivot, tm)
    if pivot.l == -1; @warn("problem is infeasible"); return tm; end
    pivot_eliminate!(tm, pivot)
    update_basis!(tm, pivot, verb)
    if verb; print(tm); end
    wsidx = findfirst(tm.B .== s)
    zsidx = findfirst(tm.B .== s+n)
    if (&)((s+n ∈ tm.B), (s ∈ tm.B))
      if (abs(tm.arr[wsidx,qidx(tm)]) <= tol) || (abs(tm.arr[zsidx,qidx(tm)]) <= tol)
        if verb; println("both ws and zs are in basis but are complementary"); end
        keepgoing = false
      end
    else
      if verb; println("one of ws or zs has dropped from basis"); end
      keepgoing = false
    end
  end
  return tm
end

"""
recover optimal solution from tableau
"""
function recover_solution(tm::TableauMatrix, lcp::LCP)
  n = size(tm)
  sol = zeros(2n)
  sol[tm.B] .= tm.arr[:,qidx(tm)]
  w = sol[1:n]
  z = sol[(n+1):(2n)]
  @assert(norm(w - lcp.M*z - lcp.q) <= tol)
  @assert(prod(w .>= -tol))
  @assert(prod(z .>= -tol))
  @assert(prod(z .* w) .<= tol)
  return (w=w, z=z)
end
