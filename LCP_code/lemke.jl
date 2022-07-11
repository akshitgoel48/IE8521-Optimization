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
        ┌─────────────────────────────────────────────────────┐
        │ w_1 │ w_2 │...│ w_n │ z_1 │ z_2 │...│ z_n │ z_0 │ q │
        ├─────────────────────────────────────────────────────┤
  arr = │  :  │  :  │   │  :  │  :  │  :  │   │  :  │  :  │ : │
        │ B_1 │ B_2 │...│ B_n │ N_1 │ N_2 │...│ N_n │ d_0 │ q │
        │  :  │  :  │   │  :  │  :  │  :  │   │  :  │  :  │ : │
        └─────────────────────────────────────────────────────┘
  =#
end
function widx(tm::TableauMatrix);  n = size(tm); return collect(1:n); end
function zidx(tm::TableauMatrix);  n = size(tm); return collect((n+1):2n); end
function z0idx(tm::TableauMatrix); n = size(tm); return 2n+1; end
function qidx(tm::TableauMatrix);  n = size(tm); return 2n+2; end

"""
pivot structure
"""
mutable struct Pivot
  l    ::  Int64  # leaving variable (index in array)
  e    ::  Int64  # entering variable (index in array)
  c    ::  Int64  # complement of leaving variable (index in array)
end
function pivot_to_name(pivot::Pivot, n::Int64)
  labs = hcat(["w$i" for i in 1:n]..., ["z$i" for i in 1:n]..., "z0")
  enter = labs[pivot.e]
  leave = labs[pivot.l]
  return enter, leave
end

"""
get complement
"""
function get_complement(vidx::Int64, tm::TableauMatrix)
  n = size(tm)
  if vidx > 2n
    @info("asking for complement of z0")
    return 0
  else
    vcidx = mod1(vidx + n, 2n)
    return vcidx
  end
end

"""
convert LCP to tableau
"""
function TableauMatrix(lcp::LCP, d::Array{T,1}) where T <: Real
  n = size(lcp)
  B0 = collect(1:n)
  arr0 = hcat(I(n), -lcp.M, -d, lcp.q)
  return TableauMatrix(B0, arr0)
end
function Base.size(tm::TableauMatrix); return length(tm.B); end
function Base.print(tm::TableauMatrix)
  n = size(tm)
  labs = hcat(["w$i" for i in 1:n]..., ["z$i" for i in 1:n]..., "z0")
  header = hcat("basis", [(i ∈ tm.B ? "*w$i" : "w$i") for i in 1:n]..., [(n+i ∈ tm.B ? "*z$i" : "z$i")  for i in 1:n]..., (z0idx(tm) ∈ tm.B ? "*z0" : "z0"), "q")
  pretty_table(hcat(labs[tm.B], map(x -> (abs(x) >= tol ? @sprintf("%0.2f",x) : "⋅"), tm.arr)), header)
end

"""
get initial pivot
"""
function get_initial_pivot(tm::TableauMatrix)
  #=
  find leaving variable via mrt: q / M_{z0}
  =#
  l = tm.B[argmin(tm.arr[:,qidx(tm)] ./ -tm.arr[:,z0idx(tm)])]

  #=
  find entering variable
  =#
  e = z0idx(tm)

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
  λ = tm.arr[findfirst(pivot.l .== tm.B), pivot.e]
  r = tm.arr[findfirst(pivot.l .== tm.B),:]
  for i in 1:n
    if i == findfirst(pivot.l .== tm.B); continue; end
    c = tm.arr[i,pivot.e]
    tm.arr[i,:] .-= (c/λ) .* r
  end
  tm.arr[findfirst(pivot.l .== tm.B),:] .= r ./ λ
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
function process_lcp(lcp::LCP, verb::Bool=true)
  tm = TableauMatrix(lcp, ones(size(lcp)))
  keepgoing = true

  ## initial basis and pivot
  if verb; print(tm); end
  pivot = get_initial_pivot(tm); println(pivot)
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
    keepgoing = z0idx(tm) ∈ tm.B
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

# """
# bimatrix game structure
# """
# mutable struct BimatrixGame
#   A0 ::  AbstractArray{Float64,2}
#   B0 ::  AbstractArray{Float64,2}
# end
# function LCP(bmg::BimatrixGame)
#   mA, nA = size(A0)
#   mB, nB = size(B0)
#
#   ## make all entries stricly positive
#   a = minimum(A0)
#   b = minimum(B0)
#   if a <= 0; a = abs(a)+1; end
#   if b <= 0; b = abs(b)+1; end
#   A = A0 .+ a
#   B = B0 .+ b
#
#   ## build M and q
#   M = [ zeros(mA, mB)  A ; B'  zeros(nB, nA) ]
#   q = ones(mA + nB)
#   return LCP(-q, -M)
# end
# bmg = BimatrixGame(A0, B0)
# lcp = LCP(bmg)
# tm = process_lcp(lcp)