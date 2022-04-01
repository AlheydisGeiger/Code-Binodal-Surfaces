# This document contains functions to investigate binodal varieties for a fixed support of the surfaces given by a (small) lattice polytope

#Dependencies: the following file contains helper functions that allow the code below to be shorter.
include("helper_functions.jl")

#Given a 3-polytope by means of a matrix filled with the lattice points,
#this function computes the generalized binodal variety of those polynomials with the given support 
#INPUT: matrix(ZZ,n,3,[...])
#OUTPUT: returns the radical ideal generating the generalized binodal variety
function general_binodal(A::MatElem)
    n = nrows(A)
    R, c =  PolynomialRing(QQ,vcat(["c$i" for i in 1:n+7]))
    x = c[n+1]
    y = c[n+2]
    z = c[n+3]
    u = c[n+4]
    v = c[n+5]
    w = c[n+6]
    t = c[n+7]
    #(x,y,z) stand for the position of the first node, (u,v,w) for the position of the second node
    #t is an additional variable to make sure that the nodes are contained in the torus
    Id = identity_matrix(ZZ,n)
    M = hcat(Id,A)
    M = hcat(M, zero_matrix(ZZ, n, 4)) #add coefficients for u,v,w,t
    h = MPolyBuildCtx(R)
    for i in 1:n
      push_term!(h,1,Vector{Int64}([M[i,j] for j in 1:ncols(M)]))
    end
    p= finish(h)
    q = evaluate(p,[x,y,z],[u,v,w])
    I = ideal(R,[p,q,derivative(p,x),derivative(p,y),derivative(p,z),derivative(q,u),derivative(q,v),derivative(q,w),1-t*x*y*z*u*v*w])
    Sat = saturation(I, ideal(R,[x-u,y-v,z-w]))
    K = eliminate(Sat,[t,x,y,z,u,v,w])
    J = gens(K)
    RR, a = PolynomialRing(QQ, ["a$i" for i in 1:n])
    V = Vector{fmpq_mpoly}(undef,length(J))
    for i in 1:length(J)
        V[i] = evaluate(J[i],vcat(a,[zero(RR) for i in 1:7]))
    end
    II = ideal(RR,V)
    Radical = radical(II)
    return Radical
end

#Given a 3-polytope by means of a matrix filled with the n lattice points,
#this function computes the binodal variety of those polynomials with the given support 
#If the dimension of the generalized binodal variety is less than n-2 the function returns an error.
#INPUT: matrix(ZZ,n,3,[...])
#OUTPUT: returns the radical ideal generating the binodal variety
function binodal(A::MatElem) 
  n = nrows(A)
  Var = general_binodal(A);
  RR = base_ring(Var)
  monomes = RR(Rational{BigInt}[1],[fill(1,n)])# this is the monomial a[1]*...*a[n]
  Mon = ideal(RR, monomes)
  dimension = dim(Var)
  @assert dimension >= n-2 "Generalized binodal variety is of too small dimension."
  Bin = helper_component_check(Var)
  return Bin
end

#Given a 3-polytope by means of a matrix filled with the n lattice points,
#this function computes the binodal variety of those polynomials with the given support 
#If the dimension of the generalized binodal variety is less than n-2 the function returns an error.
#INPUT: matrix(ZZ,n,3,[...])
#OUTPUT: returns the radical ideal generating the binodal variety
function binodal(I::MPolyIdeal{fmpq_mpoly}) 
  RR = base_ring(I)
  n = nvars(RR);
  monomes = RR(Rational{BigInt}[1],[fill(1,n)])# this is the monomial a[1]*...*a[n]
  Mon = ideal(RR, monomes)
  dimension = dim(I)
  @assert dimension >= n-2 "Generalized binodal variety is of too small dimension."
  Bin = helper_component_check(I)
  return Bin
end

#Given an ideal this function computes the degree
#INPUT: ideal
#OUTPUT: the degree of the given ideal (an integer)
function compute_degree(Path::MPolyIdeal{fmpq_mpoly})
  Gens = gens(Path)
  k = nvars(base_ring(Path))
  Genshom = Vector(undef, length(Gens))
  for i in 1:length(Gens) 
      Genshom[i] = homogenization(Gens[i],"z",k+1)
  end
  RRhom = parent(Genshom[1]) # This is RR, homogenous with additional variable e3 and weights [1, 1,..., 1]
  IIhom = ideal(RRhom, RRhom.(Genshom))
  Q, _ = quo(RRhom, IIhom)
  deg = degree(quo(RRhom, IIhom)[1]) #this is the degree of the hyperplane section of the binodal variety in the left out points
  return deg
end



#Given an ideal ths function prints the affine dimension, degree and whether the variety is empty or not of expected dimension.
#INPUT: matrix(ZZ,n,3,[...])
#OUTPUT: Pair (dimension, degree)
function investigate_binodal(I::MPolyIdeal{fmpq_mpoly})
  deg = compute_degree(I)
  dimension = dim(I)
  n = nvars(base_ring(I))
  println("Affine dimension of the binodal variety is ", dimension)
  if isone(I)
    println("Binodal variety is empty")
  elseif dimension==n-2
    println("The binodal variety is of expected affine dimension.")
    println("Degree of the binodal variet is ", deg)
  else 
    println("The binodal variety is not of expected affine dimension.")
  end
  return dimension, deg
end

#Given a 3-polytope by means of a matrix filled with the lattice points,
#this function computes the binodal variety of those polynomials with the given support 
#and prints the affine dimension, degree and whether the binodal variety is empty or not of expected dimension.
#INPUT: matrix(ZZ,n,3,[...])
#OUTPUT: Pair (dimension, degree)
function investigate_binodal(A::MatElem)
  I = binodal(A)
  deg = compute_degree(I)
  dimension = dim(I)
  n = nvars(base_ring(I))
  println("Affine dimension of the generalized variety is ", dimension)
  if isone(I)
    println("Binodal variety is empty")
  elseif dimension==n-2
    println("The binodal variety is of expected affine dimension.")
    println("Degree of the binodal variet is ", deg)
    println(A, " could be a binodal polytope.")
  else 
    println("The binodal variety is not of expected affine dimension.")
  end
  return dimension, deg
end

# This function computes the intersection of the binodal variety for a 3-polytope given as a matrix filled with the lattice points, 
# with 4 linear spaces of codimension 1 given by substituting all variables apart from those with indices in p with the values given in v ("Path Ideal")
# INPUT: The n lattice points of the polytope as matrix, the positions of the two left out lattice points as a Vector, a generic point Vector of length n-2
# OUTPUT: the Path ideal
function con_mult(A::MatElem, p::Vector{Int}, v::Vector)
  @assert p[1]!= p[2] "Left out points must me different"
  @assert length(p)==2 "This function is only meant for lattice paths with 2 left out points."
  @assert p[1] < p[2] "Left out points must be ordered from small to large."
  I = binodal(A)
  J = gens(I);
  n = nrows(A)
  o = helper_find_indices(n,p)
  V = Vector{fmpq_mpoly}(undef,length(J))# we fill this vector according to the left out point conditions p[1] and p[2] with the values from v
  for i in 1:length(J)
    V[i] = evaluate(J[i],[o[i] for i in 1:(n-2)], [v[i] for i in 1:(n-2)])
  end
  RRR, e = PolynomialRing(QQ,["e$i" for i in 1:2])
  w = helper_switch_ring(n,p,[true,true],v,RRR)
  m = length(J)
  V2 = helper_finish_switch_ring(n,m,w,V)
  II = ideal(RRR,V2) # this is the ideal of the intersection of the binodal variety with the hyperplanes given by v
  Radical = radical(II)
  Path = helper_component_check(Radical)
  return Path
end

# This function computes the intersection of the binodal variety for a 3-polytope given as a matrix filled with the lattice points, 
# with 4 linear spaces of codimension 1 given by substituting all variables apart from those with indices in p with the values given in v ("Path Ideal")
# INPUT: the binodal variety, the positions of the two left out lattice points as a Vector, a generic point Vector of length n-2
# OUTPUT: the Path ideal
function con_mult(I::MPolyIdeal{fmpq_mpoly}, p::Vector{Int}, v::Vector)
  @assert p[1]!= p[2] "Left out points must me different"
  @assert length(p)==2 "This function is only meant for lattice paths with 2 left out points."
  @assert p[1] < p[2] "Left out points must be ordered from small to large."
  J = gens(I);
  n = nvars(base_ring(I))
  o = helper_find_indices(n,p)
  V = Vector{fmpq_mpoly}(undef,length(J))# we fill this vector according to the left out point conditions p[1] and p[2] with the values from v
  for i in 1:length(J)
    V[i] = evaluate(J[i],[o[i] for i in 1:(n-2)], [v[i] for i in 1:(n-2)])
  end
  RRR, e = PolynomialRing(QQ,["e$i" for i in 1:2])
  w = helper_switch_ring(n,p,[true,true],v,RRR)
  m = length(J)
  V2 = helper_finish_switch_ring(n,m,w,V)
  II = ideal(RRR,V2) # this is the ideal of the intersection of the binodal variety with the hyperplanes given by v
  Radical = radical(II)
  Path = helper_component_check(Radical)
  return Path
end


# b[1] = true means, that the point at position p[1] is a left out point, not part of a gap.
# For the moment, p and b are only allowed exactly 2 entries! It should be possible to expand the function if needed.
# p needs to be filled according to size, i.e., the smallest value first, the largest at the end.
# INPUT: The n lattice points of the polytope as matrix, the positions of the two left out lattice points as a Vector, a generic point Vector of length n-2
# OUTPUT: prints the lattice path multiplicity if the points in v are generic enough.
function path_mult(A::MatElem, p::Vector{Int}, b::Vector{Bool}, v::Vector)
  n = nrows(A)
  @assert p[1]!= p[2] "Left out points must me different"
  @assert length(p) ==2 "This function is only meant for lattice paths with 2 left out points."
  @assert length(p) ==length(b)  "$p and $b need to have the same length"
  @assert length(v)== n-length(p) "Number of point conditions in $v is wrong."
  @assert p[1] < p[2] "Left out points must be ordered from small to large."
  @assert b[1] || b[2] "This function is only meant for lattice paths with at most one gap."
  @assert (!(p[1] == 1 || p[1] == 2) || !b[1] || b[2] || p[2]>=4) "Invalid Path given." #This path is invalid because: 
  @assert (b[1] || !(p[1]<3||p[1]==n-1)) "Invalid Path given." #This path is invalid because: 
  I = binodal(A)
  J = gens(I)
  if b[1] && b[2] #case of the connected paths (no gap in the path)
    Path = con_mult(A,p,v)
  else #disconnected multiplicty
    if b[1] && !b[2] #first left out point, then gap
      o = helper_find_indices(p[2],p)
      V = Vector{fmpq_mpoly}(undef,length(J))
      for i in 1:length(J)
          V[i] = evaluate(J[i],[o[i] for i in 1:(p[2]-2)], [v[i] for i in 1:(p[2]-2)])
        end
        println(V)
        RRR, e = PolynomialRing(QQ,["e$i" for i in 1:2])
        w = helper_switch_ring(n,p,b,v,RRR)
        V2 = helper_finish_switch_ring(n,length(J),w,V)
    elseif !b[1] && b[2] # here the gap appears first
      o = helper_find_indices(p[1],p)
      V = Vector{fmpq_mpoly}(undef,length(J))
      for i in 1:length(J)
        V[i] = evaluate(J[i],[o[i] for i in 1:(p[1]-1)], [v[i] for i in 1:(p[1]-1)])
      end
      RRR, e = PolynomialRing(QQ,["e$i" for i in 1:2])
      w = helper_switch_ring(n,p,b,v,RRR)
      V2 = helper_finish_switch_ring(n,length(J),w,V)
    end
    II = ideal(RRR,V2)
    Radical = radical(II)
    Path = helper_component_check(Radical)
  end
  deg = compute_degree(Path)
  if isone(Path)
    println("Path variety is empty")
  else 
    dimension = dim(Path)
    println("Affine dimension is ", dimension) 
    if dimension == 0 
      println("Path multiplicity is ", deg)
    else
      error("Given point conditions were not generic enough or the polytope is not binodal.")
    end
  end
  return Path, deg
end

# b[1] = true means, that the point at position p[1] is a left out point, not part of a gap.
# For the moment, p and b are only allowed exactly 2 entries! It should be possible to expand the function if needed.
# p needs to be filled according to size, i.e., the smallest value first, the largest at the end.
# INPUT: The binodal variety, the positions of the two left out lattice points as a Vector, a generic point Vector of length n-2
# OUTPUT: prints the lattice path multiplicity if the points in v are generic enough.
function path_mult(I::MPolyIdeal{fmpq_mpoly}, p::Vector{Int}, b::Vector{Bool}, v::Vector)
  n = nvars(base_ring(I))
  @assert p[1]!= p[2] "Left out points must me different"
  @assert length(p)==2 "This function is only meant for lattice paths with 2 left out points."
  @assert length(p)==length(b) "$p and $b need to have the same length"
  @assert length(v)== n-length(p) "Number of point conditions in $v is wrong."
  @assert p[1] < p[2] "Left out points must be ordered from small to large."
  @assert b[1] || b[2] "This function is only meant for lattice paths with at most one gap."
  @assert (!(p[1] == 1 || p[1] == 2) || !b[1] || b[2] || p[2]>=4) "Invalid Path given." #This path is invalid because: 
  @assert (b[1] || !(p[1]<3||p[1]==n-1)) "Invalid Path given." #This path is invalid because: 
  J = gens(I)
  if b[1] && b[2] #case of the connected paths (no gap in the path)
    Path = con_mult(I,p,v)
  else #disconnected multiplicty
    if b[1] && !b[2] #first left out point, then gap
      o = helper_find_indices(p[2],p)
      V = Vector{fmpq_mpoly}(undef,length(J))
      for i in 1:length(J)
          V[i] = evaluate(J[i],[o[i] for i in 1:(p[2]-2)], [v[i] for i in 1:(p[2]-2)])
        end
        println(V)
        RRR, e = PolynomialRing(QQ,["e$i" for i in 1:2])
        w = helper_switch_ring(n,p,b,v,RRR)
        V2 = helper_finish_switch_ring(n,length(J),w,V)
    elseif !b[1] && b[2] # here the gap appears first
      o = helper_find_indices(p[1],p)
      V = Vector{fmpq_mpoly}(undef,length(J))
      for i in 1:length(J)
        V[i] = evaluate(J[i],[o[i] for i in 1:(p[1]-1)], [v[i] for i in 1:(p[1]-1)])
      end
      RRR, e = PolynomialRing(QQ,["e$i" for i in 1:2])
      w = helper_switch_ring(n,p,b,v,RRR)
      V2 = helper_finish_switch_ring(n,length(J),w,V)
    end
    II = ideal(RRR,V2)
    Radical = radical(II)
    Path = helper_component_check(Radical)
  end
  deg = compute_degree(Path)
  if isone(Path)
    println("Path variety is empty")
  else 
    dimension = dim(Path)
    println("Affine dimension is ", dimension) 
    if dimension == 0 
      println("Path multiplicity is ", deg)
    else
      error("Given point conditions were not generic enough or the polytope is not binodal.")
    end
  end
  return Path, deg
end


# First step to check for a generic point in the binodal variety (same function as con_mult)
# This function computes the intersection of the binodal variety for a 3-polytope given as a matrix filled with the lattice points, 
# with 4 linear spaces of codimension 1 given by substituting all variables apart from those with indices in p with the values given in v ("Path Ideal")
# INPUT: The n lattice points of the polytope as matrix, the positions of the two left out lattice points as a Vector, a generic point Vector of length n-2
# OUTPUT: the Path ideal
function find_generic_point(A::MatElem, p::Vector{Int}, v::Vector)
  Path = con_mult(A,p,v);
  return Path
end


#After the output of find_generic_point(A::MatElem, p::Vector{Int}, v::Vector)  was solved for e1 and e2 by hand
#we obtain a Vector filled with the values for a[i]
#This vector induces a surface in the binodal variety. 
# The following function computes the singular locus of a surface with coefficients in v and monomials given by the matrix A
#INPUT: the monomials given as matrix(ZZ,n,3,[...]) and the ocefficients given as Vector{Rational{Int64}} 
#OUTPUT: Pair (the ideal of the singular locus, its dimension)
function singular_locus(A::MatElem,v::Vector{Rational{Int64}})
  R, x =  PolynomialRing(QQ,["x$i" for i in 1:3])
  n = nrows(A)
  h = MPolyBuildCtx(R)
  for i in 1:n
      push_term!(h,v[i],Vector{Int64}([A[i,j] for j in 1:ncols(A)]))
  end
  p= finish(h)
  I = ideal(R,[p,derivative(p,x[1]),derivative(p,x[2]),derivative(p,x[3])])
  return I, dim(I)
end

#Works for a simple NumberField K
function singular_locus(A::MatElem,v::Vector{Any},K::AnticNumberField )
  R, x =  PolynomialRing(K,["x$i" for i in 1:3])
  n = nrows(A)
  h = MPolyBuildCtx(R)
  for i in 1:n
      push_term!(h,v[i]*one(K),Vector{Int64}([A[i,j] for j in 1:ncols(A)]))
  end
  p= finish(h)
  I = ideal(R,[p,derivative(p,x[1]),derivative(p,x[2]),derivative(p,x[3])])
  return I, dim(I)
end