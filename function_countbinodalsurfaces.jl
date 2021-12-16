#Given a 3-polytope by means of a matrix filled with the lattice points,
#this function computes the binodal variety of those polynomials with the given support 
#and prints the affine dimension, degree and whether the binodal variety contains monomials

#INPUT: matrix(ZZ,n,3,[...])
#OUTPUT: returns binodal variety
function binodal(A::MatElem)
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
    #  @show Vector{Int64}([M[i,j] for j in 1:ncols(M)])
      push_term!(h,1,Vector{Int64}([M[i,j] for j in 1:ncols(M)]))
    end
    p= finish(h)
    q = evaluate(p,[x,y,z],[u,v,w])
    I = ideal(R,[p,q,derivative(p,x),derivative(p,y),derivative(p,z),derivative(q,u),derivative(q,v),derivative(q,w),1-t*x*y*z*u*v*w])
    Sat = saturation(I, ideal(R,[x-u,y-v,z-w]))
    K = eliminate(Sat,[t,x,y,z,u,v,w])
    J = gens(K)
    RR, a = PolynomialRing(QQ, ["a$i" for i in 1:n])
    V = Vector(undef,length(J))
    for i in 1:length(J)
        V[i] = evaluate(J[i],vcat(a,[zero(RR) for i in 1:7]))
    end
    ##V is currently of type Any it would be better to have it as a Vector of type fmpq_mpoly
    II = ideal(RR,V)
    Radical = radical(II)
    return Radical
end

#############################
#Given a 3-polytope by means of a matrix filled with the lattice points,
#prints the affine dimension, degree of the associated binodal variety and whether it contains monomials

#INPUT: matrix(ZZ,n,3,[...])
#OUTPUT: prints dimension, degree and whether the binodal variety contains monomials
function investigate_binodalvariety(A::MatElem) 
    n = nrows(A)
    Var = binodal(A);
    dimension = dim(Var)
    println("Affine Dimension is ", dimension)
    if dimension == n-2
        println("The binodal variety is of expected affine dimension.")
    else
        println("The binodal variety has not the expected affine dimension.")
    end
    # the degree is by definition the degree of RR/Radical, but everything must be homogenous
    RR, a = PolynomialRing(QQ, ["a$i" for i in 1:n])
    RRhom, = grade(RR) # This is RR, homogenous with weights [1, 1,..., 1]
    IIhom = ideal(RRhom, RRhom.(gens(Var)))
    Q, _ = quo(RRhom, IIhom)
    deg = degree(quo(RRhom, IIhom)[1])
    println("Degree is ", deg)
    # Soon there will be a all_monomials(polynomial_ring, degree) function,
    # but for the moment we can do the following
    # Note that we are checking more then we should, but it does not matter
    # for these small parameter values
    monomial_iterator = AbstractAlgebra.ProductIterator(a, n)
    # this is the iterator of d x d x ... x d (n times)
    found = false
    for mon in monomial_iterator
      # m = is of the form (d[2], d[1], ..., d[2])
      m = prod(mon)
      if m in Var
        found = true
        println("Contains monomials.")
        break
      end
    end
    if !found
      println("Does not contain monomials.")
      if dimension == n-2
        println(A, " could be a binodal polytope.")
      end
    end
end

###############################
# To compute the path multiplicities for a binodal polytope of a connected path
# INPUT: The n lattice points of the polytope as matrix, the positions of the two left out lattice points, a generic point Vector of length n-2
# OUTPUT: prints the lattice path multiplicity if the points are generic enough.
function con_mult(A::MatElem, s::Int, t::Int, v::Vector)
    if s == t
      println("Error: Left out points must me different")
    else
      I = binodal(A)
      J = gens(I);
      n = nrows(A)
      o = Vector{Int64}(undef,n-2) #this is auxiliary vector to help us fill the point conditions in the vector V
      for i in 1:n
        if i<s && i<t
          o[i] = i
        elseif s<i && i<t
          o[i-1] = i
        elseif i<s && t<i
          o[i-1] = i
        elseif s<i && t<i
          o[i-2] = i
        end
      end
      V = Vector(undef,length(J))# we fill this vector according to the left out point conditions s and t with the values from v
      for i in 1:length(J)
          V[i] = evaluate(J[i],[o[i] for i in 1:(n-2)], [v[i] for i in 1:(n-2)])
      end
      RRR, e = PolynomialRing(QQ,["e$i" for i in 1:2])
      w = Vector(undef,n)
      for i in 1:n
        if i!=s && i!=t
          w[i] = zero(RRR)
        elseif (s==i && s<t) || (t==i && t<s)
          w[i] = e[1]
        elseif (s==i && t<s) || (t==i && s<t)
          w[i] = e[2]
        end
      end
      V2 = Vector(undef,length(J))
      for i in 1:length(J)
        V2[i] = evaluate(V[i],[w[i] for i in 1:n])
      end
      II = ideal(RRR,V2) # this is the ideal of the intersection of the binodal variety with the hyperplanes given by v
      Radical = radical(II)
      dimension = dim(Radical)
      println("Affine Dimension is ", dimension) 
      if dimension == 0 #this means that the point conditions were general enough
        Gens = gens(Radical)
        Genshom = Vector(undef, length(Gens))
        for i in 1:length(Gens) 
          Genshom[i] = homogenization(Gens[i],"e3",3)
        end
        RRhom = parent(Genshom[1]) # This is RR, homogenous with additional variable e3 and weights [1, 1,..., 1]
        IIhom = ideal(RRhom, RRhom.(Genshom))
        Q, _ = quo(RRhom, IIhom)
        deg = degree(quo(RRhom, IIhom)[1]) #this is the degree of the hyperplane section of the binodal variety in the left out points
        println("Path Multiplicity is ", deg)
        return deg;
      else
        println("Error: Point conditions not generic enough or not binodal polytope. The intersection is empty or too large.")
      end
    end
  end