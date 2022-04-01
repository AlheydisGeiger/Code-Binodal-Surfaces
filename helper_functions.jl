#This document contains helper functions for the investigation of binodal varieties in code_binodal.jl

#Function to find fill a vector with the indices of the variables that have to be substituted by fixed values
#INPUT: integer n giving the maximal index up to which we consider the substitution, vector p stating the indices that are not to be substituted
#OUTPUT: Vector{Int64}
function helper_find_indices(n::Int,p::Vector{Int})
    o = Vector{Int64}() #this is auxiliary vector to help us fill the point conditions in the vector V
    s = p[1]
    t = p[2]
    for i in 1:n
        if i<s && i<t
            append!(o,i)
        elseif s<i && i<t
            append!(o,i)
        elseif i<s && t<i
            append!(o,i)
        elseif s<i && t<i
            append!(o,i)
        end
    end
    return o
end

# First helper function to switch from the Polynomial ring with n variables to 2 variables in the Computation of the Path multiplcities
# Depending on connected (b[1]&&b[2]) or disconnected paths (b[1] &&!b[2])||(!b[1]&&b[2]) the vector w gets filled differently
# This is explained in detail in Algorithms 5 & 6 in the PhD Thesis "Tropical Geometric Counting Problems" by the author
#INPUT: Integer n, Vector p, Vector of bools b, PolynomialRing RRR
#OUTPUT: Vector w filled with zero(RRR) for all indices not in p and with monomials in e[1], e[2] for the entries of p
function helper_switch_ring(n::Int,p::Vector{Int},b::Vector{Bool},v::Vector,RRR::FmpqMPolyRing)
    w = Vector(undef,n)
    RRR, e = PolynomialRing(base_ring(RRR),symbols(RRR))
    if b[1] && b[2]
        for i in 1:n
            if i!=p[1] && i!=p[2]
                w[i] = zero(RRR)
            elseif (p[1]==i && p[1]<p[2]) || (p[2]==i && p[2]<p[1])
                w[i] = e[1]
            elseif (p[1]==i && p[2]<p[1]) || (p[2]==i && p[1]<p[2])
            w[i] = e[2]
            end
        end
    elseif b[1] && !b[2]
        for i in 1:n
            if i<p[1] || p[1]<i<p[2]
              w[i] = zero(RRR)
            elseif i == p[1]
              w[i] = e[1]
            elseif i == p[2]
              w[i] = e[2]
            else
              w[i] = v[i-2]*e[2]
            end
        end
    elseif !b[1] && b[2]
        for i in 1:n
            if i<p[1]
              w[i] = zero(RRR)
            elseif i == p[1]
              w[i] = e[1]
            elseif i == p[2]
              w[i] = e[1]*e[2]
            elseif p[1]<i<p[2]
              w[i] = v[i-1]*e[1]
            else 
              w[i] = v[i-2]*e[1]
            end
          end
    end
    return w
end

# Second helper function to switch from the Polynomial ring with n variables to 2 variables in the Computation of the Path multiplcities
#INPUT:  Integer n for the number of inital variables, Integer m for the number of polynomials that have to be substituted, 
#        vector w that contains the polynomials in the new variables e[1], e[2] that have to be substituted in the indices complementary to those found by using helper_find_indices
#        vector V in which the entries with indices found by helper_find_indices are already substituted by fixed values
#OUTPUT: Vector V2 that consists for each variable of the original ring of fixed values or monomials in e[1],e[2]
function helper_finish_switch_ring(n::Int,m::Int,w::Vector,V::Vector)
    V2 = Vector(undef,m)
    for i in 1:m
        V2[i] = evaluate(V[i],[w[i] for i in 1:n])
    end
    return V2
end

#This function takes an ideal and computes its primary decomposition
#For each primary component it checks whether the associated prime ideal contains monomials
#The output is the radical ideal consisting of the intersection of all those prime ideals that do not contain monomials
#INPUT: Ideal I
#OUTPUT: the maximal radical ideal in I such that no primary component contains monomials
function helper_component_check(I::MPolyIdeal{fmpq_mpoly})
    RRR = base_ring(I)
    k = nvars(RRR)
    monomes = RRR(Rational{BigInt}[1],[fill(1,k)]) #e[1]*e[2]
    Mon = ideal(RRR, monomes)
    L = primary_decomposition(I)
    II = ideal(RRR, one(RRR)) 
    l = length(L)
    for i in 1:l
        if !issubset(Mon,L[i][2])
            II = intersect(II,L[i][2])
        end
    end
    return II
end


