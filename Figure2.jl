#FIGURE 3, c)

struct parameters
    X::Float64
    Y1::Float64
    Y2::Float64
    bx::Float64
    ux::Float64
    ui::Float64
    ei::Float64
    bi::Float64
    c::Float64
    βi::Float64
    t::Float64
end

#each strain has its own values of:
r1 = rand(Float64) #random number between 0-1
r2 = rand(Float64) #random number between 0-1
r3 = rand(Float64) #chose from a uniform distribution over (0, r1) for my figure4

bi = bx*r3*(1-r1*r2)
ei = bx*(1-r3)*(1-r1*r2)
ui = rand(Float64)
βi = r1 - α*bi/bx
    return β = 0
        if β =< 0
        else β = β
    end

c = #0.5 and 4 -- constant
bx = 1 #constant
ux = 0.2 #constant

X = 1
Y1 = 1
Y2 = 1

t = 0

#Equations 1,2,3
dXdt = (bx*X + e1*Y1 +e2*Y2)*(1 - (X + Y1 + Y2)/K) - ux*X - c*(B1*Y1+B2*Y2)*X
dY1dt = b1*Y1*(1 - (X + Y1 + Y2)/K) - u1*Y1 +c*B1*X*Y1
dY2dt = b2*Y2*(1 - (X + Y1 + Y2)/K) - u2*Y2 +c*B2*X*Y2

# Values of X,Y1,Y2 at next time step
X = X + dXdt
Y1 = Y1 + dY1dt
Y2 = Y2 + dY2dt

#need to introduce 1000 strains, one at a time, every 1000 host generations
