using JLD2
using Primes
using IntervalArithmetic
using Optim

setdisplay(; sigdigits=20)

# Defining constants
V_list_small_length = 10^8
π_up = interval(3.14159265358979323847)
π_low = interval(3.14159265358979323846)
π_2_up = interval(0.6601619)
π_2_low = interval(0.6601618)
li_2_up = interval(1.0451638)
li_2_low = interval(1.0451637)
V_step = 10^4
L_1 = 10^8
L_2 = 10^10
c_3 = interval(12.6244)
L_4 = 10^7

if isfile("V_list.jld2")
    @load "V_list.jld2" V_list
end
if isfile("V_list_small.jld2")
    @load "V_list_small.jld2" V_list_small
end

# Calculates f(z*) exactly
# Lemma 2.2
function f(z)
    sum = interval(0)
    product = prod(primes(3,z))
    divs = divisors(product)[2:end]
    for d in divs
        if length(factor(d)) % 2 == 1
            sum += interval(1) / (interval(2)^interval(length(factor(d))) * interval(d)^interval(2))
        else
            sum -= interval(1) / (interval(2)^interval(length(factor(d))) * interval(d)^interval(2))
        end
    end
    return sum
end

# Stores each f(z*) to avoid recalculating
# The maximum value of z* was set to 37 to avoid excessive computation
f_zs = []
for z in primes(3,37)
    push!(f_zs, f(z))
end

# Calculates upper bound on the portion of the integal containing the V(z) term: Int_startp^endp{Li(t)/(V*t^2)} dt
# Lemma 3.1
function int_V(startp, endp, V_start)
    I_end = log(log(endp)) - (Li_low(endp) - li_2_up)/endp
    I_start = log(log(startp)) - (Li_up(startp)-li_2_low)/startp
    I = I_end - I_start
    return I / V_start
end

# Calculates upper bound on the portion of the integal containing the R(x,z) term: Int_startp^endp{r(z)^2 * c_π * sqrt(t)log(t) / t^2} dt
# Lemma 3.1
function int_R(startp, endp, z_end)
    I_end = -interval(2) * (log(endp) + interval(2)) / sqrt(endp)
    I_start = -interval(2) * (log(startp) + interval(2)) / sqrt(startp)
    I = I_end - I_start
    return c_π_up(startp) * I * r_coeff_up(z_end)^interval(2)
end

# Calculates li_up(x)
# Lemma 3.1
function Li_up(x)
    return π1_up(x) + sqrt(x)*log(x)/(interval(8) * π_low) + li_2_up
end

# Calculates li_low(x)
# Lemma 3.1
function Li_low(x)
    return π1_low(x) - sqrt(x)*log(x)/(interval(8) * π_low) - li_2_up
end

# Calculates upper bound on π(x)
# Lemma 3.1
function π1_up(x)
    if inf(x) >= 4 * 10^9
        return x / log(x) * (interval(1) + interval(1) / log(x) + interval(2) / log(x)^interval(2) + interval(7.32) / log(x)^interval(3))
    else
        throw(DomainError(x, "x must be greater than or equal to 4 * 10^9"))
    end
end

# Calculates lower bound on π(x)
# Lemma 3.1
function π1_low(x)
    if inf(x) >= 88789
        return x / log(x) * (interval(1) + interval(1) / log(x) + interval(2) / log(x)^interval(2))
    else
        throw(DomainError(x, "x must be greater than or equal to 88789"))
    end
end

# Calculates the function r(z)
# Lemma 2.2
function r_coeff_up(z)
    if sup(z) < L_4
        return interval(M_list[floor(Int, z)])
    else
        return (interval(4)/π_up^interval(2) + interval(2)/(sqrt(z)-interval(1)))*z + interval(0.5) * sqrt(z)/interval(2)
    end
end

# Calculates upper bound on c_π
# Theorem 2.2
function c_π_up(x)
    return interval(3)/(interval(8) * π_low) + (interval(6) + interval(1)/interval(π_low)) / (interval(4)log(x)) + interval(6)/log(x)^interval(2) + li_2_up/(sqrt(x) * log(x))
end

# Calculates upper bound on D(L,z)
# Lemma 2.1
function D(L, z)
    return interval(3)*c_3 / sqrt(z) - interval(9)*c_3 / sqrt(L)
end

# Store the values of M(z) for z up to L_4 (set at 10^7)
function create_M_list()
    M_list = zeros(Int32, L_4)
    count = 0
    for i in 1:L_4
        if i % 2 == 1 && squarefree(i)
            count += 1
        end
        M_list[i] = count
    end
    @save "M_list.jld2" M_list
end

# Store the values of V(z) for z up to V_list_small_length (set at 10^8)
function create_V_list_small()
    V_list_small = zeros(Float64, V_list_small_length)
    sum = interval(1)
    for i in 3:V_list_small_length
        if squarefree(i) && i % 2 == 1
            prod = interval(1)
            for p in factor(i)
                prod *= interval(1) / (interval(p[1]) - interval(2))
            end
            sum += prod
        end
        V_list_small[i] = inf(sum)
    end
    @save "V_list_small.jld2" V_list_small
end

# Store the values of V(z') for z' up to x; allows for using an existing calculation as a starting point
function create_V_list(x)
    if isfile("V_list.jld2")
        @load "V_list.jld2" V_list
    else
        V_list = Float64[]
    end
    start = floor(Int, length(V_list) / V_step)

    if length(V_list) >= floor(Int, x / V_step)
        return
    end

    # Can use previous calculations as starting point
    if length(V_list) > 0
        sum = interval(V_list[end])
    else
        sum = interval(0)
    end

    V_list = vcat(V_list, zeros(Float64, floor(Int, x / V_step) - start))
    prod = interval(1)
    for i in 1 + start:floor(Int, x/V_step)
        # Save incrementally
        if (i * V_step) % (10^8) == 0
            @save "V_list.jld2" V_list
        end
        # Compute the remaining values of V(z')
        for n in (V_step * (i-1) + 1):1:(V_step*i)
            if squarefree(n) && n % 2 == 1
                if n == 1
                    sum += interval(1)
                    continue
                end
                prod = interval(1)
                for p in factor(n)
                    prod *= interval(1) / (interval(p[1]) - interval(2))
                end
                sum += prod
            end
        end
        V_list[i] = inf(sum)
    end
    @save "V_list.jld2" V_list
end

# Returns a lower bound on V(z)
# Lemma 2.1
function V_est_low(z)
    if inf(z) > L_2
        sum = V_est_low(L_2) + (log(z) - log(interval(L_2))) / (interval(2) * π_2_up) + D(interval(L_2), z)
    elseif inf(z) > L_1
        sum = interval(V_list[floor(Int, inf(z/interval(V_step)))])
    elseif inf(z) > V_step
        sum = V_low(z)
    else
        sum = interval(V_list_small[floor(Int,inf(z))])
    end
    return sum
end    

# Calculates V(z) explicitly using V(z') as a starting point
function V_low(z)
    z = interval(floor(Int, inf(z)))
    # Begin at the sum for the nearest z' below z
    if floor(Int, inf(z / interval(V_step))) > 0
        sum = interval(V_list[floor(Int, inf(z / interval(V_step)))])
    else
        sum = interval(0)
    end
    # Add the missing values to the sum
    if sup((floor(Int, inf(z / interval(V_step))) * interval(V_step))) + 1 <= inf(z)
        for i in ceil(Int, sup((interval(floor(Int, inf(z / interval(V_step))))*interval(V_step))))+1:1:floor(Int, inf(z))
            if squarefree(i) && i % 2 == 1
                prod = interval(1)
                for p in factor(i)
                    prod *= interval(1) / (interval(p[1]) - interval(2))
                end
                sum += prod
            end
        end
    end
    return sum
end

# Determines if a number is squarefree
function squarefree(x)
    for i in values(factor(x))
        if i > 1
            return false
        end
    end
    return true
end

# Calculates the bound for B(m_k, ∞)
# Lemma 3.2
function bound_remainder(x)
    return interval(16) * interval(2) * interval(π_2_up) / log(x) + interval(8) / sqrt(x)
end

# Calculates c_1
# Lemma 3.1
function calculate_c_1(startp, endp)
    I_end = log(log(endp)) - (Li_low(endp) - li_2_up)/endp
    I_start = log(log(startp)) - (Li_up(startp)-li_2_low)/startp
    return I_end - I_start
end

# Calculates c_2
# Lemma 3.1
function calculate_c_2(startp, endp)
    I_end = -interval(2) * (log(endp) + interval(2)) / sqrt(endp)
    I_start = -interval(2) * (log(startp) + interval(2)) / sqrt(startp)
    I = I_end - I_start
    return c_π_up(startp) * I
end

# Optimises B(Start, End) given an optimal z from a previous calculation of B(i, Start)
function parameter_sweep(Start, End, start_z)
    c_1 = calculate_c_1(Start, End)
    c_2 = calculate_c_2(Start, End)
    function int(z)
        return sup(c_1 / V_est_low(z) + c_2 * r_coeff_up(z)^interval(2))
    end
    # Optimises for z
    opt = optimize(int, start_z, start_z*(End/Start))
    z = Optim.minimizer(opt)
    return (sup(interval(2) * interval(int(z))), z)
end

# Creates a file storing optimal parameters found when calculating B(m_1, m_k)
function create_B_parameters(midpoints, Start_z)
    b_result = Dict()
    b_result["midpoints"] = midpoints
    zs = []
    sums = []
    for (i, m) in enumerate(midpoints)
        print(".")
        if i < length(midpoints)
            if i == 1
                params = parameter_sweep(midpoints[i], midpoints[i+1], Start_z)
            else
                params = parameter_sweep(midpoints[i], midpoints[i+1], zs[end][end])
            end
            push!(zs, params[2])
            push!(sums, params[1])
        end
    end
    b_result["zs"] = zs
    b_result["sums"] = sums
    b_result["total_sum"] = sup(sum([interval(s) for s in b_result["sums"]]))
    @save "B_result.jld2" b_result
    return sum(sums)
end

# Uses the previously calculated parameters to give an upper bound on B(a, b), where m_1<a<b<m_k
function calculate_B_from_result(a, b)
    @load "B_result.jld2" b_result
    if a >= b
        throw(DomainError((a, b), "a must be less than b"))
    elseif a < b_result["midpoints"][1] || b > b_result["midpoints"][end]
        throw(DomainError((a, b), "a must be >= m_1 and b must be <= m_k"))
    end
    if a in b_result["midpoints"]
        startp = findfirst(i -> i == a, b_result["midpoints"])
    else
        startp = findfirst(i -> i > a, b_result["midpoints"])-1
    end
    endp = findfirst(i -> i >= b, b_result["midpoints"])-1
    return sup(sum([interval(s) for s in b_result["sums"][startp:endp]]))
end


# The parameters used to calculate B as in the paper
m = 0.2
midpoints1 = [19 + i*m for i in 0:1981/m]
midpoints = union([4*BigInt(10)^18], [BigInt(10)^i for i in midpoints1])
@time create_B_parameters(midpoints, 650)

println("B: $(sup(sum(b_result["sums"]) + bound_remainder(BigInt(10)^2000) + interval(1.840518)))")