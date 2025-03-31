function get_intervals(left_border_x, right_border_x, num_intervals)::Vector{Tuple{Float64, Float64}}
    """
    Get num_intervals intervals within [left_border_x, left_border_y]
    """
    Δ = (right_border_x - left_border_x)/num_intervals
    intervals = []
    for i in 1:num_intervals
        push!(intervals, (left_border_x + (i-1) * Δ, left_border_x + i * Δ))
    end
    return intervals
end


function get_midpoints(intervals)
    """
    Return middle points of intervals
    """
    mid_points = []
    for interval in intervals
        push!(mid_points, (interval[1] + interval[2])/2)
    end
    return mid_points
end


function calc_right_border(tt)
    """
    Calculate right border function in tt points
    """
    all_roots = []
    for t in tt
        poly_to_find_roots = Polynomial([-1, 10, -6t, t*t])
        found_roots = roots(poly_to_find_roots)
        found_roots = filter(z -> imag(z) == 0, found_roots)
        push!(all_roots, real(found_roots[1]))
    end
    return all_roots
end


function calc_left_border(tt)
    """
    Calculate right border function in tt points
    """
    all_roots = []
    for t in tt
        poly_to_find_roots = Polynomial([-1, 10, 6t, t*t])
        found_roots = roots(poly_to_find_roots)
        found_roots = filter(z -> imag(z) == 0, found_roots)
        push!(all_roots, real(found_roots[1]))
    end
    return all_roots
end


function basis_function(x, interval, degree)
    """
    Return value of basis function ϕ⁽ˡ⁾(x) within interval where l = degree
    """
    return Pl(2(x-(interval[1] + interval[2])/2)/(interval[2] - interval[1]), degree)
end


function numerical_flux(a, b)
    """
    Return value of numerical flux corresponding values of u_h to the left and the right of the point 
    """
    if a <= b
        return min(a*a, b*b)/2
    end
    return max(a*a, b*b)/2
end


function init_condition(interval::Tuple{Float64, Float64})::Vector{Float64}
    """
    Return values of u₀(0) and u₁(0) within interval 
    """
    x = (interval[1] + interval[2])/2
    Δ = interval[2] - interval[1]
    u_0 = 1/Δ * (atan(interval[2]) - atan(interval[1]))
    u_1 = 6/(Δ*Δ) * ((log(1 + (interval[2])^2) - log(1 + (interval[1])^2))/2 - x * atan(interval[2]) + x * atan(interval[1]))
    return [u_0, u_1]
end


function init_condition(intervals::Vector{Tuple{Float64, Float64}})::Vector{Vector{Float64}}
    """
    Return values of u₀(0) and u₁(0)
    """
    coeffs = []
    for interval in intervals 
        push!(coeffs, init_condition(interval))
    end
    return coeffs
end

function calculate_u_h(coeff::Vector{Float64}, interval::Tuple{Float64, Float64})::Float64
    """
    Calculate uₕ(x, t) in middle point of interval using known u₀(t) and u₁(t)  
    """
    return coeff[1] + coeff[2] * basis_function((interval[1] + interval[2])/2, interval, 1)
end

function calculate_u_h(coeffs::Vector{Vector{Float64}}, intervals::Vector{Tuple{Float64, Float64}})
    """
    Calculate uₕ(x, t) in middle point of intervals using known u₀(t) and u₁(t)  
    """
    @assert length(coeffs) == length(intervals)
    values = []
    for (coeff, interval) in zip(coeffs, intervals)
        push!(values, calculate_u_h(coeff, interval))
    end
    return values
end


function calculate_u_h(coeff, interval, x)
    """
    Calculate uₕ(x, t) in point x within interval using known u₀(t) and u₁(t)  
    """
    @assert( (interval[1] <= x) & (x <= interval[2]) )
    return coeff[1] + coeff[2] * basis_function(x, interval, 1)
end


function form_intervals_coeffs(intervals, coeffs)
    """
        Gets intervals to special form for count deriative function
    """
    
    superintervals = []
    supercoeffs = []
    last_index = lastindex(intervals)
    push!(superintervals, (intervals[1], intervals[2]))
    push!(supercoeffs, (coeffs[1], coeffs[2]))
    for i in 2:(last_index-1)
        push!(superintervals, (intervals[i-1], intervals[i], intervals[i+1]))
        push!(supercoeffs, (coeffs[i-1], coeffs[i], coeffs[i+1]))
    end
    push!(superintervals, (intervals[last_index - 1], intervals[last_index]))
    push!(supercoeffs, (coeffs[last_index - 1], coeffs[last_index]))

    return superintervals, supercoeffs
end


function form_coeffs(coeffs)
    """
        Gets intervals to special form for count deriative function
    """
    
    supercoeffs = []
    last_index = lastindex(coeffs)
    push!(supercoeffs, (coeffs[1], coeffs[2]))
    for i in 2:(last_index-1)
        push!(supercoeffs, (coeffs[i-1], coeffs[i], coeffs[i+1]))
    end
    push!(supercoeffs, (coeffs[last_index - 1], coeffs[last_index]))

    return supercoeffs
end


function calculate_u_0_deriative(interval, first_flux, second_flux)
    Δ = interval[2] - interval[1]
    return -1/Δ * (second_flux - first_flux)
end


function calculate_u_1_deriative(coeffs, interval, first_flux, second_flux)
    Δ = interval[2] - interval[1]
    xⱼ = (interval[2] + interval[1])/2
    u₀ = coeffs[1]
    u₁ = coeffs[2]
    intergral_part = 3/Δ *(u₀^2 + u₀ * u₁ * ((interval[2] - xⱼ)^2  - (interval[1] - xⱼ)^2) + (u₁^2) /3 * ((interval[2] - xⱼ)^3  - (interval[1] - xⱼ)^3))
    return intergral_part -1/Δ * (second_flux + first_flux)
end


function deriative_right_part(intervals, coeffs)::Vector{Float64}
    """
        Calculate the right of deriative within interval

        intervals must have form

        [leftmost_interval, interval, rightmost_interval]
        
        coeffs must be Vector{Turple{Float64, Float64}} that have following form:

        [coeffs_in_left_interval, coeffs_within_interval, coeffs_in_right_interval]

        - coeffs_in_left_interval - u₀ and u₁ coeffs in the interval to the left of the given one
        - coeffs_within_interval - u₀ and u₁ coeffs in the given interval 
        - coeffs_in_right_interval - u₀ and u₁ coeffs in the interval to the right of the given one


        if left = true then we work with the leftmost interval, coeffs must have form [coeffs_within_interval, coeffs_in_right_interval]
        if right = true then we work with the rightmost interval, coeffs must have form [coeffs_in_left_interval, coeffs_within_interval]

        if left | right == true then we will need an auxiliar_value.
        left & right must be equal false.
    """

    @assert (length(coeffs) == 3)
    @assert (length(intervals) == 3)

    leftmost_interval = intervals[1]
    interval = intervals[2]
    rightmost_interval = intervals[3]

    first_x_point = interval[1]
    second_x_point = interval[2]

    leftmost_coeffs = coeffs[1]
    within_coeffs = coeffs[2]
    rightmost_coeffs = coeffs[3]

    first_left_value = calculate_u_h(leftmost_coeffs, leftmost_interval, first_x_point)
    first_right_value = calculate_u_h(within_coeffs, interval, first_x_point)
    second_left_value = calculate_u_h(within_coeffs, interval, second_x_point)
    second_right_value = calculate_u_h(rightmost_coeffs, rightmost_interval, second_x_point)

    first_flux = numerical_flux(first_left_value, first_right_value)
    second_flux = numerical_flux(second_left_value, second_right_value)

    u_h_der = [calculate_u_0_deriative(interval, first_flux, second_flux),
               calculate_u_1_deriative(within_coeffs, interval, first_flux, second_flux)]

    return u_h_der
end


function deriative_right_part(intervals::Tuple{Tuple{Float64, Float64}, Tuple{Float64, Float64}}, 
                              coeffs,
                              auxiliary_value::Float64; left::Bool = false)::Vector{Float64}
    """
        Calculate the right of deriative within interval

        intervals must have form

        [leftmost_interval, interval, rightmost_interval]
        
        coeffs must be Vector{Turple{Float64, Float64}} that have following form:

        [coeffs_in_left_interval, coeffs_within_interval, coeffs_in_right_interval]

        - coeffs_in_left_interval - u₀ and u₁ coeffs in the interval to the left of the given one
        - coeffs_within_interval - u₀ and u₁ coeffs in the given interval 
        - coeffs_in_right_interval - u₀ and u₁ coeffs in the interval to the right of the given one


        if left = true then we work with the leftmost interval, coeffs must have form [coeffs_within_interval, coeffs_in_right_interval]
        if right = true then we work with the rightmost interval, coeffs must have form [coeffs_in_left_interval, coeffs_within_interval]

        if left | right == true then we will need an auxiliar_value.
        left & right must be equal false.
    """

    @assert (length(coeffs) == 2)
    @assert (length(intervals) == 2)

    if left
        interval = intervals[1]
        rightmost_interval = intervals[2]

        first_x_point = interval[1]
        second_x_point = interval[2]

        within_coeffs = coeffs[1]
        rightmost_coeffs = coeffs[2]

        first_left_value = auxiliary_value
        first_right_value = calculate_u_h(within_coeffs, interval, first_x_point)
        second_left_value = calculate_u_h(within_coeffs, interval, second_x_point)
        second_right_value = calculate_u_h(rightmost_coeffs, rightmost_interval, second_x_point)

    else
        leftmost_interval = intervals[1]
        interval = intervals[2]

        first_x_point = interval[1]
        second_x_point = interval[2]
        
        leftmost_coeffs = coeffs[1]
        within_coeffs = coeffs[2]

        first_left_value = calculate_u_h(leftmost_coeffs, leftmost_interval, first_x_point)
        first_right_value = calculate_u_h(within_coeffs, interval, first_x_point)
        second_left_value = calculate_u_h(within_coeffs, interval, second_x_point)
        second_right_value = auxiliary_value
    end

    first_flux = numerical_flux(first_left_value, first_right_value)
    second_flux = numerical_flux(second_left_value, second_right_value)

    u_h_der = [calculate_u_0_deriative(interval, first_flux, second_flux),
               calculate_u_1_deriative(within_coeffs, interval, first_flux, second_flux)]

    return u_h_der
end


function form_RK_k1(intervals, coeffs, left_border_value, right_border_value)::Vector{Vector{Float64}}
    superintervals, supercoeffs = form_intervals_coeffs(intervals, coeffs)
    leftmost_intervals = superintervals[1]
    rightmost_intervals = last(superintervals)
    leftmost_coeffs = supercoeffs[1]
    rightmost_coeffs = last(supercoeffs)
    regular_intervals = superintervals[2:(lastindex(superintervals) - 1)]
    regular_coeffs = supercoeffs[2:(lastindex(supercoeffs) - 1)]
    leftmost_der_value = deriative_right_part(leftmost_intervals, leftmost_coeffs, left_border_value, left = true)
    rightmost_der_value = deriative_right_part(rightmost_intervals, rightmost_coeffs, right_border_value, left = false)
    der_values = []
    for (reg_intervals, reg_coeffs) in zip(regular_intervals, regular_coeffs)
        push!(der_values, deriative_right_part(reg_intervals, reg_coeffs))
    end
    push!(der_values, rightmost_der_value)
    pushfirst!(der_values, leftmost_der_value)
    return der_values
end


function form_RK_k2(intervals, coeffs, k1, left_border_value, right_border_value)

    k1_coeffs::Vector{Vector{Float64}} = []
    for (coeff, k_1_values) in zip(coeffs, k1)
        u₀ = coeff[1]
        u₁ = coeff[2]
        push!(k1_coeffs, [u₀ + 0.1/2 * k_1_values[1], u₀ + 0.1/2 * k_1_values[2]])
    end

    superintervals, supercoeffs = form_intervals_coeffs(intervals, k1_coeffs)
    leftmost_intervals = superintervals[1]
    rightmost_intervals = last(superintervals)
    leftmost_coeffs = supercoeffs[1]
    rightmost_coeffs = last(supercoeffs)
    regular_intervals = superintervals[2:(lastindex(superintervals) - 1)]
    regular_coeffs = supercoeffs[2:(lastindex(supercoeffs) - 1)]
    leftmost_der_value = deriative_right_part(leftmost_intervals, leftmost_coeffs, left_border_value, left = true)
    rightmost_der_value = deriative_right_part(rightmost_intervals, rightmost_coeffs, right_border_value, left = false)
    der_values = []
    for (reg_intervals, reg_coeffs) in zip(regular_intervals, regular_coeffs)
        push!(der_values, deriative_right_part(reg_intervals, reg_coeffs))
    end
    push!(der_values, rightmost_der_value)
    pushfirst!(der_values, leftmost_der_value)
    return der_values
end


function form_RK_k3(intervals, coeffs, k2, left_border_value, right_border_value)

    k2_coeffs::Vector{Vector{Float64}} = []
    for (coeff, k_2_values) in zip(coeffs, k2)
        u₀ = coeff[1]
        u₁ = coeff[2]
        push!(k2_coeffs, [u₀ + 0.1/2 * k_2_values[1], u₀ + 0.1/2 * k_2_values[2]])
    end

    superintervals, supercoeffs = form_intervals_coeffs(intervals, k2_coeffs)
    leftmost_intervals = superintervals[1]
    rightmost_intervals = last(superintervals)
    leftmost_coeffs = supercoeffs[1]
    rightmost_coeffs = last(supercoeffs)
    regular_intervals = superintervals[2:(lastindex(superintervals) - 1)]
    regular_coeffs = supercoeffs[2:(lastindex(supercoeffs) - 1)]
    leftmost_der_value = deriative_right_part(leftmost_intervals, leftmost_coeffs, left_border_value, left = true)
    rightmost_der_value = deriative_right_part(rightmost_intervals, rightmost_coeffs, right_border_value, left = false)
    der_values = []
    for (reg_intervals, reg_coeffs) in zip(regular_intervals, regular_coeffs)
        push!(der_values, deriative_right_part(reg_intervals, reg_coeffs))
    end
    push!(der_values, rightmost_der_value)
    pushfirst!(der_values, leftmost_der_value)
    return der_values
end


function form_RK_k4(intervals, coeffs, k3, left_border_value, right_border_value)

    k3_coeffs::Vector{Vector{Float64}} = []
    for (coeff, k_3_values) in zip(coeffs, k3)
        u₀ = coeff[1]
        u₁ = coeff[2]
        push!(k3_coeffs, [u₀ + 0.1 * k_3_values[1], u₀ + 0.1 * k_3_values[2]])
    end

    superintervals, supercoeffs = form_intervals_coeffs(intervals, k3_coeffs)
    leftmost_intervals = superintervals[1]
    rightmost_intervals = last(superintervals)
    leftmost_coeffs = supercoeffs[1]
    rightmost_coeffs = last(supercoeffs)
    regular_intervals = superintervals[2:(lastindex(superintervals) - 1)]
    regular_coeffs = supercoeffs[2:(lastindex(supercoeffs) - 1)]
    leftmost_der_value = deriative_right_part(leftmost_intervals, leftmost_coeffs, left_border_value, left = true)
    rightmost_der_value = deriative_right_part(rightmost_intervals, rightmost_coeffs, right_border_value, left = false)
    der_values = []
    for (reg_intervals, reg_coeffs) in zip(regular_intervals, regular_coeffs)
        push!(der_values, deriative_right_part(reg_intervals, reg_coeffs))
    end
    push!(der_values, rightmost_der_value)
    pushfirst!(der_values, leftmost_der_value)
    return der_values
end