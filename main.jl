include("./hopf_solution.jl")
using Polynomials
using LegendrePolynomials

#initialize auxiliary entities 
tt = 0.0:0.1:2.0
left_border_x = -3.0
right_border_x = 3.0
num_intervals = 3
intervals = get_intervals(left_border_x, right_border_x, num_intervals)
xx = get_midpoints(intervals)
frac_xx = unique(collect(Iterators.flatten(intervals)))
inner_frac_xx = frac_xx[2:(lastindex(frac_xx)-1)]
left_border_value = calc_left_border(tt)
right_border_value = calc_right_border(tt)
u_h_coeffs_tt = Vector{Vector{Vector{Float64}}}([])
u_h_tt = Vector{Vector{Float64}}([])

#initialize u_h(0)

u_h_0_coeffs = init_condition(intervals)
push!(u_h_coeffs_tt, u_h_0_coeffs)
u_h_0_values = calculate_u_h(u_h_0_coeffs, intervals)
push!(u_h_0_values, right_border_value[1])
pushfirst!(u_h_0_values, left_border_value[1])
push!(u_h_tt, u_h_0_values)

for t in 2:20
    coeffs = last(u_h_coeffs_tt)
    k1 = form_RK_k1(intervals, coeffs, left_border_value[t-1], right_border_value[t-1])
    k2 = form_RK_k2(intervals, coeffs, k1, left_border_value[t-1], right_border_value[t-1])
    k3 = form_RK_k3(intervals, coeffs, k2, left_border_value[t-1], right_border_value[t-1])
    k4 = form_RK_k4(intervals, coeffs, k3, left_border_value[t-1], right_border_value[t-1])
    new_coeffs = Vector{Vector{Float64}}([])
    for object in zip(coeffs, k1, k2, k3, k4)
        u₀ = object[1][1] + 0.1/6 * (object[2][1] + 2*object[3][1] + 2*object[4][1] + object[5][1])
        u₁ = object[1][2] + 0.1/6 * (object[2][2] + 2*object[3][2] + 2*object[4][2] + object[5][2])
        push!(new_coeffs, [u₀, u₁])
    end
    push!(u_h_coeffs_tt, new_coeffs)
    new_values = calculate_u_h(new_coeffs, intervals)
    push!(new_values , right_border_value[t])
    pushfirst!(new_values , left_border_value[t])
    push!(u_h_tt, new_values)
end