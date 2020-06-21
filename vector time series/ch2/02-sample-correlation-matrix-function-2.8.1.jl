#!/usr/bin/env -S julia --sysimage=../sys_vts.so --

using HypothesisTests
using DataFrames
using CSV
using ArgParse
using Statistics
using Distributions

aps = ArgParseSettings(
    description = "This program computes the sample correlation matrix for a given " *
                  "difference D^(d) of a multi-variate time series at a set of given lags. " *
                  "We suppose that the time series is" *
                  " presented in order from earlier to later," *
                  " and that the intervals between" *
                  " each time point are equal.",
)
@add_arg_table! aps begin
    "--time-column"
    help = "name of the column containing the time information, so we can ignore this column"
    arg_type = String
    default = "Period"
    
    "--p-value"
    help = "p-value for significance test whether the cross-correlations are significantly different from zero (based on t-test)"
    arg_type = Real
    default = 0.045
    
    "--diff"
    help = "number of times the difference operator D is applied to the input series"
    arg_type = Int
    default = 0
    
    "--lag-from"
    help = "minimal lag k"
    arg_type = Int
    default = 0
    
    "--lag-to"
    help = "maximal lag k"
    arg_type = Int
    default = nothing
    
    "--store"
    help = "store the sample correlation matrix function table"
    arg_type = String
    required = false
    
    "csv"
    help = "path to the CSV containing the vector time series"
    required = false
    arg_type = String
    default = "WW2a.csv"
end

args = parse_args(aps)

###

function diff_series(x, how_often::Int = 1)
    if how_often == 0
        x
    elseif how_often > 0
        x0 = diff_series(x, how_often - 1)
        (x0[2:length(x0)] - x0[1:length(x0)-1])
    else
        throw("Integration is not implemented in diff_series!")
    end
end

###

df = CSV.read(args["csv"])
time_col = Symbol(args["time-column"])
apply_D_count = args["diff"]

p_value = args["p-value"]

lag_from = args["lag-from"]
lag_to = args["lag-to"]

if lag_to === nothing
    lag_to = size(df)[1] - 3 - apply_D_count
end



components = [col for col in names(df) if !(col == time_col)]

diff_components = Dict(col => Real[] for col = components)
avg_components = Dict{Symbol,Real}(col => 0 for col = components)
centered_components = Dict(col => Real[] for col = components)

n = size(df)[1] - apply_D_count

lag_column = Symbol("  lag k")
i_column = Symbol(" j-component")

output = Dict{Symbol,Array{Any,1}}(Symbol("rho^[D^($apply_D_count)["*String(j_col)*"],j](k)") => Real[] for j_col = components)

for j_col = components
    output[Symbol("p(rho[D^($apply_D_count)["*String(j_col)*"],j](k) nonzero)")] = Real[]
    output[Symbol(" p-sgn(rho[D^($apply_D_count)["*String(j_col)*"],j](k))")] = String[]
end

output[lag_column] = Int[]
output[i_column] = String[]

gamma_hat_i_i_at_zero = Dict{Symbol,Real}(col => 0 for col = components)

for c = components
    diff_components[c] = diff_series(df[!,c],apply_D_count)
    avg_components[c] = Statistics.mean(diff_components[c])
    centered_components[c] = [x - avg_components[c] for x = diff_components[c]]
end


for lag_k = lag_from:lag_to
    println("\n cross-correlation matrix for lag = $lag_k")
    for i_col = components
        append!(output[lag_column],[lag_k])
        append!(output[i_column],["D^($apply_D_count)["*String(i_col)*"]"])
        print("  ")
        for j_col = components

            # This is on page 23, yet in the book the denominator
            # does not seem to be right (sums over all n samples instead of the used n-k samples)

            sigma_j = sqrt(mean([x^2 for x = centered_components[j_col][1:(n-lag_k)]]))
            sigma_i = sqrt(mean([x^2 for x = centered_components[i_col][(1+lag_k):n]]))
            rho_hat = mean([(centered_components[j_col][t])*(centered_components[i_col][t+lag_k])
                for t in 1:(n-lag_k)
                ])/ (sigma_i*sigma_j)
                # rho^_{i,j}(k) = gamma^_{i,j}(k) / ((gamma^_{i,i}(0)^.5)(gamma^_{j,j}(0)^.5))
            append!(output[Symbol("rho^[D^($apply_D_count)["*String(i_col)*"],j](k)")],[
                rho_hat                
            ])
            #clamp the value of rho_hat
            if rho_hat > 1
                rho_hat = 1
            elseif rho_hat < -1
                rho_hat = -1
            end
            test_value = rho_hat * sqrt((n-lag_k-2)/(1-rho_hat^2))
            degree_of_freedom = n-lag_k-2
            tdist = Distributions.TDist(degree_of_freedom)
            
            calc_p = pvalue(tdist,test_value)
            
            append!(output[Symbol("p(rho[D^($apply_D_count)["*String(j_col)*"],j](k) nonzero)")],
            [calc_p])
            if calc_p < p_value
                if rho_hat < 0
                    append!(output[Symbol(" p-sgn(rho[D^($apply_D_count)["*String(j_col)*"],j](k))")],["-"])
                    print(" -")
                else
                    append!(output[Symbol(" p-sgn(rho[D^($apply_D_count)["*String(j_col)*"],j](k))")],["+"])
                    print(" +")
                end
            else
                append!(output[Symbol(" p-sgn(rho[D^($apply_D_count)["*String(j_col)*"],j](k))")],["."])
                print(" .")
            end
        end
        println()
    end
end

df_output = DataFrame(output)

println(df_output)

store = args["store"]
if store !== nothing
    CSV.write(args["store"], df_output)
end
