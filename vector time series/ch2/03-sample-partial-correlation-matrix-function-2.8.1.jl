#!/usr/bin/env -S julia --sysimage=../sys_vts.so --

using HypothesisTests
using DataFrames
using CSV
using ArgParse
using Statistics
using Distributions
using GLM
using Plots


aps = ArgParseSettings(
    description = "This program computes the partial lag correlation matrices for a given " *
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
    default = 0.01
    
    "--diff"
    help = "number of times the difference operator D is applied to the input series"
    arg_type = Int
    default = 1
    
    "--center"
    help = "center every observation variable after differencing"
    arg_type = Bool
    default = true
    
    "--lag-from"
    help = "minimal lag k"
    arg_type = Int
    default = 1
    
    "--lag-to"
    help = "maximal lag k"
    arg_type = Int
    default = nothing
    
    "--store-matrix"
    help = "store the partial sample correlation matrix function table"
    arg_type = String
    required = false

    "--store"
    help = "store the test statistics table"
    arg_type = String
    required = false
    
    "--show-plot"
    help = "show a plot with the test statistics"
    arg_type = Bool
    default = false

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

function center_series(x)
    x .- mean(x)
end

###

"""
    (x, lag::Int)
    
    x __ data frame with value columns, rows must be in order of the time column (which must be omitted)
"""
function calculate_residual_vectors(x, lag::Int)

    # convert data frame to array ...
    x = Array{Float64}(x)

    # (1.6),(1.7) on p. 5 [beware that there is a typo in (1.7): when s=1, V_{s-1,t} = Z_t not Z_{t+1}]
    n_row,n_col = size(x)
    
    # calculate sizes
    rows = n_row - (lag+1)
            
    if (lag > 1)        
        # s>=2 involves removing linear dependence through regression
        
        X_cols = n_col * (lag-1)
        
        
        X = Array{Float64}(undef, rows, X_cols)
        
        # fill covariates
        
        for r in 1:rows
            for l in 1:(lag-1)
                X[r,1 + (l-1)*n_col : l*n_col] = x[r+l,:]
            end
        end
        
        # calculate residuals
        
        U = Array{Float64}(undef, rows, n_col)
        V = Array{Float64}(undef, rows, n_col)
        
        # each column has to be fitted for themselves
        for c in 1:n_col
                yu = x[1:rows,c]
                yv = x[1+lag:lag+rows,c]
                mu = GLM.fit(LinearModel, X, yu, true)
                mv = GLM.fit(LinearModel, X, yv, true)
                yhat_u = predict(mu,X)
                yhat_v = predict(mv,X)
                U[:,c] = yu - yhat_u
                V[:,c] = yv - yhat_v
        end
        
        return U,V
    
    elseif (lag == 1)
        # for s=1, the formula does not involve any regression at all.
        
        U = Array{Float64}(undef, rows, n_col)
        V = Array{Float64}(undef, rows, n_col)
        
        # each column has to be fitted for themselves
        for c in 1:n_col
                yu = x[1:rows,c]
                yv = x[1+lag:lag+rows,c]
                U[:,c] = yu
                V[:,c] = yv
        end
        
        return U,V
        
    end
end

###

"""
    (U,V)
    
    calculates the covariance matrix between U and V
    
    optional parameter:
      diagonal_inverted _ set all non-diagonal elements to zero insteard of calculating the expected values,
                          and calculate 1/(sqrt(mean)) for the diagonal elements.
                          e.g. calculates D_U from C_UU.
"""
function C_UV(U,V,diagonal_inverted::Bool=false)
    
    U_row,U_col = size(U)
    V_row,V_col = size(V)
    
    @assert(U_row == V_row)
    
    C = Array{Float64}(undef, U_col, V_col)
    
    for u in 1:U_col
        
        x = U[:,u]
        
        for v in 1:V_col
            if diagonal_inverted && (u != v)
                C[u,v] = 0
            else
                y = V[:,v]
                if diagonal_inverted
                    variance = mean(x .* y)
                    C[u,v] = 1/(sqrt(variance))
                else
                    C[u,v] = mean(x .* y)
                end
            end
        end
    end
    
    return C
end

"""
    (U,V)
    
    U,V  __ time series, each row corresponds to the series at a given time; must be ordered accordingly.
    
    determines the (sample) cross-correlation matrix between U and V
    
     -- as an exercise for the careful reader, prove or disprove why this computes the same matrix
        as C_UV(V,V,true) * C_UV(V,U) * C_UV(U,U,true) 
         [which still is probably more effective since we do not calculate the variances more than once there.]
"""
function sample_cor_matrix(U,V)
    
    U_row,U_col = size(U)
    V_row,V_col = size(V)
    
    @assert(U_row == V_row)
    
    R = Array{Float64}(undef, U_col, V_col)
    
    for u in 1:U_col
        for v in 1:V_col
            x = U[:,u]
            y = V[:,v]
            x_centered = x .- Statistics.mean(x)
            y_centered = y .- Statistics.mean(y)
            sigma_x = sqrt(mean([xt^2 for xt in x]))
            sigma_y = sqrt(mean([yt^2 for yt in y]))
            R[u,v] = mean(x .* y)/(sigma_x*sigma_y)
        end
    end
    
    return R
end

"""
    (P, p_value, n, lag_k)
    
    does a T-test whether the elements in the matrix P are indeed no-zero; and set the values to zero if otherwise.
    Should work for sample correlation matrices. But it's not what is done in the book.
"""
function zap_insignificant_to_zero(P, p_value, n, lag_k)
    degree_of_freedom = n-lag_k-2
    tdist = Distributions.TDist(degree_of_freedom)
    
    n_row,n_col = size(P)
    
    P_sig = Array{Float64}(0, n_row,n_col)
    
    for i in 1:n_row
        for j in 1:n_col
            rho_hat = P[i,j]
            #clamp the value of rho_hat
            if rho_hat > 1
                rho_hat = 1
            elseif rho_hat < -1
                rho_hat = -1
            end
            test_value = rho_hat * sqrt((n-lag_k-2)/(1-rho_hat^2))
        
            calc_p = pvalue(tdist,test_value)
    
            if calc_p < p_value
                P_sig[i,j] = P[i,j]
            end
        end
    end
    
    return P_sig
end

###

df = CSV.read(args["csv"])
time_col = Symbol(args["time-column"])
apply_D_count = args["diff"]

p_value = args["p-value"]

lag_from = args["lag-from"]
lag_to = args["lag-to"]

if lag_to === nothing
    lag_to = min(35,size(df)[1] - 3 - apply_D_count)
end



components = [col for col in names(df) if !(col == time_col)]

if args["center"]
    diff_components = DataFrame(Dict(col => center_series(diff_series(df[!,col],apply_D_count)) for col = components))
else
    diff_components = DataFrame(Dict(col => diff_series(df[!,col],apply_D_count) for col = components))
end


lag_column = Symbol("  lag k")
j_column = Symbol(" j-component")
i_column_names = [Symbol("P[D^($apply_D_count)["*String(j_col)*"],j](k)") for j_col in components]

output = Dict{Symbol,Array{Any,1}}(key => Real[] for key in i_column_names)

stats = Dict{Symbol,Array{Real,1}}(key => Real[] for key in [Symbol("lag"),Symbol("X(lag)"),
                            Symbol("X(lag) with insignificant entries")]) #,Symbol("XR(lag)")])

output[lag_column] = Int[]
output[j_column] = String[]

column_names = propertynames(diff_components)
n_columns, = size(column_names)

n_obs, = size(diff_components)

for lag_k = lag_from:lag_to
    global P,X,U,V,CVU,DU,DV,P0
    println("\n partial cross-correlation matrix for lag = $lag_k")
    U,V = calculate_residual_vectors(diff_components,lag_k)
    CVU = C_UV(V,U)
    DU = C_UV(U,U,true)
    DV = C_UV(V,V,true)
    
    #R = sample_cor_matrix(U,V)
    P0 = DV * CVU * DU
    X0 = sum(P0.^2)*(n_obs)
    
    P = DV * CVU * DU
    for j in 1:n_columns
        append!(output[lag_column],[lag_k])
        append!(output[j_column],["D^($apply_D_count)["*String(column_names[j])*"]"])
        for i in 1:n_columns
            append!(output[i_column_names[i]],[P[i,j]])
        end
    end
    
    
    # calculate the test statistics 
    X = sum(P.^2)*(n_obs-lag_k)
    #XR = sum(R.^2)*(n_obs-lag_k)
    append!(stats[Symbol("lag")],[lag_k])
    append!(stats[Symbol("X(lag)")],[X])
    append!(stats[Symbol("X(lag) with insignificant entries")],[X0])
    #append!(stats[Symbol("XR(lag)")],[XR])
end

df_stats = DataFrame(stats)

println(df_stats)

# optionally store the partial lag matrices
store = args["store-matrix"]
if store !== nothing
    df_output = DataFrame(output)

    CSV.write(args["store-matrix"], df_output)
end


if args["show-plot"]
    p = scatter(
        df_stats[!,Symbol("lag")],
        df_stats[!,Symbol("X(lag)")])
        
        
    display(p)
    gui(p)
    println("Showing plot, press enter to exit...")
    readline()

end