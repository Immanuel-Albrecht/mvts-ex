#!/usr/bin/env -S julia --sysimage=../sys_vts.so --

using HypothesisTests
using DataFrames
using CSV
using ArgParse

aps = ArgParseSettings(
    description = "This program performs Augmented Dickey-Fuller tests " *
                  "on each column of a given data set with increasing " *
                  "number of self-differences of the series, until " *
                  "the all differenced series' may be considered unit root free. " *
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
        help = "p-value for univariate Augmented Dickey-Fuller unit root test"
        arg_type = Real
        default = .01
    "--lag"
        help = "lag used to calculate the Augmented Dicker-Fuller unit root test"
        arg_type = Int
        default = 1
    "--store"
        help = "store the univariate p-value table"
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

function diff_series(x, how_often::Int=1)
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
target_p = args["p-value"]
lag = args["lag"]

dp = Dict(col => Real[] for col in names(df) if ! (col == time_col))
d_col = Symbol(" D^( ... )")
dp[d_col] = []

println("Accepted p-value <$target_p.")

for d in 0:(size(df)[1]-lag)
    nonzero_count = 0
    append!(dp[d_col], d)
    for col in names(df)
        if ! (col == time_col)
            y = diff_series(df[!,col],d)
            p = pvalue(ADFTest(y, :none, lag))
            significant = (p<target_p) ? " *" : ""
            append!(dp[col], p)
            if p >= target_p
                nonzero_count += 1
            end
        end
    end
    if nonzero_count < 1
        break
    end
end

df_p = DataFrame(dp)

println(df_p)

store = args["store"]
if store != nothing
    CSV.write(args["store"], df_p)
end
