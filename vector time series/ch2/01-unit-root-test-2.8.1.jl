using HypothesisTests
using DataFrames
using CSV
using ArgParse
using Plots

aps = ArgParseSettings(
    description = "This program does stuff. " * #TODO: add precise description
                  "We suppose that the time series is" *
                  " presented in order from earlier to later," *
                  " and that the intervals between" *
                  " each time point are equal.",
)
@add_arg_table! aps begin
    # "--opt1"
    #     help = "an option with an argument"
    # "--opt2", "-o"
    #     help = "another option with an argument"
    #     arg_type = Int
    #     default = 0
    # "--flag1"
    #     help = "an option without argument, i.e. a flag"
    #     action = :store_true
    "--time-column"
        help = "name of the column containing the time information"
        arg_type = String
        default = "Period"
    "--p-value"
        help = "p-value for univariate Augmented Dickey-Fuller unit root test"
        arg_type = Real
        default = .01
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



println(names(df))

print( HypothesisTests.ADFTest(df[!, :Automobile], :none, 1) )

println("done.")
