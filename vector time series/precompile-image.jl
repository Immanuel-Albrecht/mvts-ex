

using HypothesisTests
using DataFrames
using CSV
using ArgParse
using Plots
using JuliaFormatter
using StatsBase
using Pkg
using Distributions

using TimeSeries
using GLM
using MixedModels
using StatsFuns
using StatsModels
using Distances
using Clustering
using KernelDensity
using MultivariateStats
using Survival

Distributions.Normal(2,2)

p = plot(rand(5), rand(5))
display(p)
plotly()
p = plot(rand(5), rand(5))
display(p)
gui(p)

df = CSV.read("../data/WW2a.csv")

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
    "csv"
    help = "path to the CSV containing the vector time series"
    required = false
    arg_type = String
    default = "WW2a.csv"
end

args = parse_args(aps)

format_file("format-dummy.jl")

println(names(df))

print( HypothesisTests.ADFTest(df[!, :Automobile], :none, 1) )

print( HypothesisTests.ADFTest(df[!, :Automobile], :trend, 2) )
