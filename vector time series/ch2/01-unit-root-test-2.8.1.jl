using HypothesisTests
using DataFrames
using CSV
using ArgParse

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

df = CSV.read(args["csv"])

println(names(df))

println("done.")
