#!/usr/bin/env -S julia --sysimage=../sys_vts.so --

using DataFrames
using CSV
using ArgParse
using Plots

aps = ArgParseSettings(
    description = "This program plots vector time series. " *
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
    "--one-plot-per-component"
        help = "create a plot for each of the component series found in the data, as opposed to plotting" *
               " all series into a common plot"
        action = :store_true
    "--D"
        help = "number of times the time series are differenced before plotting"
        arg_type = Int
        default = 0
    "--gr"
        help = "use gr backend instead of plotly backend."
        action = :store_true
    "csv"
        help = "path to the CSV containing the vector time series"
        required = false
        arg_type = String
        default = "WW2a.csv"
end

args = parse_args(aps)

if ! args["gr"]
    plotly()
end

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

println(args)

df = CSV.read(args["csv"])

d = args["D"]
time_col = Symbol(args["time-column"])
new_plot_for_each = args["one-plot-per-component"]
t = df[!, time_col][1+d:size(df)[1]]

p = nothing


for col in names(df)
    if col != time_col
        title_s = if d > 0 
            "$d x Differenced Series"
        else 
            "Plain Series"
        end
        if new_plot_for_each 
            p = plot(t, diff_series(df[!, col],d), title=title_s * " for $col",legend=:none)
            display(p)
            gui(p)
            println("Showing $title_s for $col. Enter to continue...")
            readline()
        else
            global p
            if p == nothing
                p = plot(t, diff_series(df[!, col],d), title=title_s, label=col, legend=:outerbottomright)
            else
                plot!(p, diff_series(df[!, col],d), label=col)
            end
        end
    end
end

if ! new_plot_for_each
    display(p)
    gui(p)
    println("Done, enter to close.")
    readline()    
end
