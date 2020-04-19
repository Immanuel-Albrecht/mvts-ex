using PackageCompiler

PackageCompiler.create_sysimage(
    [:Plots, :HypothesisTests, :DataFrames, :ArgParse, :CSV],
    sysimage_path = "sys_vts.so",
    precompile_execution_file = "precompile-image.jl",
)
