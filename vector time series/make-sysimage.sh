#!/bin/bash

cd $(dirname $0)

julia --sysimage=sys_vts.so requirements.jl || julia requirements.jl
julia --sysimage=sys_vts.so do-precompile.jl || julia do-precompile.jl
