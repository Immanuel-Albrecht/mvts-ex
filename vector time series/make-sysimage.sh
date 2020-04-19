#!/bin/bash

cd $(dirname $0)

julia requirements.jl
julia do-precompile.jl
