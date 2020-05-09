#!/bin/zsh

res=$1
wfn=$2

julia anim.jl res${res}_wfn${wfn}.txt | anim
