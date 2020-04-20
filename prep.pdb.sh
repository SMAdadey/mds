#!/usr/bin/env bash

for i in $(cat covid19.acc); do echo "-c -d \'load ${i/_*/.pdb}; indicate chain ${i/*_/}; save ${i/_/_chain}.pdb\'"; done | parallel echo pymol {} | xargs -I input -P10 sh -c "echo input"
