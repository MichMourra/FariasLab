#!/bin/bash
# El siguiente comando puede ser utilizado para correr este programa y que no si interrumpa si se cierra sesion:
# nohup resto.sh &
    


python calc_elongation.py --gene d1dp4a_ --polymer ./output_input_ALL_PBPs/d1dp4a_._abs_energies.dat --output ./output_input_ALL_PBPs &

python calc_elongation.py --gene d1ewka_ --polymer ./output_input_ALL_PBPs/d1ewka_._abs_energies.dat --output ./output_input_ALL_PBPs &

python calc_elongation.py --gene d3hsya_ --polymer ./output_input_ALL_PBPs/d3hsya_._abs_energies.dat --output ./output_input_ALL_PBPs &

wait
