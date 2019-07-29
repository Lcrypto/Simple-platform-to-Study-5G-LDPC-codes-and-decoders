function [crc_bits] = nr15_crc24b_encode(input_bits)

L = 24;

crc_bits = zeros(1,L);
