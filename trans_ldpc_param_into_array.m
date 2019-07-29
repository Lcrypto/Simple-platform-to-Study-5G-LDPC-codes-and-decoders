% transform ldpc_param structure into array(for mex function)
function ldpc_param_array = trans_ldpc_param_into_array(ldpc_param)

ldpc_param_array = zeros(1,17);
ldpc_param_array(1) = ldpc_param.B;
ldpc_param_array(2) = floor(ldpc_param.code_rate*1024);
ldpc_param_array(3) = ldpc_param.K_b;
ldpc_param_array(4) = ldpc_param.K_p;
ldpc_param_array(5) = ldpc_param.K_n;
ldpc_param_array(6) = ldpc_param.C;
ldpc_param_array(7) = ldpc_param.L;
ldpc_param_array(8) = ldpc_param.iLS;
ldpc_param_array(9) = ldpc_param.BG_sel;
ldpc_param_array(10) = ldpc_param.Z_c;
ldpc_param_array(11) = ldpc_param.K;
ldpc_param_array(12) = ldpc_param.N;
ldpc_param_array(13) = size(ldpc_param.H_BG,2);
ldpc_param_array(14) = size(ldpc_param.H_BG,1);
ldpc_param_array(15) = size(ldpc_param.H_col,2) - 1;
ldpc_param_array(16) = size(ldpc_param.H_row,2) - 1;
ldpc_param_array(17) = ldpc_param.K_cb;

end