function [C,Cu,Cx] = OCP_GEN_C_Cu_Cx(in1,in2,in3,parIdx)
%OCP_GEN_C_CU_CX
%    [C,CU,CX] = OCP_GEN_C_CU_CX(IN1,IN2,IN3,PARIDX)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    06-Nov-2020 17:34:35

C = zeros(0,1);
if nargout > 1
    Cu = zeros(0,8);
end
if nargout > 2
    Cx = zeros(0,14);
end
