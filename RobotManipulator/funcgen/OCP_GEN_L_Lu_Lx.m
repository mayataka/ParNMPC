function [L,Lu,Lx] = OCP_GEN_L_Lu_Lx(in1,in2,in3)
%OCP_GEN_L_LU_LX
%    [L,LU,LX] = OCP_GEN_L_LU_LX(IN1,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    06-Nov-2020 17:34:34

p1 = in3(1,:);
p2 = in3(2,:);
p3 = in3(3,:);
p4 = in3(4,:);
p5 = in3(5,:);
p6 = in3(6,:);
p7 = in3(7,:);
u1 = in1(1,:);
u2 = in1(2,:);
u3 = in1(3,:);
u4 = in1(4,:);
u5 = in1(5,:);
u6 = in1(6,:);
u7 = in1(7,:);
u8 = in1(8,:);
x1 = in2(1,:);
x2 = in2(2,:);
x3 = in2(3,:);
x4 = in2(4,:);
x5 = in2(5,:);
x6 = in2(6,:);
x7 = in2(7,:);
x8 = in2(8,:);
x9 = in2(9,:);
x10 = in2(10,:);
x11 = in2(11,:);
x12 = in2(12,:);
x13 = in2(13,:);
x14 = in2(14,:);
t2 = p1-x1;
t4 = p2-x2;
t6 = p3-x3;
t8 = p4-x4;
t10 = p5-x5;
t12 = p6-x6;
t14 = p7-x7;
L = u1.^2./2.0e+3+u2.^2./2.0e+3+u3.^2./2.0e+3+u4.^2./2.0e+3+u5.^2./2.0e+3+u6.^2./2.0e+3+u7.^2./2.0e+3+u8.^2.*1.0e+3+x8.^2./2.0e+1+x9.^2./2.0e+1+x10.^2./2.0e+1+x11.^2./2.0e+1+x12.^2./2.0e+1+x13.^2./2.0e+1+x14.^2./2.0e+1+t2.*(p1./2.0-x1./2.0)+t4.*(p2./2.0-x2./2.0)+t6.*(p3./2.0-x3./2.0)+t8.*(p4./2.0-x4./2.0)+t10.*(p5./2.0-x5./2.0)+t12.*(p6./2.0-x6./2.0)+t14.*(p7./2.0-x7./2.0);
if nargout > 1
    Lu = [u1./1.0e+3,u2./1.0e+3,u3./1.0e+3,u4./1.0e+3,u5./1.0e+3,u6./1.0e+3,u7./1.0e+3,u8.*2.0e+3];
end
if nargout > 2
    t3 = -t2;
    t5 = -t4;
    t7 = -t6;
    t9 = -t8;
    t11 = -t10;
    t13 = -t12;
    t15 = -t14;
    Lx = [t3,t5,t7,t9,t11,t13,t15,x8./1.0e+1,x9./1.0e+1,x10./1.0e+1,x11./1.0e+1,x12./1.0e+1,x13./1.0e+1,x14./1.0e+1];
end