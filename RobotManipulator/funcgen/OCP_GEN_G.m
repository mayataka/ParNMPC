function G = OCP_GEN_G(in1,in2,in3)
%OCP_GEN_G
%    G = OCP_GEN_G(IN1,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    06-Nov-2020 17:34:34

u1 = in1(1,:);
u2 = in1(2,:);
u3 = in1(3,:);
u4 = in1(4,:);
u5 = in1(5,:);
u6 = in1(6,:);
u7 = in1(7,:);
u8 = in1(8,:);
x8 = in2(8,:);
x9 = in2(9,:);
x10 = in2(10,:);
x11 = in2(11,:);
x12 = in2(12,:);
x13 = in2(13,:);
x14 = in2(14,:);
t2 = pi./2.0;
G = [u1+1.0e+1;u2+1.0e+1;u3+1.0e+1;u4+1.0e+1;u5+1.0e+1;u6+1.0e+1;u7+1.0e+1;-u1+1.0e+1;-u2+1.0e+1;-u3+1.0e+1;-u4+1.0e+1;-u5+1.0e+1;-u6+1.0e+1;-u7+1.0e+1;u8;t2+u8+x8;t2+u8+x9;t2+u8+x10;t2+u8+x11;t2+u8+x12;t2+u8+x13;t2+u8+x14;t2+u8-x8;t2+u8-x9;t2+u8-x10;t2+u8-x11;t2+u8-x12;t2+u8-x13;t2+u8-x14];
