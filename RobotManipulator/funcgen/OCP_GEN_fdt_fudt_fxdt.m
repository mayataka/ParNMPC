function [fdt,fudt,fxdt] = OCP_GEN_fdt_fudt_fxdt(u,x,p,parIdx)
   [f,fu,fx] = f_fu_fx_Wrapper(u,x,p,parIdx);
   fdt = f*0.041667;
   fudt = fu*0.041667;
   fxdt = fx*0.041667;
end