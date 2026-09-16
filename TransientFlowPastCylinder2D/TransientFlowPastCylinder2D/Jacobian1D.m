function J = Jacobian1D(x,xi,Ordr)
   J = dNdxi1D(xi,Ordr) * x;
end