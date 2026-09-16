function tau = calculateTau(h,mu,u)

    tau = 1/sqrt( (2*u/h)^2 + ((4*mu)/(h*h))^2 );

    % tau = (h*h)/(4*mu);

end