function y = fit_cross_polarizers(A, om, phi, x)
    y = A*sin(om*x+phi)*sin(om*x+phi);
end