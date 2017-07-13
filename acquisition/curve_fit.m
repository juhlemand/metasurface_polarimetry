g = fittype(@(I0, A, phi, x) I0+A*cos((pi/180)*(x+phi)).^2);
X = p1_new(:, 1);
Y = p1_new(:, 2);
plot(X, Y)
fit_params = fit(X, Y, g);
hold on
plot(fit_params)
hold off