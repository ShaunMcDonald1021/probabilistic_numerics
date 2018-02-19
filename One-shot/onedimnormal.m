function [f, mle] = onedimnormal(theta)

f = @(t) normpdf(t,theta(1), theta(2));
mle = theta(1);




