
plot(function(x) dlaplace(x, 0, 1), -4, 4)

plot(function(x) plaplace(x, 0, 1), -4, 4)

x <- rlaplace(1e5, 5, pi)
hist(x)
hist(plaplace(x, 5, pi))
