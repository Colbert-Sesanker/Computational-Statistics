# Question 1 in Octave (runs in Matlab)
# Colbert Sesanker HW3 STA 613

# draw normal and poission samples
norm_samples = normrnd(0 ,1, 500,1);
theta = [-1,.5,2];
for i = 1:3
  pois_samples(:,i) = poissrnd(exp(theta(i)*norm_samples))
end

# plots
plot(pois_samples(:,1), norm_samples)
plot(pois_samples(:,2), norm_samples)
plot(pois_samples(:,3), norm_samples)
plot(log(pois_samples(:,1)), norm_samples)

# Newton's method to estimate theta from samples 1b
it = 20;
t = zeros(it,3);
t(1,:) = [1,1,1]; # initialize, t is the theta estimate
for i= 1:3
  y = pois_samples(:,i);
  for j=1:it-1
   y_ = exp(t(j,i)*norm_samples);
   h = norm_samples'*(y-y_); # scalar gradient
   h_prime = -(norm_samples.^2)'*exp(t(j,i)*norm_samples); # scalar hessian
   t(j+1,i) = t(j,i) - h/h_prime; 
  end
end

# 1d
for i=1:3
 plot(t(:,i)); hold on;
end
print convergence_of_estimates
%{ it takes ABOUT TEN estimates to converge, notice how it goes out then comes back in. It converges pretty close to the actual theta values
%}

