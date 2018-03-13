  phi=function(x,t){
#
# Funtion phi computes the exact solution of Burgers' equation
# for comparison with the numerical solution.  It is also used to
# define the initial and boundary conditions for the numerical 
# solution.
#
# Analytical solution
  a=(0.05/vis)*(x-0.5+4.95*t);
  b=(0.25/vis)*(x-0.5+0.75*t);
  c=( 0.5/vis)*(x-0.375);
  ea=exp(-a);
  eb=exp(-b);
  ec=exp(-c);
  ua=(0.1*ea+0.5*eb+ec)/(ea+eb+ec);
  return(c(ua));
}