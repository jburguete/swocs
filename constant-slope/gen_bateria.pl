my $L=$ARGV[0];
my $S0=$ARGV[1];
my $Q=$ARGV[2];
my $t=$ARGV[3];
my $nm=$ARGV[4]; #numerical model
my $pm=$ARGV[5]; #physical model
my $z0;

$z0=$S0*$L; #zfinal es siempre cero

printf "%g 1 0 1\n", $L; 

print <<END1;
1
1 1 1
2
END1

printf "0 %g\n", $z0; 
printf "%g 0\n", $L; 

print <<END2;
0.03
0 1 0 1
10
1
END2

printf "0 %g\n1\n0 %g\n501 1\n", $Q,$Q; 

printf "%g 0 0.9 0.01\n", $t; 

printf "%d 2\n", $nm; 
printf "%d\n\n", $pm; 

print <<END3;
(channel length) (bottom width) (wall slope) (channel depth)
(outlet type)
(friction model) (infiltration model) (diffusion model)
(points number of geometry)
(x-coordinate of the geometry point) (z-coordinate of the geometry point)
...
(friction coefficients) ...
(infiltration coefficients) ...
(diffusion coefficients) ...
(points number of the inlet water hydrogam)
(time of the inlet water hydrogam) (discharge of the inlet water hydrogram)
...
(points number of the inlet solute hydrogam)
(time of the inlet solute hydrogam) (discharge of the inlet solute hydrogram)
...
(number of mesh nodes) (initial conditions type)
(initial conditions parameters)
...
(final time) (interval time of measures) (cfl number) (minimum_depth)
(numerical model of surface flow) (numerical model of diffusion)
(model of surface flow)
END3
