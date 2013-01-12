my $B0=$ARGV[0];
my $Z=$ARGV[1];
my $S0=$ARGV[2];
my $r=$ARGV[3];
my $h1=$ARGV[4];
my $Q1=$ARGV[5];
my $h2=$ARGV[6];
my $Q2=$ARGV[7];
my $cfl=$ARGV[8];
my $nm=$ARGV[9]; #numerical model
my $pm=$ARGV[10]; #physical model
my $L;
my $z0;
my $A1;
my $A2;
my $t;
my $U;

$L=100;
$z0=$S0*$L; #zfinal es siempre cero
$A1=($B0+$Z*$h1)*$h1;
$A2=($B0+$Z*$h2)*$h2;
$U=($Q1-$Q2)/($A1-$A2);
$t=0.5*$L/$U;

printf "%.17lg %.17lg %.17lg 1\n", $L, $B0, $Z; 

print <<END1;
2 2
1 1 1
2
END1

printf "0 %.17lg\n", $z0; 
printf "%.17lg 0\n", $L; 
printf "%.17lg\n", $r;

print <<END2;
0 1 0 1
0
END2

printf "1\n0 %.17lg\n1\n0 %.17lg\n101 2\n4\n", $Q1, $Q1; 

printf "0 %.17lg %.17lg 1\n", $A1, $Q1; 
printf "%.17lg %.17lg %.17lg 1\n", 0.25*$L, $A1, $Q1, $A1; 
printf "%.17lg %.17lg %.17lg 0\n", 0.25*$L, $A2, $Q2; 
printf "%.17lg %.17lg %.17lg 0\n", $L, $A2, $Q2; 

printf "%.17lg 0 %.17lg 0.01\n", $t, $cfl; 

printf "%d 2\n", $nm; 
printf "%d\n\n", $pm; 

print <<END3;
(channel length) (bottom width) (wall slope) (channel depth)
(inlet tupe) (outlet type)
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
