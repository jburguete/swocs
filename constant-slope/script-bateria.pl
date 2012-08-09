my $S0,$Q;
my $nm,$pm;
my $L=$ARGV[0];
my $t=$ARGV[1];


$nm=1;
$pm=3;
for($S0=0.00001;$S0<1.;$S0=$S0*10.)
{
	for($Q=0.001;$Q<100.;$Q=$Q*10.)
	{
		system("perl gen_bateria.pl $L $S0 $Q $t $nm $pm > i-$nm-$pm-$t-$S0-$Q");
		printf "i-$nm-$pm-$t-$S0-$Q\n";
		system("../swocs i-$nm-$pm-$t-$S0-$Q o-$nm-$pm-$t-$S0-$Q");
	}
}

$pm=4;
for($S0=0.00001;$S0<1.;$S0=$S0*10.)
{
	for($Q=0.001;$Q<100.;$Q=$Q*10.)
	{
		system("perl gen_bateria.pl $L $S0 $Q $t $nm $pm > i-$nm-$pm-$t-$S0-$Q");	
		printf "i-$nm-$pm-$t-$S0-$Q\n";
		system("../swocs i-$nm-$pm-$t-$S0-$Q o-$nm-$pm-$t-$S0-$Q");
	}
}

$nm=2;
$pm=1;
for($S0=0.00001;$S0<1.;$S0=$S0*10.)
{
	for($Q=0.001;$Q<100.;$Q=$Q*10.)
	{
		system("perl gen_bateria.pl $L $S0 $Q $t $nm $pm > i-$nm-$pm-$t-$S0-$Q");	
		printf "i-$nm-$pm-$t-$S0-$Q\n";
		system("../swocs i-$nm-$pm-$t-$S0-$Q o-$nm-$pm-$t-$S0-$Q");
	}
}

$pm=2;
for($S0=0.00001;$S0<1.;$S0=$S0*10.)
{
	for($Q=0.001;$Q<100.;$Q=$Q*10.)
	{
		system("perl gen_bateria.pl $L $S0 $Q $t $nm $pm > i-$nm-$pm-$t-$S0-$Q");	
		printf "i-$nm-$pm-$t-$S0-$Q\n";
		system("../swocs i-$nm-$pm-$t-$S0-$Q o-$nm-$pm-$t-$S0-$Q");
	}
}
