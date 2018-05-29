package OSPath;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
$VERSION = 1.00; # Or higher
@ISA = qw(Exporter);
@EXPORT = qw(convert_cygwin2perl_win_path convert_win2perl_win_path convert_win2cygwin_path get_newline convert_cygwin2win_path); # Symbols to autoexport (:DEFAULT tag)

sub convert_cygwin2perl_win_path{
	my $path=shift;
	if($path=~m!\/cygdrive\/([a-z])!i){
		my $drive=$1;
		$path=~s!\/cygdrive\/[a-z]!$drive:!i;
	};
	return $path;
};

sub convert_cygwin2win_path{
	my $path=shift;
	if($path=~m!\/cygdrive\/([a-z])!i){
		my $drive=$1;
		$path=~s!\/cygdrive\/[a-z]!$drive:!i;
	};
	$path=~s!\/!\\!g;
	return $path;
};

sub convert_win2perl_win_path{
	my $path=shift;
	$path=~s!\\!\/!g;
	return $path;
};

sub convert_win2cygwin_path{
	my $path=shift;
	$path=~s!\\!\/!g;
	$path=~s!([a-z]):!\/cygdrive\/\1!i;
	return $path;
};

sub get_newline{
	my $ifile=shift;
	open IFILE, "< $ifile" or die "\nUnable to open input file: $ifile";
	my ($nr,$nn,$N)=(0,0,0);
	while(<IFILE>){
		if(/\r\n$/){
			$nr++;
			$nn++;
		}elsif(/\n$/){
			$nn++;
		}else{
			while(m/\r/g){$nr++;};
		};
		$N++;
	};
	close IFILE;
	if($nn==$nr){
		return "\r\n";
	}elsif($nn!=$N){
		return "\r" if $nr>0&&$N==1;
		return undef;
	};
	return $/;
};

1;