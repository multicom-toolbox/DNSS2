package Time;

# Module: 	Time.pm
# Author: 	Matt Spencer
# Made:		~4/1/13
# Last Mod:	4/16/13
#
# This module calls the time and can return it in a variety of 
# convenient forms.
#
# Functions:
#	get_localtime : gets array of all time
#	get_time : gets array of time
#	get_date : gets array of date
#	formatted_time : gives string of time
#	formatted_date : gives string of date
#	formatted_localtime: gives string of date and time
#
# Dependencies:
#	none

use strict;
use Exporter;
use vars qw(@ISA @EXPORT @EXPORT_OK);

@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw(get_localtime get_time get_date formatted_time format_time formatted_date format_date formatted_localtime difference time_to_secs secs_to_time);

##############################################################
# Name : get_localtime
# Takes: nothing
# Returns: array containing sec, min, hr, day, mon, year
##############################################################

sub get_localtime {
	my ($sec, $min, $hr, $day, $mon, $year) = localtime(time);
	$mon += 1;
	$year += 1900;

	return ($sec, $min, $hr, $day, $mon, $year);
}

##############################################################
# Name : get_time
# Takes: nothing
# Returns: array containing hour, min, sec
##############################################################

sub get_time {
	my @time = get_localtime();
	return ($time[2], $time[1], $time[0]);
}

##############################################################
# Name : formatted_time
# Takes: nothing
# Returns: string containing hour, min, sec
##############################################################

sub formatted_time {
	my ($hr, $min, $sec) = get_time();
	my $time = sprintf("%02d:%02d:%02d", $hr, $min, $sec);
	return $time;
}

##############################################################
# Name : format_time
# Takes: array time
# Returns: string containing hour, min, sec
##############################################################

sub format_time {
	my ($hr, $min, $sec) = @_;
	my $time = sprintf("%02d:%02d:%02d", $hr, $min, $sec);
	return $time;
}

##############################################################
# Name : get_date
# Takes: nothing
# Returns: array containing day, month, year
##############################################################

sub get_date {
	my @time = get_localtime();
	return ($time[3], $time[4], $time[5]);
}

##############################################################
# Name : formatted_date
# Takes: nothing
# Returns: string date
##############################################################

sub formatted_date {
	my ($day, $mon, $year) = get_date();
	my $date = sprintf("%02d/%02d/%04d", $mon, $day, $year);
	return $date;
}

##############################################################
# Name : format_date
# Takes: array date
# Returns: string date
##############################################################

sub format_date {
	my ($day, $mon, $year) = @_;
	my $date = sprintf("%02d/%02d/%04d", $mon, $day, $year);
	return $date;
}

##############################################################
# Name : formatted_localtime
# Takes: nothing
# Returns: string date and time
##############################################################

sub formatted_localtime {
	my $date = formatted_date();
	my $time = formatted_time();

	return "[$date $time] ";
}

sub difference {
	my ($time_ref1, $time_ref2) = @_;
	my @time1 = @{ $time_ref1 };
	my @time2 = @{ $time_ref2 };

	my $secs1 = time_to_secs(@time1);
	my $secs2 = time_to_secs(@time2);
	
	my $diffsecs;
	$diffsecs = $secs2-$secs1 if ($secs2>$secs1);
	$diffsecs = $secs1-$secs2 if ($secs1>=$secs2);
	my @diff = secs_to_time($diffsecs);
	return @diff;
}

sub time_to_secs {
	my @time = @_;
	
	for (my $ii=0; $ii<$#time; $ii++){
		$time[$ii+1] += $time[$ii]*60;
	}
	return $time[$#time];
}

sub secs_to_time {
	my ($secs) = @_;
	my @time = (0, 0, $secs);

	for (my $ii=2; $ii>=0; $ii--){
		if ($time[$ii]>=60){
			$time[$ii-1]+= int($time[$ii]/60);
			$time[$ii] = $time[$ii]%60;
		}
	}
	return @time;
}


1;
