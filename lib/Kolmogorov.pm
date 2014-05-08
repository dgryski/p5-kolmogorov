package Kolmogorov;

use warnings;
use strict;

use Exporter;
our @ISA= qw( Exporter );
our @EXPORT_OK= qw(
        KS
);
our @EXPORT= qw(
        KS
);

sub KS {

        my ($data1, $data2) = @_;

        my ($n1, $n2) = (scalar @$data1, scalar @$data2);

        my ($en1, $en2) = ($n1+0.0, $n2+0.0);

        my ($d, $fn1, $fn2) = (0,0,0);

        @$data1 = sort {$a <=> $b} @$data1;
        @$data2 = sort {$a <=> $b} @$data2;

        my $j1 = 0;
        my $j2 = 0;

        while ($j1 < $n1 && $j2 < $n2) {
                my $d1 = $data1->[$j1];
                my $d2 = $data2->[$j2];

                if ($d1 <= $d2) {
                        while ($j1 < $n1 && $d1 == $data1->[$j1]) {
                                $j1++;
                                $fn1 = $j1/$en1;
                        }
                }

                if ($d2 <= $d1) {
                        while ($j2 < $n2 && $d2 == $data2->[$j2]) {
                                $j2++;
                                $fn2 = $j2/$en2;
                        }
                }

                my $dt = abs($fn2-$fn1);
                if ($dt > $d) { $d = $dt };
        }

        my $en = sqrt(($n1*$n2)/($n1+$n2));
        # R and Octave don't use this approximation that NR does
        # return qks((en + 0.12 + 0.11/en) * d)
        return qks($en * $d);
}

sub qks {
        my $z = shift;

       die "bad z in qks" if $z < 0;

        return 1 if $z == 1;

        return 1 - pks($z) if $z < 1.18;

        my $x = exp(-2 * ($z*$z));

        my $x4 = ($x*$x*$x*$x);
        my $x9 = $x4*$x4*$x;

        return 2.0 * ($x - $x4 + $x9);
}


sub pks {
        my $z = shift;

        die "bad z in pks" if $z < 0;

        return 0 if $z == 0;

        if ($z < 1.18) {
                my $y = exp(-1.23370055013616983 / ($z * $z));

                my $y4 = $y*$y*$y*$y;
                my $y8 = $y4*$y4;
                my $y16 = $y8*$y8;
                my $y32 = $y16*$y16;

                my $y9 = $y8 * $y;
                my $y25 = $y16 * $y9;
                my $y49 = $y32 * $y16 * $y;

                return 2.25675833419102515 * sqrt(-log($y)) * ($y + $y9 + $y25 + $y49);
        } else {
                my $x = exp(-2 * ($z*$z));

                my $x4 = ($x*$x*$x*$x);
                my $x9 = $x4*$x4*$x;

                return 1.0 - 2.0 * ($x - $x4 + $x9);
        }
}

1;

