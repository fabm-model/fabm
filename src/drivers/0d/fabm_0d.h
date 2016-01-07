! Taken from the General Ocean Turbulence Model:

#define PATH_MAX	255

#define stderr		0
#define stdout		6

#define STDOUT write(stdout,*)
#define STDERR write(stderr,*)
#define LEVEL0 STDERR
#define LEVEL1 STDERR '   ',
#define LEVEL2 STDERR '       ',
#define LEVEL3 STDERR '           ',
#define LEVEL4 STDERR '               ',

#define LINE "------------------------------------------------------------------------"
