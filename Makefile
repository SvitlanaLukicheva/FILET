CC = gcc
CFLAGS = -O3 -Wall -lm -lgsl -lgslcblas
BASICS = miscCode/stringWrap.c miscCode/sequenceMatrix.c pgSummaryStats.c miscCode/vector.c miscCode/numerical.c miscCode/nrutil.c miscCode/bedFile.c

all: twoPopnStats_forML pgStatsBedSubpop_forML

twoPopnStats_forML:	twoPopnStats_forML.c msGeneralStats.c
		$(CC) twoPopnStats_forML.c msGeneralStats.c -o twoPopnStats_forML $(CFLAGS)

pgStatsBedSubpop_forML:	pgStatsBedSubpop_forML_printNumSites.c $(BASICS)
		$(CC) pgStatsBedSubpop_forML_printNumSites.c $(BASICS) -o pgStatsBedSubpop_forML $(CFLAGS)

pgStatsBedOnePop_forGhostML:	pgStatsBedOnePop_forGhostML_printNumSites.c $(BASICS)
		$(CC) pgStatsBedOnePop_forGhostML_printNumSites.c $(BASICS) -o pgStatsBedOnePop_forGhostML $(CFLAGS)

onePopnStats_forGhostIntroML:	onePopnStats_forGhostIntroML.c msGeneralStats.c
		$(CC) onePopnStats_forGhostIntroML.c msGeneralStats.c -o onePopnStats_forGhostIntroML $(CFLAGS) 	

clean:
	rm pgStatsBedSubpop_forML twoPopnStats_forML
