all: pageRankGSOMP pageRankGS

pageRankGSOMP:pageRankGSOMP.c
	gcc pageRankGSOMP.c -o prgsomp -fopenmp -O3

pageRankGS:pageRankGS.c
	gcc pageRankGS.c -o prgs -O3

clear:
	rm prgs prgsOMP
