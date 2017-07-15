wannaAln: edlib.cpp edlib.h kseq.h revcomp.h wannaAln.c
	g++ wannaAln.c edlib.cpp -o wannaAln -lz

clean:
	rm -f wannaAln
