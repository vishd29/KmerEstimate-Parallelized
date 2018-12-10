Compile and Run
------------------------------
Compile:


		g++ -o pkmerEst parallelKmerCountEstimate.cpp -std=c++11 -O3 -march=native

Run:

		./pkmerEst -f <seq.fa> -k  <kmerLen> -s <minHeap_Size> -c <coverage> -t <thread_count> -q <queue_size> -o <out.txt>