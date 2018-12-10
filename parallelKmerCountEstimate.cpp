#include <iostream>
#include <climits>
#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include <time.h>
#include "metrohash64.cpp"
#include <stdint.h>
#include <unordered_map>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <cassert>
#include <string.h>
#include <cstring>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <math.h>
#include <sys/time.h>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <list>
#include <stack>
#include <limits.h>
#include <map>
#include <bitset>
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <cstring>
#include <iostream>
#include <random>
#include <cinttypes>
#include <thread>
#include <atomic>
#include <queue>
#include "ntHashIterator.hpp"
#include "sparsepp/spp.h"

using spp::sparse_hash_map;
using namespace std;

#define MAX_THREADPOOL 1000

//global variables
bool eof;


unsigned trailing_zeros(uint64_t n) {
    return n ? __builtin_ctzll(n) : -1;
}

typedef sparse_hash_map<uint64_t, uint32_t> SMap;

using namespace std;
KSEQ_INIT(int, read)


void help() {

    cout << "KmerEst [options] -f <fasta/fastq> -k <k-mer length>  -s <sample size> -o <output file>" << endl
         << "  -h               help" << endl
         << "  -f <file>       Input sequence file " << endl
         << "  -k <k-mer size >        kmer size (default 31) " << endl
         << "  -s <sample size>        sample size (default 25m)" << endl
         << "  -c coverage>       coverage (default 64)" << endl
         << "  -o         	  Prefix of the Output file " << endl
         << "  -t         	  Number of threads required " << endl
         << "  -q         	  Queue size " << endl;

    exit(0);
}


void consume(vector<SMap> &vmap, queue<char *> &Q, mutex &mtx, condition_variable &cond, int k_size, int L,
             int &th, uint64_t &no_kmers) {

    int count = 0;
    uint64_t hash;
    uint8_t tz;

    while (!eof || !Q.empty()) {
        unique_lock<mutex> lck(mtx);
        while (Q.empty()) {
            cond.wait(lck);
            if (eof)
                break;
        }
        while (!Q.empty()) {
            ntHashIterator itr(Q.front(), 1, k_size);
            while (itr != itr.end()) {
                hash = (*itr)[0];
                no_kmers++;
                tz = trailing_zeros(hash);
                if (tz >= th) {
                    if (vmap[tz].find(hash) != vmap[tz].end())
                        vmap[tz][hash] += 1;
                    else { //// insert if not there
                        vmap[tz].insert(make_pair(hash, 1));
                        ++count;  // insert if not there
                        if (count == L) {
                            int cnt = vmap[th].size();
                            count = count - cnt;
                            SMap().swap(vmap[th]);
                            ++th;
                        }
                    }
                }
                ++itr;
            }
            Q.pop();
        }
        cond.notify_one();
    }
}

int main(int argc, char **argv) {
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

    if (argc == 1) {
        cout << argv[0] << " -f <seq.fa> -k  <kmerLen> -s <minHeap_Size> -c <coverage> -o <out.txt>" << endl;
        exit(0);
    }

    //making all required variable local.
    int k = 31;
    int s = 25000000;
    int thread_pool_size = 4;
    int c = 64;
    string in_file, out_file;
    int q_size = 10;



    //Read input
    for (int i = 1; i < argc; i++) {

        if (!strcmp(argv[i], "-h")) { help(); }
        else if (!strcmp(argv[i], "-k")) {
            k = atoi(argv[i + 1]);
            i++;
        } else if (!strcmp(argv[i], "-f")) {
            in_file = argv[i + 1];
            i++;
        } else if (!strcmp(argv[i], "-s")) {
            s = atoi(argv[i + 1]);
            //k = s / thread_count;
            i++;
        } else if (!strcmp(argv[i], "-c")) {
            c = atoi(argv[i + 1]);
            i++;
        } else if (!strcmp(argv[i], "-o")) {
            out_file = argv[i + 1];
            i++;
        } else if (!strcmp(argv[i], "-t")) {
            thread_pool_size = atoi(argv[i + 1]);
            i++;
        }
        else if (!strcmp(argv[i], "-q")) {
            q_size = atoi(argv[i + 1]);
            i++;
        }
    }

    if (in_file.empty() || out_file.empty()) {
        help();
    }


    kseq_t *seq;
    FILE *fp;

    fp = fopen(in_file.c_str(), "r");
    if (fp == Z_NULL) {
        cout << "File: " << in_file << " does not exist" << endl;
        exit(1);
    }
    seq = kseq_init(fileno(fp));
    eof = false;


    uint64_t total = 0;
    int th_list[thread_pool_size];
    uint64_t no_kmers_t[thread_pool_size];
    uint64_t no_kmers = 0;

    bool flag = false;

    vector<SMap> consolidated_vmap(64);
    queue<char *> Q[thread_pool_size];
    thread t_pool[thread_pool_size];
    vector<SMap> all_maps[thread_pool_size];
    mutex mtx[thread_pool_size];
    condition_variable cond[thread_pool_size];


    for (int i = 0; i < thread_pool_size; ++i) {
        no_kmers_t[i] = 0;
        th_list[i] = 0;
        all_maps[i] = vector<SMap>(64);
        t_pool[i] = thread(&consume, ref(all_maps[i]), ref(Q[i]), ref(mtx[i]), ref(cond[i]), k, s / thread_pool_size,
                      ref(th_list[i]),
                      ref(no_kmers_t[i]));
    }


    //producer logic
    int thread_no = 0;  // thread index
    while (true)
    {
        unique_lock<mutex> lck(mtx[thread_no]);

        while (!Q[thread_no].empty())
            cond[thread_no].wait(lck);

        for (int push = 0; push < q_size; push++) {
            int l = kseq_read(seq);
            if (l < 0) {
                flag = true;
                break;
            }
            ++total;
            char *st = strdup(seq->seq.s);
            Q[thread_no].push(st);
        }
        if (!Q[thread_no].empty()) {
            cond[thread_no].notify_one();
        }
        if (flag)
            break;
        thread_no = (thread_no + 1) % thread_pool_size;
    }

    eof = true;

    for (int i = 0; i < thread_pool_size; ++i) {
        cond[i].notify_one();
        t_pool[i].join();
    }


    //consolidation of threads data
    int co = 0;
    for (auto &item : all_maps) {
        no_kmers += no_kmers_t[co];
        for (int i = th_list[co]; i < 64; i++) {
            SMap::iterator it = item[i].begin();
            while (it != item[i].end()) {
                if (consolidated_vmap[i].find(it->first) == consolidated_vmap[i].end()) {
                    consolidated_vmap[i][it->first] = it->second;
                } else {
                    consolidated_vmap[i][it->first] = consolidated_vmap[i][it->first] + it->second;
                }
                ++it;
            }
        }
        co++;
    }

    int *th = max_element(th_list, th_list + thread_pool_size);
    cout << "th: " << *th << endl;
    cout << "No. of sequences: " << total << endl;
    FILE *fo = fopen(out_file.c_str(), "w");
    uint32_t csize = 0;
    for (int i = *th; i < 64; i++) csize += consolidated_vmap[i].size();
    unsigned long F0 = csize * pow(2, (*th));
    cout << "F0: " << F0 << endl;
    fprintf(fo, "F1\t%llu\n", (uint64_t) no_kmers);
    fprintf(fo, "F0\t%lu\n", F0);
    cout << endl;
    cout << "total: " << total << endl;
    cout << "no_kmer: " << no_kmers << endl;
    unsigned long *freq = new unsigned long[c];
    for (int i = 1; i <= c; i++) freq[i] = 0;
    unsigned long tot = 0;
    int xx = 0;
    for (int i = *th; i < 64; i++) {
        for (auto &p: consolidated_vmap[i]) {
            if (p.second <= c) freq[p.second]++;
        }
    }
    cout << "th: " << *th << endl;
    for (int i = 1; i <= c; i++) {
        unsigned long fff = (freq[i] * pow(2, *th));
        fprintf(fo, "f%d\t%lu\n", i, fff);
    }
    fclose(fo);



    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> dur = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    double elapsed = dur.count();
    cout << "\n\n Time taken = " << elapsed << " seconds\n" << endl;

    return 0;
}
