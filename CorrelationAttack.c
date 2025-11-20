#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
typedef unsigned long long u64;

static const double p=0.25;


double calT(double N) {     //T = (1 − 2p) + 2(Φ^(−1)(Pn))*√ p(1 − p)N .
    double sp = sqrt((p * (1.0 - p))/N); 
    return (1.0 - 2.0*p) + 2.0 * (-3.0) * (sp);
}


int calN(int l) {             //N=Φ−1(1 − Pf) − 2Φ−1(Pn)√p(1 − p)1 − 2p
    float b = sqrt(p * (1 - p));
    float a = (sqrt(l) - ((2 * (-3)) * b)) / (1 - (2 * p));
    return (ceil(a * a));
}
// double calT(int l,int N){
//     return (sqrt(l)/sqrt(N));
// }

double calC(int * geffeSeq,int* seq,int N){
    float c=0.0;
    for(int i=0;i<N;i++){
        c+=(seq[i]^geffeSeq[i]);
    }
    return 1.0 - (2.0 * c / N);
}




/* ---------- primitive polynomials ---------- */
static const uint64_t poly_deg1[]  = { 0x3 }; 
static const uint64_t poly_deg2[]  = { 0x7 };
static const uint64_t poly_deg3[]  = { 0xB, 0xD };
static const uint64_t poly_deg4[]  = { 0x13, 0x19 };
static const uint64_t poly_deg5[]  = { 0x25, 0x29, 0x2F };
static const uint64_t poly_deg6[]  = { 0x43, 0x47, 0x4D, 0x53 };
static const uint64_t poly_deg7[]  = { 0x89, 0x8F, 0xA1, 0xA7, 0xB1, 0xB7 };
static const uint64_t poly_deg8[]  = { 0x11D, 0x12B, 0x12D, 0x14D, 0x15F, 0x163 };
static const uint64_t poly_deg9[]  = { 0x211, 0x247, 0x253, 0x257, 0x283, 0x287, 0x2C9 };
static const uint64_t poly_deg10[] = {
    0x409, 0x40F, 0x413, 0x415, 0x443,
   0x445, 0x46D, 0x481, 0x48B, 0x4C3
};
static const uint64_t poly_deg11[] = {
    0x821, 0x829, 0x841, 0x853, 0x857,
    0x859, 0x85F, 0x86B, 0x871, 0x883
};
static const uint64_t poly_deg12[] = {
    0x1053, 0x1069, 0x106F, 0x108B,
    0x10A1, 0x10AF, 0x10D1, 0x10DB, 0x10ED
};

/* Pointer table */
static const uint64_t *primitive_poly[] = {
    poly_deg1, poly_deg2, poly_deg3, poly_deg4, poly_deg5,
    poly_deg6, poly_deg7, poly_deg8, poly_deg9, poly_deg10,
    poly_deg11, poly_deg12
};

/* Count table */
static const int primitive_poly_sizes[] = {
    1,1,2,2,3,4,6,6,7,10,10,9
};



//geffe sequence generator
 void geffe(int gef[], int s1[],int s2[],int s3[],int l,int period1,int period2,int period3){
    for(int i=0;i<l;i++){
        gef[i]=((s1[i%period1]&s2[i%period2])^(s2[i%period2]&s3[i%period3])^s3[i%period3]);
    }
}

//LFSR seq
void generateSequence(int* initial_state, int* tap_positions, int* seq, int len, int no_of_taps,int N) {
    int j = 0;
    unsigned lfsr = 0;

    
    for(int i = 0; i < len; i++) {
        lfsr |= (initial_state[i] << (len - i - 1));
    }

    unsigned initial_lfsr = lfsr;
int period = (1 << len) - 1;

for (int t = 0; t < N; t++) {
    seq[t] = lfsr & 1;

    unsigned feedback = 0;
    for(int i = 0; i < no_of_taps; i++)
        feedback ^= (lfsr >> (tap_positions[i] - 1)) & 1;

    lfsr = (lfsr >> 1) | (feedback << (len - 1));
}  
}



int extract_taps(uint64_t poly, int degree, int* taps)
{
    int count = 0;

    for (int bit = 0; bit < degree; bit++) {
        if (poly & (1ULL << bit)) {
            taps[count++] = bit + 1;   // 1-based index
        }
    }

    return count;
}




void print_poly(uint64_t poly) {
    printf("Standard form: ");

    int highest = 63;
    while (highest >= 0 && ((poly >> highest) & 1) == 0)
        highest--;

    int first = 1;

    for (int i = highest; i >= 0; i--) {
        if ((poly >> i) & 1) {
            if (!first) printf(" + ");

            if (i == 0)       printf("1");
            else if (i == 1)  printf("x");
            else              printf("x^%d", i);

            first = 0;
        }
    }
}


void correlationAttack(int len, int *Geffeseq, float T, int N)
{
    int L = (1 << len) - 1;

    /* initial state for the first generator */
    int initial_state[30];
    for (int i = 0; i < len; ++i)
        initial_state[i] = 1;

    int seq[5000];      /* adjust if needed */
    int taps[30];

    /* pick first primitive polynomial for this length */
    uint64_t poly = primitive_poly[len - 1][0];
    int no_of_taps = extract_taps(poly, len, taps);

    generateSequence(initial_state, taps, seq, len, no_of_taps, N);

    int c = 0;

    /* loop over all cyclic shifts */
    for (int shift = 0; shift < L; shift++) {

        int ini[30];
        for (int k = 0; k < len; k++)
            ini[k] = seq[(shift + k) % L];

        /* loop over all primitive polynomials for this degree */
        const uint64_t *polylist = primitive_poly[len - 1];
        int poly_count = primitive_poly_sizes[len - 1];

        for (int pi = 0; pi < poly_count; pi++) {

            uint64_t polyy = polylist[pi];

            int taps2[30];
            int seq2[5000];

            int no_of_taps2 = extract_taps(polyy, len, taps2);

            generateSequence(ini, taps2, seq2, len, no_of_taps2, N);

            c++;

            float C = calC(seq2, Geffeseq, N);

            if (C > T) {
                print_poly(polyy);
                printf("\ninitial state is:\n");
                for (int t = 0; t < len; ++t)
                    printf("%d ", ini[t]);
                printf("\n");
            }
        }
    }
}





int main() {

    int AttackingLFSR_Length=3;
    int N=calN(AttackingLFSR_Length);//calN(p,AttackingLFSR_Length);
    double T=calT(N);//calT(AttackingLFSR_Length,N);
    printf("N : %d", N);
    printf(" T : %f\n",T);
    //lfsr1
    int len1=3;
    int initial_state1[]={1,0,1};//[s2,s1,s0]
    int no_of_taps1=2;
    int feedback_poly1[]={1,2};//x^3 + x + 1
    int period1 = (1 << len1) - 1;
    int seq1[N];
    generateSequence(initial_state1, feedback_poly1, seq1, len1, no_of_taps1,N);
    //lfsr1 seq
    printf("LFSR-1 sequence is:\n");
    for(int i = 0; i < N; i++)
        printf("%d ", seq1[i]);
    printf("\n");

    //lfsr2
    int len2=4;
    int initial_state2[]={1,1,0,1};//[s3,s2,s1,s0]
    int no_of_taps2=2;
    int feedback_poly2[]={1,2};// x^4 + x + 1
    int period2 = (1 << len2) - 1;
    int seq2[N];
    generateSequence(initial_state2, feedback_poly2, seq2, len2, no_of_taps2,N);
    //lfsr2 seq
    printf("LFSR-2 sequence is:\n");
    for(int i = 0; i < N; i++)
        printf("%d ", seq2[i]);
    printf("\n");

    //lfsr-3
    int len3=5;
    int initial_state3[]={1,0,0,1,1};//[s4,s3,s2,s1,s0]
    int no_of_taps3=2;
    int feedback_poly3[]={1,3};// x^5 + x^2 + 1
    int period3 = (1 << len3) - 1;
    int seq3[N];
    generateSequence(initial_state3, feedback_poly3, seq3, len3, no_of_taps3,N);
    //lfsr3 seq
    printf("LFSR-3 sequence is:\n");
    for(int i = 0; i < N; i++)
        printf("%d ", seq3[i]);
    printf("\n");

    int gef[N];
    geffe(gef,seq1,seq2,seq3,N,period1,period2,period3);
    printf("Geffe seq : ");
    for(int i=0;i<N;i++){
        printf("%d ",gef[i]);
    }
    
    printf("\n");
    correlationAttack(AttackingLFSR_Length,gef,T,N);
    



    return 0;
}
