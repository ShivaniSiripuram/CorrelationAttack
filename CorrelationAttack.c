#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

typedef unsigned long long u64;

static const double p=0.25;// prob[sigma(t)!=st(t)]

double phiInvPn=0.8416212335729143;  //(Φ^(−1)(Pn))
double phiInvOneMinusPf=8.209536151601387; //Φ^(−1)(1 − Pf) 


double calT(double N){     //T = (1 − 2p) + 2(Φ^(−1)(Pn))*√ p(1 − p)N .
    double sp = sqrt((p * (1.0 - p))/N); 
    return (1.0 - (2.0*p)) + (2.0 * (phiInvPn) * (sp));
}


int calN(int l) {             //N=Φ−1^(-1)(1 − Pf) − 2Φ−1(Pn)√p(1 − p)1 − 2p
    double x=(phiInvOneMinusPf-2*(phiInvPn*sqrt(p*(1-p))))/1-2*p;
    return (x*x);
}
// double calT(int l,int N){
//     return (sqrt(l)/sqrt(N));
// }

double calC(int * geffeSeq,int* seq,int N){
    double c=0.0;
    for(int i=0;i<N;i++){
        if(seq[i]!=geffeSeq[i])
        c+=1;
    }
    return 1.0 - ((2.0 / N)*c);
}




/* ---------- primitive polynomials ---------- */

static const uint64_t poly_deg1[] = { 0X3 };
static const uint64_t poly_deg2[] = { 0X7 };
static const uint64_t poly_deg3[] = { 0XB, 0XD };
static const uint64_t poly_deg4[] = { 0X13, 0X19 };
static const uint64_t poly_deg5[] = { 0X25, 0X29, 0X2F, 0X37, 0X3B, 0X3D };
static const uint64_t poly_deg6[] = { 0X43, 0X5B, 0X61, 0X67, 0X6D, 0X73 };
static const uint64_t poly_deg7[] = { 0X83, 0X89, 0X8F, 0X91, 0X9D, 0XA7, 0XAB, 0XB9, 0XBF, 0XC1, 0XCB, 0XD3, 0XD5, 0XE5, 0XEF, 0XF1, 0XF7, 0XFD };
static const uint64_t poly_deg8[] = { 0X11D, 0X12B, 0X12D, 0X14D, 0X15F, 0X163, 0X165, 0X169, 0X171, 0X187, 0X18D, 0X1A9, 0X1C3, 0X1CF, 0X1E7, 0X1F5 };
static const uint64_t poly_deg9[] = { 0X211, 0X21B, 0X221, 0X22D, 0X233, 0X259, 0X25F, 0X269, 0X26F, 0X277, 0X27D, 0X287, 0X295, 0X2A3, 0X2A5, 0X2AF, 0X2B7, 0X2BD, 0X2CF, 0X2D1, 0X2DB, 0X2F5, 0X2F9, 0X313, 0X315, 0X31F, 0X323, 0X331, 0X33B, 0X34F, 0X35B, 0X361, 0X36B, 0X36D, 0X373, 0X37F, 0X385, 0X38F, 0X3B5, 0X3B9, 0X3C7, 0X3CB, 0X3CD, 0X3D5, 0X3D9, 0X3E3, 0X3E9, 0X3FB };
static const uint64_t poly_deg10[] = { 0X409, 0X41B, 0X427, 0X42D, 0X465, 0X46F, 0X481, 0X48B, 0X4C5, 0X4D7, 0X4E7, 0X4F3, 0X4FF, 0X50D, 0X519, 0X523, 0X531, 0X53D, 0X543, 0X557, 0X56B, 0X585, 0X58F, 0X597, 0X5A1, 0X5C7, 0X5E5, 0X5F7, 0X5FB, 0X613, 0X615, 0X625, 0X637, 0X643, 0X64F, 0X65B, 0X679, 0X67F, 0X689, 0X6B5, 0X6C1, 0X6D3, 0X6DF, 0X6FD, 0X717, 0X71D, 0X721, 0X739, 0X747, 0X74D, 0X755, 0X759, 0X763, 0X77D, 0X78D, 0X793, 0X7B1, 0X7DB, 0X7F3, 0X7F9 };
static const uint64_t poly_deg11[] = { 0X805, 0X817, 0X82B, 0X82D, 0X847, 0X863, 0X865, 0X871, 0X87B, 0X88D, 0X895, 0X89F, 0X8A9, 0X8B1, 0X8CF, 0X8D1, 0X8E1, 0X8E7, 0X8EB, 0X8F5, 0X90D, 0X913, 0X925, 0X929, 0X93B, 0X93D, 0X945, 0X949, 0X951, 0X95B, 0X973, 0X975, 0X97F, 0X983, 0X98F, 0X9AB, 0X9AD, 0X9B9, 0X9C7, 0X9D9, 0X9E5, 0X9F7, 0XA01, 0XA07, 0XA13, 0XA15, 0XA29, 0XA49, 0XA61, 0XA6D, 0XA79, 0XA7F, 0XA85, 0XA91, 0XA9D, 0XAA7, 0XAAB, 0XAB3, 0XAB5, 0XAD5, 0XADF, 0XAE9, 0XAEF, 0XAF1, 0XAFB, 0XB03, 0XB09, 0XB11, 0XB33, 0XB3F, 0XB41, 0XB4B, 0XB59, 0XB5F, 0XB65, 0XB6F, 0XB7D, 0XB87, 0XB8B, 0XB93, 0XB95, 0XBAF, 0XBB7, 0XBBD, 0XBC9, 0XBDB, 0XBDD, 0XBE7, 0XBED, 0XC0B, 0XC0D, 0XC19, 0XC1F, 0XC57, 0XC61, 0XC6B, 0XC73, 0XC85, 0XC89, 0XC97, 0XC9B, 0XC9D, 0XCB3, 0XCBF, 0XCC7, 0XCCD, 0XCD3, 0XCD5, 0XCE3, 0XCE9, 0XCF7, 0XD03, 0XD0F, 0XD1D, 0XD27, 0XD2D, 0XD41, 0XD47, 0XD55, 0XD59, 0XD63, 0XD6F, 0XD71, 0XD93, 0XD9F, 0XDA9, 0XDBB, 0XDBD, 0XDC9, 0XDD7, 0XDDB, 0XDE1, 0XDE7, 0XDF5, 0XE05, 0XE1D, 0XE21, 0XE27, 0XE2B, 0XE33, 0XE39, 0XE47, 0XE4B, 0XE55, 0XE5F, 0XE71, 0XE7B, 0XE7D, 0XE81, 0XE93, 0XE9F, 0XEA3, 0XEBB, 0XECF, 0XEDD, 0XEF3, 0XEF9, 0XF0B, 0XF19, 0XF31, 0XF37, 0XF5D, 0XF6B, 0XF6D, 0XF75, 0XF83, 0XF91, 0XF97, 0XF9B, 0XFA7, 0XFAD, 0XFB5, 0XFCD, 0XFD3, 0XFE5, 0XFE9 };
static const uint64_t poly_deg12[] = { 0X1053, 0X1069, 0X107B, 0X107D, 0X1099, 0X10D1, 0X10EB, 0X1107, 0X111F, 0X1123, 0X113B, 0X114F, 0X1157, 0X1161, 0X116B, 0X1185, 0X11B3, 0X11D9, 0X11DF, 0X120D, 0X1237, 0X123D, 0X1267, 0X1273, 0X127F, 0X12B9, 0X12C1, 0X12CB, 0X130F, 0X131D, 0X1321, 0X1339, 0X133F, 0X134D, 0X1371, 0X1399, 0X13A3, 0X13A9, 0X1407, 0X1431, 0X1437, 0X144F, 0X145D, 0X1467, 0X1475, 0X14A7, 0X14AD, 0X14D3, 0X150F, 0X151D, 0X154D, 0X1593, 0X15C5, 0X15D7, 0X15DD, 0X15EB, 0X1609, 0X1647, 0X1655, 0X1659, 0X16A5, 0X16BD, 0X1715, 0X1719, 0X1743, 0X1745, 0X1775, 0X1789, 0X17AD, 0X17B3, 0X17BF, 0X17C1, 0X1857, 0X185D, 0X1891, 0X1897, 0X18B9, 0X18EF, 0X191B, 0X1935, 0X1941, 0X1965, 0X197B, 0X198B, 0X19B1, 0X19BD, 0X19C9, 0X19CF, 0X19E7, 0X1A1B, 0X1A2B, 0X1A33, 0X1A69, 0X1A8B, 0X1AD1, 0X1AE1, 0X1AF5, 0X1B0B, 0X1B13, 0X1B1F, 0X1B57, 0X1B91, 0X1BA7, 0X1BBF, 0X1BC1, 0X1BD3, 0X1C05, 0X1C11, 0X1C17, 0X1C27, 0X1C4D, 0X1C87, 0X1C9F, 0X1CA5, 0X1CBB, 0X1CC5, 0X1CC9, 0X1CCF, 0X1CF3, 0X1D07, 0X1D23, 0X1D43, 0X1D51, 0X1D5B, 0X1D75, 0X1D85, 0X1D89, 0X1E15, 0X1E19, 0X1E2F, 0X1E45, 0X1E51, 0X1E67, 0X1E73, 0X1E8F, 0X1EE3, 0X1F11, 0X1F1B, 0X1F27, 0X1F71, 0X1F99, 0X1FBB, 0X1FBD, 0X1FC9 };

/* Pointer table */
static const uint64_t *primitive_poly[] = {
    poly_deg1, poly_deg2, poly_deg3, poly_deg4, poly_deg5,
    poly_deg6, poly_deg7, poly_deg8, poly_deg9, poly_deg10,
    poly_deg11, poly_deg12
};

/* Count table */
static const int primitive_poly_sizes[] = {
    1,1,2,2,6,6,18,16,48,60,176,144
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

 
    int initial_state[30];
    for (int i = 0; i < len; ++i)
        initial_state[i] = 1;

    int seq[5000];    
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

            double C = calC(seq2, Geffeseq, N);
            //printf("C value : %lf\n",C);
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

    int AttackingLFSR_Length=6;
    int N=calN(AttackingLFSR_Length);//calN(p,AttackingLFSR_Length);
    double T=calT(N);//calT(AttackingLFSR_Length,N);
    printf("N : %d", N);
    printf(" T : %f\n",T);
    //lfsr1
    int len1=6;
    int initial_state1[]={0,0,1,0,0,1};//[s2,s1,s0]
    int no_of_taps1=2;
    int feedback_poly1[]={1,2};
    int period1 = (1 << len1) - 1;
    int seq1[N];
    generateSequence(initial_state1, feedback_poly1, seq1, len1, no_of_taps1,N);
    //lfsr1 seq
    printf("LFSR-1 sequence is:\n");
    for(int i = 0; i < N; i++)
        printf("%d ", seq1[i]);
    printf("\n");

    //lfsr2
    int len2=10;
    int initial_state2[]={0,1,0,1,0,1,1,1,0,1};//[s3,s2,s1,s0]
    int no_of_taps2=2;
    int feedback_poly2[]={1,4};
    int period2 = (1 << len2) - 1;
    int seq2[N];
    generateSequence(initial_state2, feedback_poly2, seq2, len2, no_of_taps2,N);
    //lfsr2 seq
    printf("LFSR-2 sequence is:\n");
    for(int i = 0; i < N; i++)
        printf("%d ", seq2[i]);
    printf("\n");

    //lfsr-3
    int len3=11;
    int initial_state3[]={1,1,1,0,1,0,0,0,1,0,1};//[s4,s3,s2,s1,s0]
    int no_of_taps3=2;
    int feedback_poly3[]={1,3};
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



