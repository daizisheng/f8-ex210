#define precision 18
#define progress_mask 0x1fffff
#define maxn 2000000 \

#define radix 1000000000000000000LL
#define MODP 1000000007LL
#define STABLE_LOW 3
#define STABLE_HIGH 60
#define maxprec ((precision+17) /18)  \

#define memtyp unsigned int
#define memsize 1000000000LL
#define wtsize 10000000LL
#define oldmemsize 450000000
#define deg 9 \

#define bitsperword 8*sizeof(unsigned long long)  \

#define vert(k) (g->vertices+(k) -1) 
#define vertnum(v) ((v) -g->vertices+1)  \

#define encode(x) ((x) <0?'#':(x) <10?(x) +'0':(x) -10+'a')  \

/*3:*/
#line 112 "dynaham.w"

#include <stdio.h> 
#include <stdlib.h> 
#include "gb_graph.h"
#include "gb_save.h"
int m;
int n;
/*5:*/
#line 171 "dynaham.w"

typedef struct bignum_struct{
long long val[maxprec];
}bignum;

/*:5*/
#line 119 "dynaham.w"
;
/*6:*/
#line 176 "dynaham.w"

bignum zero;
bignum one;
bignum infty;
int prec;

/*:6*//*15:*/
#line 379 "dynaham.w"

memtyp oldp;
long long contribs;
memtyp*mem;
unsigned long long*oldmem;
memtyp memptr;
memtyp wtptr;
int q,oldq;
int code[maxn];
int oldcode[maxn];
bignum*weight,*oldweight;
bignum minweight,maxweight;
bignum count[maxn+1];
int maxdeg;
unsigned long long pack;
int spack;
int owp,omp;
long long maxmemptr,maxwtptr,maxomp;
int tmap[deg],itmap[deg],omap[deg],iomap[deg];
int tms,tmx,oms,omx;

/*:15*//*27:*/
#line 624 "dynaham.w"

int fr[maxn],ifr[maxn+1];
int ofr[maxn];
int q0;

/*:27*//*35:*/
#line 732 "dynaham.w"

int path[maxn];
int mate[maxn],oldmate[maxn];
int bmate[maxn];
int mp[maxn+1];
int imap[maxn];
int r;
int nbr[maxn];
int steps;

/*:35*/
#line 120 "dynaham.w"
;
/*8:*/
#line 187 "dynaham.w"

void add_to_bignum(bignum*x,bignum delta){
  x->val[0] = (x->val[0] + delta.val[0]) % MODP;
}


/*:8*//*9:*/
#line 203 "dynaham.w"

int bignum_comp(bignum x,bignum y){
  if (x.val[0] < y.val[0]) return -1;
  if (x.val[0] > y.val[0]) return 1;
  return 0;
}


/*:9*//*10:*/
#line 211 "dynaham.w"

void print_bignum(FILE*stream,bignum x){
  fprintf(stream, "%lld", x.val[0]);
}


/*:10*//*16:*/
#line 409 "dynaham.w"

memtyp trielookup(void){
register int j,l,k,kk;
register memtyp p,pp;
tms= 0,tmx= p= 1;
for(l= 0;l<q;l++,p= mem[pp]){
j= code[l];
if(j<=0)pp= p+j;
else if(j==tmx){
/*17:*/
#line 430 "dynaham.w"

if(tmx> maxdeg){
maxdeg= tmx;
if(tmx==deg){
fprintf(stderr,"Overflow: code digits must be less than %d!\n",
deg);
exit(-66);
}
}
tmap[tms]= tmx,itmap[tmx]= tms;
tms++,tmx++;
pp= p+tms;

/*:17*/
#line 418 "dynaham.w"
;
}else{
/*18:*/
#line 443 "dynaham.w"

k= tmap[--tms],kk= itmap[j];
tmap[kk]= k,itmap[k]= kk;
pp= p+1+kk;

/*:18*/
#line 420 "dynaham.w"
;
}
if(mem[pp]==0){
if(l+1<q)/*19:*/
#line 451 "dynaham.w"

{
register int slots;
if(l+1+tms<q){
mem[memptr]= mem[memptr+1]= 0,memptr+= 2;
slots= tms+(l+2+tms==q?0:1);
}else slots= tms;
if(memptr+slots+1>=memsize){
fprintf(stderr,"Oops: Dictionary overflow (more than %lld pointers)!\n",
memsize);
exit(-666);
}
mem[pp]= memptr-1;
for(j= 0;j<slots;j++)mem[memptr+j]= 0;
memptr+= slots;
}

/*:19*/
#line 423 "dynaham.w"

else/*20:*/
#line 468 "dynaham.w"

{
mem[pp]= ++wtptr;
if(wtptr>=wtsize){
fprintf(stderr,"Oops: Dictionary overflow (more than %lld classes)!\n",
wtsize);
exit(-6666);
}
weight[wtptr-1]= zero;
}

/*:20*/
#line 424 "dynaham.w"
;
}
}
return p-1;
}

/*:16*//*39:*/
#line 817 "dynaham.w"

void contribute(void){
register int j,k,t;
register memtyp p;
/*33:*/
#line 711 "dynaham.w"

for(t= 0,k= 1;k<=q;k++){
j= mate[k];
if(j<=0)code[k-1]= j;
else if(j> k)code[k-1]= code[j-1]= ++t;
}

/*:33*/
#line 821 "dynaham.w"
;
p= trielookup();
if (m >= STABLE_LOW && m <= STABLE_HIGH) {
    fprintf(stderr, "T m=%d old=", m);
    for (int _l = 0; _l < oldq; _l++) fputc(encode(oldcode[_l]), stderr);
    fprintf(stderr, " new=");
    for (int _l = 0; _l < q; _l++) fputc(encode(code[_l]), stderr);
    fputc('\n', stderr);
}
add_to_bignum(&weight[p],oldweight[oldp]);
contribs++;
}

/*:39*//*44:*/
#line 896 "dynaham.w"

int add_derived(int i,int j){
register int k,kk;
if(mate[i]<0||mate[j]<0)return 0;
if(!mate[i]){
if(!mate[j]){
if(i==j)goto cycle;
mate[i]= j,mate[j]= i;
return 1;
}else{
mate[i]= mate[j],mate[mate[j]]= i,mate[j]= -1;
return 1;
}
}else if(!mate[j]){
mate[j]= mate[i],mate[mate[i]]= j,mate[i]= -1;
return 1;
}else if(mate[i]!=j){
mate[mate[i]]= mate[j],mate[mate[j]]= mate[i];
mate[i]= mate[j]= -1;
return 1;
}else/*45:*/
#line 935 "dynaham.w"

{
cycle:mate[i]= mate[j]= -1;
for(k= 1;k<=q0;k++){
if(mate[k]>=0)break;
if(fr[k-1]!=m+k)return 0;
}
for(kk= k;kk<=q;kk++)
if(mate[kk])return 0;
{
register int l;
fprintf(stderr,"Class ");
for(l= 0;l<oldq;l++)fprintf(stderr,"%c",
encode(oldcode[l]));
fprintf(stderr," contributes ");
print_bignum(stderr,oldweight[oldp]);
fprintf(stderr," to a %d-cycle.\n",
m+k-1);
if (m >= STABLE_LOW && m <= STABLE_HIGH) {
    fprintf(stderr, "C m=%d cycle_at=%d old=", m, m+k-1);
    for (int _l = 0; _l < oldq; _l++) fputc(encode(oldcode[_l]), stderr);
    fputc('\n', stderr);
}
add_to_bignum(&count[m+k-1],oldweight[oldp]);
}
return 0;
}

/*:45*/
#line 916 "dynaham.w"
;
}

/*:44*//*53:*/
#line 1053 "dynaham.w"

void report_cycles(int m){
if(bignum_comp(count[m],one)>=0){
printf("There are ");
print_bignum(stdout,count[m]);
printf(" Hamiltonian %d-cycles.\n",
m);
fflush(stdout);
}
}

/*:53*/
#line 121 "dynaham.w"
;
/*21:*/
#line 490 "dynaham.w"

void compress(int l,memtyp p,int d){

register int j,k,kk,kkk;
register unsigned long long bits;
if(l==q)/*23:*/
#line 527 "dynaham.w"

{
for(k= 0;k<prec;k++)
oldweight[owp].val[k]= weight[p-1].val[k];
owp++;
}

/*:23*/
#line 495 "dynaham.w"

else{
for(j= (l+oms==q?1:-1),k= 0,bits= 0;k<d;j++,k++)
if(mem[p+j])bits+= 1LL<<k;
/*22:*/
#line 517 "dynaham.w"

if(spack+d> bitsperword){
oldmem[omp++]= pack,pack= bits,spack= d;
if(omp>=oldmemsize){
fprintf(stderr,"Oops: oldmem overflow (more than %d bytes)!\n",
oldmemsize);
exit(-666666);
}
}else pack+= bits<<spack,spack+= d;

/*:22*/
#line 499 "dynaham.w"
;
for(j= (l+oms==q?1:-1),k= 0;k<d;j++,k++)if(bits&(1LL<<k)){
if(j> 0){
if(j> oms)omap[oms]= omx,iomap[omx]= oms,oms++,omx++,kk= 0;
else kk= omap[--oms],kkk= omap[j-1],omap[j-1]= kk,iomap[kk]= j-1;
}
compress(l+1,mem[p+j],oms+(l+1+oms==q?0:
l+2+oms==q?2:3));
if(j> 0){
if(!kk)oms--,omx--;
else omap[j-1]= kkk,omap[oms]= kk,iomap[kk]= oms++;
}
}
}
}

/*:21*//*25:*/
#line 558 "dynaham.w"

void uncompress(int l,int d){
register int i,j,k,ii,kk,kkk;
register unsigned long long bits;
if(l==oldq){
if (m >= STABLE_LOW && m <= STABLE_HIGH + 1) {
    fprintf(stderr, "W m=%d oldp=%u code=", m-1, oldp);
    for (int _l = 0; _l < oldq; _l++) fputc(encode(oldcode[_l]), stderr);
    fprintf(stderr, " w=%lld\n", oldweight[oldp].val[0]);
}
/*47:*/
#line 970 "dynaham.w"

{
/*34:*/
#line 720 "dynaham.w"

for(k= 1;k+k<=oldq;k++)path[k]= 0;
for(k= 1;k<=oldq;k++){
j= oldcode[k-1];
if(j<=0)oldmate[k]= j;
else{
l= path[j];
if(!l)path[j]= k;
else oldmate[k]= l,oldmate[l]= k;
}
}

/*:34*/
#line 972 "dynaham.w"
;
/*48:*/
#line 980 "dynaham.w"

if(oldp==0)/*49:*/
#line 986 "dynaham.w"

{
fprintf(stderr,"The first one is ");
for(l= 0;l<oldq;l++)fprintf(stderr,"%c",
encode(oldcode[l]));
fprintf(stderr,"\nand its weight is ");
print_bignum(stderr,oldweight[oldp]);
fprintf(stderr,".\n");
minweight= maxweight= oldweight[oldp];
steps= contribs= 0;
fflush(stderr);
}

/*:49*/
#line 981 "dynaham.w"

else/*50:*/
#line 999 "dynaham.w"

{
if(bignum_comp(oldweight[oldp],minweight)<0){
minweight= oldweight[oldp];
fprintf(stderr,"Class ");
for(l= 0;l<oldq;l++)fprintf(stderr,"%c",
encode(oldcode[l]));
fprintf(stderr," has weight ");
print_bignum(stderr,oldweight[oldp]);
fprintf(stderr,".\n");
}
if(bignum_comp(oldweight[oldp],maxweight)> 0){
maxweight= oldweight[oldp];
fprintf(stderr,"Class ");
for(l= 0;l<oldq;l++)fprintf(stderr,"%c",
encode(oldcode[l]));
fprintf(stderr," has weight ");
print_bignum(stderr,oldweight[oldp]);
fprintf(stderr,".\n");
}
if((steps++&progress_mask)==progress_mask){
fprintf(stderr,".");
fflush(stderr);
}
}

/*:50*/
#line 982 "dynaham.w"
;

/*:48*/
#line 973 "dynaham.w"
;
/*37:*/
#line 804 "dynaham.w"

for(j= 1;j<=q0;j++)bmate[j]= mp[1+oldmate[imap[j]]];
for(;j<=q;j++)bmate[j]= 0;

/*:37*/
#line 974 "dynaham.w"
;
if(oldmate[1]<0)/*38:*/
#line 808 "dynaham.w"

{
for(j= 1;j<=q;j++)mate[j]= bmate[j];
contribute();
}

/*:38*/
#line 975 "dynaham.w"

else if(oldmate[1]==0)/*41:*/
#line 848 "dynaham.w"

{
for(i= 0;i<r;i++)for(ii= i+1;ii<r;ii++){
for(j= 1;j<=q;j++)mate[j]= bmate[j];
if(add_derived(nbr[i],nbr[ii]))contribute();
}
}

/*:41*/
#line 976 "dynaham.w"

else/*42:*/
#line 865 "dynaham.w"

{
for(i= 0;i<r;i++){
for(j= 1;j<=q;j++)mate[j]= bmate[j];
if(add_derived(mp[1+oldmate[1]],nbr[i]))contribute();
}
}

/*:42*/
#line 977 "dynaham.w"
;
}

/*:47*/
#line 563 "dynaham.w"

oldp++;
}else{
/*26:*/
#line 583 "dynaham.w"

if(spack+d> bitsperword)
pack= oldmem[++omp],spack= 0;
bits= pack>>spack,spack+= d;

/*:26*/
#line 566 "dynaham.w"
;
for(j= (l+oms==oldq?1:-1),k= 0;k<d;j++,k++)if(bits&(1LL<<k)){
if(j<=0)oldcode[l]= j;
else if(j> oms)
oldcode[l]= omx,omap[oms]= omx,iomap[omx]= oms,omx++,oms++,kk= 0;
else oldcode[l]= omap[j-1],kk= omap[--oms],
kkk= omap[j-1],omap[j-1]= kk,iomap[kk]= j-1;
uncompress(l+1,oms+(l+1+oms==oldq?0:
l+2+oms==oldq?2:3));
if(j> 0){
if(!kk)oms--,omx--;
else omap[j-1]= kkk,omap[oms]= kk,iomap[kk]= oms++;
}
}
}
}

/*:25*/
#line 122 "dynaham.w"
;
int main(int argc,char*argv[]){
register int c,d,i,j,k,l,p,t,x,ii,ik,iv;
register Graph*g;
register Arc*a;
/*4:*/
#line 139 "dynaham.w"

if(argc!=2){
fprintf(stderr,"Usage: %s foo.gb\n",argv[0]);
exit(-1);
}
g= restore_graph(argv[1]);
if(!g){
fprintf(stderr,"I couldn't reconstruct graph %s!\n",argv[1]);
exit(-2);
}
n= g->n;
if(n> maxn){
fprintf(stderr,"Recompile me: I allow at most %d vertices!\n",
maxn);
exit(-3);
}
printf("Dynamic Hamiltonian cycles of the graph %s",
g->id);
printf(" (%d vertices, %ld edges):\n",
n,g->m/2);
fflush(stdout);

/*:4*/
#line 127 "dynaham.w"
;
/*7:*/
#line 182 "dynaham.w"

one.val[0]= 1;
for(k= 0;k<maxprec;k++)infty.val[k]= radix-1;
prec= 1;

/*:7*//*14:*/
#line 369 "dynaham.w"

mem= (memtyp*)malloc(memsize*sizeof(memtyp));
oldmem= (unsigned long long*)malloc(oldmemsize);
weight= (bignum*)malloc(wtsize*sizeof(bignum));
oldweight= (bignum*)malloc(wtsize*sizeof(bignum));
if(!mem||!oldmem||!weight||!oldweight){
fprintf(stderr,"I can't allocate the big tables!\n");
exit(-6);
}

/*:14*//*28:*/
#line 629 "dynaham.w"

for(k= 0;k<n;k++)fr[k]= k+1,ifr[k+1]= k;
q= 1;
mem[0]= mem[2]= 0,mem[1]= 1,memptr= 3;
weight[0]= one,wtptr= 1;

/*:28*//*36:*/
#line 794 "dynaham.w"

mp[0]= -1;

/*:36*/
#line 128 "dynaham.w"
;
for(m= 1;;m++){
/*32:*/
#line 688 "dynaham.w"

fprintf(stderr,"\nThe frontier for %d-classes is",
m-1);
for(k= 0;k<q;k++)
fprintf(stderr," %s",
vert(fr[k])->name);
fprintf(stderr,".\n");

/*:32*/
#line 130 "dynaham.w"
;
if(wtptr==0)/*52:*/
#line 1037 "dynaham.w"

{
fprintf(stderr,"\n");
for(k= m;k<=n;k++){
if(k> m+m)break;
report_cycles(k);
}
fprintf(stderr,"\nThat's all; there are no %d-configs!\n",
m-1);
fprintf(stderr,"Storage requirements: %lld memsize, %lld omemsize,",
(long long)maxmemptr+2,(long long)maxomp+1);
fprintf(stderr," %lld wtsize, %d maxprec, %d deg.\n",
(long long)maxwtptr,18*prec,maxdeg+1);
exit(0);
}

/*:52*/
#line 131 "dynaham.w"
;
/*24:*/
#line 534 "dynaham.w"

fprintf(stderr,"There are %lld %d-classes,\n",
(long long)wtptr,m-1);
fprintf(stderr," resulting from %lld contributions and filling %lld pointers.\n",
contribs,(long long)memptr);
spack= omp= owp= oms= 0,omx= 1;
pack= 0;
oldq= q;
compress(0,1,3);
oldmem[omp]= pack;
if(omp> maxomp)maxomp= omp;

/*:24*/
#line 132 "dynaham.w"
;
/*51:*/
#line 1025 "dynaham.w"

if(memptr> maxmemptr)maxmemptr= memptr;
if(wtptr> maxwtptr)maxwtptr= wtptr;
mem[0]= mem[1]= mem[2]= 0,memptr= 3,wtptr= 0;

/*29:*/
#line 635 "dynaham.w"

for(j= 1;j<q;j++)ofr[j]= fr[j];
oldq= q;
iv= ifr[m+1];
if(iv<--q){
x= fr[q];
fr[0]= m+1,ifr[m+1]= 0;
fr[q]= m,ifr[m]= q;
fr[iv]= x,ifr[x]= iv;
}else{
fr[0]= m+1,ifr[m+1]= 0,fr[iv]= m,ifr[m]= iv;
if(iv!=q)q++;
}
q0= q;
/*30:*/
#line 656 "dynaham.w"

for(a= vert(m)->arcs;a;a= a->next){
k= vertnum(a->tip);
if(k<m)continue;
ik= ifr[k];
if(ik>=q){
x= fr[q];
fr[q]= k,ifr[k]= q,fr[ik]= x,ifr[x]= ik;
q++;
}
}

/*:30*/
#line 649 "dynaham.w"
;
/*31:*/
#line 668 "dynaham.w"

for(k= 2;k<q0;k++)if(fr[k]<fr[k-1]){
for(t= fr[k],j= k-1;;j--){
fr[j+1]= fr[j],ifr[fr[j]]= j+1;
if(j==0||fr[j-1]<t)break;
}
fr[j]= t,ifr[t]= j;
}
for(k= q0+1;k<q;k++)if(fr[k]<fr[k-1]){
for(t= fr[k],j= k-1;;j--){
fr[j+1]= fr[j],ifr[fr[j]]= j+1;
if(j==q0||fr[j-1]<t)break;
}
fr[j]= t,ifr[t]= j;
}

/*:31*/
#line 650 "dynaham.w"
;
/*40:*/
#line 833 "dynaham.w"

for(r= 0,a= vert(m)->arcs;a;a= a->next){
k= vertnum(a->tip);
if(k<m)continue;
nbr[r++]= ifr[k]+1;
}

/*:40*//*43:*/
#line 877 "dynaham.w"

imap[1]= 0;
for(j= 2;j<=oldq;j++){
mp[j+1]= 1+ifr[ofr[j-1]];
imap[mp[j+1]]= j;
}

/*:43*/
#line 651 "dynaham.w"
;

/*:29*/
#line 1030 "dynaham.w"
;

/*:51*/
#line 133 "dynaham.w"
;
/*46:*/
#line 966 "dynaham.w"

spack= oldp= omp= oms= 0,omx= 1,pack= oldmem[0];
uncompress(0,3);

/*:46*/
#line 134 "dynaham.w"
;
report_cycles(m);
}
}

/*:3*/
