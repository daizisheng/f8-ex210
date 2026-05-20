#define precision 18
#define progress_mask 0x1fffff
#define maxn 2000000 \

#define radix 1000000000000000000LL
#define MODP 1000000007LL
#define STABLE_LOW 30
#define STABLE_HIGH 41
#define maxprec ((precision+17) /18)  \

#define memsize 1000000000LL
#define wtsize 10000000LL
#define oldmemsize 450000000
#define deg 10 \

#define bitsperword 8*sizeof(unsigned long long)  \

#define vert(k) (g->vertices+(k) -1) 
#define vertnum(v) ((v) -g->vertices+1)  \

#define encode(x) ((x) <0?'#':(x) <10?(x) +'0':(x) -10+'a')  \

/*2:*/
#line 89 "dynahamp.w"

#include <stdio.h> 
#include <stdlib.h> 
#include "gb_graph.h"
#include "gb_save.h"
int m;
int n;
/*4:*/
#line 149 "dynahamp.w"

typedef struct bignum_struct{
long long val[maxprec];
}bignum;

/*:4*/
#line 96 "dynahamp.w"
;
/*5:*/
#line 154 "dynahamp.w"

bignum zero;
bignum one;
bignum infty;
int prec;

/*:5*//*14:*/
#line 356 "dynahamp.w"

long long oldp;
long long contribs;
unsigned long long*mem;
unsigned long long*oldmem;
long long memptr,maxmemptr;
long long wtptr,maxwtptr;
int q,oldq;
int code[maxn];
int oldcode[maxn];
bignum*weight,*oldweight;
bignum minweight,maxweight;
bignum count[maxn+1];
int maxdeg;
unsigned long long pack;
int spack;
long long owp,omp;
int maxomp;
int tmap[deg],itmap[deg],omap[deg],iomap[deg];
int tms,tmx,oms,omx;

/*:14*//*26:*/
#line 604 "dynahamp.w"

int fr[maxn],ifr[maxn+1];
int ofr[maxn];
int q0;

/*:26*//*34:*/
#line 713 "dynahamp.w"

int path[maxn];
int mate[maxn],oldmate[maxn];
int bmate[maxn];
int mp[maxn+1];
int imap[maxn];
int r;
int nbr[maxn];
int steps;

/*:34*/
#line 97 "dynahamp.w"
;
/*7:*/
#line 165 "dynahamp.w"

void add_to_bignum(bignum*x,bignum delta){
  x->val[0] = (x->val[0] + delta.val[0]) % MODP;
}


/*:7*//*8:*/
#line 181 "dynahamp.w"

int bignum_comp(bignum x,bignum y){
  if (x.val[0] < y.val[0]) return -1;
  if (x.val[0] > y.val[0]) return 1;
  return 0;
}


/*:8*//*9:*/
#line 189 "dynahamp.w"

void print_bignum(FILE*stream,bignum x){
  fprintf(stream, "%lld", x.val[0]);
}


/*:9*//*15:*/
#line 386 "dynahamp.w"

long long trielookup(void){
register int j,l,k,kk;
register unsigned long long p,pp;
tms= 0,tmx= p= 1;
for(l= 0;l<q;l++,p= mem[pp]){
j= code[l];
if(j<=0)pp= p+j;
else if(j==tmx){
/*16:*/
#line 407 "dynahamp.w"

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

/*:16*/
#line 395 "dynahamp.w"
;
}else{
/*17:*/
#line 420 "dynahamp.w"

k= tmap[--tms],kk= itmap[j];
tmap[kk]= k,itmap[k]= kk;
pp= p+1+kk;

/*:17*/
#line 397 "dynahamp.w"
;
}
if(mem[pp]==0){
if(l+1<q)/*18:*/
#line 428 "dynahamp.w"

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

/*:18*/
#line 400 "dynahamp.w"

else/*19:*/
#line 445 "dynahamp.w"

{
mem[pp]= ++wtptr;
if(wtptr>=wtsize){
fprintf(stderr,"Oops: Dictionary overflow (more than %lld classes)!\n",
wtsize);
exit(-6666);
}
weight[wtptr-1]= zero;
}

/*:19*/
#line 401 "dynahamp.w"
;
}
}
return p-1;
}

/*:15*//*38:*/
#line 801 "dynahamp.w"

void contribute(void){
register int j,k,t;
register long long p;
/*32:*/
#line 692 "dynahamp.w"

for(t= 0,k= 1;k<=q;k++){
j= mate[k];
if(j<=0)code[k-1]= j;
else if(j> k)code[k-1]= code[j-1]= ++t;
}

/*:32*/
#line 805 "dynahamp.w"
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

/*:38*//*44:*/
#line 885 "dynahamp.w"

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
#line 924 "dynahamp.w"

{
cycle:mate[i]= mate[j]= -1;
if(mate[1]>=0)return 0;
for(k= 2;k<=q0;k++){
if(mate[k]>=0)break;
if(fr[k-1]!=m+k-1)return 0;
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
fprintf(stderr," to a %d-path.\n",
m+k-2);
if (m >= STABLE_LOW && m <= STABLE_HIGH) {
    fprintf(stderr, "C m=%d cycle_at=%d old=", m, m+k-2);
    for (int _l = 0; _l < oldq; _l++) fputc(encode(oldcode[_l]), stderr);
    fputc('\n', stderr);
}
add_to_bignum(&count[m+k-2],oldweight[oldp]);
}
return 0;
}

/*:45*/
#line 905 "dynahamp.w"
;
}

/*:44*//*53:*/
#line 1042 "dynahamp.w"

void report_cycles(int m){
if(bignum_comp(count[m],one)>=0){
printf("There are ");
print_bignum(stdout,count[m]);
printf(" Hamiltonian %d-paths.\n",
m);
fflush(stdout);
}
}

/*:53*/
#line 98 "dynahamp.w"
;
/*20:*/
#line 467 "dynahamp.w"

void compress(int l,unsigned long long p,int d){

register int j,k,kk,kkk;
register unsigned long long bits;
if(l==q)/*22:*/
#line 504 "dynahamp.w"

{
for(k= 0;k<prec;k++)
oldweight[owp].val[k]= weight[p-1].val[k];
owp++;
}

/*:22*/
#line 472 "dynahamp.w"

else{
for(j= (l+oms==q?1:-1),k= 0,bits= 0;k<d;j++,k++)
if(mem[p+j])bits+= 1LL<<k;
/*21:*/
#line 494 "dynahamp.w"

if(spack+d> bitsperword){
oldmem[omp++]= pack,pack= bits,spack= d;
if(omp>=oldmemsize){
fprintf(stderr,"Oops: oldmem overflow (more than %d bytes)!\n",
oldmemsize);
exit(-666666);
}
}else pack+= bits<<spack,spack+= d;

/*:21*/
#line 476 "dynahamp.w"
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

/*:20*//*24:*/
#line 535 "dynahamp.w"

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
#line 960 "dynahamp.w"

{
/*33:*/
#line 701 "dynahamp.w"

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

/*:33*/
#line 962 "dynahamp.w"
;
/*48:*/
#line 970 "dynahamp.w"

if(oldp==0)/*49:*/
#line 976 "dynahamp.w"

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
#line 971 "dynahamp.w"

else/*50:*/
#line 989 "dynahamp.w"

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
#line 972 "dynahamp.w"
;

/*:48*/
#line 963 "dynahamp.w"
;
/*36:*/
#line 788 "dynahamp.w"

for(j= 1;j<=q0;j++)bmate[j]= mp[1+oldmate[imap[j]]];
for(;j<=q;j++)bmate[j]= 0;

/*:36*/
#line 964 "dynahamp.w"
;
if(oldmate[2]<0)/*37:*/
#line 792 "dynahamp.w"

{
for(j= 1;j<=q;j++)mate[j]= bmate[j];
contribute();
}

/*:37*/
#line 965 "dynahamp.w"

else if(oldmate[2]==0)/*41:*/
#line 837 "dynahamp.w"

{
for(i= 0;i<r;i++)for(ii= i+1;ii<r;ii++){
for(j= 1;j<=q;j++)mate[j]= bmate[j];
if(add_derived(nbr[i],nbr[ii]))contribute();
}
}

/*:41*/
#line 966 "dynahamp.w"

else/*42:*/
#line 854 "dynahamp.w"

{
for(i= 0;i<r;i++){
for(j= 1;j<=q;j++)mate[j]= bmate[j];
if(add_derived(mp[1+oldmate[2]],nbr[i]))contribute();
}
}

/*:42*/
#line 967 "dynahamp.w"
;
}

/*:47*/
#line 540 "dynahamp.w"

oldp++;
}else{
/*25:*/
#line 560 "dynahamp.w"

if(spack+d> bitsperword)
pack= oldmem[++omp],spack= 0;
bits= pack>>spack,spack+= d;

/*:25*/
#line 543 "dynahamp.w"
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

/*:24*/
#line 99 "dynahamp.w"
;
int main(int argc,char*argv[]){
register int c,d,i,j,k,l,p,t,x,ii,ik,iv;
register Graph*g;
register Arc*a;
/*3:*/
#line 116 "dynahamp.w"

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
if(n>=maxn){
fprintf(stderr,"Recompile me: I allow at most %d vertices!\n",
maxn-1);
exit(-3);
}
(g->vertices+n)->name= gb_save_string("!");
printf("Dynamic Hamiltonian paths of the graph %s",
g->id);
printf(" (%d vertices, %ld edges):\n",
n,g->m/2);
fflush(stdout);

/*:3*/
#line 104 "dynahamp.w"
;
/*6:*/
#line 160 "dynahamp.w"

one.val[0]= 1;
for(k= 0;k<maxprec;k++)infty.val[k]= radix-1;
prec= 1;

/*:6*//*13:*/
#line 346 "dynahamp.w"

mem= (unsigned long long*)malloc(memsize*sizeof(unsigned long long));
oldmem= (unsigned long long*)malloc(oldmemsize);
weight= (bignum*)malloc(wtsize*sizeof(bignum));
oldweight= (bignum*)malloc(wtsize*sizeof(bignum));
if(!mem||!oldmem||!weight||!oldweight){
fprintf(stderr,"I can't allocate the big tables!\n");
exit(-6);
}

/*:13*//*27:*/
#line 609 "dynahamp.w"

fr[0]= n+1,ifr[n+1]= 0;
for(k= 1;k<=n;k++)fr[k]= ifr[k]= k;
q= 2;
mem[0]= mem[2]= 0,mem[1]= 1,memptr= 3;
weight[0]= one,wtptr= 1;

/*:27*//*35:*/
#line 777 "dynahamp.w"

mp[0]= -1,mp[2]= 1;
imap[1]= 1;

/*:35*//*39:*/
#line 819 "dynahamp.w"

nbr[0]= 1;

/*:39*/
#line 105 "dynahamp.w"
;
for(m= 1;;m++){
/*31:*/
#line 669 "dynahamp.w"

fprintf(stderr,"\nThe frontier for %d-classes is",
m-1);
for(k= 0;k<q;k++)
fprintf(stderr," %s",
vert(fr[k])->name);
fprintf(stderr,".\n");

/*:31*/
#line 107 "dynahamp.w"
;
if(wtptr==0)/*52:*/
#line 1027 "dynahamp.w"

{
fprintf(stderr,"\n");
for(k= m;k<=n;k++){
if(k> m+m)break;
report_cycles(k);
}
fprintf(stderr,"\nThat's all; there are no %d-configs!\n",
m-1);
fprintf(stderr,
"Storage requirements: %lld memsize, %d deg, %lld wtsize, %d maxprec.\n",
maxmemptr,maxdeg+1,maxwtptr,18*prec);
exit(0);
}

/*:52*/
#line 108 "dynahamp.w"
;
/*23:*/
#line 511 "dynahamp.w"

fprintf(stderr,"There are %lld %d-classes,\n",
wtptr,m-1);
fprintf(stderr," resulting from %lld contributions and filling %lld pointers.\n",
contribs,memptr);
spack= omp= owp= oms= 0,omx= 1;
pack= 0;
oldq= q;
compress(0,1,3);
oldmem[omp]= pack;
if(omp> maxomp)maxomp= omp;

/*:23*/
#line 109 "dynahamp.w"
;
/*51:*/
#line 1015 "dynahamp.w"

if(memptr> maxmemptr)maxmemptr= memptr;
if(wtptr> maxwtptr)maxwtptr= wtptr;
mem[0]= mem[1]= mem[2]= 0,memptr= 3,wtptr= 0;

/*28:*/
#line 616 "dynahamp.w"

for(j= 1;j<q;j++)ofr[j]= fr[j];
oldq= q;
iv= ifr[m+1];
if(iv<--q){
x= fr[q];
fr[1]= m+1,ifr[m+1]= 1;
fr[q]= m,ifr[m]= q;
fr[iv]= x,ifr[x]= iv;
}else{
fr[1]= m+1,ifr[m+1]= 1,fr[iv]= m,ifr[m]= iv;
if(iv!=q)q++;
}
q0= q;
/*29:*/
#line 637 "dynahamp.w"

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

/*:29*/
#line 630 "dynahamp.w"
;
/*30:*/
#line 649 "dynahamp.w"

for(k= 3;k<q0;k++)if(fr[k]<fr[k-1]){
for(t= fr[k],j= k-1;;j--){
fr[j+1]= fr[j],ifr[fr[j]]= j+1;
if(j==1||fr[j-1]<t)break;
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

/*:30*/
#line 631 "dynahamp.w"
;
/*40:*/
#line 822 "dynahamp.w"

for(r= 1,a= vert(m)->arcs;a;a= a->next){
k= vertnum(a->tip);
if(k<m)continue;
nbr[r++]= ifr[k]+1;
}

/*:40*//*43:*/
#line 866 "dynahamp.w"

imap[2]= 0;
for(j= 3;j<=oldq;j++){
mp[j+1]= 1+ifr[ofr[j-1]];
imap[mp[j+1]]= j;
}

/*:43*/
#line 632 "dynahamp.w"
;

/*:28*/
#line 1020 "dynahamp.w"
;

/*:51*/
#line 110 "dynahamp.w"
;
/*46:*/
#line 956 "dynahamp.w"

spack= oldp= omp= oms= 0,omx= 1,pack= oldmem[0];
uncompress(0,3);

/*:46*/
#line 111 "dynahamp.w"
;
report_cycles(m);
}
}

/*:2*/
