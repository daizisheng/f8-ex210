/*11:*/
#line 175 "simulate_count.w"

#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include <time.h> 
#include <stdint.h> 
#include <assert.h> 

/*12:*/
#line 205 "simulate_count.w"

#ifdef _OPENMP
#include <omp.h> 
#endif
#line 209 "simulate_count.w"

/*:12*/
#line 184 "simulate_count.w"


/*13:*/
#line 214 "simulate_count.w"

#ifndef MODP
#define MODP 1000000007ULL
#endif
#line 218 "simulate_count.w"
#define MAX_LINE 1024
#define MAX_M_STEPS 64
#define MAX_STATE_LEN 64

/*:13*/
#line 186 "simulate_count.w"


/*14:*/
#line 226 "simulate_count.w"

typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

/*:14*/
#line 188 "simulate_count.w"


/*15:*/
#line 233 "simulate_count.w"

static double now_sec(void){
struct timespec ts;
clock_gettime(CLOCK_MONOTONIC,&ts);
return ts.tv_sec+ts.tv_nsec*1e-9;
}

/*:15*//*16:*/
#line 249 "simulate_count.w"

typedef struct{char*str;u32 idx;}entry;
typedef struct{
entry*tab;
int bits;
int mask;
int n;
char**all_strs;
int all_cap;
}smap;

static u64 fnv1a(const char*s){
u64 h= 0xcbf29ce484222325ULL;
while(*s){h^= (unsigned char)*s++;h*= 0x100000001b3ULL;}
return h;
}

/*:16*//*17:*/
#line 266 "simulate_count.w"

static void smap_init(smap*s){
s->bits= 4;
s->mask= (1<<s->bits)-1;
s->tab= calloc(1<<s->bits,sizeof(entry));
s->n= 0;
s->all_strs= NULL;
s->all_cap= 0;
}

static void smap_grow(smap*s){
int old_size= 1<<s->bits;
entry*old= s->tab;
s->bits++;
s->mask= (1<<s->bits)-1;
s->tab= calloc(1<<s->bits,sizeof(entry));
for(int i= 0;i<old_size;i++)if(old[i].str){
u32 h= (u32)(fnv1a(old[i].str)&s->mask);
while(s->tab[h].str)h= (h+1)&s->mask;
s->tab[h]= old[i];
}
free(old);
}

/*:17*//*18:*/
#line 293 "simulate_count.w"

static u32 smap_add(smap*s,const char*str){
if((s->n+1)*2> (1<<s->bits))smap_grow(s);
u32 h= (u32)(fnv1a(str)&s->mask);
while(s->tab[h].str){
if(strcmp(s->tab[h].str,str)==0)return s->tab[h].idx;
h= (h+1)&s->mask;
}
s->tab[h].str= strdup(str);
s->tab[h].idx= (u32)s->n;
if(s->n==s->all_cap){
s->all_cap= s->all_cap?s->all_cap*2:16;
s->all_strs= realloc(s->all_strs,s->all_cap*sizeof(char*));
}
s->all_strs[s->n]= s->tab[h].str;
s->n++;
return(u32)(s->n-1);
}

/*:18*//*19:*/
#line 312 "simulate_count.w"

static u32 smap_lookup(smap*s,const char*str){
u32 h= (u32)(fnv1a(str)&s->mask);
while(s->tab[h].str){
if(strcmp(s->tab[h].str,str)==0)return s->tab[h].idx;
h= (h+1)&s->mask;
}
return UINT32_MAX;
}

/*:19*//*20:*/
#line 326 "simulate_count.w"

static int cmpstrptr(const void*a,const void*b){
return strcmp(*(const char**)a,*(const char**)b);
}

static void smap_finalize(smap*s){
qsort(s->all_strs,s->n,sizeof(char*),cmpstrptr);
for(int i= 0;i<s->n;i++){
u32 h= (u32)(fnv1a(s->all_strs[i])&s->mask);
while(s->tab[h].str&&strcmp(s->tab[h].str,s->all_strs[i])!=0)
h= (h+1)&s->mask;
s->tab[h].idx= (u32)i;
}
}

/*:20*//*21:*/
#line 341 "simulate_count.w"

static void smap_free(smap*s){
if(s->tab){
int sz= 1<<s->bits;
for(int i= 0;i<sz;i++)free(s->tab[i].str);
free(s->tab);
}
free(s->all_strs);
s->tab= NULL;s->all_strs= NULL;s->n= 0;
}

/*:21*/
#line 190 "simulate_count.w"


/*22:*/
#line 357 "simulate_count.w"

typedef struct{
int n_src;
int n_dst;
u32*col_ptr;
u32*src;
}bin_csc;

static void bin_csc_init(bin_csc*m,int n_src,int n_dst){
m->n_src= n_src;
m->n_dst= n_dst;
m->col_ptr= calloc(n_dst+1,sizeof(u32));
m->src= NULL;
}

static void bin_csc_free(bin_csc*m){
free(m->col_ptr);
free(m->src);
m->col_ptr= NULL;m->src= NULL;
}

/*:22*//*23:*/
#line 388 "simulate_count.w"

static void bin_csc_apply(const bin_csc*m,const u32*v,u32*y){
int n= m->n_dst;
const u32*col_ptr= m->col_ptr;
const u32*src= m->src;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
#line 396 "simulate_count.w"
 for(int j= 0;j<n;j++){
u64 acc= 0;
u32 lo= col_ptr[j],hi= col_ptr[j+1];
for(u32 k= lo;k<hi;k++)
acc+= (u64)v[src[k]];
y[j]= (u32)(acc%MODP);
}
}

/*:23*/
#line 192 "simulate_count.w"


/*24:*/
#line 412 "simulate_count.w"

typedef struct{
int n;
int cap;
u32*src;
u32*mult;
}cycle_list;

static void cycle_list_init(cycle_list*c){
c->n= 0;c->cap= 0;c->src= NULL;c->mult= NULL;
}

static void cycle_list_push(cycle_list*c,u32 s,u32 m){
if(c->n==c->cap){
c->cap= c->cap?c->cap*2:8;
c->src= realloc(c->src,c->cap*sizeof(u32));
c->mult= realloc(c->mult,c->cap*sizeof(u32));
}
c->src[c->n]= s;
c->mult[c->n]= m;
c->n++;
}

static void cycle_list_free(cycle_list*c){
free(c->src);free(c->mult);
}

/*:24*/
#line 194 "simulate_count.w"


/*25:*/
#line 451 "simulate_count.w"

static u64 mod_inv(u64 a,u64 p){
u64 result= 1,base= a%p,exp= p-2;
while(exp> 0){
if(exp&1)result= result*base%p;
base= base*base%p;
exp>>= 1;
}
return result;
}

static int berlekamp_massey(const u64*s,int n,u64*c_out,u64 p){
u64*c= calloc(n+1,sizeof(u64));
u64*b= calloc(n+1,sizeof(u64));
int big_l= 0,m= 1;
u64 last_b= 1;
c[0]= 1;b[0]= 1;
int c_len= 1,b_len= 1;
for(int i= 0;i<n;i++){
u64 d= s[i]%p;
for(int j= 1;j<=big_l;j++)
d= (d+c[j]*s[i-j])%p;
if(d==0){m++;continue;}
if(2*big_l<=i){
u64*t= malloc(c_len*sizeof(u64));
memcpy(t,c,c_len*sizeof(u64));
int t_len= c_len;
u64 coef= d*mod_inv(last_b,p)%p;
while(c_len<b_len+m){c[c_len++]= 0;}
for(int j= 0;j<b_len;j++)
c[j+m]= (c[j+m]+p-coef*b[j]%p)%p;
big_l= i+1-big_l;
free(b);b= t;b_len= t_len;
last_b= d;m= 1;
}else{
u64 coef= d*mod_inv(last_b,p)%p;
while(c_len<b_len+m){c[c_len++]= 0;}
for(int j= 0;j<b_len;j++)
c[j+m]= (c[j+m]+p-coef*b[j]%p)%p;
m++;
}
}
memcpy(c_out,c,(big_l+1)*sizeof(u64));
free(c);free(b);
return big_l;
}

/*:25*/
#line 196 "simulate_count.w"


/*26:*/
#line 509 "simulate_count.w"

typedef struct{int m_offset;char old_s[MAX_STATE_LEN];char new_s[MAX_STATE_LEN];}trec;
typedef struct{int m_offset;int cycle_at;char old_s[MAX_STATE_LEN];}crec;
typedef struct{int m_offset;char code[MAX_STATE_LEN];u64 w;}wrec;

int main(int argc,char**argv){
/*27:*/
#line 538 "simulate_count.w"

if(argc<5){
fprintf(stderr,
"Usage: %s <dump.log> <max_iters> <target_cycle_at> <m_start>\n",
argv[0]);
return 1;
}
const char*path= argv[1];
int max_iters= atoi(argv[2]);
int target_cycle_at= atoi(argv[3]);
int m_start= atoi(argv[4]);
int period= 10;
if(getenv("PERIOD"))period= atoi(getenv("PERIOD"));
fprintf(stderr,
"Path: %s, max_iters: %d, target_cycle_at: %d, m_start: %d, period: %d\n",
path,max_iters,target_cycle_at,m_start,period);
#ifdef _OPENMP
fprintf(stderr,"OpenMP threads: %d\n",omp_get_max_threads());
#else
#line 557 "simulate_count.w"
 fprintf(stderr,"Build: single-threaded (compile with -fopenmp for parallel SpMV)\n");
#endif
#line 559 "simulate_count.w"

/*:27*/
#line 515 "simulate_count.w"
;
/*28:*/
#line 560 "simulate_count.w"

trec*T= NULL;int nT= 0,capT= 0;
crec*C= NULL;int nC= 0,capC= 0;
wrec*W= NULL;int nW= 0,capW= 0;
smap chained_in[MAX_M_STEPS];
for(int i= 0;i<MAX_M_STEPS;i++)smap_init(&chained_in[i]);

/*:28*/
#line 516 "simulate_count.w"
;
/*29:*/
#line 571 "simulate_count.w"

{
double t0= now_sec();
FILE*f= fopen(path,"r");
if(!f){perror("fopen dump");return 1;}
char line[MAX_LINE*4];
while(fgets(line,sizeof(line),f)){
if(line[0]=='T'&&line[1]==' ')/*30:*/
#line 588 "simulate_count.w"
{
int m_step;
char old_s[MAX_STATE_LEN],new_s[MAX_STATE_LEN];
if(sscanf(line,"T m=%d old=%63s new=%63s",&m_step,old_s,new_s)==3){
int off= m_step-m_start;
if(off>=0&&off<MAX_M_STEPS){
if(nT==capT){
capT= capT?capT*2:1024;
T= realloc(T,capT*sizeof(trec));
}
T[nT].m_offset= off;
strcpy(T[nT].old_s,old_s);
strcpy(T[nT].new_s,new_s);
nT++;
}
}
}

/*:30*/
#line 578 "simulate_count.w"

else if(line[0]=='C'&&line[1]==' ')/*31:*/
#line 606 "simulate_count.w"
{
int m_step,ca;
char old_s[MAX_STATE_LEN];
if(sscanf(line,"C m=%d cycle_at=%d old=%63s",&m_step,&ca,old_s)==3){
int off= m_step-m_start;
if(off>=1&&off<MAX_M_STEPS&&ca==target_cycle_at){
if(nC==capC){
capC= capC?capC*2:64;
C= realloc(C,capC*sizeof(crec));
}
C[nC].m_offset= off;
C[nC].cycle_at= ca;
strcpy(C[nC].old_s,old_s);
nC++;
}
}
}

/*:31*/
#line 579 "simulate_count.w"

else if(line[0]=='W'&&line[1]==' ')/*32:*/
#line 624 "simulate_count.w"
{
int m_step,oldp_dummy;
char code[MAX_STATE_LEN];
unsigned long long w;
if(sscanf(line,"W m=%d oldp=%d code=%63s w=%llu",
&m_step,&oldp_dummy,code,&w)==4){
if(m_step==m_start){
if(nW==capW){
capW= capW?capW*2:64;
W= realloc(W,capW*sizeof(wrec));
}
W[nW].m_offset= 0;
strcpy(W[nW].code,code);
W[nW].w= w;
nW++;
}
}
}

/*:32*/
#line 580 "simulate_count.w"

}
fclose(f);
fprintf(stderr,
"Parsed in %.2fs:  T=%d, C(cycle_at=%d)=%d, W(m_start)=%d\n",
now_sec()-t0,nT,target_cycle_at,nC,nW);
}

/*:29*/
#line 517 "simulate_count.w"
;
/*33:*/
#line 651 "simulate_count.w"

{
for(int t= 0;t<nT;t++){
int k= T[t].m_offset-1;
if(k<0||k+1> period)continue;
smap_add(&chained_in[k],T[t].old_s);
smap_add(&chained_in[k+1],T[t].new_s);
}
for(int k= 0;k<=period;k++)smap_finalize(&chained_in[k]);
fprintf(stderr,"chained_in sizes:");
for(int k= 0;k<=period;k++)
fprintf(stderr," [%d]=%d",k,chained_in[k].n);
fprintf(stderr,"\n");
}

/*:33*/
#line 518 "simulate_count.w"
;
/*34:*/
#line 666 "simulate_count.w"

bin_csc M[MAX_M_STEPS];
{
for(int k= 0;k<period;k++){
int n_src= chained_in[k].n;
int n_dst= chained_in[k+1].n;
bin_csc_init(&M[k],n_src,n_dst);
/*35:*/
#line 684 "simulate_count.w"

u32*col_count= calloc(n_dst,sizeof(u32));
int tmp_cap= 0,tmp_n= 0;
u32*tmp_src= NULL,*tmp_dst= NULL;
for(int t= 0;t<nT;t++){
if(T[t].m_offset!=k+1)continue;
u32 s= smap_lookup(&chained_in[k],T[t].old_s);
u32 d= smap_lookup(&chained_in[k+1],T[t].new_s);
if(s==UINT32_MAX||d==UINT32_MAX)continue;
if(tmp_n==tmp_cap){
tmp_cap= tmp_cap?tmp_cap*2:64;
tmp_src= realloc(tmp_src,tmp_cap*sizeof(u32));
tmp_dst= realloc(tmp_dst,tmp_cap*sizeof(u32));
}
tmp_src[tmp_n]= s;tmp_dst[tmp_n]= d;
tmp_n++;
col_count[d]++;
}

/*:35*/
#line 673 "simulate_count.w"
;
/*36:*/
#line 703 "simulate_count.w"

M[k].col_ptr[0]= 0;
for(int j= 0;j<n_dst;j++)
M[k].col_ptr[j+1]= M[k].col_ptr[j]+col_count[j];
M[k].src= malloc(tmp_n*sizeof(u32));
u32*cursor= calloc(n_dst,sizeof(u32));
for(int t= 0;t<tmp_n;t++){
int d= tmp_dst[t];
M[k].src[M[k].col_ptr[d]+cursor[d]]= tmp_src[t];
cursor[d]++;
}
free(cursor);free(col_count);
free(tmp_src);free(tmp_dst);

/*:36*/
#line 674 "simulate_count.w"
;
/*37:*/
#line 720 "simulate_count.w"

{
int cmp_u32(const void*a,const void*b){
u32 x= *(const u32*)a,y= *(const u32*)b;
return(x> y)-(x<y);
}
for(int j= 0;j<n_dst;j++){
u32 lo= M[k].col_ptr[j],hi= M[k].col_ptr[j+1];
if(hi-lo> 1)
qsort(&M[k].src[lo],hi-lo,sizeof(u32),cmp_u32);
}
}

/*:37*/
#line 675 "simulate_count.w"
;
}
u64 total_nnz= 0;
for(int k= 0;k<period;k++)
total_nnz+= M[k].col_ptr[M[k].n_dst];
fprintf(stderr,"Built %d matrices, total nnz = %llu\n",
period,(unsigned long long)total_nnz);
}

/*:34*/
#line 519 "simulate_count.w"
;
/*38:*/
#line 733 "simulate_count.w"

cycle_list cycles_by_offset[MAX_M_STEPS];
for(int k= 0;k<=period;k++)cycle_list_init(&cycles_by_offset[k]);
{
int cycle_count= 0;
for(int t= 0;t<nC;t++){
int off= C[t].m_offset;
if(off<1||off> period)continue;
u32 s= smap_lookup(&chained_in[off-1],C[t].old_s);
if(s==UINT32_MAX)continue;
cycle_list_push(&cycles_by_offset[off],s,1);
cycle_count++;
}
fprintf(stderr,"Cycle events captured: %d\n",cycle_count);
}

/*:38*/
#line 520 "simulate_count.w"
;
/*39:*/
#line 752 "simulate_count.w"

int dim0= chained_in[0].n;
u32*v= calloc(dim0,sizeof(u32));
for(int t= 0;t<nW;t++){
u32 i= smap_lookup(&chained_in[0],W[t].code);
if(i!=UINT32_MAX)v[i]= (u32)(W[t].w%MODP);
}
u32*remap_post_to_pre= malloc(chained_in[period].n*sizeof(u32));
for(int i= 0;i<chained_in[period].n;i++)
remap_post_to_pre[i]= UINT32_MAX;
for(int i= 0;i<chained_in[period].n;i++){
const char*s= chained_in[period].all_strs[i];
u32 pre_i= smap_lookup(&chained_in[0],s);
remap_post_to_pre[i]= pre_i;
}

/*:39*/
#line 521 "simulate_count.w"
;
/*40:*/
#line 776 "simulate_count.w"

int max_dim= 0;
for(int k= 0;k<=period;k++)
if(chained_in[k].n> max_dim)max_dim= chained_in[k].n;
u32*buf_a= calloc(max_dim,sizeof(u32));
u32*buf_b= calloc(max_dim,sizeof(u32));
memcpy(buf_a,v,dim0*sizeof(u32));
u64*sequence= malloc(max_iters*sizeof(u64));
double main_t= now_sec();
for(int iter= 0;iter<max_iters;iter++){
u64 count_iter= 0;
for(int s_idx= 0;s_idx<period;s_idx++){
int offset= s_idx+1;
/*41:*/
#line 801 "simulate_count.w"

{
cycle_list*cl= &cycles_by_offset[offset];
for(int k= 0;k<cl->n;k++){
u32 vi= buf_a[cl->src[k]];
u32 mult= cl->mult[k];
u64 contrib= (mult==1)?vi:(u64)vi*mult%MODP;
count_iter= (count_iter+contrib)%MODP;
}
}

/*:41*/
#line 789 "simulate_count.w"
;
bin_csc_apply(&M[s_idx],buf_a,buf_b);
u32*tmp= buf_a;buf_a= buf_b;buf_b= tmp;
}
sequence[iter]= count_iter;
/*42:*/
#line 812 "simulate_count.w"

{
memset(v,0,dim0*sizeof(u32));
int post_n= chained_in[period].n;
for(int i= 0;i<post_n;i++){
u32 r= remap_post_to_pre[i];
if(r!=UINT32_MAX){
u64 s= (u64)v[r]+(u64)buf_a[i];
v[r]= (u32)(s>=MODP?s-MODP:s);
}
}
memcpy(buf_a,v,dim0*sizeof(u32));
}

/*:42*/
#line 794 "simulate_count.w"
;
/*43:*/
#line 826 "simulate_count.w"

if(iter<5||iter%1000==0){
fprintf(stderr,"  iter %d: count=%llu, elapsed %.2fs\n",
iter,(unsigned long long)count_iter,now_sec()-main_t);
}

/*:43*/
#line 795 "simulate_count.w"
;
}
double elapsed= now_sec()-main_t;
fprintf(stderr,"\nMain loop: %d iters in %.2fs = %.3f ms/iter\n",
max_iters,elapsed,1000.0*elapsed/max_iters);

/*:40*/
#line 522 "simulate_count.w"
;
/*44:*/
#line 832 "simulate_count.w"

{
FILE*out= fopen("/tmp/simulate_count.out","w");
if(out){
fprintf(out,"# Mod p sequence, p=%llu, count=%d\n",
(unsigned long long)MODP,max_iters);
for(int i= 0;i<max_iters;i++)
fprintf(out,"%llu\n",(unsigned long long)sequence[i]);
fclose(out);
fprintf(stderr,"Wrote %d entries to /tmp/simulate_count.out\n",max_iters);
}
}

/*:44*/
#line 523 "simulate_count.w"
;
/*45:*/
#line 848 "simulate_count.w"

{
u64*poly= calloc(max_iters+1,sizeof(u64));
double bm_t= now_sec();
int L= berlekamp_massey(sequence,max_iters,poly,MODP);
fprintf(stderr,"BM done in %.3fs. Order L = %d, surplus = %lld\n",
now_sec()-bm_t,L,(long long)max_iters-2LL*L);
FILE*pf= fopen("/tmp/simulate_count.poly","w");
if(pf){
fprintf(pf,"# min poly mod %llu, degree %d\n",
(unsigned long long)MODP,L);
for(int i= 0;i<=L;i++)
fprintf(pf,"%d %llu\n",i,(unsigned long long)poly[i]);
fclose(pf);
fprintf(stderr,"Wrote polynomial of degree %d to /tmp/simulate_count.poly\n",L);
}
free(poly);
}

/*:45*/
#line 524 "simulate_count.w"
;
/*46:*/
#line 867 "simulate_count.w"

free(v);free(buf_a);free(buf_b);
free(sequence);free(remap_post_to_pre);
free(T);free(C);free(W);
for(int k= 0;k<period;k++)bin_csc_free(&M[k]);
for(int k= 0;k<=period;k++)cycle_list_free(&cycles_by_offset[k]);
for(int k= 0;k<=period;k++)smap_free(&chained_in[k]);

/*:46*/
#line 525 "simulate_count.w"
;
return 0;
}

/*:26*/
#line 198 "simulate_count.w"


/*:11*/
