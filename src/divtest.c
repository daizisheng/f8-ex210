/*3:*/
#line 42 "divtest.w"

#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include <flint/fmpz.h> 
#include <flint/fmpz_mod.h> 
#include <flint/fmpz_mod_poly.h> 

/*4:*/
#line 60 "divtest.w"

#ifndef MODP
#define MODP 1000000007
#endif
#line 64 "divtest.w"

/*:4*/
#line 50 "divtest.w"


/*5:*/
#line 74 "divtest.w"

static int read_poly(const char*path,fmpz_mod_poly_t out,fmpz_mod_ctx_t ctx){
FILE*f= fopen(path,"r");
if(!f){perror(path);return-1;}
char line[256];
fmpz_mod_poly_zero(out,ctx);
while(fgets(line,sizeof(line),f)){
if(line[0]=='#'||line[0]=='\n')continue;
int i;unsigned long long v;
if(sscanf(line,"%d %llu",&i,&v)==2){
fmpz_t fv;
fmpz_init_set_ui(fv,v);
fmpz_mod_poly_set_coeff_fmpz(out,i,fv,ctx);
fmpz_clear(fv);
}
}
fclose(f);
return 0;
}

/*:5*/
#line 52 "divtest.w"


/*6:*/
#line 99 "divtest.w"

int main(int argc,char**argv){
if(argc<4){
fprintf(stderr,
"Usage: %s <Q.txt> <Q+.txt> <power>\n"
"  e.g. %s Q5.txt Q5_plus.txt 3\n",
argv[0],argv[0]);
return 1;
}
int k= atoi(argv[3]);

fmpz_mod_ctx_t ctx;
fmpz_t modulus;
fmpz_init_set_ui(modulus,MODP);
fmpz_mod_ctx_init(ctx,modulus);

fmpz_mod_poly_t Q,Qplus,Q_pow,quot,rem;
fmpz_mod_poly_init(Q,ctx);
fmpz_mod_poly_init(Qplus,ctx);
fmpz_mod_poly_init(Q_pow,ctx);
fmpz_mod_poly_init(quot,ctx);
fmpz_mod_poly_init(rem,ctx);

if(read_poly(argv[1],Q,ctx)<0)return 1;
if(read_poly(argv[2],Qplus,ctx)<0)return 1;

slong dQ= fmpz_mod_poly_degree(Q,ctx);
slong dQplus= fmpz_mod_poly_degree(Qplus,ctx);
fprintf(stderr,"deg Q = %ld, deg Q+ = %ld, testing Q^%d divides Q+\n",
dQ,dQplus,k);

fmpz_mod_poly_pow(Q_pow,Q,(ulong)k,ctx);
slong dPow= fmpz_mod_poly_degree(Q_pow,ctx);
fprintf(stderr,"deg Q^%d = %ld\n",k,dPow);

fmpz_mod_poly_divrem(quot,rem,Qplus,Q_pow,ctx);
if(fmpz_mod_poly_is_zero(rem,ctx)){
printf("DIVISIBLE: Q^%d divides Q+\n",k);
}else{
slong dR= fmpz_mod_poly_degree(rem,ctx);
printf("NOT DIVISIBLE: Q^%d does not divide Q+ "
"(remainder degree %ld)\n",k,dR);
}

fmpz_mod_poly_clear(Q,ctx);
fmpz_mod_poly_clear(Qplus,ctx);
fmpz_mod_poly_clear(Q_pow,ctx);
fmpz_mod_poly_clear(quot,ctx);
fmpz_mod_poly_clear(rem,ctx);
fmpz_mod_ctx_clear(ctx);
fmpz_clear(modulus);
return 0;
}

/*:6*/
#line 54 "divtest.w"


/*:3*/
