
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#define NDEBUG 1
// #include <boost/math/special_functions/factorials.hpp>
// #include <boost/math/distributions/poisson.hpp>
// #include <boost/random.hpp>
// #include <boost/random/uniform_real_distribution.hpp>
#include <math.h>
#include <cmath>
#include <random>
using namespace Rcpp;


std::random_device rd;
std::default_random_engine generator(rd());

// [[Rcpp::export]]
float Poiss(float m, int x) {
  float dens;
  if(x==0){
    dens = exp(-m);
  }
  else if(x<30){
    float fact=1;
    for(int i=1;i<=x;i++){
      fact=fact/i;
    }
    dens = exp(-m)*fact*pow(m,x);
  }else{
    
    // float fact=0;
    // for(int i=1;i<=x;i++){
    //   fact+=log(i);
    // }
    // float ldens = -m+x*log(m)-fact;
    // dens = exp(ldens);
    float pi = 3.14159265358979323846;
    float e = exp(1);
    dens = exp(-m)*pow(m*e/x,x)/sqrt(2*pi*x);
  }
  return dens;
}

// [[Rcpp::export]]
int discsamp(arma::vec p){
  // std::random_device rd;
  // std::default_random_engine generator(rd());
  std::uniform_real_distribution<float> uniform_distribution(0.0,1.0);
  float uuu=uniform_distribution(generator);
  float sss=sum(p);
  if(sss==0){
    p=arma::ones(p.n_elem);
  }
  int x=0;
  for(unsigned int kk=0;kk<p.n_elem-1;kk++){
    p[kk+1]+=p[kk];
    if(p[kk]/sss<uuu){
      x=kk+1;
    }
  }
  return x;
}

// [[Rcpp::export]]
arma::mat Fast_Matrix_exponential(arma::mat Q,
                                  arma::vec Lambda,
                                  float Inter,
                                  int M2,
                                  int pp=10,
                                  bool msg=false) {
  //get dimension of Q
  int M = Q.n_rows;
  //find the right scale for scaling and squaring
  int scale = std::ceil(3*max(-2*Q.diag()+2*Lambda)*Inter);
  arma::mat A = (Q-diagmat(Lambda))*Inter/scale;
  arma::mat B = diagmat(Lambda*Inter/scale);
  
  arma::mat exptemp = arma::zeros(M*pp,M*(pp+1));
  exptemp.submat(0,0,M-1,2*M-1)=arma::join_rows(A,B);
  
  
  for(int ii=1;ii<=pp-1;ii++){
    exptemp.submat(M*ii,0,M+M*ii-1,M-1)=
      exptemp.submat(M*(ii-1),0,M+M*(ii-1)-1,M-1)*A;
    for(int iii=1;iii<=ii+1;iii++){
      exptemp.submat(ii*M,iii*M,ii*M-1+M,M+iii*M-1) =
        exptemp.submat((ii-1)*M,iii*M,(ii-1)*M-1+M,M+iii*M-1)*A+
        exptemp.submat((ii-1)*M,(iii-1)*M,(ii-1)*M-1+M,M+(iii-1)*M-1)*B;
    }
  }
  
  exptemp.submat(0,0,M-1,M-1) += arma::eye(M,M);
  arma::mat expest = exptemp.submat(0,0,M-1,M*(pp+1)-1);
  float fact = 1;
  for(int ii=2; ii<=pp;ii++){
    fact=fact/ii;
    expest+=exptemp.submat((ii-1)*M,0,(ii-1)*M+M-1,M*(pp+1)-1)*fact; //boost::math::factorial<double>(ii);
  }
  
  int tempdim = M2;
  if(M*(pp+1)<M2){
    tempdim=M*(pp+1);
  }
  
  arma::mat exptot = arma::zeros(2*M,M2);
  exptot.submat(0,0,M-1,tempdim-1)=expest.submat(0,0,M-1,tempdim-1);
  int tempk=0;
  int k1 = 0;
  int k2 = 0;
  int progress=0;
  for(int ii=1;ii<scale;ii++){
    if(msg==true){
      if(ii>=progress*(scale/20.0)){
        Rcpp::Rcout << "\r  Interval sampling, Initialisation:"<<(progress)*5<<"%                    ";  
        progress++;
      }
    }
    k1 = ii%2;
    k2 = (ii+1)%2;
    
    for(int j=1; j<=M2/M; j++){
      if(pp+1<j){
        tempk = pp+1;
      }else{
        tempk = j;
      }
      exptot.submat(k1*M,(j-1)*M,M+k1*M-1,M+(j-1)*M-1)=arma::zeros(M,M);
      for(int l=1;l<=tempk;l++){
        exptot.submat(k1*M,(j-1)*M,M+k1*M-1,M+(j-1)*M-1)+=
          exptot.submat(k2*M,(j-l)*M,M+k2*M-1,M+(j-l)*M-1)*expest.submat(0,(l-1)*M,M-1,(l-1)*M+M-1);
      }
    }
  }
  
  arma::mat P=exptot.submat(k1*M,0,k1*M+M-1,M2-1);
  
  //int P=1;
  
  return P;
}

// [[Rcpp::export]]
List Forward_Accumulation_cpp(arma::ivec y_0Tprime,
                              float Inter,
                              arma::mat Q,
                              arma::vec Lambda,
                              arma::vec Rho,
                              int pp=12,
                              bool messages=true) {
  // Initialize the random number generators
  // boost::random::mt19937 generator( time(0) ); 
  // boost::random::uniform_real_distribution< >  uniform_distribution(0.0,1.0);
  // std::random_device rd;
  // std::default_random_engine generator(rd());
  std::uniform_real_distribution<double> uniform_distribution(0.0,1.0);
  if(messages==true){
    Rcpp::Rcout << "\r Backwards calculation: Initialisation                    "; 
  }
  //length of input
  int t = y_0Tprime.n_elem;
  //Initialize output
  arma::vec x_0T = arma::ones(t+1);
  int M = Q.n_rows; //Dim of state space
  int c_max = y_0Tprime.max(); //maximum observation
  // Generator for auxilliary V_t process
  
  /////////////////////////////////////////////////////////////////////
  // //Part 1: Calculate marginal posteriors recursively
  /////////////////////////////////////////////////////////////////////
  
  arma::mat P= Fast_Matrix_exponential(Q,Lambda,Inter,(c_max+1)*M,pp=10);
  // A matrix for marginal probability distributions
  arma::mat A=arma::zeros(M,M*(t+1));
  // Initialize A
  A.submat(0,M*t,M-1,M*(t+1)-1)=P.submat(0,0+M*y_0Tprime[t-1],M-1,M-1+M*y_0Tprime[t-1]);
  // A-normalizer variable
  float k=arma::accu(A.submat(0,M*t,M-1,M*(t+1)-1)); //arma::ones(1,A.n_cols)* A * arma::ones(A.n_cols, 1);
  // normalize A
  A.submat(0,M*t,M-1,M*(t+1)-1)=A.submat(0,M*t,M-1,M*(t+1)-1)/k;
  //messages
  int progress=20;
  for(int i = t; i >= 1; i--) {
    
    if(messages==true){
      if(i<=progress*(t/20.0)){
        Rcpp::Rcout << "\r Backwards calculation:"<<(20-progress)*5<<"%                    ";  
        progress--;
      }
    }
    // Calculate each A element recursively
    arma::mat A_now=P.submat(0,0+M*y_0Tprime[i-1],M-1,M-1+M*y_0Tprime[i-1])*
      A.submat(0,M*i,M-1,M*(i+1)-1);
    //Normalize
    float k=arma::accu(A_now); 
    A.submat(0,M*(i-1),M-1,M*(i)-1)=A_now/k;
  }
  //Part 2: Sample forward
  if(messages==true){
    Rcpp::Rcout << "\r Forwards sampling:                          ";
  }
  //Initialize M-dimensional probability vector
  arma::vec p=(Rho%(A.submat(0,0,M-1,M-1)*arma::ones(M)));
  //Normalize p
  p=p/arma::accu(p);
  
  // boost::random::discrete_distribution<> dist(p);
  //Draw first sample
  // x_0T[0]=dist(generator)+1;
  x_0T[0]=discsamp(p)+1;
  arma::vec help_vec=arma::ones(M);
  //Sample forward
  progress=0;
  for(int i = 1; i <= t-1; i++) {
    if(messages==true){
      if(i>=progress*(t/20.0)){
        Rcpp::Rcout << "\r Forwards sampling:"<<(progress)*5<<"%                    ";  
        progress++;
      }
    }
    for(int ll = 1; ll <= M-1; ll++) {
      arma::vec help_vec=A.submat(0,0+M*i,M-1,M-1+M*i) * arma::ones(M);
      p[ll-1]=P(x_0T[i-1]-1,M*y_0Tprime[i-1]+ll-1)*
        help_vec(ll-1)/help_vec(x_0T[i-1]-1);
    }
    // boost::random::discrete_distribution<> dist(p);
    // x_0T[i]=dist(generator)+1;
    x_0T[i]=discsamp(p)+1;
  }
  //Draw last sample
  int i=t;
  for(int ll = 1; ll <= M-1; ll++) {
    arma::vec help_vec=A.submat(0,0+M*i,M-1,M-1+M*i) * arma::ones(M);
    p[ll-1]=P(x_0T[i-1]-1,M*y_0Tprime[i-1]+ll-1)*
      help_vec(ll-1)/help_vec(x_0T[i-1]-1);
  }
  // boost::random::discrete_distribution<> dist_help(p);
  // x_0T[i]=dist_help(generator)+1;
  x_0T[i]=discsamp(p)+1;
  //Return results
  return List::create(Named("x") = x_0T);
}

// [[Rcpp::export]]
List Interval_sampling_accumulation_cpp(arma::ivec x_0T,
                                        arma::ivec y_0Tprime,
                                        float Inter,
                                        arma::mat Q,
                                        arma::vec Lambda,
                                        int pp=12,
                                        bool messages=true) {
  //initialize random number generators
  // boost::random::mt19937 generator( time(0) ); 
  // boost::random::uniform_real_distribution< >  uniform_distribution(0.0,1.0);
  // std::random_device rd;
  // std::default_random_engine generator(rd());
  std::uniform_real_distribution<double> uniform_distribution(0.0,1.0);
  Rcpp::Rcout << "\r Interval sampling:Initialisation                    ";  
  
  //length of input
  int t = x_0T.n_elem;
  //dim of state space
  int M= Q.n_rows;
  //transition matrix
  arma::mat D=arma::zeros(M,M);
  for(int i=1;i<=M;i++){
    D.row(i-1)=Q.row(i-1)*Q(i-1,i-1);
    D(i-1,i-1)=0;
  }
  int c_max = y_0Tprime.max(); //maximum observation
  
  if(c_max<5){
    c_max=5;
  }
  //Accumulation interval transition matrix
  arma::mat P= Fast_Matrix_exponential(Q,Lambda,Inter,(c_max+1)*M,pp,messages);
  
  // int c_max2 = 4;
  int M2 = (c_max+2)*M;
  arma::vec dd= Q.diag()-Lambda;
  float mu = -dd.min();
  int M22=2;
  arma::mat A_mat=arma::eye<arma::mat>(M,M)+(Q-arma::diagmat(Lambda))/mu;
  arma::mat B_mat=arma::diagmat(Lambda)/mu;
  arma::mat R = arma::eye<arma::mat>(M22*M,M22*M);
  R.submat(0,0,M22*M-1,M22*M-1)+=  arma::join_cols(arma::join_rows(Q-arma::diagmat(Lambda),
                                               arma::diagmat(Lambda)),
                                               arma::join_rows(arma::zeros(M,M),
                                                               Q-arma::diagmat(Lambda)))/mu;
  
  int N2 = 2*c_max; //number of already conducted matrix calculations on R
  //storage for all R^n in efficient form
  arma::mat R_n = arma::eye<arma::mat>(M,M2*(N2+1)); 
  // Calculate R^n efficiently using the sparse block form of R
  for(int i = 1; i <= N2; i++) {
    R_n.submat(0, 0+i*M2 , M-1, M-1+ i*M2 )=
      R_n.submat(0, 0+(i-1)*M2 , M-1, M-1+ (i-1)*M2 )*A_mat;
    for(int j = 1; j <= M2/M-1; j++) {
      R_n.submat(0, i*M2+j*M , M-1, M-1+j*M+i*M2 )=
        R_n.submat(0, j*M+(i-1)*M2 , M-1, M-1+j*M+ (i-1)*M2 )*A_mat+
        R_n.submat(0, (j-1)*M+(i-1)*M2 , M-1, M-1+(j-1)*M+ (i-1)*M2 )*B_mat;
    }
  }
  
  arma::ivec x_0Tprime(10*t);
  x_0Tprime.fill(0); //return vector of sampled states (with initialized length, might be adjusted later)
  arma::vec  t_0Tprime=arma::zeros(10*t); //return vector of sampled transition times
  int pathcounter1=0; //counter of current position in return vector
  int pathcounter2=10*t; //counter of length of return vector
  int iiii=0; //counter of number of state transitions into each state during each interval
  arma::vec Count = arma::zeros(M); //count of generated events during each state
  int progress=0;
  for(int i = 0; i <= (t-2); i++) {
    if(messages==true){
      if(i>=progress*(t/20.0)){
        Rcpp::Rcout << "\r Interval sampling:"<<(progress)*5<<"%                    ";
        progress++;
      }
    }
    // rejection sampling for easy transitions (no events and same begin and end state) since faster, but unreliable
    if((x_0T[i]==x_0T[i+1])&(y_0Tprime[i]==0)){
      bool endstatecheck=FALSE; //bool variable to check if sampled path has right end state
      arma::ivec state(100);//overlong state and time storage vectors to avoid extensions
      state.fill(0);
      arma::vec time=arma::zeros(100);
      //initialize
      state[0]=x_0T[i]; 
      iiii=1; 
      //variable of transition time
      float tau;
      //rate variable
      float c;
      int limitcounter=1; //iteration limiter
      //while the trajectory does not end right, sample again
      while((endstatecheck==FALSE)&(limitcounter<100)){
        limitcounter++;
        //set c to be equal to current rate
        c=-Q(x_0T[i]-1,x_0T[i]-1)+Lambda(x_0T[i]-1);
        //  G(x_0T[i]-1,x_0T[i]-1);
        // sample an exponentially distr. transition time with rate c
        tau=-log(uniform_distribution(generator))/c;
        iiii=1;
        //while transition is inside the acc. interval, sample another step
        while(tau<Inter){
          //sample state after transition
          arma::vec pp=(D.row(x_0T[i]-1)).t();
          // boost::random::discrete_distribution<> dist(pp);
          //put sampled state in storage vector
          // state[iiii]=dist(generator)+1;
          state[iiii]=discsamp(pp)+1;
          time[iiii]=tau;
          c=-Q(state[iiii]-1,state[iiii]-1)+Lambda(state[iiii]-1);
          // G(state[iiii]-1,state[iiii]-1);
          //sample another time
          tau+=-log(uniform_distribution(generator))/c;
          iiii++;
        }
        //set checker variable to true if end state is sampled right
        if(state[iiii-1]==x_0T[i]){
          endstatecheck=true;
        }
      }
      //put sampled state in storage vector
      state[iiii]=x_0T[i+1];
      time[iiii]=Inter;
      iiii+=1;
      //check if return vector has to be extended
      //if not, add sampled path to return vector
      if(pathcounter1+iiii<pathcounter2){
        x_0Tprime.rows(pathcounter1,pathcounter1+iiii-1)=state.rows(0,iiii-1);
        t_0Tprime.rows(pathcounter1,pathcounter1+iiii-1)=time.rows(0,iiii-1)+i*Inter;
      }else{
        //if yes, extend return vector
        pathcounter2+=10*t;
        x_0Tprime.resize(pathcounter2);
        t_0Tprime.resize(pathcounter2);
        x_0Tprime.rows(pathcounter1,pathcounter1+iiii-1)=state.rows(0,iiii-1);
        t_0Tprime.rows(pathcounter1,pathcounter1+iiii-1)=time.rows(0,iiii-1)+i*Inter;
      }
      //add transitions into each state to counter variable
      arma::vec x_0Tprimee = arma::conv_to<arma::vec>::from(x_0Tprime);
      arma::vec temp=arma::diff(arma::ceil(x_0Tprimee.rows(pathcounter1,pathcounter1+iiii-1)/M));
      arma::vec temp2=x_0Tprimee.rows(pathcounter1,pathcounter1+iiii-1)-
        arma::floor((x_0Tprimee.rows(pathcounter1,pathcounter1+iiii-1)-1)/M)*M;
      for(int kk=0; kk<M;kk++){
        for(int mm=0;mm<iiii-1;mm++){
          if(temp2[mm]==kk+1){
            Count[kk]+= temp[mm];
          }
        }
      }
      pathcounter1+=iiii;
    }
    // modified rejection sampling for more complex paths
    else{
      int N=0; //number of virtual state transitions (must be at least the number of events)
      float P_ab=P(x_0T[i]-1,x_0T[i+1]+M*y_0Tprime[i]-1); //total probability of begin and end state path
      float lll=mu*Inter; //auxilliary process rate
      // boost::math::poisson_distribution< >  p_lll(lll);
      //initialize cumulated probability from formula in Nihms paper 
      float prob= Poiss(lll,N)*R_n(x_0T[i]-1,x_0T[i+1]-1+M*y_0Tprime[i]+N*M2)/P_ab; 
      // float prob= pdf(p_lll,N)*R_n(x_0T[i]-1,x_0T[i+1]-1+M*y_0Tprime[i]+N*M2)/P_ab; 
      float u = uniform_distribution(generator);
      int jjj=1; //loop limiter
      float diffcount=0;
      //test if cum. prob. is greater than u, if not increase N by one and calculate cum. prob.
      int maxit=3000;
      while(prob<u&&jjj<maxit){
        
        jjj++;
        //test if R^n storage matrix has to be extended
        if(N<N2){
          //if not increase N by one and calculate cum. prob.
          N++;
          // diffcount=pdf(p_lll,N)/(P_ab)*R_n(x_0T[i]-1,x_0T[i+1]-1+M*y_0Tprime[i]+N*M2);
          diffcount=Poiss(lll,N)/(P_ab)*R_n(x_0T[i]-1,x_0T[i+1]-1+M*y_0Tprime[i]+N*M2);
          prob+= diffcount;
          //if yes, extend and add more R^n calculations to storage matrix
        }else{
          int N_temp=N2;
          N2=N2*4;
          R_n.resize(M,M2*(N2+1));
          for(int pp = N_temp+1; pp <= N2; pp++) {
            R_n.submat(0, 0+pp*M2 , M-1, M-1+ pp*M2 )=
              R_n.submat(0, 0+(pp-1)*M2 , M-1, M-1+ (pp-1)*M2 )*A_mat;
            for(int j = 1; j <= M2/M-1; j++) {
              R_n.submat(0, pp*M2+j*M , M-1, M-1+j*M+pp*M2 )=
                R_n.submat(0, j*M+(pp-1)*M2 , M-1, M-1+j*M+ (pp-1)*M2 )*A_mat+
                R_n.submat(0, (j-1)*M+(pp-1)*M2 , M-1, M-1+(j-1)*M+ (pp-1)*M2 )*B_mat;
            }
          }
          N++;
          // diffcount=pdf(p_lll,N)/(P_ab)*R_n(x_0T[i]-1,x_0T[i+1]-1+M*y_0Tprime[i]+N*M2);
          diffcount=Poiss(lll,N)/(P_ab)*R_n(x_0T[i]-1,x_0T[i+1]-1+M*y_0Tprime[i]+N*M2);
          prob+= diffcount;
        }
        if((prob>0)&(diffcount<0.1e-10)){
          N=0;
          P_ab=prob*P_ab;
          // prob= pdf(p_lll,N)*R_n(x_0T[i]-1,x_0T[i+1]-1+M*y_0Tprime[i]+N*M2)/P_ab;
          prob= Poiss(lll,N)*R_n(x_0T[i]-1,x_0T[i+1]-1+M*y_0Tprime[i]+N*M2)/P_ab;
          diffcount=0;
          maxit=6000;
        }
        if(maxit==5990){
          Rcpp::Rcout<<"\n something weird happend (@_@)\n";
        }
        
      }
      //if N is equal to zero, path is trivial
      if(N==0){
        arma::vec time;
        time<<0<<Inter;
        arma::ivec state;
        state<<x_0T[i]<<x_0T[i+1];
        
        iiii=state.n_elem;
        //check if return vector has to be extended (like before)
        if(pathcounter1+iiii<pathcounter2){
          x_0Tprime.rows(pathcounter1,pathcounter1+iiii-1)=state.rows(0,iiii-1);
          t_0Tprime.rows(pathcounter1,pathcounter1+iiii-1)=time.rows(0,iiii-1)+i*Inter;
        }else{
          pathcounter2+=10*t;
          x_0Tprime.resize(pathcounter2);
          t_0Tprime.resize(pathcounter2);
          x_0Tprime.rows(pathcounter1,pathcounter1+iiii-1)=state.rows(0,iiii-1);
          t_0Tprime.rows(pathcounter1,pathcounter1+iiii-1)=time.rows(0,iiii-1)+i*Inter;
        }
        //add transitions into each state to counter variable
        arma::vec x_0Tprimee = arma::conv_to<arma::vec>::from(x_0Tprime);
        arma::vec temp=arma::diff(arma::ceil(x_0Tprimee.rows(pathcounter1,pathcounter1+iiii-1)/M));
        arma::vec temp2=x_0Tprimee.rows(pathcounter1,pathcounter1+iiii-1)-
          arma::floor((x_0Tprimee.rows(pathcounter1,pathcounter1+iiii-1)-1)/M)*M;
        for(int kk=0; kk<M;kk++){
          for(int mm=0;mm<iiii-1;mm++){
            if(temp2[mm]==kk+1){
              Count[kk]+= temp[mm];
            }
          }
        }
        pathcounter1+=iiii;
      }
      //if N is equal to 1, path is trivial too
      else if(N==1){
        float uu = uniform_distribution(generator)*Inter;
        arma::vec time;
        time<<0<<uu<<Inter;
        arma::ivec state;
        state<<x_0T[i]<<x_0T[i+1]+M*y_0Tprime[i]<<x_0T[i+1]+M*y_0Tprime[i];
        iiii=state.n_elem;
        //check if return vector has to be extended (like before)
        if(pathcounter1+iiii<pathcounter2){
          x_0Tprime.rows(pathcounter1,pathcounter1+iiii-1)=state.rows(0,iiii-1);
          t_0Tprime.rows(pathcounter1,pathcounter1+iiii-1)=time.rows(0,iiii-1)+i*Inter;
        }else{
          pathcounter2+=10*t;
          x_0Tprime.resize(pathcounter2);
          t_0Tprime.resize(pathcounter2);
          x_0Tprime.rows(pathcounter1,pathcounter1+iiii-1)=state.rows(0,iiii-1);
          t_0Tprime.rows(pathcounter1,pathcounter1+iiii-1)=time.rows(0,iiii-1)+i*Inter;
        }
        //add transitions into each state to counter variable
        arma::vec x_0Tprimee = arma::conv_to<arma::vec>::from(x_0Tprime);
        arma::vec temp=arma::diff(arma::ceil(x_0Tprimee.rows(pathcounter1,pathcounter1+iiii-1)/M));
        arma::vec temp2=x_0Tprimee.rows(pathcounter1,pathcounter1+iiii-1)-
          arma::floor((x_0Tprimee.rows(pathcounter1,pathcounter1+iiii-1)-1)/M)*M;
        for(int kk=0; kk<M;kk++){
          for(int mm=0;mm<iiii-1;mm++){
            if(temp2[mm]==kk+1){
              Count[kk]+= temp[mm];
            }
          }
        }
        pathcounter1+=iiii;
      }
      //if path is not trivial, sample the path
      else{
        //initialize the path storage vector with the right length and begin and end states/times
        arma::vec time=arma::zeros(N+2);
        time[N+1]=Inter;
        arma::ivec state(N+2);
        state.fill(1);
        // =arma::ones(N+2);
        state[0]=x_0T[i];
        state[N+1]=x_0T[i+1]+M*y_0Tprime[i];
        //for each N, sample a virtual state transition with formula from Nihms (and using the efficient R^n storage)
        for(int ll =0; ll<N;ll++){
          
          time[ll+1]=uniform_distribution(generator)*Inter;
          float uuu=uniform_distribution(generator);
          // Rcpp::Rcout<<"\n uuu:"<<uuu<<", state[ll]="<<state[ll]<<", x_0T[i]="<<
          //   x_0T[i]<<", x_0T[i+1]="<<x_0T[i+1]<<", y_0Tprime[i]="<<y_0Tprime[i]<<"ll"<<ll<<"N"<<N;
          arma::vec p=arma::zeros(M2); //probability vector
          for(int kk=0;kk<M2-1;kk++){
            
            if((-(kk/M)*M+x_0T[i+1]+M*y_0Tprime[i]+(N-(ll+1))*M2-1>-1)){
              
              float R_insert=0;
              // int tempstatell;
              // int kkk;
              int Mbase;
              if(state[ll]-1<=kk){
                Mbase = std::floor((state[ll]-1)/M);
              }else{
                Mbase = std::floor((kk)/M);
              }
              
              if(!((state[ll]-1-Mbase*M>=2*M)|(kk-Mbase*M>=2*M))){
                R_insert=R(state[ll]-1-Mbase*M,kk-Mbase*M);
              }
              
              p[kk]=R_insert*R_n(kk%M,-(kk/M)*M+x_0T[i+1]+M*y_0Tprime[i]+(N-(ll+1))*M2-1);
              
              // Rcpp::Rcout<<" p["<<kk<<"]"<<p[kk];
              // Rcpp::Rcout<<" R_n["<<kk%M<<","<<-(kk/M)*M+x_0T[i+1]+M*y_0Tprime[i]+(N-(ll+1))*M2-1<<"]"
              // <<R_n(kk%M,-(kk/M)*M+x_0T[i+1]+M*y_0Tprime[i]+(N-(ll+1))*M2-1);
              // Rcpp::Rcout<<" R_insert"<<R_insert;
            }else{
              p[kk]=0;
            }
          }
          float sss=sum(p);
          // int control=0;
          for(int kk=0;kk<M2-1;kk++){
            p[kk+1]+=p[kk];
            if(p[kk]/sss<uuu){
              state[ll+1]=kk+2;
              // control=1;
            }
          }
          
          // //////////////////////////////////////////////////////////
          // if((state[ll]>M)&(state[ll+1]<M)){
          //   Rcpp::Rcout << "state[ll]:"<< state[ll];
          //   Rcpp::Rcout << " state[ll+1]:"<< state[ll+1];
          //   Rcpp::Rcout << " state[N+1]:"<< state[N+1];
          //   Rcpp::Rcout << " sss:"<< p[0];
          // }
          // //////////////////////////////////////////////////////////
          // 
          // Rcpp::Rcout<<" state["<<ll+1<<"]"<<state[ll+1]<<"\n";
        }
        //check if return vector has to be extended (like before)
        iiii=state.n_elem;
        while(pathcounter1+iiii>=pathcounter2){
          pathcounter2+=2*pathcounter2;
          x_0Tprime.resize(pathcounter2);
          t_0Tprime.resize(pathcounter2);
        }
        x_0Tprime.rows(pathcounter1,pathcounter1+iiii-1)=state.rows(0,iiii-1);
        t_0Tprime.rows(pathcounter1,pathcounter1+iiii-1)=arma::sort(time.rows(0,iiii-1))+i*Inter;
        //add transitions into each state to counter variable
        arma::vec x_0Tprimee = arma::conv_to<arma::vec>::from(x_0Tprime);
        arma::vec temp=arma::diff(arma::ceil(x_0Tprimee.rows(pathcounter1,pathcounter1+iiii-1)/M));
        arma::vec temp2=x_0Tprimee.rows(pathcounter1,pathcounter1+iiii-1)-
          arma::floor((x_0Tprimee.rows(pathcounter1,pathcounter1+iiii-1)-1)/M)*M;
        for(int kk=0; kk<M;kk++){
          for(int mm=0;mm<iiii-1;mm++){
            if(temp2[mm]==kk+1){
              
              // //////////////////////////////////////////////////////////
              // if(temp[mm]<0){
              //   Rcpp::Rcout << "N:"<< N;
              //   Rcpp::Rcout << " iiii:"<< iiii;
              //   Rcpp::Rcout << " state[mm]:"<< state[mm];
              //   Rcpp::Rcout << " state[mm+1]:"<< state[mm+1];
              //   Rcpp::Rcout << " state[mm+2]:"<< state[mm+2];
              //   Rcpp::Rcout << " state[N+1]:"<< state[N+1];
              //   Rcpp::Rcout << " y_0T[i]:"<< y_0Tprime[i];
              //   Rcpp::Rcout << "Phase3 temp["<<mm<<"]:"<<temp[mm]<<"\n";
              // }
              // //////////////////////////////////////////////////////////
              
              Count[kk]+= temp[mm];
            }
          }
        }
        pathcounter1+=iiii;
      }
    }
  }
  if(messages==true){
    Rcpp::Rcout << "\r Interval sampling: Export results                    ";  
  }
  
  //cut the excess parts of the return vectors
  arma::vec x_0Tprimee = arma::conv_to<arma::vec>::from(x_0Tprime);
  arma::vec x=x_0Tprimee.rows(0,pathcounter1-1)-
    arma::floor((x_0Tprimee.rows(0,pathcounter1-1)-1)/M)*M;
  arma::vec time=t_0Tprime.rows(0,pathcounter1-1);
  return List::create(Named("x") = x,
                      Named("t") = time,
                      Named("Count") = Count);
  
}

// [[Rcpp::export]]
List Z_Sampling_cpp(arma::ivec y_0Tprime,
                    arma::ivec path_x,
                    arma::vec path_t,
                    float Inter,
                    float Lambda_Z,
                    arma::vec Lambda,
                    bool messages=true) {
  //initialize random number generator
  // boost::random::mt19937 generator2( time(0) );
  std::random_device rd_Z;
  std::default_random_engine generator_Z(rd_Z());
  std::uniform_real_distribution<float> uniform_distribution_Z(0.0,1.0);
  int t = y_0Tprime.n_elem; //number of acc. intervals
  //Initialize the index vectors that will mark the beginning and ending of acc. intervals in path_x and path_t
  arma::ivec start_index(t);
  arma::ivec end_index(t);
  start_index(0)=1;
  end_index(t-1)=path_t.n_elem;
  
  //rounding precision set to four digits
  float Ndigits=100000;
  //initialize return vector
  arma::ivec Z(t);
  Z.fill(1);  
  int dummy_var=2; //dummy variable that marks the current position in path_t
  //iterate through all acc. intervals (-1 because the ultimate starting and ending indices are known)
  int progress=0;
  for(int i = 1; i <= t-1; i++) { 
    if(messages==true){
      if(i>=progress*(t/20.0)){
        Rcpp::Rcout << "\r Latent variable sampling, stage 1:"<<(progress)*5<<"%                    ";
        progress++;
      }
    }
    //increase dummy variable while the next end/start time is not reached yet while iterating through path_t
    while(round(path_t(dummy_var-1)*Ndigits)/Ndigits<i){ 
      dummy_var++;
    }
    //when the next end/start time is reached, set the number in the index vectors to be the position in path_t
    start_index(i)=dummy_var;
    end_index(i-1)=dummy_var;
  }
  progress=0;
  //iterate through all intervals again for sampling Z
  for(int i = 1; i <= t; i++) {
    if(messages==true){
      if(i>=progress*(t/20.0)){
        Rcpp::Rcout << "\r Latent variable sampling, stage 2:"<<(progress)*5<<"%                    ";
        progress++;
      }
    }
    //if y=0, Z has to be equal to 0 too
    if(y_0Tprime[i-1]==0){
      Z[i-1]=0;
      //similarly, if y=1 --> Z=1
    }else if(y_0Tprime[i-1]==1){
      Z[i-1]=1;
    }else{
      //calculate lambda_I_k as in formula 4.6 of my master thesis:
      float Lambda_I_k = 0; //initialize to be zero
      //iterate through all entries in path_x during the current acc. interval
      for(int kk=start_index[i-1];kk<end_index[i-1];kk++){
        //add the corresponding value to lambda_I_k
        Lambda_I_k+=Lambda[path_x[kk-1]-1]*(path_t[kk]-path_t[kk-1]);
      }
      //initialize the probabilities and the poisson distribution with rate lambda_I_k
      arma::vec Probs=arma::zeros(y_0Tprime[i-1]);
      // boost::math::poisson_distribution< >  p1(Lambda_I_k+0.0000000000000001);
      
      //iterate through all possible values of Z to calculate the corresponding probabilities 
      //according to formula 4.7 of my master thesis
      for(int zz=1;zz<=y_0Tprime[i-1];zz++){
        // boost::math::poisson_distribution< >  p2(Lambda_Z*zz+0.0000000000000001);
        // Probs[zz-1]=1e100*pdf(p1,zz)*pdf(p2,y_0Tprime[i-1]-zz);
        Probs[zz-1]=Poiss(Lambda_I_k+0.0000000000000001,zz)*
        Poiss(Lambda_Z*zz+0.0000000000000001,y_0Tprime[i-1]-zz);
      }
      //sample Z from a discrete probability distribution with the corresponding weights
      arma::vec Probs2 = Probs;
      float uuu=uniform_distribution_Z(generator_Z);
      float sss=sum(Probs);
      float Ndigits=1000000;
      if(sss==0){
        // Rcpp::Rcout<<"\n \n wtf happened?!?!? \n \n";
        Probs=arma::ones(Probs.n_elem);
      }else{
        Probs=Probs/sss;
        Probs=floor(Probs*Ndigits)/Ndigits;
      }
      Z[i-1]=1;
      for(unsigned int kk=0;kk<Probs.n_elem-1;kk++){
        Probs[kk+1]+=Probs[kk];
        if(Probs[kk]<uuu){
          Z[i-1]=kk+2;
        }
      }
      // boost::random::discrete_distribution<> dist(Probs2);
      // Z[i-1]=dist(generator2)+1;
      // Z[i-1]=discsamp(Probs)+1;
      // arma::vec test=arma::ones(3);
      // int test2=discsamp(Probs)+1;
    }
  } 
  // Rcpp::Rcout<<"\n";
  return List::create(Named("z") = Z);
}

