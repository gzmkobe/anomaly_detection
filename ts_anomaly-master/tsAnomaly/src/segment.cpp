#include <Rcpp.h>
#include <iterator>
#include <vector>
#include <tuple>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

namespace std{
    namespace
    {
        template <class T>
        inline void hash_combine(std::size_t& seed, T const& v)
        {
            seed ^= hash<T>()(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        }
        template <class Tuple, size_t Index = std::tuple_size<Tuple>::value - 1>
        struct HashValueImpl
        {
          static void apply(size_t& seed, Tuple const& tuple)
          {
            HashValueImpl<Tuple, Index-1>::apply(seed, tuple);
            hash_combine(seed, get<Index>(tuple));
          }
        };
        template <class Tuple>
        struct HashValueImpl<Tuple,0>
        {
          static void apply(size_t& seed, Tuple const& tuple)
          {
            hash_combine(seed, get<0>(tuple));
          }
        };
    }

    template <typename ... TT>
    struct hash<std::tuple<TT...>> 
    {
        size_t
        operator()(std::tuple<TT...> const& tt) const
        {                                              
            size_t seed = 0;                             
            HashValueImpl<std::tuple<TT...> >::apply(seed, tt);    
            return seed;                                 
        }                                              
    };
}

template <typename ReturnType, typename... Args>
std::function<ReturnType (Args...)> memoize(std::function<ReturnType (Args...)> f) {
  //std::map<std::tuple<Args...>, ReturnType> mem;
  std::unordered_map<std::tuple<Args...>, ReturnType> mem;
  return ([=](Args... args) mutable {
    std::tuple<Args...> key(args...);
    if (mem.find(key) == mem.end())
      mem[key] = f(args...);
    return mem[key];
  });
}

// std::function<double (int, int, double*)> MU;
// double mu(int s, int t, double *X) {
//   if (s == t) {
//     return *(X + t-1);
//   }
//   else {
//     double temp = (t - 1) - s + 1;
//     temp *= MU(s, t-1, X) + *(X + t-1);
//     temp /= ((t-1)-s + 2);
//     return temp;
//   }
// }
//
// std::function<double (int, int, double*)> D;
// double d(int s,int t, double *X) {
//   if (s == t){
//     return 0.;
//   }
//   else {
//     double temp = (t - 1) - s + 1;
//     temp *= pow(MU(s, t-1, X) - MU(s, t, X), 2);
//     temp += pow( X[t-1] - MU(s, t, X), 2);
//     return D(s, t-1, X) + temp;
//   }
// }

std::function<double (int, int, double*)> MU;
double mu(int s, int t, double *X) {
  if (s == t) {
    return *(X+t-1);
  }
  else {
    NumericVector temp(t-s+1);
    for (int i = 0; i <= (t-s); i++){
      temp[i] = *(X+s+i);
    }
    return median(temp, 0);
  }
}

std::function<double (int, int, double*)> D;
double d(int s, int t, double *X) {
  if (s == t){
    return 0.;
  }
  else {
    double temp = (t - 1) - s + 1;
    temp *= abs( MU(s, t-1, X) - MU(s, t, X) );
    temp += abs( *(X+t-1) - MU(s, t, X) );
    return D(s, t-1, X) + temp;
  }
}

double schefes_score(NumericVector group_a, NumericVector group_b, double var){
  double num = mean(group_a) - mean(group_b);
  num = pow(num, 2.);
  double den = (1./group_a.size()) + (1./group_b.size());
  den *= var;
  return num/den;
}

double group_var(NumericVector x){
  return var(x)*(x.size());
}

bool is_mod_n(int X, int N){
  return (X - 1) % N == 0;
}

NumericVector every_n(NumericVector X, int N) {
  IntegerVector inter = seq_along(X);
  return X[mapply(inter, N, is_mod_n)];
}

//' Segment series into adjacent regimes
//' 
//' @param Y Numeric. Vector to be segmented
//' 
//' @param K Integer. Number of segments. Defaults to 5.
//'
//' @param alpha Numeric. Critical value for F-statistic. Defaults to 0.05. 
//'
//' @param exact Logical. Skip reduction of time series. Defaults to FALSE. 
//'
//' @return Vector of integers specifying segmentation points of the series.
//'
//' @export
//'
// [[Rcpp::export]]
NumericVector segment(NumericVector Y, int K=5, double alpha=0.05, bool exact=false){
  MU = memoize(std::function<double (int, int, double*)>(mu));
  D = memoize(std::function<double (int, int, double*)>(d));
  int pre_T = Y.size();
  NumericVector X;
  int N;
  int T;

  if (pre_T < 150 || exact == true) {
    N = 1;
    X = Y;
    T = X.size();
  } else {
    N = (pre_T/150) + 1;
    X = every_n(Y, N);
    T = X.size();
  }

  NumericMatrix cost(K + 1, T + 1);
  NumericMatrix Z(K + 1, T + 1);
  NumericMatrix That(K + 1, K + 1);
  
  //cost(1,0) = 0;
  for(int t = 1; t <= T; t++){
    cost(1 , t) = D(1, t, &X[0]);
    //Z(1,t) = 0;
  }
  for(int k=2; k <= K; k++){
    //cost(k,0) = 0
    for (int t =1; t <= T; t++){
      NumericVector vals(t-1);
      for (int s = 1; s <= t - 1; s++){
        vals[s - 1] = cost(k - 1, s) + D(s + 1, t, &X[0]);
        if(t > 1){
          Z(k, t) = which_min(vals);
          cost(k,t) = vals[Z(k, t)];
        }
      }
    }
  }
  
  for(int k=1; k <= K; k++){
    That(k,k) = T - 1;
    for(int kk=k; kk >= 2; kk--){
      That(k, kk - 1) = Z(k, That(k, kk) );
    }
  }

  That = N*That;
  for(int k=1; k <= K; k++){
    That(k,k) = pre_T - 1;
  }

  //std::cout << "\n";
  //std::cout << That;

  NumericVector group;
  for(int g=K; g >= 1; g--) {
    NumericVector subset =  That(g, _ );
    if(subset[1] == 0){
      continue;
    }

    double df1 = g - 1;
    double df2 = T - g;
    group = subset[Range(0, g)];

    double criticalF = R::qf(alpha, df1, df2, 0, 0);
    criticalF *= df1;

    List groups(group.size() - 1);
    for(int i = 0; i < group.size() - 1; i++) {
      int start;
      if( i != 0){
        start = group[i] + 1;
      } else {
        start = group[i];
      }
      NumericVector out = Y[Range(start, group[i+1])];
      groups[i] = out;
    }
    double var = sum(sapply(groups, group_var));
    var /= df2;

    IntegerVector test(g-1);
    for(int i = 0; i < (g - 1); i++){

      NumericVector group_a = groups[i];
      NumericVector group_b = groups[i+1];
      double score = schefes_score(group_a, group_b, var);

      if( score >= criticalF){
        test[i] = 1;
      } else {
        test[i] = 0;
      }

    }
    if (sum(test) == test.size()){
      break;
    }
  }
  return group;
}
