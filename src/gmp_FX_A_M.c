#include <assert.h>
#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void bincoef_array(mpz_t b[], int K) {
  // Length of b:  K
  // [bincoef(K, 1), bincoef(K, 2), ....., bincoef(K, K)]
  for (int i = 0; i < K; i++) {
    mpz_init(b[i]);
    mpz_bin_uiui(b[i], K, i + 1);
  }
}

void factK_Stirling(mpz_t S, int K, int l, mpz_t bincoef_a[]) {
  // Implements K! Stirling(K, l), K <= l of the second kind
  // from (1.6.7), H. S. Wilf, Generating functionology, 2nd ed., Elsevier, 2013
  // bincoef_a contains [bincoef(K, 1), bincoef(K, 2), ....., bincoef(K, K)]
  // modifies S.

  mpz_t f_S;
  mpz_t product, powil;

  mpz_init_set_ui(f_S, 0);
  mpz_init(product);
  mpz_init(powil);

  for (int i = 1; i <= K; ++i) {
    mpz_ui_pow_ui(powil, i, l);
    mpz_mul(product, powil, bincoef_a[i - 1]);
    mpz_mul_si(product, product, pow(-1, K - i));
    mpz_add(f_S, f_S, product);
  }
  mpz_set(S, f_S);

  mpz_clear(product);
  mpz_clear(powil);
  mpz_clear(f_S);
}

void factA_Stirling_array(mpz_t f_S_a[], int A, int m) {
  // compute A! Stirling(A, l), for A <= l <= m. Modifies f_S_a.
  mpz_t bincoef_a[A];

  bincoef_array(bincoef_a, A);
  for (int j = A; j <= m; j++) {
    mpz_init_set_ui(f_S_a[j - A], 0);
    factK_Stirling(f_S_a[j - A], A, j, bincoef_a);
  }

  for (int i = 0; i < A; i++) {
    mpz_clear(bincoef_a[i]);
  }
}

void F_array(int A, int M, int n, double pg[], int len_pg, double F_threshold) {
  // Computes P(X_{A, M}) = n and P(X_{A, M} <= t), for A + 1 <= t <= n, writing
  // the results as t, P(X_{A, M} = t, P(X_{A, M} <= t) to a csv file. It loops
  // over the values of pg, the probability of extracting a red ball. For each
  // pg, a different csv file is generated. F_threshold allows to specify an
  // early stop, if the cdf(t) > F_threshold
  mpf_t pgg, prg;

  mpz_t f_S_a[n - A];
  factA_Stirling_array(f_S_a, A, n - 1);

  mpf_t prob, sum_prob;
  mpf_init(prob);
  mpf_init(sum_prob);

  mpf_t one_minus_pg_pow, pr_pow;
  mpf_init(one_minus_pg_pow);
  mpf_init(pr_pow);
  mpf_t prod;
  mpf_init(prod);

  mpz_t bin_coef;
  mpf_t bin_coef_f;
  mpz_init(bin_coef);
  mpf_init(bin_coef_f);

  mpf_t s_f;
  mpf_init(s_f);

  char file_name[50];

  int t;

  for (int i = 0; i < len_pg; i++) {
    sprintf(file_name, "../results/pmf_cdf_A_%d_M_%d_n_%d_%d.csv", A, M, n,
            (int)(pg[i] * 1000));
    FILE *results_file = fopen(file_name, "w+");
    fprintf(results_file, "t,p,F\n");

    mpf_init_set_d(pgg, pg[i]);
    mpf_init_set_d(prg, pg[i] / M);

    mpf_set_ui(sum_prob, 0);
    t = A;
    while (t <= n && mpf_cmp_d(sum_prob, F_threshold) < 0) {
      ++t;
      mpf_set_ui(prob, 0);
      for (int l = A; l < t; l++) {
        mpf_ui_sub(one_minus_pg_pow, 1, pgg);
        mpf_pow_ui(one_minus_pg_pow, one_minus_pg_pow, t - l - 1);
        mpf_pow_ui(pr_pow, prg, l + 1);
        mpf_mul(prod, one_minus_pg_pow, pr_pow);
        mpf_mul_ui(prod, prod, A + 1);
        mpz_bin_uiui(bin_coef, t - 1, l);
        mpf_set_z(bin_coef_f, bin_coef);
        mpf_mul(prod, prod, bin_coef_f);
        mpf_set_z(s_f, f_S_a[l - A]);
        mpf_mul(prod, prod, s_f);
        mpf_add(prob, prob, prod);
      }
      mpz_bin_uiui(bin_coef, M, A + 1);
      mpf_set_z(bin_coef_f, bin_coef);
      mpf_mul(prob, prob, bin_coef_f);
      mpf_add(sum_prob, sum_prob, prob);
      fprintf(results_file, "%d,", t);
      mpf_out_str(results_file, 10, 10, prob);
      fprintf(results_file, ",");
      mpf_out_str(results_file, 10, 10, sum_prob);
      fprintf(results_file, "\n");
    }
    printf("pg: %f, F(%d) computed, value:", pg[i], t);
    mpf_out_str(stdout, 10, 10, sum_prob);
    printf("\n");

    fclose(results_file);
  }
  mpf_clear(pgg);
  mpf_clear(prg);
  for (int k = 0; k < n - A; k++) {
    mpz_clear(f_S_a[k]);
  }
  mpf_clear(prob);
  mpf_clear(sum_prob);
  mpf_clear(prod);
  mpf_clear(one_minus_pg_pow);
  mpf_clear(pr_pow);
  mpz_clear(bin_coef);
  mpf_clear(bin_coef_f);
  mpf_clear(s_f);
}

int main(int argc, char *argv[]) {
  if (argc <= 3) {
    printf("Usage: %s <number A> <number M> <number n> \n", argv[0]);
    return 2;
  }
  int A;
  int M;
  int n;

  A = atoi(argv[1]);
  assert(A >= 0);
  M = atoi(argv[2]);
  assert(M >= A);
  n = atoi(argv[3]);
  assert(n >= M);

  double pg[3] = {0.7, 0.95, 0.999};

  int len_pg = sizeof(pg) / sizeof(pg[0]);

  F_array(A, M, n, pg, len_pg, 0.96);

  printf("Done!");

  return 1;
}



