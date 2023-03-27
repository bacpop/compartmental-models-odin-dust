#include <dust/r/dust.hpp>

template <typename real_type, typename container>
__host__ __device__ real_type odin_sum1(const container x, size_t from, size_t to);
template <typename real_type, typename container>
__host__ __device__ real_type odin_sum2(const container x, int from_i, int to_i, int from_j, int to_j, int dim_x_1);
// [[dust::class(WrightFisher_nGenotypes_haploid)]]
// [[dust::param(gene_no, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(GeneFitness, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(Genotypes, has_default = FALSE, default_value = NULL, rank = 2, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(Pop_ini, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(pop_size, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(species_no, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(dt, has_default = TRUE, default_value = 1L, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
class WrightFisher_nGenotypes_haploid {
public:
  using real_type = double;
  using rng_state_type = dust::random::generator<real_type>;
  using data_type = dust::no_data;
  struct shared_type {
    int dim_FitnessMatrix;
    int dim_FitnessMatrix_1;
    int dim_FitnessMatrix_2;
    int dim_GeneFitness;
    int dim_GenotypeFitness;
    int dim_Genotypes;
    int dim_Genotypes_1;
    int dim_Genotypes_2;
    int dim_Pop;
    int dim_Pop_ini;
    int dim_Pop_ini_norm;
    int dim_probs;
    int dim_y;
    real_type dt;
    std::vector<real_type> FitnessMatrix;
    int gene_no;
    std::vector<real_type> GeneFitness;
    std::vector<real_type> GenotypeFitness;
    std::vector<real_type> Genotypes;
    std::vector<real_type> initial_Pop;
    real_type initial_time;
    std::vector<real_type> Pop_ini;
    std::vector<real_type> Pop_ini_norm;
    real_type pop_size;
    int species_no;
  };
  struct internal_type {
    std::vector<real_type> probs;
    std::vector<real_type> y;
  };
  WrightFisher_nGenotypes_haploid(const dust::pars_type<WrightFisher_nGenotypes_haploid>& pars) :
    shared(pars.shared), internal(pars.internal) {
  }
  size_t size() {
    return shared->dim_Pop + 1;
  }
  std::vector<real_type> initial(size_t step) {
    std::vector<real_type> state(shared->dim_Pop + 1);
    state[0] = shared->initial_time;
    std::copy(shared->initial_Pop.begin(), shared->initial_Pop.end(), state.begin() + 1);
    return state;
  }
  void update(size_t step, const real_type * state, rng_state_type& rng_state, real_type * state_next) {
    const real_type * Pop = state + 1;
    state_next[0] = (step + 1) * shared->dt;
    for (int i = 1; i <= shared->dim_probs; ++i) {
      internal.probs[i - 1] = (1 + shared->GenotypeFitness[i - 1]) * Pop[i - 1] / (real_type) shared->pop_size;
    }
    {
       int i = 1;
       internal.y[i - 1] = (internal.probs[0] / (real_type) odin_sum1<real_type>(internal.probs.data(), 0, shared->species_no) < 1 ? dust::random::binomial<real_type>(rng_state, shared->pop_size, internal.probs[0] / (real_type) odin_sum1<real_type>(internal.probs.data(), 0, shared->species_no)) : shared->pop_size);
    }
    for (int i = 2; i <= (shared->species_no - 1); ++i) {
      internal.y[i - 1] = (internal.probs[i - 1] / (real_type) odin_sum1<real_type>(internal.probs.data(), i - 1, shared->species_no) < 1 ? dust::random::binomial<real_type>(rng_state, shared->pop_size - odin_sum1<real_type>(internal.y.data(), 0, i - 1), internal.probs[i - 1] / (real_type) odin_sum1<real_type>(internal.probs.data(), i - 1, shared->species_no)) : shared->pop_size - odin_sum1<real_type>(internal.y.data(), 0, i - 1));
    }
    {
       int i = shared->species_no;
       internal.y[i - 1] = shared->pop_size - odin_sum1<real_type>(internal.y.data(), 0, shared->species_no - 1);
    }
    for (int i = 1; i <= shared->dim_Pop; ++i) {
      state_next[1 + i - 1] = internal.y[i - 1];
    }
  }
private:
  std::shared_ptr<const shared_type> shared;
  internal_type internal;
};
template <typename real_type, typename container>
__host__ __device__ real_type odin_sum2(const container x, int from_i, int to_i, int from_j, int to_j, int dim_x_1) {
  real_type tot = 0.0;
  for (int j = from_j; j < to_j; ++j) {
    int jj = j * dim_x_1;
    for (int i = from_i; i < to_i; ++i) {
      tot += x[i + jj];
    }
  }
  return tot;
}
#include <array>
#include <cpp11/R.hpp>
#include <cpp11/sexp.hpp>
#include <cpp11/doubles.hpp>
#include <cpp11/integers.hpp>
#include <cpp11/list.hpp>
#include <cpp11/strings.hpp>
#include <memory>
#include <vector>

template <typename T>
inline bool is_na(T x);

template <>
inline bool is_na(int x) {
  return x == NA_INTEGER;
}

template <>
inline bool is_na(double x) {
  return ISNA(x);
}

inline size_t object_length(cpp11::sexp x) {
  return ::Rf_xlength(x);
}

template <typename T>
void user_check_value(T value, const char *name, T min, T max) {
  if (is_na(value)) {
    cpp11::stop("'%s' must not be NA", name);
  }
  if (!is_na(min) && value < min) {
    cpp11::stop("Expected '%s' to be at least %g", name, (double) min);
  }
  if (!is_na(max) && value > max) {
    cpp11::stop("Expected '%s' to be at most %g", name, (double) max);
  }
}

template <typename T>
void user_check_array_value(const std::vector<T>& value, const char *name,
                            T min, T max) {
  for (auto& x : value) {
    user_check_value(x, name, min, max);
  }
}

inline size_t user_get_array_rank(cpp11::sexp x) {
  if (!::Rf_isArray(x)) {
    return 1;
  } else {
    cpp11::integers dim = cpp11::as_cpp<cpp11::integers>(x.attr("dim"));
    return dim.size();
  }
}

template <size_t N>
void user_check_array_rank(cpp11::sexp x, const char *name) {
  size_t rank = user_get_array_rank(x);
  if (rank != N) {
    if (N == 1) {
      cpp11::stop("Expected a vector for '%s'", name);
    } else if (N == 2) {
      cpp11::stop("Expected a matrix for '%s'", name);
    } else {
      cpp11::stop("Expected an array of rank %d for '%s'", N, name);
    }
  }
}

template <size_t N>
void user_check_array_dim(cpp11::sexp x, const char *name,
                          const std::array<int, N>& dim_expected) {
  cpp11::integers dim = cpp11::as_cpp<cpp11::integers>(x.attr("dim"));
  for (size_t i = 0; i < N; ++i) {
    if (dim[(int)i] != dim_expected[i]) {
      Rf_error("Incorrect size of dimension %d of '%s' (expected %d)",
               i + 1, name, dim_expected[i]);
    }
  }
}

template <>
inline void user_check_array_dim<1>(cpp11::sexp x, const char *name,
                                    const std::array<int, 1>& dim_expected) {
  if ((int)object_length(x) != dim_expected[0]) {
    cpp11::stop("Expected length %d value for '%s'", dim_expected[0], name);
  }
}

template <size_t N>
void user_set_array_dim(cpp11::sexp x, const char *name,
                        std::array<int, N>& dim) {
  cpp11::integers dim_given = cpp11::as_cpp<cpp11::integers>(x.attr("dim"));
  std::copy(dim_given.begin(), dim_given.end(), dim.begin());
}

template <>
inline void user_set_array_dim<1>(cpp11::sexp x, const char *name,
                                  std::array<int, 1>& dim) {
  dim[0] = object_length(x);
}

template <typename T>
T user_get_scalar(cpp11::list user, const char *name,
                  const T previous, T min, T max) {
  T ret = previous;
  cpp11::sexp x = user[name];
  if (x != R_NilValue) {
    if (object_length(x) != 1) {
      cpp11::stop("Expected a scalar numeric for '%s'", name);
    }
    // TODO: when we're getting out an integer this is a bit too relaxed
    if (TYPEOF(x) == REALSXP) {
      ret = cpp11::as_cpp<T>(x);
    } else if (TYPEOF(x) == INTSXP) {
      ret = cpp11::as_cpp<T>(x);
    } else {
      cpp11::stop("Expected a numeric value for %s", name);
    }
  }

  if (is_na(ret)) {
    cpp11::stop("Expected a value for '%s'", name);
  }
  user_check_value<T>(ret, name, min, max);
  return ret;
}

template <>
inline float user_get_scalar<float>(cpp11::list user, const char *name,
                                    const float previous, float min, float max) {
  double value = user_get_scalar<double>(user, name, previous, min, max);
  return static_cast<float>(value);
}

template <typename T>
std::vector<T> user_get_array_value(cpp11::sexp x, const char * name,
                                    T min, T max) {
  std::vector<T> ret = cpp11::as_cpp<std::vector<T>>(x);
  user_check_array_value<T>(ret, name, min, max);
  return ret;
}

template <typename T, size_t N>
std::vector<T> user_get_array_fixed(cpp11::list user, const char *name,
                                    const std::vector<T> previous,
                                    const std::array<int, N>& dim,
                                    T min, T max) {
  cpp11::sexp x = user[name];
  if (x == R_NilValue) {
    if (previous.size() == 0) {
      cpp11::stop("Expected a value for '%s'", name);
    }
    return previous;
  }

  user_check_array_rank<N>(x, name);
  user_check_array_dim<N>(x, name, dim);

  return user_get_array_value<T>(x, name, min, max);
}

template <typename T, size_t N>
std::vector<T> user_get_array_variable(cpp11::list user, const char *name,
                                       std::vector<T> previous,
                                       std::array<int, N>& dim,
                                       T min, T max) {
  cpp11::sexp x = user[name];
  if (x == R_NilValue) {
    if (previous.size() == 0) {
      cpp11::stop("Expected a value for '%s'", name);
    }
    return previous;
  }

  user_check_array_rank<N>(x, name);
  user_set_array_dim<N>(x, name, dim);

  return user_get_array_value<T>(x, name, min, max);
}

template <>
inline std::vector<float> user_get_array_value(cpp11::sexp x, const char * name,
                                               float min, float max) {
  // NOTE: possible under/overflow here for min/max because we've
  // downcast this.
  std::vector<double> value = user_get_array_value<double>(x, name, min, max);
  std::vector<float> ret(value.size());
  std::copy(value.begin(), value.end(), ret.begin());
  return ret;
}

// This is sum with inclusive "from", exclusive "to", following the
// same function in odin
template <typename real_type, typename container>
__host__ __device__
real_type odin_sum1(const container x, size_t from, size_t to) {
  real_type tot = 0.0;
  for (size_t i = from; i < to; ++i) {
    tot += x[i];
  }
  return tot;
}

inline cpp11::writable::integers integer_sequence(size_t from, size_t len) {
  cpp11::writable::integers ret(len);
  int* data = INTEGER(ret);
  for (size_t i = 0, j = from; i < len; ++i, ++j) {
    data[i] = j;
  }
  return ret;
}
namespace dust {
template<>
dust::pars_type<WrightFisher_nGenotypes_haploid> dust_pars<WrightFisher_nGenotypes_haploid>(cpp11::list user) {
  using real_type = typename WrightFisher_nGenotypes_haploid::real_type;
  auto shared = std::make_shared<WrightFisher_nGenotypes_haploid::shared_type>();
  WrightFisher_nGenotypes_haploid::internal_type internal;
  shared->initial_time = 0;
  shared->gene_no = NA_INTEGER;
  shared->pop_size = NA_REAL;
  shared->species_no = NA_INTEGER;
  shared->dt = 1;
  shared->dt = user_get_scalar<real_type>(user, "dt", shared->dt, NA_REAL, NA_REAL);
  shared->gene_no = user_get_scalar<int>(user, "gene_no", shared->gene_no, NA_INTEGER, NA_INTEGER);
  shared->pop_size = user_get_scalar<real_type>(user, "pop_size", shared->pop_size, NA_REAL, NA_REAL);
  shared->species_no = user_get_scalar<int>(user, "species_no", shared->species_no, NA_INTEGER, NA_INTEGER);
  shared->dim_FitnessMatrix_1 = shared->gene_no;
  shared->dim_FitnessMatrix_2 = shared->species_no;
  shared->dim_GeneFitness = shared->gene_no;
  shared->dim_GenotypeFitness = shared->species_no;
  shared->dim_Genotypes_1 = shared->gene_no;
  shared->dim_Genotypes_2 = shared->species_no;
  shared->dim_Pop = shared->species_no;
  shared->dim_Pop_ini = shared->species_no;
  shared->dim_Pop_ini_norm = shared->species_no;
  shared->dim_probs = shared->species_no;
  shared->dim_y = shared->species_no;
  shared->GenotypeFitness = std::vector<real_type>(shared->dim_GenotypeFitness);
  shared->initial_Pop = std::vector<real_type>(shared->dim_Pop);
  shared->Pop_ini_norm = std::vector<real_type>(shared->dim_Pop_ini_norm);
  internal.probs = std::vector<real_type>(shared->dim_probs);
  internal.y = std::vector<real_type>(shared->dim_y);
  shared->dim_FitnessMatrix = shared->dim_FitnessMatrix_1 * shared->dim_FitnessMatrix_2;
  shared->dim_Genotypes = shared->dim_Genotypes_1 * shared->dim_Genotypes_2;
  shared->GeneFitness = user_get_array_fixed<real_type, 1>(user, "GeneFitness", shared->GeneFitness, {shared->dim_GeneFitness}, NA_REAL, NA_REAL);
  shared->Pop_ini = user_get_array_fixed<real_type, 1>(user, "Pop_ini", shared->Pop_ini, {shared->dim_Pop_ini}, NA_REAL, NA_REAL);
  shared->FitnessMatrix = std::vector<real_type>(shared->dim_FitnessMatrix);
  shared->Genotypes = user_get_array_fixed<real_type, 2>(user, "Genotypes", shared->Genotypes, {shared->dim_Genotypes_1, shared->dim_Genotypes_2}, NA_REAL, NA_REAL);
  for (int i = 1; i <= shared->dim_Pop_ini_norm; ++i) {
    shared->Pop_ini_norm[i - 1] = shared->Pop_ini[i - 1] / (real_type) odin_sum1<real_type>(shared->Pop_ini.data(), 0, shared->species_no);
  }
  for (int i = 1; i <= shared->dim_FitnessMatrix_1; ++i) {
    for (int j = 1; j <= shared->dim_FitnessMatrix_2; ++j) {
      shared->FitnessMatrix[i - 1 + shared->dim_FitnessMatrix_1 * (j - 1)] = (shared->GeneFitness[i - 1] * shared->Genotypes[shared->dim_Genotypes_1 * (j - 1) + i - 1]);
    }
  }
  for (int i = 1; i <= shared->dim_Pop; ++i) {
    shared->initial_Pop[i - 1] = static_cast<int>(shared->pop_size * shared->Pop_ini_norm[i - 1]);
  }
  for (int i = 1; i <= shared->dim_GenotypeFitness; ++i) {
    shared->GenotypeFitness[i - 1] = odin_sum2<real_type>(shared->FitnessMatrix.data(), 0, shared->dim_FitnessMatrix_1, i - 1, i, shared->dim_FitnessMatrix_1);
  }
  return dust::pars_type<WrightFisher_nGenotypes_haploid>(shared, internal);
}
template <>
cpp11::sexp dust_info<WrightFisher_nGenotypes_haploid>(const dust::pars_type<WrightFisher_nGenotypes_haploid>& pars) {
  const std::shared_ptr<const WrightFisher_nGenotypes_haploid::shared_type> shared = pars.shared;
  cpp11::writable::strings nms({"time", "Pop"});
  cpp11::writable::list dim(2);
  dim[0] = cpp11::writable::integers({1});
  dim[1] = cpp11::writable::integers({shared->dim_Pop});
  dim.names() = nms;
  cpp11::writable::list index(2);
  index[0] = cpp11::writable::integers({1});
  index[1] = integer_sequence(2, shared->dim_Pop);
  index.names() = nms;
  size_t len = 1 + shared->dim_Pop;
  using namespace cpp11::literals;
  return cpp11::writable::list({
           "dim"_nm = dim,
           "len"_nm = len,
           "index"_nm = index});
}
}

cpp11::sexp dust_WrightFisher_nGenotypes_haploid_capabilities() {
  return dust::r::dust_capabilities<WrightFisher_nGenotypes_haploid>();
}

cpp11::sexp dust_WrightFisher_nGenotypes_haploid_gpu_info() {
  return dust::gpu::r::gpu_info();
}
using model_cpu = dust::dust_cpu<WrightFisher_nGenotypes_haploid>;

SEXP dust_cpu_WrightFisher_nGenotypes_haploid_alloc(cpp11::list r_pars, bool pars_multi, cpp11::sexp r_time,
                             cpp11::sexp r_n_particles, int n_threads,
                             cpp11::sexp r_seed, bool deterministic,
                             cpp11::sexp gpu_config, cpp11::sexp ode_control) {
  return dust::r::dust_cpu_alloc<WrightFisher_nGenotypes_haploid>(r_pars, pars_multi, r_time, r_n_particles,
                                        n_threads, r_seed, deterministic,
                                        gpu_config, ode_control);
}

SEXP dust_cpu_WrightFisher_nGenotypes_haploid_run(SEXP ptr, cpp11::sexp r_time_end) {
  return dust::r::dust_run<model_cpu>(ptr, r_time_end);
}

SEXP dust_cpu_WrightFisher_nGenotypes_haploid_simulate(SEXP ptr, cpp11::sexp r_time_end) {
  return dust::r::dust_simulate<model_cpu>(ptr, r_time_end);
}

SEXP dust_cpu_WrightFisher_nGenotypes_haploid_set_index(SEXP ptr, cpp11::sexp r_index) {
  dust::r::dust_set_index<model_cpu>(ptr, r_index);
  return R_NilValue;
}

SEXP dust_cpu_WrightFisher_nGenotypes_haploid_update_state(SEXP ptr, SEXP r_pars, SEXP r_state,
                                           SEXP r_time, SEXP r_set_initial_state, SEXP index, SEXP reset_step_size) {
  return dust::r::dust_update_state<model_cpu>(ptr, r_pars, r_state, r_time,
                                                      r_set_initial_state, index, reset_step_size);
}

SEXP dust_cpu_WrightFisher_nGenotypes_haploid_state(SEXP ptr, SEXP r_index) {
  return dust::r::dust_state<model_cpu>(ptr, r_index);
}

SEXP dust_cpu_WrightFisher_nGenotypes_haploid_time(SEXP ptr) {
  return dust::r::dust_time<model_cpu>(ptr);
}

void dust_cpu_WrightFisher_nGenotypes_haploid_reorder(SEXP ptr, cpp11::sexp r_index) {
  return dust::r::dust_reorder<model_cpu>(ptr, r_index);
}

SEXP dust_cpu_WrightFisher_nGenotypes_haploid_resample(SEXP ptr, cpp11::doubles r_weights) {
  return dust::r::dust_resample<model_cpu>(ptr, r_weights);
}

SEXP dust_cpu_WrightFisher_nGenotypes_haploid_rng_state(SEXP ptr, bool first_only, bool last_only) {
  return dust::r::dust_rng_state<model_cpu>(ptr, first_only, last_only);
}

SEXP dust_cpu_WrightFisher_nGenotypes_haploid_set_rng_state(SEXP ptr, cpp11::raws rng_state) {
  dust::r::dust_set_rng_state<model_cpu>(ptr, rng_state);
  return R_NilValue;
}

SEXP dust_cpu_WrightFisher_nGenotypes_haploid_set_data(SEXP ptr, cpp11::list data,
                                       bool shared) {
  dust::r::dust_set_data<model_cpu>(ptr, data, shared);
  return R_NilValue;
}

SEXP dust_cpu_WrightFisher_nGenotypes_haploid_compare_data(SEXP ptr) {
  return dust::r::dust_compare_data<model_cpu>(ptr);
}

SEXP dust_cpu_WrightFisher_nGenotypes_haploid_filter(SEXP ptr, SEXP time_end,
                                     bool save_trajectories,
                                     cpp11::sexp time_snapshot,
                                     cpp11::sexp min_log_likelihood) {
  return dust::r::dust_filter<model_cpu>(ptr, time_end,
                                                save_trajectories,
                                                time_snapshot,
                                                min_log_likelihood);
}

void dust_cpu_WrightFisher_nGenotypes_haploid_set_n_threads(SEXP ptr, int n_threads) {
  return dust::r::dust_set_n_threads<model_cpu>(ptr, n_threads);
}

int dust_cpu_WrightFisher_nGenotypes_haploid_n_state(SEXP ptr) {
  return dust::r::dust_n_state<model_cpu>(ptr);
}

void dust_cpu_WrightFisher_nGenotypes_haploid_set_stochastic_schedule(SEXP ptr, SEXP time) {
  dust::r::dust_set_stochastic_schedule<model_cpu>(ptr, time);
}

void dust_cpu_WrightFisher_nGenotypes_haploid_ode_statistics(SEXP ptr) {
  dust::r::dust_ode_statistics<model_cpu>(ptr);
}
