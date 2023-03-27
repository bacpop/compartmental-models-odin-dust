#include <cpp11.hpp>

[[cpp11::register]]
cpp11::sexp dust_WrightFisher_nGenotypes_haploid_capabilities();

[[cpp11::register]]
cpp11::sexp dust_WrightFisher_nGenotypes_haploid_gpu_info();
[[cpp11::register]]
SEXP dust_cpu_WrightFisher_nGenotypes_haploid_alloc(cpp11::list r_pars, bool pars_multi, cpp11::sexp r_time,
                         cpp11::sexp r_n_particles, int n_threads,
                         cpp11::sexp r_seed, bool deterministic,
                         cpp11::sexp gpu_config, cpp11::sexp ode_control);

[[cpp11::register]]
SEXP dust_cpu_WrightFisher_nGenotypes_haploid_run(SEXP ptr, cpp11::sexp r_time_end);

[[cpp11::register]]
SEXP dust_cpu_WrightFisher_nGenotypes_haploid_simulate(SEXP ptr, cpp11::sexp time_end);

[[cpp11::register]]
SEXP dust_cpu_WrightFisher_nGenotypes_haploid_set_index(SEXP ptr, cpp11::sexp r_index);

[[cpp11::register]]
SEXP dust_cpu_WrightFisher_nGenotypes_haploid_update_state(SEXP ptr, SEXP r_pars, SEXP r_state,
                                           SEXP r_time, SEXP r_set_initial_state,
                                           SEXP index, SEXP reset_step_size);

[[cpp11::register]]
SEXP dust_cpu_WrightFisher_nGenotypes_haploid_state(SEXP ptr, SEXP r_index);

[[cpp11::register]]
SEXP dust_cpu_WrightFisher_nGenotypes_haploid_time(SEXP ptr);

[[cpp11::register]]
void dust_cpu_WrightFisher_nGenotypes_haploid_reorder(SEXP ptr, cpp11::sexp r_index);

[[cpp11::register]]
SEXP dust_cpu_WrightFisher_nGenotypes_haploid_resample(SEXP ptr, cpp11::doubles r_weights);

[[cpp11::register]]
SEXP dust_cpu_WrightFisher_nGenotypes_haploid_rng_state(SEXP ptr, bool first_only, bool last_only);

[[cpp11::register]]
SEXP dust_cpu_WrightFisher_nGenotypes_haploid_set_rng_state(SEXP ptr, cpp11::raws rng_state);

[[cpp11::register]]
SEXP dust_cpu_WrightFisher_nGenotypes_haploid_set_data(SEXP ptr, cpp11::list data, bool shared);

[[cpp11::register]]
SEXP dust_cpu_WrightFisher_nGenotypes_haploid_compare_data(SEXP ptr);

[[cpp11::register]]
SEXP dust_cpu_WrightFisher_nGenotypes_haploid_filter(SEXP ptr, SEXP time_end,
                                     bool save_trajectories,
                                     cpp11::sexp time_snapshot,
                                     cpp11::sexp min_log_likelihood);

[[cpp11::register]]
void dust_cpu_WrightFisher_nGenotypes_haploid_set_n_threads(SEXP ptr, int n_threads);

[[cpp11::register]]
int dust_cpu_WrightFisher_nGenotypes_haploid_n_state(SEXP ptr);

[[cpp11::register]]
void dust_cpu_WrightFisher_nGenotypes_haploid_set_stochastic_schedule(SEXP ptr, SEXP time);

[[cpp11::register]]
void dust_cpu_WrightFisher_nGenotypes_haploid_ode_statistics(SEXP ptr);
