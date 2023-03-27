WrightFisher_nGenotypes_haploid <- R6::R6Class(
  "dust",
  cloneable = FALSE,

  private = list(
    pars_ = NULL,
    pars_multi_ = NULL,
    index_ = NULL,
    info_ = NULL,
    n_threads_ = NULL,
    n_particles_ = NULL,
    n_particles_each_ = NULL,
    shape_ = NULL,
    ptr_ = NULL,
    gpu_config_ = NULL,
    ode_control_ = NULL,
    methods_ = NULL,
    param_ = list(gene_no = list(has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE),
     GeneFitness = list(has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE),
     Genotypes = list(has_default = FALSE, default_value = NULL, rank = 2, min = -Inf, max = Inf, integer = FALSE),
     Pop_ini = list(has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE),
     pop_size = list(has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE),
     species_no = list(has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE),
     dt = list(has_default = TRUE, default_value = 1L, rank = 0, min = -Inf, max = Inf, integer = FALSE)),
    reload_ = list(path = "build", base = "dust2252a231")
  ),

  public = list(
    initialize = function(pars, time, n_particles, n_threads = 1L,
                          seed = NULL, pars_multi = FALSE,
                          deterministic = FALSE,
                          gpu_config = NULL, ode_control = NULL) {
      if (is.null(gpu_config)) {
        private$methods_ <- list(
           alloc = dust_cpu_WrightFisher_nGenotypes_haploid_alloc,
           run = dust_cpu_WrightFisher_nGenotypes_haploid_run,
           simulate = dust_cpu_WrightFisher_nGenotypes_haploid_simulate,
           set_index = dust_cpu_WrightFisher_nGenotypes_haploid_set_index,
           n_state = dust_cpu_WrightFisher_nGenotypes_haploid_n_state,
           update_state = dust_cpu_WrightFisher_nGenotypes_haploid_update_state,
           state = dust_cpu_WrightFisher_nGenotypes_haploid_state,
           time = dust_cpu_WrightFisher_nGenotypes_haploid_time,
           reorder = dust_cpu_WrightFisher_nGenotypes_haploid_reorder,
           resample = dust_cpu_WrightFisher_nGenotypes_haploid_resample,
           rng_state = dust_cpu_WrightFisher_nGenotypes_haploid_rng_state,
           set_rng_state = dust_cpu_WrightFisher_nGenotypes_haploid_set_rng_state,
           set_n_threads = dust_cpu_WrightFisher_nGenotypes_haploid_set_n_threads,
           set_data = dust_cpu_WrightFisher_nGenotypes_haploid_set_data,
           compare_data = dust_cpu_WrightFisher_nGenotypes_haploid_compare_data,
           filter = dust_cpu_WrightFisher_nGenotypes_haploid_filter,
           set_stochastic_schedule = dust_cpu_WrightFisher_nGenotypes_haploid_set_stochastic_schedule,
           ode_statistics = dust_cpu_WrightFisher_nGenotypes_haploid_ode_statistics)
      } else {
        private$methods_ <- list(alloc = function(...) {
          stop("GPU support not enabled for this object")
        })
      }
      res <- private$methods_$alloc(pars, pars_multi, time, n_particles,
                        n_threads, seed, deterministic, gpu_config, ode_control)
      private$pars_ <- pars
      private$pars_multi_ <- pars_multi
      private$n_threads_ <- n_threads
      private$ptr_ <- res[[1L]]
      private$info_ <- res[[2L]]
      private$shape_ <- res[[3L]]
      private$gpu_config_ <- res[[4L]]
      private$ode_control_ <- res[[5L]]
      private$n_particles_ <- prod(private$shape_)
      if (pars_multi) {
        private$n_particles_each_ <- private$n_particles_ / length(pars)
      } else {
        private$n_particles_each_ <- private$n_particles_
      }
    },

    name = function() {
      "WrightFisher_nGenotypes_haploid"
    },

    param = function() {
      private$param_
    },

    run = function(time_end) {
      m <- private$methods_$run(private$ptr_, time_end)
      rownames(m) <- names(private$index_)
      m
    },

    simulate = function(time_end) {
      m <- private$methods_$simulate(private$ptr_, time_end)
      rownames(m) <- names(private$index_)
      m
    },

    set_index = function(index) {
      private$methods_$set_index(private$ptr_, index)
      private$index_ <- index
      invisible()
    },

    index = function() {
      private$index_
    },

    ode_control = function() {
      private$ode_control_
    },

    ode_statistics = function() {
      private$methods_$ode_statistics(private$ptr_)
    },

    n_threads = function() {
      private$n_threads_
    },

    n_state = function() {
      private$methods_$n_state(private$ptr_)
    },

    n_particles = function() {
      private$n_particles_
    },

    n_particles_each = function() {
      private$n_particles_each_
    },

    shape = function() {
      private$shape_
    },

    update_state = function(pars = NULL, state = NULL, time = NULL,
                            set_initial_state = NULL, index = NULL,
                            reset_step_size = NULL) {
      info <- private$methods_$update_state(private$ptr_, pars, state, time,
                                          set_initial_state, index,
                                          reset_step_size)
      if (!is.null(pars)) {
        private$info_ <- info
        private$pars_ <- pars
      }
      invisible()
    },

    state = function(index = NULL) {
      m <- private$methods_$state(private$ptr_, index)
      rownames(m) <- names(index)
      m
    },

    time = function() {
      private$methods_$time(private$ptr_)
    },

    set_stochastic_schedule = function(time) {
      private$methods_$set_stochastic_schedule(private$ptr, time)
      invisible()
    },

    reorder = function(index) {
      storage.mode(index) <- "integer"
      private$methods_$reorder(private$ptr_, index)
      invisible()
    },

    resample = function(weights) {
      invisible(private$methods_$resample(private$ptr_, weights))
    },

    info = function() {
      private$info_
    },

    pars = function() {
      private$pars_
    },

    rng_state = function(first_only = FALSE, last_only = FALSE) {
      private$methods_$rng_state(private$ptr_, first_only, last_only)
    },

    set_rng_state = function(rng_state) {
      private$methods_$set_rng_state(private$ptr_, rng_state)
      invisible()
    },

    has_openmp = function() {
      dust_WrightFisher_nGenotypes_haploid_capabilities()[["openmp"]]
    },

    has_gpu_support = function(fake_gpu = FALSE) {
      if (fake_gpu) {
        FALSE
      } else {
        dust_WrightFisher_nGenotypes_haploid_capabilities()[["gpu"]]
      }
    },

    has_compare = function() {
      dust_WrightFisher_nGenotypes_haploid_capabilities()[["compare"]]
    },

    real_size = function() {
      dust_WrightFisher_nGenotypes_haploid_capabilities()[["real_size"]]
    },

    rng_algorithm = function() {
      dust_WrightFisher_nGenotypes_haploid_capabilities()[["rng_algorithm"]]
    },

    uses_gpu = function(fake_gpu = FALSE) {
      real_gpu <- private$gpu_config_$real_gpu
      !is.null(real_gpu) && (fake_gpu || real_gpu)
    },

    n_pars = function() {
      if (private$pars_multi_) length(private$pars_) else 0L
    },

    set_n_threads = function(n_threads) {
      prev <- private$n_threads_
      private$methods_$set_n_threads(private$ptr_, n_threads)
      private$n_threads_ <- n_threads
      invisible(prev)
    },

    set_data = function(data, shared = FALSE) {
      private$methods_$set_data(private$ptr_, data, shared)
    },

    compare_data = function() {
      private$methods_$compare_data(private$ptr_)
    },

    filter = function(time_end = NULL, save_trajectories = FALSE,
                      time_snapshot = NULL, min_log_likelihood = NULL) {
      private$methods_$filter(private$ptr_, time_end, save_trajectories,
                              time_snapshot, min_log_likelihood)
    },

    gpu_info = function() {
      ret <- dust_WrightFisher_nGenotypes_haploid_gpu_info()
      parent <- parent.env(environment())
      if (ret$has_cuda && exists("private", parent, inherits = FALSE)) {
        ret$config <- private$gpu_config_
      }
      ret
    }
  ))
class(WrightFisher_nGenotypes_haploid) <- c("dust_generator", class(WrightFisher_nGenotypes_haploid))
