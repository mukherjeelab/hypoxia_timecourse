# Null test for SHAP-profile clustering, shared by 03_shap.Rmd and 03b_shap_subgroup_null.Rmd.
#
# k-means always returns clusters and the silhouette is always positive for some k, so a
# silhouette means nothing without a null. The null permutes each family INDEPENDENTLY across
# genes: co-occurrence between families is destroyed, every family's marginal distribution is
# preserved.
#
# THE STATISTIC MATTERS MORE THAN THE NULL. Mean silhouette falls monotonically with k, so the
# null's best-over-k is almost always at k = 2. Comparing that against an observed best that
# happens to sit at k = 3 compares two different quantities, and it gets the answer backwards:
# on this data it returned p = 0.74 (no structure) when every k >= 3 in fact exceeds the null's
# MAXIMUM over 200 permutations. The comparison must be made WITHIN k, and only then maximised
# over k to pay for having been free to choose k. Null moments are leave-one-out so a
# permutation is not standardised against itself.
subgroup_silhouette <- function(m, k, seed = 9) {
  set.seed(seed)
  km <- stats::kmeans(m, centers = k, nstart = 25)
  if (length(unique(km$cluster)) < 2) return(NA_real_)
  mean(cluster::silhouette(km$cluster, stats::dist(m))[, 3])
}

subgroup_null <- function(profile, ks = 2:6, n_perm = 200, seed0 = 1000) {
  obs <- vapply(ks, function(k) subgroup_silhouette(profile, k), numeric(1))

  # The permutation is drawn ONCE per replicate, before the k loop. subgroup_silhouette() calls
  # set.seed() internally, so drawing inside the loop would reset the stream and make every
  # replicate identical from the second k onwards - a null with zero variance.
  null_mat <- t(vapply(seq_len(n_perm), function(i) {
    set.seed(seed0 + i)
    perm <- apply(profile, 2, sample)
    vapply(ks, function(k) subgroup_silhouette(perm, k), numeric(1))
  }, numeric(length(ks))))
  colnames(null_mat) <- as.character(ks)

  mu  <- colMeans(null_mat)
  sdv <- apply(null_mat, 2, stats::sd)
  stopifnot("null has zero variance - the permutation is not being redrawn" = all(sdv > 0))

  per_k <- data.frame(
    k          = ks,
    observed   = obs,
    null_median = apply(null_mat, 2, stats::median),
    null_min   = apply(null_mat, 2, min),
    null_max   = apply(null_mat, 2, max),
    z          = (obs - mu) / sdv,
    p_within_k = vapply(seq_along(ks),
                        function(j) (sum(null_mat[, j] >= obs[j]) + 1) / (n_perm + 1),
                        numeric(1))
  )

  z_obs  <- max(per_k$z, na.rm = TRUE)
  z_null <- vapply(seq_len(n_perm), function(i) {
    m2 <- null_mat[-i, , drop = FALSE]
    max((null_mat[i, ] - colMeans(m2)) / apply(m2, 2, stats::sd), na.rm = TRUE)
  }, numeric(1))
  p_overall <- (sum(z_null >= z_obs) + 1) / (n_perm + 1)

  # Two different "best k", deliberately kept apart.
  #   best_k_z   - where the standardised statistic peaks. Correct for TESTING, but it favours
  #                the k at which the null happens to be tightest (here k = 6), which is a
  #                property of the null's variance, not of the structure.
  #   best_k_sil - where the raw silhouette peaks. What to DESCRIBE clusters at.
  list(ks = ks, obs = obs, null_mat = null_mat, per_k = per_k,
       best_k_z   = ks[which.max(per_k$z)],
       best_k_sil = ks[which.max(obs)],
       z_obs = z_obs, z_null = z_null, p_overall = p_overall, n_perm = n_perm)
}
