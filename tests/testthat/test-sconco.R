test_that("the shipped databases and gene positions load", {
  db <- sconco_load_db()
  expect_s3_class(db, "data.frame")
  expect_true(all(c("gene", "Level1", "Level4_Abb") %in% names(db)))
  expect_true("state" %in% db$part)                       # v3.3: glioblastoma state modules
  db32 <- sconco_load_db(system.file("extdata", "DB_pancancer_human_v3.2.R", package = "scONCO"))
  expect_false("state" %in% db32$part)
  gp <- sconco_load_gene_positions()
  expect_true(all(c("gene", "chrom", "start", "end") %in% names(gp)))
  expect_gt(nrow(gp), 20000)
})

test_that("the tumour test is negative on noise and positive on a shifted clone", {
  set.seed(1)
  ref <- paste0("r", 1:200); cand <- paste0("c", 1:200)
  sc <- data.frame(cell = c(ref, cand), burden = rexp(400), coherence = 0, row.names = c(ref, cand))
  expect_false(sconco_tumour_test(sc, ref, cand)$present)
  sc[cand[1:80], "burden"] <- sc[cand[1:80], "burden"] + 20
  tt <- sconco_tumour_test(sc, ref, cand)
  expect_true(tt$present)
  expect_equal(tt$frac_above_null, 0.4, tolerance = 0.05)
})
