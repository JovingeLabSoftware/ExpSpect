library(stringr)
context("CMAP")

test_that("Values match those from CMAP website", {
  expect_vecs_equal <- function(a, b) {
    a <- as.vector(a)
    b <- as.vector(b)
    eval(bquote(expect_true(all.equal(.(a), .(b), tolerance=0.00001))))
  }
  data("cmap_test")
  ## complete test if full data matrix available (too big to include in package):
  # load("/mnt/lincs/cmap_rank_matrix.Rdata")
  # scores <- exps$cmap(cmap, list(up=up, down=down))
  ## all the following should be true
  # round(scores[4]*100000)/100000 == -0.39704
  # round(scores[7]*100000)/100000 == -0.31993
  # round(scores[10]*100000)/100000 == -0.44053
  # sum(scores[1:10] < 0) == 4
  # sum(scores[1:10] == 0) == 6

  exps <- ExpSpect$new()
  scores <- exps$cmap(cmap_test, list(up=up, down=down))
  expect_vecs_equal(scores, c(0.0, 0.0, 0.0, -0.8086089, -1.0, 0.0, -0.6515715, 0.0, 0.0, -0.8971877))
  scores <- exps$cmap(cmap_test, cmap_test[,1], threshold=1.5)
  expect_vecs_equal(scores, c(-1, -0.44022, -0.45385, -0.49065, 0.93272, 1, 0.90563, -0.12315, 0.685, 0))
})
          
          
