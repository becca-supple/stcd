test_that("dsgn class properly applied", {
  expect_equal(class(dsgn("a")), c("dsgn", "character"))
})
