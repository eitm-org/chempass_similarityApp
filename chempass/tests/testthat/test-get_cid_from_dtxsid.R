# testing get_cid_from_dtxsid

test_that("get_cid_from_dtxsid works with valid DTXSID", {
  # Example: DTXSID70860951 maps to CID 5363 (benzo[a]pyrene)
  result <- get_cid_from_dtxsid("DTXSID70860951")
  expect_type(result, "character")
  expect_false(result == "error")
  expect_true(grepl("^[0-9]+$", result))
})

test_that("get_cid_from_dtxsid returns error for invalid DTXSID", {
  result <- get_cid_from_dtxsid("DTXSID9999999")  # Likely invalid
  expect_equal(result, "error")
})

test_that("get_cid_from_dtxsid returns error for empty string", {
  result <- get_cid_from_dtxsid("")
  expect_equal(result, "error")
})

test_that("get_cid_from_dtxsid returns error for non-string input", {
  # Note: Python function expects a string; in R, send something incorrect
  result <- get_cid_from_dtxsid(12345)
  expect_equal(result, "error")
})


