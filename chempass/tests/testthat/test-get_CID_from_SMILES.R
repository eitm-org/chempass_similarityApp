# get_CID_from_SMILES

test_that("get_CID_from_SMILES works with valid SMILES", {
  # SMILES for Aspirin: CC(=O)OC1=CC=CC=C1C(=O)O
  encoded_smiles <- base64encode(charToRaw("CC(=O)OC1=CC=CC=C1C(=O)O"))
  result <- get_CID_from_SMILES(encoded_smiles)
  
  expect_type(result, "character")
  expect_false(result == "error")
  expect_equal(result, "2244")  # Known CID for aspirin
})

test_that("get_CID_from_SMILES returns error for invalid SMILES", {
  encoded_smiles <- base64encode(charToRaw("%%%"))  # Invalid SMILES
  result <- get_CID_from_SMILES(encoded_smiles)
  
  expect_equal(result, "error")
})

test_that("get_CID_from_SMILES returns error for malformed base64", {
  result <- get_CID_from_SMILES("not_base64_text")
  expect_equal(result, "error")
})

test_that("get_CID_from_SMILES returns error for empty input", {
  result <- get_CID_from_SMILES("")
  expect_equal(result, "error")
})