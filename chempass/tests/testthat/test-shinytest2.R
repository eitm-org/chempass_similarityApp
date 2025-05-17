library(shinytest2)

test_that("{shinytest2} recording: cid_input_test", {
  app <- AppDriver$new(variant = platform_variant(), name = "cid_input_test", height = 981, 
      width = 1409)
  rlang::warn(paste0("`file_upload` should be the path to the file, relative to the app's tests/testthat directory.\n", 
      "Remove this warning when the file is in the correct location."))
  app$upload_file(file_upload = "CID_test_file1.csv")
  app$click("process_file")
  app$set_inputs(fingerprint_type = "FCFP4")
  app$click("fingerprint_button")
  app$set_inputs(cutoff = 0)
  app$set_inputs(cutoff = 0.6)
  app$click("cluster")
  app$set_inputs(`plotly_hover-A` = "[{\"curveNumber\":2,\"pointNumber\":2,\"x\":-0.15295728129252464,\"y\":0.1633038703657812}]", 
      allow_no_input_binding_ = TRUE, priority_ = "event")
  app$set_inputs(`plotly_hover-A` = character(0), allow_no_input_binding_ = TRUE, 
      priority_ = "event")
  app$set_inputs(`plotly_hover-A` = "[{\"curveNumber\":3,\"pointNumber\":0,\"x\":0.40857599691199487,\"y\":0.15319323587898168}]", 
      allow_no_input_binding_ = TRUE, priority_ = "event")
  app$set_inputs(`plotly_hover-A` = character(0), allow_no_input_binding_ = TRUE, 
      priority_ = "event")
  app$set_inputs(`plotly_hover-A` = "[{\"curveNumber\":3,\"pointNumber\":1,\"x\":0.2867178101269046,\"y\":0.1906331500457626}]", 
      allow_no_input_binding_ = TRUE, priority_ = "event")
  app$set_inputs(`plotly_hover-A` = character(0), allow_no_input_binding_ = TRUE, 
      priority_ = "event")
  app$set_inputs(`plotly_hover-A` = "[{\"curveNumber\":3,\"pointNumber\":2,\"x\":0.3454651344110293,\"y\":0.055251335682748974}]", 
      allow_no_input_binding_ = TRUE, priority_ = "event")
  app$set_inputs(`plotly_hover-A` = character(0), allow_no_input_binding_ = TRUE, 
      priority_ = "event")
  app$set_inputs(`plotly_hover-A` = "[{\"curveNumber\":0,\"pointNumber\":22,\"x\":0.3454651344110293,\"y\":0.055251335682748974}]", 
      allow_no_input_binding_ = TRUE, priority_ = "event")
  app$set_inputs(`plotly_hover-A` = character(0), allow_no_input_binding_ = TRUE, 
      priority_ = "event")
  app$set_inputs(`plotly_hover-A` = "[{\"curveNumber\":3,\"pointNumber\":2,\"x\":0.3454651344110293,\"y\":0.055251335682748974}]", 
      allow_no_input_binding_ = TRUE, priority_ = "event")
  app$set_inputs(`plotly_hover-A` = character(0), allow_no_input_binding_ = TRUE, 
      priority_ = "event")
  app$expect_values()
  app$expect_screenshot()
})
