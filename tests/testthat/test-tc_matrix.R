context("tc_matrix")

test_that("a simple tc matrix is computed correctly", {
   tcs <- SEAR::tc_matrix(
	ref_fp=data.frame(
		smiles=c(
			"c1cc(c(cc1CC[NH3+])O)O",
			"c1cc2c(cc1O)c(c[nH]2)CC[NH3+]",
			"C[NH2+]C[C@@H](c1ccc(c(c1)O)O)O"),
		compound=c(
			"dopamine",
			"serotonin",
			"adrenaline")),
	query_fp=data.frame(
		smiles=c(
			"COc1cc(ccc1O)CC(=O)[O-]",
			"c1cc(c(cc1C[C@@H](C(=O)[O-])[NH3+])O)O"),
		compound=c(
			"homovanillic acid",
			"levodopa")))
   expect_equal(
       tcs,
       structure(
           list(
               ref_cid = c(
                   "serotonin", "serotonin",
                   "adrenaline", "adrenaline",
                   "dopamine", "dopamine"),
               query_cid = c(
                   "levodopa", "homovanillic acid",
                   "levodopa", "homovanillic acid",
                   "levodopa", "homovanillic acid"),
               tc = c(
                   0.166666666667, 0.1875, 0.292682926829, 
                   0.227272727273, 0.454545454545, 0.324324324324)),
           .Names = c("ref_cid", "query_cid", "tc"),
           class = c("tbl_df", "tbl", "data.frame"),
           row.names = c(NA, -6L),
           spec = structure(list(cols = structure(list(
                 ref_cid = structure(list(), class = c("collector_character", "collector")),
                 query_cid = structure(list(), class = c("collector_character", "collector")),
                 tc = structure(list(), class = c("collector_double", "collector"))),
              .Names = c("ref_cid", "query_cid", "tc")),
            default = structure(list(), class = c("collector_guess", "collector"))),
           .Names = c("cols", "default"), class = "col_spec")))
 }
