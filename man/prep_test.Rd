% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustOpt.R
\name{prep_test}
\alias{prep_test}
\title{prep_test}
\usage{
prep_test(input, subject_ids, dtype = "scRNA", test_id)
}
\arguments{
\item{input}{Seurat object}

\item{subject_ids}{Metadata field that identifies unique subjects.}

\item{dtype}{Type of data in the Seurat object "scRNA" or "CyTOF", default
is "scRNA". CyTOF data is expected to be arcsinh normalized.}

\item{test_id}{subject_id for the test sample}
}
\value{
Training data formatted for sil_score, format depends on dtype
}
\description{
Prepare test data for random forest
}
