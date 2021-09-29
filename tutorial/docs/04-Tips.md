# Helpful Tips

## Handle	nulls	in	categorical	annotation
Such nan’s in categorical annotation would cause trouble in VIP because it cannot be converted to string. Here is how to handle it, let’s call the annotation X_annotation:

PY CODE HERE

## Display	full	traceback	stack	for	debugging	in	VIP
It follows the global setting. Please set “—verbose” to launch cellxgene server. 

## Pitfall	of	using	special	characters
In the mode which allows user to create manual annotation in cellxgene, user should try to avoid using hyphen (“-”) in name label.It would cause client-side issue. Please try to use underscores.

## Potential	use	for	bulk	or	pseudo	bulk	sample	dataset
Once the data matrix is replaced by sample x gene matrix, cellxgene VIP framework can handle regular bulk / pseudobulk RNAseq datasets. Simply replace “cells” by “samples”. All plotting functions would still work.
