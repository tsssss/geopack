Most routines in this folder deal with arrays of data in multi-dimensions. The default format is [n, m1, m2, ...], where n is the number of records of array, [m1, m2, ...] is the dimension of data at each record. For example, an array of 3-d electric field is in [n, 3], an array of density is in [n].

The reason for choosing [n, m1, m2, ...] format as default follows: (1) the 1st dimension is always the record; (2) [*,m] will appreciate idl's strength of fast array calculation.
