#!/bin/bash

# TODO: command line inputs
one=work/tc2/symmetric/chksum_diag
two=work/tc2/asymmetric/chksum_diag

for chk in ${one} ${two}; do
    awk '{print $(NF-2) " " $(NF-1) " " $(NF),$0}' ${chk} | sort > ${chk}.sorted
done

cmp ${one}.sorted ${two}.sorted

#diff <(sort ${one}) <(sort ${two})

#diff <(sort $(awk '{print $(NF-2) " " $(NF-1) " " $(NF),$0}' ${one})) <(sort $(awk '{print $(NF-2) " " $(NF-1) " " $(NF),$0}' ${two}))
